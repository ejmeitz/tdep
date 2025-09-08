#include "precompilerdefinitions"
module type_canonical_configs
!! Information about an MD simulation
use konstanter, only: r8, lo_pi, lo_huge, lo_hugeint, lo_sqtol, lo_status, &
                      lo_exitcode_param, lo_bohr_to_A, lo_Hartree_to_eV, &
                      lo_exitcode_io, lo_velocity_au_to_Afs
use gottochblandat, only: open_file, tochar, walltime 
use mpi_wrappers, only:  lo_mpi_helper, lo_stop_gracefully, MPI_SUM, MPI_INTEGER, MPI_INFO_NULL
use hdf5_wrappers, only: lo_hdf5_helper, lo_h5_store_attribute, lo_h5_store_data
use type_crystalstructure, only: lo_crystalstructure
use hdf5
use mpi

use iso_fortran_env ! REMOVE LATER

implicit none
private
public :: lo_canonical_configs

!> All the energies in an simulation
type lo_canonical_configs_stat
    !> kinetic energy
    real(r8), dimension(:), allocatable :: kinetic_energy
    !> dipole-dipole energy
    real(r8), dimension(:), allocatable :: polar_potential_energy
    !> harmonic potential energy
    real(r8), dimension(:), allocatable :: secondorder_potential_energy
    !> third order potential energy
    real(r8), dimension(:), allocatable :: thirdorder_potential_energy
    !> fourth order potential energy
    real(r8), dimension(:), allocatable :: fourthorder_potential_energy
end type

!> information from the reference starting position, such as unitcell and so on
type lo_canonical_configs_extra
    !> latticevectors
    real(r8), dimension(:, :), allocatable :: unitcell_latticevectors
    real(r8), dimension(:, :), allocatable :: supercell_latticevectors
    !> reference positions
    real(r8), dimension(:, :), allocatable :: unitcell_positions
    real(r8), dimension(:, :), allocatable :: supercell_positions
    !> atomic numbers
    integer, dimension(:), allocatable :: unitcell_atomic_numbers
    integer, dimension(:), allocatable :: supercell_atomic_numbers
    !> Alloy things, if needed
    integer, dimension(:), allocatable ::      unitcell_componentcounter
    integer, dimension(:), allocatable ::      supercell_componentcounter
    integer, dimension(:, :), allocatable ::    unitcell_components
    integer, dimension(:, :), allocatable ::    supercell_components
    real(r8), dimension(:, :), allocatable :: unitcell_concentrations
    real(r8), dimension(:, :), allocatable :: supercell_concentrations
end type

!> A set of displacement/velocity configs from the canonical configurations command
type lo_canonical_configs
    !> energies, stress temperatures
    type(lo_canonical_configs_stat) :: stat
    !> general information
    type(lo_canonical_configs_extra) :: extra
    !> reference lattice for simulation
    type(lo_crystalstructure) :: crystalstructure
    !> thermostat temperature
    real(r8) :: temperature_thermostat = -lo_huge
    !> Are these alloy configurations?
    logical :: alloy = .false.
    !> number of configurations
    integer :: nt = -lo_hugeint
    !> number of atoms
    integer :: na = -lo_hugeint
    !> positions in fractional coordinates
    real(r8), allocatable, dimension(:, :, :) :: r
    !> velocities in cartesian coordinates
    real(r8), allocatable, dimension(:, :, :) :: v
    !> atomic numbers (or species? Not sure. This is unambiguous at least.)
    integer, allocatable, dimension(:) :: atomic_numbers
contains
    !> allocate and create empty object
    procedure :: init_empty
    !> set values for specific timestep 
    procedure :: set_step
    !> write to hdf5
    procedure :: write_to_hdf5

    procedure :: write_hdf5_header

    procedure :: write_hdf5_mpi
end type

contains

!> Set a specific timestep 
subroutine set_step(cc, positions, velocities, kinetic_energy, ep, e2, e3, e4, idx)
    !> md simulation
    class(lo_canonical_configs), intent(inout) :: cc
    !> positions, in fractional coordinates
    real(r8), dimension(:, :), intent(in) :: positions
    !> velocities, in Cartesian coordinates
    real(r8), dimension(:, :), intent(in) :: velocities
    !> kinetic energy (per atom)
    real(r8), intent(in) :: kinetic_energy
    !> polar component of potential energy
    real(r8), intent(in) :: ep
    !> harmonic component of potential energy
    real(r8), intent(in) :: e2
    !> third-order component of potential energy
    real(r8), intent(in) :: e3
    !> fourth-order component of potential energy
    real(r8), intent(in) :: e4
    !> index to set in storage
    integer, intent(in) :: idx

    integer :: tmax
    character(len=256) :: msg

    ! Max number of timesteps
    tmax = size(cc%r, 3)

    ! Sanity tests
    if ((idx .lt. 1) .or. (idx .gt. tmax)) then
        write(msg,'(A," idx=",I0,", tmax=",I0," (valid: 1..tmax)")') &
        'Not enough space to store timestep.', idx, tmax
        call lo_stop_gracefully([trim(msg)], lo_exitcode_param, __FILE__, __LINE__)
    end if

    ! Start storing things
    cc%r(:, :, idx) = positions
    cc%v(:, :, idx) = velocities

    ! Store energies and stuff
    cc%stat%kinetic_energy(idx) = kinetic_energy
    cc%stat%polar_potential_energy(idx) = ep
    cc%stat%secondorder_potential_energy(idx) = e2
    cc%stat%thirdorder_potential_energy(idx) = e3
    cc%stat%fourthorder_potential_energy(idx) = e4

end subroutine

!> create empty sim container, to be filled incrementally
subroutine init_empty(cc, uc, ss, nstep, temperature)
    !> md simulation
    class(lo_canonical_configs), intent(out) :: cc
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> how many steps
    integer, intent(in) :: nstep
    !> specify temperature
    real(r8), intent(in) :: temperature

    ! Set metadata
    init: block
        ! Some sanity tests
        if (ss%info%supercell .eqv. .false.) then
            call lo_stop_gracefully(['Need a supercell/unitcell pair to generate empty simulation'], lo_exitcode_param, __FILE__, __LINE__)
        end if
        if (nstep .lt. 1) then
            call lo_stop_gracefully(['Need at least one timestep to generate empty simulation'], lo_exitcode_param, __FILE__, __LINE__)
        end if

        ! Set some basic things
        cc%crystalstructure = ss
        cc%na = ss%na
        cc%nt = nstep
        cc%temperature_thermostat = temperature

        if (uc%info%alloy) then
            cc%alloy = .true.
        else
            cc%alloy = .false.
        end if

    end block init

    ! Info on the perfect structure, will come handy later
    structure: block
        integer :: i, j, k, l

        lo_allocate(cc%extra%unitcell_positions(3, uc%na))
        lo_allocate(cc%extra%supercell_positions(3, ss%na))
        lo_allocate(cc%extra%unitcell_atomic_numbers(uc%na))
        lo_allocate(cc%extra%supercell_atomic_numbers(ss%na))
        lo_allocate(cc%extra%unitcell_latticevectors(3, 3))
        lo_allocate(cc%extra%supercell_latticevectors(3, 3))
        cc%extra%unitcell_positions = uc%r
        cc%extra%supercell_positions = ss%r
        cc%extra%unitcell_atomic_numbers = uc%atomic_number
        cc%extra%supercell_atomic_numbers = ss%atomic_number
        cc%extra%unitcell_latticevectors = uc%latticevectors
        cc%extra%supercell_latticevectors = ss%latticevectors

        ! Also, if alloy store a lot of extra things
        if (cc%alloy) then
            ! Figure out the max number of components
            l = 0
            do i = 1, uc%na
                l = max(l, uc%alloyspecies(uc%species(i))%n)
            end do
            ! Space for alloy specification
            lo_allocate(cc%extra%unitcell_componentcounter(uc%na))
            lo_allocate(cc%extra%supercell_componentcounter(ss%na))
            lo_allocate(cc%extra%unitcell_components(l, uc%na))
            lo_allocate(cc%extra%supercell_components(l, ss%na))
            lo_allocate(cc%extra%unitcell_concentrations(l, uc%na))
            lo_allocate(cc%extra%supercell_concentrations(l, ss%na))
            cc%extra%unitcell_componentcounter = -1
            cc%extra%supercell_componentcounter = -1
            cc%extra%unitcell_components = -1
            cc%extra%supercell_components = -1
            cc%extra%unitcell_concentrations = 0.0_r8
            cc%extra%supercell_concentrations = 0.0_r8
            ! Store alloy specification
            do i = 1, uc%na
                j = uc%species(i)
                cc%extra%unitcell_componentcounter(i) = uc%alloyspecies(j)%n
                do k = 1, cc%extra%unitcell_componentcounter(i)
                    cc%extra%unitcell_components(k, i) = uc%alloyspecies(j)%atomic_number(k)
                    cc%extra%unitcell_concentrations(k, i) = uc%alloyspecies(j)%concentration(k)
                end do
            end do
            do i = 1, ss%na
                j = ss%species(i)
                cc%extra%supercell_componentcounter(i) = ss%alloyspecies(j)%n
                do k = 1, cc%extra%supercell_componentcounter(i)
                    cc%extra%supercell_components(k, i) = ss%alloyspecies(j)%atomic_number(k)
                    cc%extra%supercell_concentrations(k, i) = ss%alloyspecies(j)%concentration(k)
                end do
            end do
        else
            ! Allocate dummy empty arrays
            lo_allocate(cc%extra%unitcell_componentcounter(1))
            lo_allocate(cc%extra%supercell_componentcounter(1))
            lo_allocate(cc%extra%unitcell_components(1, 1))
            lo_allocate(cc%extra%supercell_components(1, 1))
            lo_allocate(cc%extra%unitcell_concentrations(1, 1))
            lo_allocate(cc%extra%supercell_concentrations(1, 1))
            cc%extra%unitcell_componentcounter = -lo_hugeint
            cc%extra%supercell_componentcounter = -lo_hugeint
            cc%extra%unitcell_components = -lo_hugeint
            cc%extra%supercell_components = -lo_hugeint
            cc%extra%unitcell_concentrations = -lo_huge
            cc%extra%supercell_concentrations = -lo_huge
        end if
  
    end block structure

    trajectories: block
        ! Space for trajectories
        lo_allocate(cc%r(3, cc%na, nstep))
        lo_allocate(cc%v(3, cc%na, nstep))
        cc%r = 0.0_r8
        cc%v = 0.0_r8

        if (cc%alloy) then
            lo_allocate(cc%atomic_numbers(cc%na))
            cc%atomic_numbers = 0
        end if
    end block trajectories

    ! And some space for energies
    energies: block
        lo_allocate(cc%stat%kinetic_energy(nstep))
        lo_allocate(cc%stat%polar_potential_energy(nstep))
        lo_allocate(cc%stat%secondorder_potential_energy(nstep))
        lo_allocate(cc%stat%thirdorder_potential_energy(nstep))
        lo_allocate(cc%stat%fourthorder_potential_energy(nstep))
        cc%stat%kinetic_energy = 0.0_r8
        cc%stat%polar_potential_energy = 0.0_r8
        cc%stat%secondorder_potential_energy = 0.0_r8
        cc%stat%thirdorder_potential_energy = 0.0_r8
        cc%stat%fourthorder_potential_energy = 0.0_r8
    end block energies
end subroutine

!> write a simulation to hdf5 (assumes single rank writing)
subroutine write_to_hdf5(cc, uc, ss, filename, verbosity)
    !> md simulation
    class(lo_canonical_configs), intent(in) :: cc
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(inout) :: ss
    !> filename
    character(len=*), intent(in) :: filename
    !> Talk a lot?
    integer, intent(in) :: verbosity

    real(r8) :: timer

    init: block
        if (verbosity .gt. 0) then
            timer = walltime()
            write (*, *) ''
            write (*, *) 'Writing simulation to "'//trim(filename)//'"'
        end if

        ! Check that the supercell really is a supercell
        call ss%classify('supercell', uc)

    end block init

    writefile: block
        type(lo_hdf5_helper) :: h5
        real(r8), dimension(:, :, :, :), allocatable :: dw
        real(r8), dimension(:, :, :), allocatable :: dr
        real(r8), dimension(:), allocatable :: ds
        ! Initialize hdf5 properly
        call h5%init(__FILE__, __LINE__)
        call h5%open_file('write', trim(filename))

        ! Store some metadata. Not sure if effective.
        call lo_h5_store_attribute(cc%na, h5%file_id, 'number_of_atoms')
        call lo_h5_store_attribute(cc%nt, h5%file_id, 'number_of_timesteps')
        call lo_h5_store_attribute(cc%temperature_thermostat, h5%file_id, 'temperature_thermostat')
        call lo_h5_store_attribute(cc%alloy, h5%file_id, 'is_simulation_alloy')

        ! Write positions
        lo_allocate(dr(3, cc%na, cc%nt))
        dr = cc%r(:, :, 1:cc%nt)
        call lo_h5_store_data(dr, h5%file_id, 'positions', enhet='fractional')
        if (verbosity .gt. 0) write (*, *) '... wrote positions'

        ! Write velocities
        dr = cc%v(:, :, 1:cc%nt)*lo_velocity_au_to_Afs
        call lo_h5_store_data(dr, h5%file_id, 'velocities', enhet='A/fs')
        lo_deallocate(dr)
        if (verbosity .gt. 0) write (*, *) '... wrote velocities'

        ! Write some alloy things?
        if (cc%alloy) then
            call lo_h5_store_data(cc%atomic_numbers, h5%file_id, 'atomic_numbers', enhet='Z')
            ! Then the alloy specification
            call lo_h5_store_data(cc%extra%unitcell_componentcounter, h5%file_id, 'unitcell_componentcounter')
            call lo_h5_store_data(cc%extra%supercell_componentcounter, h5%file_id, 'supercell_componentcounter')
            call lo_h5_store_data(cc%extra%unitcell_components, h5%file_id, 'unitcell_components')
            call lo_h5_store_data(cc%extra%supercell_components, h5%file_id, 'supercell_components')
            call lo_h5_store_data(cc%extra%unitcell_concentrations, h5%file_id, 'unitcell_concentrations')
            call lo_h5_store_data(cc%extra%supercell_concentrations, h5%file_id, 'supercell_concentrations')
        end if

        ! Write a lot of energies
        lo_allocate(ds(cc%nt))
        ds = cc%stat%polar_potential_energy(1:cc%nt)*lo_Hartree_to_eV
        call lo_h5_store_data(ds, h5%file_id, 'polar_potential_energy', enhet='eV')
        ds = cc%stat%secondorder_potential_energy(1:cc%nt)*lo_Hartree_to_eV
        call lo_h5_store_data(ds, h5%file_id, 'secondorder_potential_energy', enhet='eV')
        ds = cc%stat%thirdorder_potential_energy(1:cc%nt)*lo_Hartree_to_eV
        call lo_h5_store_data(ds, h5%file_id, 'thirdorder_potential_energy', enhet='eV')
        ds = cc%stat%fourthorder_potential_energy(1:cc%nt)*lo_Hartree_to_eV
        call lo_h5_store_data(ds, h5%file_id, 'fourthorder_potential_energy', enhet='eV')
        ds = cc%stat%kinetic_energy(1:cc%nt)*lo_Hartree_to_eV
        call lo_h5_store_data(ds, h5%file_id, 'kinetic_energy', enhet='eV')
        lo_deallocate(ds)

        ! Maybe some auxiliary stuff
        call lo_h5_store_data(uc%latticevectors*lo_bohr_to_A, h5%file_id, 'unitcell_latticevectors', enhet='A')
        call lo_h5_store_data(ss%latticevectors*lo_bohr_to_A, h5%file_id, 'supercell_latticevectors', enhet='A')
        call lo_h5_store_data(uc%r, h5%file_id, 'unitcell_positions', enhet='dimensionless')
        call lo_h5_store_data(ss%r, h5%file_id, 'supercell_positions', enhet='dimensionless')
        call lo_h5_store_data(uc%atomic_number, h5%file_id, 'unitcell_atomic_numbers', enhet='e')
        call lo_h5_store_data(ss%atomic_number, h5%file_id, 'supercell_atomic_numbers', enhet='e')
        if (verbosity .gt. 0) write (*, *) '... wrote energies and metadata'
        ! ! And close
        call h5%close_file()
        call h5%destroy(__FILE__, __LINE__)
    end block writefile

    if (verbosity .gt. 0) write (*, *) 'wrote simulation (', tochar(walltime() - timer), 's)'
end subroutine

!> write header to hdf5 (assumes single rank writing)
subroutine write_hdf5_header(cc, uc, ss, filename, total_configs, verbosity)
    !> md simulation
    class(lo_canonical_configs), intent(in) :: cc
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(inout) :: ss
    !> filename
    character(len=*), intent(in) :: filename
    !> total number of configurations 
    integer, intent(in) :: total_configs
    !> Talk a lot?
    integer, intent(in) :: verbosity

    type(lo_hdf5_helper) :: h5

    ! Initialize hdf5 properly
    call h5%init(__FILE__, __LINE__)
    call h5%open_file('write', trim(filename))

    ! Store some metadata. Not sure if effective.
    call lo_h5_store_attribute(cc%na, h5%file_id, 'number_of_atoms')
    call lo_h5_store_attribute(total_configs, h5%file_id, 'number_of_configurations')
    call lo_h5_store_attribute(cc%temperature_thermostat, h5%file_id, 'temperature_thermostat')
    call lo_h5_store_attribute(cc%alloy, h5%file_id, 'is_alloy')

    ! Write some alloy things?
    if (cc%alloy) then
        call lo_h5_store_data(cc%atomic_numbers, h5%file_id, 'atomic_numbers', enhet='Z')
        ! Then the alloy specification
        call lo_h5_store_data(cc%extra%unitcell_componentcounter, h5%file_id, 'unitcell_componentcounter')
        call lo_h5_store_data(cc%extra%supercell_componentcounter, h5%file_id, 'supercell_componentcounter')
        call lo_h5_store_data(cc%extra%unitcell_components, h5%file_id, 'unitcell_components')
        call lo_h5_store_data(cc%extra%supercell_components, h5%file_id, 'supercell_components')
        call lo_h5_store_data(cc%extra%unitcell_concentrations, h5%file_id, 'unitcell_concentrations')
        call lo_h5_store_data(cc%extra%supercell_concentrations, h5%file_id, 'supercell_concentrations')
    end if

    ! Maybe some auxiliary stuff
    call lo_h5_store_data(uc%latticevectors*lo_bohr_to_A, h5%file_id, 'unitcell_latticevectors', enhet='A')
    call lo_h5_store_data(ss%latticevectors*lo_bohr_to_A, h5%file_id, 'supercell_latticevectors', enhet='A')
    call lo_h5_store_data(uc%r, h5%file_id, 'unitcell_positions', enhet='dimensionless')
    call lo_h5_store_data(ss%r, h5%file_id, 'supercell_positions', enhet='dimensionless')
    call lo_h5_store_data(uc%atomic_number, h5%file_id, 'unitcell_atomic_numbers', enhet='e')
    call lo_h5_store_data(ss%atomic_number, h5%file_id, 'supercell_atomic_numbers', enhet='e')
    if (verbosity .gt. 0) write (*, *) '... wrote energies and metadata'       
    
    call h5%close_file()
    call h5%destroy(__FILE__, __LINE__)


end subroutine write_hdf5_header

!> Parallel write of:
!>   /header/... (small arrays; written by rank 1 only)
!>   /data/positions   : real(r8) [3, NA, NT_global]
!>   /data/velocities  : real(r8) [3, NA, NT_global]
!>   /data/energies/...  : real(r8) [NT_global] (4 vectors)
!>
!> Inputs:
!>   mw_comm  : MPI communicator (e.g., MPI_COMM_WORLD)   [integer]
!>   cc      : type(lo_canonical_configs) on each rank, holding this rank's local configs
!>   filename : character(*) HDF5 output path
!>
!> Layout/assumptions:
!>   cc%r, cc%v :: shape (3, NA, NT_local_on_this_rank)
!>   cc%stat%*   :: length NT_local_on_this_rank (same per-rank NT used for energies)
!>
subroutine write_hdf5_mpi(cc, mw, filename)

  class(lo_canonical_configs), intent(in) :: cc
  type(lo_mpi_helper), intent(inout) :: mw
  character(len=*), intent(in) :: filename


  ! Sizes (local/global) and offsets
  integer :: na, nt_local, nt_global, offset_ccs
  integer(HSIZE_T) :: dims_g3(3), dims_l3(3), start3(3), count3(3)
  integer(HSIZE_T) :: dims_g1(1), dims_l1(1), start1(1), count1(1)

  ! HDF5 handles
  integer(HID_T) :: fapl, file_id
  integer(HID_T) :: grp_header, grp_data
  integer(HID_T) :: dset_pos, dset_vel
  integer(HID_T) :: dset_ke, dset_pe_dd, dset_pe_h2, dset_pe_h3, dset_pe_h4
  integer(HID_T) :: filespace3, memspace3, filespace1, memspace1
  integer(HID_T) :: dxpl

  integer :: h5err, p 

  ! Convenience locals to header fields we’ll store (small demo subset)
  integer(HID_T) :: dset_ucell_lv, dset_scell_lv, dset_atnums
  integer(HSIZE_T) :: dims_2x(2), dims_1x(1)


  ! ===== Local sizes from the object on this rank =====
  na       = cc%na
  nt_local = cc%nt

  ! Global NT and starting offset along config dimension
  call MPI_Allreduce(nt_local, nt_global, 1, MPI_INTEGER, MPI_SUM, mw%comm, mw%error)
  call MPI_Exscan(nt_local, offset_ccs, 1, MPI_INTEGER, MPI_SUM, mw%comm, mw%error)
  if (mw%r == 0) offset_ccs = 0

  do p = 0, mw%n-1
    if (mw%r == p) then
        write(*,'(A,I0,3(A,I0))') 'rank ', mw%r, &
            '  na=', na, '  nt_local=', nt_local, '  offset=', offset_ccs
        flush(6)
    end if
    call mw%barrier()
  end do

  ! ===== HDF5 setup =====
  call h5open_f(h5err)

  call h5pcreate_f(H5P_FILE_ACCESS_F, fapl, h5err)
  call h5pset_fapl_mpio_f(fapl, mw%comm, MPI_INFO_NULL, h5err)

  call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, h5err, access_prp=fapl)

  ! Create groups collectively
  call h5gcreate_f(file_id, "header", grp_header, h5err)
  call h5gcreate_f(file_id, "data",   grp_data,   h5err)

  ! ===== Header datasets (created collectively, written by HEADER_RANK only) =====
  ! Example: unitcell/supercell lattice vectors and atomic numbers
  if (allocated(cc%extra%unitcell_latticevectors)) then
     dims_2x = (/ int(size(cc%extra%unitcell_latticevectors,1),HSIZE_T), &
                  int(size(cc%extra%unitcell_latticevectors,2),HSIZE_T) /)
     call h5screate_simple_f(2, dims_2x, filespace3, h5err)
     call h5dcreate_f(grp_header, "unitcell_latticevectors", H5T_NATIVE_DOUBLE, filespace3, dset_ucell_lv, h5err)
     call h5sclose_f(filespace3, h5err)
  end if

  if (allocated(cc%extra%supercell_latticevectors)) then
     dims_2x = (/ int(size(cc%extra%supercell_latticevectors,1),HSIZE_T), &
                  int(size(cc%extra%supercell_latticevectors,2),HSIZE_T) /)
     call h5screate_simple_f(2, dims_2x, filespace3, h5err)
     call h5dcreate_f(grp_header, "supercell_latticevectors", H5T_NATIVE_DOUBLE, filespace3, dset_scell_lv, h5err)
     call h5sclose_f(filespace3, h5err)
  end if

  if (allocated(cc%atomic_numbers)) then
     dims_1x = (/ int(size(cc%atomic_numbers,1),HSIZE_T) /)
     call h5screate_simple_f(1, dims_1x, filespace1, h5err)
     call h5dcreate_f(grp_header, "atomic_numbers", H5T_NATIVE_INTEGER, filespace1, dset_atnums, h5err)
     call h5sclose_f(filespace1, h5err)
  end if

  ! Only HEADER_RANK writes the small header payloads
  if (mw%talk) then
     if (allocated(cc%extra%unitcell_latticevectors)) then
        call h5dwrite_f(dset_ucell_lv, H5T_NATIVE_DOUBLE, cc%extra%unitcell_latticevectors, &
                        shape(cc%extra%unitcell_latticevectors, kind=HSIZE_T), h5err)
     end if
     if (allocated(cc%extra%supercell_latticevectors)) then
        call h5dwrite_f(dset_scell_lv, H5T_NATIVE_DOUBLE, cc%extra%supercell_latticevectors, &
                        shape(cc%extra%supercell_latticevectors, kind=HSIZE_T), h5err)
     end if
     if (allocated(cc%atomic_numbers)) then
        call h5dwrite_f(dset_atnums, H5T_NATIVE_INTEGER, cc%atomic_numbers, &
                        (/ int(size(cc%atomic_numbers),HSIZE_T) /), h5err)
     end if
  end if

  ! Close header datasets if they were created
  if (allocated(cc%extra%unitcell_latticevectors)) call h5dclose_f(dset_ucell_lv, h5err)
  if (allocated(cc%extra%supercell_latticevectors)) call h5dclose_f(dset_scell_lv, h5err)
  if (allocated(cc%atomic_numbers))                call h5dclose_f(dset_atnums,   h5err)

  call mw%barrier() ! ensure header is written before large dataset I/O

  ! ===== Create global 3D datasets for positions/velocities (collective) =====
  dims_g3 = (/ 3_hsize_t, int(na,HSIZE_T), int(nt_global,HSIZE_T) /)

  call h5screate_simple_f(3, dims_g3, filespace3, h5err)
  call h5dcreate_f(grp_data, "positions",  H5T_NATIVE_DOUBLE, filespace3, dset_pos, h5err)
  call h5dclose_f(dset_pos, h5err)   ! close/reopen not needed, but free filespace3 reuse clarity
  call h5sclose_f(filespace3, h5err)

  call h5screate_simple_f(3, dims_g3, filespace3, h5err)
  call h5dcreate_f(grp_data, "velocities", H5T_NATIVE_DOUBLE, filespace3, dset_vel, h5err)
  call h5dclose_f(dset_vel, h5err)
  call h5sclose_f(filespace3, h5err)

  ! Energies (4 vectors) — create now
  dims_g1 = (/ int(nt_global,HSIZE_T) /)
  call h5screate_simple_f(1, dims_g1, filespace1, h5err)
  call h5dcreate_f(grp_data, "kinetic_energy",            H5T_NATIVE_DOUBLE, filespace1, dset_ke,     h5err)
  call h5dclose_f(dset_ke, h5err)
  call h5dcreate_f(grp_data, "polar_potential_energy",    H5T_NATIVE_DOUBLE, filespace1, dset_pe_dd,  h5err)
  call h5dclose_f(dset_pe_dd, h5err)
  call h5dcreate_f(grp_data, "secondorder_potential_energy", H5T_NATIVE_DOUBLE, filespace1, dset_pe_h2, h5err)
  call h5dclose_f(dset_pe_h2, h5err)
  call h5dcreate_f(grp_data, "thirdorder_potential_energy",  H5T_NATIVE_DOUBLE, filespace1, dset_pe_h3, h5err)
  call h5dclose_f(dset_pe_h3, h5err)
  call h5dcreate_f(grp_data, "fourthorder_potential_energy", H5T_NATIVE_DOUBLE, filespace1, dset_pe_h4, h5err)
  call h5sclose_f(filespace1, h5err)
  call h5dclose_f(dset_pe_h4, h5err)

  ! ===== Parallel writes (collective) for positions & velocities =====
  ! Re-open datasets and set up collective xfer property
  call h5dopen_f(grp_data, "positions",  dset_pos, h5err)
  call h5dopen_f(grp_data, "velocities", dset_vel, h5err)

  call h5pcreate_f(H5P_DATASET_XFER_F, dxpl, h5err)
  call h5pset_dxpl_mpio_f(dxpl, H5FD_MPIO_COLLECTIVE_F, h5err)

  ! File hyperslab: (0,0,offset) with count (3, NA, nt_local)
  call h5dget_space_f(dset_pos, filespace3, h5err)
  start3 = (/ 0_hsize_t, 0_hsize_t, int(offset_ccs,HSIZE_T) /)
  count3 = (/ 3_hsize_t, int(na,HSIZE_T), int(nt_local,HSIZE_T) /)
  call h5sselect_hyperslab_f(filespace3, H5S_SELECT_SET_F, start3, count3, h5err)

  ! Memory dataspace: (3, NA, nt_local)
  dims_l3 = count3
  call h5screate_simple_f(3, dims_l3, memspace3, h5err)

  ! Parallel write positions
  if (nt_local > 0) then
    call h5dwrite_f(dset_pos, H5T_NATIVE_DOUBLE, cc%r, dims_l3, h5err, &
         file_space_id=filespace3, mem_space_id=memspace3, xfer_prp=dxpl)
  else
    ! zero-length selection is okay; still collective call
    call h5dwrite_f(dset_pos, H5T_NATIVE_DOUBLE, cc%r, dims_l3, h5err, &
         file_space_id=filespace3, mem_space_id=memspace3, xfer_prp=dxpl)
  end if

  ! Repeat for velocities
  call h5sclose_f(filespace3, h5err)
  call h5dget_space_f(dset_vel, filespace3, h5err)
  start3 = (/ 0_hsize_t, 0_hsize_t, int(offset_ccs,HSIZE_T) /)
  count3 = (/ 3_hsize_t, int(na,HSIZE_T), int(nt_local,HSIZE_T) /)
  call h5sselect_hyperslab_f(filespace3, H5S_SELECT_SET_F, start3, count3, h5err)

  ! Reuse memspace3
  if (nt_local > 0) then
    call h5dwrite_f(dset_vel, H5T_NATIVE_DOUBLE, cc%v, dims_l3, h5err, &
         file_space_id=filespace3, mem_space_id=memspace3, xfer_prp=dxpl)
  else
    call h5dwrite_f(dset_vel, H5T_NATIVE_DOUBLE, cc%v, dims_l3, h5err, &
         file_space_id=filespace3, mem_space_id=memspace3, xfer_prp=dxpl)
  end if

  call h5sclose_f(filespace3, h5err)
  call h5sclose_f(memspace3,   h5err)
  call h5dclose_f(dset_pos,    h5err)
  call h5dclose_f(dset_vel,    h5err)

  ! ===== Parallel writes for the four energy vectors =====
  ! Reopen each vector dataset and write my [offset:offset+nt_local)
  call h5dopen_f(grp_data, "kinetic_energy",                 dset_ke,    h5err)
  call h5dopen_f(grp_data, "polar_potential_energy",         dset_pe_dd, h5err)
  call h5dopen_f(grp_data, "secondorder_potential_energy",   dset_pe_h2, h5err)
  call h5dopen_f(grp_data, "thirdorder_potential_energy",    dset_pe_h3, h5err)
  call h5dopen_f(grp_data, "fourthorder_potential_energy",   dset_pe_h4, h5err)

  dims_g1 = (/ int(nt_global,HSIZE_T) /)
  dims_l1 = (/ int(nt_local,HSIZE_T) /)
  start1  = (/ int(offset_ccs,HSIZE_T) /)
  count1  = dims_l1

  ! KE
  call h5dget_space_f(dset_ke, filespace1, h5err)
  call h5sselect_hyperslab_f(filespace1, H5S_SELECT_SET_F, start1, count1, h5err)
  call h5screate_simple_f(1, dims_l1, memspace1, h5err)
  call h5dwrite_f(dset_ke, H5T_NATIVE_DOUBLE, cc%stat%kinetic_energy, dims_l1, h5err, &
       file_space_id=filespace1, mem_space_id=memspace1, xfer_prp=dxpl)
  call h5sclose_f(filespace1, h5err)
  call h5sclose_f(memspace1, h5err)

  ! PE (dipole-dipole)
  call h5dget_space_f(dset_pe_dd, filespace1, h5err)
  call h5sselect_hyperslab_f(filespace1, H5S_SELECT_SET_F, start1, count1, h5err)
  call h5screate_simple_f(1, dims_l1, memspace1, h5err)
  call h5dwrite_f(dset_pe_dd, H5T_NATIVE_DOUBLE, cc%stat%polar_potential_energy, dims_l1, h5err, &
       file_space_id=filespace1, mem_space_id=memspace1, xfer_prp=dxpl)
  call h5sclose_f(filespace1, h5err)
  call h5sclose_f(memspace1, h5err)

  ! Harmonic 2nd order
  call h5dget_space_f(dset_pe_h2, filespace1, h5err)
  call h5sselect_hyperslab_f(filespace1, H5S_SELECT_SET_F, start1, count1, h5err)
  call h5screate_simple_f(1, dims_l1, memspace1, h5err)
  call h5dwrite_f(dset_pe_h2, H5T_NATIVE_DOUBLE, cc%stat%secondorder_potential_energy, dims_l1, h5err, &
       file_space_id=filespace1, mem_space_id=memspace1, xfer_prp=dxpl)
  call h5sclose_f(filespace1, h5err)
  call h5sclose_f(memspace1, h5err)

  ! 3rd order
  call h5dget_space_f(dset_pe_h3, filespace1, h5err)
  call h5sselect_hyperslab_f(filespace1, H5S_SELECT_SET_F, start1, count1, h5err)
  call h5screate_simple_f(1, dims_l1, memspace1, h5err)
  call h5dwrite_f(dset_pe_h3, H5T_NATIVE_DOUBLE, cc%stat%thirdorder_potential_energy, dims_l1, h5err, &
       file_space_id=filespace1, mem_space_id=memspace1, xfer_prp=dxpl)
  call h5sclose_f(filespace1, h5err)
  call h5sclose_f(memspace1, h5err)

  ! 4th order
  call h5dget_space_f(dset_pe_h4, filespace1, h5err)
  call h5sselect_hyperslab_f(filespace1, H5S_SELECT_SET_F, start1, count1, h5err)
  call h5screate_simple_f(1, dims_l1, memspace1, h5err)
  call h5dwrite_f(dset_pe_h4, H5T_NATIVE_DOUBLE, cc%stat%fourthorder_potential_energy, dims_l1, h5err, &
       file_space_id=filespace1, mem_space_id=memspace1, xfer_prp=dxpl)
  call h5sclose_f(filespace1, h5err)
  call h5sclose_f(memspace1, h5err)

  call h5dclose_f(dset_ke,    h5err)
  call h5dclose_f(dset_pe_dd, h5err)
  call h5dclose_f(dset_pe_h2, h5err)
  call h5dclose_f(dset_pe_h3, h5err)
  call h5dclose_f(dset_pe_h4, h5err)

  ! ===== Tear down =====
  call h5pclose_f(dxpl, h5err)
  call h5gclose_f(grp_data,   h5err)
  call h5gclose_f(grp_header, h5err)
  call h5fclose_f(file_id,    h5err)
  call h5pclose_f(fapl,       h5err)
  call h5close_f(h5err)
end subroutine write_hdf5_mpi

! subroutine write_hdf5_mpi(cc, mw, filename)
!   class(lo_canonical_configs), intent(in) :: cc
!   type(lo_mpi_helper), intent(inout) :: mw
!   character(len=*), intent(in) :: filename

!   ! Sizes (local/global) and offsets
!   integer :: na, nt_local, nt_global, offset_ccs
!   integer(HSIZE_T) :: dims_g3(3), dims_l3(3), start3(3), count3(3)
!   integer(HSIZE_T) :: dims_g1(1), dims_l1(1), start1(1), count1(1)

!   ! HDF5 handles
!   integer(HID_T) :: fapl, file_id
!   integer(HID_T) :: grp_header, grp_data
!   integer(HID_T) :: dset_pos, dset_vel
!   integer(HID_T) :: dset_ke, dset_pe_dd, dset_pe_h2, dset_pe_h3, dset_pe_h4
!   integer(HID_T) :: filespace3, memspace3, filespace1, memspace1
!   integer(HID_T) :: dxpl

!   ! header helpers
!   integer(HID_T) :: dset_ucell_lv, dset_scell_lv, dset_atnums
!   integer(HSIZE_T) :: dims_2x(2), dims_1x(1)

!   ! DCPL / compression
!   integer(HID_T)   :: dcpl3, dcpl1
!   integer(HSIZE_T) :: chunk3(3), chunk1(1)
!   logical          :: have_deflate
!   integer          :: h5err, h5err2, p

!   ! ===== Local sizes on this rank =====
!   na       = cc%na
!   nt_local = cc%nt

!   ! Global NT and offset along config dimension
!   call MPI_Allreduce(nt_local, nt_global, 1, MPI_INTEGER, MPI_SUM, mw%comm, mw%error)
!   call MPI_Exscan(nt_local, offset_ccs, 1, MPI_INTEGER, MPI_SUM, mw%comm, mw%error)
!   if (mw%r == 0) offset_ccs = 0

!   do p = 0, mw%n-1
!     if (mw%r == p) then
!         write(*,'(A,I0,3(A,I0))') 'rank ', mw%r, &
!             '  na=', na, '  nt_local=', nt_local, '  offset=', offset_ccs
!         flush(6)
!     end if
!     call mw%barrier()
!   end do

!   ! ===== HDF5 setup =====
!   call h5open_f(h5err)
!   call h5pcreate_f(H5P_FILE_ACCESS_F, fapl, h5err)
!   call h5pset_fapl_mpio_f(fapl, mw%comm, MPI_INFO_NULL, h5err)

!   ! Helpful for parallel metadata (if available in your HDF5 build)
!   call h5pset_all_coll_metadata_ops_f(fapl, .true., h5err)
!   call h5pset_coll_metadata_write_f(fapl, .true., h5err)

!   call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, h5err, access_prp=fapl)

!   ! Create groups collectively
!   call h5gcreate_f(file_id, "header", grp_header, h5err)
!   call h5gcreate_f(file_id, "data",   grp_data,   h5err)

!   ! ===== Header datasets (created collectively, payload written by talker only) =====
!   if (allocated(cc%extra%unitcell_latticevectors)) then
!      dims_2x = (/ int(size(cc%extra%unitcell_latticevectors,1),HSIZE_T), &
!                   int(size(cc%extra%unitcell_latticevectors,2),HSIZE_T) /)
!      call h5screate_simple_f(2, dims_2x, filespace3, h5err)
!      call h5dcreate_f(grp_header, "unitcell_latticevectors", H5T_NATIVE_DOUBLE, filespace3, dset_ucell_lv, h5err)
!      call h5sclose_f(filespace3, h5err)
!   end if

!   if (allocated(cc%extra%supercell_latticevectors)) then
!      dims_2x = (/ int(size(cc%extra%supercell_latticevectors,1),HSIZE_T), &
!                   int(size(cc%extra%supercell_latticevectors,2),HSIZE_T) /)
!      call h5screate_simple_f(2, dims_2x, filespace3, h5err)
!      call h5dcreate_f(grp_header, "supercell_latticevectors", H5T_NATIVE_DOUBLE, filespace3, dset_scell_lv, h5err)
!      call h5sclose_f(filespace3, h5err)
!   end if

!   if (allocated(cc%atomic_numbers)) then
!      dims_1x = (/ int(size(cc%atomic_numbers,1),HSIZE_T) /)
!      call h5screate_simple_f(1, dims_1x, filespace1, h5err)
!      call h5dcreate_f(grp_header, "atomic_numbers", H5T_NATIVE_INTEGER, filespace1, dset_atnums, h5err)
!      call h5sclose_f(filespace1, h5err)
!   end if

!   if (mw%talk) then
!      if (allocated(cc%extra%unitcell_latticevectors)) then
!         call h5dwrite_f(dset_ucell_lv, H5T_NATIVE_DOUBLE, cc%extra%unitcell_latticevectors, &
!                         shape(cc%extra%unitcell_latticevectors, kind=HSIZE_T), h5err)
!      end if
!      if (allocated(cc%extra%supercell_latticevectors)) then
!         call h5dwrite_f(dset_scell_lv, H5T_NATIVE_DOUBLE, cc%extra%supercell_latticevectors, &
!                         shape(cc%extra%supercell_latticevectors, kind=HSIZE_T), h5err)
!      end if
!      if (allocated(cc%atomic_numbers)) then
!         call h5dwrite_f(dset_atnums, H5T_NATIVE_INTEGER, cc%atomic_numbers, &
!                         (/ int(size(cc%atomic_numbers),HSIZE_T) /), h5err)
!      end if
!   end if

!   if (allocated(cc%extra%unitcell_latticevectors)) call h5dclose_f(dset_ucell_lv, h5err)
!   if (allocated(cc%extra%supercell_latticevectors)) call h5dclose_f(dset_scell_lv, h5err)
!   if (allocated(cc%atomic_numbers))                call h5dclose_f(dset_atnums,   h5err)

!   call mw%barrier()

!   ! ================== DATASET CREATION with robust DCPL ==================
!   dims_g3 = (/ int(3,HSIZE_T), int(na,HSIZE_T), int(nt_global,HSIZE_T) /)
!   dims_g1 = (/ int(nt_global,HSIZE_T) /)

!   ! decide chunking (align chunks with your access: full 3xNA plane per time-slice)
!   if (nt_global > 0) then
!      chunk3 = (/ int(3,HSIZE_T), int(na,HSIZE_T), int(min(max(1,64), nt_global),HSIZE_T) /)
!      chunk1 = (/ int(min(max(1,4096), nt_global),HSIZE_T) /)
!   else
!      chunk3 = (/ int(3,HSIZE_T), int(na,HSIZE_T), int(1,HSIZE_T) /)
!      chunk1 = (/ int(1,HSIZE_T) /)
!   end if

!   call h5zfilter_avail_f(H5Z_FILTER_DEFLATE_F, have_deflate, h5err)

!   ! ---- positions ----
!   call h5screate_simple_f(3, dims_g3, filespace3, h5err)
!   call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl3, h5err)
!   call h5pset_chunk_f(dcpl3, 3, chunk3, h5err)
!   if (have_deflate) call h5pset_deflate_f(dcpl3, 6, h5err)
!   call h5dcreate_f(grp_data, "positions", H5T_NATIVE_DOUBLE, filespace3, dset_pos, h5err, dcpl_id=dcpl3)
!   if (h5err < 0) then
!      if (mw%r==0) write(*,*) "Create with deflate failed for positions; retrying without filters"
!      call h5pclose_f(dcpl3, h5err2)
!      call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl3, h5err)
!      call h5pset_chunk_f(dcpl3, 3, chunk3, h5err)
!      call h5dcreate_f(grp_data, "positions", H5T_NATIVE_DOUBLE, filespace3, dset_pos, h5err, dcpl_id=dcpl3)
!      if (h5err < 0) then
!         if (mw%r==0) write(*,*) "Chunked create failed for positions; falling back to contiguous"
!         call h5dcreate_f(grp_data, "positions", H5T_NATIVE_DOUBLE, filespace3, dset_pos, h5err)
!      end if
!   end if
!   call h5pclose_f(dcpl3, h5err2)
!   call h5sclose_f(filespace3, h5err)

!   ! ---- velocities ----
!   call h5screate_simple_f(3, dims_g3, filespace3, h5err)
!   call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl3, h5err)
!   call h5pset_chunk_f(dcpl3, 3, chunk3, h5err)
!   if (have_deflate) call h5pset_deflate_f(dcpl3, 6, h5err)
!   call h5dcreate_f(grp_data, "velocities", H5T_NATIVE_DOUBLE, filespace3, dset_vel, h5err, dcpl_id=dcpl3)
!   if (h5err < 0) then
!      if (mw%r==0) write(*,*) "Create with deflate failed for velocities; retrying without filters"
!      call h5pclose_f(dcpl3, h5err2)
!      call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl3, h5err)
!      call h5pset_chunk_f(dcpl3, 3, chunk3, h5err)
!      call h5dcreate_f(grp_data, "velocities", H5T_NATIVE_DOUBLE, filespace3, dset_vel, h5err, dcpl_id=dcpl3)
!      if (h5err < 0) then
!         if (mw%r==0) write(*,*) "Chunked create failed for velocities; falling back to contiguous"
!         call h5dcreate_f(grp_data, "velocities", H5T_NATIVE_DOUBLE, filespace3, dset_vel, h5err)
!      end if
!   end if
!   call h5pclose_f(dcpl3, h5err2)
!   call h5sclose_f(filespace3, h5err)

!   ! ---- energy vectors ----
!   call h5screate_simple_f(1, dims_g1, filespace1, h5err)
!   call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl1, h5err)
!   call h5pset_chunk_f(dcpl1, 1, chunk1, h5err)
!   if (have_deflate) call h5pset_deflate_f(dcpl1, 6, h5err)

!   call h5dcreate_f(grp_data, "kinetic_energy", H5T_NATIVE_DOUBLE, filespace1, dset_ke,    h5err, dcpl_id=dcpl1)
!   if (h5err < 0) call h5dcreate_f(grp_data, "kinetic_energy", H5T_NATIVE_DOUBLE, filespace1, dset_ke,    h5err)

!   call h5dcreate_f(grp_data, "polar_potential_energy", H5T_NATIVE_DOUBLE, filespace1, dset_pe_dd, h5err, dcpl_id=dcpl1)
!   if (h5err < 0) call h5dcreate_f(grp_data, "polar_potential_energy", H5T_NATIVE_DOUBLE, filespace1, dset_pe_dd, h5err)

!   call h5dcreate_f(grp_data, "secondorder_potential_energy", H5T_NATIVE_DOUBLE, filespace1, dset_pe_h2, h5err, dcpl_id=dcpl1)
!   if (h5err < 0) call h5dcreate_f(grp_data, "secondorder_potential_energy", H5T_NATIVE_DOUBLE, filespace1, dset_pe_h2, h5err)

!   call h5dcreate_f(grp_data, "thirdorder_potential_energy", H5T_NATIVE_DOUBLE, filespace1, dset_pe_h3, h5err, dcpl_id=dcpl1)
!   if (h5err < 0) call h5dcreate_f(grp_data, "thirdorder_potential_energy", H5T_NATIVE_DOUBLE, filespace1, dset_pe_h3, h5err)

!   call h5dcreate_f(grp_data, "fourthorder_potential_energy", H5T_NATIVE_DOUBLE, filespace1, dset_pe_h4, h5err, dcpl_id=dcpl1)
!   if (h5err < 0) call h5dcreate_f(grp_data, "fourthorder_potential_energy", H5T_NATIVE_DOUBLE, filespace1, dset_pe_h4, h5err)

!   call h5pclose_f(dcpl1, h5err2)
!   call h5sclose_f(filespace1, h5err)

!   ! ================== COLLECTIVE WRITES ==================
!   call h5dopen_f(grp_data, "positions",  dset_pos, h5err)
!   call h5dopen_f(grp_data, "velocities", dset_vel, h5err)

!   call h5pcreate_f(H5P_DATASET_XFER_F, dxpl, h5err)
!   call h5pset_dxpl_mpio_f(dxpl, H5FD_MPIO_COLLECTIVE_F, h5err)

!   ! File hyperslab: (0,0,offset) with count (3, NA, nt_local)
!   call h5dget_space_f(dset_pos, filespace3, h5err)
!   start3 = (/ int(0,HSIZE_T), int(0,HSIZE_T), int(offset_ccs,HSIZE_T) /)
!   count3 = (/ int(3,HSIZE_T), int(na,HSIZE_T), int(nt_local,HSIZE_T) /)
!   call h5sselect_hyperslab_f(filespace3, H5S_SELECT_SET_F, start3, count3, h5err)

!   dims_l3 = count3
!   call h5screate_simple_f(3, dims_l3, memspace3, h5err)
!   call h5dwrite_f(dset_pos, H5T_NATIVE_DOUBLE, cc%r, dims_l3, h5err, &
!        file_space_id=filespace3, mem_space_id=memspace3, xfer_prp=dxpl)

!   call h5sclose_f(filespace3, h5err)
!   call h5dget_space_f(dset_vel, filespace3, h5err)
!   call h5sselect_hyperslab_f(filespace3, H5S_SELECT_SET_F, start3, count3, h5err)
!   call h5dwrite_f(dset_vel, H5T_NATIVE_DOUBLE, cc%v, dims_l3, h5err, &
!        file_space_id=filespace3, mem_space_id=memspace3, xfer_prp=dxpl)

!   call h5sclose_f(filespace3, h5err)
!   call h5sclose_f(memspace3,   h5err)
!   call h5dclose_f(dset_pos,    h5err)
!   call h5dclose_f(dset_vel,    h5err)

!   ! Energies
!   call h5dopen_f(grp_data, "kinetic_energy",               dset_ke,    h5err)
!   call h5dopen_f(grp_data, "polar_potential_energy",       dset_pe_dd, h5err)
!   call h5dopen_f(grp_data, "secondorder_potential_energy", dset_pe_h2, h5err)
!   call h5dopen_f(grp_data, "thirdorder_potential_energy",  dset_pe_h3, h5err)
!   call h5dopen_f(grp_data, "fourthorder_potential_energy", dset_pe_h4, h5err)

!   dims_l1 = (/ int(nt_local,HSIZE_T) /)
!   start1  = (/ int(offset_ccs,HSIZE_T) /)
!   count1  = dims_l1

!   call h5dget_space_f(dset_ke, filespace1, h5err)
!   call h5sselect_hyperslab_f(filespace1, H5S_SELECT_SET_F, start1, count1, h5err)
!   call h5screate_simple_f(1, dims_l1, memspace1, h5err)
!   call h5dwrite_f(dset_ke, H5T_NATIVE_DOUBLE, cc%stat%kinetic_energy, dims_l1, h5err, &
!        file_space_id=filespace1, mem_space_id=memspace1, xfer_prp=dxpl)
!   call h5sclose_f(filespace1, h5err)
!   call h5sclose_f(memspace1, h5err)

!   call h5dget_space_f(dset_pe_dd, filespace1, h5err)
!   call h5sselect_hyperslab_f(filespace1, H5S_SELECT_SET_F, start1, count1, h5err)
!   call h5screate_simple_f(1, dims_l1, memspace1, h5err)
!   call h5dwrite_f(dset_pe_dd, H5T_NATIVE_DOUBLE, cc%stat%polar_potential_energy, dims_l1, h5err, &
!        file_space_id=filespace1, mem_space_id=memspace1, xfer_prp=dxpl)
!   call h5sclose_f(filespace1, h5err)
!   call h5sclose_f(memspace1, h5err)

!   call h5dget_space_f(dset_pe_h2, filespace1, h5err)
!   call h5sselect_hyperslab_f(filespace1, H5S_SELECT_SET_F, start1, count1, h5err)
!   call h5screate_simple_f(1, dims_l1, memspace1, h5err)
!   call h5dwrite_f(dset_pe_h2, H5T_NATIVE_DOUBLE, cc%stat%secondorder_potential_energy, dims_l1, h5err, &
!        file_space_id=filespace1, mem_space_id=memspace1, xfer_prp=dxpl)
!   call h5sclose_f(filespace1, h5err)
!   call h5sclose_f(memspace1, h5err)

!   call h5dget_space_f(dset_pe_h3, filespace1, h5err)
!   call h5sselect_hyperslab_f(filespace1, H5S_SELECT_SET_F, start1, count1, h5err)
!   call h5screate_simple_f(1, dims_l1, memspace1, h5err)
!   call h5dwrite_f(dset_pe_h3, H5T_NATIVE_DOUBLE, cc%stat%thirdorder_potential_energy, dims_l1, h5err, &
!        file_space_id=filespace1, mem_space_id=memspace1, xfer_prp=dxpl)
!   call h5sclose_f(filespace1, h5err)
!   call h5sclose_f(memspace1, h5err)

!   call h5dget_space_f(dset_pe_h4, filespace1, h5err)
!   call h5sselect_hyperslab_f(filespace1, H5S_SELECT_SET_F, start1, count1, h5err)
!   call h5screate_simple_f(1, dims_l1, memspace1, h5err)
!   call h5dwrite_f(dset_pe_h4, H5T_NATIVE_DOUBLE, cc%stat%fourthorder_potential_energy, dims_l1, h5err, &
!        file_space_id=filespace1, mem_space_id=memspace1, xfer_prp=dxpl)
!   call h5sclose_f(filespace1, h5err)
!   call h5sclose_f(memspace1, h5err)

!   call h5dclose_f(dset_ke,    h5err)
!   call h5dclose_f(dset_pe_dd, h5err)
!   call h5dclose_f(dset_pe_h2, h5err)
!   call h5dclose_f(dset_pe_h3, h5err)
!   call h5dclose_f(dset_pe_h4, h5err)

!   ! ===== Tear down =====
!   call h5pclose_f(dxpl, h5err)
!   call h5gclose_f(grp_data,   h5err)
!   call h5gclose_f(grp_header, h5err)
!   call h5fclose_f(file_id,    h5err)
!   call h5pclose_f(fapl,       h5err)
!   call h5close_f(h5err)
! end subroutine write_hdf5_mpi

end module
