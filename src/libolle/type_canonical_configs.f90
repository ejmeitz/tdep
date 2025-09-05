#include "precompilerdefinitions"
module type_canonical_configs
!! Information about an MD simulation
use konstanter, only: r8, lo_pi, lo_huge, lo_hugeint, lo_sqtol, lo_status, &
                      lo_exitcode_param, lo_bohr_to_A, lo_Hartree_to_eV, &
                      lo_exitcode_io, lo_velocity_au_to_Afs
use gottochblandat, only: open_file, tochar, walltime 
use mpi_wrappers, only:  lo_stop_gracefully
use hdf5_wrappers, only: lo_hdf5_helper, lo_h5_store_attribute, lo_h5_store_data
use type_crystalstructure, only: lo_crystalstructure

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

    ! Max number of timesteps
    tmax = size(cc%r, 3)

    ! Sanity tests
    if (idx .gt. tmax) then
        call lo_stop_gracefully(['Not enough space to store timestep. Initialize canonical_config storage with more space.'], lo_exitcode_param, __FILE__, __LINE__)
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
            !lo_allocate(cc%r_ref(3, cc%na, nstep))
            lo_allocate(cc%atomic_numbers(cc%na))
            !cc%r_ref = 0.0_r8
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

!> write a simulation to hdf5
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


end module
