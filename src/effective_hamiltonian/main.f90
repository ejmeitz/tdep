#include "precompilerdefinitions"
program effective_hamiltonian
!!{!src/effective_hamiltonian/manual.md!}
use konstanter, only: r8, lo_tol, lo_kb_hartree, lo_bohr_to_A, lo_frequency_Hartree_to_THz, lo_eV_to_Hartree, &
                      lo_exitcode_physical, lo_exitcode_param, lo_temperaturetol, lo_Hartree_to_eV
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
use gottochblandat, only: tochar, walltime, lo_stop_gracefully, open_file, lo_progressbar_init, &
                            lo_progressbar, lo_does_file_exist, lo_trueNtimes, lo_mean

use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_canonical_configs, only: lo_canonical_configs
use type_mdsim, only: lo_mdsim


use lo_epot, only: lo_energy_differences

use type_crystalstructure, only: lo_crystalstructure
use options, only: lo_opts

implicit none
type(lo_opts) :: opts
type(lo_crystalstructure) :: ss, uc
type(lo_energy_differences) :: pot
type(lo_mdsim) :: sim
type(lo_canonical_configs) :: cc


type(lo_forceconstant_secondorder) :: fc2
type(lo_forceconstant_thirdorder) :: fc3
type(lo_forceconstant_fourthorder) :: fc4
! POLAR IFCS?

type(lo_mpi_helper) :: mw
type(lo_mem_helper) :: mem
real(r8), dimension(:, :), allocatable :: ebuf

logical :: generate_configs = .false.


call opts%parse()
call mw%init()
call mem%init()

init: block

    integer :: f, i, j, l, readrank, local_nconf
    logical :: readonthisrank, mpiparallel
    real(r8) :: t0


    if (.not. mw%talk) opts%verbosity = -100

    if (opts%nconf .gt. 0) generate_configs = .true.

    if (mw%talk) then
        write (*, *) 'Recap of the parameters governing the calculation'
        write (*, '(1X,A40,L3)') 'Thirdorder contribution                 ', opts%thirdorder
        write (*, '(1X,A40,L3)') 'Fourthorder contribution                ', opts%fourthorder
        write (*, '(1X,A40,I5)') 'Stride                                  ', opts%stride
        write(*, '(1X,A40,L3)') 'Generate canonical configs               ', generate_configs
        if(generate_configs) then
            write(*, '(1X,A40,I8)') 'Num Configs                          ', opts%nconf
            write (*, '(1X,A40,L3)') 'Quantum configurations              ', opts%quantum
            write (*, '(1X,A40,F20.12)') 'Temperature                     ', opts%temperature
            write (*, '(1X,A40,L3)') 'Dump configurations                 ', opts%dumpconfigs
        end if
    end if

    ! Read structures
    if (mw%talk) write (*, *) '... reading infiles'
    call ss%readfromfile('infile.ssposcar', verbosity=opts%verbosity)
    call uc%readfromfile('infile.ucposcar', verbosity=opts%verbosity)

    ! Match the supercell to the unitcell, always a good idea
    call uc%classify('spacegroup', timereversal=.true.)
    call ss%classify('supercell', uc)

    call fc2%readfromfile(uc, 'infile.forceconstant', mem, verbosity=-1)
    if (opts%thirdorder) call fc3%readfromfile(uc, 'infile.forceconstant_thirdorder')
    if (opts%fourthorder) call fc4%readfromfile(uc, 'infile.forceconstant_fourthorder')
    if (mw%talk) write (*, *) '... read forceconstants'

    if (generate_configs) then
        if ((opts%quantum .eqv. .false.) .and. (opts%temperature .lt. lo_temperaturetol)) then
            call lo_stop_gracefully(['For classical statistics temperature has to be nonzero'], lo_exitcode_physical, __FILE__, __LINE__)
        end if
        if (mw%talk .and. opts%quantum) write (*, *) '... will generate quantum configurations'
        if (mw%talk .and. opts%quantum) write (*, *) '... will generate classical configurations'
    
        if (opts%temperature .lt. lo_temperaturetol) then
            call lo_stop_gracefully(['If nconf is passed, the temperature flag must be set'], lo_exitcode_param, __FILE__, __LINE__)
        end if

    end if

    call pot%setup(uc, ss, fc2, fc3, fc4, mw, opts%verbosity + 1)
    if (mw%talk) write (*, *) '... setup potential energy calculator'

    if (.not. generate_configs) then
        ! If packed simulation is there might as well read it
        if (lo_does_file_exist('infile.sim.hdf5')) then
            call sim%read_from_hdf5('infile.sim.hdf5', verbosity=opts%verbosity + 2, stride=opts%stride)
        else 
            call sim%read_from_file(verbosity = opts%verbosity + 2, stride = opts%stride, &
                                     magnetic=.false., dielectric=.false., nrand=-1, mw=mw)
        end if
        if (mw%talk) write (*, *) '... parsed simulation data'
    else
        if (opts%dumpconfigs) then
            ! Calculate local number of configurations for this rank
            ! Assumes round-robin parallelization
            if (mw%r == 0) then
                local_nconf = opts%nconf / mw%n
            else if (mw%r <= opts%nconf) then
                local_nconf = (opts%nconf - mw%r) / mw%n + 1
            else
                local_nconf = 0
            end if 

            call cc%init_empty(uc, ss, local_nconf, opts%temperature)

        end if
    end if

end block init




energy : block

    integer :: i, u, ierr
    real(r8), dimension(:, :), allocatable :: f2, f3, f4, fp
    real(r8) :: e2, e3, e4, ep, total_energy, to_ev_per_atom
    character(len=100) :: filename
    integer :: global_offset
    ! real(r8), dimension(:, :, :) allocatable :: r_buf, v_buf


    if (generate_configs) then

        if (mw%talk) write (*, *) '... generating canonical configurations'
        call mem%allocate(ebuf, [opts%nconf, 5], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        ebuf = 0.0_r8

        if (opts%dumpconfigs) then
            call pot%statistical_sampling(uc, ss, fc2, opts%nconf, opts%temperature, opts%quantum, ebuf, mw, mem, opts%verbosity, cc)
            call cc%write_hdf5_mpi(mw, 'outfile.canonical_configs.hdf5')
        else
            call pot%statistical_sampling(uc, ss, fc2, opts%nconf, opts%temperature, opts%quantum, ebuf, mw, mem, opts%verbosity)
        end if

    else

        call mem%allocate(ebuf, [sim%nt, 4], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        ebuf = 0.0_r8

        ! Dummy space for force
        call mem%allocate(f2, [3, ss%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(f3, [3, ss%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(f4, [3, ss%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(fp, [3, ss%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        f2 = 0.0_r8
        f3 = 0.0_r8
        f4 = 0.0_r8
        fp = 0.0_r8
        
        ! how do I add a progress bar here without race condition?
        do i = 1, sim%nt 

            if (mod(i, mw%n) .ne. mw%r) cycle

            ! Calculate the energy, e2/e3/e4/ep are zeroed inside of this call
            call pot%energies_and_forces(sim%u(:,:,i), e2, e3, e4, ep, f2, f3, f4, fp)

            ebuf(i, 1) = e2
            ebuf(i, 2) = e3
            ebuf(i, 3) = e4
            ebuf(i, 4) = ep
        end do

        call mw%allreduce('sum', ebuf)

    end if

    if (mw%talk) write (*, *) '... calculated energies'

end block energy

! Modified from anharmonic_free_energy
epotthings: block
    real(r8), dimension(3, 5) :: cumulant
    real(r8), dimension(:, :), allocatable :: ediff
    real(r8) :: inverse_kbt, U0, total_energy, hartree_to_mev, to_mev_per_atom
    integer :: i, u

    to_mev_per_atom = 1000*lo_Hartree_to_eV / real(ss%na, r8)

    if (.not. generate_configs) then

        if (mw%talk) then
            u = open_file('out', 'outfile.energies')
            write (u, '(A,A)') '# Unit:      ', 'meV/atom'
            write (u, '(A,A)') '# no. atoms: ', tochar(ss%na)
            write (u, *) '# U0 [meV/atom]: ', U0
            write (u, "(A)") '#  conf      Etotal_actual            Etotal_tdep                Epolar              &
                &Epair               Etriplet            Equartet'
            
            do i = 1, sim%nt
                total_energy = sum(ebuf(i,:))*to_mev_per_atom + U0
                write (u, "(1X,I8,6(2X,E20.12))") i, sim%stat%potential_energy(i)*to_mev_per_atom, total_energy, ebuf(i, 4)*to_mev_per_atom, &
                                                    ebuf(i, 1)*to_mev_per_atom, ebuf(i, 2)*to_mev_per_atom, &
                                                    ebuf(i, 3)*to_mev_per_atom
            end do

            ! close outfile.energies
            close (u)

            write (*, '(A)') ' ... energies writen to `outfile.energies`'
        end if
    
    ! Don't have sim%stat%potential energy when generating canonical_configs
    else
        if (mw%talk) then
            u = open_file('out', 'outfile.energies')
            write (u, '(A,A)') '# Unit:      ', 'meV/atom'
            write (u, *) '# Temperature (K) : ', opts%temperature
            write (u, '(A,A)') '# no. atoms: ', tochar(ss%na)
            write (u, "(A)") '#  conf      Ekinetic            Epolar              &
                &Epair               Etriplet            Equartet'
            
            do i = 1, opts%nconf
                write (u, "(1X,I8,5(2X,E20.12))") i, ebuf(i, 5)*to_mev_per_atom, ebuf(i, 4)*to_mev_per_atom, ebuf(i, 1)*to_mev_per_atom, &
                                                    ebuf(i, 2)*to_mev_per_atom, ebuf(i, 3)*to_mev_per_atom
            end do

            close (u)

            write (*, '(A)') ' ... energies writen to `outfile.energies`'
        end if
    end if
    

end block epotthings

call mw%destroy()

end program