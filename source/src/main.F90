!----------------------------------------------------------------------
!MODYLAS ver. 1.1.0 
!
!Copyright (c) 2014-2019 Nagoya University
!              2020-2023 The University of Tokyo
!
!Released under the MIT license.
!see https://opensource.org/licenses/MIT
!----------------------------------------------------------------------
!MODYLAS Developers:
!Yoshimichi Andoh, Kazushi Fujimoto, Tatsuya Sakashita, Noriyuki Yoshii, 
!Zhiye Tang, Jiachao Zhang, Yuta Asano, Ryo Urano, Tetsuro Nagai, 
!Atsushi Yamada, Hidekazu Kojima, Kensuke Iwahashi, Fumiyasu Mizutani, 
!Shin-ichi Ichikawa, and Susumu Okazaki.
!----------------------------------------------------------------------
!>
!! \file
!! \brief Main routine of MODYLAS
!! \author Yoshimichi Andoh
!<
!======================================================================
!===                                                                ===
!===     MODYLAS --- MOlecular DYnamics simulation software         ===
!===                 ^^        ^^                                   ===
!===                 for LArge Systems                              ===
!===                     ^^    ^                                    ===
!===   ! SI unit is used if no explicit specification               ===
!===                                                                ===
!======================================================================
program modylas
  use commandline_args
  use version
  use atom_virial
  use file_mdrun
  use md_condition
  use md_monitors
  use md_multiplestep
  use mpi_tool
  use comm_direct2_dr, only : comm_direct_2_dr
  use comm_direct3_dr, only : comm_direct_3_dr
  use trajectory_mpi
  use shake_rattle_roll
  use subcell
  use cell_shape
  use boundary
  use nve_mod
  use nvt_mod
  use opt_mod
  use npt_a_mod
  use npt_pr_mod
  use npt_z_mod
  use force_wrap, only : md_calculate_forces
  use input_mod
  use display_log_mod
  use ensemble_numbers
  use file_application
  use file_mdmntr
  use file_mdff, only : write_mdffbin_restart
#ifdef TIP4
  use tip4p, only : distribute_force_on_msite
#endif
  use dist_atom2cell, only : atom2cell_wrap
  use openmp_tool
#include "timing.h90"
  use system_dipole, only: calculating_system_dipole_in_each_process, &
&                          bSystemDipole, close_system_dipole
  implicit none
  include 'mpif.h'
  integer(4) :: i0,k0

  call mpistart()
  call initialize_openmp()

  call parse_args()
  call parse_input(.true.)
  call initialize_application()

  if(myrank==0) write(*,*) 'Done up to initialize_application'
  flush(6)

  call output_log   ! output log to stdout
  TIME_INIT

  !
  !     prepare initial state for MD calculation
  !
  call cell_edge()      
  call calc_ia2c()      
  call atom2cell_wrap() 
  !
  if(totnconst > 0) call update_shake_local()
  !
  call comm_direct_3_dr(wkxyz) 
  call apply_pbc()
  !
  MTl=maxMTl ! Initial value for 0th-step
  MTm=maxMTm ! Initial valur for 0th-step
  call md_calculate_forces()

  if(myrank.eq.mpiout) call record_current_runtime_info()

#ifdef DEBUGFCE
   call calc_hamiltonian_nve()
   if(myrank==mpiout) call record_current_monitors
   if(myrank==0)then
    write(*,*) '***********************************'
    write(*,*) 'MODYLAS stoped at 0th step, because'
    write(*,*) '-DDEBUGFCE is set.                 '
    write(*,*) '***********************************'
   endif
   call mpiend()
#endif

  !for check 0 step potential
  if(is_stopped_at_0step)then
    call calc_hamiltonian_nve()
    if(myrank==mpiout) call record_current_monitors
    if(myrank==0)then
    write(*,*) '***********************************'
    write(*,*) 'MODYLAS stoped at 0th step, because'
    write(*,*) '0step_stop=yes is set in .mddef.   '
    write(*,*) '***********************************'
    endif
    call mpiend
  endif
  
  TIME_PROFILE_START
  TIME_BARRIER(TMB_PRE_MAINLOOP)
  TIME_PROFILE_STOP
  
  if(myrank==0) write(*,*) 'Starting MD main loop'
  flush(6)

  !
  !     MD main loop
  !
  
  do while (mdstep<md_condition__howmany_steps)

     if(mdstep >=  10)  TIME_PROFILE_START
     TIME_BARRIER(TMB_MDSTEP_TOP)
     TIME_START(TM_ALL)

     if(    md_condition__ensemble==NVE) then
        TIME_START(TM_NVE)
        call nve_integrate()
        TIME_STOP(TM_NVE)
     elseif(md_condition__ensemble==NVT .or. & 
    &       md_condition__ensemble==NVLXLYLZT) then
        call nvt_integrate()
     elseif(md_condition__ensemble==NPT_A) then
        TIME_START(TM_NPT_A)
        call npt_a_integrate()
        TIME_STOP(TM_NPT_A)
     elseif(md_condition__ensemble==NPT_Z .or. &
    &       md_condition__ensemble==NPTLZ .or. &
    &       md_condition__ensemble==NLXLYPZT .or. &
    &       md_condition__ensemble==NPXPYPZT .or. &
    &       md_condition__ensemble==NLXLYLZT .or. &
    &       md_condition__ensemble==NPTLZxy) then
        call npt_z_integrate()
     elseif(md_condition__ensemble==NPT_PR) then
        call npt_pr_integrate()
     elseif(md_condition__ensemble==OPT) then
        call opt_integrate()
     endif

     TIME_STOP(TM_ALL)
     TIME_PROFILE_STOP
     
     if(    md_condition__ensemble==NVE) then
        call calc_hamiltonian_nve()
     elseif(md_condition__ensemble==NVT .or. &
          & md_condition__ensemble==NVLXLYLZT) then
        call calc_hamiltonian_nvt()
     elseif(md_condition__ensemble==NPT_A) then
        call calc_hamiltonian_npt_a()
     elseif(md_condition__ensemble==NPT_Z .or. &
          &       md_condition__ensemble==NPTLZ .or. &
          &       md_condition__ensemble==NLXLYPZT .or. &
          &       md_condition__ensemble==NPXPYPZT .or. &
          &       md_condition__ensemble==NLXLYLZT .or. &
          &       md_condition__ensemble==NPTLZxy) then
        call calc_hamiltonian_npt_z()
     elseif(md_condition__ensemble==NPT_PR) then
        call calc_hamiltonian_npt_pr()
     elseif(md_condition__ensemble==OPT) then
        call calc_potential_opt()
     endif
    if(bSystemDipole) then
      call calculating_system_dipole_in_each_process(mdstep + 1)
    endif
     
    if (volume_scaling)       call apply_volume_scaling()       ! default:off
    if (cellgeometry_scaling) call apply_cellgeometry_scaling() ! default:off
    if (velocity_scaling)     call apply_velocity_scaling()     ! default:off
     
     mdstep=mdstep+1

     if(mdstep >= 10)  TIME_PROFILE_START

     TIME_BARRIER(TMB_IO_RECCURSTATE)
     TIME_START(TM_IO_RECCURSTATE)
     call record_current_state()
     TIME_STOP(TM_IO_RECCURSTATE)
     TIME_PROFILE_STOP

  enddo

  if(bSystemDipole) then
    call close_system_dipole()
  endif

  if( allow_bond_breaking >= 0.0d0 ) then
    call write_mdffbin_restart
  endif
  
  if(myrank.eq.mpiout) call cleanup()
  
  if(md_condition__ensemble .ne. OPT) then
     if(myrank==0)then
        write(*,*) '***********************'
        write(*,*) 'MODYLAS normally ended!'
        write(*,*) '***********************'
        write(*,*) '<P_inner>=             '
        write(*,'(es20.12,a5)')  (sum_P_part(1))/dble(icount_vir),' [Pa]'
        write(*,*) '<P_outer>=             '
        write(*,'(es20.12,a5)')  (sum_P_part(2))/dble(icount_vir),' [Pa]'
        write(*,*) '<P_outermost>=         '
        write(*,'(es20.12,a5)')  (sum_P_part(3))/dble(icount_vir),' [Pa]'
        write(*,*) '***********************'
        write(*,*) '<P_total>=             '
        write(*,'(es20.12,a5)')   (sum_P)/dble(icount_vir),' [Pa]'
        write(*,*) '***********************'
     endif
  endif
  !NI-->
  if(md_condition__ensemble==NPT_PR) then
     if(myrank==0)then
        call cell_convert2(dummy(1),dummy(2),dummy(3),dummy(4),dummy(5),dummy(6),sum_box/dble(mdstep))
        write(*,*) '////////////////////////'
        write(*,*) ' Average of Cell Box [A]'
        write(*,'(1a4,1e19.12)') ' lx=', dummy(1)
        write(*,'(1a4,1e19.12)') ' ly=', dummy(2)
        write(*,'(1a4,1e19.12)') ' lz=', dummy(3)
        write(*,'(1a7,1e19.12)') ' alpha=', dummy(4)
        write(*,'(1a6,1e19.12)') ' beta=', dummy(5)
        write(*,'(1a7,1e19.12)') ' gamma=', dummy(6)
        write(*,*) '////////////////////////'
     endif
  endif
  !<--NI

  TIME_FINALIZE
  call mpiend()

end program modylas
