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
!! \brief  Subroutine to finalize timer.
!! \author Tatsuya Sakashita
!<
  subroutine timer_close()
    use mpi_tool, only : myrank
#include "timing.h90"
    implicit none

    if (myrank == 0)  write(*,'(/,a)') '**** timings'
    if (myrank == 0)  write(*, '(a)') '---------------'
    call maprof_print_time_mpi(TM_ALL, 'Total of MODYLAS: ')
    if (myrank == 0)  write(*, '(a)') '---------------'
    call maprof_print_time_mpi(TM_ENERGY_DIRECT, 'P2P:           ')
    call maprof_print_time_mpi(TM_COM_D3,       'P2P_comm3 :     ')
    call maprof_print_time_mpi(TM_COM_D2,       'P2P_comm2 :     ')
    if (myrank == 0)  write(*, '(a)') '---------------'
    call maprof_print_time_mpi(TM_FMM, 'Total of FMM: ')
    if (myrank == 0)  write(*, '(a)') '---------------'
    call maprof_print_time_mpi(TM_P2M,           'Particle2M:    ')
    call maprof_print_time_mpi(TM_M2M,           'M2M:           ')
    call maprof_print_time_mpi(TM_FMM_EWALD,     'FMM_Ewald:     ')
    call maprof_print_time_mpi(TM_M2L,           'M2L:           ')
    call maprof_print_time_mpi(TM_L2L,           'L2L:           ')
    call maprof_print_time_mpi(TM_L2P,           'L2Particle:    ')
    if (myrank == 0)  write(*, '(a)') '---------------'
    call maprof_print_time_mpi(TM_FMM_COMM_LOCAL,        'FMM_comm_local:     ')
    call maprof_print_time_mpi(TM_FMM_COMM_SUPER,        'FMM_comm_super:     ')
    call maprof_print_time_mpi(TM_FMM_COMM_GLOBAL,       'FMM_comm_global:    ')
    if (myrank == 0)  write(*, '(a)') '---------------'
    call maprof_print_time_mpi(TM_COM_BD,       'bound_comm:     ')
!debug
!   call maprof_print_time_mpi(TM_COM_BD_pre,       'bound_comm_pre:     ')
!   call maprof_print_time_mpi(TM_COM_BD_pst,       'bound_comm_pst:     ')
!
    if (myrank == 0)  write(*, '(a)') '---------------'

    call maprof_output()

  end subroutine timer_close
