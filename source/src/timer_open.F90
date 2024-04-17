!----------------------------------------------------------------------
!MODYLAS ver. 1.1.0 
!
!Copyright (c) 2014-2019 Nagoya University
!              2020-2023 The University of Tokyo
!
!Released under the MIT license.
!see https://opensource.org/licenses/MIT
!---------------------------------------------------------------------
!MODYLAS Developers:
!Yoshimichi Andoh, Kazushi Fujimoto, Tatsuya Sakashita, Noriyuki Yoshii, 
!Zhiye Tang, Jiachao Zhang, Yuta Asano, Ryo Urano, Tetsuro Nagai, 
!Atsushi Yamada, Hidekazu Kojima, Kensuke Iwahashi, Fumiyasu Mizutani, 
!Shin-ichi Ichikawa, and Susumu Okazaki.
!----------------------------------------------------------------------
!>
!! \file
!! \brief  Subroutines to open timer.
!! \author Tatsuya Sakashita
!<
  subroutine print_params(nmax)
#include "timing.h90"
    implicit none
    integer, intent(in) :: nmax

    call maprof_setup("MODYLAS", "1.1")
    call maprof_add_section("Energy_Direct", TM_ENERGY_DIRECT)
    call maprof_add_section("FMM", TM_FMM)
    call maprof_add_section("Particle2M", TM_P2M)
    call maprof_add_section("M2M", TM_M2M)
    call maprof_add_section("FMM_Ewald", TM_FMM_EWALD)
    call maprof_add_section("M2L", TM_M2L)
    call maprof_add_section("L2L", TM_L2L)
    call maprof_add_section("L2Particle", TM_L2P)
    call maprof_add_section("FMM_comm_global", TM_FMM_COMM_GLOBAL)
    call maprof_add_section("FMM_comm_local", TM_FMM_COMM_LOCAL)
    call maprof_add_section("FMM_comm_super", TM_FMM_COMM_SUPER)
    call maprof_add_section("PME_near", TM_PME_NEAR)
    call maprof_add_section("FMM_near", TM_FMM_NEAR)
#ifdef TI0
    call maprof_add_section("PME_ss_ww", TM_PME_ss_ww)
    call maprof_add_section("BG_ENE", TM_BG_ENE)
    call maprof_add_section("BG_BCG_SURF", TM_BCG_SURF)
#endif
    call maprof_add_section("PME", TM_PME)
    call maprof_add_section("All", TM_ALL)
!
    call maprof_add_section("BND_com", TM_COM_BD)
    call maprof_add_section("P2P_com3", TM_COM_D3)
    call maprof_add_section("P2P_com2", TM_COM_D2)
!debug
    call maprof_add_section("BND_com_pre", TM_COM_BD_pre)
    call maprof_add_section("BND_com_pst", TM_COM_BD_pst)

    call maprof_profile_add_problem_size("nmax", nmax)

  end subroutine print_params
