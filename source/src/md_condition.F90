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
!! \brief  Module and subroutines to controm MD calculation condition.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to store MD calculation conditions.
!! \author Kensuke Iwahashi
!<
module md_condition
  implicit none
  integer(4) :: mdstep=0
  real(8) :: dt=0.0d0,rdt=0.0d0,rdt2=0.0d0,rdtb=0.0d0,rdtb2=0.0d0
  real(8) :: dth=0.0d0
  integer(4) :: degree_of_freedom=0
  real(8) :: degree_of_freedom_inverse=0.0d0
  integer(4) :: md_condition__howmany_steps=0
  integer(4) :: md_condition__ensemble=0
  integer(4) :: md_condition__force_field=100
  integer(4) :: initial_step = 0
  logical :: velocity_scaling=.false.
  logical :: velocity_scaling_region_wise=.false.
  logical :: LJ_LRC=.true. ! default
  logical :: volume_scaling=.false. ! defalut
  logical :: cellgeometry_scaling=.false. ! defalut
!debug option
  logical:: is_stopped_at_0step=.false.

  real(8) :: deformx=0d0
  real(8) :: deformy=0d0
  real(8) :: deformz=0d0
  real(8) :: allow_bond_breaking = -1.0d0

! TZYBEGIN
  logical :: bNonEvenDihedralPotential = .true.
  ! .false. = method in https://en.wikipedia.org/wiki/Dihedral_angle
  ! .true. = cross product method
! TZYEND

contains

! TZYBEGIN
subroutine fmod_md_condition__bNonEvenDihedralPotential(ivalue)
    use mpi_tool
    implicit none
    integer(4), intent(in) :: ivalue
    if    (ivalue==1) then
       bNonEvenDihedralPotential = .true.
    elseif(ivalue==0) then
       bNonEvenDihedralPotential = .false.
    elseif(ivalue==-1) then
       if(myrank==0)then
          write(*,*) 'You did not select non_even_dihedral_potential(yes/no), so', &
               &             ' default (non_even_dihedral_potential=yes) is set.'
       endif
    endif
end subroutine

! TZYEND

subroutine fmod_md_condition__deformx(value)
  implicit none 
  real(8), intent(in) :: value
  deformx = value * dt
end subroutine fmod_md_condition__deformx

subroutine fmod_md_condition__deformy(value)
  implicit none
  real(8), intent(in) :: value
  deformy = value * dt
end subroutine fmod_md_condition__deformy

subroutine fmod_md_condition__deformz(value)
  implicit none
  real(8), intent(in) :: value
  deformz = value * dt
end subroutine fmod_md_condition__deformz

subroutine fmod_set_allow_bond_breaking(input)
  implicit none
  real(8) :: input
  allow_bond_breaking = input
end subroutine

!>
!! \brief  Subroutine to set delta t (dt).
!! \author Kensuke Iwahashi
!<
  subroutine fmod_md_condition__dt(value)
    implicit none
    real(8), intent(in) :: value
    dt = value
    dth = 0.5d0 * dt
    rdt = 1.0d0 / dt
    rdt2 = rdt * rdt
    rdtb = rdt + rdt
    rdtb2 = rdtb * rdtb
  end subroutine fmod_md_condition__dt
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to set MD steps.
!! \author Kensuke Iwahashi
!<
  subroutine fmod_howmany_steps(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    md_condition__howmany_steps = ivalue
  end subroutine fmod_howmany_steps
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to set initial MD steps.
!! \author Kensuke Iwahashi
!<
  subroutine fmod_mdstep(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    initial_step = ivalue
    mdstep = 0
  end subroutine fmod_mdstep
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to set ensemble.
!! \author Kensuke Iwahashi
!<
  subroutine fmod_md_condition__ensemble(ivalue)
    use unit_cell, only : cuboid
    use ensemble_numbers
    implicit none
    integer(4), intent(in) :: ivalue
    md_condition__ensemble = ivalue
    if (ivalue == NPT_PR)  cuboid = .false.
  end subroutine fmod_md_condition__ensemble
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to initilize degree of freedom.
!! \author Kensuke Iwahashi
!<
  subroutine init_md_condition
    use trajectory_org
    use shake_rattle_roll, only : totnconst
#ifdef TIP4
    use tip4p, only : msite
#endif
    implicit none
    integer(4) :: k

    degree_of_freedom = n * 3 - 3
#ifdef TIP4
    degree_of_freedom = degree_of_freedom - 3*msite
#endif
    if(totnconst .ne. 0) then
       degree_of_freedom = degree_of_freedom - totnconst
    endif
    degree_of_freedom_inverse = 1.0d0 / degree_of_freedom

  end subroutine init_md_condition

  subroutine fmod_velocity_scaling(ivalue, ivalue2)
    implicit none
    logical, intent(in) :: ivalue, ivalue2

    velocity_scaling =  ivalue
    velocity_scaling_region_wise = ivalue2

  end subroutine fmod_velocity_scaling
!----------------------------------------------------------------------
  subroutine fmod_volume_scaling(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    if (ivalue /= 0) then
       volume_scaling = .true.
    else
       volume_scaling = .false.
    endif
  end subroutine fmod_volume_scaling
!----------------------------------------------------------------------
  subroutine fmod_cellgeometry_scaling(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    write(*,*) 'cellgeometry_scaling called.'
    if (ivalue /= 0) then
       cellgeometry_scaling = .true.
    else
       cellgeometry_scaling = .false.
    endif
  end subroutine fmod_cellgeometry_scaling
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to set a flag for LJ long-range correction.
!! \author Kensuke Iwahashi
!<
  subroutine fmod_LJ_correction(ivalue)
    use mpi_tool
    implicit none
    integer(4), intent(in) :: ivalue
    if    (ivalue==1) then
       LJ_LRC = .true.
    elseif(ivalue==0) then
       LJ_LRC = .false.
    elseif(ivalue==-1) then
       if(myrank==0)then
          write(*,*) 'You did not select LJcorrection(yes/no), so', &
               &             ' default (LJcorrection=yes) is set.'
       endif
    endif
#ifdef CHARMMFSW
    LJ_LRC = .false.
#endif
  end subroutine fmod_LJ_correction

  subroutine fmod_set_0step_stop
  implicit none

  is_stopped_at_0step=.true.

  end subroutine fmod_set_0step_stop

end module md_condition
