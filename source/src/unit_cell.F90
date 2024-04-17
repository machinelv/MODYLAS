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
!! \brief  Module and subroutines to store unit cell geometry.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to store unit cell geometry.
!! \author Kensuke Iwahashi
!<
module unit_cell
  implicit none
  real(8) :: cellx=0.0d0, celly=0.0d0, cellz=0.0d0
  real(8) :: cellxh=0.0d0, cellyh=0.0d0, cellzh=0.0d0
  real(8) :: alpha=90.0d0, beta=90.0d0, gamma=90.0d0
  real(8) :: cellvol=0.0d0
  real(8) :: sinbeta=1.0d0, cosbeta=0.0d0
  real(8) :: singamma=1.0d0, cosgamma=0.0d0
  real(8) :: b_factor=0.0d0, g_factor=1.0d0
  real(8) :: r_singamma=1.0d0, r_g_factor=1.0d0
  logical :: cuboid=.true.

contains

  subroutine fmod_set_alpha(value)
    implicit none
    real(8), intent(in) :: value
    alpha = value
    if (abs(value-90.0) .gt. 1.0d-5)  cuboid = .false.
  end subroutine fmod_set_alpha
!----------------------------------------------------------------------
  subroutine fmod_set_beta(value)
    implicit none
    real(8), intent(in) :: value
    beta = value
    if (abs(value-90.0) .gt. 1.0d-5)  cuboid = .false.
  end subroutine fmod_set_beta
!----------------------------------------------------------------------
  subroutine fmod_set_gamma(value)
    implicit none
    real(8), intent(in) :: value
    gamma = value
    if (abs(value-90.0) .gt. 1.0d-5)  cuboid = .false.
  end subroutine fmod_set_gamma
!----------------------------------------------------------------------
  subroutine fmod_cellx(value)
    implicit none
    real(8), intent(in) :: value
    cellx = value
  end subroutine fmod_cellx
!----------------------------------------------------------------------
  subroutine fmod_celly(value)
    implicit none
    real(8), intent(in) :: value
    celly = value
  end subroutine fmod_celly
!----------------------------------------------------------------------
  subroutine fmod_cellz(value)
    implicit none
    real(8), intent(in) :: value
    cellz = value
  end subroutine fmod_cellz
!----------------------------------------------------------------------
  subroutine fmod_cellxh(value)
    implicit none
    real(8), intent(in) :: value
    cellxh = value
  end subroutine fmod_cellxh
!----------------------------------------------------------------------
  subroutine fmod_cellyh(value)
    implicit none
    real(8), intent(in) :: value
    cellyh = value
  end subroutine fmod_cellyh
!----------------------------------------------------------------------
  subroutine fmod_cellzh(value)
    implicit none
    real(8), intent(in) :: value
    cellzh = value
  end subroutine fmod_cellzh

end module unit_cell
