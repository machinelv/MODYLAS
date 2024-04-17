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
!! \brief  Module to store fmm parameter (nmax).
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to store fmm parameter (nmax).
!! \author Noriyuki Yoshii, Kensuke Iwahashi
!<
module fmm_parameters
  implicit none
  integer(4), protected :: nmax=4

contains

!>
!! \brief  Subroutine to set nmax.
!! \author Noriyuki Yoshii, Kensuke Iwahashi
!<
  subroutine fmod_set_nmax(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    nmax = ivalue
  end subroutine fmod_set_nmax

end module fmm_parameters
