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
!! \brief  Module and subroutines to set random seed.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to set random seed.
!! \author Kensuke Iwahashi
!<
module random_seed
  implicit none
  integer(4) :: randomseed=1235

contains

  subroutine fmod_randomseed(ivalue) ! only for input version=0.9.0
    implicit none
    integer(4), intent(in) :: ivalue
    randomseed = ivalue
  end subroutine fmod_randomseed

end module random_seed
