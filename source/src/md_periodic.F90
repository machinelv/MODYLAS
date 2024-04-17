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
!! \brief  Module to store periodic type (method to calculate long-range
!!         part of Coulombic interaction) information.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to store periodic type information.
!! \author Kensuke Iwahashi
!<
module md_periodic
  implicit none
  integer(4),parameter :: CLUSTER=0,SIMPLE=1,EWALD=2,PMEWALD=3,FMM=4
  integer(4) :: md_periodic__type=-1

contains

  subroutine fmod_md_periodic__type(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    md_periodic__type = ivalue
  end subroutine fmod_md_periodic__type

end module md_periodic
