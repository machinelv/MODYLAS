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
!! \brief  Module and subroutines to store the Ewald parameters.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to store the Ewald parameters.
!! \author Kensuke Iwahashi
!<
!----------------------------------------------------------------------
module ewald_variables
  implicit none
  real(8) :: ewald_alpha=3.767d+9
  logical :: ewald_sterm=.false. ! default

contains

!>
!! \brief  Subroutine to set the Ewald alpha value.
!! \author Kensuke Iwahashi
!<
  subroutine fmod_md_ewald__alpha(value)
    implicit none
    real(8), intent(in) :: value
    ewald_alpha = value
  end subroutine fmod_md_ewald__alpha
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to set the flag for the surface term.
!! \author Kensuke Iwahashi
!<
  subroutine fmod_ewald_sterm(ivalue)
    use mpi_tool
    implicit none
    integer(4), intent(in) :: ivalue
    if (ivalue /= 0) then
       ewald_sterm = .true.
    else
       ewald_sterm = .false.
    endif
  end subroutine fmod_ewald_sterm

end module ewald_variables
