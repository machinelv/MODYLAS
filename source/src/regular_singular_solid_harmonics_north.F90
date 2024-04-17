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
!! \brief  Module to calculate singular solid harmonics at north pole to be used in M2L of FMM.
!<
!----------------------------------------------------------------------
!>
!! \brief  Calculating singular solid harmonics at north pole (with respect to the basis of rotated spherical coordinates).
!! \author Tatsuya Sakashita, Noriyuki Yoshii.
!<
module regular_singular_solid_harmonics_north
  implicit none

contains

  real(8) function regular_harmonics_north(n, rad)
    use math_functions, only : factorial
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: rad

    regular_harmonics_north = rad**n / factorial(n+0)
  end function regular_harmonics_north

  real(8) function singular_harmonics_north(n, rad)
    use math_functions, only : factorial
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: rad

    singular_harmonics_north = (-1)**(n+0) * factorial(n-0) / rad**(n+1)
  end function singular_harmonics_north

end module regular_singular_solid_harmonics_north
