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
!! \brief  Module for regular and singular solid harmonics (to be used in FMM).
!<
!----------------------------------------------------------------------
!>
!! \brief  Calculating regular and singular solid harmonics in spherical coordinates.
!! \author Tatsuya Sakashita, Noriyuki Yoshii.
!<
module regular_singular_solid_harmonics
  implicit none

contains

  complex(8) function regular_harmonics_pz(n, m, rad, csthe, phi)
    use math_functions
    implicit none
    integer, intent(in) :: n, m
    real(8), intent(in) :: rad, csthe, phi
    complex(8) :: zphi

    if(n.lt.iabs(m)) then
       regular_harmonics_pz = 0
       return
    endif

    zphi=dcmplx(0.d0,m*phi)
    regular_harmonics_pz = plgndr(n,m,csthe)*exp(zphi) * rad**n / factorial(n+m)
  end function regular_harmonics_pz

  subroutine regular_harmonics_pz_reals(n, m, rad, csthe, phi, re, im)
    use math_functions
    implicit none
    integer, intent(in) :: n, m
    real(8), intent(in) :: rad, csthe, phi
    real(8), intent(out) :: re, im
    real(8) :: mphi, val

    if(n.lt.iabs(m)) then
       re = 0
       im = 0
       return
    endif

    mphi = m*phi
    val = plgndr(n,m,csthe) * rad**n / factorial(n+m)
    re = val * cos(mphi)
    im = val * sin(mphi)
  end subroutine regular_harmonics_pz_reals

  complex(8) function singular_harmonics_pz(n, m, rad, csthe, phi)
    use math_functions
    implicit none
    integer, intent(in) :: n, m
    real(8), intent(in) :: rad, csthe, phi
    complex(8) :: zphi

    if(n.lt.iabs(m)) then
       singular_harmonics_pz = 0
       return
    endif

    zphi=dcmplx(0.d0,m*phi)
    singular_harmonics_pz = plgndr(n,m,csthe)*exp(zphi) * (-1)**(n+m) * factorial(n-m) / rad**(n+1)

  end function singular_harmonics_pz
  
  subroutine singular_harmonics_pz_reals(n, m, rad, csthe, phi, re, im)
    use math_functions
    implicit none
    integer, intent(in) :: n, m
    real(8), intent(in) :: rad, csthe, phi
    real(8), intent(out) :: re, im
    real(8) :: mphi, val
    
    if(n.lt.iabs(m)) then
       re = 0
       im = 0
       return
    endif
    
    mphi = m*phi
    val = plgndr(n,m,csthe) * (-1)**(n+m) * factorial(n-m) / rad**(n+1)
    re = val * cos(mphi)
    im = val * sin(mphi)
  end subroutine singular_harmonics_pz_reals

  subroutine regular_harmonics_reals(n, m, rad, csthe, phi, re, im)
    implicit none
    integer, intent(in) :: n, m
    real(8), intent(in) :: rad, csthe, phi
    real(8), intent(out) :: re, im

    call regular_harmonics_pz_reals(n, iabs(m), rad, csthe, phi, re, im)
    if (m < 0) then
       re = re * (-1)**m
       im = - im * (-1)**m
    endif
  end subroutine regular_harmonics_reals

  subroutine singular_harmonics_reals(n, m, rad, csthe, phi, re, im)
    implicit none
    integer, intent(in) :: n, m
    real(8), intent(in) :: rad, csthe, phi
    real(8), intent(out) :: re, im

    call singular_harmonics_pz_reals(n, iabs(m), rad, csthe, phi, re, im)
    if (m < 0) then
       re = re * (-1)**m
       im = - im * (-1)**m
    endif
  end subroutine singular_harmonics_reals

  complex(8) function regular_harmonics(n, m, rad, csthe, phi)
    implicit none
    integer, intent(in) :: n, m
    real(8), intent(in) :: rad, csthe, phi

    if (m >= 0) then
       regular_harmonics = regular_harmonics_pz(n, m, rad, csthe, phi)
    else
       regular_harmonics = (-1)**m * conjg(regular_harmonics_pz(n, iabs(m), rad, csthe, phi))
    endif
  end function regular_harmonics

  complex(8) function singular_harmonics(n, m, rad, csthe, phi)
    implicit none
    integer, intent(in) :: n, m
    real(8), intent(in) :: rad, csthe, phi

    if (m >= 0) then
       singular_harmonics = singular_harmonics_pz(n, m, rad, csthe, phi)
    else
       singular_harmonics = (-1)**m * conjg(singular_harmonics_pz(n, iabs(m), rad, csthe, phi))
    endif
  end function singular_harmonics

end module regular_singular_solid_harmonics
