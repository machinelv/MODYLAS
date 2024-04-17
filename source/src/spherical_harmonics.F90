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
!! \brief  Module and subroutines which relates to the spherical harmonics.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module which relates to the spherical harmonics.
!! \author Noriyuki Yoshii, Kensuke Iwahashi
!<
module spherical_harmonics
  implicit none
  real(8), parameter :: eps=3.d-7
  integer, parameter :: itmax=100

contains

!----------------------------------------------------------------------
!>
!! \brief Subroutine which converts Cartesian coordinate into spherical coordinate
!! \author Noriyuki Yoshii, Kensuke Iwahashi
!<
  subroutine cart2angle(xta,yta,zta,rad,the,csthe,phi)
!----------------------------------------------------------------------
    use md_const, only : PI
    implicit none
    real(8), intent(in) :: xta,yta,zta
    real(8), intent(out) :: rad,the,csthe,phi
    !******* calclate angles
    rad=dsqrt(xta*xta+yta*yta+zta*zta)
    if(rad.eq.0.d0) then
       the=0.d0
       phi=0.d0
       csthe=1.d0
    elseif(zta.ne.0.d0) then
       the=dacos(zta/rad)
       csthe=zta/rad
       if(xta.ne.0.d0) then
          phi=datan(yta/xta)
          if(xta.le.0.d0)then
             if(yta.ge.0.d0) phi=phi+PI
             if(yta.lt.0.d0) phi=phi-PI
          endif
       else
          if(yta.gt.0.d0) then
             phi=0.5d0*PI
          elseif(yta.lt.0.d0) then
             phi=-0.5d0*PI
          else
             phi = 0  ! to suppress warning of undefined value
          endif
       endif
    else
       the=0.5d0*PI
       csthe=0.d0
       if(xta.ne.0.d0) then
          phi=datan(yta/xta)
          if(xta.le.0.d0)then
             if(yta.ge.0.d0) phi=phi+PI
             if(yta.lt.0.d0) phi=phi-PI
          endif
       else
          if(yta.gt.0.d0) phi=0.5d0*PI
          if(yta.lt.0.d0) phi=-0.5d0*PI
       endif
    endif

  end subroutine cart2angle
!*********************************************************************
!>
!! \file
!! \brief Subroutine which calculate Legendre polynomials
!! \author Noriyuki Yoshii
!<
  function algndr(n,m,x)
!*********************************************************************
    implicit none
    integer(4), intent(in) :: n,m
    real(8), intent(in) ::x
    real(8) :: algndr
    real(8) :: apnn2,apnn1,p=0d0
    integer(4) :: i,k

    apnn2=1.d0
    if(m.gt.0) then
       p=dsqrt(1.d0-x*x)
       do i=1,m
          apnn2=-apnn2*(2.d0*i-1.d0)*p
       enddo
    endif
    if(n.eq.m) then
       algndr=apnn2
    else
       apnn1=x*(2*m+1)*apnn2
       if(n.eq.m+1) then
          algndr=apnn1
       else
          do k=m+2,n
             p=(x*(2*k-1)*apnn1-(k+m-1)*apnn2)/(k-m)
             apnn2=apnn1
             apnn1=p
          enddo
          algndr=p
       endif
    endif
  end function algndr

!*********************************************************************
!>
!! \file
!! \brief Subroutine which calculate n factorial
!! \author Noriyuki Yoshii
!<
  function factorial(n)
!*********************************************************************
    implicit none
    integer, intent(in) :: n
    real(8) :: factorial
    integer :: i
    factorial = 1.d0
    do i=1,n
       factorial = factorial * i
    enddo
  end function factorial
!*********************************************************************
!>
!! \file
!! \brief Subroutine which calculate incomplete gamma function
!! \author Noriyuki Yoshii
!<
  function incmpgamm(a,x)
!*********************************************************************
    real(8), intent(in) :: a,x
    real(8) :: incmpgamm
    if(x.lt.a+1.d0)then
       incmpgamm=kyu(a,x)
    else
       incmpgamm=1.d0-ren(a,x)
    endif
  end function incmpgamm
!*********************************************************************
  function ren(a,x)
!*********************************************************************
    real(8), intent(in) :: a,x
    real(8) :: ren,gln
    real(8), parameter :: fpmin=1.d-30
    integer :: i
    real(8) :: an,b,c,d,del,h
    gln=lngamma(a)
    b=x+1.d0-a
    c=1.d0/fpmin
    d=1.d0/b
    h=d
    do 11 i=1,itmax
       an=-i*(i-a)
       b=b+2.d0
       d=an*d+b
       if(abs(d).lt.fpmin)d=fpmin
       c=b+an/c
       if(abs(c).lt.fpmin)c=fpmin
       d=1.d0/d
       del=d*c
       h=h*del
       if(abs(del-1.d0).lt.eps)goto 1
11     continue
1      ren=exp(-x+a*log(x)-gln)*h

  end function ren
!*********************************************************************
  function kyu(a,x)
!*********************************************************************
    real(8), intent(in) :: a,x
    real(8) :: kyu,gln
    integer :: n
    real(8) :: ap,del,sum
    gln=lngamma(a)
    if(x.le.0.d0)then
       kyu=0.d0
       return
    endif
    ap=a
    sum=1.d0/a
    del=sum
    do 11 n=1,itmax
       ap=ap+1.d0
       del=del*x/ap
       sum=sum+del
       if(abs(del).lt.abs(sum)*eps)goto 1
11     continue
1      kyu=sum*exp(-x+a*log(x)-gln)
  end function kyu
!*********************************************************************
  function lngamma(xx)
!*********************************************************************
    implicit none
    real(8), intent(in) :: xx
    real(8) :: lngamma,ser,tmp

    tmp  = xx + 4.5d0
    tmp  = tmp-(xx - 0.5d0)*dlog(tmp)
    ser = 1.000000000190015d0 &
         &     + (76.18009172947146d0   / xx) &
         &     - (86.50532032941677d0   / (xx + 1.0d0)) &
         &     + (24.01409824083091d0   / (xx + 2.0d0)) &
         &     - (1.231739572450155d0   / (xx + 3.0d0)) &
         &     + (0.1208650973866179d-2 / (xx + 4.0d0)) &
         &     - (0.5395239384953d-5    / (xx + 5.0d0))

    lngamma=dlog(2.5066282746310005d0 * ser) - tmp
  end function lngamma
!*********************************************************************

end module spherical_harmonics

