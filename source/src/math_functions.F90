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
!! \brief  Module to store mathmatical functions/subroutines.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to store mathmatical functions/subroutines.
!<
module math_functions
  implicit none

contains

  function plgndr(l,m,x)
    implicit none
    integer, intent(in) :: l,m
    real(8), intent(in) :: x
    real(8) :: plgndr
    integer i,ll
    real(8) :: fact,pll=0d0,pmm,pmmp1,somx2

    if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.d0) then
       write(0,*) 'bad arguments in plgndr'
       stop
    endif
    pmm=1.d0
    if(m.gt.0) then
       somx2=dsqrt((1.d0-x)*(1.d0+x))
       fact=1.d0
       do i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.d0
       enddo
       endif
       if(l.eq.m) then
          plgndr=pmm
       else
          pmmp1=x*(2*m+1)*pmm
          if(l.eq.m+1) then
             plgndr=pmmp1
          else
             do ll=m+2,l
                pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
                pmm=pmmp1
                pmmp1=pll
             enddo
          plgndr=pll
       endif
    endif

  end function plgndr
!**********************************************************************
  function factorial(n)
    implicit none
    integer, intent(in) :: n
    real(8) factorial
    integer :: i

    factorial = 1.d0
    do i=1,n
       factorial = factorial * i
    enddo

  end function factorial
!**********************************************************************
  subroutine cart2angle(xta,yta,zta,rad,the,csthe,phi)
    use math_const
    implicit none
    real(8), intent(in) :: xta,yta,zta
    real(8), intent(out) :: rad,the,csthe,phi
    !******** calclate angles
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
             if(yta.ge.0.d0) phi=phi+pi
             if(yta.lt.0.d0) phi=phi-pi
          endif
       else
          if(yta.gt.0.d0) phi=0.5d0*pi
          if(yta.lt.0.d0) phi=-0.5d0*pi
       endif
    else
       the=0.5d0*pi
       csthe=0.d0
       if(xta.ne.0.d0) then
          phi=datan(yta/xta)
          if(xta.le.0.d0)then
             if(yta.ge.0.d0) phi=phi+pi
             if(yta.lt.0.d0) phi=phi-pi
          endif
       else
          if(yta.gt.0.d0) phi=0.5d0*pi
          if(yta.lt.0.d0) phi=-0.5d0*pi
       endif
    endif
  end subroutine cart2angle
!**********************************************************************
  function fa(n,m)
    implicit none
    integer, intent(in) :: n,m
    real(8) fa
    fa=(-1)**n/dsqrt(factorial(n-m)*factorial(n+m))
  end function fa
!**********************************************************************
  function gammp(a,x)
    real(8), intent(in) :: a,x
    real(8) :: gammp
    ! uses gcf,gser
    real(8) gammcf,gamser,gln
    if(x.lt.0.d0.or.a.le.0.d0) then
       write(0,*) 'bad arguments in gammp'
       stop
    endif
    if(x.lt.a+1.d0)then
       call gser(gamser,a,x,gln)
       gammp=gamser
    else
       call gcf(gammcf,a,x,gln)
       gammp=1.-gammcf
    endif
  end function gammp
!**********************************************************************
  subroutine gcf(gammcf,a,x,gln)
    real(8), intent(out) :: gammcf
    real(8), intent(in) :: a,x
    real(8), intent(out) :: gln
    integer itmax
    real(8) eps,fpmin
    parameter (itmax=100,eps=3.d-7,fpmin=1.d-30)
    ! uses gammln
    integer i
    real(8) an,b,c,d,del,h
    gln=gammln(a)
    b=x+1.d0-a
    c=1.d0/fpmin
    d=1.d0/b
    h=d
    do i=1,itmax
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
    enddo
    write(0,*) 'a too large, itmax too small in gcf'
    stop
1   gammcf=exp(-x+a*log(x)-gln)*h
  end subroutine gcf
!**********************************************************************
  subroutine gser(gamser,a,x,gln)
    real(8), intent(out) :: gamser
    real(8), intent(in) :: a, x
    real(8), intent(out) :: gln
    integer itmax
    real(8) eps
    parameter (itmax=100,eps=3.d-7)
    ! uses gammln
    integer n
    real(8) ap,del,sum
    gln=gammln(a)
    if(x.le.0.d0)then
       if(x.lt.0.d0) then
          write(0,*) 'x < 0 in gser'
          stop
       end if
       gamser=0.d0
       return
    endif
    ap=a
    sum=1.d0/a
    del=sum
    do n=1,itmax
       ap=ap+1.d0
       del=del*x/ap
       sum=sum+del
       if(abs(del).lt.abs(sum)*eps)goto 1
    enddo
    write(0,*) 'a too large, itmax too small in gser'
    stop
1   gamser=sum*exp(-x+a*log(x)-gln)
  end subroutine gser
!**********************************************************************
  function gammln(xx)
    real(8), intent(in) :: xx
    real(8) :: gammln
    integer j
    double precision ser,stp,tmp,x,y,cof(6)
    save cof,stp
    data cof,stp/76.18009172947146d0,-86.50532032941677d0, &
         & 24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
         & -.5395239384953d-5,2.5066282746310005d0/
    x=xx
    y=x
    tmp=x+5.5d0
    tmp=(x+0.5d0)*log(tmp)-tmp
    ser=1.000000000190015d0
    do j=1,6
       y=y+1.d0
       ser=ser+cof(j)/y
    enddo
    gammln=tmp+log(stp*ser/x)
  end function gammln
!**********************************************************************
  function gammq(a,x)
    real(8), intent(in) :: a,x
    real(8) :: gammq
    ! uses gcf,gser
    real(8) gammcf,gamser,gln
    if(x.lt.0..or.a.le.0.) then
       write(0,*) 'bad arguments in gammq'
       stop
    endif
    if(x.lt.a+1.)then
       call gser(gamser,a,x,gln)
       gammq=1.-gamser
    else
       call gcf(gammcf,a,x,gln)
       gammq=gammcf
    endif
  end function gammq

end module math_functions
