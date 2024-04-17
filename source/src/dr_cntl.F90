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
!! \brief  Module and subroutines to control MTD method. 
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to control MTD method.
!! \author Yoshimichi Andoh, Shi-ichi Ichikawa
!<
!----------------------------------------------------------------------
module dr_cntl
    integer, parameter :: nbd=1, nbd2=2*nbd, nadj=2
    integer, parameter :: mbd =2, mbdp1 =mbd+1, mbd2 =2*mbd
    integer, parameter :: mtbd=3, mtbdp2=mtbd+2, mtbd2=2*mtbd
!
    integer(4), allocatable :: ldr_zp(:,:)
    integer(4), allocatable :: ldr_zm(:,:)
    integer(4), allocatable :: ldr_yp(:,:)
    integer(4), allocatable :: ldr_ym(:,:)
    integer(4), allocatable :: ldr_xp(:,:)
    integer(4), allocatable :: ldr_xm(:,:)
!
    integer(4), allocatable :: slist(:,:)  ! copy from edirect_dr

contains

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to initialize arrays for P2P. 
!! \author Yoshimichi Andoh, Shi-ichi Ichikawa
!<
   subroutine init_drtbl
!----------------------------------------------------------------------
   use domain, only : lxdiv, lydiv, lzdiv
   implicit none

   allocate(ldr_zp(1-nbd:lydiv+nbd, 1-nbd:lxdiv+nbd))
   allocate(ldr_zm(1-nbd:lydiv+nbd, 1-nbd:lxdiv+nbd))
   allocate(ldr_yp(1-nbd:lxdiv+nbd, -nadj:nadj))
   allocate(ldr_ym(1-nbd:lxdiv+nbd, -nadj:nadj))
   allocate(ldr_xp(1-nbd:lydiv+nbd, -nadj:nadj))
   allocate(ldr_xm(1-nbd:lydiv+nbd, -nadj:nadj))

   ldr_zm = -1
   ldr_zp = -1
   ldr_ym = -1
   ldr_yp = -1
   ldr_xm = -1
   ldr_xp = -1

   end subroutine init_drtbl
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to get number of subcells in the halo. 
!! \author Yoshimichi Andoh, Shi-ichi Ichikawa
!<
  subroutine get_halo_number(ip0,np,nl,nscell, &
     &                  nbd_m,nbd_p,mbd0)
!----------------------------------------------------------------------
  implicit none
  integer(4),intent(in) :: ip0          ! variable
  integer(4),intent(in) :: np,nl,nscell ! constant
  integer(4),intent(out) :: nbd_m,nbd_p
  integer(4),intent(out) :: mbd0
  integer(4) npow,icg0,icg1
  integer(4) ip

  if(    ip0 .gt. np-1)then
    ip=ip0-np  ! periodic boundary +
  elseif(ip0 .lt. 0   )then
    ip=ip0+np  ! periodic boundary -
  else
    ip=ip0
  endif

  npow = nl
  icg0=(nscell *  ip       )/np + 1
  icg1=(nscell * (ip+1) - 1)/np + 1
  if(npow==2)then
     nbd_m = mbdp1 - mod(icg0   ,npow)
     nbd_p = mbd   + mod(icg1   ,npow)
     mbd0  = mbd
  else
     nbd_m = mtbd   + mod(icg0+2,npow)
     nbd_p = mtbdp2 - mod(icg1+2,npow)
     mbd0  = mtbd
  endif

  end subroutine get_halo_number
!----------------------------------------------------------------------
end module dr_cntl
