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
!-DPRECISION_P2P_MIX : opens the mixed-precision of P2P subroutine.
!-DPRECISION_P2P_SP  : opens the all single-precision of P2P subroutine.
!
!-DPRECISION_M2L_MIX : opens the mixed_precision of M2L subroutine.
!-DPRECISION_M2L_SP  : opens the all single-precision of M2L subroutine.
!----------------------------------------------------------------------
!>
!! \file
!! \brief Module to control calculation precision (FP32 or FP64).
!<
!----------------------------------------------------------------------
!>
!! \brief Module to control calculation precision (FP32 or FP64).
!! \author Jiachao Zhang, Shin-ichi Ichikawa, Yoshimichi Andoh
!<
module precision
  implicit none
  integer,  public, parameter :: sp = selected_real_kind(6, 37)
  integer,  public, parameter :: dp = selected_real_kind(15, 307)

#if defined (PRECISION_P2P_MIX)
  integer,  public, parameter :: wrap= dp ! target of summation
  integer,  public, parameter :: wrp = sp ! intermediate variables/arrays
#elif defined (PRECISION_P2P_SP)
  integer,  public, parameter :: wrap= sp
  integer,  public, parameter :: wrp = sp
#else
!Default 
  integer,  public, parameter :: wrap= dp !selected_real_kind(15, 307)
  integer,  public, parameter :: wrp = dp 
#endif /** PRECISION_P2P_MIX **/

#if defined(PRECISION_PME_SP)
  integer,  public, parameter :: wrrp = sp
#else
  integer,  public, parameter :: wrrp = dp
#endif

!Detailes is discussed in merge request page.
  integer,  public, parameter :: dp_void = selected_real_kind(15, 307)
  integer,  public, parameter :: wrp_void = dp_void
  
!Default is all-dp for remove_special14lj
  integer,  public, parameter :: dp_special = selected_real_kind(15, 307)
  integer,  public, parameter :: wrp_special = dp_special

!Default is all-dp for remove_scaling14lj/cl.
  integer,  public, parameter :: dp_scale = selected_real_kind(15, 307)
  integer,  public, parameter :: wrp_scale = dp_scale
  
#if defined (PRECISION_M2L_MIX)
  integer,  public, parameter :: wvap = dp     ! for multipole moment accumulation array.
  integer,  public, parameter :: wmp = sp      ! for M2L matrix.
  integer,  public, parameter :: wvp = sp      ! for multipole moment.
#elif defined (PRECISION_M2L_SP)
  integer,  public, parameter :: wvap = sp     ! for multipole moment accumulation array.
  integer,  public, parameter :: wmp = sp      ! for M2L matrix.
  integer,  public, parameter :: wvp = sp      ! for multipole moment.
#else
  integer,  public, parameter :: wvap = dp     ! for multipole moment accumulation array.
  integer,  public, parameter :: wmp = dp      ! for M2L matrix.
  integer,  public, parameter :: wvp = dp      ! for multipole moment.
#endif

  integer,  public, parameter :: wip = kind(0)

! ff parameters
  integer,  public, parameter :: wfp = dp      ! ff parameters in real(8)
! integer,  public, parameter :: wfp = sp      ! ff parameters in real(4)

  INCLUDE 'mpif.h' ! this line is essential to refer to MPI_REAL, MPI_DOUBLE_PRECISION

! ff parameters
  integer,  public, parameter :: MPI_FF_PARAM = MPI_DOUBLE_PRECISION  ! MPI_Bcast in real(8)
! integer,  public, parameter :: MPI_FF_PARAM = MPI_REAL              ! MPI_Bcast in real(4)

! parameter for comm_FMM, must be SAME as wvap
! parameter for comm_FMM, must be SAME as wvp
#if defined (PRECISION_M2L_MIX)
  integer,  public, parameter :: MPI_PREC_MOMVEC_SUPER = MPI_DOUBLE_PRECISION  ! wvap.
  integer,  public, parameter :: MPI_PREC_MOMVEC_COM   = MPI_REAL              ! wvp
#elif defined (PRECISION_M2L_SP)
  integer,  public, parameter :: MPI_PREC_MOMVEC_SUPER = MPI_REAL              ! wvap.
  integer,  public, parameter :: MPI_PREC_MOMVEC_COM   = MPI_REAL              ! wvp
#else /** double precision (default) **/
  integer,  public, parameter :: MPI_PREC_MOMVEC_SUPER = MPI_DOUBLE_PRECISION  ! wvap.
  integer,  public, parameter :: MPI_PREC_MOMVEC_COM   = MPI_DOUBLE_PRECISION  ! wvp
#endif
  integer,  public, parameter :: MPI_PREC_CMPLX_SUPER = MPI_DOUBLE_COMPLEX
  integer,  public, parameter :: MPI_PREC_CMPLX_COM   = MPI_DOUBLE_COMPLEX 

#if defined(PRECISION_PME_SP)
  integer,  public, parameter :: MPI_PREC_PME_COM = MPI_REAL
#else
  integer,  public, parameter :: MPI_PREC_PME_COM = MPI_DOUBLE_PRECISION
#endif

! XYZ
  integer,  public, parameter :: MPI_PREC_FORCE_SUPER  = MPI_DOUBLE_PRECISION
  integer,  public, parameter :: MPI_PREC_COORD_COM    = MPI_DOUBLE_PRECISION

  real(kind=dp ),  parameter :: Cdzero=0.0d0
  real(kind=dp ),  parameter :: Cdone=1.0d0
  real(kind=dp ),  parameter :: Cdtwo=2.0d0
  real(kind=dp ),  parameter :: Cdthree=3.0d0

  real(kind=wrp),  parameter :: Czero=0.0d0
  real(kind=wrp),  parameter :: Cone=1.0d0
  real(kind=wrp),  parameter :: Ctwo=2.0d0
  real(kind=wrp),  parameter :: Cfour=4.0d0
  real(kind=wrp),  parameter :: Csix=6.0d0
  real(kind=wrp),  parameter :: Ctwelve=12.0d0
  real(kind=wrp),  parameter :: Chalf=0.5d0
  real(kind=wrp),  parameter :: Cquarter=0.25d0
! wvp
  real(kind=wvp),  parameter :: Cvzero=0.d0
  real(kind=wvp),  parameter :: Cvone=1.0d0
  real(kind=wvp),  parameter :: Cvhalf=0.5d0
! wvap
  real(kind=wvap), parameter :: Cwvzero=0.d0
  real(kind=wvap), parameter :: Cwve10=1.d10
#ifdef PCFF9
  real(kind=wrp),  parameter :: Cthree=3.0d0
  real(kind=wrp),  parameter :: Cnine=9.0d0
  real(kind=wrp),  parameter :: C18=18.0d0
#endif

end module precision
