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
! Dummy module for file_xtc.F90 because they depend on exernal library
! without -DXTC, this module is compiled istead for successful comilation
!>
!! \file
!! \brief  Dummy module for file_xtc.F90.
!<
!----------------------------------------------------------------------
!>
!! \brief  Dummy module for file_xtc.F90.
!! \author Ryo Urano
!<
module xdr
    use, intrinsic :: iso_c_binding, only: C_PTR, C_CHAR, C_FLOAT, C_DOUBLE, C_INT
    type, public :: xtcfile
          real(C_FLOAT) :: dummy
    end type xtcfile
  end module xdr
