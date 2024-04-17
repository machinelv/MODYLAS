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
!! \brief  Program to convert mdff and mdxyz from ascii to binary.
!! \author Kensuke Iwahashi
!<
program modylas_text2bin
  use commandline_args
  use input_mod
  use mpi_tool
  use file_mdff, only : write_mdffbin
  use file_mdxyz_bin, only : write_mdxyzbin
  implicit none

  call mpistart
  call parse_args()
  call parse_input(.false.)
  call write_mdffbin
  !     For debug force field.
  !     call write_memory_mdff
  call write_mdxyzbin
  call mpiend
end program modylas_text2bin
