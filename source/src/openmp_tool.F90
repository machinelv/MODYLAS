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
!! \brief  Module and subroutine to initialize/finalize OpenMP
!<
!----------------------------------------------------------------------
!>
!! \brief  Module initialize/finalize OpenMP
!! \author Tatsuya Sakashita
!>
module openmp_tool
  use omp_lib
  implicit none
  integer(4),protected :: nomp=1

contains

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to get number of OMP threads.
!! \author Tatsuya Sakashita
!<
  subroutine initialize_openmp
!----------------------------------------------------------------------
    implicit none
!$  nomp = omp_get_max_threads()
  end subroutine initialize_openmp

end module openmp_tool
