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
!! \brief  Module to show program version.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to show program version.
!! \author Kensuke Iwahashi
!<
module version
#include "../config.h"

  character(20) :: MODYLAS_version="1.1.0"
  character(10) :: input_version

contains

!>
!! \brief  Subroutine to show MODYLAS version.
!! \author Kensuke Iwahashi
!<
  subroutine show_version_exit
    use mpi_tool
    implicit none
    character(LEN=*),parameter :: str_package=PACKAGE_NAME
    character(LEN=*),parameter :: str_version=PACKAGE_VERSION
    if (myrank == mpiout) then
       write(6,'(a,a,a,a)') str_package, &
            &    ' (MOlecular DYnamics simulation software for LArge Systems)', &
            &    ' --- version ', str_version
    endif
    call modylas_abort
  end subroutine show_version_exit

end module version
