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
!! \brief  Module and subroutines to execute post processes after 
!!         reading inputs.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to execute post processes after reading inputs.
!! \author Kensuke Iwahashi
!<
module commandline_args
  use mpi_tool
  implicit none

contains

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to read command line arguments.
!! \author Kensuke Iwahashi
!<
  subroutine parse_args
    use session_name_mod, only : fmod_session_name
    use version
    implicit none
    character(LEN=1024) :: arg, session
    integer(4) :: i, narg, n, nlen, iargc
    include 'mpif.h'
    integer(4) :: ierr
    if(myrank.eq.mpiout) then
       narg = iargc()
       n = 0
       do i = 1, narg
          call getarg(i,arg)
          if (arg .eq. '-h') then
             call show_help_exit
          else if (arg .eq. '-V') then
             call show_version_exit
          else if (arg(1:1) .eq. '-') then
             call show_usage_error_abort
          else
             n = 1
             session = arg
             exit
          endif
       enddo
    endif
    call mpi_bcast(n,1,mpi_integer,mpiout,mpi_comm_world,ierr)
    call mpi_bcast(session,1024,mpi_character,mpiout, &
         & mpi_comm_world,ierr)
    if (n /= 1)  call show_usage_error_abort
    nlen = len_trim(session)
    call fmod_session_name(session, nlen)
  end subroutine parse_args
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to show usage of MODYLAS.
!! \author Kensuke Iwahashi
!<
  subroutine show_usage_error_abort
    implicit none
    if (myrank == mpiout) then
       write(0,'(a)') 'Usage: modylas [options] session-name'
       write(0,'(a)') 'for more information, type "modylas -h"'
    endif

    call modylas_abort
  end subroutine show_usage_error_abort
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to show brief usage of MODYLAS.
!! \author Kensuke Iwahashi
!<
  subroutine show_help_exit
    use mpi_tool
    implicit none
    if (myrank == mpiout) then
       write(6,'(a)') 'Usage:    modylas [options] session-name'
       write(6,'(a)') '    session-name'
       write(6,'(a)') '      session name of MD simulation'
       write(6,*)
       write(6,'(a)') '    options'
       write(6,'(a)') '      -h'
       write(6,'(a)') '          show this help message and exit'
       write(6,'(a)') '      -V'
       write(6,'(a)') '          show version number of program and exit'
    endif

    call mpiend
  end subroutine show_help_exit

end module commandline_args
