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
!! \brief  Module and subroutines to read session name.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to read and set session name.
!! \author Kensuke Iwahashi
!<
module session_name_mod
  implicit none
  character(LEN=1024) :: session_name

contains

  subroutine fmod_session_name(string, l)
    implicit none
    integer(4), intent(in) :: l
    character(LEN=l), intent(in) :: string
    session_name = string
    !  There are two variables (one is 'session_name' in Fortran90,
    !  another is 'g_main__session_name' in C) for session name.
    !      call cmod_session_name(session_name, l)
  end subroutine fmod_session_name

end module session_name_mod
