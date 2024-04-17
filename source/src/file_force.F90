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
!! \brief  Module and subroutines to control .force file output.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to control .force file output.
!! \author Yoshimichi Andoh
!<
module file_force
  implicit none
  integer(4) :: force_start=0,force_interval=1

contains

!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to open file for mean force.
!! \author Yoshimichi Andoh
!<
  subroutine open_force
    use device_numbers
    use session_name_mod
    use file_utility
    implicit none

    call create_file(f_force, trim(session_name)// '.force')
  end subroutine open_force
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to write header part of calculated mean force.
!! \author Yoshimichi Andoh
!<
  subroutine write_force_header
    use device_numbers
    implicit none

    write(f_force,'(a10,4a20,a12)') &
         &  "#     step",  &
         &  "      mean force A-B", &
         &  "       raw force A-B", &
         &  "         raw force A", &
         &  "         raw force B", &
         &  "     dist_AB"
    call flush(f_force)
  end subroutine write_force_header
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to record mean force
!! \author Yoshimichi Andoh
!<
  subroutine record_current_force
    use device_numbers
    use md_condition
    use shake_rattle_roll
    implicit none
    real(8)::aMF_AB

    if (mod((mdstep-force_start),force_interval) == 0) then
       aMF_AB=sMF_AB/dble(mdstep)
       write(f_force,'(i10,4e20.12,f12.5)')  (mdstep + initial_step), aMF_AB, MF_AB, MForceA, MForceB, drCOM*1d+10
       call flush(f_force)
    endif
  end subroutine record_current_force

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to set start step for .force file output.
!! \author Yoshimichi Andoh
!<
  subroutine fmod_force_start(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    force_start = ivalue
  end subroutine fmod_force_start

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to set interval step for .force file output.
!! \author Yoshimichi Andoh
!<
  subroutine fmod_force_interval(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    force_interval = ivalue
  end subroutine fmod_force_interval

end module file_force
