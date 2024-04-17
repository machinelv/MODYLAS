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
!! \brief  Module and subroutines to control .mdtrj file output.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to control .mdtrj file output.
!! \author Kensuke Iwahashi
!<
module file_mdtrj

  integer(4) :: trj_start=0,trj_interval=1
  logical :: trj_output=.false.  ! mdtrj file

contains
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to open file for MD trajectory.
!! \author Kensuke Iwahashi
!<
  subroutine open_trj
    use device_numbers
    use session_name_mod
    use md_condition
    use mpi_tool
    use file_utility, only : create_binary_file
    implicit none

    if(trj_output)then
       call create_binary_file(f_trj, trim(session_name)// '.mdtrj')
    endif
  end subroutine open_trj
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to write header part of MD trajectory.
!! \author Kensuke Iwahashi
!<
  subroutine write_trajectory_header
    use device_numbers
    use md_condition
    implicit none
    if(trj_output)then
    endif
  end subroutine write_trajectory_header
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to record trajectory of current step.
!! \author Kensuke Iwahashi
!<
  subroutine record_current_trajectory
    use trajectory_org
    use unit_cell
    use trajectory_org
    use thermostat
    use md_condition
    use md_condition
    use device_numbers
    implicit none

    if (mod((mdstep-trj_start),trj_interval) == 0) then
       write(f_trj) (mdstep + initial_step), (mdstep + initial_step)*dt
       !  Write coordinates and velocities of atoms.
       write(f_trj) n
       write(f_trj) xyz(1:3,:), v(1:3,:)
       !  Write positions and velocities of thermostats.
       write(f_trj) nnhc
       write(f_trj) rss, vss
       !  Write positions and velocities of barostats.
       write(f_trj) nnhc
       write(f_trj) rssb,vssb
       !  Write cell parameters (length and angles).
       write(f_trj) cellx,celly,cellz,alpha,beta,gamma,vboxg
       call flush(f_trj)
    endif

  end subroutine record_current_trajectory
!----------------------------------------------------------------------
  subroutine fmod_mdtrj_output(ivalue)
    use md_condition
    implicit none
    integer(4), intent(in) :: ivalue
    if (ivalue /= 0) then
       trj_output = .true.
    else
       trj_output = .false.
    endif
  end subroutine fmod_mdtrj_output
!----------------------------------------------------------------------
  subroutine fmod_trj_start(ivalue)
    use md_condition
    implicit none
    integer(4), intent(in) :: ivalue
    trj_start = ivalue
  end subroutine fmod_trj_start
!----------------------------------------------------------------------
  subroutine fmod_trj_interval(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    trj_interval = ivalue
  end subroutine fmod_trj_interval

end module file_mdtrj

