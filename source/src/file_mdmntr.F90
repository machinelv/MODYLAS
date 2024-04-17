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
!! \brief  Module and subroutines to control .mdmntr file output.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to control .mdmntr file output.
!! \author Kensuke Iwahashi
!<
module file_mdmntr
  implicit none
  integer(4) :: mntr_start=0,mntr_interval=1

contains

!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to open file for MD monitor.
!! \author Kensuke Iwahashi
!<
  subroutine open_mntr
    use device_numbers
    use session_name_mod
    use mpi_tool
    use file_utility, only : create_file
    implicit none

    call create_file(f_mntr, trim(session_name)// '.mdmntr')
  end subroutine open_mntr
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to write header part of MD monitor.
!! \author Kensuke Iwahashi
!<
  subroutine write_monitors_header
    use device_numbers
    use session_name_mod
    use md_periodic
    implicit none
    write(f_mntr,'(a,a,a)') &
         & '## ', trim(session_name), '.mdmntr' // &
         & ' -- monitor variables output from MD calculation by modylas'
    write(f_mntr,'(a)') '#'
    write(f_mntr,'(a)') '# datas below are formated as:'

if (    md_periodic__type == FMM) then
    write(f_mntr,'(13a)') &
         &  '#',' step     '          , &
         &      ' time               ', &
         &      ' Hamiltonian        ', &
         &      ' potential-E        ', &
         &      ' kinetic-E          ', &
         &      ' total energy       ', &
         &      ' temperature        ', &
         &      ' volume             ', &
         &      ' pressure           ', &
         &      ' box-length(x)      ', &
         &      ' box-length(y)      ', &
         &      ' box-length(z)      '
else
    write(f_mntr,'(22a)') &
         &  '#',' step     '          , &
         &      ' time               ', &
         &      ' Hamiltonian        ', &
         &      ' potential-E        ', &
         &      ' kinetic-E          ', &
         &      ' total energy       ', &
         &      ' temperature        ', &
         &      ' volume             ', &
         &      ' pressure           ', &
         &      ' box-length(x)      ', &
         &      ' box-length(y)      ', &
         &      ' box-length(z)      ', &
         &      ' alpha              ', &
         &      ' beta               ', &
         &      ' gamma              ', &
         &      ' Pxx                ', &
         &      ' Pyy                ', &
         &      ' Pzz                ', &
         &      ' Pxy                ', &
         &      ' Pxz                ', &
         &      ' Pyz                '
endif

if (    md_periodic__type == FMM) then
    write(f_mntr,'(13a)') &
         &  '#', '          '          , &
         &       '              [sec] ', &
         &       '           [J/cell] ', &
         &       '           [J/cell] ', &
         &       '           [J/cell] ', &
         &       '           [J/cell] ', &
         &       '                [K] ', &
         &       '               [m3] ', &
         &       '               [Pa] ', &
         &       '                [m] ', &
         &       '                [m] ', &
         &       '                [m] '
else
    write(f_mntr,'(22a)') &
         &  '#', '          '          , &
         &       '              [sec] ', &
         &       '           [J/cell] ', &
         &       '           [J/cell] ', &
         &       '           [J/cell] ', &
         &       '           [J/cell] ', &
         &       '                [K] ', &
         &       '               [m3] ', &
         &       '               [Pa] ', &
         &       '                [m] ', &
         &       '                [m] ', &
         &       '                [m] ', &
         &       '           [degree] ', &
         &       '           [degree] ', &
         &       '           [degree] ', &
         &       '               [Pa] ', &
         &       '               [Pa] ', &
         &       '               [Pa] ', &
         &       '               [Pa] ', &
         &       '               [Pa] ', &
         &       '               [Pa] '
endif
    write(f_mntr,'(a)') '#'
    call flush(f_mntr)
  end subroutine write_monitors_header
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to record various outputs of MD such as Hamiltonian.
!! \author Kensuke Iwahashi
!<
  subroutine record_current_monitors
    use unit_cell
    use md_condition
    use device_numbers
    use md_monitors
    use md_periodic
    implicit none

    if (mod((mdstep - mntr_start), mntr_interval) == 0) then
if (    md_periodic__type == FMM) then
       write(f_mntr,'(i10,11es20.12)') &
            &   (mdstep + initial_step), (mdstep + initial_step)*dt, hamiltonian, &
            &   p_energy, k_energy, &
            &   t_energy, temperature, &
            &   cellvol, pressure, cellx, celly, cellz
else
       write(f_mntr,'(i10,21es20.12)') &
            &   (mdstep + initial_step), (mdstep + initial_step)*dt, hamiltonian, &
            &   p_energy, k_energy, &
            &   t_energy, temperature, &
            &   cellvol, pressure, cellx, celly, cellz, &
            &   alpha, beta, gamma, &
            &   ptensor(1),ptensor(2),ptensor(3), &  ! Pxx, Pyy, Pzz
            &   ptensor(4),ptensor(5),ptensor(6)     ! Pxy, Pxz, Pyz
endif
       call flush(f_mntr)
    endif

  end subroutine record_current_monitors
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to set start step for .mdmntr file output.
!! \author Kensuke Iwahashi
!<
  subroutine fmod_mntr_start(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    mntr_start = ivalue
  end subroutine fmod_mntr_start
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to set interval step for .mdmntr file output.
!! \author Kensuke Iwahashi
!<
  subroutine fmod_mntr_interval(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    mntr_interval = ivalue
  end subroutine fmod_mntr_interval

end module file_mdmntr
