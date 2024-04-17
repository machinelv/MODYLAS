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
!! \brief  Module and subroutines to control .restart.* file I/O.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to control .restart.* file I/O.
!! \author Kensuke Iwahashi
!<
module file_restart
  implicit none
  integer(4), parameter :: f_restart_asc=15,f_restart_bin=16
  logical :: rst_output=.false.  ! restart or restart.bin file
  integer(4) :: restart_start=0,restart_interval=1

contains
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to open file to restart MD.
!! \author Kensuke Iwahashi
!<
  subroutine open_restart
    use device_numbers
    use session_name_mod
    use md_condition
    use mpi_tool
    use file_utility, only : create_binary_file
    implicit none

    if(rst_output)then
       call create_binary_file(f_restart_bin, trim(session_name)// '.restart.bin')
    endif
  end subroutine open_restart
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to record restart file as binary format.
!! \author Kensuke Iwahashi
!<
  subroutine record_restart_binary
    use trajectory_org
    use unit_cell
    use thermostat
    use md_condition
    use device_numbers
    implicit none

    if (mod((mdstep-restart_start),restart_interval) == 0) then
       !  Write coordinates and velocities of atoms.
       rewind(f_restart_bin)
       write(f_restart_bin) n
       write(f_restart_bin) xyz(1:3,:), v(1:3,:)
       !  Write positions and velocities of thermostats.
       write(f_restart_bin) nnhc
       write(f_restart_bin) rss, vss
       !  Write positions and velocities of barostats.
       write(f_restart_bin) nnhc
       write(f_restart_bin) rssb, vssb
       !  Write cell parameters (length and angles).
       write(f_restart_bin) cellx,celly,cellz, alpha,beta,gamma, vboxg
       call flush(f_restart_bin)
    endif

  end subroutine record_restart_binary
!-----------------------------------------------------------------------
  subroutine record_restart_ascii
!>
!! \brief  Subroutine to record restart file as ascii format.
!! \author Kensuke Iwahashi
!<
    use  version
    implicit none

    ! read input files depending on input version
    if(     trim(input_version) == '0.9.0') then
       call  record_restart_ascii_v0_9_0
    else if(trim(input_version) == '1.0.0') then
       call  record_restart_ascii_v1_0_0
    endif
  end subroutine record_restart_ascii
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to record restart file as ascii format for
!!         version 0.9.0.
!! \author Kensuke Iwahashi
!<
  subroutine record_restart_ascii_v0_9_0 !by input version 0.9.0
    use md_condition
    use trajectory_org
    use unit_cell
    use thermostat
    use device_numbers
    use session_name_mod
    use file_utility, only : open_file
    implicit none
    integer(4)::i,j,ii
    real(8), parameter :: m2a=1.0d+10  ! [m]->[A]

    call open_file(f_restart_asc, trim(session_name)//'.restart.asc')

    write(f_restart_asc,'(a10,i10)') '#### step=', (mdstep + initial_step)
    write(f_restart_asc,'(a6)') '<atom>'
    write(f_restart_asc,'(a8,i10)') '  natom=', n

    ! --- position ---
    write(f_restart_asc,'(a13)') '  <positions>'
    do i=1,n
       write(f_restart_asc,9000) (xyz(j,i)*m2a,j=1,3)
    enddo
9000 format(2X,es24.17,1X,es24.17,1X,es24.17)
    write(f_restart_asc,'(a14)') '  </positions>'

    ! --- velocity ---
    write(f_restart_asc,'(a14)') '  <velocities>'
    do i=1,n
       write(f_restart_asc,9100) v(1,i)*m2a, v(2,i)*m2a, v(3,i)*m2a
    enddo
9100 format(2X,es24.17,1X,es24.17,1X,es24.17)
    write(f_restart_asc,'(a15)') '  </velocities>'
    write(f_restart_asc,'(a7)') '</atom>'

    ! --- thermostat ---
    write(f_restart_asc,'(a12)') '<thermostat>'
    write(f_restart_asc,'(a14,i1)') '  nthermostat=',nnhc
    write(f_restart_asc,'(a13)') '  <positions>'
    do i=1,nnhc
       ii=i-1
       write(f_restart_asc,'(a2,E24.17)') '  ',rss(i)
    enddo
    write(f_restart_asc,'(a14)') '  </positions>'
    write(f_restart_asc,'(a14)') '  <velocities>'
    do i=1,nnhc
       ii=i-1
       write(f_restart_asc,'(a2,E24.17)') '  ',vss(i)
    enddo
    write(f_restart_asc,'(a15)') '  </velocities>'
    write(f_restart_asc,'(a13)') '</thermostat>'

    !--- barostat ---
    write(f_restart_asc,'(a10)') '<barostat>'
    write(f_restart_asc,'(a12,i1)') '  nbarostat=',nnhc
    write(f_restart_asc,'(a13)') '  <positions>'
    do i=1,nnhc
       ii=i-1
       write(f_restart_asc,'(a2,E24.17)') '  ',rssb(i)
    enddo
    write(f_restart_asc,'(a14)') '  </positions>'
    write(f_restart_asc,'(a14)') '  <velocities>'
    do i=1,nnhc
       ii=i-1
       write(f_restart_asc,'(a2,E24.17)') '  ',vssb(i)
    enddo
    write(f_restart_asc,'(a15)') '  </velocities>'
    write(f_restart_asc,'(a11)') '</barostat>'

    ! --- cell info. ---
    write(f_restart_asc,'(a10)') '<periodic>'
    write(f_restart_asc,'(a8)') '  <cell>'
    write(f_restart_asc,'(a12)') '    <length>'
    write(f_restart_asc,9200) 'x=',cellx*m2a, 'y=',celly*m2a, 'z=',cellz*m2a
9200 format(4X,a2,es20.10,2X,a2,es20.10,2X,a2,es20.10)
    write(f_restart_asc,'(a13)') '    </length>'
    write(f_restart_asc,'(a11)') '    <angle>'
    write(f_restart_asc,'(a13,es20.10)')'      alpha= ',alpha
    write(f_restart_asc,'(a13,es20.10)')'      beta = ',beta
    write(f_restart_asc,'(a13,es20.10)')'      gamma= ',gamma
    write(f_restart_asc,'(a12)') '    </angle>'
    write(f_restart_asc,'(a11)') '    <vboxg>'
    ! <<< NOTE: vboxg is a symmetric matrix >>>
    write(f_restart_asc,9300) vboxg(1,1), vboxg(1,2), vboxg(1,3)
    write(f_restart_asc,9300) vboxg(1,2), vboxg(2,2), vboxg(2,3)
    write(f_restart_asc,9300) vboxg(1,3), vboxg(2,3), vboxg(3,3)
9300 format(4X,E24.17,1X,E24.17,1X,E24.17)
    write(f_restart_asc,'(a12)') '    </vboxg>'
    write(f_restart_asc,'(a9)') '  </cell>'
    write(f_restart_asc,'(a11)') '</periodic>'

    close(f_restart_asc)
  end subroutine record_restart_ascii_v0_9_0

!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to record restart file as ascii format for
!!         version 1.0.0.
!! \author Kensuke Iwahashi
!<
  subroutine record_restart_ascii_v1_0_0 !by input version 1.0.0(2013.10.24)
    use md_condition
    use trajectory_org
    use unit_cell
    use thermostat
    use device_numbers
    use session_name_mod
    use file_utility, only : open_file
    implicit none
    integer(4)::i,j,ii
    real(8), parameter :: m2a=1.0d+10  ! [m]->[A]

    call open_file(f_restart_asc, trim(session_name)//'.restart.asc')

    write(f_restart_asc,'(a10,i10)') '#### step=', (mdstep + initial_step)
    write(f_restart_asc,'(a6)') '<atom>'
    write(f_restart_asc,'(a8,i10)') '  natom=', n

    ! --- position ---
    write(f_restart_asc,'(a13)') '  <positions>'
    do i=1,n
       write(f_restart_asc,9000) (xyz(j,i)*m2a,j=1,3)
    enddo
9000 format(2X,es24.17,1X,es24.17,1X,es24.17)
    write(f_restart_asc,'(a14)') '  </positions>'

    ! --- velocity ---
    write(f_restart_asc,'(a14)') '  <velocities>'
    do i=1,n
       write(f_restart_asc,9100) v(1,i)*m2a, v(2,i)*m2a, v(3,i)*m2a
    enddo
9100 format(2X,es24.17,1X,es24.17,1X,es24.17)
    write(f_restart_asc,'(a15)') '  </velocities>'
    write(f_restart_asc,'(a7)') '</atom>'

    ! --- thermostat ---
    write(f_restart_asc,'(a12)') '<thermostat>'
    write(f_restart_asc,'(a14,i1)') '  nthermostat=',nnhc
    write(f_restart_asc,'(a13)') '  <positions>'
    do i=1,nnhc
       ii=i-1
       write(f_restart_asc,'(a2,E24.17)') '  ',rss(i)
    enddo
    write(f_restart_asc,'(a14)') '  </positions>'
    write(f_restart_asc,'(a14)') '  <velocities>'
    do i=1,nnhc
       ii=i-1
       write(f_restart_asc,'(a2,E24.17)') '  ',vss(i)
    enddo
    write(f_restart_asc,'(a15)') '  </velocities>'
    write(f_restart_asc,'(a13)') '</thermostat>'

    !--- barostat ---
    write(f_restart_asc,'(a10)') '<barostat>'
    write(f_restart_asc,'(a12,i1)') '  nbarostat=',nnhc
    write(f_restart_asc,'(a13)') '  <positions>'
    do i=1,nnhc
       ii=i-1
       write(f_restart_asc,'(a2,E24.17)') '  ',rssb(i)
    enddo
    write(f_restart_asc,'(a14)') '  </positions>'
    write(f_restart_asc,'(a14)') '  <velocities>'
    do i=1,nnhc
       ii=i-1
       write(f_restart_asc,'(a2,E24.17)') '  ',vssb(i)
    enddo
    write(f_restart_asc,'(a15)') '  </velocities>'
    write(f_restart_asc,'(a11)') '</barostat>'

    ! --- cell info. ---
    write(f_restart_asc,'(a15)') '<periodic cell>'
    write(f_restart_asc,9200) '  <length> '
    write(f_restart_asc,9300) cellx*m2a, celly*m2a, cellz*m2a
    write(f_restart_asc,9200) '  </length>'
    write(f_restart_asc,9200) '  <angle>  '
    write(f_restart_asc,9300) alpha, beta, gamma
    write(f_restart_asc,9200) '  </angle> '
    write(f_restart_asc,9200) '  <vboxg>  '
    ! <<< NOTE: vboxg is a symmetric matrix >>>
    write(f_restart_asc,9300) vboxg(1,1), vboxg(1,2), vboxg(1,3)
    write(f_restart_asc,9300) vboxg(1,2), vboxg(2,2), vboxg(2,3)
    write(f_restart_asc,9300) vboxg(1,3), vboxg(2,3), vboxg(3,3)
    write(f_restart_asc,9200) '  </vboxg> '
9200 format(a11)
9300 format(4X,E24.17,1X,E24.17,1X,E24.17)

    write(f_restart_asc,'(a16)') '</periodic cell>'

    close(f_restart_asc)
  end subroutine record_restart_ascii_v1_0_0
!----------------------------------------------------------------------
  subroutine fmod_restart_start(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    rst_output=.true.
    restart_start = ivalue
  end subroutine fmod_restart_start
!----------------------------------------------------------------------
  subroutine fmod_restart_interval(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    restart_interval = ivalue
  end subroutine fmod_restart_interval

end module file_restart

