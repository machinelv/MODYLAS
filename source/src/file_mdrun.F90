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
!! \brief  Module and subroutines to control .mdrun file output.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to control .mdrun file output.
!! \author Kensuke Iwahashi
!<
module file_mdrun

contains
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to open file to write MD step, time and speed.
!! \author Kensuke Iwahashi
!<
  subroutine open_run
    use device_numbers
    use session_name_mod
    use mpi_tool
    use file_utility, only : create_file
    implicit none

    call create_file(f_run, trim(session_name)// '.mdrun')
  end subroutine open_run
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to show the MD step, time and speed.
!! \author Kensuke Iwahashi
!<
  subroutine record_current_runtime_info
    use md_condition
    use device_numbers
    use session_name_mod
    use mpi_tool
! TZTBEGIN
    use md_multiplestep, only : maxMTm, maxMTl
! TZYEND
    implicit none
    real(4) :: t,ta(2)
    real(4) :: etime
    real(8) :: t0
    real(8),save :: tsave
! TZYBEGIN
    real(8) :: buffer
    character(len=8) :: unit
! TZYEND
    include 'mpif.h'
    if (mdstep == 0) then
       tsave = MPI_WTIME()
    else
       t     = etime(ta)
       t0    = MPI_WTIME()
       rewind(f_run)
       write(f_run,'(3a)') '## ', trim(session_name)// '.mdrun', &
            &          ' -- run time information of MD calculation by modylas'
       write(f_run,'(a)') '#'
       write(f_run,'(a,i10)') 'step:      ', (mdstep + initial_step)
       write(f_run,'(a,f15.6,a)') 'CPU time:  ', t , ' [sec]'
       write(f_run,'(a,f15.6,a)') 'for MD:    ', t0-tsave, ' [sec]'
       write(f_run,'(a,f15.6,a)') 'time/step: ', &
            &                         (t0-tsave)/dble(mdstep),' [sec/step]'
       write(f_run,'(a)') '#'
! TZYBEGIN
!      if (mdstep == md_condition__howmany_steps) then
           ! time for 1ns
           buffer = 1.0d-9 / dt * (t0-tsave)/dble(mdstep)
           unit = ' [sec]'
           if( buffer > 60.0d0 ) then
               buffer = buffer / 60.0d0
               unit = ' [min]'
               if( buffer > 60.0d0 ) then
                   buffer = buffer / 60.0d0
                   unit = ' [hour]'
                   if( buffer > 24.0d0 ) then
                       buffer = buffer / 24.0d0
                       unit = ' [day]'
                   endif
               endif
           endif
           write(f_run,'(a,f15.6,a)') 'time/ns:   ', buffer, unit
           ! ns per day
           buffer = 86400.0d0 / ((t0-tsave)/dble(mdstep)) * dt * 1.0d9
           write(f_run,'(a,f15.6,a)') 'ns/day:    ', buffer, ' [ns/day]'
           ! md length
!          buffer = dt * md_condition__howmany_steps * 1.0d15
           buffer = dt * mdstep * 1.0d15
           unit = ' [fs]'
           if( buffer > 1000.0d0 ) then
               buffer = buffer * 0.001d0
               unit = ' [ ps]'
               if( buffer > 1000.0d0 ) then
                   buffer = buffer * 0.001d0
                   unit = ' [ ns]'
                   if( buffer > 1000.0d0 ) then
                       buffer = buffer * 0.001d0
                       unit = ' [ us]'
                       if( buffer > 1000.0d0 ) then
                           buffer = buffer * 0.001d0
                           unit = ' [ ms]'
                           if( buffer > 1000.0d0 ) then
                               buffer = buffer * 0.001d0
                               unit = ' [ s]'
                           endif
                       endif
                   endif
               endif
           endif
           write(f_run,'(a,f15.6,a)') 'md length: ', buffer, unit
!          write(f_run,'(a)')         'calculation times per step:'
!          write(f_run,'(a,i15,a)')   '    short: ', maxMTm * maxMTl, ' [times]'
!          write(f_run,'(a,i15,a)')   '    middle:', maxMTl, ' [times]'
!          write(f_run,'(a,i15,a)')   '    long:  ', 1, ' [times]'
!      else
!          write(f_run,'(a)') 'time/ns and ns/day will be provided at the final step'
!      endif
! TZYEND
       call flush(f_run)
    endif
  end subroutine record_current_runtime_info

end module file_mdrun
