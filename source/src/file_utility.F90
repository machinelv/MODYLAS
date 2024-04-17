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
!! \brief  Module and subroutines as utility tools for file I/O.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to control .dcd file output.
!! \author Kensuke Iwahasi, Zhiye Tang
!<
module file_utility
  implicit none
  logical :: ascii_output=.false.
  logical :: backup_output=.true.

contains

!>
!! \brief  Subroutine for backing up existing file
!! \author Zhiye Tang
!<
subroutine backup_existing_file(file_name)
  implicit none
  character(LEN=*), intent(in) :: file_name
  character(LEN=256) :: new_file_name
  character(LEN=16) :: new_file_index
  logical :: file_exists
  integer(4) :: suffix

  inquire(file=trim(file_name), exist=file_exists)
  if( file_exists ) then
      suffix = 1
      do while( file_exists )
          write(new_file_index,"(I10.1)") suffix
          new_file_index = adjustl(new_file_index)
          new_file_name = "#"//trim(file_name)//"_bk."//trim(new_file_index)//"#"
          inquire(file=trim(new_file_name), exist=file_exists)
          suffix = suffix + 1
      enddo
      write(*,*)"Backup "//trim(file_name)//" to "//trim(new_file_name)
      call rename(trim(file_name), trim(new_file_name))
  endif
end subroutine backup_existing_file

!>
!! \brief  Subroutine to create an output file.
!! \author Kensuke Iwahasi
!<
  subroutine create_file(f_num, file_name)
    use mpi_tool, only : modylas_abort
    integer(4), intent(in) :: f_num
    character(LEN=*), intent(in) :: file_name
    integer(4) :: io

    if(backup_output)then
      call backup_existing_file(trim(file_name))
    endif

    open(f_num, file=file_name, iostat=io, status='replace', access='sequential',form='formatted')
    if (io /= 0) then
       write(0,*) 'ERROR: Cannot create ' // trim(file_name)
       call modylas_abort()
    endif
  end subroutine create_file

!>
!! \brief  Subroutine to create an output file in binary format.
!! \author Kensuke Iwahasi
!<
  subroutine create_binary_file(f_num, file_name)
    use mpi_tool, only : modylas_abort
    integer(4), intent(in) :: f_num
    character(LEN=*), intent(in) :: file_name
    integer(4) :: io

    if(backup_output)then
      call backup_existing_file(trim(file_name))
    endif

    open(f_num, file=file_name, iostat=io, status='replace', access='sequential',form='unformatted')
    if (io /= 0) then
       write(0,*) 'ERROR: Cannot create ' // trim(file_name)
       call modylas_abort()
    endif
  end subroutine create_binary_file

!>
!! \brief  Subroutine for creating new binary file while backing up existing file
!! \author Zhiye Tang
!<
  subroutine create_new_binary_file(f_num, file_name)
    use mpi_tool, only : modylas_abort
    integer(4), intent(in) :: f_num
    character(LEN=*), intent(in) :: file_name
    integer(4) :: io
    logical :: file_exists

    if(backup_output)then
      call backup_existing_file(trim(file_name))
    endif

    open(f_num, file=file_name, iostat=io, status='new', access='sequential',form='unformatted')
    if (io /= 0) then
       write(0,*) 'ERROR: Cannot create ' // trim(file_name)
       inquire(file=trim(file_name), exist=file_exists)
       if (file_exists) then
          write(0,*) 'This file already exists.'
       endif
       call modylas_abort()
    endif
  end subroutine create_new_binary_file

!>
!! \brief  Subroutine to open an input file.
!! \author Kensuke Iwahasi
!<
  subroutine open_file(f_num, file_name)
    use mpi_tool, only : modylas_abort
    integer(4), intent(in) :: f_num
    character(LEN=*), intent(in) :: file_name
    integer(4) :: io

    open(f_num, file=file_name, iostat=io, access='sequential',form='formatted')
    if (io /= 0) then
       write(0,*) 'ERROR: Cannot open ' // trim(file_name)
       call modylas_abort()
    endif
  end subroutine open_file

!>
!! \brief  Subroutine to open an input file in binary format.
!! \author Kensuke Iwahasi
!<
  subroutine open_binary_file(f_num, file_name)
    use mpi_tool, only : modylas_abort
    integer(4), intent(in) :: f_num
    character(LEN=*), intent(in) :: file_name
    integer(4) :: io

    open(f_num, file=file_name, iostat=io, status='old', access='sequential',form='unformatted')
    if (io /= 0) then
       write(0,*) 'ERROR: Cannot open ' // trim(file_name)
       call modylas_abort()
    endif
  end subroutine open_binary_file

!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to check whether file is binary.
!! \author Kensuke Iwahashi
!<
  logical function can_read_binary_file(fname)
    implicit none
    character(len=*),intent(in) :: fname
    integer(4) :: io
    open(99, file=fname, &
         &         iostat=io, status='old', &
         &         access='sequential',form='unformatted')
    if (io == 0) then
       can_read_binary_file = .true.
       close(99)
    else
       can_read_binary_file = .false.
    endif
  end function can_read_binary_file
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to set a flag to output .restart.asc
!! \author Kensuke Iwahasi
!<
  subroutine fmod_ascii_output(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    if (ivalue /= 0) then
       ascii_output = .true.
    else
       ascii_output = .false.
    endif
  end subroutine fmod_ascii_output
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to set a flag to output #***_bk# files
!! \author Zhiye Tang
!<
  subroutine fmod_backup_output(ivalue)
    use md_condition
    implicit none
    integer(4), intent(in) :: ivalue
    if (ivalue /= 0) then
       backup_output = .true.
    else
       backup_output = .false.
    endif
  end subroutine fmod_backup_output
!----------------------------------------------------------------------

end module file_utility
