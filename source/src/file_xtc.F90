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
!! \brief  Module and subroutines to control .xtc file output.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to control .xtc file output.
!! \author Ryo Urano
!<
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This subroutine require external library for xtc file formats because native their file is written
!! by C/C++. Due to the difference of binary style, this requires such library.
!! To install the library, see source/external/README
!! you need path and -l option for Makefile or in cmake command 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module file_xtc
!  1. Use the xdr interface
   use xdr, only: xtcfile
  implicit none
  logical :: xtc_output=.false.  ! xtc file
  integer(4) :: xtc_start=0,xtc_interval=1000000
 ! 2. Declare a variable of type xtcfile
  type(xtcfile) :: xtcf
  type(xtcfile) :: xtc_out

#ifdef XTC
contains
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to open xtc file.
!! \author Ryo Urano 
!<
  subroutine open_xtc
    use device_numbers
    use session_name_mod
    use md_condition
    use mpi_tool
    use file_utility, only : create_binary_file
    implicit none
    if(xtc_output)then
       call xtc_out % init(trim(session_name)// '.xtc', 'w')
    endif
  end subroutine open_xtc

  subroutine fmod_xtc_output(ivalue)
    use md_condition
    implicit none
    integer(4), intent(in) :: ivalue
    if (ivalue /= 0) then
       xtc_output = .true.
    else
       xtc_output = .false.
    endif
  end subroutine fmod_xtc_output
!----------------------------------------------------------------------
  subroutine fmod_xtc_start(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    xtc_start = ivalue
  end subroutine fmod_xtc_start
!----------------------------------------------------------------------
  subroutine fmod_xtc_interval(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    xtc_interval = ivalue
  end subroutine fmod_xtc_interval
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to record trajectory of current step to xtc file.
!! \author Ryo Urano
!<
  subroutine record_current_trajectory_xtc
    use trajectory_org
    use unit_cell
    use cell_shape
    use mpi_tool
    use md_condition
    use device_numbers
    implicit none
    real(8)::H(3,3)

    call cell_convert1(cellx,celly,cellz,alpha,beta,gamma,    H)
    
    xtcf % time= (mdstep + initial_step)*dt *1e+12
    xtcf % STEP= (mdstep + initial_step)-xtc_start
    xtcf % box(1,1)=cellx*1e+9
    xtcf % box(2,2)=celly*1e+9
    xtcf % box(3,3)=cellz*1e+9
    xtcf % prec=1000.000000

    call xtc_out % write(n, xtcf % STEP, xtcf % time, xtcf % box, real(xyz)*1e+9, xtcf % prec)

  end subroutine record_current_trajectory_xtc

#endif

!----------------------------------------------------------------------
end module file_xtc
