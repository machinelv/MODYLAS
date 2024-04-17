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
!! \brief  Module and subroutines to control .dcd file output.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to control .dcd file output.
!! \author Yoshimichi Andoh
!<
module file_dcd
  implicit none
  logical :: dcd_output=.false.  ! dcd file
  integer(4) :: dcd_start=0,dcd_interval=1

contains
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to open .dcd file.
!! \author Yoshimichi Andoh
!<
  subroutine open_dcd
    use device_numbers
    use session_name_mod
    use md_condition
    use mpi_tool
    use file_utility, only : create_binary_file
    implicit none
    if(dcd_output)then
       call create_binary_file(f_dcd, trim(session_name)// '.dcd')
    endif
  end subroutine open_dcd
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to write header part of .dcd file.
!! \author Yoshimichi Andoh
!<
  subroutine write_dcd_header
    use device_numbers
    use md_condition
    use trajectory_org
    implicit none
    character(len=4)::Aname
    integer(4) :: i,nstr
    integer(4),dimension(20)::icntrl

    if(dcd_output)then
       Aname='CORD'
       icntrl=0
       nstr=0

       icntrl(1)=md_condition__howmany_steps/dcd_interval
       icntrl(2)=1
       icntrl(3)=1
       icntrl(4)=md_condition__howmany_steps
       icntrl(8)=n*3
       icntrl(10)=981668463
       icntrl(11)=1
!#if defined(CHARMM36) || defined(CHARMMFSW)
       icntrl(20)=36
!#else
!     icntrl(20)=27
!#endif

       write(f_dcd) Aname,icntrl
       write(f_dcd) nstr
       write(f_dcd) n

    endif
  end subroutine write_dcd_header
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to record trajectory of current step to .dcd file.
!! \author Yoshimichi Andoh
!<
  subroutine record_current_trajectory_dcd
    use trajectory_org
    use unit_cell
    use cell_shape
    use mpi_tool
    use md_condition
    use device_numbers
    implicit none
    integer(4)::i
    real(8)::H(3,3)
    real(8),dimension(6)::cellstr
    real(4),dimension(n)::flpx,flpy,flpz

    call cell_convert1(cellx,celly,cellz,alpha,beta,gamma,H)
    cellstr(1)=cellx*1d+10 !A
    cellstr(2)=gamma       !degree
    cellstr(3)=celly*1d+10 !A
    cellstr(4)=beta        !degree
    cellstr(5)=alpha       !degree
    cellstr(6)=cellz*1d+10 !A

    do i=1,n
       flpx(i)=real(xyz(1,i))
       flpy(i)=real(xyz(2,i))
       flpz(i)=real(xyz(3,i))
    enddo

    write(f_dcd) cellstr
    write(f_dcd) (flpx(i)*1e+10,i=1,n)
    write(f_dcd) (flpy(i)*1e+10,i=1,n)
    write(f_dcd) (flpz(i)*1e+10,i=1,n)
    call flush(f_dcd)

  end subroutine record_current_trajectory_dcd
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to set a flag for .dcd file output.
!! \author Yoshimichi Andoh
!<
  subroutine fmod_dcd_output(ivalue)
    use md_condition
    implicit none
    integer(4), intent(in) :: ivalue
    if (ivalue /= 0) then
       dcd_output = .true.
    else
       dcd_output = .false.
    endif
  end subroutine fmod_dcd_output
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to set start step for .dcd file output.
!! \author Yoshimichi Andoh
!<
  subroutine fmod_dcd_start(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    dcd_start = ivalue
  end subroutine fmod_dcd_start
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to set interval step for .dcd file output.
!! \author Yoshimichi Andoh
!<
  subroutine fmod_dcd_interval(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    dcd_interval = ivalue
  end subroutine fmod_dcd_interval

end module file_dcd

