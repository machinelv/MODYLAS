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
!! \author Yoshimichi Andoh, Yohei Yamaguchi
!<
module file_dcd_divided

#ifdef DIVIDEDCD

contains
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to open dcd.header file.
!! \author Youhei Yamaguchi
!<
  subroutine open_dcd_info
    use file_dcd
    use device_numbers
    use session_name_mod
    use md_condition
    use mpi_tool
    use file_utility, only : create_binary_file
    implicit none

    if(dcd_output)then
       call create_binary_file(f_dcdxyz, trim(session_name)// '.dcd.header')
    endif

  end subroutine open_dcd_info
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to open dcd.xyz file. 
!! \author Youhei Yamaguchi
!<
  subroutine open_dcd_xyz
    use file_dcd
    use device_numbers 
    use session_name_mod
    use md_condition
    use mpi_tool
    use file_utility, only : create_binary_file
    implicit none
    character(20) :: i

    if(dcd_output)then
      write(i,'(i0.4)') myrank
      call create_binary_file(f_dcdxyz+myrank+1, trim(session_name)// '.dcd.'//trim(adjustl(i)))
    end if

  end subroutine open_dcd_xyz
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to write header part of dcd file.
!! \author Youhei Yamaguchi
!<
  subroutine write_dcd_info_header
    use file_dcd
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
!   icntrl(20)=27
!#endif
      write(f_dcdxyz) Aname,icntrl
      write(f_dcdxyz) nstr
      write(f_dcdxyz) n
      call flush(f_dcdxyz)
    endif

  end subroutine write_dcd_info_header
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to record info of current step to dcd file.
!! \author Youhei Yamaguchi
!<
  subroutine record_current_trajectory_dcd_info
    use device_numbers
    use unit_cell
    use cell_shape
    use trajectory_org
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

    write(f_dcdxyz) cellstr
    call flush(f_dcdxyz)

  end subroutine record_current_trajectory_dcd_info
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to record trajectory of current step to dcd file.
!! \author Youhei Yamaguchi
!<
  subroutine record_current_trajectory_dcd_xyz
    use device_numbers
    use trajectory_org
    use trajectory_mpi
    use subcell
    implicit none
    integer(4) :: i, k, n_selatom
    real(8)::H(3,3)
    real(8),dimension(6)::cellstr
    real(4),dimension(n)::flpx,flpy,flpz
    integer(4),dimension(n) :: atom_number
    integer(4) :: i0, k0

    k = 0
    do k0=1,nselfseg
    do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
      k = k + 1
      i = m2i(i0)
      atom_number(k) = i
      flpx(k)=real(wkxyz(1,i0))
      flpy(k)=real(wkxyz(2,i0))
      flpz(k)=real(wkxyz(3,i0))
    end do ! i0
    end do ! k0
    n_selatom = k

    write(f_dcdxyz+myrank+1)  n_selatom
    write(f_dcdxyz+myrank+1) (atom_number(i),i=1,n_selatom)
    write(f_dcdxyz+myrank+1) (flpx(i)*1e+10,i=1,n_selatom)
    write(f_dcdxyz+myrank+1) (flpy(i)*1e+10,i=1,n_selatom)
    write(f_dcdxyz+myrank+1) (flpz(i)*1e+10,i=1,n_selatom)
    call flush(f_dcdxyz+myrank+1)

  end subroutine record_current_trajectory_dcd_xyz
!-----------------------------------------------------------------------    
#endif
end module file_dcd_divided
