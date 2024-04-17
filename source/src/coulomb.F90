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
!! \brief  Module to store partial charge information.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to store partial charge infomation.
!! \author Kensuke Iwahashi, Kazushi Fujimoto
!<
module coulomb_mod
  use omp_lib
  implicit none
  real(8),allocatable :: chgv(:)
  real(8),allocatable :: chg_list(:,:), chgfactor_list(:)

contains

!>
!! \brief  Subroutine to allocate charge arrays.
!! \author Kensuke Iwahashi, Kazushi Fujimoto
!<
  subroutine fmod_alloc_chg_list(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    integer(4) :: nomp=1
!$  nomp = omp_get_max_threads()
    allocate(chg_list(ivalue,0:nomp-1))
  end subroutine fmod_alloc_chg_list
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to allocate charge arrays.
!! \author Kensuke Iwahashi, Kazushi Fujimoto
!<
  subroutine fmod_alloc_chgfactor_list(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(chgfactor_list(ivalue))
  end subroutine fmod_alloc_chgfactor_list
!---------------------------------------------------------------------
!>
!! \brief  Subroutine to set chgv.
!! \author Kensuke Iwahashi, Kazushi Fujimoto
!<
  subroutine fmod_set_mol_chgv(value)
    implicit none
    real(8), intent(in) :: value
    integer(4),save :: i=0

    i = i + 1
    chgv(i) = value

  end subroutine fmod_set_mol_chgv
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to read chgv.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_chgv
    use device_numbers, only : f_mdff
    use param
    use mpi_tool
    implicit none
    include 'mpif.h'
    integer(4) :: ierr

    if(myrank==0) read(f_mdff) chgv
    call MPI_Bcast(chgv, npara, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  end subroutine read_chgv
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write chgv in mdff.bin File.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_chgv
    use device_numbers, only : f_mdff
    implicit none
    write(f_mdff) chgv
  end subroutine write_chgv
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write chgv.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_chgv
    implicit none
    write(*,*) '[write_chgv]'
    write(*,*) chgv
  end subroutine write_memory_chgv

end module coulomb_mod
