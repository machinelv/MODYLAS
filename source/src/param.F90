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
!! \brief  Module and subroutines to set parameter number.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to set parameter number.
!! \author Kazushi FUJIMOTO
!<
module param
  implicit none
  integer(4) :: npara
  integer(4),allocatable :: paranum(:)

contains

!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write parameter number mdff.bin
!! \author Kazushi FUJIMOTO
!<
  subroutine write_parameter_number
    use device_numbers, only : f_mdff
    implicit none
    write(f_mdff) npara
    write(f_mdff) paranum(:)
  end subroutine write_parameter_number
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write parameter number.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_parameter_number
    implicit none
    write(*,*) '[write_parameter_number]'
    write(*,*) npara
    write(*,*) paranum(:)
  end subroutine write_memory_parameter_number

  subroutine fmod_set_param_num(np,nmol,nmolmol,natommol)
    use mpi_tool
    use trajectory_org
    use molecules, only : fmod_alloc_moltop, fmod_howmany_molecules, fmod_alloc_mol_natoms
    implicit none
    integer(4), intent(in) :: np, nmol, nmolmol(nmol), natommol(nmol)
    integer(4) :: i, j, k, atomnum, molnum
    integer(4) :: natomin, howmany_mols

    allocate(paranum(n))
    molnum=0
    atomnum=0
    howmany_mols=0
    do i = 1, nmol
       if(i==1) then
          molnum=1
       else
          molnum = molnum + natommol(i-1)
       endif
       howmany_mols=howmany_mols+nmolmol(i)
       do j = 1, nmolmol(i)
          do k = 1, natommol(i)
             atomnum = atomnum + 1
             paranum(atomnum) = molnum + (k-1)
          enddo
       enddo
    enddo

    if(n.ne.atomnum)then
       write(0,*) 'ERROR: natomin (by pdb) .ne. atomnum (by mdff)', &
            &             natomin,atomnum
       call modylas_abort
    endif
    nmols=howmany_mols

    call fmod_howmany_molecules(howmany_mols)
    call fmod_alloc_mol_natoms(howmany_mols)
    call fmod_alloc_moltop(howmany_mols)

  end subroutine fmod_set_param_num

end module param
