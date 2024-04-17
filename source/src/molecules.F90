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
!! \brief  Module and subroutines to store molecule information.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to store molecule information.
!! \author Kazushi FUJIMOTO
!<
module molecules
  implicit none
  integer(4),allocatable :: atom2mol(:),mol2atom(:),moltop(:)
  integer(4) :: md_periodic__howmany_molecules=-1
  integer(4),allocatable :: molecules__howmany_atoms(:)

contains

  subroutine fmod_alloc_atom2mol(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(atom2mol(ivalue))
  end subroutine fmod_alloc_atom2mol
!----------------------------------------------------------------------
  subroutine fmod_alloc_mol2atom(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(mol2atom(ivalue))
  end subroutine fmod_alloc_mol2atom
!----------------------------------------------------------------------
  subroutine fmod_alloc_moltop(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(moltop(ivalue))
  end subroutine fmod_alloc_moltop

  subroutine fmod_set_mol2atom(imol0,j0,ivalue0,line_number)
    use trajectory_org
    use md_periodic
    use segments, only : nsegments, irearrange, seg_ready
    implicit none
    integer(4), intent(in) :: imol0, j0, ivalue0, line_number
    integer(4) :: imol, j, ivalue
    integer(4),save :: ic=0

    imol = imol0 + 1
    ivalue = ivalue0 + 1
    j = j0 + 1
    if (j > molecules__howmany_atoms(imol)) then
       write(0,'(a)')  'ERROR : ', &
            &    'The number of md_periodic.molecules.[1].atoms[] is ' // &
            &    'out of bounds.  '// &
            &    'It must be less than ', molecules__howmany_atoms(imol), '.'
       call modylas_abort()
    endif
    if(j == 1) then
       moltop(imol) = ivalue0
    endif
    ic = ivalue
    if(nsegments < 0) then
       mol2atom(ic) = ivalue0
    else
       mol2atom(ic) = irearrange(ivalue)
    endif

    atom2mol(irearrange(ivalue)) = imol

    if (nsegments < 0 .and. ic == n) then
       seg_ready = .TRUE.
    endif
  end subroutine fmod_set_mol2atom
!----------------------------------------------------------------------
  subroutine fmod_howmany_molecules(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    md_periodic__howmany_molecules = ivalue
  end subroutine fmod_howmany_molecules
!----------------------------------------------------------------------
  subroutine fmod_alloc_mol_natoms(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(molecules__howmany_atoms(ivalue))
    molecules__howmany_atoms = -1
  end subroutine fmod_alloc_mol_natoms
!----------------------------------------------------------------------
  subroutine fmod_set_mol_natoms(imol0, ivalue, line_number)
    use mpi_tool
    implicit none
    integer(4), intent(in) :: imol0, ivalue, line_number
    integer(4) :: imol
    imol = imol0 + 1
    if (imol > md_periodic__howmany_molecules) then
       write(0,'(a,i0)')  'ERROR (line ', line_number, ') : '
       write(0,'(a,i0,a)')  &
            &    'The number of md_periodic.molecules[].howmany_atoms ' // &
            &    'is out of bounds.  '// &
            &    'It must be less than ', md_periodic__howmany_molecules, '.'
       call modylas_abort()
    endif
    molecules__howmany_atoms(imol) = ivalue
  end subroutine fmod_set_mol_natoms

end module molecules
