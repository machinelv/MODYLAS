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
module mol_info
  implicit none
  integer(4),allocatable :: molinfo_nmoles(:), molinfo_natoms(:)
  integer(4),allocatable :: molinfo_start(:), molinfo_end(:)
  integer(4),allocatable :: para_start(:), para_end(:)
  integer(4) :: nID, id, nmol
  integer(4) :: nmseg, segstart=0
  integer(4) :: nmbond, bondstart=0
  integer(4) :: nmangle, anglestart=0
  integer(4) :: nmub, ubstart=0
  integer(4) :: nmdihd, dihdstart=0
  integer(4) :: nmitor, itorstart=0

contains

!----------------------------------------------------------------------
!K. FUJIMOTO 2013.08.23
!--------------------------------------------------------
  subroutine fmod_alloc_molinfo(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    nID=ivalue
    allocate(molinfo_nmoles(ivalue))
    allocate(molinfo_natoms(ivalue))
    allocate(molinfo_start(ivalue))
    allocate(molinfo_end(ivalue))
    allocate(para_start(ivalue))
    allocate(para_end(ivalue))
  end subroutine fmod_alloc_molinfo
!----------------------------------------------------------------------
  subroutine fmod_set_molinfo(ivalue0, ivalue2, ivalue3)
    implicit none
    integer(4), intent(in) :: ivalue0, ivalue2, ivalue3
    integer(4) :: ivalue
    ivalue = ivalue0 + 1
    molinfo_nmoles(ivalue) = ivalue2
    molinfo_natoms(ivalue) = ivalue3
    if(ivalue==1) then
       molinfo_start(ivalue) = 1
       molinfo_end(ivalue) = molinfo_start(ivalue) + ivalue2 * &
            &ivalue3 -1
    else
       molinfo_start(ivalue) = molinfo_end(ivalue-1) + 1
       molinfo_end(ivalue) = molinfo_start(ivalue) + ivalue2 * &
            &ivalue3 -1
    endif

    if(ivalue==1) then
       para_start(ivalue) = 1
       para_end(ivalue) =  ivalue3
    else
       para_start(ivalue) = para_end(ivalue-1) + 1
       para_end(ivalue) = para_start(ivalue) + ivalue3 - 1
    endif
  end subroutine fmod_set_molinfo
!-----------------------------------------------------------------------
  subroutine  fmod_set_mol_seg_num(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    if(id==1) then
       segstart = 0
    else
       segstart = segstart + nmseg  * molinfo_nmoles(id-1)
    endif
    nmseg = ivalue
  end subroutine fmod_set_mol_seg_num
!-----------------------------------------------------------------------
  subroutine fmod_set_mol_seg_natom(ivalue1,ivalue2)
    use segments, only : fmod_set_seg_natoms
    implicit none
    integer(4), intent(in) :: ivalue1, ivalue2
    integer(4) ::  temp1, temp2
    integer(4) :: i

    do i = 1, molinfo_nmoles(id)
       temp1 = segstart + ivalue1 + nmseg * (i -1)
       call  fmod_set_seg_natoms(temp1, ivalue2, 1)
    enddo
  end subroutine fmod_set_mol_seg_natom
!----------------------------------------------------------------------
  subroutine fmod_set_mol_seg2atom(ivalue1,ivalue2,ivalue3)
    use segments, only : fmod_set_seg2atom
    implicit none
    integer(4), intent(in) :: ivalue1, ivalue2, ivalue3
    integer(4) :: temp1, temp2
    integer(4) :: i

    do i = 1, molinfo_nmoles(id)
       temp1 = segstart + ivalue1 + nmseg * (i-1)
       temp2 = ivalue3 + molinfo_natoms(id) * (i-1) + &
            & (molinfo_start(id) - 1)
       call fmod_set_seg2atom(temp1,ivalue2,temp2,1)
    enddo
  end subroutine fmod_set_mol_seg2atom
!----------------------------------------------------------------------
  subroutine fmod_set_mol_mol2atom
    use molecules, only : fmod_set_mol2atom
    implicit none
    integer(4), save :: atomnum=-1, molnum2=-1
    integer(4) :: i,j

    do i =  1,  molinfo_nmoles(id)
       molnum2 = molnum2 + 1
       do  j = 1, molinfo_natoms(id)
          atomnum = atomnum + 1
          call fmod_set_mol2atom(molnum2,j - 1, atomnum, 1)
       enddo
    enddo
  end subroutine fmod_set_mol_mol2atom
!----------------------------------------------------------------------
  subroutine fmod_set_mol_ID(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    id = ivalue + 1
  end subroutine fmod_set_mol_ID
!-----------------------------------------------------------------------
  subroutine fmod_set_mol_mol_natoms
    use molecules, only : fmod_set_mol_natoms
    implicit none
    integer(4) :: temp1
    integer(4) :: i
    integer(4) :: molnum = -1

    do i = 1, molinfo_nmoles(id)
       molnum = molnum + 1
       temp1 = molinfo_start(id) -1 +   (i-1)
       call fmod_set_mol_natoms(molnum, molinfo_natoms(id), 1)
    enddo
  end subroutine fmod_set_mol_mol_natoms

  subroutine add_to_mol_shake_list(num,ivalue0, ivalue, value)
    use parse_shake, only : add_to_shake_list
    use mpi_tool
    implicit none
    integer(4), intent(in) :: num, ivalue0, ivalue
    real(8), intent(in) :: value
    integer(4) :: temp1, temp2
    integer(4) :: i, zure

    zure = 0
    do i = 1, num-1
       zure = zure + molinfo_natoms(i)
    end do

    temp1 = ivalue0 + zure
    temp2 = ivalue  + zure

    call add_to_shake_list(temp1, temp2, value, 1)
  end subroutine add_to_mol_shake_list

end module mol_info
