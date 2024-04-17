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
!! \brief  Module and subroutines to parse parameters for 
!!         1-4 special interactions.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to parse 1-4 special interaction parameters.
!! \author Kensuke Iwahashi
!<
module parse_special
  ! LJ
  type yyparse_lj
     integer(4) :: atom1,atom2
     real(8) :: eps1,eps2,r1,r2
     type(yyparse_lj),pointer :: next
  end type yyparse_lj

  type(yyparse_lj),pointer :: yyparse_lj_top(:)
  integer(4),allocatable :: yyparse_lj_n(:)
  !^^^special^^^
  integer(4) :: yyparse_lj_total=0

  ! Coulomb
  type yyparse_coulomb
     integer(4) :: atom1,atom2
     real(8) :: charge1,charge2
     logical :: is_defined1,is_defined2
     type(yyparse_coulomb),pointer :: next
  end type yyparse_coulomb
  type(yyparse_coulomb),pointer :: yyparse_coulomb_top(:)
  integer(4),allocatable :: yyparse_coulomb_n(:)
  !^^^special^^^
  integer(4) :: yyparse_coulomb_total=0

contains

!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to keep LJ pair in liner list temporally.
!! \author Kensuke Iwahashi
!<
  subroutine add_to_lj_list(atom01,atom02,eps,r, line_number)
    use trajectory_org
    use segments
    implicit none
    integer(4) :: atom01, atom02
    integer(4) :: atom1,atom2,line_number
    real(8) :: eps,r
    type(yyparse_lj),pointer :: new1,new2,p,next
    if (seg_ready) then
       atom1 = irearrange(atom01 + 1)
       atom2 = irearrange(atom02 + 1)
    else
       call segment_is_not_ready
       return
    endif
    if (atom1 > n) then
       write(0,'(a,i0)')  'ERROR (line ', line_number, ') : '
       write(0,'(a,i0,a)')  &
            &    'The number of ljs is out of bounds.  '// &
            &    'It must be less than ', n, '.'
       call modylas_abort()
    endif
    if (atom2 > n) then
       write(0,'(a,i0)')  'ERROR (line ', line_number, ') : '
       write(0,'(a,i0,a)')  &
            &    'The number of ljs is out of bounds.  '// &
            &    'It must be less than ', n, '.'
       call modylas_abort()
    endif
    allocate(new1)
    nullify(new1%next)
    new1%atom1 = atom1
    new1%atom2 = atom2
    new1%eps1 = eps
    new1%r1 = r
    new1%eps2 = -1.0d0
    p => yyparse_lj_top(atom1)
    !     insert to list, sorting from smaller
    do while (.true.)
       next => p%next
       if (.not. associated(next)) then
          p%next => new1
          yyparse_lj_n(atom1) = yyparse_lj_n(atom1) + 1
          yyparse_lj_total = yyparse_lj_total + 1
          exit
       endif
       !       doubled must be merged
       if (atom2 == next%atom2) then
          next%eps1 = new1%eps1
          next%r1 = new1%r1
          deallocate(new1)
          exit
       else if (atom2 < next%atom2) then
          new1%next => p%next
          p%next => new1
          yyparse_lj_n(atom1) = yyparse_lj_n(atom1) + 1
          yyparse_lj_total = yyparse_lj_total + 1
          exit
       endif
       p => next
    enddo
    allocate(new2)
    nullify(new2%next)
    new2%atom1 = atom2
    new2%atom2 = atom1
    new2%eps2 = eps
    new2%r2 = r
    new2%eps1 = -1.0d0
    p => yyparse_lj_top(atom2)
    !     insert to list, sorting from smaller
    do while (.true.)
       next => p%next
       if (.not. associated(next)) then
          p%next => new2
          yyparse_lj_n(atom2) = yyparse_lj_n(atom2) + 1
          yyparse_lj_total = yyparse_lj_total + 1
          exit
       endif
       !       doubled must be merged
       if (atom1 == next%atom2) then
          next%eps2 = new2%eps2
          next%r2 = new2%r2
          deallocate(new2)
          exit
       else if (atom1 < next%atom2) then
          new2%next => p%next
          p%next => new2
          yyparse_lj_n(atom2) = yyparse_lj_n(atom2) + 1
          yyparse_lj_total = yyparse_lj_total + 1
          exit
       endif
       p => next
    enddo
  end subroutine add_to_lj_list

  subroutine add_to_mol_lj_list(ivalue1, atom1, atom2,  eps, r)
    use mol_info, only : para_start
    implicit none
    real(8), intent(in) :: eps, r
    integer(4), intent(in) :: ivalue1
    integer(4), intent(inout) :: atom1, atom2

    atom1 = atom1 + para_start(ivalue1)
    atom2 = atom2 + para_start(ivalue1)
    call add_to_lj_list(atom1-1, atom2-1, eps, r, 1)
  end subroutine add_to_mol_lj_list

!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to keep coulomb pair in liner list temporally.
!! \author Kensuke Iwahashi
!<
  subroutine add_to_coulomb_list(atom01,atom02,charge,line_number)
    use trajectory_org
    use segments
    implicit none
    integer(4) :: atom01, atom02
    integer(4) :: atom1,atom2,line_number
    real(8) :: charge
    type(yyparse_coulomb),pointer :: new1,new2,p,next
    if (seg_ready) then
       atom1 = irearrange(atom01 + 1)
       atom2 = irearrange(atom02 + 1)
    else
       call segment_is_not_ready
       return
    endif
    if (atom1 > n) then
       write(0,'(a,i0)')  'ERROR (line ', line_number, ') : '
       write(0,'(a,i0,a)')  &
            &    'The number of coulombs is out of bounds.  '// &
            &    'It must be less than ', n, '.'
       call modylas_abort()
    endif
    if (atom2 > n) then
       write(0,'(a,i0)')  'ERROR (line ', line_number, ') : '
       write(0,'(a,i0,a)')  &
            &    'The number of coulombs is out of bounds.  '// &
            &    'It must be less than ', n, '.'
       call modylas_abort()
    endif
    !	  making a new element of list
    !	  it must be normalized but must keep information to which atom
    !	  the charge belongs
    allocate(new1)
    nullify(new1%next)
    new1%atom1 = atom1
    new1%atom2 = atom2
    new1%charge1 = charge
    new1%is_defined1 = .true.
    new1%is_defined2 = .false.
    !     insert to list, sorting from smaller
    p => yyparse_coulomb_top(atom1)
    do while (.true.)
       next => p%next
       if (.not. associated(next)) then
          p%next => new1
          yyparse_coulomb_n(atom1)=yyparse_coulomb_n(atom1)+1
          yyparse_coulomb_total = yyparse_coulomb_total + 1
          exit
       endif
       !       doubled must be merged
       if (atom2 == next%atom2) then
          next%charge1 = charge
          next%is_defined1 = .true.
          deallocate(new1)
          exit
       else if (atom2 < next%atom2) then
          new1%next => p%next
          p%next => new1
          yyparse_coulomb_n(atom1)=yyparse_coulomb_n(atom1)+1
          yyparse_coulomb_total = yyparse_coulomb_total + 1
          exit
       endif
       p => next
    enddo
    !
    allocate(new2)
    nullify(new2%next)
    new2%atom1 = atom2
    new2%atom2 = atom1
    new2%charge2 = charge
    new2%is_defined1 = .false.
    new2%is_defined2 = .true.
    !     insert to list, sorting from smaller
    p => yyparse_coulomb_top(atom2)
    do while (.true.)
       next => p%next
       if (.not. associated(next)) then
          p%next => new2
          yyparse_coulomb_n(atom2)=yyparse_coulomb_n(atom2)+1
          yyparse_coulomb_total = yyparse_coulomb_total + 1
          exit
       endif
       !       doubled must be merged
       if (atom1 == next%atom2) then
          next%charge2 = charge
          next%is_defined2 = .true.
          deallocate(new2)
          exit
       else if (atom1 < next%atom2) then
          new2%next => p%next
          p%next => new2
          yyparse_coulomb_n(atom2)=yyparse_coulomb_n(atom2)+1
          yyparse_coulomb_total = yyparse_coulomb_total + 1
          exit
       endif
       p => next
    enddo
  end subroutine add_to_coulomb_list
!-----------------------------------------------------------------------
  subroutine add_to_mol_coulomb_list(ivalue1,atom1,atom2, charge)
    use mol_info
    implicit none
    real(8), intent(in)  :: charge
    integer(4), intent(in) :: ivalue1
    integer(4), intent(out) :: atom1,atom2

    atom1 = atom1 + para_start(ivalue1)
    atom2 = atom2 + para_start(ivalue1)

    call add_to_coulomb_list(atom1-1,atom2-1,charge,1)
  end subroutine add_to_mol_coulomb_list

end module parse_special
