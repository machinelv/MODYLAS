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
!! \brief  Module and subroutines to calculate dihedral potential.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to prepare index for dihedral table.
!! \author Kazushi Fujimoto, Zhiye Tang
!<
#ifdef DIHEDRAL_TABLE
module mod_dihedral_table_index
  implicit none

  type dihedral_table_index
     integer(4) :: atom1, atom2, atom3, atom4
     integer(4) :: group
  end type dihedral_table_index

  type dihedral_table_index_vector
     type(dihedral_table_index), allocatable :: items(:)
     integer(4) :: size, alloc
   contains
     procedure :: init      => init_index_vector
     procedure :: resize    => resize_index_vector
     procedure :: push_back => push_back_index_vector
  end type dihedral_table_index_vector

contains
  subroutine init_index_vector(v)
    implicit none
    class(dihedral_table_index_vector), intent(inout) :: v

    if( allocated(v%items) ) deallocate(v%items)
    v%size  = 0
    v%alloc = 0
  end subroutine init_index_vector

  subroutine resize_index_vector(v, size)
    implicit none
    class(dihedral_table_index_vector), intent(inout) :: v
    integer(4)                        , intent(in)    :: size
    type(dihedral_table_index), allocatable :: items_tmp(:)
    integer(4) :: n

    if( size <= v%alloc ) then
       v%size = size
       return
    end if

    if( v%alloc == 0 ) then
       v%alloc = 1
    end if
    do while(v%alloc < size)
       v%alloc = v%alloc * 2
    end do

    if( 0 == v%size ) then
       allocate(v%items(v%alloc))
    else
       allocate(items_tmp(v%size))
       items_tmp(1:v%size) = v%items(1:v%size)
       deallocate(v%items)
       allocate(v%items(v%alloc))
       v%items(1:v%size) = items_tmp(1:v%size)
       deallocate(items_tmp)
    end if

    v%size = size
  end subroutine resize_index_vector

  subroutine push_back_index_vector(v, index)
    implicit none
    class(dihedral_table_index_vector), intent(inout) :: v
    type(dihedral_table_index)        , intent(in)    :: index

    call v%resize(v%size + 1)

    v%items(v%size) = index

  end subroutine push_back_index_vector

  function is_same_index_atoms(index1, index2) result(same)
    implicit none
    type(dihedral_table_index), intent(in) :: index1, index2
    logical :: same
    
    same = .false.
    if( index1%atom1 /= index2%atom1 ) return
    if( index1%atom2 /= index2%atom2 ) return
    if( index1%atom3 /= index2%atom3 ) return
    if( index1%atom4 /= index2%atom4 ) return
    same = .true.
    
  end function is_same_index_atoms
  
end module mod_dihedral_table_index
!>
!! \brief  Module to input data into dihedral table.
!! \author Kazushi Fujimoto, Zhiye Tang
!<
module mod_dihedral_table_value
  use precision, only : wfp
  implicit none
  
  type dihedral_table_value
     real(kind=wfp) :: Kd, nd, delta
  end type dihedral_table_value

  type dihedral_table_value_vector
     type(dihedral_table_value), allocatable :: items(:)
     integer(4) :: size, alloc
   contains
     procedure :: init      => init_value_vector
     procedure :: resize    => resize_value_vector
     procedure :: push_back => push_back_value_vector
  end type dihedral_table_value_vector

contains
  subroutine init_value_vector(v)
    implicit none
    class(dihedral_table_value_vector), intent(inout) :: v

    if( allocated(v%items) ) deallocate(v%items)
    v%size  = 0
    v%alloc = 0
  end subroutine init_value_vector

  subroutine resize_value_vector(v, size)
    implicit none
    class(dihedral_table_value_vector), intent(inout) :: v
    integer(4)                        , intent(in)    :: size
    type(dihedral_table_value), allocatable :: items_tmp(:)
    integer(4) :: n

    if( size <= v%alloc ) then
       v%size = size
       return
    end if

    if( v%alloc == 0 ) then
       v%alloc = 1
    end if
    do while(v%alloc < size)
       v%alloc = v%alloc * 2
    end do

    if( 0 == v%size ) then
       allocate(v%items(v%alloc))
    else
       allocate(items_tmp(v%size))
       items_tmp(1:v%size) = v%items(1:v%size)
       deallocate(v%items)
       allocate(v%items(v%alloc))
       v%items(1:v%size) = items_tmp(1:v%size)
       deallocate(items_tmp)
    end if

    v%size = size
  end subroutine resize_value_vector

  subroutine push_back_value_vector(v, val)
    implicit none
    class(dihedral_table_value_vector), intent(inout) :: v
    type(dihedral_table_value)        , intent(in)    :: val

    call v%resize(v%size + 1)

    v%items(v%size) = val

  end subroutine push_back_value_vector

  function is_same_value(val1, val2) result(same)
    implicit none
    type(dihedral_table_value), intent(in) :: val1, val2
    logical :: same
    
    same = .false.
    if( val1%Kd    /= val2%Kd    ) return
    if( val1%nd    /= val2%nd    ) return
    if( val1%delta /= val2%delta ) return
    same = .true.
    
  end function is_same_value

  function is_same_value_group(values1, values2) result(same)
    implicit none
    type(dihedral_table_value), intent(in) :: values1(:), values2(:)
    logical :: same
    type(dihedral_table_value) :: val1, val2
    integer(4) :: i

    same = .false.

    if( size(values1) /= size(values2) ) return
       
    do i = 1, size(values1)
       val1 = values1(i)
       val2 = values2(i)
       if( .not. is_same_value(val1, val2) ) return
    end do

    same = .true.
    
  end function is_same_value_group

end module mod_dihedral_table_value
#endif /*DIHEDRAL_TABLE*/
!>
!! \brief  Module to parse parameters of dihedral potential.
!! \author Kazushi Fujimoto, Zhiye Tang
!<
module parse_dihedral
  use mpi_tool
  use precision, only : wfp
#ifdef DIHEDRAL_TABLE
  use mod_dihedral_table_index
  use mod_dihedral_table_value
#endif
  implicit none

  type yyparse_dihedralA
     integer(4) :: atom2, atom3, atom4
     real(kind=wfp) :: Kd, nd, delta
     type(yyparse_dihedralA),pointer :: next
  end type yyparse_dihedralA

  type(yyparse_dihedralA),pointer :: yyparse_dihedralA_top(:)
  integer(4),allocatable :: yyparse_dihedralA_n(:)
  integer(4) :: yyparse_dihedralA_total=0

  type yyparse_dihedralB
     integer(4) :: atom1, atom3, atom4
     real(kind=wfp) :: Kd, nd, delta
     type(yyparse_dihedralB),pointer :: next
  end type yyparse_dihedralB

  type(yyparse_dihedralB),pointer :: yyparse_dihedralB_top(:)
  integer(4),allocatable :: yyparse_dihedralB_n(:)
  integer(4) :: yyparse_dihedralB_total=0

#ifdef DIHEDRAL_TABLE
  type(dihedral_table_index_vector) :: v_index
  type(dihedral_table_value_vector) :: v_value
  type(dihedral_table_value_vector) :: v_value_tmp
  integer(4),allocatable            :: v_group_begin(:)

#ifdef DIHEDRAL_TABLE_DEBUG
  type dihedral_raw_data
     integer(4) :: atom1, atom2, atom3, atom4
     real(kind=wfp) :: Kd, nd, delta
  end type dihedral_raw_data

  type(dihedral_raw_data), allocatable :: raw_data(:)
#endif
#endif

contains
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to keep dihedral pair in liner list temporally.
!! \author Kensuke Iwahashi
!<
  subroutine add_to_dihedralA_list(atom1,atom2,atom3,atom4,Kd,nd,delta)
    use trajectory_org
    use param
    implicit none
    integer(4), intent(in) :: atom1,atom2,atom3,atom4
    real(kind=wfp), intent(in) :: Kd,nd,delta
    type(yyparse_dihedralA),pointer :: new1,p,next
    if (atom1 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of dihedral is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
    endif
    if (atom2 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of dihedral is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
    endif
    if (atom3 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of dihedral is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
    endif
    if (atom4 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of dihedral is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
    endif
    allocate(new1)
    nullify(new1%next)
    new1%atom2 = atom2
    new1%atom3 = atom3
    new1%atom4 = atom4
    new1%Kd = Kd
    new1%nd = nd
    new1%delta = delta
    p => yyparse_dihedralA_top(atom1)
    !     insert to list, sorting from smaller
    do while (.true.)
       next => p%next
       if (.not. associated(next)) then
          p%next => new1
          yyparse_dihedralA_n(atom1) = yyparse_dihedralA_n(atom1) + 1
          yyparse_dihedralA_total = yyparse_dihedralA_total + 1
          exit
       endif
! disable sorting for dihedral table debug
#ifndef DIHEDRAL_TABLE_DEBUG
       if (atom2 <= next%atom2 .and. atom3 <= next%atom3 .and. atom4 <= next%atom4) then
          new1%next => p%next
          p%next => new1
          yyparse_dihedralA_n(atom1) = yyparse_dihedralA_n(atom1) + 1
          yyparse_dihedralA_total = yyparse_dihedralA_total + 1
          exit
        endif
#endif
        p => next
     enddo
   end subroutine add_to_dihedralA_list
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to keep dihedral pair in liner list temporally.
!! \author Kensuke Iwahashi
!<
   subroutine add_to_dihedralB_list(atom2,atom1,atom3,atom4,Kd,nd,delta)
     use trajectory_org
     use param
     implicit none
     integer(4), intent(in) :: atom1,atom2,atom3,atom4
     real(kind=wfp), intent(in) :: Kd,nd,delta
     type(yyparse_dihedralB),pointer :: new1,p,next
     if (atom1 > npara) then
        write(0,'(a,i0,a)')  'ERROR: '// &
             &    'The number of dihedral is out of bounds.  '// &
             &    'It must be less than ', npara, '.'
     endif
     if (atom2 > npara) then
        write(0,'(a,i0,a)')  'ERROR: '// &
             &    'The number of dihedral is out of bounds.  '// &
             &    'It must be less than ', npara, '.'
     endif
     if (atom3 > npara) then
        write(0,'(a,i0,a)')  'ERROR: '// &
              &    'The number of dihedral is out of bounds.  '// &
              &    'It must be less than ', npara, '.'
      endif
      if (atom4 > npara) then
         write(0,'(a,i0,a)')  'ERROR: '// &
              &    'The number of dihedral is out of bounds.  '// &
              &    'It must be less than ', npara, '.'
      endif
      allocate(new1)
      nullify(new1%next)
      new1%atom1 = atom1
      new1%atom3 = atom3
      new1%atom4 = atom4
      new1%Kd = Kd
      new1%nd = nd
      new1%delta = delta
      p => yyparse_dihedralB_top(atom2)
      !     insert to list, sorting from smaller
      do while (.true.)
         next => p%next
         if (.not. associated(next)) then
            p%next => new1
            yyparse_dihedralB_n(atom2) = yyparse_dihedralB_n(atom2) + 1
            yyparse_dihedralB_total = yyparse_dihedralB_total + 1
          exit
       endif
! disable sorting for dihedral table debug
#ifndef DIHEDRAL_TABLE_DEBUG
       if (atom1 <= next%atom1 .and. atom3 <= next%atom3 .and. atom4 < next%atom4) then
          new1%next => p%next
          p%next => new1
          yyparse_dihedralB_n(atom2) = yyparse_dihedralB_n(atom2) + 1
          yyparse_dihedralB_total = yyparse_dihedralB_total + 1
          exit
       endif
#endif
       p => next
    enddo
  end subroutine add_to_dihedralB_list

end module parse_dihedral
!----------------------------------------------------------------------
!>
!! \brief  Module to calculate dihedral potential.
!! \author Kensuke Iwahashi, Kazushi Fujimoto, Zhiye Tang
!<
!----------------------------------------------------------------------
module dihedral
  use omp_lib
  use parse_dihedral
  use subcell
  use precision, only : wfp, MPI_FF_PARAM
  implicit none
  integer(4) :: ndihedrals=0
  integer(4) :: ndihedralAs=0, ndihedralBs=0
  integer(4), allocatable :: atom1(:), atom2(:), atom3(:), atom4(:)
  real(kind=wfp), allocatable :: Kd(:), nd(:), delta(:)
  integer(4), allocatable :: atom2A(:,:), atom3A(:,:), atom4A(:,:)
  real(kind=wfp), allocatable :: KdA(:,:), ndA(:,:), deltaA(:,:)
  integer(4), allocatable :: atom1B(:,:), atom3B(:,:), atom4B(:,:)
  real(kind=wfp), allocatable :: KdB(:,:), ndB(:,:), deltaB(:,:)
  integer(4), allocatable :: topA(:), nA(:)
  integer(4), allocatable :: topB(:), nB(:)
  integer(4), allocatable :: nArm(:),atom4Arm(:,:) !! AMBER/OPLS => all case
  integer(4), allocatable :: a4A_tmp(:) !! temporary array
    real(8), allocatable, save :: fi0_2(:,:),fi0_3(:,:,:) 
    real(8), allocatable, save :: w3_f1(:,:,:),w3_f2(:,:,:),w3_f3(:,:,:) 
    integer(4) :: jsize,jsize2,kblk  
    integer(4), allocatable, save :: jlist(:,:) 
    integer(4), allocatable, save :: i0a_2(:,:),i0c_2(:,:),i0d_2(:,:),i0_2(:,:) 
    integer(4), allocatable :: jmax(:) 
    integer(4) :: jmax_2,k0_w,k0_s,k0_e,k0_n,k0_b 
#ifdef DIHEDRAL_TABLE
  integer(4) :: ndihedralA_t
  integer(4), allocatable :: atom2A_t(:,:), atom3A_t(:,:), atom4A_t(:,:), indexA_t(:,:), nA_t(:)
  real(kind=wfp), allocatable :: KdA_t(:), ndA_t(:), deltaA_t(:)
  logical, allocatable :: continue_flagA_t(:)

  integer(4) :: ndihedralB_t
  integer(4), allocatable :: atom1B_t(:,:), atom3B_t(:,:), atom4B_t(:,:), indexB_t(:,:), nB_t(:)
  real(kind=wfp), allocatable :: KdB_t(:), ndB_t(:), deltaB_t(:)
  logical, allocatable :: continue_flagB_t(:)
#endif

contains
  !subroutine init_dihedrals
  !end subroutine init_dihedrals

  !subroutine intra_dihedrals
  !end subroutine intra_dihedrals

!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to check whether j-,k-,l-atoms are communicated correctly
!! \author Yoshimichi Andoh
!<
  subroutine check_atombound_dihedral()
!-----------------------------------------------------------------------
    use trajectory_mpi
    use param
    use md_condition
    use mpi_tool
    implicit none
    integer(4) :: k0,i0,j,i0a,i0b,i0c,i0d
    integer(4) :: iA,ipar
    integer(4) :: icheck

    icheck=0
#ifdef DIHEDRAL_TABLE
    if(allow_bond_breaking<=0.0d0) then
      
!$omp parallel default(none) &
!$omp& private(k0,i0,i0a,i0b,i0c,i0d,j,iA,ipar) &
!$omp& shared(nselfseg,lsegtop,lseg_natoms,m2i,i2m) &
!$omp& shared(paranum,nA_t,nB_t) &
!$omp& shared(atom1B_t,atom3B_t,atom4B_t,atom2A_t,atom3A_t,atom4A_t) &
!$omp& reduction(+:icheck)
!$omp do
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          iA   = m2i(i0)
          ipar = paranum(iA)
          do j = 1, nB_t(ipar)
             i0a = i2m(atom1B_t(ipar,j)+(iA-ipar))
             i0c = i2m(atom3B_t(ipar,j)+(iA-ipar))
             i0d = i2m(atom4B_t(ipar,j)+(iA-ipar))
             if(i0a.eq.-1) icheck = icheck -1
             if(i0c.eq.-1) icheck = icheck -1
             if(i0d.eq.-1) icheck = icheck -1
          enddo
       enddo
    enddo
!$omp end do
!$omp end parallel
   else
!$omp parallel default(none) &
!$omp& private(k0,i0,i0a,i0b,i0c,i0d,j,iA,ipar) &
!$omp& shared(nselfseg,lsegtop,lseg_natoms,m2i,i2m) &
!$omp& shared(paranum,nA_t,nB_t) &
!$omp& shared(atom1B_t,atom3B_t,atom4B_t,atom2A_t,atom3A_t,atom4A_t) &
!$omp& reduction(+:icheck) &
!$omp& shared(KdA_t, KdB_t,indexA_t,indexB_t)
!$omp do
      do k0=1,nselfseg
         do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
            iA   = m2i(i0)
            ipar = paranum(iA)
            do j = 1, nB_t(ipar)
               i0a = i2m(atom1B_t(ipar,j)+(iA-ipar))
               i0c = i2m(atom3B_t(ipar,j)+(iA-ipar))
               i0d = i2m(atom4B_t(ipar,j)+(iA-ipar))
               if(indexB_t(ipar,j)>=0 .and. i0a.eq.-1) icheck = icheck -1
               if(indexB_t(ipar,j)>=0 .and. i0c.eq.-1) icheck = icheck -1
               if(indexB_t(ipar,j)>=0 .and. i0d.eq.-1) icheck = icheck -1
            enddo
         enddo
      enddo
!$omp end do
!$omp end parallel
   endif

#else /* DIHEDRAL_TABLE */
    if(allow_bond_breaking<=0.0d0) then
      
!$omp parallel default(none) &
!$omp& private(k0,i0,i0a,i0b,i0c,i0d,j,iA,ipar) &
!$omp& shared(nselfseg,lsegtop,lseg_natoms,m2i,i2m) &
!$omp& shared(paranum,nA,nB) &
!$omp& shared(atom1B,atom3B,atom4B,atom2A,atom3A,atom4A) &
!$omp& reduction(+:icheck)
!$omp do
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          iA   = m2i(i0)
          ipar = paranum(iA)
          do j = 1, nB(ipar)
             i0a = i2m(atom1B(ipar,j)+(iA-ipar))
             i0c = i2m(atom3B(ipar,j)+(iA-ipar))
             i0d = i2m(atom4B(ipar,j)+(iA-ipar))
             if(i0a.eq.-1) icheck = icheck -1
             if(i0c.eq.-1) icheck = icheck -1
             if(i0d.eq.-1) icheck = icheck -1
          enddo
       enddo
    enddo
!$omp end do
!$omp end parallel
   else
!$omp parallel default(none) &
!$omp& private(k0,i0,i0a,i0b,i0c,i0d,j,iA,ipar) &
!$omp& shared(nselfseg,lsegtop,lseg_natoms,m2i,i2m) &
!$omp& shared(paranum,nA,nB) &
!$omp& shared(atom1B,atom3B,atom4B,atom2A,atom3A,atom4A) &
!$omp& reduction(+:icheck) &
!$omp& shared(KdA, KdB)
!$omp do
      do k0=1,nselfseg
         do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
            iA   = m2i(i0)
            ipar = paranum(iA)
            do j = 1, nB(ipar)
               i0a = i2m(atom1B(ipar,j)+(iA-ipar))
               i0c = i2m(atom3B(ipar,j)+(iA-ipar))
               i0d = i2m(atom4B(ipar,j)+(iA-ipar))
               if(KdB(ipar,j).ne.0.0d0 .and. i0a.eq.-1) icheck = icheck -1
               if(KdB(ipar,j).ne.0.0d0 .and. i0c.eq.-1) icheck = icheck -1
               if(KdB(ipar,j).ne.0.0d0 .and. i0d.eq.-1) icheck = icheck -1
            enddo
         enddo
      enddo
!$omp end do
!$omp end parallel
   endif
#endif /* DIHEDRAL_TABLE */

    if(icheck.le.-1) then
       write(0,*) 'ERROR[add_charmm_dihedral_a] &
            &   There is a particle outside the area.'
       write(0,*) 'Myrank,mdstep=',myrank,mdstep
       call modylas_abort()
    end if
  end subroutine check_atombound_dihedral

!----------------------------------------------------------------------
  subroutine fmod_alloc_yyparse_dihedA_top(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    integer(4) :: i
    type(yyparse_dihedralA),pointer :: top
    allocate(yyparse_dihedralA_top(ivalue))
    do i=1,ivalue
       allocate(top)
       top%atom2 = -1
       top%atom3 = -1
       top%atom4 = -1
       nullify(top%next)
       yyparse_dihedralA_top(i) = top
    enddo
  end subroutine fmod_alloc_yyparse_dihedA_top
  !----------------------------------------------------------------------
  subroutine fmod_alloc_yyparse_dihedA_n(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(yyparse_dihedralA_n(ivalue))
    yyparse_dihedralA_n = 0
  end subroutine fmod_alloc_yyparse_dihedA_n
  !-----------------------------------------------------------------------
  subroutine fmod_alloc_yyparse_dihedB_top(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    integer(4) :: i
    type(yyparse_dihedralB),pointer :: top
    allocate(yyparse_dihedralB_top(ivalue))
    do i=1,ivalue
       allocate(top)
       top%atom1 = -1
       top%atom3 = -1
       top%atom4 = -1
       nullify(top%next)
       yyparse_dihedralB_top(i) = top
    enddo
  end subroutine fmod_alloc_yyparse_dihedB_top
  !----------------------------------------------------------------------
  subroutine fmod_alloc_yyparse_dihedB_n(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(yyparse_dihedralB_n(ivalue))
    yyparse_dihedralB_n = 0
  end subroutine fmod_alloc_yyparse_dihedB_n
  !-----------------------------------------------------------------------
  subroutine add_to_mol_dihedral_list &
       &           (ivalue1, atom1, atom2, atom3, atom4, Kd, nd, delta)
    use mol_info, only : para_start
    implicit none
    integer(4), intent(in) :: ivalue1
    integer(4), intent(inout) :: atom1, atom2, atom3, atom4
    real(8), intent(in) :: Kd, nd, delta

    real(kind=wfp) :: Kd_t, nd_t, delta_t
    
    Kd_t=Kd; nd_t=nd; delta_t=delta
    

    atom1 = atom1 + para_start(ivalue1)
    atom2 = atom2 + para_start(ivalue1)
    atom3 = atom3 + para_start(ivalue1)
    atom4 = atom4 + para_start(ivalue1)

#ifdef DIHEDRAL_TABLE
    call add_to_dihedral_table(atom1, atom2, atom3, atom4, &
         &                     Kd_t, nd_t, delta_t)
#ifndef DIHEDRAL_TABLE_DEBUG
    return
#endif
#endif

    call add_to_dihedralA_list(atom1, atom2, atom3, atom4, &
         &                     Kd_t, nd_t, delta_t)
    call add_to_dihedralA_list(atom4, atom3, atom2, atom1, &
         &                     Kd_t, nd_t, delta_t)
    call add_to_dihedralB_list(atom2, atom1, atom3, atom4, &
         &                     Kd_t, nd_t, delta_t)
    call add_to_dihedralB_list(atom3, atom4, atom2, atom1, &
         &                     Kd_t, nd_t, delta_t)
  end subroutine add_to_mol_dihedral_list
!-----------------------------------------------------------------------
  subroutine set_md_charmm_dihedralA
    use param
    use force_field_numbers
    use md_condition
    use mpi_tool ! debug
    implicit none
    integer(4) :: i,j,k,maxi(1)
    integer(4) :: iA,jA,jj,jA_tmp
!for debug of cyclic motief
    integer(4) :: ia1,ia2,ia3,ia4
!
    type(yyparse_dihedralA),pointer :: p,freed
    ndihedralAs= yyparse_dihedralA_total
    if (yyparse_dihedralA_total == 0)  return
    maxi = maxloc(yyparse_dihedralA_n)
    maxi = yyparse_dihedralA_n(maxi(1))
    allocate(atom2A(npara,maxi(1)))
    allocate(atom3A(npara,maxi(1)))
    allocate(atom4A(npara,maxi(1)))
    allocate(KdA(npara,maxi(1)))
    allocate(ndA(npara,maxi(1)))
    allocate(deltaA(npara,maxi(1)))
    allocate(nA(npara))
    !     copying
    do i=1, npara
       nA(i) = yyparse_dihedralA_n(i)
       p => yyparse_dihedralA_top(i)
       do j=1,yyparse_dihedralA_n(i)
          p => p%next
          atom2A(i,j) = p%atom2
          atom3A(i,j) = p%atom3
          atom4A(i,j) = p%atom4
          KdA(i,j) = p%Kd
          ndA(i,j) = p%nd
          deltaA(i,j) = p%delta
       enddo
    enddo

!! set  nArm for scaling 1-4 in OPLS/AMBER
!   if(md_condition__force_field == OPLSAA .or. &
!        &   md_condition__force_field == AMBER  )then
       allocate(nArm(npara)) !; write(*,*) 'nArm allocated,',npara
       allocate(atom4Arm(npara,maxi(1)))
       allocate(a4A_tmp(maxi(1)))
       nArm=0; atom4Arm=0 ! initialize to 0
       jA_tmp=0; a4A_tmp=0

       do i=1, npara
          jA_tmp=0
          do j=1,nA(i)
             if(jA_tmp==0)then
                jA_tmp=1
                a4A_tmp(jA_tmp) = atom4A(i,j)
             else
                do jj=1,jA_tmp
                   if(a4A_tmp(jj)==atom4A(i,j))goto 9
                enddo
                jA_tmp=jA_tmp+1 !! not overlapped
                a4A_tmp(jA_tmp) = atom4A(i,j)
9               continue
             endif
          enddo ! j
!!            if(myrank==0) write(*,*) 'i,jA_tmp=',i-1,jA_tmp
! check cyclic motief
          jA=0
          do j=1,jA_tmp
          do jj=1,nA(i)
              if(a4A_tmp(j) == atom3A(i,jj))then  !! 5-ring case
!!                if(myrank==0) write(*,*) '4 and 3 overlap', a4A_tmp(j)-1
                  goto 8
              else
                  continue
              endif
          enddo ! jj
          jA=jA+1
          atom4Arm(i,jA)=a4A_tmp(j)
8         continue
          enddo ! j
!
          nArm(i)=jA
          !debug
#ifdef OPLS_DEBUG
          if(myrank==0) write(*,*) '##checking 1-4 pair list(debug)###'
          if(myrank==0) write(*,*) 'i,nArm(i)=',i-1,nArm(i) !! i-1 correspoinds to id in .mdff
          do j=1,jA
             if(myrank==0)write(*,*) i-1,atom4Arm(i,j)-1
          enddo
#endif
          !debug
       enddo ! i
!   endif

    deallocate(yyparse_dihedralA_top, yyparse_dihedralA_n)
  end subroutine set_md_charmm_dihedralA
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to read dihedral parameter B.
!! \author Kazushi FUJIMOTO
!<
  subroutine set_md_charmm_dihedralB
    use param
    implicit none
    integer(4) :: i,j,k,maxi(1)
    type(yyparse_dihedralB),pointer :: p,freed
    ndihedralBs= yyparse_dihedralB_total
    if (yyparse_dihedralB_total == 0)  return
    maxi = maxloc(yyparse_dihedralB_n)
    maxi = yyparse_dihedralB_n(maxi(1))
    allocate(atom1B(npara,maxi(1)))
    allocate(atom3B(npara,maxi(1)))
    allocate(atom4B(npara,maxi(1)))
    allocate(KdB(npara,maxi(1)))
    allocate(ndB(npara,maxi(1)))
    allocate(deltaB(npara,maxi(1)))
    allocate(nB(npara))
    !     copying
    do i=1, npara
       nB(i) = yyparse_dihedralB_n(i)
       p => yyparse_dihedralB_top(i)
       do j=1,yyparse_dihedralB_n(i)
          p => p%next
          atom1B(i,j) = p%atom1
          atom3B(i,j) = p%atom3
          atom4B(i,j) = p%atom4
          KdB(i,j) = p%Kd
          ndB(i,j) = p%nd
          deltaB(i,j) = p%delta
       enddo
    enddo

    deallocate(yyparse_dihedralB_top, yyparse_dihedralB_n)
  end subroutine set_md_charmm_dihedralB
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to read dihedral parameter A.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_md_charmm_dihedralA
    use trajectory_org
    use device_numbers, only : f_mdff
    use param
    use mpi_tool
    use force_field_numbers
    use md_condition
    implicit none
    include 'mpif.h'
    integer(4) :: ierr
    integer(4) :: num
    integer(4) :: i,j,jA,jj,jA_tmp

#ifndef DIHEDRAL_TABLE
    if(myrank==0) then
       read(f_mdff) ndihedralAs
    endif
    call MPI_Bcast(ndihedralAs,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    if(ndihedralAs.eq.0) return

    if(myrank==0) read(f_mdff) num
    call MPI_Bcast(num,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)

    allocate(nA(npara))
    allocate(atom2A(npara,num))
    allocate(atom3A(npara,num))
    allocate(atom4A(npara,num))
    allocate(KdA(npara,num))
    allocate(ndA(npara,num))
    allocate(deltaA(npara,num))

    if(myrank==0) then
       !      read(f_mdff) topA
       read(f_mdff) nA
       read(f_mdff) atom2A
       read(f_mdff) atom3A
       read(f_mdff) atom4A
       read(f_mdff) KdA
       read(f_mdff) ndA
       read(f_mdff) deltaA
    endif
    call MPI_Bcast(nA,npara,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom2A,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom3A,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom4A,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(KdA,npara*num,MPI_FF_PARAM,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(ndA,npara*num,MPI_FF_PARAM,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(deltaA,npara*num,MPI_FF_PARAM,0,MPI_COMM_WORLD,ierr)

!! set  nArm for scaling 1-4 in OPLS/AMBER
       allocate(nArm(npara)) !; write(*,*) 'nArm allocated,',npara
       allocate(atom4Arm(npara,num))
       allocate(a4A_tmp(num))

! copied from set_md_charmm_dihedralA
       do i=1, npara
          jA_tmp=0
          do j=1,nA(i)
             if(jA_tmp==0)then
                jA_tmp=1
                a4A_tmp(jA_tmp) = atom4A(i,j)
             else
                do jj=1,jA_tmp
                   if(a4A_tmp(jj)==atom4A(i,j))goto 9
                enddo
                jA_tmp=jA_tmp+1 !! not overlapped
                a4A_tmp(jA_tmp) = atom4A(i,j)
9               continue
             endif
          enddo ! j
!!        if(myrank==0) write(*,*) 'i,jA_tmp=',i-1,jA_tmp
! check cyclic motief
          jA=0
          do j=1,jA_tmp
          do jj=1,nA(i)
              if(a4A_tmp(j) == atom3A(i,jj))then  !! 5-ring case
!!                if(myrank==0) write(*,*) '4 and 3 overlap', a4A_tmp(j)-1
                  goto 8
              else
                  continue
              endif
          enddo ! jj
          jA=jA+1
          atom4Arm(i,jA)=a4A_tmp(j)
8         continue
          enddo ! j
!
          nArm(i)=jA
#ifdef OPLS_DEBUG
          if(myrank==0) write(*,*) '##checking 1-4 pair list(debug)###'
          if(myrank==0) write(*,*) 'i,nArm(i)=',i-1,nArm(i) !! i-1 correspoinds to id in .mdff
          do j=1,jA
             if(myrank==0)write(*,*) i-1,atom4Arm(i,j)-1
          enddo
#endif
       enddo ! i


#else /* DIHEDRAL_TABLE*/
   if(myrank==0) then
      read(f_mdff) ndihedralAs
   endif
   call MPI_Bcast(ndihedralAs,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
   if(ndihedralAs.eq.0) return

   if(myrank==0) then
      read(f_mdff) ndihedralA_t
   endif
   call MPI_Bcast(ndihedralA_t,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)

   if(myrank==0) read(f_mdff) num
   call MPI_Bcast(num,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)

    ! index
    allocate(nA_t(npara))
    allocate(atom2A_t(npara,num))
    allocate(atom3A_t(npara,num))
    allocate(atom4A_t(npara,num))
    allocate(indexA_t(npara,num))
    ! table
    allocate(KdA_t(ndihedralA_t))
    allocate(ndA_t(ndihedralA_t))
    allocate(deltaA_t(ndihedralA_t))
    allocate(continue_flagA_t(ndihedralA_t))

    if(myrank==0) then
       !      read(f_mdff) topB
       read(f_mdff) nA_t
       read(f_mdff) atom2A_t
       read(f_mdff) atom3A_t
       read(f_mdff) atom4A_t
       read(f_mdff) indexA_t
       read(f_mdff) KdA_t
       read(f_mdff) ndA_t
       read(f_mdff) deltaA_t
       read(f_mdff) continue_flagA_t
    endif
    call MPI_Bcast(nA_t,npara,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom2A_t,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom3A_t,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom4A_t,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(indexA_t,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(KdA_t,ndihedralA_t,MPI_FF_PARAM,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(ndA_t,ndihedralA_t,MPI_FF_PARAM,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(deltaA_t,ndihedralA_t,MPI_FF_PARAM,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(continue_flagA_t,ndihedralA_t,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

!! set  nArm for scaling 1-4 in OPLS/AMBER
allocate(nArm(npara)) !; write(*,*) 'nArm allocated,',npara
allocate(atom4Arm(npara,num))
allocate(a4A_tmp(num))

! copied from set_md_charmm_dihedralA
do i=1, npara
  jA_tmp=0
  do j=1,nA_t(i)
     if(jA_tmp==0)then
        jA_tmp=1
        a4A_tmp(jA_tmp) = atom4A_t(i,j)
     else
        do jj=1,jA_tmp
           if(a4A_tmp(jj)==atom4A_t(i,j))goto 9
        enddo
        jA_tmp=jA_tmp+1 !! not overlapped
        a4A_tmp(jA_tmp) = atom4A_t(i,j)
9               continue
     endif
  enddo ! j
!!        if(myrank==0) write(*,*) 'i,jA_tmp=',i-1,jA_tmp
! check cyclic motief
  jA=0
  do j=1,jA_tmp
  do jj=1,nA_t(i)
      if(a4A_tmp(j) == atom3A_t(i,jj))then  !! 5-ring case
!!                if(myrank==0) write(*,*) '4 and 3 overlap', a4A_tmp(j)-1
          goto 8
      else
          continue
      endif
  enddo ! jj
  jA=jA+1
  atom4Arm(i,jA)=a4A_tmp(j)
8         continue
  enddo ! j
!
  nArm(i)=jA
  !debug
#ifdef OPLS_DEBUG
  if(myrank==0) write(*,*) '##checking 1-4 pair list(debug)###'
  if(myrank==0) write(*,*) 'i,nArm(i)=',i-1,nArm(i) !! i-1 correspoinds to id in .mdff
  do j=1,jA
     if(myrank==0)write(*,*) i-1,atom4Arm(i,j)-1
  enddo
#endif
  !debug
enddo ! i

#endif /* DIHEDRAL_TABLE*/

  end subroutine read_md_charmm_dihedralA
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write dihedral parameter A in mdff.bin File.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_md_charmm_dihedralA
    use device_numbers, only : f_mdff
    implicit none
    integer(4) :: num
#ifdef DIHEDRAL_TABLE
    write(f_mdff) ndihedralAs
    if(ndihedralAs.eq.0) return
    write(f_mdff) ndihedralA_t
    num = size(atom2A_t(:,:), 2)
    write(f_mdff) num
    write(f_mdff) nA_t
    write(f_mdff) atom2A_t
    write(f_mdff) atom3A_t
    write(f_mdff) atom4A_t
    write(f_mdff) indexA_t
    write(f_mdff) KdA_t
    write(f_mdff) ndA_t
    write(f_mdff) deltaA_t
    write(f_mdff) continue_flagA_t
#else /* DIHEDRAL_TABLE */
    write(f_mdff) ndihedralAs
    if(ndihedralAs.eq.0) return
    num = size(atom2A(:,:), 2)
    write(f_mdff) num
    write(f_mdff) nA
    write(f_mdff) atom2A
    write(f_mdff) atom3A
    write(f_mdff) atom4A
    write(f_mdff) KdA
    write(f_mdff) ndA
    write(f_mdff) deltaA
#endif /* DIHEDRAL_TABLE */
  end subroutine write_md_charmm_dihedralA
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write dihedral parameter A.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_md_charmm_dihedralA
    implicit none
    integer(4) :: num
    write(*,*) '[write_md_charmm_dihedralA]'
    write(*,*) ndihedralAs
    if(ndihedralAs.eq.0) return
    num = size(atom2A(:,:), 2)
    write(*,*) num
    write(*,*) nA
    write(*,*) atom2A
    write(*,*) atom3A
    write(*,*) atom4A
    write(*,*) KdA
    write(*,*) ndA
    write(*,*) deltaA
  end subroutine write_memory_md_charmm_dihedralA
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to read dihedral parameter B.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_md_charmm_dihedralB
    use trajectory_org
    use device_numbers, only : f_mdff
    use param
    use mpi_tool
    implicit none
    include 'mpif.h'
    integer(4) :: ierr
    integer(4) :: num

#ifndef DIHEDRAL_TABLE
    if(myrank==0) then
       read(f_mdff) ndihedralBs
    endif
    call MPI_Bcast(ndihedralBs,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    if(ndihedralBs.eq.0) return

    if(myrank==0) read(f_mdff) num
    call MPI_Bcast(num,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)

    allocate(nB(npara))
    allocate(atom1B(npara,num))
    allocate(atom3B(npara,num))
    allocate(atom4B(npara,num))
    allocate(KdB(npara,num))
    allocate(ndB(npara,num))
    allocate(deltaB(npara,num))

    if(myrank==0) then
       !      read(f_mdff) topB
       read(f_mdff) nB
       read(f_mdff) atom1B
       read(f_mdff) atom3B
       read(f_mdff) atom4B
       read(f_mdff) KdB
       read(f_mdff) ndB
       read(f_mdff) deltaB
    endif
    call MPI_Bcast(nB,npara,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom1B,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom3B,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom4B,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(KdB,npara*num,MPI_FF_PARAM,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(ndB,npara*num,MPI_FF_PARAM,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(deltaB,npara*num,MPI_FF_PARAM,0,MPI_COMM_WORLD,ierr)
#else /* DIHERAL_TABLE*/
   if(myrank==0) then
      read(f_mdff) ndihedralBs
   endif
   call MPI_Bcast(ndihedralBs,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
   if(ndihedralBs.eq.0) return

   if(myrank==0) then
      read(f_mdff) ndihedralB_t
   endif
   call MPI_Bcast(ndihedralB_t,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)

   if(myrank==0) read(f_mdff) num
   call MPI_Bcast(num,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)

    ! index
    allocate(nB_t(npara))
    allocate(atom1B_t(npara,num))
    allocate(atom3B_t(npara,num))
    allocate(atom4B_t(npara,num))
    allocate(indexB_t(npara,num))
    ! table
    allocate(KdB_t(ndihedralB_t))
    allocate(ndB_t(ndihedralB_t))
    allocate(deltaB_t(ndihedralB_t))
    allocate(continue_flagB_t(ndihedralB_t))

    if(myrank==0) then
       !      read(f_mdff) topB
       read(f_mdff) nB_t
       read(f_mdff) atom1B_t
       read(f_mdff) atom3B_t
       read(f_mdff) atom4B_t
       read(f_mdff) indexB_t
       read(f_mdff) KdB_t
       read(f_mdff) ndB_t
       read(f_mdff) deltaB_t
       read(f_mdff) continue_flagB_t
    endif
    call MPI_Bcast(nB_t,npara,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom1B_t,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom3B_t,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom4B_t,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(indexB_t,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(KdB_t,ndihedralB_t,MPI_FF_PARAM,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(ndB_t,ndihedralB_t,MPI_FF_PARAM,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(deltaB_t,ndihedralB_t,MPI_FF_PARAM,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(continue_flagB_t,ndihedralB_t,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

#endif /* DIHERAL_TABLE*/

  end subroutine read_md_charmm_dihedralB
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write dihedral parameter B in mdff.bin File.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_md_charmm_dihedralB
    use device_numbers, only : f_mdff
    implicit none
    integer(4) :: num
#ifdef DIHEDRAL_TABLE
    write(f_mdff) ndihedralBs
    if(ndihedralBs.eq.0) return
    write(f_mdff) ndihedralB_t
    num = size(atom1B_t(:,:), 2)
    write(f_mdff) num
    write(f_mdff) nB_t
    write(f_mdff) atom1B_t
    write(f_mdff) atom3B_t
    write(f_mdff) atom4B_t
    write(f_mdff) indexB_t
    write(f_mdff) KdB_t
    write(f_mdff) ndB_t
    write(f_mdff) deltaB_t
    write(f_mdff) continue_flagB_t
#else /* DIHEDRAL_TABLE */
    write(f_mdff) ndihedralBs
    if(ndihedralBs.eq.0) return
    num = size(atom1B(:,:), 2)
    write(f_mdff) num
    write(f_mdff) nB
    write(f_mdff) atom1B
    write(f_mdff) atom3B
    write(f_mdff) atom4B
    write(f_mdff) KdB
    write(f_mdff) ndB
    write(f_mdff) deltaB
#endif /* DIHEDRAL_TABLE */
  end subroutine write_md_charmm_dihedralB
!------------------------------------------------------------------------
!>
!! \brief  Subroutine to write dihedral parameter B.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_md_charmm_dihedralB
    implicit none
    integer(4) :: num
    write(*,*) '[write_md_charmm_dihedralB]'
    write(*,*) ndihedralBs
    if(ndihedralBs.eq.0) return
    num = size(atom1B(:,:), 2)
    write(*,*) num
    write(*,*) nB
    write(*,*) atom1B
    write(*,*) atom3B
    write(*,*) atom4B
    write(*,*) KdB
    write(*,*) ndB
    write(*,*) deltaB
  end subroutine write_memory_md_charmm_dihedralB
!-------------------------------------------------------------------------
  subroutine read_md_charmm_dihedral
    !fujimoto
    use device_numbers, only : f_mdff
    implicit none

    read(f_mdff) ndihedrals
    if(ndihedrals.eq.0) return
    allocate(atom1(ndihedrals))
    allocate(atom2(ndihedrals))
    allocate(atom3(ndihedrals))
    allocate(atom4(ndihedrals))
    allocate(Kd(ndihedrals))
    allocate(nd(ndihedrals))
    allocate(delta(ndihedrals))
    read(f_mdff) atom1
    read(f_mdff) atom2
    read(f_mdff) atom3
    read(f_mdff) atom4
    read(f_mdff) Kd
    read(f_mdff) nd
    read(f_mdff) delta
  end subroutine read_md_charmm_dihedral
!-----------------------------------------------------------------------
!>
!! \brief  Wrapper subroutine of dihedral_a
!! \author Zhiye Tang
!<
  subroutine add_charmm_dihedral_a()
!-----------------------------------------------------------------------
   use md_condition, only : allow_bond_breaking, bNonEvenDihedralPotential

   if(bNonEvenDihedralPotential) then
      if(allow_bond_breaking<=0.0d0) then
         call add_charmm_dihedral_a_no_bond_breaking()
      else
         call add_charmm_dihedral_a_bond_breaking()
      endif
   else
      if(allow_bond_breaking<=0.0d0) then
         call add_charmm_dihedral_a_no_bond_breaking_stable()
      else
         call add_charmm_dihedral_a_bond_breaking_stable()
      endif
   endif

  end subroutine add_charmm_dihedral_a
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate energy and force of dihedral.
!! \author Kensuke Iwahashi, Zhiye Tang
!<
  subroutine add_charmm_dihedral_a_no_bond_breaking()
!-----------------------------------------------------------------------
!
!   /* V = Kd(1 + cos[n.chi - delta])
!    *
!    *            D
!    *           /   chi = angle between plain ABC and plain BCD
!    *     B----C
!    *   /
!    *  A
!    *
!    */
!
!   /* aliases
!    * delta => d, chi => c */
!
    use trajectory_org
    use trajectory_mpi
    use atom_virial

    use forces
    use md_monitors
    use md_const
    use param
    use mpi_tool
    use md_condition, only : allow_bond_breaking
#ifdef ONEPROC_AXIS
    use boundary, only : pbc_pair
#endif
    use openmp_tool, only : nomp
    use comm_direct2_dr, only : comm_direct_2_dr
    implicit none
    integer(4) :: i, j, ip, ia, ib, i0a, i0b, i0c, i0d, i0, ipar, k0
    real(8) :: a, b, c, coef, tmp_angle
    real(8) :: tmp__x, tmp__y, tmp__z, i_norm
    real(8) :: tmp1__x, tmp1__y, tmp1__z, tmp2__x, tmp2__y, tmp2__z
    real(8) :: vBA__x, vBA__y, vBA__z, vBC__x, vBC__y, vBC__z
    real(8) :: vCD__x, vCD__y, vCD__z, vCA__x, vCA__y, vCA__z
    real(8) :: v1__x, v1__y, v1__z, v2__x, v2__y, v2__z
    real(8) :: v3__x, v3__y, v3__z
    real(8) :: Fa__x, Fa__y, Fa__z, Fb__x, Fb__y, Fb__z
    real(8) :: Fc__x, Fc__y, Fc__z, Fd__x, Fd__y, Fd__z
    real(8) :: Udihedral,oUd
    real(8) :: v11,v22,v33
    real(8) :: wk_v11, wk_v22, wk_v33
    real(8) :: sfact=1.0d0  ! scaling factor to potE and virial (default)
    real(8) :: wkvd(6), vd(6) ! debug
    real(8) :: v21,v31,v32
    real(8) :: wk_v21, wk_v31, wk_v32
    real(8) :: ss1,ss2,ss3,r1,r2,r3,c0,c1,c2,s1,s2,s12,s
    real(8) :: df,a11,a22,a33,a12,a13,a23,sx2,sy2,sz2
    real(8) :: vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb4x,vb4y,vb4z,vb1xvb2,vb3xvb2,vAx,vAy,vAz,vBx,vBy,vBz,rvA2,rvB2
#ifdef DIHEDRAL_TABLE
    integer(4) :: k
#endif
    real(8) :: fi0(3)
!
    integer(4) :: iam
    include 'mpif.h'

!default: MTD
    if( nselfseg .eq. 0 ) then
        k0_w=1
        k0_n=0
    else
        k0_w=(nselfseg-1)/nomp + 1  ! note: nselfseg is variable with MD cycle
        k0_n=(nselfseg-1)/k0_w + 1  ! note: nselfseg is variable with MD cycle
    endif

    call check_atombound_dihedral()

    iam = 0
    Udihedral = 0.0d0
    wk_v11=0d0;wk_v22=0d0;wk_v33=0d0
    wk_v21=0d0;wk_v31=0d0;wk_v32=0d0

!$omp parallel default(none) &
!$omp& private(v11,v22,v33) &
!$omp& reduction(+:wk_v11,wk_v22,wk_v33) &
!$omp& private(v21,v31,v32) &
!$omp& reduction(+:wk_v21,wk_v31,wk_v32) &
!$omp& private(k0,i,iam,j,i0,ipar,jsize) &
!$omp& private(fi0,tmp_angle) &
!$omp& private(ip,ia,ib,i0a,i0b,i0c,i0d,vBA__x,vBA__y,vBA__z) &
!$omp& private(vBC__x,vBC__y,vBC__z,vCD__x,vCD__y,vCD__z) &
!$omp& private(vCA__x,vCA__y,vCA__z,v1__x,v1__y,v1__z) &
!$omp& private(v2__x,v2__y,v2__z,c,tmp__x,tmp__y,tmp__z) &
!$omp& private(coef,v3__x,v3__y,v3__z,a,b,i_norm) &
!$omp& private(Fa__x,Fa__y,Fa__z,tmp1__x,tmp1__y,tmp1__z) &
!$omp& private(tmp2__x,tmp2__y,tmp2__z) &
!$omp& private(Fb__x,Fb__y,Fb__z,Fc__x,Fc__y,Fc__z) &
!$omp& private(Fd__x,Fd__y,Fd__z) &
!$omp& private(jmax_2,k0_b,k0_s,k0_e) &
!$omp& shared(i0a_2,i0c_2,i0d_2,i0_2,fi0_2,fi0_3) &
!$omp& shared(w3_f1,w3_f2,w3_f3,jlist,jmax) &
#ifdef DIHEDRAL_TABLE
!$omp& shared(m2i,i2m,nA_t,topA,atom2A_t,atom3A_t,atom4A_t,indexA_t, continue_flagA_t) &
!$omp& shared(wkxyz,w3_f,KdA_t,ndA_t,deltaA_t,wk_vir2,ndihedralB_t,ndihedralA_t) &
!$omp& shared(nB_t,topB,atom1B_t,atom3B_t,atom4B_t,indexB_t,KdB_t,continue_flagB_t,ndB_t,deltaB_t,paranum) &
#else
!$omp& shared(m2i,i2m,nA,topA,atom2A,atom3A,atom4A) &
!$omp& shared(wkxyz,w3_f,KdA,ndA,deltaA,wk_vir2) &
!$omp& shared(nB,topB,atom1B,atom3B,atom4B,KdB,ndB,deltaB,paranum) &
#endif
!$omp& shared(tag,na_per_cell) &
!$omp& shared(lsegtop,lseg_natoms,nselfseg) &
!$omp& shared(sfact,k0_n,k0_w) &
#ifdef DIHEDRAL_TABLE
!$omp& private(k) &
#endif
!$omp& shared(myrank) &
!$omp& private(ss1,ss2,ss3,r1,r2,r3,c0,c1,c2,s1,s2,s12,s) &
!$omp& private(df,a11,a22,a33,a12,a13,a23,sx2,sy2,sz2) &
!$omp& private(vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb4x,vb4y,vb4z,vb1xvb2,vb3xvb2,vAx,vAy,vAz,vBx,vBy,vBz,rvA2,rvB2) &
!$omp& reduction(+:Udihedral)
!$  iam = omp_get_thread_num()
!$omp do
!default: MTD
    do k0_b=1,k0_n
        k0_s=(k0_b-1)*k0_w + 1 
        k0_e=k0_s + k0_w - 1 
        k0_e=min(k0_e,nselfseg) 
        do k0=k0_s,k0_e 
            do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
                ia   = m2i(i0)
                ipar = paranum(ia)
                fi0(1:3)=0d0
!default: MTD
                fi0_2(1:3,i0)=0d0
!default: MTD
            enddo ! enddo i0
        enddo ! endo k0
        jsize=0 
        do k0=k0_s,k0_e 
            do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
                ia   = m2i(i0)
                ipar = paranum(ia)
#ifdef DIHEDRAL_TABLE
                do j = 1, nB_t(ipar)
                    if(ia < atom3B_t(ipar,j)+(ia-ipar)) cycle
                    jsize=jsize+1
                    i0a = i2m(atom1B_t(ipar,j)+(ia-ipar))
                    i0c = i2m(atom3B_t(ipar,j)+(ia-ipar))
                    i0d = i2m(atom4B_t(ipar,j)+(ia-ipar))
                    i0a_2(jsize,k0_b)= i0a
                    i0c_2(jsize,k0_b)= i0c
                    i0d_2(jsize,k0_b)= i0d
                    i0_2(jsize,k0_b) = i0
                    jlist(jsize,k0_b)=j
                enddo
#else /* DIHEDRAL_TABLE */
                do j = 1, nB(ipar)
                    if(ia < atom3B(ipar,j)+(ia-ipar)) cycle
                    jsize=jsize+1
                    i0a = i2m(atom1B(ipar,j)+(ia-ipar))
                    i0c = i2m(atom3B(ipar,j)+(ia-ipar))
                    i0d = i2m(atom4B(ipar,j)+(ia-ipar))
                    i0a_2(jsize,k0_b)= i0a
                    i0c_2(jsize,k0_b)= i0c
                    i0d_2(jsize,k0_b)= i0d
                    i0_2(jsize,k0_b) = i0
                    jlist(jsize,k0_b)=j
                enddo
#endif /* DIHEDRAL_TABLE */
            enddo ! enddo i0
        enddo ! enddo k0
        jmax_2=jsize

!default: MTD
! in this branch w/ and w/o DIHEDRAL_TABLE it does the same thing
                do jsize=1,jmax_2 
                    j=jlist(jsize,k0_b) 

!default: MTD
                    i0a = i0a_2(jsize,k0_b) 
                    i0c = i0c_2(jsize,k0_b) 
                    i0d = i0d_2(jsize,k0_b) 
                    i0  = i0_2(jsize,k0_b) 
                    ia   = m2i(i0)
                    ipar = paranum(ia)
#ifdef DIHEDRAL_TABLE
                    k = indexB_t(ipar,j)
                     fi0_3(1,jsize,k0_b) = 0.0d0
                     fi0_3(2,jsize,k0_b) = 0.0d0
                     fi0_3(3,jsize,k0_b) = 0.0d0
                     w3_f1(1,jsize,iam) = 0.0d0
                     w3_f1(2,jsize,iam) = 0.0d0
                     w3_f1(3,jsize,iam) = 0.0d0
                     w3_f2(1,jsize,iam) = 0.0d0
                     w3_f2(2,jsize,iam) = 0.0d0
                     w3_f2(3,jsize,iam) = 0.0d0
                     w3_f3(1,jsize,iam) = 0.0d0
                     w3_f3(2,jsize,iam) = 0.0d0
                     w3_f3(3,jsize,iam) = 0.0d0
                     
                    do while ( .true. )
#endif /* DIHEDRAL_TABLE */
                        vb1x = wkxyz(1,i0a) - wkxyz(1,i0 )
                        vb1y = wkxyz(2,i0a) - wkxyz(2,i0 )
                        vb1z = wkxyz(3,i0a) - wkxyz(3,i0 )
                        vb2x = wkxyz(1,i0c) - wkxyz(1,i0 )
                        vb2y = wkxyz(2,i0c) - wkxyz(2,i0 )
                        vb2z = wkxyz(3,i0c) - wkxyz(3,i0 )
                        vb3x = wkxyz(1,i0d) - wkxyz(1,i0c)
                        vb3y = wkxyz(2,i0d) - wkxyz(2,i0c)
                        vb3z = wkxyz(3,i0d) - wkxyz(3,i0c)
                        vb4x = wkxyz(1,i0a) - wkxyz(1,i0c)
                        vb4y = wkxyz(2,i0a) - wkxyz(2,i0c)
                        vb4z = wkxyz(3,i0a) - wkxyz(3,i0c)
             
#ifdef ONEPROC_AXIS

                        call pbc_pair(vb1x,vb1y,vb1z)
                        call pbc_pair(vb2x,vb2y,vb2z)
                        call pbc_pair(vb3x,vb3y,vb3z)
                        call pbc_pair(vb4x,vb4y,vb4z)
#endif

                        ss1 = 1.0d0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z)
                        ss2 = 1.0d0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z)
                        ss3 = 1.0d0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z)

                        r1 = sqrt(ss1);
                        r2 = sqrt(ss2);
                        r3 = sqrt(ss3);

                        vb1xvb2 = vb1x * vb2x + vb1y * vb2y + vb1z * vb2z
                        vb3xvb2 = vb3x * vb2x + vb3y * vb2y + vb3z * vb2z

                        c0 = (vb1x * vb3x + vb1y * vb3y + vb1z * vb3z) * r1 * r3
                        c1 = vb1xvb2 * r1 * r2
                        c2 = -vb3xvb2 * r3 * r2

                        s1 = 1.0d0 - c1*c1;
                        s1 = max(s1, 0.001d0)
                        s1 = 1.0d0 / s1

                        s2 = 1.0d0 - c2*c2;
                        s2 = max(s2, 0.001d0)
                        s2 = 1.0d0 / s2

                        s12 = sqrt(s1*s2);
                        c = (c1*c2 + c0) * s12;
                        c=min(+1d0,c)
                        c=max(-1d0,c)

                        ! vb1 cross vb3 * vb2 > 0 => +
                        ! vb1 cross vb3 * vb2 <=0 => -
                        c = acos(c)
                        tmp1__x = vb1y * vb3z - vb1z * vb3y
                        tmp1__y = vb1z * vb3x - vb1x * vb3z
                        tmp1__z = vb1x * vb3y - vb1y * vb3x
                        c= sign(c, tmp1__x * vb2x + tmp1__y * vb2y + tmp1__z * vb2z)
                        ! at this point, c means the dihedral angle theta

#ifdef DIHEDRAL_TABLE
                        !         force coefficient
                        tmp_angle = ndB_t(k)*c - deltaB_t(k)
                        coef = KdB_t(k) * ndB_t(k) * sin(tmp_angle)
#else /* DIHEDRAL_TABLE */
                        !         force coefficient
                        tmp_angle = ndB(ipar,j)*c - deltaB(ipar,j)
                        coef = KdB(ipar,j) * ndB(ipar,j) * sin(tmp_angle)
#endif /* DIHEDRAL_TABLE */


                        vAx = -vb1y * vb2z + vb1z * vb2y;
                        vAy = -vb1z * vb2x + vb1x * vb2z;
                        vAz = -vb1x * vb2y + vb1y * vb2x;
                        vBx = -vb3y * vb2z + vb3z * vb2y;
                        vBy = -vb3z * vb2x + vb3x * vb2z;
                        vBz = -vb3x * vb2y + vb3y * vb2x;
                        rvA2 = 1.0e0 / (vAx * vAx + vAy * vAy + vAz * vAz);
                        rvB2 = 1.0e0 / (vBx * vBx + vBy * vBy + vBz * vBz);

                        a11 = (-coef / r2 * rvA2 );
                        a22 = (-vb1xvb2 * rvA2 * r2) * coef;
                        a33 = ( vb3xvb2 * rvB2 * r2) * coef;
                        a13 = (coef / r2 * rvB2);

!default: MTD
                        Fa__x = vAx * a11;
                        Fa__y = vAy * a11;
                        Fa__z = vAz * a11;
                        Fb__x = - Fa__x + vAx * a22 + vBx * a33;
                        Fb__y = - Fa__y + vAy * a22 + vBy * a33;
                        Fb__z = - Fa__z + vAz * a22 + vBz * a33;

             !         calculating & adding forces
#ifndef DIHEDRAL_TABLE
                        fi0_3(1,jsize,k0_b)=Fb__x
                        fi0_3(2,jsize,k0_b)=Fb__y
                        fi0_3(3,jsize,k0_b)=Fb__z
#else
                        fi0_3(1,jsize,k0_b)=fi0_3(1,jsize,k0_b)+Fb__x
                        fi0_3(2,jsize,k0_b)=fi0_3(2,jsize,k0_b)+Fb__y
                        fi0_3(3,jsize,k0_b)=fi0_3(3,jsize,k0_b)+Fb__z
#endif
                        Fd__x = vBx * a13;
                        Fd__y = vBy * a13;
                        Fd__z = vBz * a13;

                        Fc__x = - Fd__x - vAx * a22 - vBx * a33;
                        Fc__y = - Fd__y - vAy * a22 - vBy * a33;
                        Fc__z = - Fd__z - vAz * a22 - vBz * a33;

                        v11 = (vb4x*Fa__x - vb2x*Fb__x + vb3x*Fd__x)
                        v22 = (vb4y*Fa__y - vb2y*Fb__y + vb3y*Fd__y)
                        v33 = (vb4z*Fa__z - vb2z*Fb__z + vb3z*Fd__z)
                        v21 = (vb4y*Fa__x - vb2y*Fb__x + vb3y*Fd__x)
                        v31 = (vb4z*Fa__x - vb2z*Fb__x + vb3z*Fd__x)
                        v32 = (vb4z*Fa__y - vb2z*Fb__y + vb3z*Fd__y)
                        wk_v11 = wk_v11 + v11
                        wk_v22 = wk_v22 + v22
                        wk_v33 = wk_v33 + v33
                        wk_v21 = wk_v21 + v21
                        wk_v31 = wk_v31 + v31
                        wk_v32 = wk_v32 + v32
#ifndef DIHEDRAL_TABLE
                        w3_f1(1,jsize,iam) = Fa__x
                        w3_f1(2,jsize,iam) = Fa__y
                        w3_f1(3,jsize,iam) = Fa__z
                        w3_f2(1,jsize,iam) = Fc__x
                        w3_f2(2,jsize,iam) = Fc__y
                        w3_f2(3,jsize,iam) = Fc__z
                        w3_f3(1,jsize,iam) = Fd__x
                        w3_f3(2,jsize,iam) = Fd__y
                        w3_f3(3,jsize,iam) = Fd__z
#else
                        w3_f1(1,jsize,iam) = w3_f1(1,jsize,iam) + Fa__x
                        w3_f1(2,jsize,iam) = w3_f1(2,jsize,iam) + Fa__y
                        w3_f1(3,jsize,iam) = w3_f1(3,jsize,iam) + Fa__z
                        w3_f2(1,jsize,iam) = w3_f2(1,jsize,iam) + Fc__x
                        w3_f2(2,jsize,iam) = w3_f2(2,jsize,iam) + Fc__y
                        w3_f2(3,jsize,iam) = w3_f2(3,jsize,iam) + Fc__z
                        w3_f3(1,jsize,iam) = w3_f3(1,jsize,iam) + Fd__x
                        w3_f3(2,jsize,iam) = w3_f3(2,jsize,iam) + Fd__y
                        w3_f3(3,jsize,iam) = w3_f3(3,jsize,iam) + Fd__z
#endif

#ifdef DIHEDRAL_TABLE
                        Udihedral=Udihedral+KdB_t(k)*(1d0+cos(tmp_angle))
#else /* DIHEDRAL_TABLE */
                        Udihedral=Udihedral+KdB(ipar,j)*(1d0+cos(tmp_angle))
#endif /* DIHEDRAL_TABLE */
             !         Potential and virial are calculated in loop A.
#ifdef DIHEDRAL_TABLE
                        if( .not. continue_flagB_t(k) ) exit
                        k=k+1
                        if( k > ndihedralB_t ) exit ! exit loop if not continue the same group
                    enddo ! enddo while (true)
#endif /* DIHEDRAL_TABLE */
                enddo ! enddo j
!default: MTD
                do jsize=1,jmax_2 
                    i0a = i0a_2(jsize,k0_b) 
                    i0c = i0c_2(jsize,k0_b) 
                    i0d = i0d_2(jsize,k0_b) 
                    i0  = i0_2(jsize,k0_b) 
                    w3_f(1,i0a,iam) = w3_f(1,i0a,iam) + w3_f1(1,jsize,iam)  
                    w3_f(2,i0a,iam) = w3_f(2,i0a,iam) + w3_f1(2,jsize,iam)   
                    w3_f(3,i0a,iam) = w3_f(3,i0a,iam) + w3_f1(3,jsize,iam)  
                    w3_f(1,i0c,iam) = w3_f(1,i0c,iam) + w3_f2(1,jsize,iam)   
                    w3_f(2,i0c,iam) = w3_f(2,i0c,iam) + w3_f2(2,jsize,iam)   
                    w3_f(3,i0c,iam) = w3_f(3,i0c,iam) + w3_f2(3,jsize,iam)  
                    w3_f(1,i0d,iam) = w3_f(1,i0d,iam) + w3_f3(1,jsize,iam)   
                    w3_f(2,i0d,iam) = w3_f(2,i0d,iam) + w3_f3(2,jsize,iam)   
                    w3_f(3,i0d,iam) = w3_f(3,i0d,iam) + w3_f3(3,jsize,iam)   
                    fi0_2(1,i0) = fi0_2(1,i0) + fi0_3(1,jsize,k0_b)  
                    fi0_2(2,i0) = fi0_2(2,i0) + fi0_3(2,jsize,k0_b)  
                    fi0_2(3,i0) = fi0_2(3,i0) + fi0_3(3,jsize,k0_b)  
                enddo ! enddo jsize
                do k0=k0_s,k0_e 
                    do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
                        w3_f(1,i0,iam) = w3_f(1,i0,iam) + fi0_2(1,i0)
                        w3_f(2,i0,iam) = w3_f(2,i0,iam) + fi0_2(2,i0)
                        w3_f(3,i0,iam) = w3_f(3,i0,iam) + fi0_2(3,i0)
                    enddo 
                enddo 
        enddo ! k0, or k0_b
!$omp end do nowait
!$omp end parallel
   
    wk_p_energy = wk_p_energy + Udihedral * sfact 
    wk_vir2(1,0) = wk_vir2(1,0) + wk_v11 *sfact 
    wk_vir2(2,0) = wk_vir2(2,0) + wk_v22 *sfact
    wk_vir2(3,0) = wk_vir2(3,0) + wk_v33 *sfact
    wk_vir2(4,0) = wk_vir2(4,0) + wk_v21 *sfact
    wk_vir2(5,0) = wk_vir2(5,0) + wk_v31 *sfact
    wk_vir2(6,0) = wk_vir2(6,0) + wk_v32 *sfact
#ifdef DEBUGFCE
    oUd=0d0
    call mpi_allreduce(Udihedral*sfact,oUd,1, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,ipar)
#ifdef KCAL
    if(myrank==0)write(*,*)'Pot(dihed)=',oUd *kJ_mol/4.184d0,'[kcal/mol]'
#else
    if(myrank==0)write(*,*)'Pot(dihed)=',oUd *kJ_mol,'[kJ/mol]'
#endif
#endif

  end subroutine add_charmm_dihedral_a_no_bond_breaking

!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate energy and force of dihedral \n
!!         with events of bond-breaking.
!! \author Zhiye Tang
!<
  subroutine add_charmm_dihedral_a_bond_breaking()
!-----------------------------------------------------------------------
!
!   /* V = Kd(1 + cos[n.chi - delta])
!    *
!    *            D
!    *           /   chi = angle between plain ABC and plain BCD
!    *     B----C
!    *   /
!    *  A
!    *
!    */
!
!   /* aliases
!    * delta => d, chi => c */
!
    use trajectory_org
    use trajectory_mpi
    use atom_virial

    use forces
    use md_monitors
    use md_const
    use param
    use mpi_tool
    use md_condition, only : allow_bond_breaking
#ifdef ONEPROC_AXIS
    use boundary, only : pbc_pair
#endif
    use openmp_tool, only : nomp
    use comm_direct2_dr, only : comm_direct_2_dr
    implicit none
    integer(4) :: i, j, ip, ia, ib, i0a, i0b, i0c, i0d, i0, ipar, k0
    real(8) :: a, b, c, coef, tmp_angle
    real(8) :: tmp__x, tmp__y, tmp__z, i_norm
    real(8) :: tmp1__x, tmp1__y, tmp1__z, tmp2__x, tmp2__y, tmp2__z
    real(8) :: vBA__x, vBA__y, vBA__z, vBC__x, vBC__y, vBC__z
    real(8) :: vCD__x, vCD__y, vCD__z, vCA__x, vCA__y, vCA__z
    real(8) :: v1__x, v1__y, v1__z, v2__x, v2__y, v2__z
    real(8) :: v3__x, v3__y, v3__z
    real(8) :: Fa__x, Fa__y, Fa__z, Fb__x, Fb__y, Fb__z
    real(8) :: Fc__x, Fc__y, Fc__z, Fd__x, Fd__y, Fd__z
    real(8) :: Udihedral,oUd
    real(8) :: v11,v22,v33
    real(8) :: wk_v11, wk_v22, wk_v33
    real(8) :: sfact=1.0d0  ! scaling factor to potE and virial (default)
    real(8) :: wkvd(6), vd(6) ! debug
    real(8) :: v21,v31,v32
    real(8) :: wk_v21, wk_v31, wk_v32
    real(8) :: ss1,ss2,ss3,r1,r2,r3,c0,c1,c2,s1,s2,s12,s
    real(8) :: df,a11,a22,a33,a12,a13,a23,sx2,sy2,sz2
    real(8) :: vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb4x,vb4y,vb4z,vb1xvb2,vb3xvb2,vAx,vAy,vAz,vBx,vBy,vBz,rvA2,rvB2
#ifdef DIHEDRAL_TABLE
    integer(4) :: k
#endif
    real(8) :: fi0(3)
!
    integer(4) :: iam
    include 'mpif.h'

!default: MTD
    if( nselfseg .eq. 0 ) then
        k0_w=1
        k0_n=0
    else
        k0_w=(nselfseg-1)/nomp + 1  ! note: nselfseg is variable with MD cycle
        k0_n=(nselfseg-1)/k0_w + 1  ! note: nselfseg is variable with MD cycle
    endif

    call check_atombound_dihedral()

    iam = 0
    Udihedral = 0.0d0
    wk_v11=0d0;wk_v22=0d0;wk_v33=0d0
    wk_v21=0d0;wk_v31=0d0;wk_v32=0d0

!$omp parallel default(none) &
!$omp& private(v11,v22,v33) &
!$omp& reduction(+:wk_v11,wk_v22,wk_v33) &
!$omp& private(v21,v31,v32) &
!$omp& reduction(+:wk_v21,wk_v31,wk_v32) &
!$omp& private(k0,i,iam,j,i0,ipar,jsize) &
!$omp& private(fi0,tmp_angle) &
!$omp& private(ip,ia,ib,i0a,i0b,i0c,i0d,vBA__x,vBA__y,vBA__z) &
!$omp& private(vBC__x,vBC__y,vBC__z,vCD__x,vCD__y,vCD__z) &
!$omp& private(vCA__x,vCA__y,vCA__z,v1__x,v1__y,v1__z) &
!$omp& private(v2__x,v2__y,v2__z,c,tmp__x,tmp__y,tmp__z) &
!$omp& private(coef,v3__x,v3__y,v3__z,a,b,i_norm) &
!$omp& private(Fa__x,Fa__y,Fa__z,tmp1__x,tmp1__y,tmp1__z) &
!$omp& private(tmp2__x,tmp2__y,tmp2__z) &
!$omp& private(Fb__x,Fb__y,Fb__z,Fc__x,Fc__y,Fc__z) &
!$omp& private(Fd__x,Fd__y,Fd__z) &
!$omp& private(jmax_2,k0_b,k0_s,k0_e) &
!$omp& shared(i0a_2,i0c_2,i0d_2,i0_2,fi0_2,fi0_3) &
!$omp& shared(w3_f1,w3_f2,w3_f3,jlist,jmax) &
#ifdef DIHEDRAL_TABLE
!$omp& shared(m2i,i2m,nA_t,topA,atom2A_t,atom3A_t,atom4A_t,indexA_t, continue_flagA_t) &
!$omp& shared(wkxyz,w3_f,KdA_t,ndA_t,deltaA_t,wk_vir2,ndihedralB_t,ndihedralA_t) &
!$omp& shared(nB_t,topB,atom1B_t,atom3B_t,atom4B_t,indexB_t,KdB_t,continue_flagB_t,ndB_t,deltaB_t,paranum) &
#else
!$omp& shared(m2i,i2m,nA,topA,atom2A,atom3A,atom4A) &
!$omp& shared(wkxyz,w3_f,KdA,ndA,deltaA,wk_vir2) &
!$omp& shared(nB,topB,atom1B,atom3B,atom4B,KdB,ndB,deltaB,paranum) &
#endif
!$omp& shared(tag,na_per_cell) &
!$omp& shared(lsegtop,lseg_natoms,nselfseg) &
!$omp& shared(sfact,k0_n,k0_w) &
#ifdef DIHEDRAL_TABLE
!$omp& private(k) &
#endif
!$omp& private(ss1,ss2,ss3,r1,r2,r3,c0,c1,c2,s1,s2,s12,s,myrank) &
!$omp& private(df,a11,a22,a33,a12,a13,a23,sx2,sy2,sz2) &
!$omp& private(vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb4x,vb4y,vb4z,vb1xvb2,vb3xvb2,vAx,vAy,vAz,vBx,vBy,vBz,rvA2,rvB2) &
!$omp& reduction(+:Udihedral)
!$  iam = omp_get_thread_num()
!$omp do
!default: MTD
    do k0_b=1,k0_n
      k0_s=(k0_b-1)*k0_w + 1 
      k0_e=k0_s + k0_w - 1 
      k0_e=min(k0_e,nselfseg) 
     do k0=k0_s,k0_e 
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          ia   = m2i(i0)
          ipar = paranum(ia)
          fi0(1:3)=0d0
!default: MTD
          fi0_2(1:3,i0)=0d0
!default: MTD
       enddo
     enddo
     jsize=0 
     do k0=k0_s,k0_e 
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          ia   = m2i(i0)
          ipar = paranum(ia)
#ifdef DIHEDRAL_TABLE
          do j = 1, nB_t(ipar)
             if(ia < atom3B_t(ipar,j)+(ia-ipar))cycle
             jsize=jsize+1
             i0a = i2m(atom1B_t(ipar,j)+(ia-ipar))
             i0c = i2m(atom3B_t(ipar,j)+(ia-ipar))
             i0d = i2m(atom4B_t(ipar,j)+(ia-ipar))
             i0a_2(jsize,k0_b)= i0a
             i0c_2(jsize,k0_b)= i0c
             i0d_2(jsize,k0_b)= i0d
             i0_2(jsize,k0_b) = i0
             jlist(jsize,k0_b)=j
          enddo
#else /* DIHEDRAL_TABLE */
          do j = 1, nB(ipar)
             if(ia < atom3B(ipar,j)+(ia-ipar))cycle
             jsize=jsize+1
             i0a = i2m(atom1B(ipar,j)+(ia-ipar))
             i0c = i2m(atom3B(ipar,j)+(ia-ipar))
             i0d = i2m(atom4B(ipar,j)+(ia-ipar))
             i0a_2(jsize,k0_b)= i0a
             i0c_2(jsize,k0_b)= i0c
             i0d_2(jsize,k0_b)= i0d
             i0_2(jsize,k0_b) = i0
             jlist(jsize,k0_b)=j
          enddo
#endif /* DIHEDRAL_TABLE */
       enddo
     enddo
     jmax_2=jsize

!default: MTD
! in this branch w/ and w/o DIHEDRAL_TABLE it does the same thing
          do jsize=1,jmax_2 
             j=jlist(jsize,k0_b) 

!default: MTD
             i0a = i0a_2(jsize,k0_b) 
             i0c = i0c_2(jsize,k0_b) 
             i0d = i0d_2(jsize,k0_b) 
             i0  = i0_2(jsize,k0_b) 
             ia   = m2i(i0)
             ipar = paranum(ia)
#ifdef DIHEDRAL_TABLE
            fi0_3(1,jsize,k0_b) = 0.0d0
            fi0_3(2,jsize,k0_b) = 0.0d0
            fi0_3(3,jsize,k0_b) = 0.0d0
            w3_f1(1,jsize,iam) = 0.0d0
            w3_f1(2,jsize,iam) = 0.0d0
            w3_f1(3,jsize,iam) = 0.0d0
            w3_f2(1,jsize,iam) = 0.0d0
            w3_f2(2,jsize,iam) = 0.0d0
            w3_f2(3,jsize,iam) = 0.0d0
            w3_f3(1,jsize,iam) = 0.0d0
            w3_f3(2,jsize,iam) = 0.0d0
            w3_f3(3,jsize,iam) = 0.0d0
            if(indexB_t(ipar,j)<0) then
               cycle
            endif

            k = indexB_t(ipar,j)
            do while ( .true. )
#endif /* DIHEDRAL_TABLE */
             vb1x = wkxyz(1,i0a) - wkxyz(1,i0 )
             vb1y = wkxyz(2,i0a) - wkxyz(2,i0 )
             vb1z = wkxyz(3,i0a) - wkxyz(3,i0 )
             vb2x = wkxyz(1,i0c) - wkxyz(1,i0 )
             vb2y = wkxyz(2,i0c) - wkxyz(2,i0 )
             vb2z = wkxyz(3,i0c) - wkxyz(3,i0 )
             vb3x = wkxyz(1,i0d) - wkxyz(1,i0c)
             vb3y = wkxyz(2,i0d) - wkxyz(2,i0c)
             vb3z = wkxyz(3,i0d) - wkxyz(3,i0c)
             vb4x = wkxyz(1,i0a) - wkxyz(1,i0c)
             vb4y = wkxyz(2,i0a) - wkxyz(2,i0c)
             vb4z = wkxyz(3,i0a) - wkxyz(3,i0c)
#ifdef ONEPROC_AXIS
             call pbc_pair(vb1x,vb1y,vb1z)
             call pbc_pair(vb2x,vb2y,vb2z)
             call pbc_pair(vb3x,vb3y,vb3z)
             call pbc_pair(vb4x,vb4y,vb4z)
#endif
            ss1 = 1.0d0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z)
            ss2 = 1.0d0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z)
            ss3 = 1.0d0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z)

            r1 = sqrt(ss1);
            r2 = sqrt(ss2);
            r3 = sqrt(ss3);

            vb1xvb2 = vb1x * vb2x + vb1y * vb2y + vb1z * vb2z
            vb3xvb2 = vb3x * vb2x + vb3y * vb2y + vb3z * vb2z

            c0 = (vb1x * vb3x + vb1y * vb3y + vb1z * vb3z) * r1 * r3
            c1 = vb1xvb2 * r1 * r2
            c2 = -vb3xvb2 * r3 * r2

            s1 = 1.0d0 - c1*c1;
            s1 = max(s1, 0.001d0)
            s1 = 1.0d0 / s1

            s2 = 1.0d0 - c2*c2;
            s2 = max(s2, 0.001d0)
            s2 = 1.0d0 / s2

            s12 = sqrt(s1*s2);
            c = (c1*c2 + c0) * s12;
            c=min(+1d0,c)
            c=max(-1d0,c)

            ! vb1 cross vb3 * vb2 > 0 => +
            ! vb1 cross vb3 * vb2 <=0 => -
            c = acos(c)
            tmp1__x = vb1y * vb3z - vb1z * vb3y
            tmp1__y = vb1z * vb3x - vb1x * vb3z
            tmp1__z = vb1x * vb3y - vb1y * vb3x
            c= sign(c, tmp1__x * vb2x + tmp1__y * vb2y + tmp1__z * vb2z)
            ! at this point, c means the dihedral angle theta

#ifdef DIHEDRAL_TABLE
             !         force coefficient
             tmp_angle = ndB_t(k)*c - deltaB_t(k)
             coef = KdB_t(k) * ndB_t(k) * sin(tmp_angle)
#else /* DIHEDRAL_TABLE */
             !         force coefficient
             tmp_angle = ndB(ipar,j)*c - deltaB(ipar,j)
             coef = KdB(ipar,j) * ndB(ipar,j) * sin(tmp_angle)
#endif /* DIHEDRAL_TABLE */

             !         force coefficient

            vAx = -vb1y * vb2z + vb1z * vb2y;
            vAy = -vb1z * vb2x + vb1x * vb2z;
            vAz = -vb1x * vb2y + vb1y * vb2x;
            vBx = -vb3y * vb2z + vb3z * vb2y;
            vBy = -vb3z * vb2x + vb3x * vb2z;
            vBz = -vb3x * vb2y + vb3y * vb2x;
            rvA2 = 1.0e0 / (vAx * vAx + vAy * vAy + vAz * vAz);
            rvB2 = 1.0e0 / (vBx * vBx + vBy * vBy + vBz * vBz);

            a11 = (-coef / r2 * rvA2 );
            a22 = (-vb1xvb2 * rvA2 * r2) * coef;
            a33 = ( vb3xvb2 * rvB2 * r2) * coef;
            a13 = (coef / r2 * rvB2);

!default: MTD
            Fa__x = vAx * a11;
            Fa__y = vAy * a11;
            Fa__z = vAz * a11;
            Fb__x = - Fa__x + vAx * a22 + vBx * a33;
            Fb__y = - Fa__y + vAy * a22 + vBy * a33;
            Fb__z = - Fa__z + vAz * a22 + vBz * a33;

             !         calculating & adding forces
#ifndef DIHEDRAL_TABLE
             fi0_3(1,jsize,k0_b)=Fb__x
             fi0_3(2,jsize,k0_b)=Fb__y
             fi0_3(3,jsize,k0_b)=Fb__z
#else
             fi0_3(1,jsize,k0_b)=fi0_3(1,jsize,k0_b)+Fb__x
             fi0_3(2,jsize,k0_b)=fi0_3(2,jsize,k0_b)+Fb__y
             fi0_3(3,jsize,k0_b)=fi0_3(3,jsize,k0_b)+Fb__z
#endif
            Fd__x = vBx * a13;
            Fd__y = vBy * a13;
            Fd__z = vBz * a13;

            Fc__x = - Fd__x - vAx * a22 - vBx * a33;
            Fc__y = - Fd__y - vAy * a22 - vBy * a33;
            Fc__z = - Fd__z - vAz * a22 - vBz * a33;
            
             v11 = (vb4x*Fa__x - vb2x*Fb__x + vb3x*Fd__x)
             v22 = (vb4y*Fa__y - vb2y*Fb__y + vb3y*Fd__y)
             v33 = (vb4z*Fa__z - vb2z*Fb__z + vb3z*Fd__z)
             v21 = (vb4y*Fa__x - vb2y*Fb__x + vb3y*Fd__x)
             v31 = (vb4z*Fa__x - vb2z*Fb__x + vb3z*Fd__x)
             v32 = (vb4z*Fa__y - vb2z*Fb__y + vb3z*Fd__y)
             wk_v11 = wk_v11 + v11
             wk_v22 = wk_v22 + v22
             wk_v33 = wk_v33 + v33
             wk_v21 = wk_v21 + v21
             wk_v31 = wk_v31 + v31
             wk_v32 = wk_v32 + v32
#ifndef DIHEDRAL_TABLE
                        w3_f1(1,jsize,iam) = Fa__x
                        w3_f1(2,jsize,iam) = Fa__y
                        w3_f1(3,jsize,iam) = Fa__z
                        w3_f2(1,jsize,iam) = Fc__x
                        w3_f2(2,jsize,iam) = Fc__y
                        w3_f2(3,jsize,iam) = Fc__z
                        w3_f3(1,jsize,iam) = Fd__x
                        w3_f3(2,jsize,iam) = Fd__y
                        w3_f3(3,jsize,iam) = Fd__z
#else
                        w3_f1(1,jsize,iam) = w3_f1(1,jsize,iam) + Fa__x
                        w3_f1(2,jsize,iam) = w3_f1(2,jsize,iam) + Fa__y
                        w3_f1(3,jsize,iam) = w3_f1(3,jsize,iam) + Fa__z
                        w3_f2(1,jsize,iam) = w3_f2(1,jsize,iam) + Fc__x
                        w3_f2(2,jsize,iam) = w3_f2(2,jsize,iam) + Fc__y
                        w3_f2(3,jsize,iam) = w3_f2(3,jsize,iam) + Fc__z
                        w3_f3(1,jsize,iam) = w3_f3(1,jsize,iam) + Fd__x
                        w3_f3(2,jsize,iam) = w3_f3(2,jsize,iam) + Fd__y
                        w3_f3(3,jsize,iam) = w3_f3(3,jsize,iam) + Fd__z
#endif
#ifdef DIHEDRAL_TABLE
             Udihedral=Udihedral+KdB_t(k)*(1d0+cos(tmp_angle))
#else /* DIHEDRAL_TABLE */
             Udihedral=Udihedral+KdB(ipar,j)*(1d0+cos(tmp_angle))
#endif /* DIHEDRAL_TABLE */
             !         Potential and virial are calculated in loop A.
#ifdef DIHEDRAL_TABLE
             if( .not. continue_flagB_t(k) ) exit
             k=k+1
             if( k > ndihedralB_t ) exit ! exit loop if not continue the same group
            enddo
#endif /* DIHEDRAL_TABLE */
          enddo
!default: MTD
          do jsize=1,jmax_2 
             i0a = i0a_2(jsize,k0_b) 
             i0c = i0c_2(jsize,k0_b) 
             i0d = i0d_2(jsize,k0_b) 
             i0  = i0_2(jsize,k0_b) 
             w3_f(1,i0a,iam) = w3_f(1,i0a,iam) + w3_f1(1,jsize,iam)  
             w3_f(2,i0a,iam) = w3_f(2,i0a,iam) + w3_f1(2,jsize,iam)   
             w3_f(3,i0a,iam) = w3_f(3,i0a,iam) + w3_f1(3,jsize,iam)  
             w3_f(1,i0c,iam) = w3_f(1,i0c,iam) + w3_f2(1,jsize,iam)   
             w3_f(2,i0c,iam) = w3_f(2,i0c,iam) + w3_f2(2,jsize,iam)   
             w3_f(3,i0c,iam) = w3_f(3,i0c,iam) + w3_f2(3,jsize,iam)  
             w3_f(1,i0d,iam) = w3_f(1,i0d,iam) + w3_f3(1,jsize,iam)   
             w3_f(2,i0d,iam) = w3_f(2,i0d,iam) + w3_f3(2,jsize,iam)   
             w3_f(3,i0d,iam) = w3_f(3,i0d,iam) + w3_f3(3,jsize,iam)   
             fi0_2(1,i0) = fi0_2(1,i0) + fi0_3(1,jsize,k0_b)  
             fi0_2(2,i0) = fi0_2(2,i0) + fi0_3(2,jsize,k0_b)  
             fi0_2(3,i0) = fi0_2(3,i0) + fi0_3(3,jsize,k0_b)  
          enddo 
     do k0=k0_s,k0_e 
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          w3_f(1,i0,iam) = w3_f(1,i0,iam) + fi0_2(1,i0)
          w3_f(2,i0,iam) = w3_f(2,i0,iam) + fi0_2(2,i0)
          w3_f(3,i0,iam) = w3_f(3,i0,iam) + fi0_2(3,i0)
       enddo 
     enddo 
    enddo ! k0
!$omp end do nowait
!$omp end parallel

    wk_p_energy = wk_p_energy + Udihedral * sfact 
    wk_vir2(1,0) = wk_vir2(1,0) + wk_v11 *sfact 
    wk_vir2(2,0) = wk_vir2(2,0) + wk_v22 *sfact 
    wk_vir2(3,0) = wk_vir2(3,0) + wk_v33 *sfact 
    wk_vir2(4,0) = wk_vir2(4,0) + wk_v21 *sfact 
    wk_vir2(5,0) = wk_vir2(5,0) + wk_v31 *sfact 
    wk_vir2(6,0) = wk_vir2(6,0) + wk_v32 *sfact 
#ifdef DEBUGFCE
    oUd=0d0
    call mpi_allreduce(Udihedral*sfact,oUd,1, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,ipar)
#ifdef KCAL
    if(myrank==0)write(*,*)'Pot(dihed)=',oUd *kJ_mol/4.184d0,'[kcal/mol]'
#else
    if(myrank==0)write(*,*)'Pot(dihed)=',oUd *kJ_mol,'[kJ/mol]'
#endif
#endif

  end subroutine add_charmm_dihedral_a_bond_breaking

! ======================================================================================================
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate energy and force of dihedral.
!! \author Kensuke Iwahashi, Zhiye Tang
!<
  subroutine add_charmm_dihedral_a_no_bond_breaking_stable()
!-----------------------------------------------------------------------
!
!   /* V = Kd(1 + cos[n.chi - delta])
!    *
!    *            D
!    *           /   chi = angle between plain ABC and plain BCD
!    *     B----C
!    *   /
!    *  A
!    *
!    */
!
!   /* aliases
!    * delta => d, chi => c */
!
    use trajectory_org
    use trajectory_mpi
    use atom_virial

    use forces
    use md_monitors
    use md_const
    use param
    use mpi_tool
    use md_condition, only : allow_bond_breaking
#ifdef ONEPROC_AXIS
    use boundary, only : pbc_pair
#endif
    use openmp_tool, only : nomp
    use comm_direct2_dr, only : comm_direct_2_dr
    implicit none
    integer(4) :: i, j, ip, ia, ib, i0a, i0b, i0c, i0d, i0, ipar, k0
    real(8) :: a, b, c, coef, tmp_angle
    real(8) :: tmp__x, tmp__y, tmp__z, i_norm
    real(8) :: tmp1__x, tmp1__y, tmp1__z, tmp2__x, tmp2__y, tmp2__z
    real(8) :: vBA__x, vBA__y, vBA__z, vBC__x, vBC__y, vBC__z
    real(8) :: vCD__x, vCD__y, vCD__z, vCA__x, vCA__y, vCA__z
    real(8) :: v1__x, v1__y, v1__z, v2__x, v2__y, v2__z
    real(8) :: v3__x, v3__y, v3__z
    real(8) :: Fa__x, Fa__y, Fa__z, Fb__x, Fb__y, Fb__z
    real(8) :: Fc__x, Fc__y, Fc__z, Fd__x, Fd__y, Fd__z
    real(8) :: Udihedral,oUd
    real(8) :: v11,v22,v33
    real(8) :: wk_v11, wk_v22, wk_v33
    real(8) :: sfact=1.0d0  ! scaling factor to potE and virial (default)
    real(8) :: wkvd(6), vd(6) ! debug
    real(8) :: v21,v31,v32
    real(8) :: wk_v21, wk_v31, wk_v32
    real(8) :: ss1,ss2,ss3,r1,r2,r3,c0,c1,c2,s1,s2,s12,s
    real(8) :: df,a11,a22,a33,a12,a13,a23,sx2,sy2,sz2
    real(8) :: vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb4x,vb4y,vb4z
#ifdef DIHEDRAL_TABLE
    integer(4) :: k
#endif
    real(8) :: fi0(3)
!
    integer(4) :: iam
    include 'mpif.h'

!default: MTD
    if( nselfseg .eq. 0 ) then
        k0_w=1
        k0_n=0
    else
        k0_w=(nselfseg-1)/nomp + 1  ! note: nselfseg is variable with MD cycle
        k0_n=(nselfseg-1)/k0_w + 1  ! note: nselfseg is variable with MD cycle
    endif

    call check_atombound_dihedral()

    iam = 0
    Udihedral = 0.0d0
    wk_v11=0d0;wk_v22=0d0;wk_v33=0d0
    wk_v21=0d0;wk_v31=0d0;wk_v32=0d0

!$omp parallel default(none) &
!$omp& private(v11,v22,v33) &
!$omp& reduction(+:wk_v11,wk_v22,wk_v33) &
!$omp& private(v21,v31,v32) &
!$omp& reduction(+:wk_v21,wk_v31,wk_v32) &
!$omp& private(k0,i,iam,j,i0,ipar,jsize) &
!$omp& private(fi0,tmp_angle) &
!$omp& private(ip,ia,ib,i0a,i0b,i0c,i0d,vBA__x,vBA__y,vBA__z) &
!$omp& private(vBC__x,vBC__y,vBC__z,vCD__x,vCD__y,vCD__z) &
!$omp& private(vCA__x,vCA__y,vCA__z,v1__x,v1__y,v1__z) &
!$omp& private(v2__x,v2__y,v2__z,c,tmp__x,tmp__y,tmp__z) &
!$omp& private(coef,v3__x,v3__y,v3__z,a,b,i_norm) &
!$omp& private(Fa__x,Fa__y,Fa__z,tmp1__x,tmp1__y,tmp1__z) &
!$omp& private(tmp2__x,tmp2__y,tmp2__z) &
!$omp& private(Fb__x,Fb__y,Fb__z,Fc__x,Fc__y,Fc__z) &
!$omp& private(Fd__x,Fd__y,Fd__z) &
!$omp& private(jmax_2,k0_b,k0_s,k0_e) &
!$omp& shared(i0a_2,i0c_2,i0d_2,i0_2,fi0_2,fi0_3) &
!$omp& shared(w3_f1,w3_f2,w3_f3,jlist,jmax) &
#ifdef DIHEDRAL_TABLE
!$omp& shared(m2i,i2m,nA_t,topA,atom2A_t,atom3A_t,atom4A_t,indexA_t, continue_flagA_t) &
!$omp& shared(wkxyz,w3_f,KdA_t,ndA_t,deltaA_t,wk_vir2,ndihedralB_t,ndihedralA_t) &
!$omp& shared(nB_t,topB,atom1B_t,atom3B_t,atom4B_t,indexB_t,KdB_t,continue_flagB_t,ndB_t,deltaB_t,paranum) &
#else
!$omp& shared(m2i,i2m,nA,topA,atom2A,atom3A,atom4A) &
!$omp& shared(wkxyz,w3_f,KdA,ndA,deltaA,wk_vir2) &
!$omp& shared(nB,topB,atom1B,atom3B,atom4B,KdB,ndB,deltaB,paranum) &
#endif
!$omp& shared(tag,na_per_cell) &
!$omp& shared(lsegtop,lseg_natoms,nselfseg) &
!$omp& shared(sfact,k0_n,k0_w) &
#ifdef DIHEDRAL_TABLE
!$omp& private(k) &
#endif
!$omp& shared(myrank) &
!$omp& private(ss1,ss2,ss3,r1,r2,r3,c0,c1,c2,s1,s2,s12,s) &
!$omp& private(df,a11,a22,a33,a12,a13,a23,sx2,sy2,sz2) &
!$omp& private(vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb4x,vb4y,vb4z) &
!$omp& reduction(+:Udihedral)
!$  iam = omp_get_thread_num()
!$omp do
!default: MTD
    do k0_b=1,k0_n
        k0_s=(k0_b-1)*k0_w + 1 
        k0_e=k0_s + k0_w - 1 
        k0_e=min(k0_e,nselfseg) 
        do k0=k0_s,k0_e 
            do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
                ia   = m2i(i0)
                ipar = paranum(ia)
                fi0(1:3)=0d0
!default: MTD
                fi0_2(1:3,i0)=0d0
!default: MTD
            enddo ! enddo i0
        enddo ! endo k0
        jsize=0 
        do k0=k0_s,k0_e 
            do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
                ia   = m2i(i0)
                ipar = paranum(ia)
#ifdef DIHEDRAL_TABLE
                do j = 1, nB_t(ipar)
                    if(ia < atom3B_t(ipar,j)+(ia-ipar)) cycle
                    jsize=jsize+1
                    i0a = i2m(atom1B_t(ipar,j)+(ia-ipar))
                    i0c = i2m(atom3B_t(ipar,j)+(ia-ipar))
                    i0d = i2m(atom4B_t(ipar,j)+(ia-ipar))
                    i0a_2(jsize,k0_b)= i0a
                    i0c_2(jsize,k0_b)= i0c
                    i0d_2(jsize,k0_b)= i0d
                    i0_2(jsize,k0_b) = i0
                    jlist(jsize,k0_b)=j
                enddo
#else /* DIHEDRAL_TABLE */
                do j = 1, nB(ipar)
                    if(ia < atom3B(ipar,j)+(ia-ipar)) cycle
                    jsize=jsize+1
                    i0a = i2m(atom1B(ipar,j)+(ia-ipar))
                    i0c = i2m(atom3B(ipar,j)+(ia-ipar))
                    i0d = i2m(atom4B(ipar,j)+(ia-ipar))
                    i0a_2(jsize,k0_b)= i0a
                    i0c_2(jsize,k0_b)= i0c
                    i0d_2(jsize,k0_b)= i0d
                    i0_2(jsize,k0_b) = i0
                    jlist(jsize,k0_b)=j
                enddo
#endif /* DIHEDRAL_TABLE */
            enddo ! enddo i0
        enddo ! enddo k0
        jmax_2=jsize

!default: MTD
! in this branch w/ and w/o DIHEDRAL_TABLE it does the same thing
                do jsize=1,jmax_2 
                    j=jlist(jsize,k0_b) 

!default: MTD
                    i0a = i0a_2(jsize,k0_b) 
                    i0c = i0c_2(jsize,k0_b) 
                    i0d = i0d_2(jsize,k0_b) 
                    i0  = i0_2(jsize,k0_b) 
                    ia   = m2i(i0)
                    ipar = paranum(ia)
#ifdef DIHEDRAL_TABLE
                    k = indexB_t(ipar,j)
                     fi0_3(1,jsize,k0_b) = 0.0d0
                     fi0_3(2,jsize,k0_b) = 0.0d0
                     fi0_3(3,jsize,k0_b) = 0.0d0
                     w3_f1(1,jsize,iam) = 0.0d0
                     w3_f1(2,jsize,iam) = 0.0d0
                     w3_f1(3,jsize,iam) = 0.0d0
                     w3_f2(1,jsize,iam) = 0.0d0
                     w3_f2(2,jsize,iam) = 0.0d0
                     w3_f2(3,jsize,iam) = 0.0d0
                     w3_f3(1,jsize,iam) = 0.0d0
                     w3_f3(2,jsize,iam) = 0.0d0
                     w3_f3(3,jsize,iam) = 0.0d0
                     
                    do while ( .true. )
#endif /* DIHEDRAL_TABLE */
                        vb1x = wkxyz(1,i0a) - wkxyz(1,i0 )
                        vb1y = wkxyz(2,i0a) - wkxyz(2,i0 )
                        vb1z = wkxyz(3,i0a) - wkxyz(3,i0 )
                        vb2x = wkxyz(1,i0c) - wkxyz(1,i0 )
                        vb2y = wkxyz(2,i0c) - wkxyz(2,i0 )
                        vb2z = wkxyz(3,i0c) - wkxyz(3,i0 )
                        vb3x = wkxyz(1,i0d) - wkxyz(1,i0c)
                        vb3y = wkxyz(2,i0d) - wkxyz(2,i0c)
                        vb3z = wkxyz(3,i0d) - wkxyz(3,i0c)
                        vb4x = wkxyz(1,i0a) - wkxyz(1,i0c)
                        vb4y = wkxyz(2,i0a) - wkxyz(2,i0c)
                        vb4z = wkxyz(3,i0a) - wkxyz(3,i0c)
             
#ifdef ONEPROC_AXIS

                        call pbc_pair(vb1x,vb1y,vb1z)
                        call pbc_pair(vb2x,vb2y,vb2z)
                        call pbc_pair(vb3x,vb3y,vb3z)
                        call pbc_pair(vb4x,vb4y,vb4z)
#endif

                        ss1 = 1.0d0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z)
                        ss2 = 1.0d0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z)
                        ss3 = 1.0d0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z)

                        r1 = sqrt(ss1);
                        r2 = sqrt(ss2);
                        r3 = sqrt(ss3);

                        c0 = (vb1x * vb3x + vb1y * vb3y + vb1z * vb3z) * r1 * r3;
                        c1 = (vb1x * vb2x + vb1y * vb2y + vb1z * vb2z) * r1 * r2;
                        c2 = -(vb3x * vb2x + vb3y * vb2y + vb3z * vb2z) * r3 * r2;

                        s1 = 1.0d0 - c1*c1;
                        s1 = max(s1, 0.001d0)
                        s1 = 1.0d0 / s1

                        s2 = 1.0d0 - c2*c2;
                        s2 = max(s2, 0.001d0)
                        s2 = 1.0d0 / s2

                        s12 = sqrt(s1*s2);
                        c = (c1*c2 + c0) * s12;
                        c=min(+1d0,c)
                        c=max(-1d0,c)

                        s = sqrt(1.0d0 - c*c)
                        s = max(s, 0.001d0)

#ifdef DIHEDRAL_TABLE
                        !         force coefficient
                        tmp_angle = ndB_t(k)*acos(c) - deltaB_t(k)
                        coef = KdB_t(k) * ndB_t(k) * sin(tmp_angle)
#else /* DIHEDRAL_TABLE */
                        !         force coefficient
                        tmp_angle = ndB(ipar,j)*acos(c) - deltaB(ipar,j)
                        coef = KdB(ipar,j) * ndB(ipar,j) * sin(tmp_angle)
#endif /* DIHEDRAL_TABLE */

                        coef = coef / s;
                        c = c * coef;
                        s12 = s12 * coef;
                        a11 = c*ss1*s1;
                        a12 = -r1*r2*(c1*c*s1 + c2*s12);
                        a13 = -r1*r3*s12;
                        a22 = -ss2 * (2.0*c0*s12 - c*(s1+s2));
                        a33 = c*ss3*s2;
                        a23 = r2*r3*(c2*c*s2 + c1*s12);
                        sx2  = a22*vb2x + a23*vb3x + a12*vb1x;
                        sy2  = a22*vb2y + a23*vb3y + a12*vb1y;
                        sz2  = a22*vb2z + a23*vb3z + a12*vb1z;


!default: MTD
                        Fa__x = a12*vb2x + a13*vb3x + a11*vb1x
                        Fa__y = a12*vb2y + a13*vb3y + a11*vb1y
                        Fa__z = a12*vb2z + a13*vb3z + a11*vb1z
                        Fb__x = -sx2 - Fa__x;
                        Fb__y = -sy2 - Fa__y;
                        Fb__z = -sz2 - Fa__z;

             !         calculating & adding forces
#ifndef DIHEDRAL_TABLE
                        fi0_3(1,jsize,k0_b)=Fb__x
                        fi0_3(2,jsize,k0_b)=Fb__y
                        fi0_3(3,jsize,k0_b)=Fb__z
#else
                        fi0_3(1,jsize,k0_b)=fi0_3(1,jsize,k0_b)+Fb__x
                        fi0_3(2,jsize,k0_b)=fi0_3(2,jsize,k0_b)+Fb__y
                        fi0_3(3,jsize,k0_b)=fi0_3(3,jsize,k0_b)+Fb__z
#endif
                        Fd__x = a23*vb2x + a33*vb3x + a13*vb1x;
                        Fd__y = a23*vb2y + a33*vb3y + a13*vb1y;
                        Fd__z = a23*vb2z + a33*vb3z + a13*vb1z;

                        Fc__x = sx2 - Fd__x;
                        Fc__y = sy2 - Fd__y;
                        Fc__z = sz2 - Fd__z;

                        v11 = (vb4x*Fa__x - vb2x*Fb__x + vb3x*Fd__x)
                        v22 = (vb4y*Fa__y - vb2y*Fb__y + vb3y*Fd__y)
                        v33 = (vb4z*Fa__z - vb2z*Fb__z + vb3z*Fd__z)
                        v21 = (vb4y*Fa__x - vb2y*Fb__x + vb3y*Fd__x)
                        v31 = (vb4z*Fa__x - vb2z*Fb__x + vb3z*Fd__x)
                        v32 = (vb4z*Fa__y - vb2z*Fb__y + vb3z*Fd__y)
                        wk_v11 = wk_v11 + v11
                        wk_v22 = wk_v22 + v22
                        wk_v33 = wk_v33 + v33
                        wk_v21 = wk_v21 + v21
                        wk_v31 = wk_v31 + v31
                        wk_v32 = wk_v32 + v32
#ifndef DIHEDRAL_TABLE
                        w3_f1(1,jsize,iam) = Fa__x
                        w3_f1(2,jsize,iam) = Fa__y
                        w3_f1(3,jsize,iam) = Fa__z
                        w3_f2(1,jsize,iam) = Fc__x
                        w3_f2(2,jsize,iam) = Fc__y
                        w3_f2(3,jsize,iam) = Fc__z
                        w3_f3(1,jsize,iam) = Fd__x
                        w3_f3(2,jsize,iam) = Fd__y
                        w3_f3(3,jsize,iam) = Fd__z
#else
                        w3_f1(1,jsize,iam) = w3_f1(1,jsize,iam) + Fa__x
                        w3_f1(2,jsize,iam) = w3_f1(2,jsize,iam) + Fa__y
                        w3_f1(3,jsize,iam) = w3_f1(3,jsize,iam) + Fa__z
                        w3_f2(1,jsize,iam) = w3_f2(1,jsize,iam) + Fc__x
                        w3_f2(2,jsize,iam) = w3_f2(2,jsize,iam) + Fc__y
                        w3_f2(3,jsize,iam) = w3_f2(3,jsize,iam) + Fc__z
                        w3_f3(1,jsize,iam) = w3_f3(1,jsize,iam) + Fd__x
                        w3_f3(2,jsize,iam) = w3_f3(2,jsize,iam) + Fd__y
                        w3_f3(3,jsize,iam) = w3_f3(3,jsize,iam) + Fd__z
#endif

#ifdef DIHEDRAL_TABLE
                        Udihedral=Udihedral+KdB_t(k)*(1d0+cos(tmp_angle))
#else /* DIHEDRAL_TABLE */
                        Udihedral=Udihedral+KdB(ipar,j)*(1d0+cos(tmp_angle))
#endif /* DIHEDRAL_TABLE */
             !         Potential and virial are calculated in loop A.
#ifdef DIHEDRAL_TABLE
                        if( .not. continue_flagB_t(k) ) exit
                        k=k+1
                        if( k > ndihedralB_t ) exit ! exit loop if not continue the same group
                    enddo ! enddo while (true)
#endif /* DIHEDRAL_TABLE */
                enddo ! enddo j
!default: MTD
                do jsize=1,jmax_2 
                    i0a = i0a_2(jsize,k0_b) 
                    i0c = i0c_2(jsize,k0_b) 
                    i0d = i0d_2(jsize,k0_b) 
                    i0  = i0_2(jsize,k0_b) 
                    w3_f(1,i0a,iam) = w3_f(1,i0a,iam) + w3_f1(1,jsize,iam)  
                    w3_f(2,i0a,iam) = w3_f(2,i0a,iam) + w3_f1(2,jsize,iam)   
                    w3_f(3,i0a,iam) = w3_f(3,i0a,iam) + w3_f1(3,jsize,iam)  
                    w3_f(1,i0c,iam) = w3_f(1,i0c,iam) + w3_f2(1,jsize,iam)   
                    w3_f(2,i0c,iam) = w3_f(2,i0c,iam) + w3_f2(2,jsize,iam)   
                    w3_f(3,i0c,iam) = w3_f(3,i0c,iam) + w3_f2(3,jsize,iam)  
                    w3_f(1,i0d,iam) = w3_f(1,i0d,iam) + w3_f3(1,jsize,iam)   
                    w3_f(2,i0d,iam) = w3_f(2,i0d,iam) + w3_f3(2,jsize,iam)   
                    w3_f(3,i0d,iam) = w3_f(3,i0d,iam) + w3_f3(3,jsize,iam)   
                    fi0_2(1,i0) = fi0_2(1,i0) + fi0_3(1,jsize,k0_b)  
                    fi0_2(2,i0) = fi0_2(2,i0) + fi0_3(2,jsize,k0_b)  
                    fi0_2(3,i0) = fi0_2(3,i0) + fi0_3(3,jsize,k0_b)  
                enddo ! enddo jsize
                do k0=k0_s,k0_e 
                    do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
                        w3_f(1,i0,iam) = w3_f(1,i0,iam) + fi0_2(1,i0)
                        w3_f(2,i0,iam) = w3_f(2,i0,iam) + fi0_2(2,i0)
                        w3_f(3,i0,iam) = w3_f(3,i0,iam) + fi0_2(3,i0)
                    enddo 
                enddo 
        enddo ! k0, or k0_b
!$omp end do nowait
!$omp end parallel
   
    wk_p_energy = wk_p_energy + Udihedral * sfact 
    wk_vir2(1,0) = wk_vir2(1,0) + wk_v11 *sfact 
    wk_vir2(2,0) = wk_vir2(2,0) + wk_v22 *sfact
    wk_vir2(3,0) = wk_vir2(3,0) + wk_v33 *sfact
    wk_vir2(4,0) = wk_vir2(4,0) + wk_v21 *sfact
    wk_vir2(5,0) = wk_vir2(5,0) + wk_v31 *sfact
    wk_vir2(6,0) = wk_vir2(6,0) + wk_v32 *sfact
#ifdef DEBUGFCE
    oUd=0d0
    call mpi_allreduce(Udihedral*sfact,oUd,1, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,ipar)
#ifdef KCAL
    if(myrank==0)write(*,*)'Pot(dihed)=',oUd *kJ_mol/4.184d0,'[kcal/mol]'
#else
    if(myrank==0)write(*,*)'Pot(dihed)=',oUd *kJ_mol,'[kJ/mol]'
#endif
#endif

  end subroutine add_charmm_dihedral_a_no_bond_breaking_stable

!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate energy and force of dihedral \n
!!         with events of bond-breaking.
!! \author Zhiye Tang
!<
  subroutine add_charmm_dihedral_a_bond_breaking_stable()
!-----------------------------------------------------------------------
!
!   /* V = Kd(1 + cos[n.chi - delta])
!    *
!    *            D
!    *           /   chi = angle between plain ABC and plain BCD
!    *     B----C
!    *   /
!    *  A
!    *
!    */
!
!   /* aliases
!    * delta => d, chi => c */
!
    use trajectory_org
    use trajectory_mpi
    use atom_virial

    use forces
    use md_monitors
    use md_const
    use param
    use mpi_tool
    use md_condition, only : allow_bond_breaking
#ifdef ONEPROC_AXIS
    use boundary, only : pbc_pair
#endif
    use openmp_tool, only : nomp
    use comm_direct2_dr, only : comm_direct_2_dr
    implicit none
    integer(4) :: i, j, ip, ia, ib, i0a, i0b, i0c, i0d, i0, ipar, k0
    real(8) :: a, b, c, coef, tmp_angle
    real(8) :: tmp__x, tmp__y, tmp__z, i_norm
    real(8) :: tmp1__x, tmp1__y, tmp1__z, tmp2__x, tmp2__y, tmp2__z
    real(8) :: vBA__x, vBA__y, vBA__z, vBC__x, vBC__y, vBC__z
    real(8) :: vCD__x, vCD__y, vCD__z, vCA__x, vCA__y, vCA__z
    real(8) :: v1__x, v1__y, v1__z, v2__x, v2__y, v2__z
    real(8) :: v3__x, v3__y, v3__z
    real(8) :: Fa__x, Fa__y, Fa__z, Fb__x, Fb__y, Fb__z
    real(8) :: Fc__x, Fc__y, Fc__z, Fd__x, Fd__y, Fd__z
    real(8) :: Udihedral,oUd
    real(8) :: v11,v22,v33
    real(8) :: wk_v11, wk_v22, wk_v33
    real(8) :: sfact=1.0d0  ! scaling factor to potE and virial (default)
    real(8) :: wkvd(6), vd(6) ! debug
    real(8) :: v21,v31,v32
    real(8) :: wk_v21, wk_v31, wk_v32
    real(8) :: ss1,ss2,ss3,r1,r2,r3,c0,c1,c2,s1,s2,s12,s
    real(8) :: df,a11,a22,a33,a12,a13,a23,sx2,sy2,sz2
    real(8) :: vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb4x,vb4y,vb4z
#ifdef DIHEDRAL_TABLE
    integer(4) :: k
#endif
    real(8) :: fi0(3)
!
    integer(4) :: iam
    include 'mpif.h'

!default: MTD
    if( nselfseg .eq. 0 ) then
        k0_w=1
        k0_n=0
    else
        k0_w=(nselfseg-1)/nomp + 1  ! note: nselfseg is variable with MD cycle
        k0_n=(nselfseg-1)/k0_w + 1  ! note: nselfseg is variable with MD cycle
    endif

    call check_atombound_dihedral()

    iam = 0
    Udihedral = 0.0d0
    wk_v11=0d0;wk_v22=0d0;wk_v33=0d0
    wk_v21=0d0;wk_v31=0d0;wk_v32=0d0

!$omp parallel default(none) &
!$omp& private(v11,v22,v33) &
!$omp& reduction(+:wk_v11,wk_v22,wk_v33) &
!$omp& private(v21,v31,v32) &
!$omp& reduction(+:wk_v21,wk_v31,wk_v32) &
!$omp& private(k0,i,iam,j,i0,ipar,jsize) &
!$omp& private(fi0,tmp_angle) &
!$omp& private(ip,ia,ib,i0a,i0b,i0c,i0d,vBA__x,vBA__y,vBA__z) &
!$omp& private(vBC__x,vBC__y,vBC__z,vCD__x,vCD__y,vCD__z) &
!$omp& private(vCA__x,vCA__y,vCA__z,v1__x,v1__y,v1__z) &
!$omp& private(v2__x,v2__y,v2__z,c,tmp__x,tmp__y,tmp__z) &
!$omp& private(coef,v3__x,v3__y,v3__z,a,b,i_norm) &
!$omp& private(Fa__x,Fa__y,Fa__z,tmp1__x,tmp1__y,tmp1__z) &
!$omp& private(tmp2__x,tmp2__y,tmp2__z) &
!$omp& private(Fb__x,Fb__y,Fb__z,Fc__x,Fc__y,Fc__z) &
!$omp& private(Fd__x,Fd__y,Fd__z) &
!$omp& private(jmax_2,k0_b,k0_s,k0_e) &
!$omp& shared(i0a_2,i0c_2,i0d_2,i0_2,fi0_2,fi0_3) &
!$omp& shared(w3_f1,w3_f2,w3_f3,jlist,jmax) &
#ifdef DIHEDRAL_TABLE
!$omp& shared(m2i,i2m,nA_t,topA,atom2A_t,atom3A_t,atom4A_t,indexA_t, continue_flagA_t) &
!$omp& shared(wkxyz,w3_f,KdA_t,ndA_t,deltaA_t,wk_vir2,ndihedralB_t,ndihedralA_t) &
!$omp& shared(nB_t,topB,atom1B_t,atom3B_t,atom4B_t,indexB_t,KdB_t,continue_flagB_t,ndB_t,deltaB_t,paranum) &
#else
!$omp& shared(m2i,i2m,nA,topA,atom2A,atom3A,atom4A) &
!$omp& shared(wkxyz,w3_f,KdA,ndA,deltaA,wk_vir2) &
!$omp& shared(nB,topB,atom1B,atom3B,atom4B,KdB,ndB,deltaB,paranum) &
#endif
!$omp& shared(tag,na_per_cell) &
!$omp& shared(lsegtop,lseg_natoms,nselfseg) &
!$omp& shared(sfact,k0_n,k0_w) &
#ifdef DIHEDRAL_TABLE
!$omp& private(k) &
#endif
!$omp& private(ss1,ss2,ss3,r1,r2,r3,c0,c1,c2,s1,s2,s12,s,myrank) &
!$omp& private(df,a11,a22,a33,a12,a13,a23,sx2,sy2,sz2) &
!$omp& private(vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb4x,vb4y,vb4z) &
!$omp& reduction(+:Udihedral)
!$  iam = omp_get_thread_num()
!$omp do
!default: MTD
    do k0_b=1,k0_n
      k0_s=(k0_b-1)*k0_w + 1 
      k0_e=k0_s + k0_w - 1 
      k0_e=min(k0_e,nselfseg) 
     do k0=k0_s,k0_e 
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          ia   = m2i(i0)
          ipar = paranum(ia)
          fi0(1:3)=0d0
!default: MTD
          fi0_2(1:3,i0)=0d0
!default: MTD
       enddo
     enddo
     jsize=0 
     do k0=k0_s,k0_e 
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          ia   = m2i(i0)
          ipar = paranum(ia)
#ifdef DIHEDRAL_TABLE
          do j = 1, nB_t(ipar)
             if(ia < atom3B_t(ipar,j)+(ia-ipar))cycle
             jsize=jsize+1
             i0a = i2m(atom1B_t(ipar,j)+(ia-ipar))
             i0c = i2m(atom3B_t(ipar,j)+(ia-ipar))
             i0d = i2m(atom4B_t(ipar,j)+(ia-ipar))
             i0a_2(jsize,k0_b)= i0a
             i0c_2(jsize,k0_b)= i0c
             i0d_2(jsize,k0_b)= i0d
             i0_2(jsize,k0_b) = i0
             jlist(jsize,k0_b)=j
          enddo
#else /* DIHEDRAL_TABLE */
          do j = 1, nB(ipar)
             if(ia < atom3B(ipar,j)+(ia-ipar))cycle
             jsize=jsize+1
             i0a = i2m(atom1B(ipar,j)+(ia-ipar))
             i0c = i2m(atom3B(ipar,j)+(ia-ipar))
             i0d = i2m(atom4B(ipar,j)+(ia-ipar))
             i0a_2(jsize,k0_b)= i0a
             i0c_2(jsize,k0_b)= i0c
             i0d_2(jsize,k0_b)= i0d
             i0_2(jsize,k0_b) = i0
             jlist(jsize,k0_b)=j
          enddo
#endif /* DIHEDRAL_TABLE */
       enddo
     enddo
     jmax_2=jsize

!default: MTD
! in this branch w/ and w/o DIHEDRAL_TABLE it does the same thing
          do jsize=1,jmax_2 
             j=jlist(jsize,k0_b) 

!default: MTD
             i0a = i0a_2(jsize,k0_b) 
             i0c = i0c_2(jsize,k0_b) 
             i0d = i0d_2(jsize,k0_b) 
             i0  = i0_2(jsize,k0_b) 
             ia   = m2i(i0)
             ipar = paranum(ia)
#ifdef DIHEDRAL_TABLE
            fi0_3(1,jsize,k0_b) = 0.0d0
            fi0_3(2,jsize,k0_b) = 0.0d0
            fi0_3(3,jsize,k0_b) = 0.0d0
            w3_f1(1,jsize,iam) = 0.0d0
            w3_f1(2,jsize,iam) = 0.0d0
            w3_f1(3,jsize,iam) = 0.0d0
            w3_f2(1,jsize,iam) = 0.0d0
            w3_f2(2,jsize,iam) = 0.0d0
            w3_f2(3,jsize,iam) = 0.0d0
            w3_f3(1,jsize,iam) = 0.0d0
            w3_f3(2,jsize,iam) = 0.0d0
            w3_f3(3,jsize,iam) = 0.0d0
            if(indexB_t(ipar,j)<0) then
               cycle
            endif

            k = indexB_t(ipar,j)
            do while ( .true. )
#endif /* DIHEDRAL_TABLE */
             vb1x = wkxyz(1,i0a) - wkxyz(1,i0 )
             vb1y = wkxyz(2,i0a) - wkxyz(2,i0 )
             vb1z = wkxyz(3,i0a) - wkxyz(3,i0 )
             vb2x = wkxyz(1,i0c) - wkxyz(1,i0 )
             vb2y = wkxyz(2,i0c) - wkxyz(2,i0 )
             vb2z = wkxyz(3,i0c) - wkxyz(3,i0 )
             vb3x = wkxyz(1,i0d) - wkxyz(1,i0c)
             vb3y = wkxyz(2,i0d) - wkxyz(2,i0c)
             vb3z = wkxyz(3,i0d) - wkxyz(3,i0c)
             vb4x = wkxyz(1,i0a) - wkxyz(1,i0c)
             vb4y = wkxyz(2,i0a) - wkxyz(2,i0c)
             vb4z = wkxyz(3,i0a) - wkxyz(3,i0c)
#ifdef ONEPROC_AXIS
             call pbc_pair(vb1x,vb1y,vb1z)
             call pbc_pair(vb2x,vb2y,vb2z)
             call pbc_pair(vb3x,vb3y,vb3z)
             call pbc_pair(vb4x,vb4y,vb4z)
#endif
            ss1 = 1.0d0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z)
            ss2 = 1.0d0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z)
            ss3 = 1.0d0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z)

            r1 = sqrt(ss1);
            r2 = sqrt(ss2);
            r3 = sqrt(ss3);

            c0 = (vb1x * vb3x + vb1y * vb3y + vb1z * vb3z) * r1 * r3;
            c1 = (vb1x * vb2x + vb1y * vb2y + vb1z * vb2z) * r1 * r2;
            c2 = -(vb3x * vb2x + vb3y * vb2y + vb3z * vb2z) * r3 * r2;

            s1 = 1.0d0 - c1*c1;
            s1 = max(s1, 0.001d0)
            s1 = 1.0d0 / s1

            s2 = 1.0d0 - c2*c2;
            s2 = max(s2, 0.001d0)
            s2 = 1.0d0 / s2

            s12 = sqrt(s1*s2);
            c = (c1*c2 + c0) * s12;
            c=min(+1d0,c)
            c=max(-1d0,c)

            s = sqrt(1.0d0 - c*c)
            s = max(s, 0.001d0)

#ifdef DIHEDRAL_TABLE
             !         force coefficient
             tmp_angle = ndB_t(k)*acos(c) - deltaB_t(k)
             coef = KdB_t(k) * ndB_t(k) * sin(tmp_angle)
#else /* DIHEDRAL_TABLE */
             !         force coefficient
             tmp_angle = ndB(ipar,j)*acos(c) - deltaB(ipar,j)
             coef = KdB(ipar,j) * ndB(ipar,j) * sin(tmp_angle)
#endif /* DIHEDRAL_TABLE */

             !         force coefficient

            coef = coef / s;
            c = c * coef;
            s12 = s12 * coef;
            a11 = c*ss1*s1;
            a12 = -r1*r2*(c1*c*s1 + c2*s12);
            a13 = -r1*r3*s12;
            a22 = -ss2 * (2.0*c0*s12 - c*(s1+s2));
            a33 = c*ss3*s2;
            a23 = r2*r3*(c2*c*s2 + c1*s12);
            sx2  = a22*vb2x + a23*vb3x + a12*vb1x;
            sy2  = a22*vb2y + a23*vb3y + a12*vb1y;
            sz2  = a22*vb2z + a23*vb3z + a12*vb1z;

!default: MTD
            Fa__x = a12*vb2x + a13*vb3x + a11*vb1x
            Fa__y = a12*vb2y + a13*vb3y + a11*vb1y
            Fa__z = a12*vb2z + a13*vb3z + a11*vb1z
            Fb__x = -sx2 - Fa__x;
            Fb__y = -sy2 - Fa__y;
            Fb__z = -sz2 - Fa__z;

             !         calculating & adding forces
#ifndef DIHEDRAL_TABLE
             fi0_3(1,jsize,k0_b)=Fb__x
             fi0_3(2,jsize,k0_b)=Fb__y
             fi0_3(3,jsize,k0_b)=Fb__z
#else
             fi0_3(1,jsize,k0_b)=fi0_3(1,jsize,k0_b)+Fb__x
             fi0_3(2,jsize,k0_b)=fi0_3(2,jsize,k0_b)+Fb__y
             fi0_3(3,jsize,k0_b)=fi0_3(3,jsize,k0_b)+Fb__z
#endif
            Fd__x = a23*vb2x + a33*vb3x + a13*vb1x;
            Fd__y = a23*vb2y + a33*vb3y + a13*vb1y;
            Fd__z = a23*vb2z + a33*vb3z + a13*vb1z;

            Fc__x = sx2 - Fd__x;
            Fc__y = sy2 - Fd__y;
            Fc__z = sz2 - Fd__z;
            
             v11 = (vb4x*Fa__x - vb2x*Fb__x + vb3x*Fd__x)
             v22 = (vb4y*Fa__y - vb2y*Fb__y + vb3y*Fd__y)
             v33 = (vb4z*Fa__z - vb2z*Fb__z + vb3z*Fd__z)
             v21 = (vb4y*Fa__x - vb2y*Fb__x + vb3y*Fd__x)
             v31 = (vb4z*Fa__x - vb2z*Fb__x + vb3z*Fd__x)
             v32 = (vb4z*Fa__y - vb2z*Fb__y + vb3z*Fd__y)
             wk_v11 = wk_v11 + v11
             wk_v22 = wk_v22 + v22
             wk_v33 = wk_v33 + v33
             wk_v21 = wk_v21 + v21
             wk_v31 = wk_v31 + v31
             wk_v32 = wk_v32 + v32
#ifndef DIHEDRAL_TABLE
                        w3_f1(1,jsize,iam) = Fa__x
                        w3_f1(2,jsize,iam) = Fa__y
                        w3_f1(3,jsize,iam) = Fa__z
                        w3_f2(1,jsize,iam) = Fc__x
                        w3_f2(2,jsize,iam) = Fc__y
                        w3_f2(3,jsize,iam) = Fc__z
                        w3_f3(1,jsize,iam) = Fd__x
                        w3_f3(2,jsize,iam) = Fd__y
                        w3_f3(3,jsize,iam) = Fd__z
#else
                        w3_f1(1,jsize,iam) = w3_f1(1,jsize,iam) + Fa__x
                        w3_f1(2,jsize,iam) = w3_f1(2,jsize,iam) + Fa__y
                        w3_f1(3,jsize,iam) = w3_f1(3,jsize,iam) + Fa__z
                        w3_f2(1,jsize,iam) = w3_f2(1,jsize,iam) + Fc__x
                        w3_f2(2,jsize,iam) = w3_f2(2,jsize,iam) + Fc__y
                        w3_f2(3,jsize,iam) = w3_f2(3,jsize,iam) + Fc__z
                        w3_f3(1,jsize,iam) = w3_f3(1,jsize,iam) + Fd__x
                        w3_f3(2,jsize,iam) = w3_f3(2,jsize,iam) + Fd__y
                        w3_f3(3,jsize,iam) = w3_f3(3,jsize,iam) + Fd__z
#endif
#ifdef DIHEDRAL_TABLE
             Udihedral=Udihedral+KdB_t(k)*(1d0+cos(tmp_angle))
#else /* DIHEDRAL_TABLE */
             Udihedral=Udihedral+KdB(ipar,j)*(1d0+cos(tmp_angle))
#endif /* DIHEDRAL_TABLE */
             !         Potential and virial are calculated in loop A.
#ifdef DIHEDRAL_TABLE
             if( .not. continue_flagB_t(k) ) exit
             k=k+1
             if( k > ndihedralB_t ) exit ! exit loop if not continue the same group
            enddo
#endif /* DIHEDRAL_TABLE */
          enddo
!default: MTD
          do jsize=1,jmax_2 
             i0a = i0a_2(jsize,k0_b) 
             i0c = i0c_2(jsize,k0_b) 
             i0d = i0d_2(jsize,k0_b) 
             i0  = i0_2(jsize,k0_b) 
             w3_f(1,i0a,iam) = w3_f(1,i0a,iam) + w3_f1(1,jsize,iam)  
             w3_f(2,i0a,iam) = w3_f(2,i0a,iam) + w3_f1(2,jsize,iam)   
             w3_f(3,i0a,iam) = w3_f(3,i0a,iam) + w3_f1(3,jsize,iam)  
             w3_f(1,i0c,iam) = w3_f(1,i0c,iam) + w3_f2(1,jsize,iam)   
             w3_f(2,i0c,iam) = w3_f(2,i0c,iam) + w3_f2(2,jsize,iam)   
             w3_f(3,i0c,iam) = w3_f(3,i0c,iam) + w3_f2(3,jsize,iam)  
             w3_f(1,i0d,iam) = w3_f(1,i0d,iam) + w3_f3(1,jsize,iam)   
             w3_f(2,i0d,iam) = w3_f(2,i0d,iam) + w3_f3(2,jsize,iam)   
             w3_f(3,i0d,iam) = w3_f(3,i0d,iam) + w3_f3(3,jsize,iam)   
             fi0_2(1,i0) = fi0_2(1,i0) + fi0_3(1,jsize,k0_b)  
             fi0_2(2,i0) = fi0_2(2,i0) + fi0_3(2,jsize,k0_b)  
             fi0_2(3,i0) = fi0_2(3,i0) + fi0_3(3,jsize,k0_b)  
          enddo 
     do k0=k0_s,k0_e 
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          w3_f(1,i0,iam) = w3_f(1,i0,iam) + fi0_2(1,i0)
          w3_f(2,i0,iam) = w3_f(2,i0,iam) + fi0_2(2,i0)
          w3_f(3,i0,iam) = w3_f(3,i0,iam) + fi0_2(3,i0)
       enddo 
     enddo 
    enddo ! k0
!$omp end do nowait
!$omp end parallel

    wk_p_energy = wk_p_energy + Udihedral * sfact 
    wk_vir2(1,0) = wk_vir2(1,0) + wk_v11 *sfact 
    wk_vir2(2,0) = wk_vir2(2,0) + wk_v22 *sfact 
    wk_vir2(3,0) = wk_vir2(3,0) + wk_v33 *sfact 
    wk_vir2(4,0) = wk_vir2(4,0) + wk_v21 *sfact 
    wk_vir2(5,0) = wk_vir2(5,0) + wk_v31 *sfact 
    wk_vir2(6,0) = wk_vir2(6,0) + wk_v32 *sfact 
#ifdef DEBUGFCE
    oUd=0d0
    call mpi_allreduce(Udihedral*sfact,oUd,1, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,ipar)
#ifdef KCAL
    if(myrank==0)write(*,*)'Pot(dihed)=',oUd *kJ_mol/4.184d0,'[kcal/mol]'
#else
    if(myrank==0)write(*,*)'Pot(dihed)=',oUd *kJ_mol,'[kJ/mol]'
#endif
#endif

  end subroutine add_charmm_dihedral_a_bond_breaking_stable


! ======================================================================================================

#ifdef DIHEDRAL_TABLE
  function search_value_group() result(group)
    implicit none
    integer(4) :: group, i, j, ibegin, iend

    do group = 1, size(v_group_begin) - 1
       ibegin = v_group_begin(group)
       iend   = v_group_begin(group+1) - 1
       if( is_same_value_group(             &
            v_value    %items(ibegin:iend), &
            v_value_tmp%items(1:v_value_tmp%size))) return
    end do
    group = 0

  end function search_value_group

  subroutine cat_dihedral_value_tmp()
    implicit none
    integer(4) :: group, i

    if( v_value_tmp%size == 0 ) then
       return
    end if

    group = search_value_group()

    if( group == 0 ) then
       group = size(v_group_begin)
       do i = 1, v_value_tmp%size
          call v_value%push_back(v_value_tmp%items(i))
       end do
       v_group_begin = [v_group_begin, v_value%size + 1]
    end if

    call v_value_tmp%resize(0)
    v_index%items(v_index%size)%group = group

  end subroutine cat_dihedral_value_tmp

  subroutine add_to_dihedral_table(atom1, atom2, atom3, atom4, Kd, nd, delta)
    use precision, only : wfp
    implicit none
    integer(4), intent(in) :: atom1, atom2, atom3, atom4
    real(kind=wfp)   , intent(in) :: Kd, nd, delta
    type(dihedral_table_index) :: index
    type(dihedral_table_value) :: val
    integer(4) :: group, i
#ifdef DIHEDRAL_TABLE_DEBUG
    type(dihedral_raw_data) :: raw
    type(dihedral_raw_data), allocatable :: raw_data_tmp(:)
    integer(4) :: nraw_data
#endif

    index%atom1 = atom1
    index%atom2 = atom2
    index%atom3 = atom3
    index%atom4 = atom4
    index%group = 0

    val%Kd    = Kd
    val%nd    = nd
    val%delta = delta

#ifdef DIHEDRAL_TABLE_DEBUG
    raw%atom1 = atom1
    raw%atom2 = atom2
    raw%atom3 = atom3
    raw%atom4 = atom4
    raw%Kd    = Kd
    raw%nd    = nd
    raw%delta = delta
    if( allocated(raw_data) ) then
       nraw_data = size(raw_data)
       allocate(raw_data_tmp(nraw_data))
       raw_data_tmp(:) = raw_data(:)
       deallocate(raw_data)
       allocate(raw_data(nraw_data+1))
       raw_data(:) = raw_data_tmp(:)
       raw_data(nraw_data+1) = raw
    else
       allocate(raw_data(1))
       raw_data(1) = raw
    end if
#endif /*DIHEDRAL_TABLE_DEBUG*/

    if( 0 < v_index%size ) then
       if( is_same_index_atoms(index, v_index%items(v_index%size)) ) then
          if( Kd /= 0.0 ) then
             call v_value_tmp%push_back(val)
          end if
          return
       end if

       call cat_dihedral_value_tmp()
    end if

    call v_index%push_back(index)
    if( Kd /= 0.0 ) then
       call v_value_tmp%push_back(val)
    end if

  end subroutine add_to_dihedral_table

  subroutine init_dihedral_table()
    use mod_dihedral_table_index
    use mod_dihedral_table_value
    use param
    implicit none

    call v_index%init()
    call v_value%init()
    call v_value_tmp%init()
    v_group_begin = (/ 1 /)

  end subroutine init_dihedral_table

  subroutine finalize_dihedral_table()
    use parse_dihedral
    implicit none
    integer(4) :: i, j, k, group
#ifdef DIHEDRAL_TABLE_DEBUG
    integer(4) :: nraw_data, nraw_data1, nsize
    type(dihedral_raw_data), allocatable :: raw_data1(:)
    type(dihedral_raw_data) :: raw
#endif

    call cat_dihedral_value_tmp()

#ifdef DIHEDRAL_TABLE_DEBUG
    if( myrank == 0 ) then

       ! for debug check consistency
       nraw_data = size(raw_data)
       write(*,*) "nraw_data       ", nraw_data

#if 0
       write(*,*) "raw_data"
       do i = 1, nraw_data
          write(*,'(5i5,3e15.5)') i, raw_data(i)
       end do
#endif

       ! remove Kd = 0 items
       j = 1
       do i = 1, nraw_data
          if( raw_data(i)%Kd /= 0.0d0 ) then
             raw_data(j) = raw_data(i)
             j = j + 1
          end if
       end do
       nraw_data = j - 1
       write(*,*) "nraw_data(Kd/=0)", nraw_data

#if 0
       write(*,*) "raw_data(Kd/=0)"
       do i = 1, nraw_data
          write(*,'(5i5,3e15.5)') i, raw_data(i)
       end do
#endif

       ! reconstruct raw data from the tables
       nraw_data1 = 0
       do i = 1, v_index%size
          group = v_index%items(i)%group
          if( group == 0 ) cycle
          nraw_data1 = nraw_data1 + (v_group_begin(group+1) - v_group_begin(group))
       end do

       write(*,*) "nraw_data1      ", nraw_data1

       allocate(raw_data1(nraw_data1))
       k = 1
       do i = 1, v_index%size
          group = v_index%items(i)%group
          if( group == 0 ) cycle
          raw%atom1 = v_index%items(i)%atom1
          raw%atom2 = v_index%items(i)%atom2
          raw%atom3 = v_index%items(i)%atom3
          raw%atom4 = v_index%items(i)%atom4
          do j = v_group_begin(group), v_group_begin(group+1) - 1
             raw%Kd    = v_value%items(j)%Kd
             raw%nd    = v_value%items(j)%nd
             raw%delta = v_value%items(j)%delta
             raw_data1(k) = raw
             k = k + 1
          end do
       end do

#if 0
       write(*,*) "nraw_data1"
       do i = 1, nraw_data1
          write(*,'(5i5,3e15.5)') i, raw_data1(i)
       end do
#endif

       if( nraw_data /= nraw_data1 ) then
          write(*,*) "ERROR: nrawdata /= nraw_data1"
          stop
       end if

       do i = 1, nraw_data
          if(  raw_data(i)%atom1 /= raw_data1(i)%atom1 .or. &
               raw_data(i)%atom2 /= raw_data1(i)%atom2 .or. &
               raw_data(i)%atom3 /= raw_data1(i)%atom3 .or. &
               raw_data(i)%atom4 /= raw_data1(i)%atom4 .or. &
               raw_data(i)%Kd    /= raw_data1(i)%Kd    .or. &
               raw_data(i)%nd    /= raw_data1(i)%nd    .or. &
               raw_data(i)%delta /= raw_data1(i)%delta ) then
             write(*,*) "ERROR: ", i, " th parameters are inconsistent."
          end if
       end do
       write(*,*) "Dihedral Table Parse Check Passed."
    end if
#endif /*DIHEDRAL_TABLE_DEBUG*/

    call set_md_charmm_dihedralA_table()
    call set_md_charmm_dihedralB_table()

    call v_index%init()
    call v_value%init()
    call v_value_tmp%init()
    if( allocated(v_group_begin) ) deallocate(v_group_begin)

#ifdef DIHEDRAL_TABLE_DEBUG
    if( myrank == 0 ) then
       nsize = &
            (size(atom2A) + size(atom3A) + size(atom4A)) * 4   + &
            (size(KdA)    + size(ndA)    + size(deltaA)) * wfp + &
            (size(atom1B) + size(atom3B) + size(atom4B)) * 4   + &
            (size(KdB)    + size(ndB)    + size(deltaB)) * wfp + &
            (size(nA)     + size(nB)                   ) * 4   + &
            (size(nArm)   + size(atom4Arm)             ) * 4
       write(*,*) "memory size(original) ", nsize, " Bytes"
       nsize = &
            (size(atom2A_t) + size(atom3A_t) + size(atom4A_t)) * 4   + &
            (size(KdA_t)    + size(ndA_t)    + size(deltaA_t)) * wfp + &
            (size(atom1B_t) + size(atom3B_t) + size(atom4B_t)) * 4   + &
            (size(KdB_t)    + size(ndB_t)    + size(deltaB_t)) * wfp + &
            (size(nA_t)     + size(nB_t)                     ) * 4   + &
            (size(nArm)     + size(atom4Arm)                 ) * 4
       write(*,*) "memory size(table)    ", nsize, " Bytes"
    end if
#endif
  end subroutine finalize_dihedral_table

  subroutine set_md_charmm_dihedralA_table
    use param
    implicit none
    integer(4) :: i,j,k,maxv,group,atom1,atom2,atom3,atom4,index,jA_tmp,jj,jA
    real(kind=wfp) :: Kd, nd, delta
#ifdef DIHEDRAL_TABLE_DEBUG
    integer(4), allocatable :: atom2A1(:,:), atom3A1(:,:), atom4A1(:,:), nA1(:)
    real(kind=wfp), allocatable :: KdA1(:,:), ndA1(:,:), deltaA1(:,:)
    integer(4), allocatable :: atom2A2(:,:), atom3A2(:,:), atom4A2(:,:), nA2(:)
    real(kind=wfp), allocatable :: KdA2(:,:), ndA2(:,:), deltaA2(:,:)
    integer(4), allocatable :: nArm_(:), atom4Arm_(:,:)
#endif
    
    allocate(nA_t(npara))
    nA_t(:) = 0

    do i = 1, v_index%size
       atom1 = v_index%items(i)%atom1
       atom4 = v_index%items(i)%atom4
       nA_t(atom1) = nA_t(atom1) + 1
       nA_t(atom4) = nA_t(atom4) + 1
    end do
    
    maxv = maxval(nA_t)

    allocate(atom2A_t(npara,maxv))
    allocate(atom3A_t(npara,maxv))
    allocate(atom4A_t(npara,maxv))
    allocate(indexA_t(npara,maxv))
    allocate(KdA_t           (v_value%size))
    allocate(ndA_t           (v_value%size))
    allocate(deltaA_t        (v_value%size))
    allocate(continue_flagA_t(v_value%size))
    ndihedralA_t = v_value%size
    ndihedralAs = v_index%size

    nA_t(:) = 0

    do i = 1, v_index%size
       atom1 = v_index%items(i)%atom1
       atom2 = v_index%items(i)%atom2
       atom3 = v_index%items(i)%atom3
       atom4 = v_index%items(i)%atom4
       group = v_index%items(i)%group

       k = nA_t(atom1) + 1
       atom2A_t(atom1,k) = atom2
       atom3A_t(atom1,k) = atom3
       atom4A_t(atom1,k) = atom4
       if( group == 0 ) then
          indexA_t(atom1,k) = 0
       else
          indexA_t(atom1,k) = v_group_begin(group)
       end if
       nA_t(atom1) = k 

       k = nA_t(atom4) + 1
       atom2A_t(atom4,k) = atom3
       atom3A_t(atom4,k) = atom2
       atom4A_t(atom4,k) = atom1
       if( group == 0 ) then
          indexA_t(atom4,k) = 0
       else
          indexA_t(atom4,k) = v_group_begin(group)
       end if
       nA_t(atom4) = k
    end do

    KdA_t   (1:v_value%size) = v_value%items(1:v_value%size)%Kd
    ndA_t   (1:v_value%size) = v_value%items(1:v_value%size)%nd
    deltaA_t(1:v_value%size) = v_value%items(1:v_value%size)%delta

    continue_flagA_t(:) = .true.
    do i = 2, size(v_group_begin)
       continue_flagA_t(v_group_begin(i)-1) = .false.
    end do

!! set  nArm for scaling 1-4 in OPLS/AMBER
!   if(md_condition__force_field == OPLSAA .or. &
!        &   md_condition__force_field == AMBER  )then
       if( allocated(atom4Arm) ) then
#ifdef DIHEDRAL_TABLE_DEBUG
          nArm_     = nArm
          atom4Arm_ = atom4Arm
#endif
          deallocate(nArm)
          deallocate(atom4Arm)
          deallocate(a4A_tmp)
       end if
       allocate(nArm(npara)) !; write(*,*) 'nArm allocated,',npara
       allocate(atom4Arm(npara,maxv))
       allocate(a4A_tmp(maxv))
       nArm=0; atom4Arm=0 ! initialize to 0
       jA_tmp=0; a4A_tmp=0

       do i=1, npara
          jA_tmp=0
          do j=1,nA_t(i)
             if(jA_tmp==0)then
                jA_tmp=1
                a4A_tmp(jA_tmp) = atom4A_t(i,j)
             else
                do jj=1,jA_tmp
                   if(a4A_tmp(jj)==atom4A_t(i,j))goto 9
                enddo
                jA_tmp=jA_tmp+1 !! not overlapped
                a4A_tmp(jA_tmp) = atom4A_t(i,j)
9               continue
             endif
          enddo ! j
!!            if(myrank==0) write(*,*) 'i,jA_tmp=',i-1,jA_tmp
! check cyclic motief
          jA=0
          do j=1,jA_tmp
          do jj=1,nA_t(i)
              if(a4A_tmp(j) == atom3A_t(i,jj))then  !! 5-ring case
!!                if(myrank==0) write(*,*) '4 and 3 overlap', a4A_tmp(j)-1
                  goto 8
              else
                  continue
              endif
          enddo ! jj
          jA=jA+1
          atom4Arm(i,jA)=a4A_tmp(j)
8         continue
          enddo ! j
!
          nArm(i)=jA
          !debug
#ifdef OPLS_DEBUG
          if(myrank==0) write(*,*) '##checking 1-4 pair list(debug)###'
          if(myrank==0) write(*,*) 'i,nArm(i)=',i-1,nArm(i) !! i-1 correspoinds to id in .mdff
          do j=1,jA
             if(myrank==0)write(*,*) i-1,atom4Arm(i,j)-1
          enddo
#endif
          !debug
       enddo ! i
!   endif

#ifdef DIHEDRAL_TABLE_DEBUG
    if( myrank == 0 ) then
       atom2A1 = atom2A
       atom3A1 = atom3A
       atom4A1 = atom4A
       KdA1    = KdA
       ndA1    = ndA
       deltaA1 = deltaA
       nA1     = nA

       atom2A1 = 0
       atom3A1 = 0
       atom4A1 = 0
       KdA1    = 0.0d0
       ndA1    = 0.0d0
       deltaA1 = 0.0d0
       nA1     = 0

       do atom1 = 1, npara
          do i = 1, nA_t(atom1)
             atom2 = atom2A_t(atom1,i)
             atom3 = atom3A_t(atom1,i)
             atom4 = atom4A_t(atom1,i)
             index = indexA_t(atom1,i)
             if( index == 0 ) cycle
             do while(.true.)
                Kd    = KdA_t   (index)
                nd    = ndA_t   (index)
                delta = deltaA_t(index)
                k = nA1(atom1) + 1
                atom2A1(atom1,k) = atom2
                atom3A1(atom1,k) = atom3
                atom4A1(atom1,k) = atom4
                KdA1   (atom1,k) = Kd
                ndA1   (atom1,k) = nd
                deltaA1(atom1,k) = delta
                nA1    (atom1)   = k
                if( .not. continue_flagA_t(index) ) exit
                index = index + 1
             end do
          end do
       end do

       ! Remove Kd==0
       atom2A2 = atom2A
       atom3A2 = atom3A
       atom4A2 = atom4A
       KdA2    = KdA
       ndA2    = ndA
       deltaA2 = deltaA
       nA2     = nA
       do i = 1, npara
          k = 0
          do j = 1, nA2(i)
             if (KdA2(i,j) == 0.0d0) then
                cycle
             end if
             k = k + 1
             atom2A2(i,k) = atom2A2(i,j)
             atom3A2(i,k) = atom3A2(i,j)
             atom4A2(i,k) = atom4A2(i,j)
             KdA2   (i,k) = KdA2   (i,j)
             ndA2   (i,k) = ndA2   (i,j)
             deltaA2(i,k) = deltaA2(i,j)
          end do
          nA2(i) = k
       end do

       do i = 1, npara
          if( nA1(i) /= nA2(i)) then
             write(*,*) "ERROR: ", i, " th parameters of dihedralA are inconsistent."
             stop
          end if
          do j = 1, nA1(i)
             if(  atom2A1(i,j) /= atom2A2(i,j) .or. &
                  atom3A1(i,j) /= atom3A2(i,j) .or. &
                  atom4A1(i,j) /= atom4A2(i,j) .or. &
                  KdA1   (i,j) /= KdA2   (i,j) .or. &
                  ndA1   (i,j) /= ndA2   (i,j) .or. &
                  deltaA1(i,j) /= deltaA2(i,j) ) then
                write(*,*) "ERROR: (", i, j, ") th parameters of dihedralA are inconsistent."
                write(*,*) atom2A1(i,j), atom2A2(i,j)
                write(*,*) atom3A1(i,j), atom3A2(i,j)
                write(*,*) atom4A1(i,j), atom4A2(i,j)
                write(*,*) KdA1   (i,j), KdA2   (i,j)
                write(*,*) ndA1   (i,j), ndA2   (i,j)
                write(*,*) deltaA1(i,j), deltaA2(i,j)
                stop
             end if
          end do
       end do
       write(*,*) "Reproduce dihedralA Check Passed."

       do i = 1, npara
          if( nArm(i) /= nArm_(i)) then
             write(*,*) "ERROR: ", i, " th nArm is inconsistent."
             stop
          end if
          do j = 1, nArm(i)
             if( atom4Arm(i,j) /= atom4Arm_(i,j) ) then
                write(*,*) "ERROR: (", i, j, ") th atom4Arm is inconsistent."
                write(*,*) atom4Arm(i,j), atom4Arm_(i,j)
                stop
             end if
          end do
       end do
       write(*,*) "Reproduce atom4Arm Check Passed."

    endif
    deallocate(nArm_)
    deallocate(atom4Arm_)
#endif /*DIHEDRAL_TABLE_DEBUG*/

  end subroutine set_md_charmm_dihedralA_table

  subroutine set_md_charmm_dihedralB_table
    use param
    implicit none
    integer(4) :: i,j,k,maxv,group,atom1,atom2,atom3,atom4,index
    real(kind=wfp) :: Kd, nd, delta
#ifdef DIHEDRAL_TABLE_DEBUG
    integer(4), allocatable :: atom1B1(:,:), atom3B1(:,:), atom4B1(:,:), nB1(:)
    real(kind=wfp), allocatable :: KdB1(:,:), ndB1(:,:), deltaB1(:,:)
    integer(4), allocatable :: atom1B2(:,:), atom3B2(:,:), atom4B2(:,:), nB2(:)
    real(kind=wfp), allocatable :: KdB2(:,:), ndB2(:,:), deltaB2(:,:)
#endif
    
    allocate(nB_t(npara))
    nB_t(:) = 0

    do i = 1, v_index%size
       atom2 = v_index%items(i)%atom2
       atom3 = v_index%items(i)%atom3
       nB_t(atom2) = nB_t(atom2) + 1
       nB_t(atom3) = nB_t(atom3) + 1
    end do
    
    maxv = maxval(nB_t)

    allocate(atom1B_t(npara,maxv))
    allocate(atom3B_t(npara,maxv))
    allocate(atom4B_t(npara,maxv))
    allocate(indexB_t(npara,maxv))
    allocate(KdB_t           (v_value%size))
    allocate(ndB_t           (v_value%size))
    allocate(deltaB_t        (v_value%size))
    allocate(continue_flagB_t(v_value%size))
    ndihedralB_t = v_value%size
    ndihedralBs = v_index%size

    nB_t(:) = 0

    do i = 1, v_index%size
       atom1 = v_index%items(i)%atom1
       atom2 = v_index%items(i)%atom2
       atom3 = v_index%items(i)%atom3
       atom4 = v_index%items(i)%atom4
       group = v_index%items(i)%group

       k = nB_t(atom2) + 1
       atom1B_t(atom2,k) = atom1
       atom3B_t(atom2,k) = atom3
       atom4B_t(atom2,k) = atom4
       if( group == 0 ) then
          indexB_t(atom2,k) = 0
       else
          indexB_t(atom2,k) = v_group_begin(group)
       end if
       nB_t(atom2) = k 

       k = nB_t(atom3) + 1
       atom1B_t(atom3,k) = atom4
       atom3B_t(atom3,k) = atom2
       atom4B_t(atom3,k) = atom1
       if( group == 0 ) then
          indexB_t(atom3,k) = 0
       else
          indexB_t(atom3,k) = v_group_begin(group)
       end if
       nB_t(atom3) = k
    end do

    KdB_t   (1:v_value%size) = v_value%items(1:v_value%size)%Kd
    ndB_t   (1:v_value%size) = v_value%items(1:v_value%size)%nd
    deltaB_t(1:v_value%size) = v_value%items(1:v_value%size)%delta

    continue_flagB_t(:) = .true.
    do i = 2, size(v_group_begin)
       continue_flagB_t(v_group_begin(i)-1) = .false.
    end do

#ifdef DIHEDRAL_TABLE_DEBUG
    if( myrank == 0 ) then
       atom1B1 = atom1B
       atom3B1 = atom3B
       atom4B1 = atom4B
       KdB1    = KdB
       ndB1    = ndB
       deltaB1 = deltaB
       nB1     = nB

       atom1B1 = 0
       atom3B1 = 0
       atom4B1 = 0
       KdB1    = 0.0d0
       ndB1    = 0.0d0
       deltaB1 = 0.0d0
       nB1     = 0

       do atom2 = 1, npara
          do i = 1, nB_t(atom2)
             atom1 = atom1B_t(atom2,i)
             atom3 = atom3B_t(atom2,i)
             atom4 = atom4B_t(atom2,i)
             index = indexB_t(atom2,i)
             if( index == 0 ) cycle
             do while(.true.)
                Kd    = KdB_t   (index)
                nd    = ndB_t   (index)
                delta = deltaB_t(index)
                k = nB1(atom2) + 1
                atom1B1(atom2,k) = atom1
                atom3B1(atom2,k) = atom3
                atom4B1(atom2,k) = atom4
                KdB1   (atom2,k) = Kd
                ndB1   (atom2,k) = nd
                deltaB1(atom2,k) = delta
                nB1    (atom2)   = k
                if( .not. continue_flagB_t(index) ) exit
                index = index + 1
             end do
          end do
       end do

       ! Remove Kd==0
       atom1B2 = atom1B
       atom3B2 = atom3B
       atom4B2 = atom4B
       KdB2    = KdB
       ndB2    = ndB
       deltaB2 = deltaB
       nB2     = nB
       do i = 1, npara
          k = 0
          do j = 1, nB2(i)
             if (KdB2(i,j) == 0.0d0) then
                cycle
             end if
             k = k + 1
             atom1B2(i,k) = atom1B2(i,j)
             atom3B2(i,k) = atom3B2(i,j)
             atom4B2(i,k) = atom4B2(i,j)
             KdB2   (i,k) = KdB2   (i,j)
             ndB2   (i,k) = ndB2   (i,j)
             deltaB2(i,k) = deltaB2(i,j)
          end do
          nB2(i) = k
       end do

       do i = 1, npara
          if( nB1(i) /= nB2(i)) then
             write(*,*) "ERROR: ", i, " th parameters of dihedralB are inconsistent."
             stop
          end if
          do j = 1, nB1(i)
             if(  atom1B1(i,j) /= atom1B2(i,j) .or. &
                  atom3B1(i,j) /= atom3B2(i,j) .or. &
                  atom4B1(i,j) /= atom4B2(i,j) .or. &
                  KdB1   (i,j) /= KdB2   (i,j) .or. &
                  ndB1   (i,j) /= ndB2   (i,j) .or. &
                  deltaB1(i,j) /= deltaB2(i,j) ) then
                write(*,*) "ERROR: (", i, j, ") th parameters of dihedralB are inconsistent."
                write(*,*) atom1B1(i,j), atom1B2(i,j)
                write(*,*) atom3B1(i,j), atom3B2(i,j)
                write(*,*) atom4B1(i,j), atom4B2(i,j)
                write(*,*) KdB1   (i,j), KdB2   (i,j)
                write(*,*) ndB1   (i,j), ndB2   (i,j)
                write(*,*) deltaB1(i,j), deltaB2(i,j)
                stop
             end if
          end do
       end do
       write(*,*) "Reproduce dihedralB Check Passed."

    endif
#endif /*DIHEDRAL_TABLE_DEBUG*/

  end subroutine set_md_charmm_dihedralB_table
#endif /*DIHEDRAL_TABLE*/

end module dihedral
