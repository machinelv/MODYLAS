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
!! \brief  Module and subroutines to void 1-2, 1-3 interactions.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to parse void pairs.
!! \author Kensuke Iwahashi
!<
module parse_void123
  implicit none
  ! LJ
  type yyparse_ljvoid
     integer(4) :: atom1,atom2
     type(yyparse_ljvoid),pointer :: next
  end type yyparse_ljvoid
  type(yyparse_ljvoid),pointer :: yyparse_ljvoid_top(:)
  integer(4),allocatable :: yyparse_ljvoid_n(:)
  integer(4) :: yyparse_ljvoid_total=0
  ! Coulomb
  type yyparse_coulombvoid
     integer(4) :: atom1,atom2
     type(yyparse_coulombvoid),pointer :: next
  end type yyparse_coulombvoid
  type(yyparse_coulombvoid),pointer :: yyparse_coulombvoid_top(:)
  integer(4),allocatable :: yyparse_coulombvoid_n(:)
  integer(4) :: yyparse_coulombvoid_total=0

contains

  subroutine fmod_alloc_yyparse_lj_top(ivalue)
    use parse_special
    implicit none
    integer(4), intent(in) :: ivalue
    integer(4) :: i
    type(yyparse_lj),pointer :: top
    type(yyparse_ljvoid),pointer :: topv
    allocate(yyparse_lj_top(ivalue))
    allocate(yyparse_ljvoid_top(ivalue))
    do i=1,ivalue
       allocate(top)
       allocate(topv)
       top%atom1 = -1
       top%atom2 = -1
       topv%atom1 = -1
       topv%atom2 = -1
       nullify(top%next)
       nullify(topv%next)
       yyparse_lj_top(i) = top
       yyparse_ljvoid_top(i) = topv
    enddo
  end subroutine fmod_alloc_yyparse_lj_top
!----------------------------------------------------------------------
  subroutine fmod_alloc_yyparse_lj_n(ivalue)
    use parse_special
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(yyparse_lj_n(ivalue))     ! special
    allocate(yyparse_ljvoid_n(ivalue)) ! void
    yyparse_lj_n = 0
    yyparse_ljvoid_n = 0
  end subroutine fmod_alloc_yyparse_lj_n
!----------------------------------------------------------------------
  subroutine fmod_alloc_yyparse_coulomb_top(ivalue)
    use parse_special
    implicit none
    integer(4), intent(in) :: ivalue
    integer(4) :: i
    type(yyparse_coulomb),pointer :: top
    type(yyparse_coulombvoid),pointer :: topv
    allocate(yyparse_coulomb_top(ivalue))     ! special
    allocate(yyparse_coulombvoid_top(ivalue)) ! void
    do i=1,ivalue
       allocate(top)
       allocate(topv)
       top%atom1 = -1
       top%atom2 = -1
       topv%atom1 = -1
       topv%atom2 = -1
       nullify(top%next)
       nullify(topv%next)
       yyparse_coulomb_top(i) = top
       yyparse_coulombvoid_top(i) = topv
    enddo
  end subroutine fmod_alloc_yyparse_coulomb_top
!----------------------------------------------------------------------
  subroutine fmod_alloc_yyparse_coulomb_n(ivalue)
    use parse_special
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(yyparse_coulomb_n(ivalue))     ! special
    allocate(yyparse_coulombvoid_n(ivalue)) ! void
    yyparse_coulomb_n = 0
    yyparse_coulombvoid_n = 0
  end subroutine fmod_alloc_yyparse_coulomb_n
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to keep LJ void pair in liner list temporally.
!! \author Kensuke Iwahashi
!<
  subroutine add_to_lj_void(atom01,atom02,line_number)
    use trajectory_org
    use segments
    use param
    implicit none
    integer(4) :: atom01, atom02
    integer(4) :: atom1,atom2,line_number
    type(yyparse_ljvoid),pointer :: new1,new2,p,next
    if (seg_ready) then
       atom1 = irearrange(atom01 + 1)
       atom2 = irearrange(atom02 + 1)
    else
       call segment_is_not_ready
       return
    endif
    if (atom1 > npara) then
       write(0,'(a,i0)')  'ERROR:'! line ', line_number
       write(0,'(a,i0,a)') &
            &    'The number of ljs is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom2 > npara) then
       write(0,'(a,i0)')  'ERROR:'! line ', line_number
       write(0,'(a,i0,a)') &
            &    'The number of ljs is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    allocate(new1)
    nullify(new1%next)
    new1%atom1 = atom1
    new1%atom2 = atom2
    p => yyparse_ljvoid_top(atom1)
    do while (.true.)
       next => p%next
       if (.not. associated(next)) then
          p%next => new1
          yyparse_ljvoid_n(atom1) = yyparse_ljvoid_n(atom1) + 1
          yyparse_ljvoid_total = yyparse_ljvoid_total + 1
          exit
       endif
       !       doubled must be merged
       if (atom2 == next%atom2) then
          deallocate(new1)
          exit
       else if (atom2 < next%atom2) then
          new1%next => p%next
          p%next => new1
          yyparse_ljvoid_n(atom1) = yyparse_ljvoid_n(atom1) + 1
          yyparse_ljvoid_total = yyparse_ljvoid_total + 1
          exit
       endif
       p => next
    enddo
    allocate(new2)
    nullify(new2%next)
    new2%atom1 = atom2
    new2%atom2 = atom1
    p => yyparse_ljvoid_top(atom2)
    !     insert to list, sorting from smaller
    do while (.true.)
       next => p%next
       if (.not. associated(next)) then
          p%next => new2
          yyparse_ljvoid_n(atom2) = yyparse_ljvoid_n(atom2) + 1
          yyparse_ljvoid_total = yyparse_ljvoid_total + 1
          exit
       endif
       !       doubled must be merged
       if (atom1 == next%atom2) then
          deallocate(new2)
          exit
       else if (atom1 < next%atom2) then
          new2%next => p%next
          p%next => new2
          yyparse_ljvoid_n(atom2) = yyparse_ljvoid_n(atom2) + 1
          yyparse_ljvoid_total = yyparse_ljvoid_total + 1
          exit
       endif
       p => next
    enddo
  end subroutine add_to_lj_void
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to keep coulomb void pair in liner list temporally.
!! \author Kensuke Iwahashi
!<
  subroutine add_to_coulomb_void(atom01,atom02,line_number)
    use trajectory_org
    use segments
    use param
    implicit none
    integer(4) :: atom01, atom02
    integer(4) :: atom1,atom2,line_number
    type(yyparse_coulombvoid),pointer :: new1,new2,p,next

    if (seg_ready) then
       atom1 = irearrange(atom01 + 1)
       atom2 = irearrange(atom02 + 1)
    else
       call segment_is_not_ready
       return
    endif
    if (atom1 > npara) then
       write(0,'(a,i0)')  'ERROR:'! line ', line_number
       write(0,'(a,i0,a)') &
            &    'The number of coulombs is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom2 > npara) then
       write(0,'(a,i0)')  'ERROR:'! line ', line_number
       write(0,'(a,i0,a)') &
            &    'The number of coulombs is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    !	  making a new element of list
    !	  it must be normalized but must keep information to which atom
    !	  the charge belongs
    allocate(new1)
    nullify(new1%next)
    new1%atom1 = atom1
    new1%atom2 = atom2
    !     insert to list, sorting from smaller
    p => yyparse_coulombvoid_top(atom1)
    !     p => yyparse_coulomb_top(atom1)
    do while (.true.)
       next => p%next
       if (.not. associated(next)) then
          p%next => new1
          yyparse_coulombvoid_n(atom1)=yyparse_coulombvoid_n(atom1)+1
          yyparse_coulombvoid_total = yyparse_coulombvoid_total + 1
          exit
       endif
       !       doubled must be merged
       if (atom2 == next%atom2) then
          !         next%charge1 = charge
          !         next%is_defined1 = .true.
          deallocate(new1)
          exit
       else if (atom2 < next%atom2) then
          new1%next => p%next
          p%next => new1
          yyparse_coulombvoid_n(atom1)=yyparse_coulombvoid_n(atom1)+1
          yyparse_coulombvoid_total = yyparse_coulombvoid_total + 1
          !         yyparse_coulomb_n(atom1)=yyparse_coulomb_n(atom1)+1
          !         yyparse_coulomb_total = yyparse_coulomb_total + 1
          exit
       endif
       p => next
    enddo
    !
    allocate(new2)
    nullify(new2%next)
    new2%atom1 = atom2
    new2%atom2 = atom1
    !     new2%charge2 = charge
    !     new2%is_defined1 = .false.
    !     new2%is_defined2 = .true.
    !     insert to list, sorting from smaller
    p => yyparse_coulombvoid_top(atom2)
    !     p => yyparse_coulomb_top(atom2)
    do while (.true.)
       next => p%next
       if (.not. associated(next)) then
          p%next => new2
          yyparse_coulombvoid_n(atom2)=yyparse_coulombvoid_n(atom2)+1
          yyparse_coulombvoid_total = yyparse_coulombvoid_total + 1
          !         yyparse_coulomb_n(atom2)=yyparse_coulomb_n(atom2)+1
          !         yyparse_coulomb_total = yyparse_coulomb_total + 1
          exit
       endif
       !       doubled must be merged
       if (atom1 == next%atom2) then
          !         next%charge2 = charge
          !         next%is_defined2 = .true.
          deallocate(new2)
          exit
       else if (atom1 < next%atom2) then
          new2%next => p%next
          p%next => new2
          yyparse_coulombvoid_n(atom2)=yyparse_coulombvoid_n(atom2)+1
          yyparse_coulombvoid_total = yyparse_coulombvoid_total + 1
          !         yyparse_coulomb_n(atom2)=yyparse_coulomb_n(atom2)+1
          !         yyparse_coulomb_total = yyparse_coulomb_total + 1
          exit
       endif
       p => next
    enddo
  end subroutine add_to_coulomb_void
!-----------------------------------------------------------------------
  subroutine add_to_mol_coulomb_void(ivalue1,atom1,atom2)
    use mol_info, only : para_start
    implicit none
    integer(4), intent(in) :: ivalue1
    integer(4), intent(inout) :: atom1,atom2

    atom1 = atom1 + para_start(ivalue1)
    atom2 = atom2 + para_start(ivalue1)

    call add_to_coulomb_void(atom1-1,atom2-1,1)

  end subroutine add_to_mol_coulomb_void
!-----------------------------------------------------------------------
  subroutine add_to_mol_lj_void(ivalue1, atom1, atom2)
    use mol_info, only : para_start
    implicit none
    integer(4), intent(in) :: ivalue1
    integer(4), intent(inout) :: atom1, atom2

    atom1 = atom1 + para_start(ivalue1)
    atom2 = atom2 + para_start(ivalue1)
    call add_to_lj_void(atom1-1, atom2-1, 1)
  end subroutine add_to_mol_lj_void

end module parse_void123
!----------------------------------------------------------------------
!>
!! \brief  Module to void 1-2, 1-3 interactions.
!! \author Kensuke Iwahashi, Yoshimichi Andoh, Kazushi Fujimoto
!<
!----------------------------------------------------------------------
module void123
  use omp_lib
  use parse_void123
  use mpi_tool
  use domain, only : lxdiv, lydiv, lzdiv
  use subcell
  implicit none
  integer(4) :: nvoid=0
  integer(4), allocatable :: void_atom1(:,:),void_atom2(:,:)
  integer(4), allocatable :: void_n(:)
  integer(4) :: nljvoid=0,ncoulombvoid=0
  integer(4), allocatable :: ljvoid_atom1(:,:), ljvoid_atom2(:,:)
  integer(4), allocatable :: ljvoid_n(:)
  real(8), allocatable :: ljvoid_epsilon(:,:), ljvoid_R(:,:)
  integer(4) maxj0
  integer(4),parameter::maxj0_per_i0=3
  integer(4), allocatable :: void_numL(:),void_topL(:)
  integer(4), allocatable :: void_list_j0(:)

contains

!----------------------------------------------------------------------
!>
!! \brief  Subroutines to remove 1-2, 1-3 LJ interactions.
!! \author Yoshimichi Andoh
!<
  subroutine remove_void123lj()
!----------------------------------------------------------------------
    use precision !to uniform the writing style of precision.
    use trajectory_mpi
    use lj_mod
    use param
    use md_monitors
    use md_const
    use forces
    use atom_virial
#ifdef ONEPROC_AXIS
    use boundary, only : pbc_pair
#endif
    use table_functions, only : &
  &       table_density_cut2, &
  &       lj12_table,lj12_table_delta, &
  &       lj6_table,lj6_table_delta, &
  &       lj12_force_table,lj12_force_table_delta, &
  &       lj6_force_table,lj6_force_table_delta
    use dr_cntl, only : nbd
    implicit none
    integer(4) :: ii,i,i0,j,j0,ipar,jpar,k0
    integer(4) :: isp
    integer(4) :: iam
    integer(4) :: icx0,icy0,icz0
    
    real(kind=wrp_void) :: epsi,Ri,eij,Rij
    real(kind=wrp_void) :: xli,yli,zli,rx,ry,rz,r2
    real(kind=wrp_void) :: r2_r, Rr6, Rr12, rUlj6,rUlj12, coef_lj
    real(kind=wrp_void) :: R3,R6,R12,r_03,r_06,r_12
    real(kind=wrp_void) :: tlx,tly,tlz,rc,rc_r
    real(kind=wrp_void) :: void_scale=0.5d0 !Chalf !! double-counted (void)    
    
    real(kind=dp_void) :: stlx, stly, stlz
    real(8) :: v11,v22,v33
    real(8) :: v21,v31,v32
    real(8) :: wk_v11, wk_v22, wk_v33
    real(8) :: wk_v21, wk_v31, wk_v32
    integer(4) :: table_point
    real(8) :: table_point_real, interpolate_distance
    real(8) :: lj12_force_coef, lj6_force_coef
    real(8) :: elj12,elj6

    if(nvoid==0)return

    rUlj6 = 0d0 !Czero
    rUlj12= 0d0 !Czero
    iam = 0
    wk_v11=0d0;wk_v22=0d0;wk_v33=0d0
    wk_v21=0d0;wk_v31=0d0;wk_v32=0d0

!$omp parallel default(none) &
!$omp& private(ii,iam,i0,j0,isp) &
!$omp& private(icx0,icy0,icz0) &
!$omp& private(i,ipar,epsi,Ri) &
!$omp& private(j,jpar,eij,Rij) &
!$omp& private(xli,yli,zli,rx,ry,rz) &
!$omp& private(r2,r2_r,Rr6,Rr12,coef_lj) &
!$omp& private(tlx,tly,tlz,stlx,stly,stlz) &
!$omp& private(R3,R6,R12) &
!$omp& shared(void_n,void_atom2) &
!$omp& shared(m2i,i2m,paranum,epsilon_sqrt,R_half) &
!$omp& shared(wkxyz,w3_f) &
!$omp& shared(tag,na_per_cell) &
!$omp& shared(lxdiv,lydiv,lzdiv) &
#ifdef CHARMMFSW
!$omp&  private(r_03,r_06,r_12,k_f12,k_f6,Sr,rc,rc_r) &
!$omp&  private(dV2_12,dV2_6) &
!$omp&  shared(c6a,c6b,c6c,c6d,iron2,iRoff3,iRoff6,Ron,Roff,ib3) &
!$omp&  shared(c12a,c12b,Ron2,Roff2) &
#endif
!$omp& private(v11,v22,v33,v21,v31,v32) &
!$omp& reduction(+:wk_v11,wk_v22,wk_v33,wk_v21,wk_v31,wk_v32) &
!$omp&  private(table_point_real,table_point) &
!$omp&  private(interpolate_distance) &
!$omp&  private(lj12_force_coef,lj6_force_coef) &
!$omp&  private(elj12,elj6) &
!$omp&  shared(table_density_cut2) &
!$omp&  shared(lj12_table,lj6_table) &
!$omp&  shared(lj12_table_delta,lj6_table_delta) &
!$omp&  shared(lj12_force_table,lj6_force_table) &
!$omp&  shared(lj12_force_table_delta,lj6_force_table_delta) &
!$omp& reduction(+:rUlj12,rUlj6)
!$    iam = omp_get_thread_num()
    do ii=1,lxdiv*lydiv*lzdiv
!default: MTD
       icz0=mod(ii-1,lzdiv)     +nbd
       icy0=mod(ii-1,lzdiv*lydiv)
       icy0=icy0/lzdiv          +nbd
       icx0=(ii-1)/(lzdiv*lydiv)+nbd
!$omp do
       do i0=tag(icz0,icy0,icx0), &
            &      tag(icz0,icy0,icx0)+na_per_cell(icz0,icy0,icx0)-1
          i   = m2i(i0)
          ipar= paranum(i)
          epsi= epsilon_sqrt(ipar)
          Ri  = R_half(ipar)
          xli = wkxyz(1,i0)
          yli = wkxyz(2,i0)
          zli = wkxyz(3,i0)
          stlx=0d0 !Cdzero
          stly=0d0 !Cdzero
          stlz=0d0 !Cdzero          
          do isp=1,void_n(ipar)
             j    = void_atom2(ipar,isp)+(i-ipar)
             j0   = i2m(j)
             jpar = paranum(j)
             eij = epsi*epsilon_sqrt(jpar)
#if defined(OPLSAMBER) && !defined(GAFF)
             Rij = sqrt(4d0*Ri*R_half(jpar))
#else /** CHARMM **/
             Rij = Ri +R_half(jpar)
#endif
             rx = xli - wkxyz(1,j0)
             ry = yli - wkxyz(2,j0)
             rz = zli - wkxyz(3,j0)
#ifdef ONEPROC_AXIS
             call pbc_pair(rx,ry,rz)
#endif
             !^^^ LJ ^^^
             r2   = rx**2+ry**2+rz**2
             r2_r = 1d0/r2
#ifndef NOTABLE
             R3 =Rij**3
             R6 =R3 **2
             R12=R6 **2
      table_point_real     = table_density_cut2*r2_r
      table_point          = int(table_point_real)
      interpolate_distance = table_point_real - dble(table_point)
      ! ^^^ potential ^^^ 
      elj12 = lj12_table(table_point) &
     &      + interpolate_distance    &
     &      * lj12_table_delta(table_point)
      elj6  = lj6_table(table_point)  &
     &      + interpolate_distance    &
     &      * lj6_table_delta(table_point)
      elj12 =        elj12* eij * R12
      elj6  = -2d0 * elj6 * eij * R6
      rUlj12=rUlj12+elj12
      rUlj6 =rUlj6 +elj6 
      ! ^^^ force ^^^ 
      lj12_force_coef = lj12_force_table(table_point) &
     &                   + interpolate_distance       &
     &                   * lj12_force_table_delta(table_point)
      lj6_force_coef  = lj6_force_table(table_point)  &
     &                   + interpolate_distance       &
     &                   * lj6_force_table_delta(table_point)
!
      lj12_force_coef = lj12_force_coef * R12
      lj6_force_coef  = lj6_force_coef * R6
      coef_lj = 12d0 * eij * (lj12_force_coef - lj6_force_coef)
#else
             Rr6  = Rij*Rij*r2_r
             Rr6  = Rr6*Rr6*Rr6
             Rr12 = Rr6*Rr6
             rUlj12 = rUlj12 + eij*Rr12
             rUlj6  = rUlj6  + eij*(-2d0*Rr6)
             coef_lj = 12d0*eij*r2_r*(Rr12-Rr6)
#endif /**NOTABLE**/
             tlx=coef_lj*rx
             tly=coef_lj*ry
             tlz=coef_lj*rz
             stlx=stlx+tlx
             stly=stly+tly
             stlz=stlz+tlz
             v11 = rx* tlx
             v22 = ry* tly
             v33 = rz* tlz
             v21 = ry* tlx
             v31 = rz* tlx
             v32 = rz* tly
             wk_v11 = wk_v11 + v11
             wk_v22 = wk_v22 + v22
             wk_v33 = wk_v33 + v33
             wk_v21 = wk_v21 + v21
             wk_v31 = wk_v31 + v31
             wk_v32 = wk_v32 + v32
          enddo ! isp
          w3_f(1,i0,0)=w3_f(1,i0,0)-stlx
          w3_f(2,i0,0)=w3_f(2,i0,0)-stly
          w3_f(3,i0,0)=w3_f(3,i0,0)-stlz
       enddo ! i0
!$omp end do nowait
    enddo ! ii
!$omp end parallel

    rUlj6 =void_scale*rUlj6
    rUlj12=void_scale*rUlj12

    wk_p_energy = wk_p_energy - (rUlj6+rUlj12)
    wk_vir2(1,0)=wk_vir2(1,0)-wk_v11*void_scale
    wk_vir2(2,0)=wk_vir2(2,0)-wk_v22*void_scale
    wk_vir2(3,0)=wk_vir2(3,0)-wk_v33*void_scale
    wk_vir2(4,0)=wk_vir2(4,0)-wk_v21*void_scale
    wk_vir2(5,0)=wk_vir2(5,0)-wk_v31*void_scale
    wk_vir2(6,0)=wk_vir2(6,0)-wk_v32*void_scale

#ifdef DEBUGFCE
    wplj=wplj - (rUlj6+rUlj12)
#endif

  end subroutine remove_void123lj
!----------------------------------------------------------------------
!>
!! \brief  Subroutines to remove 1-2, 1-3 Coulomb interactions.
!! \author Yoshimichi Andoh
!<
  subroutine remove_void123cl()
!----------------------------------------------------------------------
    use precision !to uniform the writing style of precision.
    use trajectory_mpi
    use param
    use md_monitors
    use md_const
    use coulomb_mod
    use forces
    use atom_virial
#ifdef ONEPROC_AXIS
    use boundary, only : pbc_pair
#endif
    use dr_cntl, only : nbd
    implicit none
    integer(4) :: ii,i,i0,j,j0,ipar,jpar
    integer(4) :: isp
    integer(4) :: iam
    integer(4) :: icx0,icy0,icz0
    real(kind=wrp_void) :: Ci,Cij
    real(kind=wrp_void) :: xci,yci,zci,rcx,rcy,rcz
    real(kind=wrp_void) :: rc_r, rc2_r, rUcl, coef_cl
    real(kind=wrp_void) :: tcx,tcy,tcz,tmp_mdQQ4PiE
    real(kind=wrp_void) :: void_scale=0.5d0!Chalf !! double-counted (void)
    real(kind=dp_void) :: stcx, stcy, stcz
    real(8) :: v11,v22,v33
    real(8) :: v21,v31,v32
    real(8) :: wk_v11, wk_v22, wk_v33
    real(8) :: wk_v21, wk_v31, wk_v32

    if(nvoid==0)return

    tmp_mdQQ4PiE=md_QQ_4PiE
    rUcl = 0d0
    iam = 0
    wk_v11=0d0;wk_v22=0d0;wk_v33=0d0
    wk_v21=0d0;wk_v31=0d0;wk_v32=0d0

!$omp parallel default(none) &
!$omp& private(ii,iam,i0,j0,isp) &
!$omp& private(icx0,icy0,icz0) &
!$omp& private(i,ipar,Ci) &
!$omp& private(j,jpar,Cij) &
!$omp& private(xci,yci,zci,rcx,rcy,rcz) &
!$omp& private(rc_r,rc2_r,coef_cl) &
!$omp& private(tcx,tcy,tcz,stcx,stcy,stcz) &
!$omp& shared(ndatm,void_n,void_atom2) &
!$omp& shared(m2i,i2m,paranum,chgv,tmp_mdQQ4PiE) &
!$omp& shared(wkxyz,w3_f) &
!$omp& shared(lxdiv,lydiv,lzdiv) &
!$omp& private(v11,v22,v33,v21,v31,v32) &
!$omp& reduction(+:wk_v11,wk_v22,wk_v33,wk_v21,wk_v31,wk_v32) &
!$omp& shared(tag,na_per_cell) &
!$omp& reduction(+:rUcl)
!$    iam = omp_get_thread_num()
    do ii=1,lxdiv*lydiv*lzdiv
!default: MTD
       icz0=mod(ii-1,lzdiv)     +nbd
       icy0=mod(ii-1,lzdiv*lydiv)
       icy0=icy0/lzdiv          +nbd
       icx0=(ii-1)/(lzdiv*lydiv)+nbd
!$omp do
       do i0=tag(icz0,icy0,icx0), &
            &      tag(icz0,icy0,icx0)+na_per_cell(icz0,icy0,icx0)-1
          i   = m2i(i0)
          ipar= paranum(i)
          Ci  = tmp_mdQQ4PiE*chgv(ipar)
          xci = wkxyz(1,i0)
          yci = wkxyz(2,i0)
          zci = wkxyz(3,i0)
          stcx=0d0
          stcy=0d0
          stcz=0d0
          do isp=1,void_n(ipar)
             j    = void_atom2(ipar,isp)+(i-ipar)
             j0   = i2m(j)
             jpar = paranum(j)
             Cij = Ci *chgv(jpar)
             rcx = xci - wkxyz(1,j0)
             rcy = yci - wkxyz(2,j0)
             rcz = zci - wkxyz(3,j0)
#ifdef ONEPROC_AXIS
             call pbc_pair(rcx,rcy,rcz)
#endif
             !^^^ Coulomb ^^^
             rc2_r = 1d0/(rcx**2+rcy**2+rcz**2)
             rc_r  = sqrt(rc2_r)
             rUcl = rUcl + Cij*rc_r
             coef_cl = Cij*rc_r*rc2_r
             tcx=coef_cl*rcx
             tcy=coef_cl*rcy
             tcz=coef_cl*rcz
             stcx=stcx+tcx
             stcy=stcy+tcy
             stcz=stcz+tcz
             v11 = rcx* tcx
             v22 = rcy* tcy
             v33 = rcz* tcz
             v21 = rcy* tcx
             v31 = rcz* tcx
             v32 = rcz* tcy
             wk_v11 = wk_v11 + v11
             wk_v22 = wk_v22 + v22
             wk_v33 = wk_v33 + v33
             wk_v21 = wk_v21 + v21
             wk_v31 = wk_v31 + v31
             wk_v32 = wk_v32 + v32
          enddo ! isp
          w3_f(1,i0,0)=w3_f(1,i0,0)-stcx
          w3_f(2,i0,0)=w3_f(2,i0,0)-stcy
          w3_f(3,i0,0)=w3_f(3,i0,0)-stcz
       enddo ! i0
!$omp end do nowait
    enddo ! ii
!$omp end parallel

    rUcl  =void_scale*rUcl

    wk_p_energy = wk_p_energy - (rUcl)
    wk_vir2(1,0)=wk_vir2(1,0)-wk_v11*void_scale
    wk_vir2(2,0)=wk_vir2(2,0)-wk_v22*void_scale
    wk_vir2(3,0)=wk_vir2(3,0)-wk_v33*void_scale
    wk_vir2(4,0)=wk_vir2(4,0)-wk_v21*void_scale
    wk_vir2(5,0)=wk_vir2(5,0)-wk_v31*void_scale
    wk_vir2(6,0)=wk_vir2(6,0)-wk_v32*void_scale

#ifdef DEBUGFCE
    wpcl=wpcl - (rUcl)
#endif

  end subroutine remove_void123cl
    
!----------------------------------------------------------------------
  subroutine fmod_alloc_void_atom1(ivalue1, ivalue2)
    implicit none
    integer(4), intent(in) :: ivalue1, ivalue2
    allocate(void_atom1(ivalue1, ivalue2))
  end subroutine fmod_alloc_void_atom1
!----------------------------------------------------------------------
  subroutine fmod_alloc_void_atom2(ivalue1, ivalue2)
    implicit none
    integer(4), intent(in) :: ivalue1, ivalue2
    allocate(void_atom2(ivalue1, ivalue2))
  end subroutine fmod_alloc_void_atom2
!----------------------------------------------------------------------
  subroutine fmod_alloc_void_n(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(void_n(ivalue))
  end subroutine fmod_alloc_void_n
!-----------------------------------------------------------------------
  subroutine set_md_void
    use md_condition
    use param
    implicit none
    integer(4) :: i,j,k, maxi(1)
    type(yyparse_ljvoid),pointer :: p,freed

    nljvoid      = yyparse_ljvoid_total
    ncoulombvoid = yyparse_coulombvoid_total
    if (nljvoid.ne.ncoulombvoid)then
       write(0,*) 'ERROR:  nljvoid .ne. ncoulombvoid'
       write(0,*) '  check 1-2, -3 void pair in mdff'
       call modylas_abort()
    endif
    nvoid=nljvoid

    if (nvoid == 0)  return

    maxi    =  maxloc(yyparse_ljvoid_n)
    maxi(1) =  yyparse_ljvoid_n(maxi(1))
    call fmod_alloc_void_atom1(npara, maxi(1))
    call fmod_alloc_void_atom2(npara, maxi(1))
    call fmod_alloc_void_n(npara)
    void_atom1 = 0 ! void_atom1(npara,maxi(1))
    void_atom2 = 0 ! void_atom2(npara,maxi(1))
    void_n     = 0 ! void_n(npara)
    !
    !     copying void pair entries
    !
    k = 0
    do i = 1, npara
       void_n(i) = yyparse_ljvoid_n(i)
       p => yyparse_ljvoid_top(i)
       do j = 1,yyparse_ljvoid_n(i)
          p => p%next
          void_atom1(i,j) = p%atom1
          void_atom2(i,j) = p%atom2
       enddo
    enddo
  end subroutine set_md_void
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to read void atoms.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_md_void
    use trajectory_org
    use device_numbers, only : f_mdff
    use param
    use mpi_tool
    implicit none
    integer(4) :: num
    include 'mpif.h'
    integer(4) :: ierr

    if(myrank==0) then
       read(f_mdff) nvoid
       !       write(*,*) 'input,nvoid=',nvoid
    endif
    call MPI_Bcast(nvoid,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    if(nvoid==0) return

    if(myrank==0) read(f_mdff) num
    call MPI_Bcast(num,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)

    allocate(void_n(npara))
    allocate(void_atom1(npara,num))
    allocate(void_atom2(npara,num))

    if(myrank==0) then
       read(f_mdff) void_n
       read(f_mdff) void_atom1
       read(f_mdff) void_atom2
    endif
    call MPI_Bcast(void_n,npara,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(void_atom1,npara*num,MPI_INTEGER4,0, &
         &               MPI_COMM_WORLD,ierr)
    call MPI_Bcast(void_atom2,npara*num,MPI_INTEGER4,0, &
         &               MPI_COMM_WORLD,ierr)
  end subroutine read_md_void
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write void atoms.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_md_void
    use device_numbers, only : f_mdff
    implicit none
    integer(4) :: num
    write(f_mdff) nvoid
    if(nvoid==0) return
    num = size(void_atom1(:,:), 2)
    write(f_mdff) num
    write(f_mdff) void_n
    write(f_mdff) void_atom1
    write(f_mdff) void_atom2
  end subroutine write_md_void
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_md_void
    implicit none
    integer(4) :: num
    write(*,*) '[write_md_void]'
    write(*,*) nvoid
    if(nvoid==0) return
    num = size(void_atom1(:,:), 2)
    write(*,*) num
    write(*,*) void_n
    write(*,*) void_atom1
    write(*,*) void_atom2
  end subroutine write_memory_md_void

end module void123
