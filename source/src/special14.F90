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
!! \brief  Module and subroutines to calculate 1-4 special interaction.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to calculate 1-4 special interaction.
!! \author Kensuke Iwahashi
!<
module special14
  use omp_lib
  use parse_special
  use mpi_tool
  use domain, only : lxdiv, lydiv, lzdiv
  use subcell
#ifdef ONEPROC_AXIS
    use boundary, only : pbc_pair
#endif
  implicit none
  ! LJ
  integer(4) :: nljsp=0
  integer(4), allocatable :: ljsp_atom1(:,:), ljsp_atom2(:,:)
  integer(4), allocatable :: ljsp_n(:)
  real(8), allocatable :: ljsp_epsilon(:,:), ljsp_R(:,:)

  ! Coulomb
  integer(4) :: ncoulombsp=0
  integer(4),allocatable :: coulombsp_atom1(:,:) &
       &                        , coulombsp_atom2(:,:)
  real(8),allocatable :: coulombsp_charge_product(:,:)
  integer(4),allocatable :: coulombsp_top(:),coulombsp_n(:)
#ifdef DEBUGFCE
  real(8) :: Ulj14=0d0
  real(8) :: Ucl14=0d0
#endif

contains

  subroutine fmod_alloc_ljsp_atom1(ivalue1, ivalue2)
    implicit none
    integer(4), intent(in) :: ivalue1, ivalue2
    allocate(ljsp_atom1(ivalue1, ivalue2))
  end subroutine fmod_alloc_ljsp_atom1
!----------------------------------------------------------------------
  subroutine fmod_alloc_ljsp_atom2(ivalue1, ivalue2)
    implicit none
    integer(4), intent(in) :: ivalue1, ivalue2
    allocate(ljsp_atom2(ivalue1, ivalue2))
  end subroutine fmod_alloc_ljsp_atom2
!----------------------------------------------------------------------
  subroutine fmod_alloc_ljsp_n(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(ljsp_n(ivalue))
    ljsp_n = 0
  end subroutine fmod_alloc_ljsp_n
!----------------------------------------------------------------------
  subroutine fmod_alloc_ljsp_epsilon(ivalue1, ivalue2)
    implicit none
    integer(4), intent(in) :: ivalue1, ivalue2
    allocate(ljsp_epsilon(ivalue1, ivalue2))
  end subroutine fmod_alloc_ljsp_epsilon
!----------------------------------------------------------------------
  subroutine fmod_alloc_ljsp_r(ivalue1, ivalue2)
    implicit none
    integer(4), intent(in) :: ivalue1, ivalue2
    allocate(ljsp_R(ivalue1, ivalue2))
  end subroutine fmod_alloc_ljsp_r
!----------------------------------------------------------------------
  subroutine set_md_LJ_special
    use lj_mod
    use force_field_numbers
    use md_condition
    use md_oplsaa_special_divide
    use param
    implicit none
    integer(4) :: i,j,k, maxi(1)
    type(yyparse_lj),pointer :: p,freed

    nljsp = yyparse_lj_total

    if (nljsp==0)  return

    maxi    =  maxloc(yyparse_lj_n)
    maxi(1) =  yyparse_lj_n(maxi(1))
    call fmod_alloc_ljsp_atom1(npara, maxi(1))
    call fmod_alloc_ljsp_atom2(npara, maxi(1))
    call fmod_alloc_ljsp_epsilon(npara, maxi(1))
    call fmod_alloc_ljsp_r(npara, maxi(1))
    call fmod_alloc_ljsp_n(npara)
    ljsp_atom1   = 0
    ljsp_atom2   = 0
    ljsp_epsilon = 0.0d0
    ljsp_r       = 0.0d0
    ljsp_n       = 0
    !
    !     copying special pair entries
    !
    k = 0
    do i = 1, npara
       ljsp_n(i) = yyparse_lj_n(i)
       p => yyparse_lj_top(i)
       do j = 1,yyparse_lj_n(i)
          p => p%next
          ljsp_atom1(i,j) = p%atom1
          ljsp_atom2(i,j) = p%atom2
          if (p%eps1 < 0.0) then
             p%eps1 = epsilon_sqrt(p%atom1)*epsilon_sqrt(p%atom1)
             p%r1   = R_half(p%atom1) * 2.0d0
          endif
          if (p%eps2 < 0.0) then
             p%eps2 = epsilon_sqrt(p%atom2)*epsilon_sqrt(p%atom2)
             p%r2   = R_half(p%atom2) * 2.0d0
          endif
          if(     md_condition__force_field == CHARMM .or. &
               &            md_condition__force_field == KREMER )then
             ljsp_epsilon(i,j) = dsqrt(p%eps1 * p%eps2)
             ljsp_r(i,j) = (p%r1 + p%r2) * 0.5d0
          else if(md_condition__force_field == OPLSAA .or. &
               &            md_condition__force_field == AMBER  )then
             ljsp_epsilon(i,j) = dsqrt(p%eps1 * p%eps2) &
                  &                        / md_oplsaa_lj_sp_divide
             ljsp_r(i,j) = dsqrt(p%r1 * p%r2)
          endif
       enddo
    enddo
  end subroutine set_md_LJ_special
!----------------------------------------------------------------------
  subroutine fmod_howmany_coulombsp(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    ncoulombsp = ivalue
  end subroutine fmod_howmany_coulombsp
!----------------------------------------------------------------------
  subroutine fmod_alloc_coulombsp_atom1(ivalue1,ivalue2)
    implicit none
    integer(4), intent(in) :: ivalue1, ivalue2
    allocate(coulombsp_atom1(ivalue1,ivalue2))
  end subroutine fmod_alloc_coulombsp_atom1
!----------------------------------------------------------------------
  subroutine fmod_alloc_coulombsp_atom2(ivalue1,ivalue2)
    implicit none
    integer(4), intent(in) :: ivalue1,ivalue2
    allocate(coulombsp_atom2(ivalue1,ivalue2))
  end subroutine fmod_alloc_coulombsp_atom2
!----------------------------------------------------------------------
  subroutine fmod_alloc_coulombsp_charge(ivalue1,ivalue2)
    implicit none
    integer(4), intent(in) :: ivalue1,ivalue2
    allocate(coulombsp_charge_product(ivalue1,ivalue2))
  end subroutine fmod_alloc_coulombsp_charge
!----------------------------------------------------------------------
  subroutine fmod_alloc_coulombsp_n(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(coulombsp_n(ivalue))
    coulombsp_n = 0
  end subroutine fmod_alloc_coulombsp_n
!-----------------------------------------------------------------------
  subroutine set_md_coulomb_special
    use coulomb_mod
    use force_field_numbers
    use md_condition
    use md_const
    use md_oplsaa_special_divide
    use param
    implicit none
    integer(4) :: i, j, k, maxi(1)
    type(yyparse_coulomb),pointer :: p,freed

    nclsp = yyparse_coulomb_total

    if (nclsp==0)  return

    maxi = maxloc(yyparse_coulomb_n)
    maxi(1) = yyparse_coulomb_n(maxi(1))
    call fmod_alloc_coulombsp_atom1(npara,maxi(1))
    call fmod_alloc_coulombsp_atom2(npara, maxi(1))
    call fmod_alloc_coulombsp_charge(npara,maxi(1))
    call fmod_alloc_coulombsp_n(npara)
    coulombsp_atom1 = 0
    coulombsp_atom2 = 0
    coulombsp_charge_product = 0.0d0
    coulombsp_n = 0

    !     registering special pair entries

    k = 0
    do i=1, npara
       coulombsp_n(i) = yyparse_coulomb_n(i)
       p => yyparse_coulomb_top(i)
       do j=1,yyparse_coulomb_n(i)
          k = k + 1
          p => p%next
          coulombsp_atom1(i,j) = p%atom1
          coulombsp_atom2(i,j) = p%atom2
          if (.not. p%is_defined1) then
             p%charge1 = chgv(p%atom1)
          endif
          if (.not. p%is_defined2) then
             p%charge2 = chgv(p%atom2)
          endif
          if(     md_condition__force_field == CHARMM .or. &
               &            md_condition__force_field == KREMER)then
             !! these lines not necessary
             !           coulombsp_charge_product(i,j) =
             !    &                        p%charge1 * p%charge2 * md_QQ_4PiE
          else if(md_condition__force_field == OPLSAA .or. &
               &            md_condition__force_field == AMBER )then
             !! these lines not necessary not necessary
             !           coulombsp_charge_product(i,j) =
             !    &                        p%charge1 * p%charge2 * md_QQ_4PiE
             !    &                        / md_oplsaa_coulomb_sp_divide
          endif
       enddo
    enddo

#ifndef DEALLOCBUG
    !     freeing temporary buffer
    do i=1, npara
       p => yyparse_coulomb_top(i)%next
       do j=1,yyparse_coulomb_n(i)
          freed => p
          p => p%next
          deallocate(freed)
       enddo
    enddo
    deallocate(yyparse_coulomb_top)
#endif
  end subroutine set_md_coulomb_special

!----------------------------------------------------------------------
  subroutine remove_special14lj()
!----------------------------------------------------------------------
    use precision
    use trajectory_mpi
    use lj_mod
    use param
    use md_monitors
    use md_const
    use coulomb_mod
    use forces
    use atom_virial
    use dr_cntl, only : nbd
    implicit none
    integer(4) :: iam
    integer(4) :: ii,i,i0,j,j0,ipar,jpar
    integer(4) :: isp
    integer(4) :: icx0,icy0,icz0

    real(kind=wrp_special) :: epsi0,Ri0,eij0,Rij0
    real(kind=wrp_special) :: eij,Rij
    real(kind=wrp_special) :: xli,yli,zli,rx,ry,rz,r2
    real(kind=wrp_special) :: r2_r, Rr6, Rr12
    real(kind=wrp_special) :: R3,R6,R12,r_03,r_06,r_12
    real(kind=wrp_special) :: r0Ulj6,r0Ulj12, coef_lj0
    real(kind=wrp_special) :: rUlj6,rUlj12, coef_lj
    real(kind=wrp_special) :: tlx,tly,tlz,rc,rc_r
    real(kind=wrp_special) :: special_scale=Chalf !! double-counted (special)

    real(8) :: v11,v22,v33
    real(8) :: v21,v31,v32
    real(8) :: wk_v11, wk_v22, wk_v33
    real(8) :: wk_v21, wk_v31, wk_v32

#ifdef DEBUGFCE
    real(8) :: wk_Usp14=0d0
#endif

    if(nljsp==0)return

    r0Ulj6 = 0.0d0 !Czero, since Czero is declared by wrp, and wrp_special may be different with wrp, so it is better to not use Czero.
    r0Ulj12 = 0.0d0 !Czero
    rUlj6  = 0.0d0 !Czero
    rUlj12 = 0.0d0 !Czero
    iam = 0
    wk_v11=0d0;wk_v22=0d0;wk_v33=0d0
    wk_v21=0d0;wk_v31=0d0;wk_v32=0d0

!$omp parallel default(none) &
!$omp& private(ii,iam,i0,j0,isp) &
!$omp& private(icx0,icy0,icz0) &
!$omp& private(i,ipar,epsi0,Ri0) &
!$omp& private(j,jpar,eij0,Rij0,eij,Rij) &
!$omp& private(xli,yli,zli,rx,ry,rz) &
!$omp& private(r2,r2_r,Rr6,Rr12,coef_lj0,coef_lj) &
!$omp& private(tlx,tly,tlz) &
!$omp& shared(ljsp_n) &
!$omp& shared(m2i,i2m,paranum,epsilon_sqrt,R_half) &
!$omp& shared(ljsp_atom2,ljsp_epsilon,ljsp_R) &
!$omp& shared(wkxyz,w3_f) &
!$omp& shared(tag,na_per_cell) &
!$omp& shared(lxdiv,lydiv,lzdiv) &
#ifdef CHARMMFSW
!$omp&  private(R3,R6,R12,r_03,r_06,r_12,k_f12,k_f6,Sr,rc,rc_r) &
!$omp&  private(dV2_12,dV2_6) &
!$omp&  shared(c6a,c6b,c6c,c6d,iron2,iRoff3,iRoff6,Ron,Roff,ib3) &
!$omp&  shared(c12a,c12b,Ron2,Roff2) &
#endif
!$omp& private(v11,v22,v33,v21,v31,v32) &
!$omp& reduction(+:wk_v11,wk_v22,wk_v33,wk_v21,wk_v31,wk_v32) &
!$omp& reduction(+:r0Ulj12,r0Ulj6) &
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
          i   =m2i(i0)
          ipar=paranum(i)
          !^^^normal^^^!
          epsi0=epsilon_sqrt(ipar)
          Ri0  =R_half(ipar)
          !^^^^^^^^^^^^!
          xli = wkxyz(1,i0)
          yli = wkxyz(2,i0)
          zli = wkxyz(3,i0)
          do isp=1,ljsp_n(ipar)
             j    = ljsp_atom2(ipar,isp)+(i-ipar)
             j0   = i2m(j)
             jpar = paranum(j)
             !^^^normal^^^!
             eij0 = epsi0*epsilon_sqrt(jpar)
             Rij0 = Ri0 + R_half(jpar)

             !^^^special^^^!
             eij = ljsp_epsilon(ipar,isp)
             Rij = ljsp_R(ipar,isp)
             !^^^^^^^^^^^^!
             rx = xli - wkxyz(1,j0)
             ry = yli - wkxyz(2,j0)
             rz = zli - wkxyz(3,j0)

#ifdef ONEPROC_AXIS
             call pbc_pair(rx,ry,rz)
#endif
             r2  = rx**2+ry**2+rz**2

             r2_r = 1.0d0/r2 !Cone

             !^^^normal^^^!
             Rr6  = Rij0*Rij0*r2_r
             Rr6  = Rr6*Rr6*Rr6
             Rr12 = Rr6*Rr6
             r0Ulj12 = r0Ulj12 + eij0* Rr12
             r0Ulj6  = r0Ulj6  + eij0*(-2d0*Rr6)
             coef_lj0= 12d0*r2_r*eij0*(Rr12-Rr6)
#ifdef CHARMMFSW
             !! Steinbach & Brooks, J.Comput.Chem.,15,667(1994). [eq(10)-(13)]
             Sr=1d0
             rc   = sqrt( r2 )
             rc_r = sqrt( r2_r )
             R3 =Rij0**3
             R6 =R3 **2
             R12=R6 **2
             if(rc.lt.Ron)then ! normal
                !^^ 12-power ^^!
                dV2_12=   -eij0*R12*c12a  !! not depend on rij
                r0Ulj12 = r0Ulj12 + dV2_12
                !^^  6-power  ^^!
                dV2_6 =2d0*eij0*R6 *c6a   !! not depend on rij
                r0Ulj6  = r0Ulj6  + dV2_6
             else ! with switching func.
                r_03=r2_r*rc_r
                r_06=r_03**2
                r_12=r_06**2
                !^^ 12-power  ^^!
                k_f12=     eij0*R12*c12b    !! not depend on rij
                r0Ulj12 = k_f12*(r_06-iRoff6)**2 !! depend on rij
                !^^  6-power  ^^!
                k_f6 =-2d0*eij0*R6 *c6b     !! not depend on rij
                r0Ulj6  = k_f6 *(r_03-iRoff3)**2 !! depend on rij
                Sr=(Roff2-r2)**2*(Roff2+2d0*r2-3d0*Ron2)*ib3 !! depend on rij
             endif
             coef_lj0=coef_lj0*Sr
#endif /**CHARMMFSW**/
             !^^^special^^^!
             Rr6  = Rij *Rij *r2_r
             Rr6  = Rr6*Rr6*Rr6
             Rr12 = Rr6*Rr6
             rUlj12  = rUlj12  + eij * Rr12

             rUlj6   = rUlj6   + eij *(-2.0d0*Rr6) !Ctwo
             coef_lj = 12.0d0*r2_r*eij *(Rr12-Rr6) !Ctwelve
#ifdef CHARMMFSW
             !! Steinbach & Brooks, J.Comput.Chem.,15,667(1994). [eq(10)-(13)]
             Sr=1d0
             R3 =Rij**3
             R6 =R3 **2
             R12=R6 **2
             if(rc.lt.Ron)then ! normal
                !^^ 12-power ^^!
                dV2_12=   -eij*R12*c12a  !! not depend on rij
                rUlj12 = rUlj12 + dV2_12
                !^^  6-power  ^^!
                dV2_6 =2d0*eij*R6 *c6a   !! not depend on rij
                rUlj6  = rUlj6  + dV2_6
             else ! with switching func.
                r_03=r2_r*rc_r
                r_06=r_03**2
                r_12=r_06**2
                !^^ 12-power  ^^!
                k_f12=     eij*R12*c12b    !! not depend on rij
                rUlj12 = k_f12*(r_06-iRoff6)**2 !! depend on rij
                !^^  6-power  ^^!
                k_f6 =-2d0*eij*R6 *c6b     !! not depend on rij
                rUlj6  = k_f6 *(r_03-iRoff3)**2 !! depend on rij
                Sr=(Roff2-r2)**2*(Roff2+2d0*r2-3d0*Ron2)*ib3 !! depend on rij
             endif
             coef_lj=coef_lj*Sr
#endif /**CHARMMFSW**/
             !^^^^^^^^^^^!
             tlx=(coef_lj0-coef_lj)*rx
             tly=(coef_lj0-coef_lj)*ry
             tlz=(coef_lj0-coef_lj)*rz
             w3_f(1,i0,iam)=w3_f(1,i0,iam)-tlx
             w3_f(2,i0,iam)=w3_f(2,i0,iam)-tly
             w3_f(3,i0,iam)=w3_f(3,i0,iam)-tlz
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
       enddo ! i0
!$omp end do
    enddo ! ii
!$omp end parallel
! The all pairs are doubly counted, so *0.5 applied on
    r0Ulj12=special_scale*r0Ulj12
    r0Ulj6 =special_scale*r0Ulj6
    rUlj12 =special_scale*rUlj12
    rUlj6  =special_scale*rUlj6

    wk_p_energy = wk_p_energy - (r0Ulj6-rUlj6) - (r0Ulj12-rUlj12)
    wk_vir2(1,0)=wk_vir2(1,0)-wk_v11*special_scale
    wk_vir2(2,0)=wk_vir2(2,0)-wk_v22*special_scale
    wk_vir2(3,0)=wk_vir2(3,0)-wk_v33*special_scale
    wk_vir2(4,0)=wk_vir2(4,0)-wk_v21*special_scale
    wk_vir2(5,0)=wk_vir2(5,0)-wk_v31*special_scale
    wk_vir2(6,0)=wk_vir2(6,0)-wk_v32*special_scale
#ifdef DEBUGFCE
    wplj=wplj - (r0Ulj6-rUlj6) - (r0Ulj12-rUlj12)
    wk_Usp14=rUlj12+rUlj6
    call mpi_allreduce(wk_Usp14,Ulj14,1, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,i)
    if(myrank==0) write(*,*) & 
#ifdef KCAL
    'Pot(LJ14) =',Ulj14*kJ_mol/4.184d0,'[kcal/mol] Note: CHARMM14sp'
#else
    'Pot(LJ14) =',Ulj14*kJ_mol,'[kJ/mol] Note: CHARMM14sp'
#endif
#endif
  end subroutine remove_special14lj
!----------------------------------------------------------------------
  subroutine remove_special14kremer()  ! for Kremer's potential
!----------------------------------------------------------------------
    use trajectory_mpi
    use lj_mod
    use param
    use md_monitors
    use md_const
    use coulomb_mod
    use forces
    use atom_virial
    use dr_cntl, only : nbd
    implicit none
    integer(4) :: iam
    integer(4) :: ii,i,i0,j,j0,ipar,jpar
    integer(4) :: isp
    integer(4) :: icx0,icy0,icz0
    real(8) :: eij,Rij
    real(8) :: xli,yli,zli,rx,ry,rz,r2
    real(8) :: r2_r, Rr6, Rr12
    real(8) :: R3,R6,R12,r_03,r_06,r_12
    real(8) :: rUlj6,rUlj12, coef_lj
    real(8) :: Cij, rc_r, rUcl, coef_cl, tmp_mdQQ4PiE
    real(8) :: tlx,tly,tlz,tcx,tcy,tcz
    real(8) :: special_scale=0.5d0 !! double-counted (special)
    real(8) :: v11,v22,v33
    real(8) :: v21,v31,v32
    real(8) :: wk_v11, wk_v22, wk_v33
    real(8) :: wk_v21, wk_v31, wk_v32

#ifdef CHARMMFSW
    if(myrank==0)then
       write(0,*) 'ERROR: -DKREMER and -DCHARMMFSW incompatible'
       call modylas_abort()
    endif
#endif

    if(nljsp==0)return

    tmp_mdQQ4PiE=md_QQ_4PiE
    rUlj6  =0d0
    rUlj12 =0d0
    rUcl  =0d0
    iam = 0
    wk_v11=0d0;wk_v22=0d0;wk_v33=0d0
    wk_v21=0d0;wk_v31=0d0;wk_v32=0d0

!$omp parallel default(none) &
!$omp& private(ii,iam,i0,j0,isp) &
!$omp& private(icx0,icy0,icz0) &
!$omp& private(i,ipar) &
!$omp& private(j,jpar,eij,Rij) &
!$omp& private(xli,yli,zli,rx,ry,rz) &
!$omp& private(r2,r2_r,Rr6,Rr12,coef_lj) &
!$omp& private(Cij,rc_r,coef_cl) &
!$omp& private(tlx,tly,tlz,tcx,tcy,tcz) &
!$omp& shared(ljsp_n,tmp_mdQQ4PiE) &
!$omp& shared(m2i,i2m,paranum,epsilon_sqrt,R_half) &
!$omp& shared(ljsp_atom2,ljsp_epsilon,ljsp_R,chgv) &
!$omp& shared(wkxyz,w3_f) &
!$omp& shared(tag,na_per_cell) &
!$omp& shared(lxdiv,lydiv,lzdiv) &
!$omp& private(v11,v22,v33,v21,v31,v32) &
!$omp& reduction(+:wk_v11,wk_v22,wk_v33,wk_v21,wk_v31,wk_v32) &
!$omp& reduction(+:rUlj12,rUlj6,rUcl)
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
          i   =m2i(i0)
          ipar=paranum(i)
          !^^^^^^^^^^^^!
          xli = wkxyz(1,i0)
          yli = wkxyz(2,i0)
          zli = wkxyz(3,i0)
          do isp=1,ljsp_n(ipar)
             j    = ljsp_atom2(ipar,isp)+(i-ipar)
             j0   = i2m(j)
             jpar = paranum(j)
             !^^^special^^^!
             eij = ljsp_epsilon(ipar,isp)    ! special
             IF(eij/=0) CYCLE
             !! 1-4 interaction are not removed
             !! (completely swithced on)
             !! else, removed (swithched off)
             eij = epsilon_sqrt(ipar)*epsilon_sqrt(jpar) ! original
             Rij = R_half(ipar) + R_half(jpar)           ! original
             Cij = tmp_mdQQ4PiE*chgv(ipar)*chgv(jpar)    ! original
             !^^^^^^^^^^^^!
             rx = xli - wkxyz(1,j0)
             ry = yli - wkxyz(2,j0)
             rz = zli - wkxyz(3,j0)
#ifdef ONEPROC_AXIS
             call pbc_pair(rx,ry,rz)
#endif
             r2  = rx**2+ry**2+rz**2
             r2_r = 1d0/r2
             rc_r = sqrt(r2_r)
             !^^^special LJ ^^^!
             Rr6  = Rij *Rij *r2_r
             Rr6  = Rr6*Rr6*Rr6
             Rr12 = Rr6*Rr6
             rUlj12  = rUlj12  + eij * Rr12
             rUlj6   = rUlj6   + eij *(-2d0*Rr6)
             coef_lj = 12d0*r2_r*eij *(Rr12-Rr6)
             !^^^special CL ^^^!
             rUcl    = rUcl    + Cij*rc_r
             coef_cl = Cij*rc_r*r2_r
             !^^^^^^^^^^^!
             tlx=coef_lj*rx
             tly=coef_lj*ry
             tlz=coef_lj*rz
             tcx=coef_cl*rx
             tcy=coef_cl*ry
             tcz=coef_cl*rz
             w3_f(1,i0,0)=w3_f(1,i0,0)-tlx-tcx
             w3_f(2,i0,0)=w3_f(2,i0,0)-tly-tcy
             w3_f(3,i0,0)=w3_f(3,i0,0)-tlz-tcz
             v11 = rx* (tlx+tcx)
             v22 = ry* (tly+tcy)
             v33 = rz* (tlz+tcz)
             v21 = ry* (tlx+tcx)
             v31 = rz* (tlx+tcx)
             v32 = rz* (tly+tcy)
             wk_v11 = wk_v11 + v11
             wk_v22 = wk_v22 + v22
             wk_v33 = wk_v33 + v33
             wk_v21 = wk_v21 + v21
             wk_v31 = wk_v31 + v31
             wk_v32 = wk_v32 + v32
          enddo ! isp
       enddo ! i0
!$omp end do
    enddo ! ii
!$omp end parallel

    rUlj12 =special_scale*rUlj12
    rUlj6  =special_scale*rUlj6
    rUcl   =special_scale*rUcl

    wk_p_energy = wk_p_energy - (rUlj6+rUlj12) - rUcl
    wk_vir2(1,0)=wk_vir2(1,0)-wk_v11*special_scale
    wk_vir2(2,0)=wk_vir2(2,0)-wk_v22*special_scale
    wk_vir2(3,0)=wk_vir2(3,0)-wk_v33*special_scale
    wk_vir2(4,0)=wk_vir2(4,0)-wk_v21*special_scale
    wk_vir2(5,0)=wk_vir2(5,0)-wk_v31*special_scale
    wk_vir2(6,0)=wk_vir2(6,0)-wk_v32*special_scale
#ifdef DEBUGFCE
    wplj=wplj - (rUlj6+rUlj12)
    wpcl=wpcl -  rUcl
#endif
  end subroutine remove_special14kremer
!----------------------------------------------------------------------
  subroutine remove_scaling14lj() !! for AMBER/OPLS
!----------------------------------------------------------------------
    use precision
    use trajectory_mpi
    use param
    use md_monitors
    use md_const
    use lj_mod
    use forces
    use atom_virial
    use md_oplsaa_special_divide
    use dihedral
    use md_condition
    use dr_cntl, only : nbd
#ifdef ONEPROC_AXIS
    use boundary, only : pbc_pair
#endif
    implicit none
    integer(4) :: ii,i,i0,j,j0,ipar,jpar,k0
    integer(4) :: isp
    integer(4) :: iam
    integer(4) :: icx0,icy0,icz0

    real(kind=wrp_scale) :: epsi,Ri,eij,Rij
    real(kind=wrp_scale) :: xli,yli,zli,rx,ry,rz,r2
    real(kind=wrp_scale) :: r2_r, Rr6, Rr12, rUlj6, rUlj12, coef_lj
    real(kind=wrp_scale) :: tlx,tly,tlz, oUr
    real(kind=wrp_scale) :: void_scale=0.5d0 !! double-counted
    
    real(8) :: v11,v22,v33
    real(8) :: v21,v31,v32
    real(8) :: wk_v11, wk_v22, wk_v33
    real(8) :: wk_v21, wk_v31, wk_v32

    real(kind=wrp_scale) :: special_scale    !! 1 - (scaling_factor)

#ifdef DEBUGFCE
    real(8) :: wk_rUljsum
    real(8) :: oscale     !! scaling_factor/(1-scaling_factor)
#endif

    if(ndihedralAs==0) return

    special_scale=1d0-1d0/md_oplsaa_lj_sp_divide
#ifdef DEBUGFCE
    oscale=(1d0/md_oplsaa_lj_sp_divide)/special_scale
#endif

    rUlj6  = 0.0d0 !Czero
    rUlj12 = 0.0d0 !Czero
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
!$omp& private(tlx,tly,tlz) &
!$omp& shared(special_scale) &
!$omp& shared(m2i,i2m,paranum,epsilon_sqrt,R_half) &
!$omp& shared(nArm,atom4Arm) &
!$omp& shared(wkxyz,w3_f,myrank) &
!$omp& shared(tag,na_per_cell) &
!$omp& shared(lxdiv,lydiv,lzdiv) &
!$omp& private(v11,v22,v33,v21,v31,v32) &
!$omp& reduction(+:wk_v11,wk_v22,wk_v33,wk_v21,wk_v31,wk_v32) &
!$omp& reduction(+:rUlj6,rUlj12)
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
          do j=1,nArm(ipar)
             j0   = i2m(atom4Arm(ipar,j)+(i-ipar))
             jpar = paranum(m2i(j0))

             eij = epsi*epsilon_sqrt(jpar)

#if defined(OPLSAMBER) && !defined(GAFF)
             Rij = sqrt(4d0*Ri*R_half(jpar))
#else /** CHARMM **/
             Rij = Ri +R_half(jpar)
#endif
             rx  = xli - wkxyz(1,j0)
             ry  = yli - wkxyz(2,j0)
             rz  = zli - wkxyz(3,j0)

#ifdef ONEPROC_AXIS
             call pbc_pair(rx,ry,rz)
#endif
             !^^^ LJ ^^^
             r2   = rx**2+ry**2+rz**2

             r2_r = 1d0/r2
             Rr6  = Rij*Rij*r2_r
             Rr6  = Rr6*Rr6*Rr6
             Rr12 = Rr6*Rr6
             eij=eij*special_scale
             rUlj12 = rUlj12 + eij*Rr12

             rUlj6  = rUlj6  + eij*(-2d0*Rr6)

             coef_lj = 12d0*eij*r2_r*(Rr12-Rr6)

             tlx=coef_lj*rx
             tly=coef_lj*ry
             tlz=coef_lj*rz
             w3_f(1,i0,iam)=w3_f(1,i0,iam)-tlx
             w3_f(2,i0,iam)=w3_f(2,i0,iam)-tly
             w3_f(3,i0,iam)=w3_f(3,i0,iam)-tlz
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
       enddo ! i0
!$omp end do
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
    wk_rUljsum=(rUlj6+rUlj12)*oscale
    wplj=wplj - (rUlj6+rUlj12) - wk_rUljsum
    call mpi_allreduce(wk_rUljsum,Ulj14,1, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,i)
    if(myrank==0) write(*,*) 'Pot(LJ14) =', &
#ifdef KCAL
 &        Ulj14*kJ_mol/4.184d0,'[kcal/mol] Note: OPLS/GAFF LJ14'
#else
 &        Ulj14*kJ_mol,      '[kJ/mol] Note: OPLS/GAFF LJ14'
#endif
#endif

  end subroutine remove_scaling14lj
!----------------------------------------------------------------------
  subroutine remove_scaling14cl() !! for AMBER/OPLS
!----------------------------------------------------------------------
    use precision
    use trajectory_mpi
    use param
    use md_monitors
    use md_const
    use coulomb_mod
    use forces
    use atom_virial
    use md_oplsaa_special_divide
    use dihedral
    use md_condition
    use dr_cntl, only : nbd
#ifdef ONEPROC_AXIS
    use boundary, only : pbc_pair
#endif
    implicit none
    integer(4) :: ii,i,i0,j,j0,ipar,jpar
    integer(4) :: isp
    integer(4) :: iam
    integer(4) :: icx0,icy0,icz0

    real(kind=wrp_scale) :: Ci,Cij
    real(kind=wrp_scale) :: xci,yci,zci,rcx,rcy,rcz
    real(kind=wrp_scale) :: rc_r, rc2_r, rUcl, coef_cl, oUr
    real(kind=wrp_scale) :: tcx,tcy,tcz,tmp_mdQQ4PiE
    real(kind=wrp_scale) :: void_scale=0.5d0 !! double-counted

    real(8) :: v11,v22,v33
    real(8) :: v21,v31,v32
    real(8) :: wk_v11, wk_v22, wk_v33
    real(8) :: wk_v21, wk_v31, wk_v32

    real(kind=wrp_scale) :: special_scale    !! 1 - (scaling_factor)

#ifdef DEBUGFCE
    real(8) :: wk_rUclsum !, rUclsum
    real(8) :: oscale     !! scaling_factor/(1-scaling_factor)
#endif

    if(ndihedralAs==0) return

    special_scale=1d0-1d0/md_oplsaa_coulomb_sp_divide
#ifdef DEBUGFCE
    oscale=(1d0/md_oplsaa_coulomb_sp_divide)/special_scale
#endif

    tmp_mdQQ4PiE=md_QQ_4PiE
    rUcl  =0d0
    iam = 0
    wk_v11=0d0;wk_v22=0d0;wk_v33=0d0
    wk_v21=0d0;wk_v31=0d0;wk_v32=0d0

!$omp parallel default(none) &
!$omp& private(ii,iam,i0,j0,isp) &
!$omp& private(icx0,icy0,icz0) &
!$omp& private(i,ipar,Ci) &
!$omp& private(j,jpar,Cij) &
!$omp& private(xci,yci,zci,rcx,rcy,rcz) &
!$omp& private(rc2_r,rc_r,coef_cl) &
!$omp& private(tcx,tcy,tcz) &
!$omp& shared(special_scale,tmp_mdQQ4PiE) &
!$omp& shared(m2i,i2m,paranum,chgv) &
!$omp& shared(nArm,atom4Arm) &
!$omp& shared(wkxyz,w3_f,myrank) &
!$omp& shared(tag,na_per_cell) &
!$omp& shared(lxdiv,lydiv,lzdiv) &
!$omp& private(v11,v22,v33,v21,v31,v32) &
!$omp& reduction(+:wk_v11,wk_v22,wk_v33,wk_v21,wk_v31,wk_v32) &
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
          do j=1,nArm(ipar)
             j0   = i2m(atom4Arm(ipar,j)+(i-ipar))
             jpar = paranum(m2i(j0))

             Cij  = Ci *chgv(jpar) * special_scale
             rcx  = xci - wkxyz(1,j0)
             rcy  = yci - wkxyz(2,j0)
             rcz  = zci - wkxyz(3,j0)

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
             w3_f(1,i0,iam)=w3_f(1,i0,iam)-tcx
             w3_f(2,i0,iam)=w3_f(2,i0,iam)-tcy
             w3_f(3,i0,iam)=w3_f(3,i0,iam)-tcz
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
       enddo ! i0
!$omp end do
    enddo ! ii
!$omp end parallel

    rUcl  =void_scale*rUcl

    wk_p_energy =wk_p_energy -(rUcl)
    wk_vir2(1,0)=wk_vir2(1,0)-wk_v11*void_scale
    wk_vir2(2,0)=wk_vir2(2,0)-wk_v22*void_scale
    wk_vir2(3,0)=wk_vir2(3,0)-wk_v33*void_scale
    wk_vir2(4,0)=wk_vir2(4,0)-wk_v21*void_scale
    wk_vir2(5,0)=wk_vir2(5,0)-wk_v31*void_scale
    wk_vir2(6,0)=wk_vir2(6,0)-wk_v32*void_scale
#ifdef DEBUGFCE
    wk_rUclsum=rUcl * oscale
    wpcl=wpcl -  rUcl - wk_rUclsum
    call mpi_allreduce(wk_rUclsum,Ucl14,1, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,i)
    if(myrank==0) write(*,*) 'Pot(CL14) =', &
#ifdef KCAL
  &   Ucl14*kJ_mol/4.184d0,'[kcal/mol] Note: OPLS/GAFF CL14'
#else
  &   Ucl14*kJ_mol,        '[kJ/mol] Note: OPLS/GAFF CL14'
#endif
#endif

  end subroutine remove_scaling14cl
!-------------------------------------------------------------------------
#ifdef DEBUGFCE
!>
!! \brief  Subroutine to calculate 1-4 coulomb interaction (CHARMM), 
!!         whereas this routine is used to check CL14 with type=charmm.
!!         Calculated values are not used in main body of MD simulation.
!! \author Yoshimichi ANDOH
!<
  subroutine dummy_calc14cl() !! for CHARMM
!-------------------------------------------------------------------------
    use precision
    use trajectory_mpi
    use dihedral
    use param
    use md_const
    use coulomb_mod
    use forces
    use atom_virial
    use dihedral
    use dr_cntl, only : nbd
#ifdef ONEPROC_AXIS
    use boundary, only : pbc_pair
#endif
    implicit none
    integer(4) :: ii,i,i0,j,j0,ipar,jpar
    integer(4) :: isp
    integer(4) :: iam
    integer(4) :: icx0,icy0,icz0

    real(kind=wrp_scale) :: Ci,Cij
    real(kind=wrp_scale) :: xci,yci,zci,rcx,rcy,rcz
    real(kind=wrp_scale) :: rc_r, rc2_r, rUcl, coef_cl, oUr
    real(kind=wrp_scale) :: tcx,tcy,tcz,tmp_mdQQ4PiE
    real(kind=wrp_scale) :: void_scale=0.5d0 !! double-counted

    real(8) :: v11,v22,v33
    real(8) :: v21,v31,v32
    real(8) :: wk_v11, wk_v22, wk_v33
    real(8) :: wk_v21, wk_v31, wk_v32
    real(kind=wrp_scale) :: special_scale    !! 1 - (scaling_factor)

    if(ndihedralAs==0) return

    special_scale=1d0 ! always 1

    tmp_mdQQ4PiE=md_QQ_4PiE
    rUcl  =0d0
    iam = 0
    wk_v11=0d0;wk_v22=0d0;wk_v33=0d0
    wk_v21=0d0;wk_v31=0d0;wk_v32=0d0

    do ii=1,lxdiv*lydiv*lzdiv
!default: MTD
       icz0=mod(ii-1,lzdiv)     +nbd
       icy0=mod(ii-1,lzdiv*lydiv)
       icy0=icy0/lzdiv          +nbd
       icx0=(ii-1)/(lzdiv*lydiv)+nbd
       do i0=tag(icz0,icy0,icx0), &
            &      tag(icz0,icy0,icx0)+na_per_cell(icz0,icy0,icx0)-1
          i   = m2i(i0)
          ipar= paranum(i)
          Ci  = tmp_mdQQ4PiE*chgv(ipar)
          xci = wkxyz(1,i0)
          yci = wkxyz(2,i0)
          zci = wkxyz(3,i0)
          do j=1,nArm(ipar)  
             j0   = i2m(atom4Arm(ipar,j)+(i-ipar))
             jpar = paranum(m2i(j0))

             Cij  = Ci *chgv(jpar) * special_scale
             rcx  = xci - wkxyz(1,j0)
             rcy  = yci - wkxyz(2,j0)
             rcz  = zci - wkxyz(3,j0)

#ifdef ONEPROC_AXIS
             call pbc_pair(rcx,rcy,rcz)
#endif
             !^^^ Coulomb ^^^
             rc2_r = 1d0/(rcx**2+rcy**2+rcz**2)

             rc_r  = sqrt(rc2_r)
             rUcl = rUcl + Cij*rc_r
!!
!! The following lines are commented out, since only potential energy
!! is a target for interest.
!!
!            coef_cl = Cij*rc_r*rc2_r
!            tcx=coef_cl*rcx
!            tcy=coef_cl*rcy
!            tcz=coef_cl*rcz
!            w3_f(1,i0,iam)=w3_f(1,i0,iam)-tcx
!            w3_f(2,i0,iam)=w3_f(2,i0,iam)-tcy
!            w3_f(3,i0,iam)=w3_f(3,i0,iam)-tcz
!            v11 = rcx* tcx
!            v22 = rcy* tcy
!            v33 = rcz* tcz
!            v21 = rcy* tcx
!            v31 = rcz* tcx
!            v32 = rcz* tcy
!            wk_v11 = wk_v11 + v11
!            wk_v22 = wk_v22 + v22
!            wk_v33 = wk_v33 + v33
!            wk_v21 = wk_v21 + v21
!            wk_v31 = wk_v31 + v31
!            wk_v32 = wk_v32 + v32
          enddo ! isp
       enddo ! i0
    enddo ! ii

    rUcl  =void_scale*rUcl  !! note: Every i-j pair double counted

#ifdef DEBUGFCE
    call mpi_allreduce(rUcl,Ucl14,1, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,i)
    if(myrank==0) &
#ifdef KCAL
    write(*,*) 'Pot(CL14) =',Ucl14*kJ_mol/4.184d0,'[kcal/mol] Note: CHARMM14'
#else
    write(*,*) 'Pot(CL14) =',Ucl14*kJ_mol,'[kJ/mol] Note: CHARMM14'
#endif
#endif

  end subroutine dummy_calc14cl
#endif
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to read LJ special atoms.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_md_LJ_special
    use trajectory_org
    use device_numbers, only : f_mdff
    use param
    use mpi_tool
    implicit none
    integer(4) :: num
    include 'mpif.h'
    integer(4) :: ierr

    if(myrank==0) then
       read(f_mdff) nljsp
    endif
    call MPI_Bcast(nljsp,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    if(nljsp==0) return

    if(myrank==0) read(f_mdff) num
    call MPI_Bcast(num,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)

    allocate(ljsp_n(npara))
    allocate(ljsp_atom1(npara,num))
    allocate(ljsp_atom2(npara,num))
    allocate(ljsp_epsilon(npara,num))
    allocate(ljsp_r(npara,num))

    if(myrank==0) then
       read(f_mdff) ljsp_n
       read(f_mdff) ljsp_atom1
       read(f_mdff) ljsp_atom2
       read(f_mdff) ljsp_epsilon
       read(f_mdff) ljsp_r
    endif
    call MPI_Bcast(ljsp_n,npara,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(ljsp_atom1,npara*num,MPI_INTEGER4,0, &
         &               MPI_COMM_WORLD,ierr)
    call MPI_Bcast(ljsp_atom2,npara*num,MPI_INTEGER4,0, &
         &               MPI_COMM_WORLD,ierr)
    call MPI_Bcast(ljsp_epsilon,npara*num,MPI_REAL8,0, &
         &               MPI_COMM_WORLD,ierr)
    call MPI_Bcast(ljsp_r,npara*num,MPI_REAL8,0, &
         &               MPI_COMM_WORLD,ierr)
  end subroutine read_md_LJ_special
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write LJ special atoms in mdff.bin File.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_md_LJ_special
    use device_numbers, only : f_mdff
    implicit none
    integer(4) :: num
    write(f_mdff) nljsp
    if(nljsp==0) return
    num = size(ljsp_atom1(:,:), 2)
    write(f_mdff) num
    write(f_mdff) ljsp_n
    write(f_mdff) ljsp_atom1
    write(f_mdff) ljsp_atom2
    write(f_mdff) ljsp_epsilon
    write(f_mdff) ljsp_r
  end subroutine write_md_LJ_special
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write LJ special atoms.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_md_LJ_special
    implicit none
    integer(4) :: num
    write(*,*) '[write_md_LJ_special]'
    write(*,*) nljsp
    if(nljsp==0) return
    num = size(ljsp_atom1(:,:), 2)
    write(*,*) num
    write(*,*) ljsp_n
    write(*,*) ljsp_atom1
    write(*,*) ljsp_atom2
    write(*,*) ljsp_epsilon
    write(*,*) ljsp_r
  end subroutine write_memory_md_LJ_special
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to read coulomb special.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_md_coulomb_special
    use trajectory_org
    use device_numbers, only : f_mdff
    use param
    use mpi_tool
    implicit none
    integer(4) :: num
    include 'mpif.h'
    integer(4) :: ierr

    if(myrank==0) then
       read(f_mdff) ncoulombsp
       !       write(*,*) 'input,nclsp=',ncoulombsp
    endif
    call MPI_Bcast(ncoulombsp,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    if(ncoulombsp==0) return

    if(myrank==0) read(f_mdff) num
    call MPI_Bcast(num,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)

    allocate(coulombsp_n(npara))
    allocate(coulombsp_atom1(npara,num))
    allocate(coulombsp_atom2(npara,num))
    allocate(coulombsp_charge_product(npara,num))

    if(myrank==0) then
       read(f_mdff) coulombsp_n
       read(f_mdff) coulombsp_atom1
       read(f_mdff) coulombsp_atom2
       read(f_mdff) coulombsp_charge_product
    endif
    call MPI_Bcast(coulombsp_n,npara,MPI_INTEGER4,0, &
         &               MPI_COMM_WORLD,ierr)
    call MPI_Bcast(coulombsp_atom1,npara*num,MPI_INTEGER4,0, &
         &               MPI_COMM_WORLD,ierr)
    call MPI_Bcast(coulombsp_atom2,npara*num,MPI_INTEGER4,0, &
         &               MPI_COMM_WORLD,ierr)
    call MPI_Bcast(coulombsp_charge_product,npara*num,MPI_REAL8,0, &
         &               MPI_COMM_WORLD,ierr)
  end subroutine read_md_coulomb_special
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write coulomb special in mdff.bin File.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_md_coulomb_special
    use device_numbers, only : f_mdff
    implicit none
    integer(4) :: num
    write(f_mdff) ncoulombsp
    if(ncoulombsp==0) return
    num = size(coulombsp_atom1(:,:), 2)
    write(f_mdff) num
    write(f_mdff) coulombsp_n
    write(f_mdff) coulombsp_atom1
    write(f_mdff) coulombsp_atom2
    write(f_mdff) coulombsp_charge_product
  end subroutine write_md_coulomb_special
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write coulomb special.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_md_coulomb_special
    implicit none
    integer(4) :: num
    write(*,*) '[write_md_coulomb_special]'
    write(*,*) ncoulombsp
    if(ncoulombsp==0) return
    num = size(coulombsp_atom1(:,:), 2)
    write(*,*) num
    write(*,*) coulombsp_n
    write(*,*) coulombsp_atom1
    write(*,*) coulombsp_atom2
    write(*,*) coulombsp_charge_product
  end subroutine write_memory_md_coulomb_special

end module special14

