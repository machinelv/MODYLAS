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
!! \brief Module and subroutines to calculate near part of FMM.
!<
!----------------------------------------------------------------------
!>
!! \brief Module to calculate far part of FMM.
!! \author Yoshimichi Andoh, Jiachao Zhang, Tatsuya Sakashita, 
!!         Noriyuki Yoshii, Shin-ichi Ichikawa
!<
module fmm_near_dr
    use precision
    use dr_cntl, only : nbd, nbd2, nadj, slist, &
                        ldr_zm, ldr_zp, ldr_ym, ldr_yp, ldr_xm, ldr_xp
    use subcell, only : tag, na_per_cell, m2i
    use domain,  only : lxdiv,lydiv,lzdiv, &
                        ncellx, ncelly, ncellz
    use trajectory_mpi, only : na1cell, nadirect

    use cutoff_radius
    use md_forces
    use lj_mod
    use coulomb_mod
    use md_monitors
    use md_const
    use param
    use atom_virial
    use mpi_tool, only : myrank, mpiend
    use openmp_tool, only : nomp
    use md_condition, only : mdstep

    implicit none
    integer(4), allocatable :: plist(:)    ! add for pme local pair-list
    !In order to remove non unit stride memory access.
    real(kind=wrp), allocatable :: wkx(:), wky(:), wkz(:)
    real(kind=wrap), allocatable :: w3_fx(:), w3_fy(:), w3_fz(:)
    real(kind=wrp), allocatable :: flg(:)

contains

!----------------------------------------------------------------------
!>
!! \brief Subroutine to allocate arrays .
!! \author Yoshimichi Andoh, Jiachao Zhang, Shin-ichi Ichikawa
!<
  subroutine init_energy_direct_dr()
!----------------------------------------------------------------------
   implicit none

   allocate(slist(1:na1cell*((2*nadj+1)*nadj+nadj+1), 0:nomp-1))
   allocate(plist(na1cell*(lzdiv+nbd2)) )
   slist=0
   plist=0

   allocate(wkx(na1cell*(lzdiv+nbd2)) )
   allocate(wky(na1cell*(lzdiv+nbd2)) )
   allocate(wkz(na1cell*(lzdiv+nbd2)) )
   wkx=Czero
   wky=Czero
   wkz=Czero
   allocate(flg(na1cell*(lzdiv+nbd2)) )
   allocate(w3_fx(na1cell*(lzdiv+nbd2)) )
   allocate(w3_fy(na1cell*(lzdiv+nbd2)) )
   allocate(w3_fz(na1cell*(lzdiv+nbd2)) )
   w3_fx=Cdzero
   w3_fy=Cdzero
   w3_fz=Cdzero

  end subroutine init_energy_direct_dr

!----------------------------------------------------------------------
!>
!! \brief Subroutine to allocate arrays .
!! \author Yoshimichi Andoh, Jiachao Zhang, Shin-ichi Ichikawa
!<
  subroutine energy_direct_dr(wkxyz,w3_f)  
!---------------------------------------------------------------------
#ifdef TABLE
    use table_functions, only :  &
  &       table_density_cut2, &
  &       cl_table,cl_table_delta, &
  &       cl_force_table,cl_force_table_delta, &
  &       lj12_table,lj12_table_delta, &
  &       lj6_table,lj6_table_delta, &
  &       lj12_force_table,lj12_force_table_delta, &
  &       lj6_force_table,lj6_force_table_delta
#endif
    use void123, only : void_topL,void_numL,void_list_j0
    use subcell, only : i2m
#ifdef ONEPROC_AXIS
    use boundary, only : pbc_pair_dr
#endif
    implicit none

    ! always DP
    real(kind=dp) wkxyz(3,nadirect)
    real(kind=dp) w3_f(3,nadirect,0:nomp-1)
    !Intermediate values for main computing in single-precision
    real(kind=wrp) :: R, Rr6, Rr12, Cij, coef, eps
    real(kind=wrp) :: xi, yi, zi
    real(kind=wrp) :: rx, ry, rz, r2, r2_r
    real(kind=wrp) :: rc, rc_r, rc2_r
    real(kind=wrp) :: tlx, tly, tlz, tcx, tcy, tcz, tmp
    real(kind=wrap) :: stlcx, stlcy, stlcz
    real(kind=wrp) :: chgv_ia0
    real(kind=wrp) :: epsilon_sqrt_ia0, R_half_ia0
    real(kind=wrp) :: Ulj12, Ulj6, Ucoulomb
    real(kind=wrp) :: tmp_QQ4PiE, tmp_cutrad2
    !Double precision for accumulation.
    real(kind=wrp) :: fvarx, fvary, fvarz    
    real(kind=wrap) :: sUlj12, sUlj6
    real(kind=wrap) :: sUcoulomb

    integer(4) :: ldr
    integer(4) :: icx0, icy0, icz0
    integer(4) :: jcx1, jcy1, jcz1
    integer(4) :: icx0s, icx0e, icy0s, icy0e
    integer(4) :: jcz1s, jcz1e
    integer(4) :: jcx1t, jcy1t, icz0t
    integer(4) :: ia0, ja1, ia, ja, jas, jabof, jas_eqyx
    integer(4) :: isp, jsp
    integer(4) :: nia(0:nomp-1), niat
    integer(4) :: lc, lcs, lce, lcc, lct, ncyc
    integer(4) :: ncarea, ncline
    integer(4) :: iam=0
    integer(4) :: ierr
    include 'mpif.h'

#ifdef TABLE
    real(kind=wrp) :: R3,R6,R12
    integer(kind=wip) :: table_point
    real(kind=wrp) :: table_point_real, interpolate_distance
    real(kind=wrp) :: cl_coef, cl_force_coef
    real(kind=wrp) :: lj12_force_coef, lj6_force_coef
#endif

    call update_local_voidlist()

    tmp_QQ4PiE=md_QQ_4PiE ! type is converted implicitly
    tmp_cutrad2=cutrad2   ! type is converted implicitly

    sUlj12 = 0d0
    sUlj6  = 0d0
    sUcoulomb = 0d0
    iam = 0

    ncline = lzdiv + nbd2
    ncarea = ncline*(lydiv + nbd2)
    lct = ncarea*(lxdiv + nbd)

!$omp parallel default(none) &
#ifdef TABLE
!$omp&  private(R3,R6,R12) &
!$omp&  private(table_point_real,table_point) &
!$omp&  private(interpolate_distance) &
!$omp&  private(cl_coef,cl_force_coef) &
!$omp&  private(lj12_force_coef,lj6_force_coef) &
!$omp&  shared(table_density_cut2) &
!$omp&  shared(lj12_table,lj6_table) &
!$omp&  shared(lj12_table_delta,lj6_table_delta) &
!$omp&  shared(lj12_force_table,lj6_force_table) &
!$omp&  shared(lj12_force_table_delta,lj6_force_table_delta) &
!$omp&  shared(cl_table,cl_table_delta) &
!$omp&  shared(cl_force_table,cl_force_table_delta) &
#endif
!$omp&  shared(tag,na_per_cell) &
!$omp&  shared(lxdiv,lydiv,lzdiv) &
!$omp&  shared(ldr_zm,ldr_zp,ldr_ym,ldr_yp,ldr_xm,ldr_xp) &
!$omp&  shared(m2i,paranum) &
!$omp&  shared(epsilon_sqrt_table,R_half_table,chgv_table) &
!$omp&  shared(epsilon_sqrt,R_half,chgv) &
!$omp&  shared(tmp_cutrad2,tmp_QQ4PiE) &
!$omp&  shared(wkxyz,w3_f) &
!$omp&  shared(slist,lct,nia,ncline,ncarea,nomp) &
!$omp&  private(jcz1,jcy1,jcx1,icz0,icy0,icx0) &
!$omp&  private(jcz1s,jcz1e,icy0s,icy0e,icx0s,icx0e) &
!$omp&  private(jcy1t,jcx1t,icz0t) &
!$omp&  private(ia0,ja1,ia,ja,jas,jabof,jas_eqyx) &
!$omp&  private(ldr) &
!$omp&  private(xi,yi,zi) &
!$omp&  private(wkx, wky, wkz) &          ! add for prior data type conversion & No AOS access.
!$omp&  private(w3_fx, w3_fy, w3_fz) &    ! add for No AOS access.
!$omp&  private(epsilon_sqrt_ia0,R_half_ia0,chgv_ia0,isp,jsp) &
!$omp&  private(stlcx,stlcy,stlcz,eps,R,rx,ry,rz,r2,rc2_r) &
!$omp&  private(r2_r,Rr6,Rr12,coef,tlx,tly,tlz,Ulj6,Ulj12) &
!$omp&  private(fvarx,fvary,fvarz) &
!$omp&  private(rc,rc_r,Cij,tmp,tcx,tcy,tcz) &
!$omp&  private(Ucoulomb) &
!$omp&  private(lc,lcs,lce,lcc,ncyc,niat) &
!$omp&  private(iam) &
!$omp&  shared(void_topL,void_numL,void_list_j0) &
!$omp&  private(flg) &
!$omp&  reduction(+:sUlj6,sUcoulomb) &
!$omp&  reduction(+:sUlj12)
!$  iam=omp_get_thread_num()

    ncyc = 0
    lcs = 0
    lce = 0
    do jcx1 = 1, lxdiv + nbd
      do jcy1 = 1 - nbd, lydiv + nbd

          jabof = tag(1-nbd,jcy1,jcx1) - 1
         !++++++++++++  Make table for [eps,R,chgv] start
          do ja1 = jabof + 1, tag(lzdiv+nbd,jcy1,jcx1) + &
            &                 na_per_cell(lzdiv+nbd,jcy1,jcx1) - 1
            jsp = m2i(ja1)
            ja = ja1 - jabof
            epsilon_sqrt_table(ja,iam) = epsilon_sqrt(paranum(jsp))
            R_half_table(ja,iam)       = R_half(paranum(jsp))
            chgv_table(ja,iam)         = tmp_QQ4PiE*chgv(paranum(jsp))
            !wkxyz->wkx,wky,wkz
            wkx(ja) = wkxyz(1,ja1)    ! add for prior data type conversion & No AOS access.
            wky(ja) = wkxyz(2,ja1)
            wkz(ja) = wkxyz(3,ja1)
            w3_fx(ja) = Cdzero       ! add for No Array-Of-Structure(AOS) access.
            w3_fy(ja) = Cdzero
            w3_fz(ja) = Cdzero
          enddo
         !++++++++++++  Make table for [eps,R,chgv] ended

        do jcz1 = 1 - nbd, lzdiv + nbd

          icz0 = jcz1
          lcc = jcz1 + nbd + ncline*(jcy1 + nbd - 1) + ncarea*(jcx1 - 1)
        if(lcc > lce) then
          lcs = ncyc*nomp + 1
          lce = (ncyc + 1)*nomp
          if(lce > lct) lce = lct
          ncyc = ncyc + 1
!$omp do
        do lc = lcs, lce
          jcx1t = (lc - 1)/ncarea + 1
          jcy1t = (lc - ncarea*(jcx1t - 1) - 1)/ncline + 1 - nbd
          icz0t = lc - ncarea*(jcx1t - 1) - ncline*(jcy1t + nbd - 1) - nbd

          niat = 0
          icx0s = jcx1t - nadj
          if(icx0s < 1-nbd) icx0s = 1 - nbd
        ! exclude boundary-boundary.
          if(jcx1t == 1-nbd .and. icx0s == 1-nbd  ) icx0s = 1
          icx0e = jcx1t
        ! exclude boundary-boundary.
          if(jcx1t == lxdiv+nbd .and. icx0e == lxdiv+nbd) icx0e = lxdiv
        do icx0 = icx0s, icx0e
          icy0s = jcy1t - nadj
          if(icy0s < 1-nbd) icy0s = 1 - nbd
        ! exclude boundary-boundary.
          if(jcy1t == 1-nbd .and. icy0s == 1-nbd  ) icy0s = 1
          icy0e = jcy1t + nadj
          if(icy0e > lydiv+nbd) icy0e = lydiv + nbd
          if(icx0 == jcx1t ) icy0e = jcy1t
        ! exclude boundary-boundary.
          if(jcy1t == lydiv+nbd .and. icy0e == lydiv+nbd) icy0e = lydiv
        do icy0 = icy0s, icy0e

        ! yx-boundary arbitration.
        ! ldr_ym(1-nbd:lxdiv+nbd,-nadj:nadj)  ; -1 -> down process, +1 -> upper process.
        ! ldr_yp(1-nbd:lxdiv+nbd,-nadj:nadj)  ; -1 -> down process, +1 -> upper process.
        ! ldr_xm(1-nbd:lydiv+nbd,-nadj:nadj)  ; -1 -> down process, +1 -> upper process.
        ! ldr_xp(1-nbd:lydiv+nbd,-nadj:nadj)  ; -1 -> down process, +1 -> upper process.

          if(jcx1t == 1 .and. icx0 == 1-nbd) then
            ldr = ldr_xm(icy0,icy0-jcy1t)
            if(ldr /=  1) cycle
          endif
          if(jcx1t == lxdiv+nbd .and. icx0 == lxdiv) then
            ldr = ldr_xp(jcy1t,jcy1t-icy0)
            if(ldr /= -1) cycle
          endif
          if(jcy1t == 1-nbd .and. icy0 == 1) then
            ldr = ldr_ym(jcx1t,jcx1t-icx0)
            if(ldr /=  1) cycle
          endif
          if(jcy1t == 1 .and. icy0 == 1-nbd) then
            ldr = ldr_ym(icx0,icx0-jcx1t)
            if(ldr /=  1) cycle
          endif
          if(jcy1t == lydiv+nbd .and. icy0 == lydiv) then
            ldr = ldr_yp(jcx1t,jcx1t-icx0)
            if(ldr /= -1) cycle
          endif
          if(jcy1t == lydiv .and. icy0 == lydiv+nbd) then
            ldr = ldr_yp(icx0,icx0-jcx1t)
            if(ldr /= -1) cycle
          endif

        ! self atom list.
          do ia0 = tag(icz0t,icy0,icx0), &
                   tag(icz0t,icy0,icx0)+na_per_cell(icz0t,icy0,icx0)-1
            niat = niat + 1
            slist(niat,lc - lcs) = ia0
          end do  ! ia0.

        end do  ! icy0.
        end do  ! icx0.
          nia(lc - lcs) = niat
        end do  ! lc.
!$omp end do
        end if  ! lcc > lce.

          jas_eqyx = tag(jcz1,jcy1,jcx1)
          jcz1s = jcz1 - nadj
          if(jcz1s < 1-nbd) jcz1s = 1 - nbd
          jas = tag(jcz1s,jcy1,jcx1)
          jcz1e = jcz1 + nadj
          if(jcz1e > lzdiv+nbd) jcz1e = lzdiv + nbd
        ! z-boundary arbitration.
        ! ldr_zm(1-nbd:lydiv+nbd,1-nbd:lxdiv+nbd) ; -1 -> down process, +1 -> upper process.
        ! ldr_zp(1-nbd:lydiv+nbd,1-nbd:lxdiv+nbd) ; -1 -> down process, +1 -> upper process.
          
          if(icz0 == 1-nbd) then
            ldr = ldr_zm(jcy1,jcx1)
            if(ldr ==  1) jcz1s = 1
            if(ldr == -1) jcz1s = 1+nbd
            jas = tag(jcz1s,jcy1,jcx1)
          endif
          if(icz0 == 1) then
            ldr = ldr_zm(jcy1,jcx1)
            if(ldr ==  1) jcz1s = 1-nbd
            if(ldr == -1) jcz1s = 1
            jas = tag(jcz1s,jcy1,jcx1)
          endif
          if(icz0 == lzdiv+nbd) then
            ldr = ldr_zp(jcy1,jcx1)
            if(ldr == -1) jcz1e = lzdiv
            if(ldr ==  1) jcz1e = lzdiv-nbd
          endif
          if(icz0 == lzdiv) then
            ldr = ldr_zp(jcy1,jcx1)
            if(ldr == -1) jcz1e = lzdiv+nbd
            if(ldr ==  1) jcz1e = lzdiv
          endif

       !++++++++++++  Main calculation start
       !
       !+++++++++ Lennard-Jones and Coulomb force calculation start ++++++++++
       !
!$omp do schedule(static,1)
        do ia = 1, nia(lcc - lcs)
          ia0 = slist(ia,lcc - lcs)
          isp=m2i(ia0)
          
          epsilon_sqrt_ia0=epsilon_sqrt(paranum(isp))
          R_half_ia0      =R_half(paranum(isp))
          chgv_ia0        =chgv(paranum(isp))
          xi=wkxyz(1,ia0)
          yi=wkxyz(2,ia0)
          zi=wkxyz(3,ia0)

          stlcx=Cdzero
          stlcy=Cdzero
          stlcz=Cdzero

          if(jas_eqyx <= ia0 .and. jcz1 > 0) jas = ia0 + 1

! ya+: create flg for 1-2, 1-3 interactions, common for LJ and Coulomb
          do ja1 = jas, tag(jcz1e,jcy1,jcx1)+na_per_cell(jcz1e,jcy1,jcx1)-1
            ja = ja1 - jabof      ! local number for temporary array
            flg(ja)=Cone          ! set to initial value 1
          enddo
!      overwrite flg
          do jsp=void_topL(ia0),void_topL(ia0)+void_numL(ia0)-1
            ja1 = void_list_j0(jsp)
            if(ja1.ge.jas .and.  &  ! within access range of the next ja1 loop
               ja1.le.tag(jcz1e,jcy1,jcx1)+na_per_cell(jcz1e,jcy1,jcx1)-1)then
              ja = ja1 - jabof
              flg(ja)=Czero
            endif
          enddo

#if defined (PRECISION_P2P_SP)
!$omp simd simdlen(16)
#else
!$omp simd simdlen(8)
#endif
          do ja1 = jas, tag(jcz1e,jcy1,jcx1) + &
            &           na_per_cell(jcz1e,jcy1,jcx1)-1
             rx=xi-wkx(ja1-jabof) !! wkxyz is always double.
             ry=yi-wky(ja1-jabof)
             rz=zi-wkz(ja1-jabof)
#ifdef ONEPROC_AXIS
             call  pbc_pair_dr(rx,ry,rz)
#endif
             r2=rx*rx+ry*ry+rz*rz
             r2_r=Cone/r2
#ifdef TABLE
             table_point_real     = table_density_cut2*r2_r
             table_point          = int(table_point_real)
             interpolate_distance = table_point_real - dble(table_point)
#endif
             !------ Lennard-Jones part start
             ! ^^^ spherical cut-off ^^^
             if(r2<=tmp_cutrad2) then
               eps=epsilon_sqrt_ia0* &
                &   epsilon_sqrt_table(ja1-jabof,iam)
             else
                eps=Czero
             endif

             eps=eps*flg(ja1-jabof)

#if defined(OPLSAMBER) && !defined(GAFF)
             R=sqrt(Cfour*R_half_ia0*R_half_table(ja1-jabof,iam))
#else /** CHARMM, GAFF **/
             R=R_half_ia0+R_half_table(ja1-jabof,iam)
#endif
#ifdef TABLE
             R3=R**3
             R6=R3**2
             R12=R6**2
             ! ^^^ potential ^^^   NOTE: CHARMMFSW is included in tables
             Ulj12 = lj12_table(table_point) &
                  &      + interpolate_distance    &
                  &      * lj12_table_delta(table_point)
             Ulj6  = lj6_table(table_point)  &
                  &      + interpolate_distance    &
                  &      * lj6_table_delta(table_point)
             Ulj12 =        Ulj12* eps * R12
             Ulj6  = -Ctwo * Ulj6 * eps * R6
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
             coef = Ctwelve * eps * (lj12_force_coef - lj6_force_coef)
#else /**not TABLE**/
             Rr6=R * R * r2_r
             Rr6=Rr6 * Rr6 * Rr6
             Rr12=Rr6 * Rr6
             coef=Ctwelve * eps * r2_r * (Rr12-Rr6)
             ! ^^^ potential ^^^
             Ulj12=eps*Rr12
             Ulj6 =-Ctwo*eps*Rr6
#endif /**TABLE**/
             tlx=coef*rx
             tly=coef*ry
             tlz=coef*rz
             sUlj12 =sUlj12 +Ulj12
             sUlj6  =sUlj6  +Ulj6
             !------ Coulomb part start
             rc_r=sqrt(r2_r)
             Cij=chgv_ia0*chgv_table(ja1-jabof,iam)
             Cij=Cij*flg(ja1-jabof)
             Cij=Cij*rc_r
             tmp=Cij*r2_r
             tcx=tmp*rx
             tcy=tmp*ry
             tcz=tmp*rz
             fvarx = tlx + tcx ! add to perform intra-atom addition with single precision.
             fvary = tly + tcy
             fvarz = tlz + tcz
             Ucoulomb=Cij
             sUcoulomb=sUcoulomb+Ucoulomb

             stlcx=stlcx+fvarx
             stlcy=stlcy+fvary
             stlcz=stlcz+fvarz
             w3_fx(ja1-jabof)=w3_fx(ja1-jabof) - fvarx
             w3_fy(ja1-jabof)=w3_fy(ja1-jabof) - fvary
             w3_fz(ja1-jabof)=w3_fz(ja1-jabof) - fvarz
          end do  ! ja1.
          w3_f(1,ia0,iam)=w3_f(1,ia0,iam)+stlcx
          w3_f(2,ia0,iam)=w3_f(2,ia0,iam)+stlcy
          w3_f(3,ia0,iam)=w3_f(3,ia0,iam)+stlcz
        end do  ! ia0.
!$omp end do
       !+++++++++ Lennard-Jones and Coulomb force calculation ended ++++++++++
       !++++++++++++  Main calculation ended

        end do  ! jcz1.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        jabof = tag(1-nbd,jcy1,jcx1) - 1
        do ja1 = jabof + 1, tag(lzdiv+nbd,jcy1,jcx1) + &
            &                 na_per_cell(lzdiv+nbd,jcy1,jcx1) - 1
           w3_f(1,ja1,iam)=w3_f(1,ja1,iam)+w3_fx(ja1-jabof)
           w3_f(2,ja1,iam)=w3_f(2,ja1,iam)+w3_fy(ja1-jabof)
           w3_f(3,ja1,iam)=w3_f(3,ja1,iam)+w3_fz(ja1-jabof)
        end do  ! ja1
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do  ! jcy1.
    end do  ! jcx1.
!$omp end parallel

    wk_p_energy = wk_p_energy + sUlj12 + sUlj6 + sUcoulomb

    wk_vir2(1,0)=wk_vir2(1,0)+(Ctwelve*sUlj12+Csix*sUlj6+sUcoulomb)/Cdthree
    wk_vir2(2,0)=wk_vir2(2,0)+(Ctwelve*sUlj12+Csix*sUlj6+sUcoulomb)/Cdthree
    wk_vir2(3,0)=wk_vir2(3,0)+(Ctwelve*sUlj12+Csix*sUlj6+sUcoulomb)/Cdthree

#ifdef DEBUGFCE
    wplj=sUlj12+sUlj6
    wpcl=sUcoulomb
#endif
  
  end subroutine energy_direct_dr

  
!----------------------------------------------------------------------
!>
!! \brief Subroutine to make local voidlist of j0 for each i0 atom.
!! \author Yoshimichi Andoh
!<
  subroutine update_local_voidlist()
!----------------------------------------------------------------------
  use subcell, only : i2m,m2i,tag,na_per_cell
  use void123, only : void_topL,void_numL,void_list_j0, maxj0,  &
                      void_n, void_atom2, nvoid, maxj0_per_i0
  use mpi_tool, only : myrank, mpiend, modylas_abort
  use trajectory_mpi, only : na1cell, narea, naline
  use openmp_tool, only : nomp
  implicit none
  integer(4) icx0,icy0,icz0,icx1,icy1,icz1
  integer(4) i,ipar,j,jpar,i0,j0,isp
  integer(4) li,mi,mi0,nv,ichk(0:nomp-1),iam

  if(nvoid == 0) return  !! note void_numL=0 (default)

  iam=0
  ichk=0

!$omp  parallel default(none) &
!$omp& private(iam,icx0,icy0,icz0,icx1,icy1,icz1,i0,isp) &
!$omp& private(mi,mi0,nv,i,ipar,j,j0) &
!$omp& shared(lxdiv,lydiv,lzdiv) &
!$omp& shared(narea,naline,na1cell) &
!$omp& shared(tag,na_per_cell,m2i,i2m,paranum) &
!$omp& shared(void_n,void_atom2) &
!$omp& shared(void_numL,void_topL,void_list_j0) &
!$omp& shared(ichk)
!$ iam=omp_get_thread_num()
!$omp do collapse(2)
  do icx0=1-nbd,lxdiv+nbd
  do icy0=1-nbd,lydiv+nbd
    icx1=icx0-(1-nbd)
    icy1=icy0-(1-nbd)
    icz1=0
    mi=maxj0_per_i0*(narea*icx1+naline*icy1+na1cell*icz1)+1
    mi0=mi
    do i0=tag(1-nbd,icy0,icx0), &
   &      tag(lzdiv+nbd,icy0,icx0)+na_per_cell(lzdiv+nbd,icy0,icx0)-1
      i   =m2i(i0)
      ipar=paranum(i)
      void_topL(i0)=mi
!
      nv=0
      do isp=1,void_n(ipar)
        j    = void_atom2(ipar,isp)+(i-ipar)
        j0   = i2m(j)
        if(j0 > 0)then
          void_list_j0(mi)=j0
          mi=mi+1
          nv=nv+1
        endif
      enddo ! isp
!
      void_numL(i0)=nv
!
    enddo ! i0
    if(mi-mi0.gt.maxj0_per_i0*naline) ichk(iam)=+1
  enddo ! icy0
  enddo ! icx0
!$omp end do nowait
!$omp end parallel

  !check for overflow [essential for correct PotE result]
  do iam=0,nomp-1
    if(ichk(iam) /= 0)then
    if(myrank==0) write(0,*) 'ERROR: mi > maxj0. Set larger maxj0_per_i0.'
    call modylas_abort
    endif
  enddo

  end subroutine update_local_voidlist

end module fmm_near_dr
