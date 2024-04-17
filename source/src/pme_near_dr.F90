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
!! \brief Module and subroutines to calculate near part of PME.
!<
!----------------------------------------------------------------------
!>
!! \brief Module to calculate near part of PME.
!! \author Yoshimichi Andoh, Jiachao Zhang, Tatsuya Sakashita, 
!!         Noriyuki Yoshii, Shin-ichi Ichikawa
!<
module pme_near_dr

    use precision  !precision module is loaded by default to declare P2P subroutine
    use dr_cntl, only : nbd, nbd2, nadj, slist
    use dr_cntl, only : ldr_zm, ldr_zp, ldr_ym, ldr_yp, ldr_xm, ldr_xp
    use subcell, only : tag, na_per_cell, m2i
    use domain,  only : lxdiv,lydiv,lzdiv
    use domain,  only : ncellx, ncelly, ncellz
    use trajectory_mpi, only : na1cell, nadirect
    use comm_direct2_dr !, only : buffp, rbuffp, buffm, rbuffm
    use fmm_near_dr, only : plist
    use fmm_near_dr, only : wkx, wky, wkz, w3_fx, w3_fy, w3_fz

    use cutoff_radius
    use md_forces
    use lj_mod
    use coulomb_mod
    use md_monitors
    use md_const
    use param
    use atom_virial
    use mpi_tool, only : myrank, mpiend
    use omp_lib ! YA
    use openmp_tool, only : nomp
    use md_condition, only : mdstep

contains

!----------------------------------------------------------------------
  subroutine energy_direct_pme_dr(wkxyz,w3_f)  
!---------------------------------------------------------------------
    use cutoff_radius, only : cutrad,cutrad2
    use ewald_variables
    use table_functions, only :  &
  &       table_density_cut2, &
  &       cl_table,cl_table_delta, &
  &       cl_force_table,cl_force_table_delta, &
  &       lj12_table,lj12_table_delta, &
  &       lj6_table,lj6_table_delta, &
  &       lj12_force_table,lj12_force_table_delta, &
  &       lj6_force_table,lj6_force_table_delta
    implicit none

    real(kind=dp) wkxyz(3,nadirect)
    real(kind=dp) w3_f(3,nadirect,0:nomp-1)

    real(kind=wrp) :: R, Rr6, Rr12, Cij, coef, eps
    real(kind=wrp) :: xi, yi, zi
    real(kind=wrp) :: rx, ry, rz, r2, r2_r
    real(kind=wrp) :: rc, rc_r, rc2_r
    real(kind=wrp) :: tlx, tly, tlz, tcx, tcy, tcz, tmp
    real(kind=dp ) :: stlcx, stlcy, stlcz
    real(kind=wrp) :: chgv_ia0
    real(kind=wrp) :: epsilon_sqrt_ia0, R_half_ia0
    real(kind=dp ) :: fvarx, fvary, fvarz
    real(kind=wrp) :: Ulj12, Ulj6, Ucoulomb
    real(kind=dp ) :: sUlj, sUlj12, sUlj6
    real(kind=dp ) :: sUcoulomb
    integer(4) :: ldr
    integer(4) :: icx0, icy0, icz0
    integer(4) :: jcx1, jcy1, jcz1
    integer(4) :: icx0s, icx0e, icy0s, icy0e
    integer(4) :: jcz1s, jcz1e
    integer(4) :: jcx1t, jcy1t, icz0t
    integer(4) :: ia0, ja1, ia, ja, jas, jabof, jas_eqyx
    integer(4) :: nja
    integer(4) :: isp, jsp
    integer(4) :: nia(0:nomp-1), niat
    integer(4) :: lc, lcs, lce, lcc, lct, ncyc
    integer(4) :: ncarea, ncline
    integer(4) :: iam=0
    integer(4) :: ierr
    include 'mpif.h'
    real(kind=wrp) :: v11,v22,v33
    real(kind=dp ) :: wk_v11, wk_v22, wk_v33
    real(kind=wrp) :: v21,v31,v32
    real(kind=dp ) :: wk_v21, wk_v31, wk_v32
    real(kind=wrp) :: vtr ! for debug
    real(kind=wrp) :: a2,alpsp2,wktrcvir,derfc
    real(kind=wrp) :: R3,R6,R12
    integer(4) :: table_point
    real(kind=wrp) :: table_point_real, interpolate_distance
    real(kind=wrp) :: cl_coef, cl_force_coef
    real(kind=wrp) :: lj12_force_coef, lj6_force_coef

    sUlj   = Cdzero
    sUlj12 = Cdzero
    sUlj6  = Cdzero
    sUcoulomb = Cdzero
    iam = 0
    wk_v11=Cdzero;wk_v22=Cdzero;wk_v33=Cdzero
    wk_v21=Cdzero;wk_v31=Cdzero;wk_v32=Cdzero

    ncline = lzdiv + nbd2
    ncarea = ncline*(lydiv + nbd2)
    lct = ncarea*(lxdiv + nbd)

    a2    = ewald_alpha * ewald_alpha
    alpsp2 = 2.0d0 * ewald_alpha * r_PI_sqrt

!$omp parallel default(none) &
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
!$omp&  shared(tag,na_per_cell) &
!$omp&  shared(lxdiv,lydiv,lzdiv) &
!$omp&  shared(ldr_zm,ldr_zp,ldr_ym,ldr_yp,ldr_xm,ldr_xp) &
!$omp&  shared(m2i,paranum) &
!$omp&  shared(epsilon_sqrt_table,R_half_table,chgv_table) &
!$omp&  shared(epsilon_sqrt,R_half,chgv) &
!$omp&  shared(cutrad2) &
!$omp&  shared(wkxyz,w3_f) &
!$omp&  shared(slist,lct,nia,ncline,ncarea,nomp) &
!$omp&  shared(ewald_alpha,a2,alpsp2) &
!$omp&  private(jcz1,jcy1,jcx1,icz0,icy0,icx0) &
!$omp&  private(jcz1s,jcz1e,icy0s,icy0e,icx0s,icx0e) &
!$omp&  private(jcy1t,jcx1t,icz0t) &
!$omp&  private(ia0,ja1,ia,ja,jas,jabof,jas_eqyx) &
!$omp&  private(ldr,plist,nja) &
!$omp&  private(xi,yi,zi) &
!$omp&  private(epsilon_sqrt_ia0,R_half_ia0,chgv_ia0,isp,jsp) &
!$omp&  private(stlcx,stlcy,stlcz,eps,R,rx,ry,rz,r2,rc2_r) &
!$omp&  private(r2_r,Rr6,Rr12,coef,tlx,tly,tlz,Ulj6,Ulj12) &
!$omp&  private(fvarx,fvary,fvarz) &
!$omp&  private(rc,rc_r,Cij,tmp,tcx,tcy,tcz) &
!$omp&  private(Ucoulomb) &
!$omp&  private(lc,lcs,lce,lcc,ncyc,niat) &
!$omp&  private(iam) &
!$omp&  private(v11,v22,v33) &
!$omp&  private(wkx,wky,wkz,w3_fx,w3_fy,w3_fz) &
!$omp&  reduction(+:wk_v11,wk_v22,wk_v33) &
!$omp&  private(v21,v31,v32) &
!$omp&  reduction(+:wk_v21,wk_v31,wk_v32) &
!$omp&  reduction(+:sUlj,sUlj6,sUlj12,sUcoulomb)
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
            R_half_table(ja,iam) = R_half(paranum(jsp))
            chgv_table(ja,iam) = md_QQ_4PiE*chgv(paranum(jsp))
            wkx(ja) = wkxyz(1,ja1)  ! add for no AOS access in inner-most loop.
            wky(ja) = wkxyz(2,ja1)
            wkz(ja) = wkxyz(3,ja1)
            w3_fx(ja) = Cdzero      ! add for no AOS access in inner-most loop.
            w3_fy(ja) = Cdzero      ! Cdzero is a zero value parameter
            w3_fz(ja) = Cdzero      ! in precision_mod.
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
          R_half_ia0=R_half(paranum(isp))
          chgv_ia0=chgv(paranum(isp))
          xi=wkxyz(1,ia0)
          yi=wkxyz(2,ia0)
          zi=wkxyz(3,ia0)

          stlcx=Cdzero
          stlcy=Cdzero
          stlcz=Cdzero

          if(jas_eqyx <= ia0 .and. jcz1 > 0) jas = ia0 + 1

! local pair-list.
          nja=0
          do ja1 = jas, tag(jcz1e,jcy1,jcx1) + &
            &           na_per_cell(jcz1e,jcy1,jcx1)-1
             rx=xi-wkx(ja1-jabof)
             ry=yi-wky(ja1-jabof)
             rz=zi-wkz(ja1-jabof)
             r2=rx*rx+ry*ry+rz*rz
             if(r2 > cutrad2)then
               cycle
             endif
             nja=nja+1
             plist(nja)=ja1
          end do

!!ocl striping
!         do ja1 = jas, tag(jcz1e,jcy1,jcx1) + &
!            &           na_per_cell(jcz1e,jcy1,jcx1)-1
!ocl norecurrence
          do ja = 1, nja
             ja1=plist(ja)
             rx=xi-wkx(ja1-jabof)
             ry=yi-wky(ja1-jabof)
             rz=zi-wkz(ja1-jabof)
             r2=rx*rx+ry*ry+rz*rz
             r2_r=Cone/r2
             table_point_real     = table_density_cut2*r2_r
             table_point          = int(table_point_real)
             interpolate_distance = table_point_real - dble(table_point)
             !------ Lennard-Jones part start
             ! ^^^ spherical cut-off ^^^
             fvarx=w3_fx(ja1-jabof)
             fvary=w3_fy(ja1-jabof)
             fvarz=w3_fz(ja1-jabof)
!
             eps=epsilon_sqrt_ia0* &
             &   epsilon_sqrt_table(ja1-jabof,iam)
#if defined(OPLSAMBER) && !defined(GAFF)
             R=sqrt(4d0*R_half_ia0*R_half_table(ja1-jabof,iam))
#else /** CHARMM **/
             R=R_half_ia0+R_half_table(ja1-jabof,iam)
#endif
#ifndef NOTABLE
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
#else /** -DNOTABLE **/
             Rr6=R * R * r2_r
             Rr6=Rr6 * Rr6 * Rr6
             Rr12=Rr6 * Rr6
             coef=Ctwelve * eps * r2_r * (Rr12-Rr6)
             ! ^^^ potential ^^^
             Ulj12=eps*Rr12
             Ulj6 =-Ctwo*eps*Rr6
#endif
             tlx=coef*rx
             tly=coef*ry
             tlz=coef*rz
             sUlj  =sUlj  +(Ulj12+Ulj6)
             sUlj12=sUlj12+Ulj12
             sUlj6 =sUlj6 +Ulj6
             ! ^^^ force ^^^
             stlcx=stlcx+tlx
             stlcy=stlcy+tly
             stlcz=stlcz+tlz
             fvarx = fvarx - tlx
             fvary = fvary - tly
             fvarz = fvarz - tlz
             !------ Coulomb part start
             Cij=chgv_ia0*chgv_table(ja1-jabof,iam)
#ifndef NOTABLE
             !potential
             Ucoulomb = Cij &
                  &         * ( cl_table(table_point) &
                  &            +interpolate_distance  &
                  &            *cl_table_delta(table_point))
             !force
             cl_force_coef = cl_force_table(table_point) &
                  &            +interpolate_distance     &
                  &            *cl_force_table_delta(table_point)
             cl_force_coef = cl_force_coef * Cij
             tcx = cl_force_coef * rx
             tcy = cl_force_coef * ry
             tcz = cl_force_coef * rz
#else /** -DNOTABLE **/
             rc =sqrt(r2)
             rc_r=sqrt(r2_r)
             Ucoulomb=Cij*derfc(ewald_alpha*rc)*rc_r
             tmp=Cij*alpsp2*exp(-a2*r2)+Ucoulomb
             tmp=tmp*r2_r
             tcx=tmp*rx
             tcy=tmp*ry
             tcz=tmp*rz
#endif 
             sUcoulomb=sUcoulomb+Ucoulomb

             stlcx=stlcx+tcx
             stlcy=stlcy+tcy
             stlcz=stlcz+tcz
             w3_fx(ja1-jabof)=fvarx - tcx
             w3_fy(ja1-jabof)=fvary - tcy
             w3_fz(ja1-jabof)=fvarz - tcz
             v11 = rx*(tlx+tcx)
             v22 = ry*(tly+tcy)
             v33 = rz*(tlz+tcz)
             wk_v11 = wk_v11 + v11
             wk_v22 = wk_v22 + v22
             wk_v33 = wk_v33 + v33
             v21 = ry*(tlx+tcx)
             v31 = rz*(tlx+tcx)
             v32 = rz*(tly+tcy)
             wk_v21 = wk_v21 + v21
             wk_v31 = wk_v31 + v31
             wk_v32 = wk_v32 + v32
          end do  ! ja1.
          w3_f(1,ia0,iam)=w3_f(1,ia0,iam)+stlcx
          w3_f(2,ia0,iam)=w3_f(2,ia0,iam)+stlcy
          w3_f(3,ia0,iam)=w3_f(3,ia0,iam)+stlcz
        end do  ! ia0.
!$omp end do
       !
       !+++++++++ Lennard-Jones and Coulomb force calculation ended ++++++++++
       !++++++++++++  Main calculation ended

        end do  ! jcz1.
        jabof = tag(1-nbd,jcy1,jcx1) - 1
        do ja1 = jabof + 1, tag(lzdiv+nbd,jcy1,jcx1) + &
            &                 na_per_cell(lzdiv+nbd,jcy1,jcx1) - 1
           w3_f(1,ja1,iam)=w3_f(1,ja1,iam)+w3_fx(ja1-jabof)
           w3_f(2,ja1,iam)=w3_f(2,ja1,iam)+w3_fy(ja1-jabof)
           w3_f(3,ja1,iam)=w3_f(3,ja1,iam)+w3_fz(ja1-jabof)
        end do  ! ja1
      end do  ! jcy1.
    end do  ! jcx1.
!$omp end parallel

    wk_p_energy = wk_p_energy + sUlj + sUcoulomb
    wk_vir2(1,0)=wk_vir2(1,0)+wk_v11
    wk_vir2(2,0)=wk_vir2(2,0)+wk_v22
    wk_vir2(3,0)=wk_vir2(3,0)+wk_v33
    wk_vir2(4,0)=wk_vir2(4,0)+wk_v21
    wk_vir2(5,0)=wk_vir2(5,0)+wk_v31
    wk_vir2(6,0)=wk_vir2(6,0)+wk_v32
#ifdef DEBUGFCE
    wplj=sUlj
    wpcl=sUcoulomb
#endif
  end subroutine energy_direct_pme_dr

end module pme_near_dr
