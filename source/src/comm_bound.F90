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
!! \brief  Module and subroutines which relate to MPI communications 
!!         for passing/receiving atom coordinates across process boundary.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module which relates to MPI communications for atom coordinates.
!! \author Yoshimichi Andoh, Shin-ichi Ichikawa
!<
module comm_bound_mod
  use omp_lib
  use subcell
  use mpi_tool
  implicit none

  integer nczdiv, ncydiv, ncxdiv
  integer max_mvatom      ! max number of atoms moving to the cell.
  integer max_mvseg       ! max number of segments moving to the cell.
  integer max_cellcbd     ! max number of cells on communication buffer.
  real(8),allocatable :: abucket(:,:,:,:,:)
  integer,allocatable :: iabucket(:,:,:,:)
  integer,allocatable :: isbucket(:,:,:,:)
  integer,allocatable :: ncseg(:,:,:)
  integer,allocatable :: ncatom(:,:,:)
  real(8),allocatable :: buffp(:,:)
  real(8),allocatable :: buffm(:,:)
  integer,allocatable :: ibuffp(:)
  integer,allocatable :: ibuffm(:)
  integer,allocatable :: isbufp(:)
  integer,allocatable :: isbufm(:)
  real(8),allocatable :: rbuff_p(:,:)
  real(8),allocatable :: rbuff_m(:,:)
  integer,allocatable :: irbuff_p(:)
  integer,allocatable :: irbuff_m(:)
  integer,allocatable :: irsbuf_p(:)
  integer,allocatable :: irsbuf_m(:)
  integer,allocatable :: ncatmw(:,:,:,:)
  integer :: MSG_SEND,MSG_RECV

contains

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to enter zero into work arrays
!! \author Yoshimichi Andoh, Shin-ichi Ichikawa
!<
  subroutine init_comm_buffer()
!----------------------------------------------------------------------
    use trajectory_org        ! n
    use trajectory_mpi        ! n
    use atom_virial      ! vir
    use forces    ! f, wk_k
    use md_monitors  ! p_energy, wk_p_energy
    use domain, only : lxdiv, lydiv, lzdiv
    use dr_cntl, only : nbd
    implicit none
    integer(4) :: ii,i0,iam,nomp
    integer(4) :: icx0,icy0,icz0

    nomp = 1
!$  nomp = omp_get_max_threads()

    wk_vir  = 0d0
    wk_p_energy = 0d0

!$omp parallel default(shared) &
!$omp& private(iam,i0,ii) &
!$omp& private(icx0,icy0,icz0)
!$omp do
    do iam = 0,nomp-1
       wk_vir2(:,iam) = 0d0
!default: MTD
       do icx0=1-nbd,lxdiv+nbd; do icy0=1-nbd,lydiv+nbd
       do i0=tag(    1-nbd,icy0,icx0),   &
     &       tag(lzdiv+nbd,icy0,icx0)+na_per_cell(lzdiv+nbd,icy0,icx0)-1
          w3_f(1,i0,iam) = 0.0d0
          w3_f(2,i0,iam) = 0.0d0
          w3_f(3,i0,iam) = 0.0d0
       enddo ! i0
       enddo; enddo ! icy0,icx0
    end do ! iam
!$omp end do
!$omp end parallel

  end subroutine init_comm_buffer
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to initialize comm_bound
!! \author Yoshimichi Andoh, Shin-ichi Ichikawa
!<
  subroutine init_comm_bound()
!----------------------------------------------------------------------
    use mpi_3d_grid
    use trajectory_mpi
    use domain, only : lxdiv, lydiv, lzdiv
    implicit none
    integer(4) :: itmp

    nczdiv = lzdiv
    ncydiv = lydiv
    ncxdiv = lxdiv

    max_mvatom = max(na1cell*3,50)
    max_mvseg = max_mvatom

    max_cellcbd = max(nczdiv*ncydiv, ncydiv*ncxdiv)
    max_cellcbd = max(max_cellcbd, nczdiv*ncxdiv)
    itmp = max(ncydiv+ncxdiv, ncxdiv+nczdiv)
    itmp = max(itmp, nczdiv+ncydiv)
    max_cellcbd = max_cellcbd + 2*itmp + 4

    allocate(abucket(6, max_mvatom, nczdiv+2, ncydiv+2, ncxdiv+2))
    allocate(  iabucket(max_mvatom, nczdiv+2, ncydiv+2, ncxdiv+2))
    allocate(  isbucket(max_mvseg,  nczdiv+2, ncydiv+2, ncxdiv+2))
    allocate(                 ncseg(nczdiv+2, ncydiv+2, ncxdiv+2))
    allocate(                ncatom(nczdiv+2, ncydiv+2, ncxdiv+2))
    allocate(buffp   (6,max_cellcbd*max_mvatom))
    allocate(buffm   (6,max_cellcbd*max_mvatom))
    allocate(ibuffp  (  max_cellcbd*max_mvatom))
    allocate(ibuffm  (  max_cellcbd*max_mvatom))
    allocate(isbufp  (2*max_cellcbd + 1 + max_cellcbd*max_mvseg))
    allocate(isbufm  (2*max_cellcbd + 1 + max_cellcbd*max_mvseg))
    allocate(rbuff_p (6,max_cellcbd*max_mvatom))
    allocate(rbuff_m (6,max_cellcbd*max_mvatom))
    allocate(irbuff_p(  max_cellcbd*max_mvatom))
    allocate(irbuff_m(  max_cellcbd*max_mvatom))
    allocate(irsbuf_p(2*max_cellcbd + 1 + max_cellcbd*max_mvseg))
    allocate(irsbuf_m(2*max_cellcbd + 1 + max_cellcbd*max_mvseg))
    allocate( ncatmw(32, nczdiv+2, ncydiv+2, ncxdiv+2) )
!
    abucket(:,:,:,:,:) = 0.0d0
    iabucket(:,:,:,:)  = 0
    isbucket(:,:,:,:)  = 0
    ncseg(:,:,:)       = 0
    ncatom(:,:,:)      = 0
    buffp(:,:)         = 0.0d0
    buffm(:,:)         = 0.0d0
    ibuffp(:)          = 0
    ibuffm(:)          = 0
    isbufp(:)          = 0
    isbufm(:)          = 0
    rbuff_p(:,:)       = 0.0d0
    rbuff_m(:,:)       = 0.0d0
    irbuff_p(:)        = 0
    irbuff_m(:)        = 0
    irsbuf_p(:)        = 0
    irsbuf_m(:)        = 0
    ncatmw(:,:,:,:)    = 0

  end subroutine init_comm_bound
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to communicate atoms sticked out from subcells
!! \author Yoshimichi Andoh, Shin-ichi Ichikawa
!<
  subroutine comm_bound()
!----------------------------------------------------------------------
    use trajectory_org
    use trajectory_mpi
    use subcell, only : wseg_cx, wseg_cy, wseg_cz
    use domain, only : ncellx, ncelly, ncellz
    use mpi_3d_grid
    use param
    use unit_cell
    use cell_shape
    use md_condition !! fapp
    use dr_cntl, only : nbd, nbd2
#include "timing.h90"
    implicit none
    INCLUDE 'mpif.h'
    integer nbase, nbase2
    real(8) rdcellz, rdcelly, rdcellx
    integer iczg0pr, icyg0pr, icxg0pr
    real(8) z0, y0, x0
    integer ipz_dest, ipy_dest, ipx_dest
    integer ipz_src, ipy_src, ipx_src
    integer nsa
    integer lcse(2,3,9)
    integer iczg, icyg, icxg
    integer icz, icy, icx
    integer icz0, icz1
    integer icy0, icy1
    integer icx0, icx1
    integer i,i0,k0
    integer ncc, ncs, ncsr, ics
    integer nca, ncar, ica, ncarp, ncarm
    integer nccp, ncsp, ncap
    integer nccm, ncsm, ncam
    integer ncc_p, ncs_p, nca_p
    integer ncc_m, ncs_m, nca_m
    integer ldcell
    integer loc_csbound
    !
    real(8) :: xfrac,yfrac,zfrac
    real(8) :: box(3,3),av(3),bv(3),cv(3)
! ...  effective only when OpenMP compile.
!$    integer ncthread(2, 2, 2), itz, ity, itx, iam, nth
!$    integer omp_get_thread_num
!$    external omp_get_thread_num
!$    integer omp_get_num_threads
!$    external omp_get_num_threads
!
    integer istatus(mpi_status_size), ierr

    !
    rdcellx=dble(ncellx)/cellx
    rdcelly=dble(ncelly)/celly
    rdcellz=dble(ncellz)/cellz
    !
    call cell_convert1(cellx,celly,cellz,alpha,beta,gamma,box)
    av(:)=box(:,1)
    bv(:)=box(:,2)
    cv(:)=box(:,3)
    call cell_convert3(cellx,celly,cellz,alpha,beta,gamma,cellxh,cellyh, &
         &      cellzh,cellvol,sinbeta,cosbeta, &
         &  singamma,cosgamma,b_factor,g_factor)
    !
    ! Clear i2m for all local atoms, which is essential
    ! to detect bugs in communications
    !
!$omp parallel do default(shared) &
!$omp& private(ncc,icz,icy,icx,i0,i)
!defult: MTD
    do ncc=1,(lxdiv+nbd2)*(lydiv+nbd2)*(lzdiv+nbd2)
       icz=mod(ncc-1,lzdiv+nbd2)                   ! Ichi
       icy=mod(ncc-1,(lzdiv+nbd2)*(lydiv+nbd2))    ! Ichi
       icy=icy/(lzdiv+nbd2)                        ! Ichi
       icx=(ncc-1)/((lzdiv+nbd2)*(lydiv+nbd2))   ! Ichi
       !       if(icx.ge.3.and.icx.le.ncxdiv+2 .and. &
       !    &     icy.ge.3.and.icy.le.ncydiv+2 .and. &
       !    &     icz.ge.3.and.icz.le.nczdiv+2) cycle
       do i0=tag(icz,icy,icx),tag(icz,icy,icx) &
            &        +na_per_cell(icz,icy,icx)-1
          i=m2i(i0)
!   if(i.lt.0) stop 'ERROR:: m2i' ! debug YA
          i2m(i)=-1
       end do ! ica
    end do ! ncc

    ! ---- 3D rank order rule. ----
    !     ipx=mod(myrank,npx)
    !     ipy=mod((myrank-ipx)/npx,npy)
    !     ipz=mod((myrank-ipx-ipy*npx)/(npx*npy),npz)
    !
    !-----  bucket preparation. -----
    !
    ! assumption ; system boundary condition is NOT applied to wseg_cx,
    !              wseg_cy, wseg_cz.
    !
    !
    ncseg = 0
    ncatom = 0
    ncatmw = 0

    ! cell index of the end of previous rank.
    iczg0pr = ncellz * ipz / npz
    icyg0pr = ncelly * ipy / npy
    icxg0pr = ncellx * ipx / npx
    !
!$omp parallel default(none) &
!$omp& private(x0,y0,z0) &
!$omp& private(icx,icy,icz) &
!$omp& private(icxg,icyg,iczg) &
!$omp& private(itx,ity,itz) &
!$omp& private(ncs,nsa,nca) &
!$omp& private(ics,ica) &
!$omp& private(ncthread) &
!$omp& private(iam) &
!$omp& private(nth) &
!$omp& private(xfrac,yfrac,zfrac) &
!$omp& shared(ncellz,ncelly,ncellx) &
!$omp& shared(nselfseg) &
!$omp& shared(wseg_cx,wseg_cy,wseg_cz) &
!$omp& shared(cellxh,cellyh,cellzh) &
!$omp& shared(cellx,celly,cellz) &
!$omp& shared(icxg0pr,icyg0pr,iczg0pr) &
!$omp& shared(cosbeta,singamma,cosgamma) &
!$omp& shared(g_factor,b_factor) &
!$omp& shared(rdcellx,rdcelly,rdcellz) &
!$omp& shared(ncxdiv,ncydiv,nczdiv) &
!$omp& shared(isbucket) &
!$omp& shared(iabucket) &
!$omp& shared(abucket) &
!$omp& shared(lsegtop,lseg_natoms) &
!$omp& shared(ncseg) &
!$omp& shared(wkxyz,wkv,m2i) &
!$omp& shared(ncatmw) &
!$omp& shared(av,bv,cv) &
!$omp& shared(myrank)
!$    iam = omp_get_thread_num()
!$    nth = omp_get_num_threads()
!$    ncthread = 0
!$    if(nth == 2) then
!$       ncthread(1,1,1) = 0
!$       ncthread(2,1,1) = 0
!$       ncthread(1,2,1) = 0
!$       ncthread(2,2,1) = 0
!$       ncthread(1,1,2) = 1
!$       ncthread(2,1,2) = 1
!$       ncthread(1,2,2) = 1
!$       ncthread(2,2,2) = 1
!$    elseif(nth == 4) then
! 1 to the 4th fractional subspace boundary.
! thread id for each subspace/domain. for FX-1.
!$       ncthread(1,1,1) = 0
!$       ncthread(2,1,1) = 0
!$       ncthread(1,2,1) = 1
!$       ncthread(2,2,1) = 1
!$       ncthread(1,1,2) = 2
!$       ncthread(2,1,2) = 2
!$       ncthread(1,2,2) = 3
!$       ncthread(2,2,2) = 3
!$    elseif(nth >= 8) then
! 1 to the 8th fractional subspace boundary.
! thread id for each subspace/domain. for K computer.
!$       ncthread(1,1,1) = 0
!$       ncthread(2,1,1) = 1
!$       ncthread(1,2,1) = 2
!$       ncthread(2,2,1) = 3
!$       ncthread(1,1,2) = 4
!$       ncthread(2,1,2) = 5
!$       ncthread(1,2,2) = 6
!$       ncthread(2,2,2) = 7
!$    endif
!
!... segment loop.
    do ics = 1, nselfseg
       zfrac=wseg_cz(ics)/g_factor
       yfrac=(wseg_cy(ics)-zfrac*b_factor)/singamma
       xfrac=wseg_cx(ics)-yfrac*cosgamma-zfrac*cosbeta
       x0 = xfrac+cellxh
       y0 = yfrac+cellyh
       z0 = zfrac+cellzh
       !        x0 = wseg_cx(ics) + cellxh
       !        y0 = wseg_cy(ics) + cellyh
       !        z0 = wseg_cz(ics) + cellzh
       !
       ! one extra cell outside process and beginning address 1 for cell index
       ! is to be considered.
       ! segment address is negative in the negative direction system boundary.
       ! since intrinsic function INT truncates to the 0-direction nearest
       ! whole number, addition of 1.0 before truncation is necessary.
       ! to avoid the addition intrinsic function FLOOR which truncate to the
       ! minus direction largest whole number is preferable.
       !
       ! global address.
       icxg = floor(x0*rdcellx) + 1
       icyg = floor(y0*rdcelly) + 1
       iczg = floor(z0*rdcellz) + 1
       ! global address to local address. one cell for boundary.
       icx = icxg - icxg0pr + 1
       icy = icyg - icyg0pr + 1
       icz = iczg - iczg0pr + 1
! sort cell for thread.
! ...  effective only when OpenMP compile.
!$       itx = 2
!$       if(icx <= ncxdiv/2+1) itx = 1
!$       ity = 2
!$       if(icy <= ncydiv/2+1) ity = 1
!$       itz = 2
!$       if(icz <= nczdiv/2+1) itz = 1

! ...  effective only when OpenMP compile.
!$       if(ncthread(itz,ity,itx) /= iam) cycle
       !
       ! periodic bondary condition to coordinate wkxyz.
       !
       if(icxg  <= 0) then
          do ica = 1, lseg_natoms(ics)
             nsa = lsegtop(ics) + ica - 1
             !              wkxyz(1,nsa) = wkxyz(1,nsa) + cellx
             wkxyz(:,nsa) = wkxyz(:,nsa) + av(:)
          end do
       endif
       if(icxg  > ncellx) then
          do ica = 1, lseg_natoms(ics)
             nsa = lsegtop(ics) + ica - 1
             !              wkxyz(1,nsa) = wkxyz(1,nsa) - cellx
             wkxyz(:,nsa) = wkxyz(:,nsa) - av(:)
          end do
       endif
       if(icyg  <= 0) then
          do ica = 1, lseg_natoms(ics)
             nsa = lsegtop(ics) + ica - 1
             !              wkxyz(2,nsa) = wkxyz(2,nsa) + celly
             wkxyz(:,nsa) = wkxyz(:,nsa) + bv(:)
          end do
       endif
       if(icyg  > ncelly) then
          do ica = 1, lseg_natoms(ics)
             nsa = lsegtop(ics) + ica - 1
             !              wkxyz(2,nsa) = wkxyz(2,nsa) - celly
             wkxyz(:,nsa) = wkxyz(:,nsa) - bv(:)
          end do
       endif
       if(iczg  <= 0) then
          do ica = 1, lseg_natoms(ics)
             nsa = lsegtop(ics) + ica - 1
             !              wkxyz(3,nsa) = wkxyz(3,nsa) + cellz
             wkxyz(:,nsa) = wkxyz(:,nsa) + cv(:)
          end do
       endif
       if(iczg  > ncellz) then
          do ica = 1, lseg_natoms(ics)
             nsa = lsegtop(ics) + ica - 1
             !              wkxyz(3,nsa) = wkxyz(3,nsa) - cellz
             wkxyz(:,nsa) = wkxyz(:,nsa) - cv(:)
          end do
       endif
       !
       ! segment data.
       ncs = ncatmw(2,icz,icy,icx) + 1
       isbucket(ncs,icz,icy,icx) = lseg_natoms(ics)
       ncatmw(2,icz,icy,icx) = ncs
       ! atom data.
       !     In the bounding process, all of the atom belonging to the same
       !     segment must be treated in the same way.
       nca = ncatmw(1,icz,icy,icx)
       do ica = 1, lseg_natoms(ics)
          nsa = lsegtop(ics) + ica - 1
          nca = nca + 1
          abucket(1,nca,icz,icy,icx) = wkxyz(1,nsa)
          abucket(2,nca,icz,icy,icx) = wkxyz(2,nsa)
          abucket(3,nca,icz,icy,icx) = wkxyz(3,nsa)
          abucket(4,nca,icz,icy,icx) = wkv(1,nsa)
          abucket(5,nca,icz,icy,icx) = wkv(2,nsa)
          abucket(6,nca,icz,icy,icx) = wkv(3,nsa)
          iabucket(nca,icz,icy,icx) = m2i(nsa)
!   if(m2i(nsa).lt.0) stop 'ERROR:: m2i(nsa)' !debug YA
       end do
       ncatmw(1,icz,icy,icx) = nca
    end do
!$OMP ENDPARALLEL

!$omp parallel default(shared) &
!$omp private(ncc,icx,icy,icz)
!$omp do
    do ncc=1,(ncxdiv+2)*(ncydiv+2)*(nczdiv+2)
       icz=mod(ncc-1,nczdiv+2)            +1
       icy=mod(ncc-1,(nczdiv+2)*(ncydiv+2))
       icy=icy/(nczdiv+2)                 +1
       icx=(ncc-1)/((nczdiv+2)*(ncydiv+2))+1
             ncatom(icz,icy,icx) = ncatmw(1,icz,icy,icx)
             ncseg(icz,icy,icx) = ncatmw(2,icz,icy,icx)
    end do ! ncc
!$omp end do
!$omp end parallel
!
!      if(mdstep >= 10)  TIME_STOP(TM_COM_BD_pre)
!
    !
    !-----  boundary communication code starts here. ------
    !
    ! ( bound +Z )
    ! nearest neighbor.
    icz = nczdiv + 2
    icy0 = 2
    icy1 = ncydiv + 1
    icx0 = 2
    icx1 = ncxdiv + 1
    ldcell = 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy0
    lcse(2,2,ldcell) = icy1
    lcse(1,3,ldcell) = icx0
    lcse(2,3,ldcell) = icx1
    ! 2d-diagonal adjacent(1)
    icz = nczdiv + 2
    icy = 1
    icx0 = 2
    icx1 = ncxdiv + 1
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx0
    lcse(2,3,ldcell) = icx1
    ! 2d-diagonal adjacent(2)
    icz = nczdiv + 2
    icy = ncydiv + 2
    icx0 = 2
    icx1 = ncxdiv + 1
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx0
    lcse(2,3,ldcell) = icx1
    ! 2d-diagonal adjacent(3)
    icz = nczdiv + 2
    icy0 = 2
    icy1 = ncydiv + 1
    icx = 1
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy0
    lcse(2,2,ldcell) = icy1
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    ! 2d-diagonal adjacent(4)
    icz = nczdiv + 2
    icy0 = 2
    icy1 = ncydiv + 1
    icx = ncxdiv + 2
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy0
    lcse(2,2,ldcell) = icy1
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    ! 3d-diagonal adjacent(5)
    icz = nczdiv + 2
    icy = 1
    icx = 1
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    ! 3d-diagonal adjacent(6)
    icz = nczdiv + 2
    icy = ncydiv + 2
    icx = 1
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    ! 3d-diagonal adjacent(7)
    icz = nczdiv + 2
    icy = 1
    icx = ncxdiv + 2
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    ! 3d-diagonal adjacent(8)
    icz = nczdiv + 2
    icy = ncydiv + 2
    icx = ncxdiv + 2
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    !
    loc_csbound = 2*(ncydiv*ncxdiv + 2*ncxdiv + 2*ncydiv + 4) + 1
    ncc = 0
    ncs = loc_csbound
    nca = 0
    !
    call add_buffer(abucket,iabucket,isbucket,buffp,ibuffp,isbufp, &
         &                max_mvatom,max_mvseg,max_cellcbd, &
         &                nczdiv,ncydiv,ncxdiv, &
         &                lcse,ldcell,ncatom,ncseg,ncc,ncs,nca)
    !
    ipz_dest = mod(ipz+1-1/npz+npz,npz)*npx*npy + ipy*npx + ipx
    ipz_src  = mod(ipz-1+1/npz+npz,npz)*npx*npy + ipy*npx + ipx
    !
    call mpi_sendrecv(ncs + 1, 1, MPI_INTEGER, &
         &             ipz_dest, myrank, &
         &             ncsr, 1, MPI_INTEGER, &
         &             ipz_src, ipz_src, &
         &             mpi_comm_world, istatus, ierr )

    isbufp(ncs+1) = nca
    call mpi_sendrecv(isbufp, ncs + 1, MPI_INTEGER, &
         &             ipz_dest, myrank, &
         &             irsbuf_p, ncsr, MPI_INTEGER, &
         &             ipz_src, ipz_src, &
         &             mpi_comm_world, istatus, ierr )

    ncar = irsbuf_p(ncsr)
MSG_SEND=6*nca
MSG_RECV=6*ncar
    call mpi_sendrecv(buffp, MSG_SEND, MPI_DOUBLE_PRECISION, &
         &             ipz_dest, myrank, &
         &             rbuff_p, MSG_RECV, MPI_DOUBLE_PRECISION, &
         &             ipz_src, ipz_src, &
         &             mpi_comm_world, istatus, ierr )

    call mpi_sendrecv(ibuffp, nca, MPI_INTEGER, &
         &             ipz_dest, myrank, &
         &             irbuff_p, ncar, MPI_INTEGER, &
         &             ipz_src, ipz_src, &
         &             mpi_comm_world, istatus, ierr )

    !
    ! ( bound -Z )
    !
    ! ( bound -Z )
    ! nearest neighbor.
    icz = 1
    icy0 = 2
    icy1 = ncydiv + 1
    icx0 = 2
    icx1 = ncxdiv + 1
    ldcell = 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy0
    lcse(2,2,ldcell) = icy1
    lcse(1,3,ldcell) = icx0
    lcse(2,3,ldcell) = icx1
    ! 2d-diagonal adjacent(13)
    icz = 1
    icy = 1
    icx0 = 2
    icx1 = ncxdiv + 1
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx0
    lcse(2,3,ldcell) = icx1
    ! 2d-diagonal adjacent(14)
    icz = 1
    icy = ncydiv + 2
    icx0 = 2
    icx1 = ncxdiv + 1
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx0
    lcse(2,3,ldcell) = icx1
    ! 2d-diagonal adjacent(15)
    icz = 1
    icy0 = 2
    icy1 = ncydiv + 1
    icx = 1
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy0
    lcse(2,2,ldcell) = icy1
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    ! 2d-diagonal adjacent(16)
    icz = 1
    icy0 = 2
    icy1 = ncydiv + 1
    icx = ncxdiv + 2
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy0
    lcse(2,2,ldcell) = icy1
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    ! 3d-diagonal adjacent(17)
    icz = 1
    icy = 1
    icx = 1
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    ! 3d-diagonal adjacent(18)
    icz = 1
    icy = ncydiv + 2
    icx = 1
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    ! 3d-diagonal adjacent(19)
    icz = 1
    icy = 1
    icx = ncxdiv + 2
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    ! 3d-diagonal adjacent(20)
    icz = 1
    icy = ncydiv + 2
    icx = ncxdiv + 2
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    !
    loc_csbound = 2*(ncydiv*ncxdiv + 2*ncxdiv + 2*ncydiv + 4) + 1
    ncc = 0
    ncs = loc_csbound
    nca = 0
    !
    call add_buffer(abucket,iabucket,isbucket,buffm,ibuffm,isbufm, &
         &                max_mvatom,max_mvseg,max_cellcbd, &
         &                nczdiv,ncydiv,ncxdiv, &
         &                lcse,ldcell,ncatom,ncseg,ncc,ncs,nca)
    !
    ipz_dest = mod(ipz-1+1/npz+npz,npz)*npx*npy + ipy*npx + ipx
    ipz_src  = mod(ipz+1-1/npz+npz,npz)*npx*npy + ipy*npx + ipx
    !
    call mpi_sendrecv(ncs + 1, 1, MPI_INTEGER, &
         &             ipz_dest, myrank, &
         &             ncsr, 1, MPI_INTEGER, &
         &             ipz_src, ipz_src, &
         &             mpi_comm_world, istatus, ierr )

    isbufm(ncs+1) = nca
    call mpi_sendrecv(isbufm, ncs + 1, MPI_INTEGER, &
         &             ipz_dest, myrank, &
         &             irsbuf_m, ncsr, MPI_INTEGER, &
         &             ipz_src, ipz_src, &
         &             mpi_comm_world, istatus, ierr )

    ncar = irsbuf_m(ncsr)
MSG_SEND=6*nca
MSG_RECV=6*ncar
    call mpi_sendrecv(buffm, MSG_SEND, MPI_DOUBLE_PRECISION, &
         &             ipz_dest, myrank, &
         &             rbuff_m, MSG_RECV, MPI_DOUBLE_PRECISION, &
         &             ipz_src, ipz_src, &
         &             mpi_comm_world, istatus, ierr )

    call mpi_sendrecv(ibuffm, nca, MPI_INTEGER, &
         &             ipz_dest, myrank, &
         &             irbuff_m, ncar, MPI_INTEGER, &
         &             ipz_src, ipz_src, &
         &             mpi_comm_world, istatus, ierr )

    !
    ! merge source nearest neighbors(+Z) on receive buffer to local bucket.
    icx0 = 2
    icx1 = ncxdiv + 1
    icy0 = 2
    icy1 = ncydiv + 1
    icz = 2
    ldcell = 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy0
    lcse(2,2,ldcell) = icy1
    lcse(1,3,ldcell) = icx0
    lcse(2,3,ldcell) = icx1
    !
    !  Merge part of 2d-diagonal adjacent on z-direction receive buffer
    !  to cell bucket corresponding to y-directional nearest neighbor
    !  communication.
    ! source rank +Z 2d-diagonal adjacent(1)
    icx0 = 2
    icx1 = ncxdiv + 1
    icy = 1
    icz = 2
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx0
    lcse(2,3,ldcell) = icx1
    ! source rank +Z 2d-diagonal adjacent(2)
    icx0 = 2
    icx1 = ncxdiv + 1
    icy = ncydiv + 2
    icz = 2
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx0
    lcse(2,3,ldcell) = icx1
    ! source rank +Z 2d-diagonal adjacent(3)
    icx = 1
    icy0 = 2
    icy1 = ncydiv + 1
    icz = 2
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy0
    lcse(2,2,ldcell) = icy1
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    ! source rank +Z 2d-diagonal adjacent(4)
    icx = ncxdiv+ 2
    icy0 = 2
    icy1 = ncydiv + 1
    icz = 2
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy0
    lcse(2,2,ldcell) = icy1
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    !
    loc_csbound = 2*(ncydiv*ncxdiv + 2*ncxdiv + 2*ncydiv + 4) + 1
    ncc_p = 0
    ncs_p = loc_csbound
    nca_p = 0
    !
    call add_bucket(abucket,iabucket,isbucket, &
         &                rbuff_p,irbuff_p,irsbuf_p, &
         &                max_mvatom,max_mvseg,max_cellcbd, &
         &                nczdiv,ncydiv,ncxdiv, &
         &                lcse,ldcell,ncatom,ncseg, &
         &                ncc_p,ncs_p,nca_p)
    !
    ! merge source rank nearest neighbors(-Z) on receive buffer to local bucket.
    icx0 = 2
    icx1 = ncxdiv + 1
    icy0 = 2
    icy1 = ncydiv + 1
    icz = nczdiv + 1
    ldcell = 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy0
    lcse(2,2,ldcell) = icy1
    lcse(1,3,ldcell) = icx0
    lcse(2,3,ldcell) = icx1
    !
    !  Merge part of 2d-diagonal adjacent on z-direction receive buffer
    !  to cell bucket corresponding to y-directional nearest neighbor
    !  communication.
    ! source rank -Z 2d-diagonal adjacent(13)
    icx0 = 2
    icx1 = ncxdiv + 1
    icy = 1
    icz = nczdiv + 1
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx0
    lcse(2,3,ldcell) = icx1
    ! source rank -Z 2d-diagonal adjacent(14)
    icx0 = 2
    icx1 = ncxdiv + 1
    icy = ncydiv + 2
    icz = nczdiv + 1
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx0
    lcse(2,3,ldcell) = icx1
    ! source rank -Z 2d-diagonal adjacent(15)
    icx = 1
    icy0 = 2
    icy1 = ncydiv + 1
    icz = nczdiv + 1
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy0
    lcse(2,2,ldcell) = icy1
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    ! source rank -Z 2d-diagonal adjacent(16)
    icx = ncxdiv + 2
    icy0 = 2
    icy1 = ncydiv + 1
    icz = nczdiv + 1
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy0
    lcse(2,2,ldcell) = icy1
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    !
    loc_csbound = 2*(ncydiv*ncxdiv + 2*ncxdiv + 2*ncydiv + 4) + 1
    ncc_m = 0
    ncs_m = loc_csbound
    nca_m = 0
    !
    call add_bucket(abucket,iabucket,isbucket, &
         &                rbuff_m,irbuff_m,irsbuf_m, &
         &                max_mvatom,max_mvseg,max_cellcbd, &
         &                nczdiv,ncydiv,ncxdiv, &
         &                lcse,ldcell,ncatom,ncseg, &
         &                ncc_m,ncs_m,nca_m)

    !
    ! ( bound +Y & bound -Y )
    ! ( +Y buffer )
    ! nearest neighbor.
    icz0 = 2
    icz1 = nczdiv + 1
    icy = ncydiv + 2
    icx0 = 2
    icx1 = ncxdiv + 1
    ldcell = 1
    lcse(1,1,ldcell) = icz0
    lcse(2,1,ldcell) = icz1
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx0
    lcse(2,3,ldcell) = icx1
    ! 2d-diagonal adjacent(10)
    icz0 = 2
    icz1 = nczdiv + 1
    icy = ncydiv + 2
    icx = 1
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz0
    lcse(2,1,ldcell) = icz1
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    ! 2d-diagonal adjacent(12)
    icz0 = 2
    icz1 = nczdiv + 1
    icy = ncydiv + 2
    icx = ncxdiv + 2
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz0
    lcse(2,1,ldcell) = icz1
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    !
    loc_csbound = 2*(nczdiv*ncxdiv + 2*nczdiv + 4) + 1
    nccp = 0
    ncsp = loc_csbound
    ncap = 0
    !
    call add_buffer(abucket,iabucket,isbucket,buffp,ibuffp,isbufp, &
         &                max_mvatom,max_mvseg,max_cellcbd, &
         &                nczdiv,ncydiv,ncxdiv, &
         &                lcse,ldcell,ncatom,ncseg, &
         &                nccp,ncsp,ncap)
    !
    ! ( -Y buffer )
    ! nearest neighbor.
    icz0 = 2
    icz1 = nczdiv + 1
    icy = 1
    icx0 = 2
    icx1 = ncxdiv + 1
    ldcell = 1
    lcse(1,1,ldcell) = icz0
    lcse(2,1,ldcell) = icz1
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx0
    lcse(2,3,ldcell) = icx1
    ! 2d-diagonal adjacent(9)
    icz0 = 2
    icz1 = nczdiv + 1
    icy = 1
    icx = 1
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz0
    lcse(2,1,ldcell) = icz1
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    ! 2d-diagonal adjacent(11)
    icz0 = 2
    icz1 = nczdiv + 1
    icy = 1
    icx = ncxdiv + 2
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz0
    lcse(2,1,ldcell) = icz1
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    !
    loc_csbound = 2*(nczdiv*ncxdiv +  2*nczdiv + 4) + 1
    nccm = 0
    ncsm = loc_csbound
    ncam = 0
    !
    call add_buffer(abucket,iabucket,isbucket,buffm,ibuffm,isbufm, &
         &                max_mvatom,max_mvseg,max_cellcbd, &
         &                nczdiv,ncydiv,ncxdiv, &
         &                lcse,ldcell,ncatom,ncseg,nccm,ncsm,ncam)
    !
    ! ( -Y buffer ) +Z source rank 3d-diagonal adjacent(5)
    !
    ldcell = 1
    call add_buffb(rbuff_p,irbuff_p,irsbuf_p,buffm,ibuffm,isbufm, &
         &               max_mvatom,max_mvseg,max_cellcbd, &
         &               ncc_p,ncs_p,nca_p,nccm,ncsm,ncam,ldcell)
    !
    ! ( +Y buffer ) +Z source rank 3d-diagonal adjacent(6)
    !
    ldcell = 1
    call add_buffb(rbuff_p,irbuff_p,irsbuf_p,buffp,ibuffp,isbufp, &
         &               max_mvatom,max_mvseg,max_cellcbd, &
         &               ncc_p,ncs_p,nca_p,nccp,ncsp,ncap,ldcell)
    !
    ! ( -Y buffer ) +Z source rank 3d-diagonal adjacent(7)
    !
    ldcell = 1
    call add_buffb(rbuff_p,irbuff_p,irsbuf_p,buffm,ibuffm,isbufm, &
         &               max_mvatom,max_mvseg,max_cellcbd, &
         &               ncc_p,ncs_p,nca_p,nccm,ncsm,ncam,ldcell)
    !
    ! ( +Y buffer ) +Z source rank 3d-diagonal adjacent(8)
    !
    ldcell = 1
    call add_buffb(rbuff_p,irbuff_p,irsbuf_p,buffp,ibuffp,isbufp, &
         &               max_mvatom,max_mvseg,max_cellcbd, &
         &               ncc_p,ncs_p,nca_p,nccp,ncsp,ncap,ldcell)
    !
    ! ( -Y buffer ) -Z source rank 3d-diagonal adjacent(17)
    !
    ldcell = 1
    call add_buffb(rbuff_m,irbuff_m,irsbuf_m,buffm,ibuffm,isbufm, &
         &               max_mvatom,max_mvseg,max_cellcbd, &
         &               ncc_m,ncs_m,nca_m,nccm,ncsm,ncam,ldcell)
    !
    ! ( +Y buffer ) -Z source rank 3d-diagonal adjacent(18)
    !
    ldcell = 1
    call add_buffb(rbuff_m,irbuff_m,irsbuf_m,buffp,ibuffp,isbufp, &
         &               max_mvatom,max_mvseg,max_cellcbd, &
         &               ncc_m,ncs_m,nca_m,nccp,ncsp,ncap,ldcell)
    !
    ! ( -Y buffer ) -Z source rank 3d-diagonal adjacent(19)
    !
    ldcell = 1
    call add_buffb(rbuff_m,irbuff_m,irsbuf_m,buffm,ibuffm,isbufm, &
         &               max_mvatom,max_mvseg,max_cellcbd, &
         &               ncc_m,ncs_m,nca_m,nccm,ncsm,ncam,ldcell)
    !
    ! ( +Y buffer ) -Z source rank 3d-diagonal adjacent(20)
    !
    ldcell = 1
    call add_buffb(rbuff_m,irbuff_m,irsbuf_m,buffp,ibuffp,isbufp, &
         &               max_mvatom,max_mvseg,max_cellcbd, &
         &               ncc_m,ncs_m,nca_m,nccp,ncsp,ncap,ldcell)
    !
    !
    ! +Y comm.
    ipy_dest  = ipz*npx*npy + mod(ipy+1-1/npy+npy,npy)*npx + ipx
    ipy_src   = ipz*npx*npy + mod(ipy-1+1/npy+npy,npy)*npx + ipx

    !
    call mpi_sendrecv(ncsp + 1, 1, MPI_INTEGER, &
         &             ipy_dest, myrank, &
         &             ncsr, 1, MPI_INTEGER, &
         &             ipy_src, ipy_src, &
         &             mpi_comm_world, istatus, ierr )

    isbufp(ncsp+1) = ncap
    call mpi_sendrecv(isbufp, ncsp + 1, MPI_INTEGER, &
         &             ipy_dest, myrank, &
         &             irsbuf_p, ncsr, MPI_INTEGER, &
         &             ipy_src, ipy_src, &
         &             mpi_comm_world, istatus, ierr )

    ncarp = irsbuf_p(ncsr)
MSG_SEND=6*ncap
MSG_RECV=6*ncarp
    call mpi_sendrecv(buffp, MSG_SEND, MPI_DOUBLE_PRECISION, &
         &             ipy_dest, myrank, &
         &             rbuff_p, MSG_RECV, MPI_DOUBLE_PRECISION, &
         &             ipy_src, ipy_src, &
         &             mpi_comm_world, istatus, ierr )

    call mpi_sendrecv(ibuffp, ncap, MPI_INTEGER, &
         &             ipy_dest, myrank, &
         &             irbuff_p, ncarp, MPI_INTEGER, &
         &             ipy_src, ipy_src, &
         &             mpi_comm_world, istatus, ierr )

    !
    ! -Y comm.
    !
    ipy_dest  = ipz*npx*npy + mod(ipy-1+1/npy+npy,npy)*npx + ipx
    ipy_src   = ipz*npx*npy + mod(ipy+1-1/npy+npy,npy)*npx + ipx
    !
    call mpi_sendrecv(ncsm + 1, 1, MPI_INTEGER, &
         &             ipy_dest, myrank, &
         &             ncsr, 1, MPI_INTEGER, &
         &             ipy_src, ipy_src, &
         &             mpi_comm_world, istatus, ierr )

    isbufm(ncsm+1) = ncam
    call mpi_sendrecv(isbufm, ncsm + 1, MPI_INTEGER, &
         &             ipy_dest, myrank, &
         &             irsbuf_m, ncsr, MPI_INTEGER, &
         &             ipy_src, ipy_src, &
         &             mpi_comm_world, istatus, ierr )

    ncarm = irsbuf_m(ncsr)
MSG_SEND=6*ncam
MSG_RECV=6*ncarm
    call mpi_sendrecv(buffm, MSG_SEND, MPI_DOUBLE_PRECISION, &
         &             ipy_dest, myrank, &
         &             rbuff_m, MSG_RECV, MPI_DOUBLE_PRECISION, &
         &             ipy_src, ipy_src, &
         &             mpi_comm_world, istatus, ierr )

    call mpi_sendrecv(ibuffm, ncam, MPI_INTEGER, &
         &             ipy_dest, myrank, &
         &             irbuff_m, ncarm, MPI_INTEGER, &
         &             ipy_src, ipy_src, &
         &             mpi_comm_world, istatus, ierr )

    !
    ! +Y receive buffer.
    ! merge source rank nearest neighbors(+Y) on receive buffer to local bucket.
    icx0 = 2
    icx1 = ncxdiv + 1
    icy = 2
    icz0 = 2
    icz1 = nczdiv + 1
    ldcell = 1
    lcse(1,1,ldcell) = icz0
    lcse(2,1,ldcell) = icz1
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx0
    lcse(2,3,ldcell) = icx1
    ! source rank +Y 2d-diagonal adjacent(10)
    icx = 1
    icy = 2
    icz0 = 2
    icz1 = nczdiv + 1
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz0
    lcse(2,1,ldcell) = icz1
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    ! source rank +Y 2d-diagonal adjacent(12)
    icx = ncxdiv+ 2
    icy = 2
    icz0 = 2
    icz1 = nczdiv + 1
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz0
    lcse(2,1,ldcell) = icz1
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    ! source rank +Z 3d-diagonal adjacent(6)
    icx = 1
    icy = 2
    icz = 2
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    ! source rank +Z 3d-diagonal adjacent(8)
    icx = ncxdiv + 2
    icy = 2
    icz = 2
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    ! source rank -Z 3d-diagonal adjacent(18)
    icx = 1
    icy = 2
    icz = nczdiv + 1
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    ! source rank -Z 3d-diagonal adjacent(20)
    icx = ncxdiv + 2
    icy = 2
    icz = nczdiv + 1
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    !
    loc_csbound = 2*(nczdiv*ncxdiv + 2*nczdiv + 4) + 1
    ncc_p = 0
    ncs_p = loc_csbound
    nca_p = 0
    !
    call add_bucket(abucket,iabucket,isbucket, &
         &                rbuff_p,irbuff_p,irsbuf_p, &
         &                max_mvatom,max_mvseg,max_cellcbd, &
         &                nczdiv,ncydiv,ncxdiv, &
         &                lcse,ldcell,ncatom,ncseg, &
         &                ncc_p,ncs_p,nca_p)
    !
    ! -Y receive buffer.
    ! merge source rank nearest neighbors(-Y) on receive buffer to local bucket.
    icx0 = 2
    icx1 = ncxdiv + 1
    icy = ncydiv + 1
    icz0 = 2
    icz1 = nczdiv + 1
    ldcell = 1
    lcse(1,1,ldcell) = icz0
    lcse(2,1,ldcell) = icz1
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx0
    lcse(2,3,ldcell) = icx1
    ! source rank -Y 2d-diagonal adjacent(9)
    icx = 1
    icy = ncydiv + 1
    icz0 = 2
    icz1 = nczdiv + 1
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz0
    lcse(2,1,ldcell) = icz1
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    ! source rank -Y 2d-diagonal adjacent(11)
    icx = ncxdiv+ 2
    icy = ncydiv + 1
    icz0 = 2
    icz1 = nczdiv + 1
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz0
    lcse(2,1,ldcell) = icz1
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    ! source rank +Z 3d-diagonal adjacent(5)
    icx = 1
    icy = ncydiv + 1
    icz = 2
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    ! source rank +Z 3d-diagonal adjacent(7)
    icx = ncxdiv + 2
    icy = ncydiv + 1
    icz = 2
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    ! source rank -Z 3d-diagonal adjacent(17)
    icx = 1
    icy = ncydiv + 1
    icz = nczdiv + 1
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    ! source rank -Z 3d-diagonal adjacent(19)
    icx = ncxdiv + 2
    icy = ncydiv + 1
    icz = nczdiv + 1
    ldcell = ldcell + 1
    lcse(1,1,ldcell) = icz
    lcse(2,1,ldcell) = icz
    lcse(1,2,ldcell) = icy
    lcse(2,2,ldcell) = icy
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    !
    loc_csbound = 2*(nczdiv*ncxdiv + 2*nczdiv + 4) + 1
    ncc_m = 0
    ncs_m = loc_csbound
    nca_m = 0
    !
    call add_bucket(abucket,iabucket,isbucket, &
         &                rbuff_m,irbuff_m,irsbuf_m, &
         &                max_mvatom,max_mvseg,max_cellcbd, &
         &                nczdiv,ncydiv,ncxdiv, &
         &                lcse,ldcell,ncatom,ncseg, &
         &                ncc_m,ncs_m,nca_m)

    !
    ! ( bound +X )
    !
    ! nearest neighbor.
    icz0 = 2
    icz1 = nczdiv + 1
    icy0 = 2
    icy1 = ncydiv + 1
    icx = ncxdiv + 2
    ldcell = 1
    lcse(1,1,ldcell) = icz0
    lcse(2,1,ldcell) = icz1
    lcse(1,2,ldcell) = icy0
    lcse(2,2,ldcell) = icy1
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    !
    loc_csbound = 2*nczdiv*ncydiv + 1
    ncc = 0
    ncs = loc_csbound
    nca = 0
    !
    call add_buffer(abucket,iabucket,isbucket,buffp,ibuffp,isbufp, &
         &                max_mvatom,max_mvseg,max_cellcbd, &
         &                nczdiv,ncydiv,ncxdiv, &
         &                lcse,ldcell,ncatom,ncseg,ncc,ncs,nca)
    !
    ipx_dest  = ipz*npx*npy + ipy*npx + mod(ipx+1-1/npx+npx,npx)
    ipx_src   = ipz*npx*npy + ipy*npx + mod(ipx-1+1/npx+npx,npx)
    !
    call mpi_sendrecv(ncs + 1, 1, MPI_INTEGER, &
         &             ipx_dest, myrank, &
         &             ncsr, 1, MPI_INTEGER, &
         &             ipx_src, ipx_src, &
         &             mpi_comm_world, istatus, ierr )

    isbufp(ncs+1) = nca
    call mpi_sendrecv(isbufp, ncs + 1, MPI_INTEGER, &
         &             ipx_dest, myrank, &
         &             irsbuf_p, ncsr, MPI_INTEGER, &
         &             ipx_src, ipx_src, &
         &             mpi_comm_world, istatus, ierr )

    ncar = irsbuf_p(ncsr)
MSG_SEND=6*nca
MSG_RECV=6*ncar
    call mpi_sendrecv(buffp, MSG_SEND, MPI_DOUBLE_PRECISION, &
         &             ipx_dest, myrank, &
         &             rbuff_p, MSG_RECV, MPI_DOUBLE_PRECISION, &
         &             ipx_src, ipx_src, &
         &             mpi_comm_world, istatus, ierr )

    call mpi_sendrecv(ibuffp, nca, MPI_INTEGER, &
         &             ipx_dest, myrank, &
         &             irbuff_p, ncar, MPI_INTEGER, &
         &             ipx_src, ipx_src, &
         &             mpi_comm_world, istatus, ierr )

    !
    ! merge nearest neighbors(+X) on receive buffer to local bucket.
    icx = 2
    icy0 = 2
    icy1 = ncydiv + 1
    icz0 = 2
    icz1 = nczdiv + 1
    ldcell = 1
    lcse(1,1,ldcell) = icz0
    lcse(2,1,ldcell) = icz1
    lcse(1,2,ldcell) = icy0
    lcse(2,2,ldcell) = icy1
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    !
    loc_csbound = 2*nczdiv*ncydiv + 1
    ncc_p = 0
    ncs_p = loc_csbound
    nca_p = 0
    !
    call add_bucket(abucket,iabucket,isbucket, &
         &                rbuff_p,irbuff_p,irsbuf_p, &
         &                max_mvatom,max_mvseg,max_cellcbd, &
         &                nczdiv,ncydiv,ncxdiv, &
         &                lcse,ldcell,ncatom,ncseg, &
         &                ncc_p,ncs_p,nca_p)

    !
    ! ( bound -X )
    !
    ! nearest neighbor.
    icz0 = 2
    icz1 = nczdiv + 1
    icy0 = 2
    icy1 = ncydiv + 1
    icx = 1
    ldcell = 1
    lcse(1,1,ldcell) = icz0
    lcse(2,1,ldcell) = icz1
    lcse(1,2,ldcell) = icy0
    lcse(2,2,ldcell) = icy1
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    !
    loc_csbound = 2*nczdiv*ncydiv + 1
    ncc = 0
    ncs = loc_csbound
    nca = 0
    !
    call add_buffer(abucket,iabucket,isbucket,buffm,ibuffm,isbufm, &
         &                max_mvatom,max_mvseg,max_cellcbd, &
         &                nczdiv,ncydiv,ncxdiv, &
         &                lcse,ldcell,ncatom,ncseg,ncc,ncs,nca)
    !
    ipx_dest  = ipz*npx*npy + ipy*npx + mod(ipx-1+1/npx+npx,npx)
    ipx_src   = ipz*npx*npy + ipy*npx + mod(ipx+1-1/npx+npx,npx)

    !
    call mpi_sendrecv(ncs + 1, 1, MPI_INTEGER, &
         &             ipx_dest, myrank, &
         &             ncsr, 1, MPI_INTEGER, &
         &             ipx_src, ipx_src, &
         &             mpi_comm_world, istatus, ierr )

    isbufm(ncs+1) = nca
    call mpi_sendrecv(isbufm, ncs + 1, MPI_INTEGER, &
         &             ipx_dest, myrank, &
         &             irsbuf_m, ncsr, MPI_INTEGER, &
         &             ipx_src, ipx_src, &
         &             mpi_comm_world, istatus, ierr )

    ncar = irsbuf_m(ncsr)
MSG_SEND=6*nca
MSG_RECV=6*ncar
    call mpi_sendrecv(buffm, MSG_SEND, MPI_DOUBLE_PRECISION, &
         &             ipx_dest, myrank, &
         &             rbuff_m, MSG_RECV, MPI_DOUBLE_PRECISION, &
         &             ipx_src, ipx_src, &
         &             mpi_comm_world, istatus, ierr )

    call mpi_sendrecv(ibuffm, nca, MPI_INTEGER, &
         &             ipx_dest, myrank, &
         &             irbuff_m, ncar, MPI_INTEGER, &
         &             ipx_src, ipx_src, &
         &             mpi_comm_world, istatus, ierr )

    !
    ! merge nearest neighbors(-X) on receive buffer to local bucket.
    icx = ncxdiv + 1
    icy0 = 2
    icy1 = ncydiv + 1
    icz0 = 2
    icz1 = nczdiv + 1
    ldcell = 1
    lcse(1,1,ldcell) = icz0
    lcse(2,1,ldcell) = icz1
    lcse(1,2,ldcell) = icy0
    lcse(2,2,ldcell) = icy1
    lcse(1,3,ldcell) = icx
    lcse(2,3,ldcell) = icx
    !
    loc_csbound = 2*nczdiv*ncydiv + 1
    ncc_m = 0
    ncs_m = loc_csbound
    nca_m = 0
    !
    call add_bucket(abucket,iabucket,isbucket, &
         &                rbuff_m,irbuff_m,irsbuf_m, &
         &                max_mvatom,max_mvseg,max_cellcbd, &
         &                nczdiv,ncydiv,ncxdiv, &
         &                lcse,ldcell,ncatom,ncseg, &
         &                ncc_m,ncs_m,nca_m)
    ! ------- create new cell meta-data and its entity ---------------
    !         create segment meta-data and its entity,
    !         ie. the number of atoms per segment.
    !
    ! copy atom data from bucket to new meta-data structure.
    ! also, setup meta-data "tag" and "na_per_cell",
    ! and also, setup segment meta-data "lsegtop" and "lseg_natoms".
    !
!default: MTD
    narea  = na1cell * (lzdiv+nbd2) * (lydiv+nbd2)
    naline = na1cell * (lzdiv+nbd2)
    !
!$omp parallel default(none) &
!$omp& private(ncc,icx,icy,icz,ica) &
!$omp& private(nca,nbase,nbase2) &
!$omp& shared(ncxdiv,ncydiv,nczdiv,na1cell) &
!$omp& shared(narea,naline) &
!$omp& shared(ncatom) &
!$omp& shared(wkxyz,wkv,m2i,abucket,iabucket)
!$omp do
    do ncc = 1, ncxdiv*ncydiv
!default: MTD
       icx = (ncc - 1)/ncydiv + 2
       icy = mod(ncc - 1, ncydiv) + 2
       nbase = (icx-2+nbd)*narea
       nbase2 = nbase + (icy-2+nbd)*naline
       nca = nbase2 + (2-2+nbd)*na1cell
       do icz = 2, nczdiv + 1
          do ica = 1, ncatom(icz,icy,icx)
             nca = nca + 1
             wkxyz(1,nca) = abucket(1,ica,icz,icy,icx)
             wkxyz(2,nca) = abucket(2,ica,icz,icy,icx)
             wkxyz(3,nca) = abucket(3,ica,icz,icy,icx)
             wkv(1,nca) = abucket(4,ica,icz,icy,icx)
             wkv(2,nca) = abucket(5,ica,icz,icy,icx)
             wkv(3,nca) = abucket(6,ica,icz,icy,icx)
             m2i(nca) = iabucket(ica,icz,icy,icx)
!   if(m2i(nca).lt.0) stop 'ERROR:: m2i(nca)'   !debug YA
          end do
       end do
       !        end do
    end do
!$OMP ENDPARALLEL

    !
    ! segment meta-data.
    ncs = 0

!default: MTD
    nbase = 0
    do icx = 2, ncxdiv + 1
       nbase = nbase + narea
!default: MTD
       nbase2 = nbase
       do icy = 2, ncydiv + 1
          nbase2 = nbase2 + naline
!default: MTD
          nca = nbase2 + 1*na1cell    ! Ichi
          do icz = 2, nczdiv + 1
             do ics = 1, ncseg(icz,icy,icx)
                ncs = ncs + 1
                lseg_natoms(ncs) = isbucket(ics,icz,icy,icx)
                lsegtop(ncs) = nca + 1
                nca = nca + lseg_natoms(ncs)
             end do
          end do
       end do
    end do
    ! number of segment per process.
    nselfseg = ncs

    ! cell meta-data.
!default: MTD
    nbase = 0
    do icx = 2, ncxdiv + 1
       nbase = nbase + narea
!default: MTD
       nbase2 = nbase
       do icy = 2, ncydiv + 1
          nbase2 = nbase2 + naline
!default: MTD
          nca = nbase2 + 1*na1cell   ! Ichi
          do icz = 2, nczdiv + 1
!default: MTD
             tag(icz-2+nbd, icy-2+nbd, icx-2+nbd) = nca + 1
             nca = nca + ncatom(icz,icy,icx)
!default: MTD
             na_per_cell(icz-2+nbd, icy-2+nbd, icx-2+nbd) = ncatom(icz,icy,icx)
          end do
       end do
    end do

    !
    !update of i2m
    !
       nselfatm=0
!$omp parallel default(shared) &
!$omp& private(i,i0,k0) &
!$omp& reduction(+:nselfatm)
!$omp do
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          i=m2i(i0)
          i2m(i)=i0
          nselfatm=nselfatm+1
       enddo
    enddo
!$omp end do
!$omp end parallel

  end subroutine comm_bound
!----------------------------------------------------------------------
  subroutine add_buffer(abucket,iabucket,isbucket,buff,ibuff,isbuf, &
       &                      max_mvatom,max_mvseg,max_cell, &
       &                      nczdiv,ncydiv,ncxdiv, &
       &                      lcse,ldcell,ncatom,ncseg,ncc,ncs,nca)
    !----------------------------------------------------------------------
    !
    implicit none
    integer max_mvatom
    integer max_mvseg
    integer max_cell        ! max number of cells on communication buffer.
    integer nczdiv
    integer ncydiv
    integer ncxdiv
    real(8) abucket(6,max_mvatom,nczdiv+2,ncydiv+2,ncxdiv+2)
    integer iabucket(max_mvatom,nczdiv+2,ncydiv+2,ncxdiv+2)
    integer isbucket(max_mvseg,nczdiv+2,ncydiv+2,ncxdiv+2)
    real(8) buff(6,max_cell*max_mvatom)
    integer ibuff(max_cell*max_mvatom)
    integer isbuf(2*max_cell + 1 + max_cell*max_mvseg)
    !
    integer ldcell, ldc
    integer lcse(2,3,ldcell)
    integer ncatom(nczdiv+2,ncydiv+2,ncxdiv+2)
    integer ncseg(nczdiv+2,ncydiv+2,ncxdiv+2)
    integer ncc
    integer ncs
    integer nca
    integer ics, ica
    integer icz, icy, icx
    integer icz0, icy0, icx0
    integer icz1, icy1, icx1
    !
    do ldc = 1, ldcell
       icz0 = lcse(1,1,ldc)
       icz1 = lcse(2,1,ldc)
       icy0 = lcse(1,2,ldc)
       icy1 = lcse(2,2,ldc)
       icx0 = lcse(1,3,ldc)
       icx1 = lcse(2,3,ldc)
       do icx = icx0, icx1
          do icy = icy0, icy1
             do icz = icz0, icz1
                ncc = ncc + 1
                isbuf(ncc) = ncatom(icz,icy,icx)
                ncc = ncc + 1
                isbuf(ncc) = ncseg(icz,icy,icx)
                do ics = 1, ncseg(icz,icy,icx)
                   ncs = ncs + 1
                   isbuf(ncs) = isbucket(ics,icz,icy,icx)
                end do
                do ica = 1, ncatom(icz,icy,icx)
                   nca = nca + 1
                   buff(1,nca) = abucket(1,ica,icz,icy,icx)
                   buff(2,nca) = abucket(2,ica,icz,icy,icx)
                   buff(3,nca) = abucket(3,ica,icz,icy,icx)
                   buff(4,nca) = abucket(4,ica,icz,icy,icx)
                   buff(5,nca) = abucket(5,ica,icz,icy,icx)
                   buff(6,nca) = abucket(6,ica,icz,icy,icx)
                   ibuff(nca)  = iabucket(ica,icz,icy,icx)
                end do
             end do
          end do
       end do
    end do
    !
  end subroutine add_buffer
!
!----------------------------------------------------------------------
  subroutine add_buffb(rbuff,irbuff,irsbuf,buff,ibuff,isbuf, &
     &                    max_mvatom,max_mvseg,max_cell, &
     &                    nccr,ncsr,ncar,ncc,ncs,nca,ldcell)
!----------------------------------------------------------------------
!
    implicit none
    integer max_mvatom
    integer max_mvseg
    integer max_cell        ! max number of cells on communication buffer.

    real(8) rbuff (6,max_cell*max_mvatom)
    integer irbuff(max_cell*max_mvatom)
    integer irsbuf(2*max_cell + 1 + max_cell*max_mvseg)
    real(8) buff  (6,max_cell*max_mvatom)
    integer ibuff (max_cell*max_mvatom)
    integer isbuf (2*max_cell + 1 + max_cell*max_mvseg)

    integer ncc, nccr
    integer ncs, ncsr
    integer nca, ncar
    integer ldcell, ldc
    integer ics, ica
    !
    do ldc = 1, ldcell
       ncc = ncc + 2
       nccr = nccr + 2
       isbuf(ncc - 1) = irsbuf(nccr - 1)
       isbuf(ncc)     = irsbuf(nccr)
       ! segment data.
       do ics = 1, irsbuf(nccr)
          ncs = ncs + 1
          ncsr = ncsr + 1
          isbuf(ncs) = irsbuf(ncsr)
       end do
       ! atom data.
       do ica = 1, irsbuf(nccr - 1)
          nca = nca + 1
          ncar = ncar + 1
          buff(1,nca) = rbuff(1,ncar)
          buff(2,nca) = rbuff(2,ncar)
          buff(3,nca) = rbuff(3,ncar)
          buff(4,nca) = rbuff(4,ncar)
          buff(5,nca) = rbuff(5,ncar)
          buff(6,nca) = rbuff(6,ncar)
          ibuff(nca)  = irbuff(ncar)
       end do
    end do
    !
  end subroutine add_buffb
!
!----------------------------------------------------------------------
  subroutine add_bucket(abucket,iabucket,isbucket, &
       &                      rbuff,irbuff,irsbuf, &
       &                      max_mvatom,max_mvseg,max_cell, &
       &                      nczdiv,ncydiv,ncxdiv, &
       &                      lcse,ldcell,ncatom,ncseg, &
       &                      ncc_b,ncs_b,nca_b)
!----------------------------------------------------------------------
!
    implicit none
    integer max_mvatom
    integer max_mvseg
    integer max_cell        ! max number of cells on communication buffer.
    integer nczdiv
    integer ncydiv
    integer ncxdiv
    real(8) abucket(6,max_mvatom,nczdiv+2,ncydiv+2,ncxdiv+2)
    integer iabucket(max_mvatom,nczdiv+2,ncydiv+2,ncxdiv+2)
    integer isbucket(max_mvseg,nczdiv+2,ncydiv+2,ncxdiv+2)
    real(8) rbuff(6,max_cell*max_mvatom)
    integer irbuff(max_cell*max_mvatom)
    integer irsbuf(2*max_cell + 1 + max_cell*max_mvseg)

    integer ldcell, ldc
    integer lcse(2,3,ldcell)
    integer ncatom(nczdiv+2,ncydiv+2,ncxdiv+2)
    integer ncseg(nczdiv+2,ncydiv+2,ncxdiv+2)
    integer ncc_b
    integer ncs_b, ncs, ics
    integer nca_b, nca, ica
    integer icz, icy, icx
    integer icz0, icy0, icx0
    integer icz1, icy1, icx1
    !
    do ldc = 1, ldcell
       icz0 = lcse(1,1,ldc)
       icz1 = lcse(2,1,ldc)
       icy0 = lcse(1,2,ldc)
       icy1 = lcse(2,2,ldc)
       icx0 = lcse(1,3,ldc)
       icx1 = lcse(2,3,ldc)
       do icx = icx0, icx1
          do icy = icy0, icy1
             do icz = icz0, icz1
                ! segment data.
                ncc_b = ncc_b + 2                     ! ncc_b
                ncs = ncseg(icz,icy,icx)
                do ics = 1, irsbuf(ncc_b)
                   ncs_b = ncs_b + 1                  ! ncs_b
                   ncs = ncs + 1
                   isbucket(ncs,icz,icy,icx) = irsbuf(ncs_b)
                end do
                ncseg(icz,icy,icx) = ncs
                ! atom data.
                nca = ncatom(icz,icy,icx)
                do ica = 1, irsbuf(ncc_b - 1)
                   nca_b = nca_b + 1                  ! nca_b
                   nca = nca + 1
                   abucket(1,nca,icz,icy,icx) = rbuff(1,nca_b)
                   abucket(2,nca,icz,icy,icx) = rbuff(2,nca_b)
                   abucket(3,nca,icz,icy,icx) = rbuff(3,nca_b)
                   abucket(4,nca,icz,icy,icx) = rbuff(4,nca_b)
                   abucket(5,nca,icz,icy,icx) = rbuff(5,nca_b)
                   abucket(6,nca,icz,icy,icx) = rbuff(6,nca_b)
                   iabucket(nca,icz,icy,icx) = irbuff(nca_b)
                end do
                ncatom(icz,icy,icx) = nca
             end do
          end do
       end do
    end do
    !
  end subroutine add_bucket
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to gather distributed coordinates and velocitiess
!! \author Yoshimichi Andoh, Kensuke Iwahashi
!<
  subroutine pre_record_data
!----------------------------------------------------------------------
    use trajectory_mpi
    use trajectory_org
    use domain, only : lxdiv, lydiv, lzdiv
    use boundary
    implicit none
    integer(4) :: i,nsum
    integer(4) :: i0,k0,i00
    include 'mpif.h'
    integer(4) :: ierr
    real(8),allocatable :: snd(:,:),rcv(:,:)
    integer(4),allocatable :: natmlist(:),natmdisp(:)
    integer(4),allocatable :: nrearrange(:)
    integer(4) :: m2i_tmp(na1cell*lxdiv*lydiv*lzdiv)

    if(nprocs.eq.1) then
!$omp parallel do default(shared) &
!$omp& private(i,i0)
       do i = 1,n
          i0=i2m(i)
          xyz(1,i) = wkxyz(1,i0)
          xyz(2,i) = wkxyz(2,i0)
          xyz(3,i) = wkxyz(3,i0)
          v(1:3,i) = wkv(1:3,i0)
       end do
    else
       allocate(nrearrange(n))
       allocate(snd(6,n))
       allocate(rcv(6,n))
       allocate(natmlist(nprocs),natmdisp(nprocs))

       nselfatm=0
       do k0=1,nselfseg
          do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
             nselfatm=nselfatm+1
             m2i_tmp(nselfatm)=m2i(i0)
          enddo ! i0
       enddo ! k0

       if(nselfatm.gt.na1cell*lxdiv*lydiv*lzdiv)then
          write(0,*) 'ERROR: nselfatm overflowed!', myrank,nselfatm
          call modylas_abort()
       endif

       call mpi_allgather(nselfatm,1,mpi_integer, &
            &                     natmlist,1,mpi_integer, &
            &                     mpi_comm_world,ierr)
       !
       call mpi_barrier(mpi_comm_world,ierr)
       !
       natmdisp(1) = 0
       nsum = natmlist(1)
       do i = 2,nprocs
          natmdisp(i) = natmdisp(i-1)+natmlist(i-1)
          nsum        = nsum + natmlist(i)
       end do
       !
       call mpi_gatherv(m2i_tmp,nselfatm,mpi_integer, &
            &       nrearrange,natmlist,natmdisp,mpi_integer, &
            &                 mpiout,mpi_comm_world,ierr)
       !
!$omp parallel do default(shared) &
!$omp& private(i)
       do i = 1,nprocs
          natmlist(i) = natmlist(i)*6
          natmdisp(i) = natmdisp(i)*6
       end do
       !
       i00=0
       do k0=1,nselfseg
          do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
             i00=i00+1
             snd(1,i00) = wkxyz(1,i0)
             snd(2,i00) = wkxyz(2,i0)
             snd(3,i00) = wkxyz(3,i0)
             snd(4,i00) = wkv(1,i0)
             snd(5,i00) = wkv(2,i0)
             snd(6,i00) = wkv(3,i0)
          end do ! i0
       end do ! k0

       MSG_SEND=nselfatm*6
       call mpi_gatherv(snd,MSG_SEND,mpi_double_precision, &
            &            rcv,natmlist,natmdisp,mpi_double_precision, &
            &            mpiout,mpi_comm_world,ierr)
       !
       if(myrank.eq.mpiout) then
!$omp parallel do default(none) &
!$omp& private(i,i0) &
!$omp& shared(xyz,v,rcv,n,nrearrange)
          do i = 1,n
             i0=nrearrange(i)
             xyz(1,i0) = rcv(1,i)
             xyz(2,i0) = rcv(2,i)
             xyz(3,i0) = rcv(3,i)
             v(1,i0) = rcv(4,i)
             v(2,i0) = rcv(5,i)
             v(3,i0) = rcv(6,i)
          end do
       end if
       deallocate(snd,rcv,nrearrange)
       deallocate(natmlist,natmdisp)
    end if

    call cell_edge()
  end subroutine pre_record_data

end module comm_bound_mod
