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
!!         of grid charge for the PME method.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module which relate to MPI communications of grid charge.
!! \author Yoshimichi Andoh, Kazushi Fujimoto, Shin-ichi Ichikawa
!<
module comm_pme_mod
  !precision module is loaded by default to declare PME subroutine
  use precision, only : wrrp, MPI_PREC_PME_COM
  use margin_sizes
  use mpi_tool
  implicit none

#ifndef PME_XYDIV
!.. for yz-div of PME resulting in small number of cells along z-direction, hence shorter 
!.. inner-most loop length, in real space calc by MTD method.
  integer(4) :: icommx, icommy_ffte, icommz_ffte
  integer(4) :: iminpy,jminpz
  integer(4) :: myrankx,nrankx
#else
!.. for xy-div of PME resulting in maximum number of cells along z-direction, hence maximum
!.. inner-most loop length, in real space calc by MTD method.
  integer(4) :: icommz, icommy_ffte, icommx_ffte
  integer(4) :: kminpy,jminpx
  integer(4) :: myrankz,nrankz
#endif
  integer(4) :: imin,jmin,kmin

  integer(4),parameter :: maxorder=8
  integer(4) :: ngrid=0
  integer(4) :: nfft1=100,nfft2=100,nfft3=100,bsorder=4

  integer(4) :: ngxdiv,ngydiv,ngzdiv,ngxmin,ngymin,ngzmin
#ifndef PME_XYDIV
  integer(4) :: ngxdivmax
  integer(4) :: ngxdivy,ngydivz
#else
  integer(4) :: ngzdivmax
  integer(4) :: ngzdivy,ngydivx
#endif
  real(kind=wrrp),allocatable :: qgrid(:,:,:)

  real(kind=wrrp),allocatable :: buff(:,:)
  real(kind=wrrp),allocatable :: buff2(:,:) ! ya add
  real(kind=wrrp),allocatable :: mpitmp(:,:,:,:)

  real(kind=wrrp),allocatable :: work0(:,:,:)
  real(kind=wrrp),allocatable :: work0_omp(:,:,:,:)
  real(kind=wrrp),allocatable :: work1(:,:,:,:)
  real(kind=wrrp),allocatable :: work2(:,:,:,:)
#ifdef DBG_PME_WNCOMP
  real(kind=wrrp),allocatable :: work_dbg(:,:,:,:)
#endif

contains

!<
!========================================================================
!
!.... create communicator                                          ....
!
!========================================================================
!>
!! \brief  Subroutine to initialize comm_qgrid
!! \author Yoshimichi Andoh, Kazushi Fujimoto, Shin-ichi Ichikawa
!<
  SUBROUTINE comm_create()
    use mpi_3d_grid
    implicit none
    INCLUDE 'mpif.h'
    integer(4) :: ierr
    integer :: mpi_comm, icomm_0, color, key
    integer :: myrank_0, nrank_0

    mpi_comm = MPI_COMM_WORLD
    !
    !.... Create communicator in the X-direction
    !
#ifndef PME_XYDIV
    color = ipy * 100000 + ipz
#else
    color = ipy * 100000 + ipx
#endif
    key = myrank
#ifndef PME_XYDIV
    CALL MPI_COMM_SPLIT(mpi_comm, color, key, icommx, ierr)
    IF (ipx == 0) THEN
#else
    CALL MPI_COMM_SPLIT(mpi_comm, color, key, icommz, ierr)
    IF (ipz == 0) THEN
#endif
       color = 1
    ELSE
       color = MPI_UNDEFINED
    END IF

    key = myrank
    CALL MPI_COMM_SPLIT(mpi_comm, color, key, icomm_0, ierr)
#ifndef PME_XYDIV
    IF (ipx == 0) THEN
#else
    IF (ipz == 0) THEN
#endif
       CALL MPI_COMM_RANK(icomm_0, myrank_0, ierr)
       CALL MPI_COMM_SIZE(icomm_0, nrank_0, ierr)
       !
       !.... Create communicator in the Y-direction for FFTE
       !
#ifndef PME_XYDIV
       color = myrank_0/npy
#else
       color = ipx
#endif
       key   = myrank_0
       CALL MPI_COMM_SPLIT(icomm_0, color, key, icommy_ffte, ierr)
       !
#ifndef PME_XYDIV
       !.... Create communicator in the Z-direction for FFTE
       color = mod(myrank_0,npy)
       key   = myrank_0
       CALL MPI_COMM_SPLIT(icomm_0, color, key, icommz_ffte, ierr)
#else
       !.... Create communicator in the X-direction for FFTE
       color = ipy
       key   = myrank_0
       CALL MPI_COMM_SPLIT(icomm_0, color, key, icommx_ffte, ierr)
#endif
       !
       !
       !.... FFTE Initializing
       !
#ifndef PME_XYDIV
       CALL PZFFT3DV(work1, work2, nfft1, nfft2, nfft3, & !! ya-modify
            &                 icommy_ffte, icommz_ffte, npy ,npz, 0)
#else
       CALL PZFFT3DV(work1, work2, nfft3, nfft2, nfft1, &
            &                 icommy_ffte, icommx_ffte, npy ,npx, 0)
#endif
    END IF
    !     CALL MPI_COMM_FREE(icomm_0, ierr)  !! ya-modify

#ifndef PME_XYDIV
    CALL MPI_COMM_RANK(icommx, myrankx, ierr)
    CALL MPI_COMM_SIZE(icommx, nrankx, ierr)
#else
    CALL MPI_COMM_RANK(icommz, myrankz, ierr)
    CALL MPI_COMM_SIZE(icommz, nrankz, ierr)
#endif

  END SUBROUTINE comm_create
!========================================================================
!
!.... qgrid communication subroutine.                               ....
!
!========================================================================
!>
!! \brief  Subroutine to communicate grid data
!! \author Yoshimichi Andoh, Kazushi Fujimoto, Shin-ichi Ichikawa
!<
  SUBROUTINE comm_qgrid(ng1divmax,ng1div,ng2div,ng3div,ng1min,ng2min,ng3min,icomm_scat)

    use mpi_3d_grid
    implicit none
    INCLUDE 'mpif.h'
    integer(4) :: ierr
    integer :: ip3_dest, ip2_dest, ip1_dest
    integer :: ip3_src, ip2_src, ip1_src
    integer :: remgrid, pregrid, inds, indr
    integer :: itr, nitr, iadd
    integer :: ig1,ig2,ig3
    integer :: ng1divmax
    integer :: ng1div, ng2div, ng3div
    integer :: ng1min, ng2min, ng3min
    integer :: icomm_scat

    integer :: ista, ista_r
    integer :: mpi_comm, istatus(MPI_STATUS_SIZE)
    integer :: itag1=10, itag2=11, itag3=12
#ifdef _TIMER_
    real(8) :: tm_gather0, tm_scatter0, tm_fftd0, tm_ffti0
    real(8) :: tm_gather1, tm_scatter1, tm_fftd1, tm_ffti1
    real(8) :: tm_gatherb0, tm_scatterb0
    real(8) :: tm_gatherb1, tm_scatterb1
    real(8) :: tm_x0, tm_y0, tm_z0, tm_xb0, tm_tmp
    real(8) :: tm_x1, tm_y1, tm_z1, tm_xb1
    real(8) :: tm_xs0, tm_ys0, tm_zs0, tm_yb0
    real(8) :: tm_xs1, tm_ys1, tm_zs1, tm_yb1
#endif
    integer(4) :: bso

    bso=bsorder+bso_mgn

    mpi_comm = MPI_COMM_WORLD

#ifdef _TIMER_
    call mpi_barrier(mpi_comm,ierr)
    tm_scatter0 = MPI_WTIME()
#endif
    CALL MPI_SCATTER(mpitmp, ng1divmax*ng2div*ng3div, &
         &                 MPI_PREC_PME_COM, &
         &                 work0,  ng1divmax*ng2div*ng3div, &
         &                 MPI_PREC_PME_COM, 0, icomm_scat, ierr)

#ifdef _TIMER_
    tm_scatter1 = MPI_WTIME() - tm_scatter0
    tm_scatterb0 = MPI_WTIME()
#endif
    DO ig3 = 1, ng3div
       DO ig2 = 1, ng2div
          DO ig1 = 1, ng1divmax
             qgrid(bso+ig1,bso+ig2,bso+ig3) = work0(ig1,ig2,ig3)
          END DO
       END DO
    END DO
#ifdef _TIMER_
    tm_scatterb1 = MPI_WTIME() - tm_scatterb0
#endif

    !
    !.... Coordinate +X/+Z
    !

#ifdef _TIMER_
    tm_xb1 = 0.0d0
    tm_xs1 = 0.0d0
    call mpi_barrier(mpi_comm, ierr)
    tm_x0 = MPI_WTIME()
#endif
#ifndef PME_XYDIV
    ip1_dest  = ipz*npx*npy + ipy*npx + mod(ipx+1-1/npx+npx,npx)
    ip1_src   = ipz*npx*npy + ipy*npx + mod(ipx-1+1/npx+npx,npx)
#else
    ip1_dest = mod(ipz+1+npz,npz)*npx*npy + ipy*npx + ipx
    ip1_src  = mod(ipz-1+npz,npz)*npx*npy + ipy*npx + ipx
#endif

    nitr = (bso - 1)/ng1min + 1
    remgrid = bso
    ngrid = 0
    DO itr = 1, nitr
       pregrid = ngrid
       IF (remgrid > ng1min) THEN
          ngrid = ng1min
       ELSE
          ngrid = remgrid  ! = 4
       END IF
       remgrid = remgrid - ngrid  ! =0
       inds = mod(itr+1,2) + 1
       indr = mod(itr  ,2) + 1
       IF (itr == 1) THEN
#ifdef _TIMER_D_
          tm_xb0 = MPI_WTIME()
#endif
!OCL NORECURRENCE
          DO ig3 = 1, ng3div
             DO ig2 = 1, ng2div
                DO ig1 = 1, ngrid
                   iadd = ig1+ngrid*(ig2-1)+ngrid*ng2div*(ig3-1) ! 1,2,...
                   buff(iadd,inds) = qgrid(bso+ng1div-ngrid+ig1,bso+ig2,bso+ig3)
                END DO
             END DO
          END DO
#ifdef _TIMER_D_
          tm_xb1 = tm_xb1 + MPI_WTIME() - tm_xb0
          call mpi_barrier(mpi_comm, ierr)
          tm_xs0 = MPI_WTIME()
#endif
          CALL mpi_sendrecv(buff(1,inds), ngrid*ng2div*ng3div, &
               &                        MPI_PREC_PME_COM, ip1_dest, myrank, &
               &                        buff(1,indr), ngrid*ng2div*ng3div, &
               &                        MPI_PREC_PME_COM, ip1_src, ip1_src, &
               &                        mpi_comm, istatus, ierr )
#ifdef _TIMER_D_
          tm_xs1 = tm_xs1 + MPI_WTIME() - tm_xs0
#endif
          DO ig3 = 1, ng3div
             DO ig2 = 1, ng2div
                DO ig1 = 1, ngrid
                   iadd = ig1+ngrid*(ig2-1)+ngrid*ng2div*(ig3-1)
                   qgrid(remgrid+ig1,bso+ig2,bso+ig3) = buff(iadd,indr)
                END DO
             END DO
          END DO
       ELSE
          IF (((nfft1/ng1min) /= ng1div).or.(pregrid /= ngrid)) then
             !                       write(*,*) 'IF +X eneterd'
#ifdef _TIMER_D_
             tm_xb0 = MPI_WTIME()
#endif
!OCL NORECURRENCE
             DO ig3 = 1, ng3div
                DO ig2 = 1, ng2div
                   DO ig1 = 1, ngrid
                      iadd = ig1+ngrid*(ig2-1)+ngrid*ng2div*(ig3-1)
                      buff(iadd,inds) = qgrid(bso+ng1div-ng1min*(itr-1)-ngrid+ig1, bso+ig2,bso+ig3)
                   END DO
                END DO
             END DO
#ifdef _TIMER_D_
             tm_xb1 = tm_xb1 + MPI_WTIME() - tm_xb0
#endif
          END IF
#ifdef _TIMER_D_
          call mpi_barrier(mpi_comm, ierr)
          tm_xs0 = MPI_WTIME()
#endif
          CALL mpi_sendrecv(buff(1,inds), ngrid*ng2div*ng3div, &
               &                        MPI_PREC_PME_COM, ip1_dest, myrank, &
               &                        buff(1,indr), ngrid*ng2div*ng3div, &
               &                        MPI_PREC_PME_COM, ip1_src, ip1_src, &
               &                        mpi_comm, istatus, ierr )
#ifdef _TIMER_
          tm_xs1 = tm_xs1 + MPI_WTIME() - tm_xs0
#endif
          DO ig3 = 1, ng3div
             DO ig2 = 1, ng2div
                DO ig1 = 1, ngrid
                   iadd = ig1+ngrid*(ig2-1)+ngrid*ng2div*(ig3-1)
                   qgrid(remgrid+ig1,bso+ig2,bso+ig3) = buff(iadd,indr)
                END DO
             END DO
          END DO
       END IF
    END DO
#ifdef _TIMER_
    tm_x1 = MPI_WTIME() - tm_x0
#endif

    !
    !.... Coordinate -X  (ya+) / -Z
    !

#ifndef PME_XYDIV
    ip1_src   = ipz*npx*npy + ipy*npx + mod(ipx+1-1/npx+npx,npx)
    ip1_dest  = ipz*npx*npy + ipy*npx + mod(ipx-1+1/npx+npx,npx)
#else
    ip1_dest = mod(ipz-1+npz,npz)*npx*npy + ipy*npx + ipx
    ip1_src  = mod(ipz+1+npz,npz)*npx*npy + ipy*npx + ipx
#endif

    nitr = (bso - 1)/ng1min + 1
    remgrid = bso
    ngrid = 0
    DO itr = 1, nitr
       pregrid = ngrid
       IF (remgrid > ng1min) THEN
          ngrid = ng1min
       ELSE
          ngrid = remgrid  ! = 4
       END IF
       remgrid = remgrid - ngrid  ! =0
       inds = mod(itr+1,2) + 1
       indr = mod(itr  ,2) + 1
       IF (itr == 1) THEN
#ifdef _TIMER_D_
#endif
!OCL NORECURRENCE
          DO ig3 = 1, ng3div
             DO ig2 = 1, ng2div
                DO ig1 = 1, ngrid
                   iadd = ig1+ngrid*(ig2-1)+ngrid*ng2div*(ig3-1) ! 1,2,...
                   buff2(iadd,inds) =  qgrid(bso+remgrid+ig1,bso+ig2,bso+ig3)
                END DO
             END DO
          END DO
#ifdef _TIMER_D_
#endif
          CALL mpi_sendrecv(buff2(1,inds), ngrid*ng2div*ng3div, &
               &                        MPI_PREC_PME_COM, ip1_dest, myrank, &
               &                        buff2(1,indr), ngrid*ng2div*ng3div, &
               &                        MPI_PREC_PME_COM, ip1_src, ip1_src, &
               &                        mpi_comm, istatus, ierr )
#ifdef _TIMER_D_
#endif
          DO ig3 = 1, ng3div
             DO ig2 = 1, ng2div
                DO ig1 = 1, ngrid
                   iadd = ig1+ngrid*(ig2-1)+ngrid*ng2div*(ig3-1)
                   qgrid(bso+ng1div+remgrid+ig1,bso+ig2,bso+ig3) = buff2(iadd,indr)
                END DO
             END DO
          END DO
       ELSE
          IF (((nfft1/ng1min) /= ng1div).or.(pregrid /= ngrid)) then
#ifdef _TIMER_D_
#endif
!OCL NORECURRENCE
             DO ig3 = 1, ng3div
                DO ig2 = 1, ng2div
                   DO ig1 = 1, ngrid
                      iadd = ig1+ngrid*(ig2-1)+ngrid*ng2div*(ig3-1)
                      buff2(iadd,inds) = qgrid(bso+ng1min*(itr-1)+ig1,bso+ig2,bso+ig3)
                   END DO
                END DO
             END DO
#ifdef _TIMER_D_
#endif
          END IF
#ifdef _TIMER_D_
#endif
          CALL mpi_sendrecv(buff2(1,inds), ngrid*ng2div*ng3div, &
               &                        MPI_PREC_PME_COM, ip1_dest, myrank, &
               &                        buff2(1,indr), ngrid*ng2div*ng3div, &
               &                        MPI_PREC_PME_COM, ip1_src, ip1_src, &
               &                        mpi_comm, istatus, ierr )
#ifdef _TIMER_
#endif
          DO ig3 = 1, ng3div
             DO ig2 = 1, ng2div
                DO ig1 = 1, ngrid
                   iadd = ig1+ngrid*(ig2-1)+ngrid*ng2div*(ig3-1)
                   qgrid(bso+ng1div+remgrid+ig1,bso+ig2,bso+ig3) = buff2(iadd,indr)
                END DO
             END DO
          END DO
       END IF
    END DO
#ifdef _TIMER_
#endif

    !
    !     coordinate +Y
    !

#ifdef _TIMER_
    tm_ys1 = 0.0d0
    tm_yb1 = 0.0d0
    call mpi_barrier(mpi_comm, ierr)
    tm_y0 = MPI_WTIME()
#endif
    ip2_dest = ipz*npx*npy + mod(ipy+1-1/npy+npy,npy)*npx + ipx
    ip2_src  = ipz*npx*npy + mod(ipy-1+1/npy+npy,npy)*npx + ipx

    nitr = (bso - 1)/ng2min + 1
    remgrid  = bso
    ngrid  = 0

    DO itr = 1, nitr
       pregrid = ngrid
       IF (remgrid > ng2min) THEN
          ngrid = ng2min
       ELSE
          ngrid = remgrid
       END IF
       remgrid = remgrid - ngrid
       inds = mod(itr+1,2) + 1
       indr = mod(itr  ,2) + 1

       IF (itr == 1) THEN
#ifdef _TIMER_D_
          tm_yb0 = MPI_WTIME()
#endif
!OCL NORECURRENCE
          DO ig3 = 1, ng3div
             DO ig2 = 1, ngrid
                DO ig1 = 1, bso+ng1div+bso
                   iadd = ig1+(bso+ng1div+bso)*(ig2-1) + (bso+ng1div+bso)*ngrid*(ig3-1)
                   buff(iadd,inds) = qgrid(ig1,bso+ng2div-ngrid+ig2,bso+ig3)
                END DO
             END DO
          END DO
#ifdef _TIMER_D_
          tm_yb1 = tm_yb1 + MPI_WTIME() - tm_yb0
          call mpi_barrier(mpi_comm, ierr)
          tm_ys0 = MPI_WTIME()
#endif
          CALL MPI_SENDRECV(buff(1,inds), &
               &                        (bso+ng1div+bso)*ngrid*ng3div, &
               &                        MPI_PREC_PME_COM,ip2_dest, myrank, &
               &                        buff(1,indr), &
               &                        (bso+ng1div+bso)*ngrid*ng3div, &
               &                        MPI_PREC_PME_COM,ip2_src, ip2_src, &
               &                        mpi_comm, istatus, ierr )
#ifdef _TIMER_D_
          tm_tmp = MPI_WTIME()
          tm_ys1 = tm_ys1 + tm_tmp - tm_ys0
          tm_yb0 = tm_tmp
#endif
          DO ig3 = 1, ng3div
             DO ig2 = 1, ngrid
                DO ig1 = 1, bso+ng1div+bso
                   iadd = ig1+(bso+ng1div+bso)*(ig2-1) + (bso+ng1div+bso)*ngrid*(ig3-1)
                   qgrid(ig1,remgrid+ig2,bso+ig3) = buff(iadd,indr)
                END DO
             END DO
          END DO
#ifdef _TIMER_D_
          tm_yb1 = tm_yb1 + MPI_WTIME() - tm_yb0
#endif
       ELSE
          !                       write(*,*) 'IF +Y eneterd'
#ifdef _TIMER_D_
          tm_yb0 = MPI_WTIME()
#endif
          IF (((nfft2/npy) /= ng2div).or.(pregrid /= ngrid)) then
!OCL NORECURRENCE
             DO ig3 = 1, ng3div
                DO ig2 = 1, ngrid
                   DO ig1 = 1, bso+ng1div+bso
                      iadd = ig1+(bso+ng1div+bso)*(ig2-1) + (bso+ng1div+bso)*ngrid*(ig3-1)
                      buff(iadd,inds) = qgrid(ig1, bso+ng2div-ng2min*(itr-1)-ngrid+ig2,bso+ig3)
                   END DO
                END DO
             END DO
          END IF
#ifdef _TIMER_D_
          tm_yb1 = tm_yb1 + MPI_WTIME() - tm_yb0
          call mpi_barrier(mpi_comm, ierr)
          tm_ys0 = MPI_WTIME()
#endif
          CALL MPI_SENDRECV(buff(1,inds), &
               &                        (bso+ng1div+bso)*ngrid*ng3div, &
               &                        MPI_PREC_PME_COM,ip2_dest, myrank, &
               &                        buff(1,indr), &
               &                        (bso+ng1div+bso)*ngrid*ng3div, &
               &                        MPI_PREC_PME_COM,ip2_src, ip2_src, &
               &                        mpi_comm, istatus, ierr )
#ifdef _TIMER_D_
          tm_tmp = MPI_WTIME()
          tm_ys1 = tm_ys1 + tm_tmp - tm_ys0
          tm_yb0 = tm_tmp
#endif
          DO ig3 = 1, ng3div
             DO ig2 = 1, ngrid
                DO ig1 = 1, bso+ng1div+bso
                   iadd = ig1+(bso+ng1div+bso)*(ig2-1) + (bso+ng1div+bso)*ngrid*(ig3-1)
                   qgrid(ig1,remgrid+ig2,bso+ig3) = buff(iadd,indr)
                END DO
             END DO
          END DO
#ifdef _TIMER_D_
          tm_yb1 = tm_yb1 + MPI_WTIME() - tm_yb0
#endif
       END IF

    END DO
#ifdef _TIMER_
    tm_y1 = MPI_WTIME() - tm_y0
#endif

    !
    !     coordinate -Y (ya)
    !

#ifdef _TIMER_
#endif
    ip2_dest = ipz*npx*npy + mod(ipy-1+1/npy+npy,npy)*npx + ipx
    ip2_src  = ipz*npx*npy + mod(ipy+1-1/npy+npy,npy)*npx + ipx
    nitr = (bso - 1)/ng2min + 1
    remgrid  = bso
    ngrid  = 0

    DO itr = 1, nitr
       pregrid = ngrid ! not used
       IF (remgrid > ng2min) THEN
          ngrid = ng2min
       ELSE
          ngrid = remgrid
       END IF
       remgrid = remgrid - ngrid
       inds = mod(itr+1,2) + 1
       indr = mod(itr  ,2) + 1

       IF (itr == 1) THEN
#ifdef _TIMER_D_
#endif
!OCL NORECURRENCE
          DO ig3 = 1, ng3div
             DO ig2 = 1, ngrid
                DO ig1 = 1, bso+ng1div+bso
                   iadd = ig1+(bso+ng1div+bso)*(ig2-1) + (bso+ng1div+bso)*ngrid*(ig3-1)
                   buff2(iadd,inds) = qgrid(ig1,bso+remgrid+ig2,bso+ig3)
                END DO
             END DO
          END DO
#ifdef _TIMER_D_
#endif
          CALL MPI_SENDRECV(buff2(1,inds), &
               &                        (bso+ng1div+bso)*ngrid*ng3div, &
               &                        MPI_PREC_PME_COM,ip2_dest, myrank, &
               &                        buff2(1,indr), &
               &                        (bso+ng1div+bso)*ngrid*ng3div, &
               &                        MPI_PREC_PME_COM,ip2_src, ip2_src, &
               &                        mpi_comm, istatus, ierr )
#ifdef _TIMER_D_
#endif
          DO ig3 = 1, ng3div
             DO ig2 = 1, ngrid
                DO ig1 = 1, bso+ng1div+bso
                   iadd = ig1+(bso+ng1div+bso)*(ig2-1) + (bso+ng1div+bso)*ngrid*(ig3-1)
                   qgrid(ig1,bso+ng2div+remgrid+ig2,bso+ig3) = buff2(iadd,indr)
                END DO
             END DO
          END DO
#ifdef _TIMER_D_
#endif
       ELSE
#ifdef _TIMER_D_
#endif
          IF (((nfft2/npy) /= ng2div).or.(pregrid /= ngrid)) then
!OCL NORECURRENCE
             DO ig3 = 1, ng3div
                DO ig2 = 1, ngrid
                   DO ig1 = 1, bso+ng1div+bso
                      iadd = ig1+(bso+ng1div+bso)*(ig2-1) + (bso+ng1div+bso)*ngrid*(ig3-1)
                      buff2(iadd,inds) = qgrid(ig1,bso+ng2min*(itr-1)+ig2,bso+ig3)
                   END DO
                END DO
             END DO
          END IF
#ifdef _TIMER_D_
#endif
          CALL MPI_SENDRECV(buff2(1,inds), &
               &                        (bso+ng1div+bso)*ngrid*ng3div, &
               &                        MPI_PREC_PME_COM,ip2_dest, myrank, &
               &                        buff2(1,indr), &
               &                        (bso+ng1div+bso)*ngrid*ng3div, &
               &                        MPI_PREC_PME_COM,ip2_src, ip2_src, &
               &                        mpi_comm, istatus, ierr )
#ifdef _TIMER_D_
#endif
          DO ig3 = 1, ng3div
             DO ig2 = 1, ngrid
                DO ig1 = 1, bso+ng1div+bso
                   iadd = ig1+(bso+ng1div+bso)*(ig2-1) + (bso+ng1div+bso)*ngrid*(ig3-1)
                   qgrid(ig1,bso+ng2div+remgrid+ig2,bso+ig3) = buff2(iadd,indr)
                END DO
             END DO
          END DO
#ifdef _TIMER_D_
#endif
       END IF

    END DO
#ifdef _TIMER_
#endif

    !
    !     coordinate +Z / +X
    !

#ifdef _TIMER_
    tm_zs1 = 0.0d0
    call mpi_barrier(mpi_comm, ierr)
    tm_z0 = MPI_WTIME()
#endif
#ifndef PME_XYDIV
    ip3_dest = mod(ipz+1+npz,npz)*npx*npy + ipy*npx + ipx
    ip3_src  = mod(ipz-1+npz,npz)*npx*npy + ipy*npx + ipx
#else
    ip3_dest  = ipz*npx*npy + ipy*npx + mod(ipx+1-1/npx+npx,npx)
    ip3_src   = ipz*npx*npy + ipy*npx + mod(ipx-1+1/npx+npx,npx)
#endif
    nitr = (bso - 1)/ng3min + 1
    remgrid = bso
    ngrid = 0
    DO itr = 1, nitr

       IF (remgrid > ng3min) THEN
          ngrid = ng3min
       ELSE
          ngrid = remgrid
       END IF
       remgrid = remgrid - ngrid

       IF (itr == 1) THEN
          ista   = bso + ng3div - ngrid + 1
          ista_r = bso          - ngrid + 1
       ELSE
          ista   = ista   - ngrid
          ista_r = ista_r - ngrid
          !                       write(*,*) 'IF +Z eneterd'
       END IF
#ifdef _TIMER_D_
       call mpi_barrier(mpi_comm, ierr)
       tm_zs0 = MPI_WTIME()
#endif
       CALL MPI_SENDRECV(qgrid(1,1,ista), &
            &                     (bso+ng1div+bso)*(bso+ng2div+bso)*ngrid, &
            &                     MPI_PREC_PME_COM,ip3_dest, myrank, &
            &                     qgrid(1,1,ista_r), &
            &                     (bso+ng1div+bso)*(bso+ng2div+bso)*ngrid, &
            &                     MPI_PREC_PME_COM,ip3_src, ip3_src, &
            &                     mpi_comm, istatus, ierr )
#ifdef _TIMER_D_
       tm_zs1 = tm_zs1 + MPI_WTIME() - tm_zs0
#endif
    END DO
#ifdef _TIMER_
    tm_z1 = MPI_WTIME() - tm_z0
#endif

    !
    !     coordinate -Z / -X
    !

#ifdef _TIMER_
#endif
#ifndef PME_XYDIV
    ip3_dest = mod(ipz-1+npz,npz)*npx*npy + ipy*npx + ipx
    ip3_src  = mod(ipz+1+npz,npz)*npx*npy + ipy*npx + ipx
#else
    ip3_src   = ipz*npx*npy + ipy*npx + mod(ipx+1-1/npx+npx,npx)
    ip3_dest  = ipz*npx*npy + ipy*npx + mod(ipx-1+1/npx+npx,npx)
#endif
    nitr = (bso - 1)/ng3min + 1
    remgrid = bso
    ngrid = 0
    DO itr = 1, nitr

       IF (remgrid > ng3min) THEN
          ngrid = ng3min
       ELSE
          ngrid = remgrid
       END IF
       remgrid = remgrid - ngrid

       IF (itr == 1) THEN
          ista   = bso          + 1
          ista_r = bso + ng3div + 1
       ELSE
          ista   = ista   + ngrid
          ista_r = ista_r + ngrid
       END IF
#ifdef _TIMER_D_
#endif
       CALL MPI_SENDRECV(qgrid(1,1,ista), &
            &                     (bso+ng1div+bso)*(bso+ng2div+bso)*ngrid, &
            &                     MPI_PREC_PME_COM,ip3_dest, myrank, &
            &                     qgrid(1,1,ista_r), &
            &                     (bso+ng1div+bso)*(bso+ng2div+bso)*ngrid, &
            &                     MPI_PREC_PME_COM,ip3_src, ip3_src, &
            &                     mpi_comm, istatus, ierr )
#ifdef _TIMER_D_
#endif
    END DO
#ifdef _TIMER_
#endif


#ifdef _TIMER_
    tm_val( 1,jjj) = tm_fftd1
    tm_val( 2,jjj) = tm_ffti1
    tm_val( 3,jjj) = tm_gatherb1
    tm_val( 4,jjj) = tm_gather1
    tm_val( 5,jjj) = tm_scatter1
    tm_val( 6,jjj) = tm_scatterb1
    tm_val( 7,jjj) = tm_xb1
    tm_val( 8,jjj) = tm_xs1
    tm_val( 9,jjj) = tm_ys1
    tm_val(10,jjj) = tm_zs1
    tm_val(11,jjj) = tm_x1
    tm_val(12,jjj) = tm_y1
    tm_val(13,jjj) = tm_z1
    tm_val(14,jjj) = tm_yb1
#endif
  END SUBROUTINE comm_qgrid
!----------------------------------------------------------------------
  subroutine fmod_pmewald_bsorder(value)
    use mpi_tool
    implicit none
    real(8), intent(in) :: value
    bsorder = value
    if(bsorder.gt.maxorder)then
       write(0,*) 'ERROR: bsorder > maxorder'
       call modylas_abort()
    endif
  end subroutine fmod_pmewald_bsorder
!----------------------------------------------------------------------
  subroutine fmod_pmewald_nfft1(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    nfft1 = ivalue
  end subroutine fmod_pmewald_nfft1

  subroutine fmod_pmewald_nfft2(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    nfft2 = ivalue
  end subroutine fmod_pmewald_nfft2
!----------------------------------------------------------------------
  subroutine fmod_pmewald_nfft3(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    nfft3 = ivalue
  end subroutine fmod_pmewald_nfft3

end module comm_pme_mod
