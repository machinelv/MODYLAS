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
!!         of atom forces in the MTD method. \n
!<
!! [Ref] Y.Andoh, S.Ichikawa, T Sakashita, N Yoshii, S Okazaki,
!!       J. Comput. Chem., 42, 1073-1087 (2021).
!----------------------------------------------------------------------
!>
!! \brief  Module which relate to MPI communications of atom forces.
!! \author Yoshimichi Andoh, Shin-ichi Ichikawa
!<
module comm_direct2_dr
  use dr_cntl, only : nbd, nbd2, nadj !, nomp
  use domain, only  : lxdiv, lydiv, lzdiv
  use subcell, only : tag, na_per_cell
  use trajectory_mpi, only : naline, narea, na1cell, nadirect
  use mpi_3d_grid
#ifndef FJ_RDMA
  use mpi_tool, only : myrank
#else
  use mpi_tool
#endif
  use omp_lib ! YA
#ifdef PROFILE_COMM
  use profile_comm
#endif
  real(8),allocatable,dimension(:,:) :: buffp
  real(8),allocatable,dimension(:,:) :: rbuffp
  real(8),allocatable,dimension(:,:) :: buffm
  real(8),allocatable,dimension(:,:) :: rbuffm
  integer(4) :: max_cell_2
  integer :: MSG_SEND,MSG_RECV
#ifdef FJ_RDMA
  integer(8) :: buffp_laddr
  integer(8) :: buffm_laddr
  integer(8), pointer :: rbuffp_raddr(:)
  integer(8), pointer :: rbuffm_raddr(:)
#endif
contains

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to initialize arrays used in MPI communications.
!! \author Yoshimichi Andoh, Shin-ichi Ichikawa
!<
  subroutine init_comm_direct_2_dr()
!----------------------------------------------------------------------
    implicit none
#ifdef FJ_RDMA
    type(c_ptr) :: rbuffp_cptr
    type(c_ptr) :: rbuffm_cptr
#endif
    max_cell_2 = lydiv*lxdiv
    max_cell_2 = max(max_cell_2, (lzdiv+nbd2))
    max_cell_2 = max(max_cell_2, (lzdiv+nbd2)*(lydiv+nbd2))
    max_cell_2 = max(max_cell_2, (lzdiv+nbd2)*(lxdiv+nbd2))
    allocate(buffp (3, max_cell_2*nbd*na1cell))
    allocate(rbuffp(3, max_cell_2*nbd*na1cell))
    allocate(buffm (3, max_cell_2*nbd*na1cell))
    allocate(rbuffm(3, max_cell_2*nbd*na1cell))

#ifdef FJ_RDMA
    call rdma_register_addr(buffp,  (3*max_cell_2*nbd*na1cell)*8)
    call rdma_register_addr(rbuffp, (3*max_cell_2*nbd*na1cell)*8)
    call rdma_register_addr(buffm,  (3*max_cell_2*nbd*na1cell)*8)
    call rdma_register_addr(rbuffm, (3*max_cell_2*nbd*na1cell)*8)
    !
    buffp_laddr = rdma_get_laddr(buffp)
    buffm_laddr = rdma_get_laddr(buffm)
    !
    rbuffp_cptr = rdma_get_raddr(rbuffp)
    rbuffm_cptr = rdma_get_raddr(rbuffm)
    call c_f_pointer(rbuffp_cptr, fptr=rbuffp_raddr, shape=[nprocs])
    call c_f_pointer(rbuffm_cptr, fptr=rbuffm_raddr, shape=[nprocs])
#endif

    buffp(:,:) =0d0
    rbuffp(:,:)=0d0
    buffm(:,:) =0d0
    rbuffm(:,:)=0d0

  end subroutine init_comm_direct_2_dr
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to perform MPI communications of atom forces.
!! \author Yoshimichi Andoh, Shin-ichi Ichikawa
!<
  subroutine comm_direct_2_dr(wk_f) !! without dimension for nomp
!----------------------------------------------------------------------
    use md_condition, only : mdstep
#include "timing.h90"
    implicit none

    INCLUDE 'mpif.h'
    real(8) wk_f(3,nadirect)
    integer ipz_pdest, ipy_pdest, ipx_pdest
    integer ipz_psrc, ipy_psrc, ipx_psrc
    integer ipz_mdest, ipy_mdest, ipx_mdest
    integer ipz_msrc, ipy_msrc, ipx_msrc
    integer iczb, icyb, icz, icy, icx
    integer icz0, icz1
    integer icy0
    integer iczp
    integer iczm
    integer icyp
    integer icym
    integer icxp
    integer icxm
    integer iczbp
    integer iczbm
    integer icybp
    integer icybm
    integer icxbp
    integer icxbm
    integer ica
    integer icasp, icarp
    integer icasm, icarm
    integer nca
    integer ncap, ncarp
    integer ncam, ncarm
    integer istatus(mpi_status_size, 4), ierr

#ifndef SYNC_COM
    integer,dimension(4) :: irq
    integer nrq
#endif

#ifdef FJ_RDMA
    type(c_ptr) :: wk_f_cptr
    integer(8) :: wk_f_laddr
    integer(8),pointer :: wk_f_raddr(:)
    
    wk_f_laddr = rdma_get_laddr(wk_f)
    wk_f_cptr  = rdma_get_raddr(wk_f)
    call c_f_pointer(wk_f_cptr, fptr=wk_f_raddr, shape=[nprocs])
#endif

    !----- common parameters for coordinate communication. -----
    !    ipx=mod(myrank,npx)
    !    ipy=mod((myrank-ipx)/npx,npy)
    !    ipz=mod((myrank-ipx-ipy*npx)/(npx*npy),npz)
    !    narea = na1cell * (ncellz/npz + 4) * (ncelly/npy + 4)
    !    naline = na1cell * (ncellz/npz + 4)

    !
    !-----  dr method force communication code starts here. ------
    !
    !     force +X
    ipx_pdest  = ipz*npx*npy + ipy*npx + mod(ipx+1-1/npx+npx,npx)
    ipx_psrc   = ipz*npx*npy + ipy*npx + mod(ipx-1+1/npx+npx,npx)
    !     force -X
    ipx_mdest  = ipz*npx*npy + ipy*npx + mod(ipx-1+1/npx+npx,npx)
    ipx_msrc   = ipz*npx*npy + ipy*npx + mod(ipx+1-1/npx+npx,npx)

    icy0 = 1 - nbd

    icxp = lxdiv + 1
    icxbp = icxp - lxdiv
    ncap = narea
    icasp = tag(1,icy0,icxp) - na1cell
    icarp = tag(1,icy0,icxbp) - na1cell
    icxm = 1 - nbd
    icxbm = icxm + lxdiv
    ncam = narea
    icasm = tag(1,icy0,icxm) - na1cell
    icarm = tag(1,icy0,icxbm) - na1cell

#ifdef PROFILE_COMM
    call mpi_barrier(mpi_comm_world, ierr)
    call start_profile(comm_direct2_dr_id)
#endif
#ifdef SYNC_COM
    MSG_SEND=3*ncap
    MSG_RECV=3*ncap
    call mpi_sendrecv(wk_f(1,icasp), MSG_SEND, MPI_DOUBLE_PRECISION, &
         &            ipx_pdest, myrank, &
         &            rbuffp, MSG_RECV, MPI_DOUBLE_PRECISION, &
         &            ipx_psrc, ipx_psrc, mpi_comm_world, &
         &            istatus, ierr )
    MSG_SEND=3*ncam
    MSG_RECV=3*ncam
    call mpi_sendrecv(wk_f(1,icasm), MSG_SEND, MPI_DOUBLE_PRECISION, &
         &            ipx_mdest, myrank, &
         &            rbuffm, MSG_RECV, MPI_DOUBLE_PRECISION, &
         &            ipx_msrc, ipx_msrc, mpi_comm_world, &
         &            istatus, ierr )
#elif FJ_RDMA
    call rdma_post(ipx_psrc)
    if (ipx_psrc .ne. ipx_msrc) then
       call rdma_post(ipx_msrc)
       call rdma_wait(ipx_mdest)
    endif
    call rdma_wait(ipx_pdest)
    !
    if (ipx_pdest .ne. ipx_mdest) then
       MSG_SEND=3*ncap
       call rdma_put_post(ipx_pdest, rbuffp_raddr(ipx_pdest+1), wk_f_laddr+(icasp-1)*3*8, MSG_SEND*8)
       MSG_SEND=3*ncam
       call rdma_put_post(ipx_mdest, rbuffm_raddr(ipx_mdest+1), wk_f_laddr+(icasm-1)*3*8, MSG_SEND*8)
       call rdma_wait(ipx_psrc)
       call rdma_wait(ipx_msrc)
    else
       MSG_SEND=3*ncap
       call rdma_put(ipx_pdest, rbuffp_raddr(ipx_pdest+1),      wk_f_laddr+(icasp-1)*3*8, MSG_SEND*8)
       MSG_SEND=3*ncam
       call rdma_put_post(ipx_mdest, rbuffm_raddr(ipx_mdest+1), wk_f_laddr+(icasm-1)*3*8, MSG_SEND*8)
       call rdma_wait(ipx_msrc)
    endif
#else
    MSG_RECV=3*ncap
    MSG_SEND=3*ncap
    call mpi_irecv(rbuffp, MSG_RECV, MPI_DOUBLE_PRECISION, &
         &            ipx_psrc,       1, mpi_comm_world, &
         &            irq(1), ierr)
    call mpi_isend(wk_f(1,icasp), MSG_SEND, MPI_DOUBLE_PRECISION, &
         &            ipx_pdest,      1, mpi_comm_world, &
         &            irq(2), ierr)
    MSG_RECV=3*ncam
    MSG_SEND=3*ncam
    call mpi_irecv(rbuffm, MSG_RECV, MPI_DOUBLE_PRECISION, &
         &            ipx_msrc,       2, mpi_comm_world, &
         &            irq(3), ierr)
    call mpi_isend(wk_f(1,icasm), MSG_SEND, MPI_DOUBLE_PRECISION, &
         &            ipx_mdest,      2, mpi_comm_world, &
         &            irq(4), ierr)

    nrq = 4
    call mpi_waitall(nrq, irq, istatus, ierr)
#endif
#ifdef PROFILE_COMM
    call end_profile(comm_direct2_dr_id)
#endif

    nca = 0
    DO ica = icarp, icarp + ncap - 1
       nca = nca + 1
       wk_f(1,ica)= wk_f(1,ica) + rbuffp(1,nca)
       wk_f(2,ica)= wk_f(2,ica) + rbuffp(2,nca)
       wk_f(3,ica)= wk_f(3,ica) + rbuffp(3,nca)
    END DO
    nca = 0
    DO ica = icarm, icarm + ncam - 1
       nca = nca + 1
       wk_f(1,ica)= wk_f(1,ica) + rbuffm(1,nca)
       wk_f(2,ica)= wk_f(2,ica) + rbuffm(2,nca)
       wk_f(3,ica)= wk_f(3,ica) + rbuffm(3,nca)
    END DO

    !     force +Y
    ipy_pdest  = ipz*npx*npy + mod(ipy+1-1/npy+npy,npy)*npx + ipx
    ipy_psrc   = ipz*npx*npy + mod(ipy-1+1/npy+npy,npy)*npx + ipx
    !     force -Y
    ipy_mdest  = ipz*npx*npy + mod(ipy-1+1/npy+npy,npy)*npx + ipx
    ipy_msrc   = ipz*npx*npy + mod(ipy+1-1/npy+npy,npy)*npx + ipx

    icz0 = 1 - nbd
    icz1 = lzdiv + nbd

#ifndef YCOMM_NOBUFFERING
    icyp = lydiv + 1
    icybp = icyp - lydiv
    icym = 1 - nbd
    icybm = icym + lydiv

    ncarp = 0
    DO icx = 1, lxdiv
          DO icz = icz0, icz1
             ncarp = ncarp + na_per_cell(icz,icybp,icx)
          END DO
    END DO

    ncap = 0
    DO icx = 1, lxdiv
          DO ica = tag(icz0, icyp, icx), tag(icz1, icyp, icx) &
             &                 + na_per_cell(icz1, icyp, icx) - 1
             ncap = ncap + 1
             buffp(1,ncap) = wk_f(1,ica)
             buffp(2,ncap) = wk_f(2,ica)
             buffp(3,ncap) = wk_f(3,ica)
          END DO
    END DO

    ncarm = 0
    DO icx = 1, lxdiv
          DO icz = icz0, icz1
             ncarm = ncarm + na_per_cell(icz,icybm,icx)
          END DO
    END DO

    ncam = 0
    DO icx = 1, lxdiv
          DO ica = tag(icz0, icym, icx), tag(icz1, icym, icx) &
             &                 + na_per_cell(icz1, icym, icx) - 1
             ncam = ncam + 1
             buffm(1,ncam) = wk_f(1,ica)
             buffm(2,ncam) = wk_f(2,ica)
             buffm(3,ncam) = wk_f(3,ica)
          END DO
    END DO

#ifdef PROFILE_COMM
    call mpi_barrier(mpi_comm_world, ierr)
    call start_profile(comm_direct2_dr_id)
#endif
#ifdef SYNC_COM
    MSG_SEND=3*ncap
    MSG_RECV=3*ncarp
    call mpi_sendrecv(buffp, MSG_SEND, MPI_DOUBLE_PRECISION, &
         &            ipy_pdest, myrank,&
         &            rbuffp, MSG_RECV, MPI_DOUBLE_PRECISION, &
         &            ipy_psrc, ipy_psrc, mpi_comm_world, &
         &            istatus, ierr )

    MSG_SEND=3*ncam
    MSG_RECV=3*ncarm
    call mpi_sendrecv(buffm, MSG_SEND, MPI_DOUBLE_PRECISION, &
         &            ipy_mdest, myrank,&
         &            rbuffm, MSG_RECV, MPI_DOUBLE_PRECISION, &
         &            ipy_msrc, ipy_msrc, mpi_comm_world, &
         &            istatus, ierr )
#elif FJ_RDMA
    call rdma_post(ipy_psrc)
    if (ipy_psrc .ne. ipy_msrc) then
       call rdma_post(ipy_msrc)
       call rdma_wait(ipy_mdest)
    endif
    call rdma_wait(ipy_pdest)
    !
    if (ipy_pdest .ne. ipy_mdest) then
       MSG_SEND=3*ncap
       call rdma_put_post(ipy_pdest, rbuffp_raddr(ipy_pdest+1), buffp_laddr, MSG_SEND*8)
       MSG_SEND=3*ncam
       call rdma_put_post(ipy_mdest, rbuffm_raddr(ipy_mdest+1), buffm_laddr, MSG_SEND*8)
       call rdma_wait(ipy_psrc)
       call rdma_wait(ipy_msrc)
    else
       MSG_SEND=3*ncap
       call rdma_put(ipy_pdest, rbuffp_raddr(ipy_pdest+1), buffp_laddr, MSG_SEND*8)
       MSG_SEND=3*ncam
       call rdma_put_post(ipy_mdest, rbuffm_raddr(ipy_mdest+1), buffm_laddr, MSG_SEND*8)
       call rdma_wait(ipy_msrc)
    endif
#else
    MSG_RECV=3*ncarp
    MSG_SEND=3*ncap
    call mpi_irecv(rbuffp, MSG_RECV, MPI_DOUBLE_PRECISION, &
         &            ipy_psrc,       1, mpi_comm_world, &
         &            irq(1), ierr)
    call mpi_isend(buffp, MSG_SEND, MPI_DOUBLE_PRECISION, &
         &            ipy_pdest,      1, mpi_comm_world, &
         &            irq(2), ierr)

    MSG_RECV=3*ncarm
    MSG_SEND=3*ncam
    call mpi_irecv(rbuffm, MSG_RECV, MPI_DOUBLE_PRECISION, &
         &            ipy_msrc,       2, mpi_comm_world, &
         &            irq(3), ierr)
    call mpi_isend(buffm, MSG_SEND, MPI_DOUBLE_PRECISION, &
         &            ipy_mdest,      2, mpi_comm_world, &
         &            irq(4), ierr)

    nrq = 4
    call mpi_waitall(nrq, irq, istatus, ierr)
#endif
#ifdef PROFILE_COMM
    call end_profile(comm_direct2_dr_id)
#endif

    nca = 0
    DO icx = 1, lxdiv
          DO ica = tag(icz0, icybp, icx), tag(icz1, icybp, icx) &
                 &                 + na_per_cell(icz1, icybp, icx) - 1
             nca = nca + 1
             wk_f(1,ica)= wk_f(1,ica) + rbuffp(1,nca)
             wk_f(2,ica)= wk_f(2,ica) + rbuffp(2,nca)
             wk_f(3,ica)= wk_f(3,ica) + rbuffp(3,nca)
          END DO
    END DO
    nca = 0
    DO icx = 1, lxdiv
          DO ica = tag(icz0, icybm, icx), tag(icz1, icybm, icx) &
                 &                 + na_per_cell(icz1, icybm, icx) - 1
             nca = nca + 1
             wk_f(1,ica)= wk_f(1,ica) + rbuffm(1,nca)
             wk_f(2,ica)= wk_f(2,ica) + rbuffm(2,nca)
             wk_f(3,ica)= wk_f(3,ica) + rbuffm(3,nca)
          END DO
    END DO

#else

    DO icx = 1, lxdiv
       icyp = lydiv + 1
       icybp = icyp - lydiv
       icym = 1 - nbd
       icybm = icym + lydiv
       ncap = naline
       icasp = tag(1,icyp,icx) - na1cell
       icarp = tag(1,icybp,icx) - na1cell
       ncam = naline
       icasm = tag(1,icym,icx) - na1cell
       icarm = tag(1,icybm,icx) - na1cell

#ifdef SYNC_COM
    MSG_SEND=3*ncap
    MSG_RECV=3*ncap
       call mpi_sendrecv(wk_f(1,icasp), MSG_SEND, MPI_DOUBLE_PRECISION, &
            &            ipy_pdest, myrank,&
            &            rbuffp, MSG_RECV, MPI_DOUBLE_PRECISION, &
            &            ipy_psrc, ipy_psrc, mpi_comm_world, &
            &            istatus, ierr )

    MSG_SEND=3*ncam
    MSG_RECV=3*ncam
       call mpi_sendrecv(wk_f(1,icasm), MSG_SEND, MPI_DOUBLE_PRECISION, &
            &            ipy_mdest, myrank,&
            &            rbuffm, MSG_RECV, MPI_DOUBLE_PRECISION, &
            &            ipy_msrc, ipy_msrc, mpi_comm_world, &
            &            istatus, ierr )
#else
    MSG_RECV=3*ncap
    MSG_SEND=3*ncap
       call mpi_irecv(rbuffp, MSG_RECV, MPI_DOUBLE_PRECISION, &
            &            ipy_psrc,    1, mpi_comm_world, &
            &            irq(1), ierr)
       call mpi_isend(wk_f(1,icasp), MSG_SEND, MPI_DOUBLE_PRECISION, &
            &            ipy_pdest,   1, mpi_comm_world, &
            &            irq(2), ierr)

    MSG_RECV=3*ncam
    MSG_SEND=3*ncam
       call mpi_irecv(rbuffm, MSG_RECV, MPI_DOUBLE_PRECISION, &
            &            ipy_msrc,    2, mpi_comm_world, &
            &            irq(3), ierr)
       call mpi_isend(wk_f(1,icasm), MSG_SEND, MPI_DOUBLE_PRECISION, &
            &            ipy_mdest,   2, mpi_comm_world, &
            &            irq(4), ierr)

       nrq = 4
       call mpi_waitall(nrq, irq, istatus, ierr)
#endif

       nca = 0
       DO ica = icarp, icarp + ncap - 1
          nca = nca + 1
          wk_f(1,ica)= wk_f(1,ica) + rbuffp(1,nca)
          wk_f(2,ica)= wk_f(2,ica) + rbuffp(2,nca)
          wk_f(3,ica)= wk_f(3,ica) + rbuffp(3,nca)
       END DO
       nca = 0
       DO ica = icarm, icarm + ncam - 1
          nca = nca + 1
          wk_f(1,ica)= wk_f(1,ica) + rbuffm(1,nca)
          wk_f(2,ica)= wk_f(2,ica) + rbuffm(2,nca)
          wk_f(3,ica)= wk_f(3,ica) + rbuffm(3,nca)
       END DO

    END DO
#endif

    !     force +Z
    ipz_pdest = mod(ipz+1-1/npz+npz,npz)*npx*npy + ipy*npx + ipx
    ipz_psrc  = mod(ipz-1+1/npz+npz,npz)*npx*npy + ipy*npx + ipx
    !     force -Z
    ipz_mdest = mod(ipz-1+1/npz+npz,npz)*npx*npy + ipy*npx + ipx
    ipz_msrc  = mod(ipz+1-1/npz+npz,npz)*npx*npy + ipy*npx + ipx

    iczp = lzdiv + 1
    iczm = 1 - nbd

    iczbp = iczp - lzdiv
    iczbm = iczm + lzdiv

    ! number of receiving atoms on x-y plane.
    ! tag and na_per_cell for extra process boundary cells
    ! should be kept up-to-date. this is guaranteed by executing comm_dir_3_dr
    ! in advance. force calculation in boundary is performaed based on extra cell
    ! atoms by carrying out comm_dir_3_dr.
    ncarp = 0
    DO icx = 1, lxdiv
       DO icy = 1, lydiv
             ncarp = ncarp + na_per_cell(iczbp,icy,icx)
       END DO
    END DO

    ncap = 0
    DO icx = 1, lxdiv
       DO icy = 1, lydiv
          DO ica = tag(iczp, icy, icx), tag(iczp, icy, icx) &
             &                 + na_per_cell(iczp, icy, icx)-1
             ncap = ncap + 1
             buffp(1,ncap) = wk_f(1,ica)
             buffp(2,ncap) = wk_f(2,ica)
             buffp(3,ncap) = wk_f(3,ica)
          END DO
       END DO
    END DO

    ncarm = 0
    DO icx = 1, lxdiv
       DO icy = 1, lydiv
             ncarm = ncarm + na_per_cell(iczbm,icy,icx)
       END DO
    END DO

    ncam = 0
    DO icx = 1, lxdiv
       DO icy = 1, lydiv
          DO ica = tag(iczm, icy, icx), tag(iczm, icy, icx) &
             &                     + na_per_cell(iczm, icy, icx)-1
             ncam = ncam + 1
             buffm(1,ncam) = wk_f(1,ica)
             buffm(2,ncam) = wk_f(2,ica)
             buffm(3,ncam) = wk_f(3,ica)
          END DO
       END DO
    END DO

#ifdef PROFILE_COMM
    call mpi_barrier(mpi_comm_world, ierr)
    call start_profile(comm_direct2_dr_id)
#endif
#ifdef SYNC_COM
    MSG_SEND=3*ncap
    MSG_RECV=3*ncarp
    call mpi_sendrecv(buffp, MSG_SEND, MPI_DOUBLE_PRECISION, &
         &            ipz_pdest, myrank, &
         &            rbuffp, MSG_RECV, MPI_DOUBLE_PRECISION, &
         &            ipz_psrc, ipz_psrc, mpi_comm_world, &
         &            istatus, ierr )

    MSG_SEND=3*ncam
    MSG_RECV=3*ncarm
    call mpi_sendrecv(buffm, MSG_SEND, MPI_DOUBLE_PRECISION, &
         &            ipz_mdest, myrank, &
         &            rbuffm, MSG_RECV, MPI_DOUBLE_PRECISION, &
         &            ipz_msrc, ipz_msrc, mpi_comm_world, &
         &            istatus, ierr )
#elif FJ_RDMA
    call rdma_post(ipz_psrc)
    if (ipz_psrc .ne. ipz_msrc) then
       call rdma_post(ipz_msrc)
       call rdma_wait(ipz_mdest)
    endif
    call rdma_wait(ipz_pdest)
    !
    if (ipz_pdest .ne. ipz_mdest) then
       MSG_SEND=3*ncap
       call rdma_put_post(ipz_pdest, rbuffp_raddr(ipz_pdest+1), buffp_laddr, MSG_SEND*8)
       MSG_SEND=3*ncam
       call rdma_put_post(ipz_mdest, rbuffm_raddr(ipz_mdest+1), buffm_laddr, MSG_SEND*8)
       call rdma_wait(ipz_psrc)
       call rdma_wait(ipz_msrc)
    else
       MSG_SEND=3*ncap
       call rdma_put(ipz_pdest, rbuffp_raddr(ipz_pdest+1), buffp_laddr, MSG_SEND*8)
       MSG_SEND=3*ncam
       call rdma_put_post(ipz_mdest, rbuffm_raddr(ipz_mdest+1), buffm_laddr, MSG_SEND*8)
       call rdma_wait(ipz_msrc)
    endif
#else
    MSG_RECV=3*ncarp
    MSG_SEND=3*ncap
    call mpi_irecv(rbuffp, MSG_RECV, MPI_DOUBLE_PRECISION, &
         &            ipz_psrc,       1, mpi_comm_world, &
         &            irq(1), ierr)
    call mpi_isend(buffp, MSG_SEND, MPI_DOUBLE_PRECISION, &
         &            ipz_pdest,      1, mpi_comm_world, &
         &            irq(2), ierr)

    MSG_RECV=3*ncarm
    MSG_SEND=3*ncam
    call mpi_irecv(rbuffm, MSG_RECV, MPI_DOUBLE_PRECISION, &
         &            ipz_msrc,       2, mpi_comm_world, &
         &            irq(3), ierr)
    call mpi_isend(buffm, MSG_SEND, MPI_DOUBLE_PRECISION, &
         &            ipz_mdest,      2, mpi_comm_world, &
         &            irq(4), ierr)

    nrq = 4
    call mpi_waitall(nrq, irq, istatus, ierr)
#endif
#ifdef PROFILE_COMM
    call end_profile(comm_direct2_dr_id)
#endif

    nca = 0
    DO icx = 1, lxdiv
       DO icy = 1, lydiv
          DO ica = tag(iczbp, icy, icx), &
             &     tag(iczbp, icy, icx) &
             &   + na_per_cell(iczbp, icy, icx)-1
             nca = nca + 1
             wk_f(1,ica) = wk_f(1,ica)+rbuffp(1,nca)
             wk_f(2,ica) = wk_f(2,ica)+rbuffp(2,nca)
             wk_f(3,ica) = wk_f(3,ica)+rbuffp(3,nca)
          END DO
       END DO
    END DO
    nca = 0
    DO icx = 1, lxdiv
       DO icy = 1, lydiv
          DO ica = tag(iczbm, icy, icx), &
             &     tag(iczbm, icy, icx) &
             &   + na_per_cell(iczbm, icy, icx)-1
             nca = nca + 1
             wk_f(1,ica) = wk_f(1,ica)+rbuffm(1,nca)
             wk_f(2,ica) = wk_f(2,ica)+rbuffm(2,nca)
             wk_f(3,ica) = wk_f(3,ica)+rbuffm(3,nca)
          END DO
       END DO
    END DO

  end subroutine comm_direct_2_dr

end module comm_direct2_dr
