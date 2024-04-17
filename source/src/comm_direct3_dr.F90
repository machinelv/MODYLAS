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
!!         of atom coordinates in the MTD method.
!<
!! Ref: Y.Andoh, S.Ichikawa, T Sakashita, N Yoshii, S Okazaki,
!!      J. Comput. Chem., 42, 1073-1087 (2021).
!----------------------------------------------------------------------
!>
!! \brief  Module which relate to MPI communications of atom coordinates.
!! \author Yoshimichi Andoh, Shin-ichi Ichikawa
!<
module comm_direct3_dr
  use dr_cntl, only : nbd, nbd2, nadj !, nomp
  use domain, only  : lxdiv, lydiv, lzdiv
  use subcell, only : tag, na_per_cell, m2i, i2m, nselfatm, ndatm
  use trajectory_mpi, only : naline, narea, na1cell, nadirect
  use mpi_3d_grid
  use mpi_tool
#ifdef PROFILE_COMM
  use profile_comm
#endif
  implicit none
  integer,allocatable,dimension(:) :: icbufp
  integer,allocatable,dimension(:) :: ircbufp
  integer,allocatable,dimension(:) :: icbufm
  integer,allocatable,dimension(:) :: ircbufm
  integer,allocatable,dimension(:) :: ibuffp
  integer,allocatable,dimension(:) :: irbuffp
  integer,allocatable,dimension(:) :: ibuffm
  integer,allocatable,dimension(:) :: irbuffm
  real(8),allocatable,dimension(:,:) :: buffp
  real(8),allocatable,dimension(:,:) :: rbuffp
  real(8),allocatable,dimension(:,:) :: buffm
  real(8),allocatable,dimension(:,:) :: rbuffm
  integer :: MSG_SEND,MSG_RECV
#ifdef FJ_RDMA
  integer(8) :: icbufp_laddr
  integer(8) :: icbufm_laddr
  integer(8) :: buffp_laddr
  integer(8) :: ibuffp_laddr
  integer(8) :: buffm_laddr
  integer(8) :: ibuffm_laddr
  integer(8) :: na_per_cell_laddr
  integer(8) :: m2i_laddr
  integer(8),pointer :: ircbufp_raddr(:)
  integer(8),pointer :: ircbufm_raddr(:)
  integer(8),pointer :: rbuffp_raddr(:)
  integer(8),pointer :: irbuffp_raddr(:)
  integer(8),pointer :: rbuffm_raddr(:)
  integer(8),pointer :: irbuffm_raddr(:)
  integer(8),pointer :: na_per_cell_raddr(:)
  integer(8),pointer :: m2i_raddr(:)
#endif
contains

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to initialize arrays used in MPI communications.
!! \author Yoshimichi Andoh, Shin-ichi Ichikawa
!<
  subroutine init_comm_direct_3_dr()
!----------------------------------------------------------------------
    implicit none

#ifndef YCOMM_NOBUFFERING
    integer(4) max_cell, max_cell2
#ifdef FJ_RDMA
    type(c_ptr) :: ircbufp_cptr
    type(c_ptr) :: ircbufm_cptr 
    type(c_ptr) :: rbuffp_cptr
    type(c_ptr) :: irbuffp_cptr
    type(c_ptr) :: rbuffm_cptr
    type(c_ptr) :: irbuffm_cptr
    type(c_ptr) :: na_per_cell_cptr
    type(c_ptr) :: m2i_cptr
#endif

    max_cell = lydiv*lxdiv
    max_cell2 = (lzdiv + nbd2)*lxdiv
    if(max_cell < max_cell2) max_cell = max_cell2

    allocate(icbufp  (max_cell*nbd))
    allocate(ircbufp (max_cell*nbd))
    allocate(buffp (3,max_cell*nbd*na1cell))
    allocate(rbuffp(3,max_cell*nbd*na1cell))
    allocate(ibuffp (max_cell*nbd*na1cell))
    allocate(irbuffp(max_cell*nbd*na1cell))
!
    allocate(icbufm (max_cell*nbd))
    allocate(ircbufm(max_cell*nbd))
    allocate(buffm (3,max_cell*nbd*na1cell))
    allocate(rbuffm(3,max_cell*nbd*na1cell))
    allocate(ibuffm (max_cell*nbd*na1cell))
    allocate(irbuffm(max_cell*nbd*na1cell))

#ifdef FJ_RDMA
    call rdma_register_addr(ircbufp,     (max_cell*nbd)*4)
    call rdma_register_addr(icbufp,      (max_cell*nbd)*4)
    call rdma_register_addr(ircbufm,     (max_cell*nbd)*4)
    call rdma_register_addr(icbufm,      (max_cell*nbd)*4)
    call rdma_register_addr(buffp,       (3*max_cell*nbd*na1cell)*8)
    call rdma_register_addr(rbuffp,      (3*max_cell*nbd*na1cell)*8)
    call rdma_register_addr(ibuffp,      (max_cell*nbd*na1cell)*4)
    call rdma_register_addr(irbuffp,     (max_cell*nbd*na1cell)*4)
    call rdma_register_addr(buffm,       (3*max_cell*nbd*na1cell)*8)
    call rdma_register_addr(rbuffm,      (3*max_cell*nbd*na1cell)*8)
    call rdma_register_addr(ibuffm,      (max_cell*nbd*na1cell)*4)
    call rdma_register_addr(irbuffm ,    (max_cell*nbd*na1cell)*4)
    call rdma_register_addr(na_per_cell, (lzdiv+nbd+1)*(lydiv+nbd+1)*(lxdiv+nbd+1)*4)
    call rdma_register_addr(m2i,         nadirect*4)
    !
    icbufp_laddr      = rdma_get_laddr(icbufp)
    icbufm_laddr      = rdma_get_laddr(icbufm)
    buffp_laddr       = rdma_get_laddr(buffp)
    ibuffp_laddr      = rdma_get_laddr(ibuffp)
    buffm_laddr       = rdma_get_laddr(buffm)
    ibuffm_laddr      = rdma_get_laddr(ibuffm)
    na_per_cell_laddr = rdma_get_laddr(na_per_cell)
    m2i_laddr         = rdma_get_laddr(m2i)
    !
    ircbufp_cptr     = rdma_get_raddr(ircbufp)
    ircbufm_cptr     = rdma_get_raddr(ircbufm)
    rbuffp_cptr      = rdma_get_raddr(rbuffp)
    irbuffp_cptr     = rdma_get_raddr(irbuffp)
    rbuffm_cptr      = rdma_get_raddr(rbuffm)
    irbuffm_cptr     = rdma_get_raddr(irbuffm)
    na_per_cell_cptr = rdma_get_raddr(na_per_cell)
    m2i_cptr         = rdma_get_raddr(m2i)
    call c_f_pointer(ircbufp_cptr,     fptr=ircbufp_raddr,     shape=[nprocs])
    call c_f_pointer(ircbufm_cptr,     fptr=ircbufm_raddr,     shape=[nprocs])
    call c_f_pointer(rbuffp_cptr,      fptr=rbuffp_raddr,      shape=[nprocs])
    call c_f_pointer(irbuffp_cptr,     fptr=irbuffp_raddr,     shape=[nprocs])
    call c_f_pointer(rbuffm_cptr,      fptr=rbuffm_raddr,      shape=[nprocs])
    call c_f_pointer(irbuffm_cptr,     fptr=irbuffm_raddr,     shape=[nprocs])
    call c_f_pointer(na_per_cell_cptr, fptr=na_per_cell_raddr, shape=[nprocs])
    call c_f_pointer(m2i_cptr,         fptr=m2i_raddr,         shape=[nprocs])
#endif

#else
    allocate(icbufp  (lydiv*lxdiv*nbd))
    allocate(ircbufp (lydiv*lxdiv*nbd))
    allocate(buffp (3,lydiv*lxdiv*nbd*na1cell))
    allocate(rbuffp(3,lydiv*lxdiv*nbd*na1cell))
    allocate(ibuffp (lydiv*lxdiv*nbd*na1cell))
    allocate(irbuffp(lydiv*lxdiv*nbd*na1cell))
!
    allocate(icbufm (lydiv*lxdiv*nbd))
    allocate(ircbufm(lydiv*lxdiv*nbd))
    allocate(buffm (3,lydiv*lxdiv*nbd*na1cell))
    allocate(rbuffm(3,lydiv*lxdiv*nbd*na1cell))
    allocate(ibuffm (lydiv*lxdiv*nbd*na1cell))
    allocate(irbuffm(lydiv*lxdiv*nbd*na1cell))
#endif
    icbufp(:)  =0
    ircbufp(:) =0
    ibuffp(:)  =0
    irbuffp(:) =0
    icbufm(:)  =0
    ircbufm(:) =0
    ibuffm(:)  =0
    irbuffm(:) =0
    buffp(:,:) =0d0
    rbuffp(:,:)=0d0
    buffm(:,:) =0d0
    rbuffm(:,:)=0d0
  end subroutine init_comm_direct_3_dr
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to perform MPI communications of atom coordinates.
!! \author Yoshimichi Andoh, Shin-ichi Ichikawa
!<
  subroutine comm_direct_3_dr(wkxyz)
!----------------------------------------------------------------------
    use md_condition, only : mdstep
#include "timing.h90"
    implicit none

    INCLUDE 'mpif.h'
    real(8) wkxyz(3,nadirect)
    integer ipz_pdest, ipy_pdest, ipx_pdest
    integer ipz_psrc, ipy_psrc, ipx_psrc
    integer ipz_mdest, ipy_mdest, ipx_mdest
    integer ipz_msrc, ipy_msrc, ipx_msrc
    integer icz, icy, icx
    integer icz0, icz1
    integer icy0, icy1
    integer iczp
    integer icyp
    integer icxp
    integer iczm
    integer icym
    integer icxm
    integer iczbp
    integer iczbm
    integer icybp
    integer icybm
    integer icxbp
    integer icxbm
    integer ncc
    integer nccp
    integer nccm
    integer ica, icag
    integer icasp, icarp
    integer icasm, icarm
    integer nca
    integer ncap, ncarp
    integer ncam, ncarm
    integer nbase, nbase2, nbase3
    integer ntmp
    integer istatus(mpi_status_size, 8), ierr
#ifndef SYNC_COM
    integer,dimension(8) :: irq
    integer nrq
#endif

    !----- common parameters for coordinate communication. -----
    !     ipx=mod(myrank,npx)
    !     ipy=mod((myrank-ipx)/npx,npy)
    !     ipz=mod((myrank-ipx-ipy*npx)/(npx*npy),npz)

    !     lzdiv = (ncellz - 1)/npz + 1
    !     lydiv = (ncelly - 1)/npy + 1
    !     lxdiv = (ncellx - 1)/npx + 1

    !     narea  = na1cell * (lzdiv + 4) * (lydiv + 4)
    !     naline = na1cell * (lzdiv + 4)

    !
    !-----  coordinate communication code starts here. ------
    !
    !     coordinate +Z
    ipz_pdest = mod(ipz+1-1/npz+npz,npz)*npx*npy + ipy*npx + ipx
    ipz_psrc  = mod(ipz-1+1/npz+npz,npz)*npx*npy + ipy*npx + ipx
    !     coordinate -Z
    ipz_mdest = mod(ipz-1+1/npz+npz,npz)*npx*npy + ipy*npx + ipx
    ipz_msrc  = mod(ipz+1-1/npz+npz,npz)*npx*npy + ipy*npx + ipx

    iczp = lzdiv - nbd + 1
    nccp = 0
    DO icx = 1, lxdiv
       DO icy = 1, lydiv
          nccp = nccp + 1
          icbufp(nccp) = na_per_cell( iczp, icy, icx )
       END DO
    END DO

    iczm = 1
    nccm = 0
    DO icx = 1, lxdiv
       DO icy = 1, lydiv
          nccm = nccm + 1
          icbufm(nccm) = na_per_cell( iczm, icy, icx )
       END DO
    END DO

#ifdef PROFILE_COMM
    call mpi_barrier(mpi_comm_world, ierr)
    call start_profile(comm_direct3_dr_id)
#endif
#ifdef SYNC_COM
    call mpi_sendrecv(icbufp, nccp, MPI_INTEGER, &
         &            ipz_pdest, myrank, &
         &            ircbufp, nccp, MPI_INTEGER, &
         &            ipz_psrc, ipz_psrc, mpi_comm_world, &
         &            istatus, ierr )
    call mpi_sendrecv(icbufm, nccm, MPI_INTEGER, &
         &            ipz_mdest, myrank, &
         &            ircbufm, nccm, MPI_INTEGER, &
         &            ipz_msrc, ipz_msrc, mpi_comm_world, &
         &            istatus, ierr )
#elif FJ_RDMA
    if (ipz_pdest .ne. ipz_mdest) then
       call rdma_put_post(ipz_pdest, ircbufp_raddr(ipz_pdest+1), icbufp_laddr, nccp*4)
       call rdma_put_post(ipz_mdest, ircbufm_raddr(ipz_mdest+1), icbufm_laddr, nccm*4)
       call rdma_wait(ipz_psrc)
       call rdma_wait(ipz_msrc)
    else
       call rdma_put(ipz_pdest, ircbufp_raddr(ipz_pdest+1), icbufp_laddr, nccp*4)
       call rdma_put_post(ipz_mdest, ircbufm_raddr(ipz_mdest+1), icbufm_laddr, nccm*4)
       call rdma_wait(ipz_msrc)
    endif
#else
    call mpi_irecv(ircbufp, nccp, MPI_INTEGER, &
         &            ipz_psrc,       1, mpi_comm_world, &
         &            irq(1), ierr)
    call mpi_isend(icbufp, nccp, MPI_INTEGER, &
         &            ipz_pdest,      1, mpi_comm_world, &
         &            irq(2), ierr)
    call mpi_irecv(ircbufm, nccm, MPI_INTEGER, &
         &            ipz_msrc,       2, mpi_comm_world, &
         &            irq(3), ierr)
    call mpi_isend(icbufm, nccm, MPI_INTEGER, &
         &            ipz_mdest,      2, mpi_comm_world, &
         &            irq(4), ierr)
    nrq = 4
    call mpi_waitall(nrq, irq, istatus, ierr)
#endif
#ifdef PROFILE_COMM
    call end_profile(comm_direct3_dr_id)
#endif

    iczbp = iczp - lzdiv
    iczbm = iczm + lzdiv

    ncc = 0
    ncarp = 0
    DO icx = 1, lxdiv
       DO icy = 1, lydiv
          nca = tag(iczbp+1,icy,icx) - ircbufp(ncc+1)
          ncc = ncc + 1
          na_per_cell(iczbp,icy,icx) = ircbufp(ncc)
          tag(iczbp,icy,icx) = nca
          nca = nca + na_per_cell(iczbp,icy,icx)
          ncarp = ncarp + na_per_cell(iczbp,icy,icx)
       END DO
    END DO

    ncap = 0
    DO icx = 1, lxdiv
       DO icy = 1, lydiv
          DO ica = tag(iczp, icy, icx), tag(iczp, icy, icx) &
               &                 + na_per_cell(iczp, icy, icx)-1
             ncap = ncap + 1
             buffp(1,ncap) = wkxyz(1,ica)
             buffp(2,ncap) = wkxyz(2,ica)
             buffp(3,ncap) = wkxyz(3,ica)
             ibuffp(ncap) = m2i(ica)
          END DO
       END DO
    END DO

    ncc = 0
    ncarm = 0
    DO icx = 1, lxdiv
       DO icy = 1, lydiv
          nca = tag(iczbm-1,icy,icx) + na_per_cell(iczbm-1,icy,icx)
          ncc = ncc + 1
          na_per_cell(iczbm,icy,icx) = ircbufm(ncc)
          tag(iczbm,icy,icx) = nca
          nca = nca + na_per_cell(iczbm,icy,icx)
          ncarm = ncarm + na_per_cell(iczbm,icy,icx)
       END DO
    END DO

    ncam = 0
    DO icx = 1, lxdiv
       DO icy = 1, lydiv
          DO ica = tag(iczm, icy, icx), tag(iczm, icy, icx) &
               &                     + na_per_cell(iczm, icy, icx)-1
             ncam = ncam + 1
             buffm(1,ncam) = wkxyz(1,ica)
             buffm(2,ncam) = wkxyz(2,ica)
             buffm(3,ncam) = wkxyz(3,ica)
             ibuffm(ncam) = m2i(ica)
          END DO
       END DO
    END DO

#ifdef PROFILE_COMM
    call mpi_barrier(mpi_comm_world, ierr)
    call start_profile(comm_direct3_dr_id)
#endif
#ifdef SYNC_COM
    MSG_SEND=3*ncap
    MSG_RECV=3*ncarp
    call mpi_sendrecv(buffp,MSG_SEND, MPI_DOUBLE_PRECISION, &
         &            ipz_pdest, myrank, &
         &            rbuffp, MSG_RECV, MPI_DOUBLE_PRECISION, &
         &            ipz_psrc, ipz_psrc, mpi_comm_world, &
         &            istatus, ierr )
    call mpi_sendrecv(ibuffp, ncap, MPI_INTEGER, &
         &            ipz_pdest, myrank, &
         &            irbuffp, ncarp, MPI_INTEGER, &
         &            ipz_psrc, ipz_psrc, mpi_comm_world, &
         &            istatus, ierr )
    MSG_SEND=3*ncam
    MSG_RECV=3*ncarm
    call mpi_sendrecv(buffm,MSG_SEND, MPI_DOUBLE_PRECISION, &
         &            ipz_mdest, myrank, &
         &            rbuffm, MSG_RECV, MPI_DOUBLE_PRECISION, &
         &            ipz_msrc, ipz_msrc, mpi_comm_world, &
         &            istatus, ierr )
    call mpi_sendrecv(ibuffm, ncam, MPI_INTEGER, &
         &            ipz_mdest, myrank,  &
         &            irbuffm, ncarm, MPI_INTEGER, &
         &            ipz_msrc, ipz_msrc, mpi_comm_world, &
         &            istatus, ierr )
#elif FJ_RDMA
    if (ipz_pdest .ne. ipz_mdest) then
       MSG_SEND=3*ncap
       call rdma_put(ipz_pdest, rbuffp_raddr(ipz_pdest+1),  buffp_laddr,  MSG_SEND*8)
       call rdma_put_post(ipz_pdest, irbuffp_raddr(ipz_pdest+1), ibuffp_laddr, ncap*4)
       MSG_SEND=3*ncam
       call rdma_put(ipz_mdest, rbuffm_raddr(ipz_mdest+1),  buffm_laddr,  MSG_SEND*8)
       call rdma_put_post(ipz_mdest, irbuffm_raddr(ipz_mdest+1), ibuffm_laddr, ncam*4)
       !
       call rdma_wait(ipz_psrc)
       call rdma_wait(ipz_msrc)
    else
       MSG_SEND=3*ncap
       call rdma_put(ipz_pdest, rbuffp_raddr(ipz_pdest+1),  buffp_laddr,  MSG_SEND*8)
       call rdma_put(ipz_pdest, irbuffp_raddr(ipz_pdest+1), ibuffp_laddr, ncap*4)
       MSG_SEND=3*ncam
       call rdma_put(ipz_mdest, rbuffm_raddr(ipz_mdest+1),  buffm_laddr,  MSG_SEND*8)
       call rdma_put_post(ipz_mdest, irbuffm_raddr(ipz_mdest+1), ibuffm_laddr, ncam*4)
       !
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
    call mpi_irecv(irbuffp, ncarp, MPI_INTEGER, &
         &            ipz_psrc,       2, mpi_comm_world, &
         &            irq(3), ierr)
    call mpi_isend(ibuffp, ncap, MPI_INTEGER, &
         &            ipz_pdest,      2, mpi_comm_world, &
         &            irq(4), ierr)

    MSG_RECV=3*ncarm
    MSG_SEND=3*ncam
    call mpi_irecv(rbuffm, MSG_RECV, MPI_DOUBLE_PRECISION, &
         &            ipz_msrc,       3, mpi_comm_world, &
         &            irq(5), ierr)
    call mpi_isend(buffm, MSG_SEND, MPI_DOUBLE_PRECISION, &
         &            ipz_mdest,      3, mpi_comm_world, &
         &            irq(6), ierr)
    call mpi_irecv(irbuffm, ncarm, MPI_INTEGER, &
         &            ipz_msrc,       4, mpi_comm_world, &
         &            irq(7), ierr)
    call mpi_isend(ibuffm, ncam, MPI_INTEGER, &
         &            ipz_mdest,      4, mpi_comm_world, &
         &            irq(8), ierr)

    nrq = 8
    call mpi_waitall(nrq, irq, istatus, ierr)
#endif
#ifdef PROFILE_COMM
    call end_profile(comm_direct3_dr_id)
#endif

    nca = 0
    DO icx = 1, lxdiv
       DO icy = 1, lydiv
          DO ica = tag(iczbp, icy, icx), tag(iczbp, icy, icx) &
               &                      + na_per_cell(iczbp, icy, icx)-1
             nca = nca + 1
             wkxyz(1,ica) = rbuffp(1,nca)
             wkxyz(2,ica) = rbuffp(2,nca)
             wkxyz(3,ica) = rbuffp(3,nca)
             m2i(ica) = irbuffp(nca)
          END DO
       END DO
    END DO

    nca = 0
    DO icx = 1, lxdiv
       DO icy = 1, lydiv
          DO ica = tag(iczbm, icy, icx), tag(iczbm, icy, icx) &
               &                 + na_per_cell(iczbm, icy, icx)-1
             nca = nca + 1
             wkxyz(1,ica) = rbuffm(1,nca)
             wkxyz(2,ica) = rbuffm(2,nca)
             wkxyz(3,ica) = rbuffm(3,nca)
             m2i(ica) = irbuffm(nca)
          END DO
       END DO
    END DO

    !     coordinate +Y
    ipy_pdest  = ipz*npx*npy + mod(ipy+1-1/npy+npy,npy)*npx + ipx
    ipy_psrc   = ipz*npx*npy + mod(ipy-1+1/npy+npy,npy)*npx + ipx
    !     coordinate -Y
    ipy_mdest  = ipz*npx*npy + mod(ipy-1+1/npy+npy,npy)*npx + ipx
    ipy_msrc   = ipz*npx*npy + mod(ipy+1-1/npy+npy,npy)*npx + ipx

    icz0 = 1 - nbd
    icz1 = lzdiv + nbd

#ifndef YCOMM_NOBUFFERING
    icyp = lydiv - nbd + 1
    icym = 1
    nccp = 0
    DO icx = 1, lxdiv
       DO icz = icz0, icz1
          nccp = nccp + 1
          icbufp(nccp) = na_per_cell( icz, icyp, icx )
       END DO
    END DO
    nccm = 0
    DO icx = 1, lxdiv
       DO icz = icz0, icz1
          nccm = nccm + 1
          icbufm(nccm) = na_per_cell( icz, icym, icx )
       END DO
    END DO

#ifdef PROFILE_COMM
    call mpi_barrier(mpi_comm_world, ierr)
    call start_profile(comm_direct3_dr_id)
#endif
#ifdef SYNC_COM
    call mpi_sendrecv(icbufp, nccp, MPI_INTEGER, &
         &            ipy_pdest, myrank, &
         &            ircbufp, nccp, MPI_INTEGER, &
         &            ipy_psrc, ipy_psrc, mpi_comm_world, &
         &            istatus, ierr )
    call mpi_sendrecv(icbufm, nccm, MPI_INTEGER, &
         &            ipy_mdest, myrank, &
         &            ircbufm, nccm, MPI_INTEGER, &
         &            ipy_msrc, ipy_msrc, mpi_comm_world, &
         &            istatus, ierr )
#elif FJ_RDMA
    call rdma_post(ipy_psrc)
    if (ipy_psrc .ne. ipy_msrc) then
       call rdma_post(ipy_msrc)
       call rdma_wait(ipy_mdest)
    endif
    call rdma_wait(ipy_pdest)
    
    if (ipy_pdest .ne. ipy_mdest) then
       call rdma_put_post(ipy_pdest, ircbufp_raddr(ipy_pdest+1), icbufp_laddr, nccp*4)
       call rdma_put_post(ipy_mdest, ircbufm_raddr(ipy_mdest+1), icbufm_laddr, nccm*4)
       call rdma_wait(ipy_psrc)
       call rdma_wait(ipy_msrc)
    else
       call rdma_put(ipy_pdest, ircbufp_raddr(ipy_pdest+1), icbufp_laddr, nccp*4)
       call rdma_put_post(ipy_mdest, ircbufm_raddr(ipy_mdest+1), icbufm_laddr, nccm*4)
       call rdma_wait(ipy_msrc)
    end if
#else
    call mpi_irecv(ircbufp, nccp, MPI_INTEGER, &
         &            ipy_psrc,       1, mpi_comm_world, &
         &            irq(1), ierr)
    call mpi_isend(icbufp, nccp, MPI_INTEGER, &
         &            ipy_pdest,      1, mpi_comm_world, &
         &            irq(2), ierr)
    call mpi_irecv(ircbufm, nccm, MPI_INTEGER, &
         &            ipy_msrc,       2, mpi_comm_world, &
         &            irq(3), ierr)
    call mpi_isend(icbufm, nccm, MPI_INTEGER, &
         &            ipy_mdest,      2, mpi_comm_world, &
         &            irq(4), ierr)
    nrq = 4
    call mpi_waitall(nrq, irq, istatus, ierr)
#endif
#ifdef PROFILE_COMM
    call end_profile(comm_direct3_dr_id)
#endif
    icybp = icyp - lydiv
    icybm = icym + lydiv

    ncc = 0
    ncarp = 0
    DO icx = 1, lxdiv
       nbase = tag(1, 1, icx) - naline
       nca = nbase - ircbufp(ncc+1)
       DO icz = icz0, icz1
          ncc = ncc + 1
          na_per_cell(icz,icybp,icx) = ircbufp(ncc)
          tag(icz,icybp,icx) = nca
          nca = nca + ircbufp(ncc)
          ncarp = ncarp + ircbufp(ncc)
       END DO
    END DO

    ncc = 0
    ncarm = 0
    DO icx = 1, lxdiv
       nbase = tag(1,lydiv,icx) + naline
       nca = nbase - ircbufm(ncc+1)
       DO icz = icz0, icz1
          ncc = ncc + 1
          tag(icz,icybm,icx) = nca
          na_per_cell(icz,icybm,icx) = ircbufm(ncc)
          nca = nca + ircbufm(ncc)
          ncarm = ncarm + ircbufm(ncc)
       END DO
    END DO

    ncap = 0
    DO icx = 1, lxdiv
       DO ica = tag(icz0, icyp, icx), tag(icz1, icyp, icx) &
            &                 + na_per_cell(icz1, icyp, icx)-1
          ncap = ncap + 1
          buffp(1,ncap) = wkxyz(1,ica)
          buffp(2,ncap) = wkxyz(2,ica)
          buffp(3,ncap) = wkxyz(3,ica)
          ibuffp(ncap) = m2i(ica)
       END DO
    END DO
    ncam = 0
    DO icx = 1, lxdiv
       DO ica = tag(icz0, icym, icx), tag(icz1, icym, icx) &
            &                 + na_per_cell(icz1, icym, icx)-1
          ncam = ncam + 1
          buffm(1,ncam) = wkxyz(1,ica)
          buffm(2,ncam) = wkxyz(2,ica)
          buffm(3,ncam) = wkxyz(3,ica)
          ibuffm(ncam) = m2i(ica)
       END DO
    END DO

#ifdef PROFILE_COMM
    call mpi_barrier(mpi_comm_world, ierr)
    call start_profile(comm_direct3_dr_id)
#endif
#ifdef SYNC_COM
    MSG_SEND=3*ncap
    MSG_RECV=3*ncarp
    call mpi_sendrecv(buffp,MSG_SEND, MPI_DOUBLE_PRECISION, &
         &            ipy_pdest, myrank, &
         &            rbuffp, MSG_RECV, MPI_DOUBLE_PRECISION, &
         &            ipy_psrc, ipy_psrc, mpi_comm_world, &
         &            istatus, ierr )
    call mpi_sendrecv(ibuffp, ncap, MPI_INTEGER, &
         &            ipy_pdest, myrank, &
         &            irbuffp, ncarp, MPI_INTEGER, &
         &            ipy_psrc, ipy_psrc, mpi_comm_world, &
         &            istatus, ierr )
    MSG_SEND=3*ncam
    MSG_RECV=3*ncarm
    call mpi_sendrecv(buffm, MSG_SEND, MPI_DOUBLE_PRECISION, &
         &            ipy_mdest, myrank, &
         &            rbuffm, MSG_RECV, MPI_DOUBLE_PRECISION, &
         &            ipy_msrc, ipy_msrc, mpi_comm_world, &
         &            istatus, ierr )
    call mpi_sendrecv(ibuffm, ncam, MPI_INTEGER, &
         &            ipy_mdest, myrank, &
         &            irbuffm, ncarm, MPI_INTEGER, &
         &            ipy_msrc, ipy_msrc, mpi_comm_world, &
         &            istatus, ierr )
#elif FJ_RDMA
    if (ipy_pdest .ne. ipy_mdest) then
       MSG_SEND=3*ncap
       call rdma_put(ipy_pdest, rbuffp_raddr(ipy_pdest+1),  buffp_laddr,  MSG_SEND*8)
       call rdma_put_post(ipy_pdest, irbuffp_raddr(ipy_pdest+1), ibuffp_laddr, ncap*4)
       MSG_SEND=3*ncam
       call rdma_put(ipy_mdest, rbuffm_raddr(ipy_mdest+1),  buffm_laddr,  MSG_SEND*8)
       call rdma_put_post(ipy_mdest, irbuffm_raddr(ipy_mdest+1), ibuffm_laddr, ncam*4)
       !
       call rdma_wait(ipy_psrc)
       call rdma_wait(ipy_msrc)
    else
       MSG_SEND=3*ncap
       call rdma_put(ipy_pdest, rbuffp_raddr(ipy_pdest+1),  buffp_laddr,  MSG_SEND*8)
       call rdma_put(ipy_pdest, irbuffp_raddr(ipy_pdest+1), ibuffp_laddr, ncap*4)
       MSG_SEND=3*ncam
       call rdma_put(ipy_mdest, rbuffm_raddr(ipy_mdest+1),  buffm_laddr,  MSG_SEND*8)
       call rdma_put_post(ipy_mdest, irbuffm_raddr(ipy_mdest+1), ibuffm_laddr, ncam*4)
       !
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
    call mpi_irecv(irbuffp, ncarp, MPI_INTEGER, &
         &            ipy_psrc,       2, mpi_comm_world, &
         &            irq(3), ierr)
    call mpi_isend(ibuffp, ncap, MPI_INTEGER, &
         &            ipy_pdest,      2, mpi_comm_world, &
         &            irq(4), ierr)

    MSG_RECV=3*ncarm
    MSG_SEND=3*ncam
    call mpi_irecv(rbuffm, MSG_RECV, MPI_DOUBLE_PRECISION, &
         &            ipy_msrc,       3, mpi_comm_world, &
         &            irq(5), ierr)
    call mpi_isend(buffm, MSG_SEND, MPI_DOUBLE_PRECISION, &
         &            ipy_mdest,      3, mpi_comm_world, &
         &            irq(6), ierr)
    call mpi_irecv(irbuffm, ncarm, MPI_INTEGER, &
         &            ipy_msrc,       4, mpi_comm_world, &
         &            irq(7), ierr)
    call mpi_isend(ibuffm, ncam, MPI_INTEGER, &
         &            ipy_mdest,      4, mpi_comm_world, &
         &            irq(8), ierr)
    nrq = 8
    call mpi_waitall(nrq, irq, istatus, ierr)
#endif
#ifdef PROFILE_COMM
    call end_profile(comm_direct3_dr_id)
#endif

    nca = 0
    DO icx = 1, lxdiv
       DO ica = tag(icz0,icybp,icx), tag(icz1,icybp,icx) &
            &                   + na_per_cell(icz1,icybp,icx) - 1
          nca = nca + 1
          wkxyz(1,ica) = rbuffp(1,nca)
          wkxyz(2,ica) = rbuffp(2,nca)
          wkxyz(3,ica) = rbuffp(3,nca)
          m2i(ica) = irbuffp(nca)
       END DO
    END DO

    nca = 0
    DO icx = 1, lxdiv
       DO ica = tag(icz0,icybm,icx), tag(icz1,icybm,icx) &
            &                   + na_per_cell(icz1,icybm,icx) - 1
          nca = nca + 1
          wkxyz(1,ica) = rbuffm(1,nca)
          wkxyz(2,ica) = rbuffm(2,nca)
          wkxyz(3,ica) = rbuffm(3,nca)
          m2i(ica) = irbuffm(nca)
       END DO
    END DO

#else

    DO icx = 1, lxdiv
       icyp = lydiv - nbd + 1
       nccp = (icz1 - icz0 + 1)
       icybp = icyp - lydiv
       icym = 1
       nccm = (icz1 - icz0 + 1)
       icybm = icym + lydiv
#ifdef PROFILE_COMM
    call mpi_barrier(mpi_comm_world, ierr)
    call start_profile(comm_direct3_dr_id)
#endif
#ifdef SYNC_COM
       call mpi_sendrecv(na_per_cell(icz0,icyp,icx), nccp, MPI_INTEGER, &
            &            ipy_pdest, myrank, &
            &            na_per_cell(icz0,icybp,icx), nccp, MPI_INTEGER, &
            &            ipy_psrc, ipy_psrc, mpi_comm_world, &
            &            istatus, ierr )
       call mpi_sendrecv(na_per_cell(icz0,icym,icx), nccm, MPI_INTEGER, &
            &            ipy_mdest, myrank, &
            &            na_per_cell(icz0,icybm,icx), nccm, MPI_INTEGER, &
            &            ipy_msrc, ipy_msrc, mpi_comm_world, &
            &            istatus, ierr )
#else
       call mpi_irecv(na_per_cell(icz0,icybp,icx), nccp, MPI_INTEGER, &
            &            ipy_psrc,    1, mpi_comm_world, &
            &            irq(1), ierr)
       call mpi_isend(na_per_cell(icz0,icyp,icx), nccp, MPI_INTEGER, &
            &            ipy_pdest,   1, mpi_comm_world, &
            &            irq(2), ierr)
       call mpi_irecv(na_per_cell(icz0,icybm,icx), nccm, MPI_INTEGER, &
            &            ipy_msrc,    2, mpi_comm_world, &
            &            irq(3), ierr)
       call mpi_isend(na_per_cell(icz0,icym,icx), nccm, MPI_INTEGER, &
            &            ipy_mdest,   2, mpi_comm_world, &
            &            irq(4), ierr)
       nrq = 4
       call mpi_waitall(nrq, irq, istatus, ierr)
#endif
#ifdef PROFILE_COMM
    call end_profile(comm_direct3_dr_id)
#endif

       nbase = tag(1, 1, icx) - nbd*naline
       nca = nbase - na_per_cell(icz0,icybp,icx)
       DO icz = icz0, icz1
          tag(icz,icybp,icx) = nca
          nca = nca + na_per_cell(icz,icybp,icx)
       END DO

       nbase = tag(1,lydiv,icx) + naline
       nca = nbase - na_per_cell(icz0,icybm,icx)
       DO icz = icz0, icz1
          tag(icz,icybm,icx) = nca
          nca = nca + na_per_cell(icz,icybm,icx)
       END DO

       ncap = naline
       icasp = tag(1,icyp,icx) - nbd*na1cell
       icarp = tag(1,icybp,icx) - nbd*na1cell
       ncam = naline
       icasm = tag(1,icym,icx) - nbd*na1cell
       icarm = tag(1,icybm,icx) - nbd*na1cell

#ifdef PROFILE_COMM
    call mpi_barrier(mpi_comm_world, ierr)
    call start_profile(comm_direct3_dr_id)
#endif
#ifdef SYNC_COM
    MSG_SEND=3*ncap
    MSG_RECV=3*ncap
       call mpi_sendrecv(wkxyz(1,icasp),MSG_SEND, MPI_DOUBLE_PRECISION, &
            &            ipy_pdest, myrank, &
            &            wkxyz(1,icarp), MSG_RECV, MPI_DOUBLE_PRECISION, &
            &            ipy_psrc, ipy_psrc, mpi_comm_world, &
            &            istatus, ierr )
       call mpi_sendrecv(m2i(icasp), ncap, MPI_INTEGER, &
            &            ipy_pdest, myrank, &
            &            m2i(icarp), ncap, MPI_INTEGER, &
            &            ipy_psrc, ipy_psrc, mpi_comm_world, &
            &            istatus, ierr )
    MSG_SEND=3*ncam
    MSG_RECV=3*ncam
       call mpi_sendrecv(wkxyz(1,icasm), MSG_SEND, MPI_DOUBLE_PRECISION, &
            &            ipy_mdest, myrank, &
            &            wkxyz(1,icarm), MSG_RECV, MPI_DOUBLE_PRECISION, &
            &            ipy_msrc, ipy_msrc, mpi_comm_world, &
            &            istatus, ierr )
       call mpi_sendrecv(m2i(icasm), ncam, MPI_INTEGER, &
            &            ipy_mdest, myrank, &
            &            m2i(icarm), ncam, MPI_INTEGER, &
            &            ipy_msrc, ipy_msrc, mpi_comm_world, &
            &            istatus, ierr )
#else
    MSG_RECV=3*ncap
    MSG_SEND=3*ncap
       call mpi_irecv(wkxyz(1,icarp), MSG_RECV, MPI_DOUBLE_PRECISION, &
            &            ipy_psrc,    1, mpi_comm_world, &
            &            irq(1), ierr)
       call mpi_isend(wkxyz(1,icasp), MSG_SEND, MPI_DOUBLE_PRECISION, &
            &            ipy_pdest,   1, mpi_comm_world, &
            &            irq(2), ierr)
       call mpi_irecv(m2i(icarp), ncap, MPI_INTEGER, &
            &            ipy_psrc,    2, mpi_comm_world, &
            &            irq(3), ierr)
       call mpi_isend(m2i(icasp), ncap, MPI_INTEGER, &
            &            ipy_pdest,   2, mpi_comm_world, &
            &            irq(4), ierr)

    MSG_RECV=3*ncap
    MSG_SEND=3*ncap
       call mpi_irecv(wkxyz(1,icarm), MSG_RECV, MPI_DOUBLE_PRECISION, &
            &            ipy_msrc,    3, mpi_comm_world, &
            &            irq(5), ierr)
       call mpi_isend(wkxyz(1,icasm), MSG_SEND, MPI_DOUBLE_PRECISION, &
            &            ipy_mdest,   3, mpi_comm_world, &
            &            irq(6), ierr)
       call mpi_irecv(m2i(icarm), ncam, MPI_INTEGER, &
            &            ipy_msrc,    4, mpi_comm_world, &
            &            irq(7), ierr)
       call mpi_isend(m2i(icasm), ncam, MPI_INTEGER, &
            &            ipy_mdest,   4, mpi_comm_world, &
            &            irq(8), ierr)
       nrq = 8
       call mpi_waitall(nrq, irq, istatus, ierr)
#endif
#ifdef PROFILE_COMM
    call end_profile(comm_direct3_dr_id)
#endif

    END DO
#endif
    
    !     coordinate +X
    ipx_pdest  = ipz*npx*npy + ipy*npx + mod(ipx+1-1/npx+npx,npx)
    ipx_psrc   = ipz*npx*npy + ipy*npx + mod(ipx-1+1/npx+npx,npx)
    !     coordinate -X
    ipx_mdest  = ipz*npx*npy + ipy*npx + mod(ipx-1+1/npx+npx,npx)
    ipx_msrc   = ipz*npx*npy + ipy*npx + mod(ipx+1-1/npx+npx,npx)
    icz0 = 1 - nbd
    icz1 = lzdiv + nbd
    icy0 = 1 - nbd
    icy1 = lydiv + nbd

    icxp = lxdiv - nbd + 1
    nccp = (icz1 - icz0 +1)*(icy1 - icy0 +1)
    icxbp = icxp - lxdiv

    icxm = 1
    nccm = (icz1 - icz0 +1)*(icy1 - icy0 +1)
    icxbm = icxm + lxdiv

#ifdef PROFILE_COMM
    call mpi_barrier(mpi_comm_world, ierr)
    call start_profile(comm_direct3_dr_id)
#endif
#ifdef SYNC_COM
    call mpi_sendrecv(na_per_cell(icz0,icy0,icxp), nccp, MPI_INTEGER, &
         &            ipx_pdest, myrank, &
         &            na_per_cell(icz0,icy0,icxbp), nccp, MPI_INTEGER,&
         &            ipx_psrc, ipx_psrc, mpi_comm_world, &
         &            istatus, ierr )
    call mpi_sendrecv(na_per_cell(icz0,icy0,icxm), nccm, MPI_INTEGER, &
         &            ipx_mdest, myrank, &
         &            na_per_cell(icz0,icy0,icxbm), nccm, MPI_INTEGER,&
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
       call rdma_put_post(ipx_pdest, &
            na_per_cell_raddr(ipx_pdest+1) + (icz0 + icy0*(lzdiv+nbd+1) + icxbp*((lzdiv+nbd+1)*(lydiv+nbd+1)))*4, &
            na_per_cell_laddr              + (icz0 + icy0*(lzdiv+nbd+1) + icxp *((lzdiv+nbd+1)*(lydiv+nbd+1)))*4, &
            nccp*4)
       call rdma_put_post(ipx_mdest, &
            na_per_cell_raddr(ipx_mdest+1) + (icz0 + icy0*(lzdiv+nbd+1) + icxbm*((lzdiv+nbd+1)*(lydiv+nbd+1)))*4, &
            na_per_cell_laddr              + (icz0 + icy0*(lzdiv+nbd+1) + icxm *((lzdiv+nbd+1)*(lydiv+nbd+1)))*4,  &
            nccm*4)
       !
       call rdma_wait(ipx_psrc)
       call rdma_wait(ipx_msrc)
    else
       call rdma_put(ipx_pdest, &
            na_per_cell_raddr(ipx_pdest+1) + (icz0 + icy0*(lzdiv+nbd+1) + icxbp*((lzdiv+nbd+1)*(lydiv+nbd+1)))*4, &
            na_per_cell_laddr              + (icz0 + icy0*(lzdiv+nbd+1) + icxp *((lzdiv+nbd+1)*(lydiv+nbd+1)))*4, &
            nccp*4)
       call rdma_put_post(ipx_mdest, &
            na_per_cell_raddr(ipx_mdest+1) + (icz0 + icy0*(lzdiv+nbd+1) + icxbm*((lzdiv+nbd+1)*(lydiv+nbd+1)))*4, &
            na_per_cell_laddr              + (icz0 + icy0*(lzdiv+nbd+1) + icxm *((lzdiv+nbd+1)*(lydiv+nbd+1)))*4,  &
            nccm*4)
       !
       call rdma_wait(ipx_msrc)
    endif
#else
    call mpi_irecv(na_per_cell(icz0,icy0,icxbp), nccp, MPI_INTEGER, &
         &            ipx_psrc,       1, mpi_comm_world, &
         &            irq(1), ierr)
    call mpi_isend(na_per_cell(icz0,icy0,icxp), nccp, MPI_INTEGER, &
         &            ipx_pdest,      1, mpi_comm_world, &
         &            irq(2), ierr)
    call mpi_irecv(na_per_cell(icz0,icy0,icxbm), nccm, MPI_INTEGER, &
         &            ipx_msrc,       2, mpi_comm_world, &
         &            irq(3), ierr)
    call mpi_isend(na_per_cell(icz0,icy0,icxm), nccm, MPI_INTEGER, &
         &            ipx_mdest,      2, mpi_comm_world, &
         &            irq(4), ierr)
    nrq = 4
    call mpi_waitall(nrq, irq, istatus, ierr)
#endif
#ifdef PROFILE_COMM
    call end_profile(comm_direct3_dr_id)
#endif

    nbase3 = tag(1,1,1) - nbd*narea - (nbd+1)*naline
    nbase2 = nbase3
    DO icy = icy0, icy1
       nbase2 = nbase2 + naline
       nca = nbase2 - na_per_cell(icz0,icy,icxbp)
       DO icz = icz0, icz1
          tag(icz,icy,icxbp) = nca
          nca = nca + na_per_cell(icz,icy,icxbp)
       END DO
    END DO

    nbase3 = tag(1,1,lxdiv) + narea - (nbd+1)*naline
    nbase2 = nbase3
    DO icy = icy0, icy1
       nbase2 = nbase2 + naline
       nca = nbase2 - na_per_cell(icz0,icy,icxbm)
       DO icz = icz0, icz1
          tag(icz,icy,icxbm) = nca
          nca = nca + na_per_cell(icz,icy,icxbm)
       END DO
    END DO

    ncap = narea
    icasp = tag(1,icy0,icxp) - nbd*na1cell
    icarp = tag(1,icy0,icxbp) - nbd*na1cell
    ncam = narea
    icasm = tag(1,icy0,icxm) - nbd*na1cell
    icarm = tag(1,icy0,icxbm) - nbd*na1cell

#ifdef PROFILE_COMM
    call mpi_barrier(mpi_comm_world, ierr)
    call start_profile(comm_direct3_dr_id)
#endif
#ifdef SYNC_COM
    MSG_SEND=3*ncap
    MSG_RECV=3*ncap
    call mpi_sendrecv(wkxyz(1,icasp), MSG_SEND, MPI_DOUBLE_PRECISION, &
         &            ipx_pdest, myrank, &
         &            wkxyz(1,icarp), MSG_RECV, MPI_DOUBLE_PRECISION, &
         &            ipx_psrc, ipx_psrc, mpi_comm_world, &
         &            istatus, ierr )
    call mpi_sendrecv(m2i(icasp), ncap, MPI_INTEGER, &
         &            ipx_pdest, myrank, &
         &            m2i(icarp), ncap, MPI_INTEGER, &
         &            ipx_psrc, ipx_psrc, mpi_comm_world, &
         &            istatus, ierr )
    MSG_SEND=3*ncam
    MSG_RECV=3*ncam
    call mpi_sendrecv(wkxyz(1,icasm), MSG_SEND, MPI_DOUBLE_PRECISION, &
         &            ipx_mdest, myrank, &
         &            wkxyz(1,icarm), MSG_RECV, MPI_DOUBLE_PRECISION, &
         &            ipx_msrc, ipx_msrc, mpi_comm_world, &
         &            istatus, ierr )
    call mpi_sendrecv(m2i(icasm), ncam, MPI_INTEGER, &
         &            ipx_mdest, myrank, &
         &            m2i(icarm), ncam, MPI_INTEGER, &
         &            ipx_msrc, ipx_msrc, mpi_comm_world, &
         &            istatus, ierr )
#elif FJ_RDMA
    if (ipx_pdest .ne. ipx_mdest) then
       MSG_SEND=3*ncap
       call rdma_put(ipx_pdest, wkxyz_raddr(ipx_pdest+1) + (3*(icarp-1))*8, &
            wkxyz_laddr + (3*(icasp-1))*8, MSG_SEND*8)
       call rdma_put_post(ipx_pdest, m2i_raddr(ipx_pdest+1) + (icarp-1)*4, &
            m2i_laddr + (icasp-1)*4, ncap*4)
       MSG_SEND=3*ncam
       call rdma_put(ipx_mdest, wkxyz_raddr(ipx_mdest+1) + (3*(icarm-1))*8, &
            wkxyz_laddr + (3*(icasm-1))*8, MSG_SEND*8)
       call rdma_put_post(ipx_mdest, m2i_raddr(ipx_mdest+1) + (icarm-1)*4, &
            m2i_laddr + (icasm-1)*4, ncam*4)
       !
       call rdma_wait(ipx_psrc)
       call rdma_wait(ipx_msrc)
    else
       MSG_SEND=3*ncap
       call rdma_put(ipx_pdest, wkxyz_raddr(ipx_pdest+1) + (3*(icarp-1))*8, &
            wkxyz_laddr + (3*(icasp-1))*8, MSG_SEND*8)
       call rdma_put(ipx_pdest, m2i_raddr(ipx_pdest+1) + (icarp-1)*4, &
            m2i_laddr + (icasp-1)*4, ncap*4)
       MSG_SEND=3*ncam
       call rdma_put(ipx_mdest, wkxyz_raddr(ipx_mdest+1) + (3*(icarm-1))*8, &
            wkxyz_laddr + (3*(icasm-1))*8, MSG_SEND*8)
       call rdma_put_post(ipx_mdest, m2i_raddr(ipx_mdest+1) + (icarm-1)*4, &
            m2i_laddr + (icasm-1)*4, ncam*4)
       !
       call rdma_wait(ipx_msrc)
    endif
#else
    MSG_RECV=3*ncap
    MSG_SEND=3*ncap
    call mpi_irecv(wkxyz(1,icarp), MSG_RECV, MPI_DOUBLE_PRECISION, &
         &            ipx_psrc,       1, mpi_comm_world, &
         &            irq(1), ierr)
    call mpi_isend(wkxyz(1,icasp), MSG_SEND, MPI_DOUBLE_PRECISION, &
         &            ipx_pdest,      1, mpi_comm_world, &
         &            irq(2), ierr)
    call mpi_irecv(m2i(icarp), ncap, MPI_INTEGER, &
         &            ipx_psrc,       2, mpi_comm_world, &
         &            irq(3), ierr)
    call mpi_isend(m2i(icasp), ncap, MPI_INTEGER, &
         &            ipx_pdest,      2, mpi_comm_world, &
         &            irq(4), ierr)

    MSG_RECV=3*ncam
    MSG_SEND=3*ncam
    call mpi_irecv(wkxyz(1,icarm), MSG_RECV, MPI_DOUBLE_PRECISION, &
         &            ipx_msrc,       3, mpi_comm_world, &
         &            irq(5), ierr)
    call mpi_isend(wkxyz(1,icasm), MSG_SEND, MPI_DOUBLE_PRECISION, &
         &            ipx_mdest,      3, mpi_comm_world, &
         &            irq(6), ierr)
    call mpi_irecv(m2i(icarm), ncam, MPI_INTEGER, &
         &            ipx_msrc,       4, mpi_comm_world, &
         &            irq(7), ierr)
    call mpi_isend(m2i(icasm), ncam, MPI_INTEGER, &
         &            ipx_mdest,      4, mpi_comm_world, &
         &            irq(8), ierr)

    nrq = 8
    call mpi_waitall(nrq, irq, istatus, ierr)
#endif
#ifdef PROFILE_COMM
    call end_profile(comm_direct3_dr_id)
#endif

    !=== create i2m ===!
    ntmp=0
!$omp parallel default(none) &
!$omp& private(ncc,icx,icy,icz,ica,icag) &
!$omp& shared(lxdiv,lydiv,lzdiv,tag,na_per_cell) &
!$omp& shared(m2i,i2m) &
!$omp& reduction(+:ntmp)
!$omp do
    do icx=0,lxdiv+nbd
    do icy=0,lydiv+nbd
    do icz=0,lzdiv+nbd

       if(icx.ge.1 .and. icx.le.lxdiv .and. &
    &     icy.ge.1 .and. icy.le.lydiv .and. &
    &     icz.ge.1 .and. icz.le.lzdiv) cycle
       do ica=tag(icz,icy,icx), &
    &         tag(icz,icy,icx)+na_per_cell(icz,icy,icx)-1
          icag=m2i(ica)
          i2m(icag)=ica  !! bug without -DSYNC_COM
          ntmp=ntmp+1
          !         if(icag .le. -1) cycle
       end do ! ica
    end do ! icz
    end do ! icy
    end do ! icx
!$omp end do
!$omp end parallel
    ndatm=nselfatm+ntmp

  end subroutine comm_direct_3_dr

end module comm_direct3_dr
