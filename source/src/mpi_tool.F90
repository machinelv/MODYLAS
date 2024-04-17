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
!! \brief  Module and subroutines to initialize/finalize MPI
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to initialize/finalize MPI
!! \author Kensuke Iwahashi, Yoshimichi Andoh
!<
module mpi_tool
#ifdef FJ_RDMA
  use fj_rdma
#endif
  integer(4) :: myrank=0, nprocs=1, mpiout=0

  ! for mpi_wtime
  real(8) :: time0,time1,time2
  integer :: on !=10

contains

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to initialize MPI
!! \author Kensuke Iwahashi, Yoshimichi Andoh
!<
  subroutine mpistart
!----------------------------------------------------------------------
    implicit none
!#ifdef MPIPARA
    include 'mpif.h'
    integer(4) :: ierr

    call mpi_init(ierr)
    call mpi_comm_size(mpi_comm_world,nprocs,ierr)
    call mpi_comm_rank(mpi_comm_world,myrank,ierr)
#ifdef FJ_RDMA
    call rdma_init(mpi_comm_world)
#endif
    mpiout=0

  end subroutine mpistart

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to finalize MPI
!! \author Kensuke Iwahashi, Yoshimichi Andoh
!<
  subroutine mpiend
!----------------------------------------------------------------------
    implicit none
    include 'mpif.h'
    integer ierr
#ifdef FJ_RDMA
    call rdma_finalize
#endif
    call mpi_finalize(ierr)
    stop
  end subroutine mpiend

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to abort modylas execution.
!! \author Kensuke Iwahashi
!<
  subroutine modylas_abort
!----------------------------------------------------------------------
    implicit none
    include 'mpif.h'
    integer ierr
    call flush(0)  !! output buffer contents to stderr !!
    if(myrank==0) then
    write(*,*) '***********************************************************'
    write(*,*) '* Program was forced to stop by calling modylas_abort.    *'
    write(*,*) '* About its reason, see also standard error file.         *'
    write(*,*) '* On Linux, the ERROR message may be output to a console. *'
    write(*,*) '* On Supercomputers, check ./*.sh.eJOBID, ./stderr, etc.  *'
    write(*,*) '***********************************************************'
    endif
    call mpi_abort(mpi_comm_world,ierr)
  end subroutine modylas_abort

end module mpi_tool
