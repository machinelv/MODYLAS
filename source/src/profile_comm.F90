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
!! \brief  Module and subroutines to profile MPI communications.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to profile MPI communications.
!! \author Kensuke Iwahashi
!<
module profile_comm
  use mpi_tool, only : myrank
#define NUM_PROFILES 5
  integer(4),parameter :: comm_direct3_dr_id = 1
  integer(4),parameter :: comm_direct2_dr_id = 2
  integer(4),parameter :: comm_fmm_local_super_id = 3
  integer(4),parameter :: comm_fmm_local_multi_dr_id = 4
  integer(4),parameter :: comm_fmm_local_top_id = 5
  real(8) :: profile_t_start(NUM_PROFILES)
  real(8) :: profile_t_end(NUM_PROFILES)
  real(8) :: profile_t(NUM_PROFILES) = 0
contains
  subroutine start_profile(id)
    implicit none
    include 'mpif.h'
    integer(4) :: id

    profile_t_start(id) = MPI_Wtime()
  end subroutine start_profile
  !
  subroutine end_profile(id)
    implicit none
    include 'mpif.h'
    integer(4) :: id

    profile_t_end(id) = MPI_Wtime()
    profile_t(id) = profile_t(id) + (profile_t_end(id) - profile_t_start(id))
  end subroutine end_profile
  !
  subroutine print_profile()
    implicit none
    include 'mpif.h'
    integer(4) :: size, ierr
    real(8) :: t_max(NUM_PROFILES), t_min(NUM_PROFILES), t_ave(NUM_PROFILES)
#ifdef FJ_RDMA
    character(len=5) :: comm_kind = "RDMA "
#else
    character(len=4) :: comm_kind = "MPI "
#endif


    call mpi_comm_size(mpi_comm_world, size, ierr)
    call mpi_allreduce(profile_t, t_max, NUM_PROFILES, mpi_double_precision, mpi_max, mpi_comm_world, ierr)
    call mpi_allreduce(profile_t, t_min, NUM_PROFILES, mpi_double_precision, mpi_min, mpi_comm_world, ierr)
    call mpi_allreduce(profile_t, t_ave, NUM_PROFILES, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
    t_ave(:) = t_ave(:) / size

    if(myrank == 0) then
       write(*,'(a,a,f10.4,a)') comm_kind, "comm_direct3_dr(MAX)", t_max(comm_direct3_dr_id), " sec."
       write(*,'(a,a,f10.4,a)') comm_kind, "comm_direct3_dr(AVE)", t_ave(comm_direct3_dr_id), " sec."
       write(*,'(a,a,f10.4,a)') comm_kind, "comm_direct3_dr(MIN)", t_min(comm_direct3_dr_id), " sec."
       write(*,*) "---"
       write(*,'(a,a,f10.4,a)') comm_kind, "comm_direct2_dr(MAX)", t_max(comm_direct2_dr_id), " sec."
       write(*,'(a,a,f10.4,a)') comm_kind, "comm_direct2_dr(AVE)", t_ave(comm_direct2_dr_id), " sec."
       write(*,'(a,a,f10.4,a)') comm_kind, "comm_direct2_dr(MIN)", t_min(comm_direct2_dr_id), " sec."
       write(*,*) "---"
       write(*,'(a,a,f10.4,a)') comm_kind, "comm_fmm_local_super(MAX)", t_max(comm_fmm_local_super_id), " sec."
       write(*,'(a,a,f10.4,a)') comm_kind, "comm_fmm_local_super(AVE)", t_ave(comm_fmm_local_super_id), " sec."
       write(*,'(a,a,f10.4,a)') comm_kind, "comm_fmm_local_super(MIN)", t_min(comm_fmm_local_super_id), " sec."
       write(*,*) "---"
       write(*,'(a,a,f10.4,a)') comm_kind, "comm_fmm_local_multi_dr(MAX)", t_max(comm_fmm_local_multi_dr_id), " sec."
       write(*,'(a,a,f10.4,a)') comm_kind, "comm_fmm_local_multi_dr(AVE)", t_ave(comm_fmm_local_multi_dr_id), " sec."
       write(*,'(a,a,f10.4,a)') comm_kind, "comm_fmm_local_multi_dr(MIN)", t_min(comm_fmm_local_multi_dr_id), " sec."
       write(*,*) "---" !ya 20191126
       write(*,'(a,a,f10.4,a)') comm_kind, "comm_fmm_local_top(MAX)", t_max(comm_fmm_local_top_id), " sec."
       write(*,'(a,a,f10.4,a)') comm_kind, "comm_fmm_local_top(AVE)", t_ave(comm_fmm_local_top_id), " sec."
       write(*,'(a,a,f10.4,a)') comm_kind, "comm_fmm_local_top(MIN)", t_min(comm_fmm_local_top_id), " sec."
    end if
  end subroutine print_profile
end module profile_comm
