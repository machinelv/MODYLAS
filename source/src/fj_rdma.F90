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
!! \brief  Module and subroutines which relate to RDMA communication (Fujitsu).
!<
!----------------------------------------------------------------------
!>
!! \brief Module which relates to RDMA communication (Fujitsu).
!! \author Shin-ichi Ichikawa
!<
module fj_rdma
  use iso_c_binding
  implicit none
  interface rdma_register_addr
     subroutine rdma_register_addr_i1(addr, size)
       integer(4) :: addr(:)
       integer(4), value :: size
     end subroutine rdma_register_addr_i1
     !
     subroutine rdma_register_addr_i3(addr, size)
       integer(4) :: addr(:,:,:)
       integer(4), value :: size
     end subroutine rdma_register_addr_i3
     !
     subroutine rdma_register_addr_r2(addr, size)
       real(8) :: addr(:,:)
       integer(4), value :: size
     end subroutine rdma_register_addr_r2
     !
     subroutine rdma_register_addr_c1(addr, size)
       complex(8) :: addr(:)
       integer(4), value :: size
     end subroutine rdma_register_addr_c1
     !
     subroutine rdma_register_addr_c4(addr, size)
       complex(8) :: addr(:,:,:,:)
       integer(4), value :: size
     end subroutine rdma_register_addr_c4
     subroutine rdma_register_addr_s5(addr, size)
       real(kind=selected_real_kind(6, 37)) :: addr(:,:,:,:,:)
       integer(4), value :: size
     end subroutine rdma_register_addr_s5
     !
     subroutine rdma_register_addr_r5(addr, size)
       real(kind=selected_real_kind(15, 307)) :: addr(:,:,:,:,:)
       integer(4), value :: size
     end subroutine rdma_register_addr_r5
  end interface rdma_register_addr
  !
  interface rdma_get_raddr
     function rdma_get_raddr_i1(addr)
       use,intrinsic :: iso_c_binding
       integer(4) :: addr(:)
       type(c_ptr) :: rdma_get_raddr_i1
     end function rdma_get_raddr_i1
     !
     function rdma_get_raddr_i3(addr)
       use,intrinsic :: iso_c_binding
       integer(4) :: addr(:,:,:)
       type(c_ptr) :: rdma_get_raddr_i3
     end function rdma_get_raddr_i3
     !
     function rdma_get_raddr_r2(addr)
       use,intrinsic :: iso_c_binding
       real(8) :: addr(:,:)
       type(c_ptr) :: rdma_get_raddr_r2
     end function rdma_get_raddr_r2
     !
     function rdma_get_raddr_c1(addr)
       use,intrinsic :: iso_c_binding
       complex(8) :: addr(:)
       type(c_ptr) :: rdma_get_raddr_c1
     end function rdma_get_raddr_c1
     !
     function rdma_get_raddr_c4(addr)
       use,intrinsic :: iso_c_binding
       complex(8) :: addr(:,:,:,:)
       type(c_ptr) :: rdma_get_raddr_c4
     end function rdma_get_raddr_c4
     function rdma_get_raddr_s5(addr)
       use,intrinsic :: iso_c_binding
       real(kind=selected_real_kind(6, 37)) :: addr(:,:,:,:,:)
       type(c_ptr) :: rdma_get_raddr_s5
     end function rdma_get_raddr_s5
     !
     function rdma_get_raddr_r5(addr)
       use,intrinsic :: iso_c_binding
       real(kind=selected_real_kind(15, 307)) :: addr(:,:,:,:,:)
       type(c_ptr) :: rdma_get_raddr_r5
     end function rdma_get_raddr_r5
  end interface rdma_get_raddr
  !
  interface rdma_get_laddr
     function rdma_get_laddr_i1(addr)
       use,intrinsic :: iso_c_binding
       integer(4) :: addr(:)
       integer(8) :: rdma_get_laddr_i1
     end function rdma_get_laddr_i1
     !
     function rdma_get_laddr_i3(addr)
       use,intrinsic :: iso_c_binding
       integer(4) :: addr(:,:,:)
       integer(8) :: rdma_get_laddr_i3
     end function rdma_get_laddr_i3
     !
     function rdma_get_laddr_r2(addr)
       use,intrinsic :: iso_c_binding
       real(8) :: addr(:,:)
       integer(8) :: rdma_get_laddr_r2
     end function rdma_get_laddr_r2
     !
     function rdma_get_laddr_c1(addr)
       use,intrinsic :: iso_c_binding
       complex(8) :: addr(:)
       integer(8) :: rdma_get_laddr_c1
     end function rdma_get_laddr_c1
     !
     function rdma_get_laddr_c4(addr)
       use,intrinsic :: iso_c_binding
       complex(8) :: addr(:,:,:,:)
       integer(8) :: rdma_get_laddr_c4
     end function rdma_get_laddr_c4
     function rdma_get_laddr_s5(addr)
       use,intrinsic :: iso_c_binding
       real(kind=selected_real_kind(6, 37)) :: addr(:,:,:,:,:)
       integer(8) :: rdma_get_laddr_s5
     end function rdma_get_laddr_s5
     !
     function rdma_get_laddr_r5(addr)
       use,intrinsic :: iso_c_binding
       real(kind=selected_real_kind(15, 307)) :: addr(:,:,:,:,:)
       integer(8) :: rdma_get_laddr_r5
     end function rdma_get_laddr_r5
  end interface rdma_get_laddr
  !
  interface
     subroutine rdma_init(fcomm)
       integer(4), value :: fcomm
     end subroutine rdma_init
     !
     subroutine rdma_finalize()
     end subroutine rdma_finalize
     !
     subroutine rdma_put(rank, raddr, laddr, length)
       integer(4), value :: rank, length
       integer(8), value :: raddr, laddr
     end subroutine rdma_put
     !
     subroutine rdma_put_post(rank, raddr, laddr, length)
       integer(4), value :: rank, length
       integer(8), value :: raddr, laddr
     end subroutine rdma_put_post
     !
     subroutine rdma_post(rank)
       integer(4), value :: rank
     end subroutine rdma_post
     !
     subroutine rdma_wait(rank)
       integer(4), value :: rank
     end subroutine rdma_wait
     !
     subroutine rdma_sync_all()
     end subroutine rdma_sync_all
     !
     subroutine rdma_get_remote_addr(memid, raddr)
       integer(4), value :: memid
       integer(8) :: raddr(:)
     end subroutine rdma_get_remote_addr
     !
     function rdma_allocate_coarray(length, offset)
       use, intrinsic :: iso_c_binding, only: c_ptr
       integer(4), value  :: length
       integer(8) :: offset
       TYPE(C_PTR) :: rdma_allocate_coarray
     end function rdma_allocate_coarray
     !
     subroutine rdma_deallocate_coarray(num)
       integer(4), value :: num
     end subroutine rdma_deallocate_coarray
     !
     function rdma_get_heap_laddr()
       integer(8)  :: rdma_get_heap_laddr
     end function rdma_get_heap_laddr
     !
     function rdma_get_heap_raddr()
       use,intrinsic :: iso_c_binding, only: c_ptr
       type(c_ptr) :: rdma_get_heap_raddr
     end function rdma_get_heap_raddr
  end interface
  
  integer(8) :: wkxyz_laddr
  integer(8),pointer :: wkxyz_raddr(:)
end module fj_rdma
