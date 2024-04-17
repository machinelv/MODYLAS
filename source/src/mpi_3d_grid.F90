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
!! \brief  Module and subroutines to set number of MPI process
!!         along each axis.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to set MPI process number along each axis.
!! \author Shin-ichi Ichikawa, Yoshimichi Andoh
!<
module mpi_3d_grid
  implicit none
  integer, protected :: npz,npy,npx
  integer, protected :: ipz,ipy,ipx
  logical :: mpi_manual_division_flg=.false.

contains

!>
!! \brief  Subroutine to wrap set_npx,y,z routines.
!! \author Shin-ichi Ichikawa
!<
  subroutine set_np_x_y_z(npx_in, npy_in, npz_in)
    integer(4), intent(in) :: npx_in, npy_in, npz_in
    call set_npx(npx_in)
    call set_npy(npy_in)
    call set_npz(npz_in)
  end subroutine set_np_x_y_z

!>
!! \brief  Subroutine to wrap set_ipx,y,z routines.
!! \author Shin-ichi Ichikawa
!<
  subroutine set_ip_x_y_z(ipx_in, ipy_in, ipz_in)
    integer(4), intent(in) :: ipx_in, ipy_in, ipz_in
    ipx = ipx_in
    ipy = ipy_in
    ipz = ipz_in
  end subroutine set_ip_x_y_z

  subroutine set_npx(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    npx = ivalue
  end subroutine set_npx
!----------------------------------------------------------------------
  subroutine set_npy(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    npy = ivalue
  end subroutine set_npy
!----------------------------------------------------------------------
  subroutine set_npz(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    npz = ivalue
  end subroutine set_npz

  integer(4) function get_ip(ipx_in, ipy_in, ipz_in)
    implicit none
    integer(4), intent(in) :: ipx_in, ipy_in, ipz_in
    get_ip = (ipz_in-1)*npx*npy + (ipy_in-1)*npx + (ipx_in - 1)
  end function get_ip

!>
!! \brief  Subroutine to get npx, npy, npz automatically.
!! \author Shin-ichi Ichikawa
!<
  subroutine init_np_x_y_z()
    use mpi_tool, only : nprocs
    implicit none
    integer(4) :: nprocx,nprocy,nprocz
    integer(4) :: maxdiv
    integer(4) :: mindiv !,mxdiv,mydiv,mzdiv

    nprocx=1
    nprocy=1
    nprocz=1

    maxdiv = nprocs

    DO WHILE (mod(maxdiv,3).eq.0)
       !!== 3-power ==!!
       nprocx = nprocx*3
       maxdiv = maxdiv/3
       if(mod(maxdiv,3).eq.0) then
          nprocy = nprocy*3
          maxdiv = maxdiv/3
       end if
       if(mod(maxdiv,3).eq.0) then
          nprocz = nprocz*3
          maxdiv = maxdiv/3
       end if
    END DO

    !!== 2-power ==!!
    DO WHILE (mod(maxdiv,2).eq.0)
       mindiv=min(nprocx,nprocy,nprocz)
       if(    mindiv==nprocx)then
          nprocx = nprocx*2
          maxdiv = maxdiv/2
       elseif(mindiv==nprocy)then
          nprocy = nprocy*2
          maxdiv = maxdiv/2
       else
          nprocz = nprocz*2
          maxdiv = maxdiv/2
       endif
    END DO
    
    call set_np_x_y_z(nprocx, nprocy, nprocz)
  end subroutine init_np_x_y_z

!>
!! \brief  Subroutine to get ipx, ipy, ipz automatically.
!! \author Shin-ichi Ichikawa
!<
  subroutine init_ip_x_y_z()
    use mpi_tool, only : myrank
    implicit none
    integer(4) :: iprocx,iprocy,iprocz
    
    iprocx = mod(myrank,npx)
    iprocy = mod((myrank-iprocx)/npx,npy)
    iprocz = mod((myrank-iprocx-iprocy*npx)/(npx*npy),npz)
    call set_ip_x_y_z(iprocx, iprocy, iprocz)
  end subroutine init_ip_x_y_z

!>
!! \brief  Subroutine to initialize ipx, ipy, ipz.
!! \author Shin-ichi Ichikawa
!<
  subroutine init_mpi_3d_grid()
    implicit none

    if (.not. mpi_manual_division_flg)  call init_np_x_y_z()

    call init_ip_x_y_z()
  end subroutine init_mpi_3d_grid

end module mpi_3d_grid
