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
!! \brief  Module and subroutines to store atom coordinates 
!!         and velocities distributed to each MPI process.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to store atom coordinates and velocities
!!         distributed to each MPI process.
!! \author Yoshimichi Andoh and Shin-ichi Ichikawa
!<
module trajectory_mpi
  implicit none
  real(8),allocatable :: wkxyz(:,:)
  real(8),allocatable :: wkv(:,:)
  integer(4) :: na1cell,na5cell,nadirect, nalongcell
  integer(4) :: na1cell_init = 250
  integer(4) :: naline,narea
#if defined(PRECISION_P2P_SP)
  real(4),allocatable :: sp_wkxyz(:,:)
#endif
contains
!----------------------------------------------------------------------
  subroutine fmod_set_na1cell(ivalue)
    use mpi_tool
    implicit none
    integer(4), intent(in) :: ivalue
    na1cell_init = ivalue
    if(myrank==0) write(*,*) 'WARNING: na1cell is set to ', ivalue
  end subroutine fmod_set_na1cell
!----------------------------------------------------------------------
end module trajectory_mpi
