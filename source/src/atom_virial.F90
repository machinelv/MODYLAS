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
!! \brief  Module to store virial information.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to store virial infomation.
!! \author Kensuke Iwahashi, Yoshimichi ANDOH
!<
module atom_virial
  implicit none
  real(8) :: wk_vir(6)
  real(8),allocatable :: wk_vir2(:,:)  !! openmp threads
  real(8) :: virialSca=0d0
  real(8) :: wkvirialSca=0d0
  real(8),allocatable :: wkvirScaS(:)
  real(8),allocatable :: wkvirScaM(:)
  real(8),allocatable :: wkvirScaL(:)
  real(8) :: vir_part(3)=0d0     !! 3 means (short,middle,long)
  real(8) :: wkvir_part(3)=0d0   !! 3 means (short,middle,long)
  real(8) :: virialTen(6)
  real(8) :: wkvirialTen(6)
  real(8),allocatable :: wkvirTenS(:,:)
  real(8),allocatable :: wkvirTenM(:,:)
  real(8),allocatable :: wkvirTenL(:,:)
  real(8) :: vir_partTen(6,3)    !! 3 means (short,middle,long)
  real(8) :: wkvir_partTen(6,3)  !! 3 means (short,middle,long)
  real(8) :: sum_vir_part(3)=0d0 !! 3 means (short,middle,long)
  integer(4) :: icount_vir=0
#ifdef DEBUGFCE
  real(8)::wplj,plj,wpcl,pcl
  real(8)::wljv(6),ljv(6)
  real(8),parameter :: kJ_mol=6.02214129d+23*1.0d-3
#endif

contains

!>
!! \brief  Subroutine to allocate virial arrays.
!! \author Kensuke Iwahashi, Yoshimichi ANDOH
!<
  subroutine fmod_alloc_atom_virial
    use omp_lib
    use md_multiplestep
    use trajectory_mpi, only : nadirect
    implicit none
    integer(4) :: nomp
    nomp = 1
!$  nomp = omp_get_max_threads()
    allocate(wk_vir2(6,0:nomp-1)) ! for interaction
    allocate(wkvirScaS(2*maxMTm*maxMTl))
    allocate(wkvirScaM(2*maxMTm*maxMTl))
    allocate(wkvirScaL(2*maxMTm*maxMTl))
    wkvirScaS=0d0;wkvirScaM=0d0;wkvirScaL=0d0
    allocate(wkvirTenS(6,2*maxMTm*maxMTl))
    allocate(wkvirTenM(6,2*maxMTm*maxMTl))
    allocate(wkvirTenL(6,2*maxMTm*maxMTl))
    wkvirTenS=0d0;wkvirTenM=0d0;wkvirTenL=0d0
  end subroutine fmod_alloc_atom_virial

end module atom_virial
