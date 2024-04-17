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
!! \brief Module and Subroutine to calculate short-ranged interactions.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to calculate short-ranged interactions
!! \author Yoshimichi Andoh
!<
module force_short
  use omp_lib
  implicit none

contains

!======================================================================
!>
!! \brief  Subroutine to calculate short-ranged interactions
!! \author Yoshimichi Andoh
!<
  subroutine md_calculate_forces_short(ftmp,virtmp,enetmp)
!======================================================================
    use trajectory_mpi
    use atom_virial
    use lj_mod
    use forces
    use md_monitors
    use md_condition
    use mpi_tool
    use comm_bound_mod
    use unit_cell
    use bond, only : nbond_as, add_charmm_bond_a
    use bond_morse, only: nbondmorse_as, add_bond_morse
    use bond_morse2, only: nbondmorse2_as, add_bond_morse2
    use angle, only : nangleAs, add_charmm_angle_a
    use angle_morse, only : nanglemorseAs, add_angle_morse
    use UB, only : nub_as, add_charmm_ub_a
    use dihedral, only : ndihedralAs, add_charmm_dihedral_a
    use CMAP, only : ncmap, add_charmm_cmap_a
    use improper_torsion, only : nitorsionAs, add_charmm_itorsion_a
    use position_constrain
    use openmp_tool, only : nomp
#include "timing.h90"
    use comm_direct2_dr, only : comm_direct_2_dr
    use dr_cntl, only : nbd
    implicit none
    integer(4) :: i,i0,k0
    integer(4) :: iam
    integer(4) :: iy,ix
    real(8) :: ftmp(3,nadirect)
    real(8) :: virtmp(6),enetmp
    include 'mpif.h'

#ifdef DEBUGFCE
    if(myrank==0) write(*,*) '###  Intra-molecule   potE ###'
#endif

    call init_comm_buffer()

    TIME_START(TM_CHARMM_BOND)
    if(nbond_as.ne.0) then
       call add_charmm_bond_a()
    endif
    TIME_STOP(TM_CHARMM_BOND)
    if(nbondmorse_as.ne.0) then
       call add_bond_morse()
    endif
    if(nbondmorse2_as.ne.0) then
       call add_bond_morse2()
    endif
    TIME_START(TM_CHARMM_ANGLE)
    if(nangleAs.ne.0) then
       call add_charmm_angle_a()
    endif
    TIME_STOP(TM_CHARMM_ANGLE)
    if(nanglemorseAs.ne.0) then
       call add_angle_morse()
    endif
    if(nub_as.ne.0) then
       call add_charmm_ub_a()
    endif
    TIME_START(TM_CHARMM_DIHEDRAL)
    if(ndihedralAs.ne.0) then
       call add_charmm_dihedral_a()
    endif
    TIME_STOP(TM_CHARMM_DIHEDRAL)
    if(ncmap.ne.0) then
       call add_charmm_CMAP_a()
    endif
    if(nitorsionAs.ne.0) then
       call add_charmm_itorsion_a()
    endif

    call add_position_constrain_forces()

!$omp parallel default(none) &
!$omp& private(iy,ix) &
!$omp& shared(lxdiv,lydiv,lzdiv) &
!$omp& private(k0,i0,iam) &
!$omp& shared(tag,na_per_cell) &
!$omp& shared(lsegtop,lseg_natoms,nselfseg) &
!$omp& shared(nomp,ftmp,w3_f,nadirect)
    do ix=1-nbd,lxdiv+nbd
    do iy=1-nbd,lydiv+nbd
!$omp do
       do i0=tag(    1-nbd,iy,ix),   &
     &      tag(lzdiv+nbd,iy,ix)+na_per_cell(lzdiv+nbd,iy,ix)-1
          ftmp(1:3,i0)=0d0
          do iam = 0,nomp-1
             ftmp(1,i0) = ftmp(1,i0) + w3_f(1,i0,iam)
             ftmp(2,i0) = ftmp(2,i0) + w3_f(2,i0,iam)
             ftmp(3,i0) = ftmp(3,i0) + w3_f(3,i0,iam)
          end do ! iam
       end do ! i0
!$omp end do nowait
    end do ! k0, or iy
    end do ! ix
!$omp end parallel

    wk_vir=0d0
    do iam = 0,nomp-1
       do i = 1,6
          wk_vir(i) = wk_vir(i) + wk_vir2(i,iam)
       end do
    end do

    virtmp=wk_vir
    enetmp=wk_p_energy

#ifdef DEBUGFCE
    call comm_direct_2_dr(ftmp)
    open(myrank+10000)
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          write(myrank+10000,'(i7.7,3es23.15)') m2i(i0),ftmp(1:3,i0)
       enddo ! i0
    enddo ! k0
    close(myrank+10000); call flush(myrank+10000)
#endif
  end subroutine md_calculate_forces_short

end module force_short
