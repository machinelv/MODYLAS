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
!! \brief Module and Subroutine to calculate long-ranged interactions.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to calculate long-ranged interactions
!! \author Yoshimichi Andoh
!<
module force_long
  use omp_lib

contains

!======================================================================
!>
!! \brief  Subroutine to calculate long-ranged interactions
!! \author Yoshimichi Andoh
!<
  subroutine md_calculate_forces_long(ftmp,virtmp,enetmp)
!======================================================================
    use trajectory_mpi
    use atom_virial
    use lj_mod
    use forces
    use md_monitors
    use md_periodic
    use ewald_variables
    use ewald_mod
    use fmm_far
    use pme_far
    use mpi_tool
    use comm_bound_mod
    use unit_cell
    use ensemble_numbers
    use md_condition, only : md_condition__ensemble
    use surface_term

    use physics_const
    use nonneutral

#include "timing.h90"
    use openmp_tool, only : nomp
    use comm_direct2_dr, only : comm_direct_2_dr
    use dr_cntl, only : nbd
    implicit none
    integer(4) :: i,i0,k0
    integer(4) :: iam
    integer(4) :: iy,ix
    real(8) :: ftmp(3,nadirect)
    real(8) :: virtmp(6),enetmp

!   include 'mpif.h'   !! module contain this line 
    integer(4)::ierr

#ifdef DEBUGFCE
    if(myrank==0) write(*,*) '###  Two-body (Long)  potE ###'
#endif

    TIME_START(TM_FLONG_PRE)

    call init_comm_buffer()

    TIME_STOP(TM_FLONG_PRE)

    if (md_periodic__type == FMM) then

       if(md_condition__ensemble==NPT_PR) then
          if(myrank==0) write(0,*)'ERROR: NPT_PR with FMM not supported'
          call modylas_abort()
       endif

       call calc_mm               ! P2M, M2M
       call energy_fmm()          ! fmm Ewald, M2L, L2L, L2P
       if(.not. ewald_sterm) then ! If sterm is off (default in FMM)

          TIME_START(TM_RM_EWALD_SURF)
          call calc_system_dipole
          call remove_ewald_surface_term
          TIME_STOP(TM_RM_EWALD_SURF)
          TIME_START(TM_BCG_SURF)
          if(bNonNeutrality)  call remove_charged_ewaldsurfaceterm
          TIME_STOP(TM_BCG_SURF)

       endif
    elseif (md_periodic__type == EWALD) then

       call md_add_coulomb_ewald_const()
       call md_add_ewald_lattice()

    elseif (md_periodic__type == PMEWALD) then

       call md_add_coulomb_ewald_const()
       TIME_START(TM_PME)
       call md_add_pmewald()
       TIME_STOP(TM_PME)

    endif

    TIME_START(TM_FLONG_POST)
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

    TIME_STOP(TM_FLONG_POST)

#ifdef DEBUGFCE
    call comm_direct_2_dr(ftmp)
    open(myrank+30000)
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          write(myrank+30000,'(i7.7,3es23.15)') m2i(i0),ftmp(1:3,i0)
       enddo ! i0
    enddo ! k0
    call flush(myrank+30000) !;call mpiend
#endif /** DEBUGFCE **/
  end subroutine md_calculate_forces_long

end module force_long
