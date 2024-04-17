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
!! \brief Module and Subroutine to calculate middle-ranged interactions.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to calculate middle-ranged interactions
!! \author Yoshimichi Andoh
!<
module force_middle
  use omp_lib
  use fmm_near_dr, only : energy_direct_dr
  use pme_near_dr, only : energy_direct_pme_dr

contains

!======================================================================
!>
!! \brief  Subroutine to calculate middle-ranged interactions
!! \author Yoshimichi Andoh
!<
  subroutine md_calculate_forces_middle(ftmp,virtmp,enetmp)
!======================================================================
    use trajectory_mpi
    use atom_virial
    use lj_mod
    use forces
    use md_monitors
    use md_periodic
    use force_field_numbers
    use md_condition
    use mpi_tool
    use comm_bound_mod
    use unit_cell
    use void123
    use special14
    use md_condition
    use openmp_tool, only : nomp
#ifdef DEBUGFCE
    use comm_direct2_dr, only : comm_direct_2_dr
#endif
    use dr_cntl, only : nbd

#include "timing.h90"
    implicit none
    integer(4) :: i,i0,k0
    integer(4) :: iam
    integer(4) :: iy,ix
    real(8) :: ftmp(3,nadirect)
    real(8) :: virtmp(6),enetmp
    include 'mpif.h'

#ifdef DEBUGFCE
    if(myrank==0) write(*,*) '###  Two-body (short) potE ###'
#endif
    TIME_START(TM_FORCES_MIDDLE_BUFINIT)

    call init_comm_buffer()

    TIME_STOP(TM_FORCES_MIDDLE_BUFINIT)

    if (    md_periodic__type == FMM) then

       TIME_START(TM_ENERGY_DIRECT)
       call energy_direct_dr(wkxyz,w3_f)
       TIME_STOP(TM_ENERGY_DIRECT)

    elseif (md_periodic__type == EWALD    .or.     &
  &         md_periodic__type == PMEWALD) then

       TIME_START(TM_PME_NEAR)
       call energy_direct_pme_dr(wkxyz,w3_f)
       TIME_STOP(TM_PME_NEAR)

    else

       write(0,*) 'Warning: This modylas is FMM/PME/Ewald version.'
       write(0,*) '       : Only FMM/PME/Ewald is supported.'
       call modylas_abort()

    endif

    TIME_START(TM_RMVOID_ETC)

    if (    md_periodic__type /= FMM) then
        call remove_void123lj()  ! note: #ifdef TABLE presumed.
        call remove_void123cl()  ! note: #ifdef TABLE presumed.
    endif

    if(     md_condition__force_field == CHARMM)then
       call remove_special14lj() !! for CHARMM
#ifdef DEBUGFCE
       call dummy_calc14cl()     !! for CHARMM
#endif
    else if(md_condition__force_field == OPLSAA .or. &
         &  md_condition__force_field == AMBER  .or. &
         &  md_condition__force_field == GAmbFF   )then
       call remove_scaling14lj() !! for AMBER/OPLS
       call remove_scaling14cl() !! for AMBER/OPLS
    else if(md_condition__force_field == KREMER )then
       call remove_special14kremer() !! for Kremer's pot
    endif
    TIME_STOP(TM_RMVOID_ETC)
    TIME_START(TM_RMVOID_ETC)

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
     &       tag(lzdiv+nbd,iy,ix)+na_per_cell(lzdiv+nbd,iy,ix)-1
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

    TIME_STOP(TM_RMVOID_ETC)

#ifdef DEBUGFCE
    call comm_direct_2_dr(ftmp)
    plj=0d0;pcl=0d0
    call mpi_allreduce(wplj,plj,1, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,i)
    call mpi_allreduce(wpcl,pcl,1, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,i)
#ifdef KCAL
    if(myrank==0) write(*,*) 'Pot(LJ_sr)=', (plj)*kJ_mol/4.184d0,'[kcal/mol] Note: LJ14 not included'
#else
    if(myrank==0) write(*,*) 'Pot(LJ_sr)=', (plj)*kJ_mol,'[kJ/mol] Note: LJ14 not included'
#endif
#ifdef KCAL
    if(myrank==0) write(*,*) 'Pot(CL_sr)=', (pcl)*kJ_mol/4.184d0,'[kcal/mol] Note: CL14 not included'
#else
    if(myrank==0) write(*,*) 'Pot(CL_sr)=', (pcl)*kJ_mol,'[kJ/mol] Note: CL14 not included'
#endif
       if(LJ_LRC)then
       if(myrank==0) write(*,*) 'Pot(LJcor)=',&
#ifdef KCAL
      &                corrector_potential/cellvol*kJ_mol/4.184d0,'[kcal/mol]'
#else
      &                corrector_potential/cellvol*kJ_mol,'[kJ/mol]'
#endif
       endif

    open(myrank+20000)
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          write(myrank+20000,'(i7.7,3es23.15)') m2i(i0),ftmp(1:3,i0)
       enddo ! i0
    enddo ! k0
    close(myrank+20000); call flush(myrank+20000)
#endif
  end subroutine md_calculate_forces_middle

end module force_middle
