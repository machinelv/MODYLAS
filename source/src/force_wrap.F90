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
!! \brief Module and Subroutine to wrap a series of force routines.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to wrap a series of force routines.
!! \author Yoshimichi Andoh
!<
module force_wrap
  use omp_lib
  use mpi_tool, only : mpiend
#include "timing.h90"
  implicit none

contains

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to wrap force routines, being called by 
!!         numerical integration code, except for npt ensemble code.
!! \author Yoshimichi Andoh
!<
subroutine md_calculate_forces()
  use forces,       only : wk_f
  use force_short,  only : md_calculate_forces_short
  use force_middle, only : md_calculate_forces_middle
  use force_long,   only : md_calculate_forces_long
  use comm_direct2_dr, only : comm_direct_2_dr
  use comm_direct3_dr, only : comm_direct_3_dr
  use md_monitors,     only : wk_p_energy
  use md_multiplestep, only : MTm, MTl, maxMTm, maxMTl, scaleM, scaleL, &
                              fshort,fmiddle,flong, &
                              virshort,virmiddle,virlong, &
                              eneshort,enemiddle,enelong
  use md_periodic
#ifdef TIP4
  use tip4p, only : distribute_force_on_msite
#endif
  use trajectory_mpi,  only : nadirect
  use subcell, only : nselfseg, lsegtop, lseg_natoms, m2i
  use mpi_tool, only : myrank
#ifdef DEBUGFCE
  use unit_cell, only : cellvol
#endif
  implicit none
  integer(4) i0,k0
#ifdef DEBUGFCE
    real(8)::Us,Um,Ul,Ut
    real(8)::virs(6),virm(6),virl(6),virt(6),scavir
    integer(4)::ierr
    real(8),parameter :: kJ_mol=6.02214129d+23*1.0d-3
    real(8),parameter :: kcal_mol=6.02214129d+23*1.0d-3/4.184d0
    include 'mpif.h'

    if(myrank==0) write(*,*) '------------------------------', &
    &                '--------------------------------------'
    if(myrank==0) write(*,*) 'Information for debugging:    ', &
    &                ' (These only appear with -DDEBUGFCE)  '
    if(myrank==0) write(*,*) '------------------------------', &
    &                '--------------------------------------'
#endif

          TIME_START(TM_FSHORT)
  !!    ^^^ short ^^^
  call md_calculate_forces_short(fshort ,virshort ,eneshort)
          TIME_STOP(TM_FSHORT)

  !!    ^^^ middle ^^^
  if(MTm==maxMTm)then
    scaleM=1d0
    call md_calculate_forces_middle(fmiddle,virmiddle,enemiddle)
  else
    scaleM=0d0
  endif
  !!    ^^^ long ^^^
  if(MTl==maxMTl.and.MTm==maxMTm)then
    scaleL=1d0
    call md_calculate_forces_long(flong  ,virlong  ,enelong)
  else
    scaleL=0d0
  endif

    TIME_START(TM_WRAP_FORCES_LOCAL)
  !!    ^^^ sum short, middle, and long ^^^
  wk_p_energy=eneshort+enemiddle+enelong

!$omp parallel do default(shared) &
!$omp& private(k0,i0)
!default: MTD
     do i0=1,nadirect
        wk_f(:,i0)=fshort(:,i0)+scaleM*fmiddle(:,i0)*maxMTm &
    &                          +scaleL*flong(:,i0)  *maxMTm*maxMTl
     enddo
    TIME_STOP(TM_WRAP_FORCES_LOCAL)

!default: MTD
    TIME_BARRIER(TMB_COMM_DIRECT2_FSHT_MDL_LNG)
    TIME_START(TM_COMM_DIRECT2_FSHT_MDL_LNG)

  call comm_direct_2_dr(wk_f) 

    TIME_STOP(TM_COMM_DIRECT2_FSHT_MDL_LNG)

#ifdef TIP4
  call distribute_force_on_msite(wk_f)
#endif
 
#ifdef DEBUGFCE
!potential energy
    call mpi_allreduce(eneshort,Us,1, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    call mpi_allreduce(enemiddle,Um,1, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    call mpi_allreduce(enelong,Ul,1, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    call mpi_allreduce(wk_p_energy,Ut,1, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    if(myrank==0) write(*,*) '------------------------------', &
    &                '--------------------------------------'
    if(myrank==0) write(*,*) 'Pot(tot),      Pot(short),    ', &
    &                        'Pot(middle),   Pot(long)      '
    if(myrank==0) write(*,*) '------------------------------', &
    &                '--------------------------------------'
    if(myrank==0) write(*,'(4es15.7,a)') &
#ifdef KCAL
    &     Ut*kcal_mol,Us*kcal_mol,Um*kcal_mol,Ul*kcal_mol,' [kcal/mol]'
#else
    &     Ut*kJ_mol,Us*kJ_mol,Um*kJ_mol,Ul*kJ_mol,' [kJ/mol]'
#endif
    if(myrank==0) write(*,*) '------------------------------', &
    &                '--------------------------------------'
    if(myrank==0) write(*,*) 'Pot(tot)   : Total sum of potential energy'
    if(myrank==0) write(*,*) 'Pot(short) : Sum of intra-molecule   terms'
    if(myrank==0) write(*,*) 'Pot(middle): Sum of two-body (short) terms'
    if(myrank==0) write(*,*) 'Pot(long)  : Sum of two-body (long)  terms'
!virial
    call mpi_allreduce(virshort, virs,6, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    call mpi_allreduce(virmiddle,virm,6, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    call mpi_allreduce(virlong,  virl,6, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    virt(1:6)=virs(1:6)+virm(1:6)+virl(1:6)
    scavir=virt(1)+virt(2)+virt(3)
    if(myrank==0) write(*,*) '------------------------------', &
    &                '--------------------------------------'
#ifdef ATM
    if(myrank==0) write(*,*) 'Pstatic [atm]'
    if(myrank==0) write(*,'(es15.7)')  scavir/(3d0*cellvol)/101325d0
#else
    if(myrank==0) write(*,*) 'Pstatic [Pa]'
    if(myrank==0) write(*,'(es15.7)')  scavir/(3d0*cellvol)
#endif
if (    md_periodic__type /= FMM) then
#ifdef ATM
    if(myrank==0) write(*,*) 'Ptensor [atm]'
    if(myrank==0) write(*,'(3es15.7,a)') virt(1:3)/cellvol/101325d0, '   Pxx, Pyy, Pzz'
    if(myrank==0) write(*,'(3es15.7,a)') virt(4:6)/cellvol/101325d0, '   Pxy, Pxz, Pyz'
#else
    if(myrank==0) write(*,*) 'Ptensor [Pa]'
    if(myrank==0) write(*,'(3es15.7,a)') virt(1:3)/cellvol, '   Pxx, Pyy, Pzz'
    if(myrank==0) write(*,'(3es15.7,a)') virt(4:6)/cellvol, '   Pxy, Pxz, Pyz'
#endif
endif
    if(myrank==0) write(*,*) 'Note: These P values are derived from potentian energy only.'
    if(myrank==0) write(*,*) '      i.e., virial from constraint forces are not included. '
    if(myrank==0) write(*,*) '------------------------------', &
    &                '--------------------------------------'
!force
    open(myrank+40000)
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          write(myrank+40000,'(i7.7,3es23.15)') m2i(i0), &
          fshort(1:3,i0)+fmiddle(1:3,i0)+flong(1:3,i0)
       enddo ! i0
    enddo ! k0
    close(myrank+40000); call flush(myrank+40000)
!  call mpiend()
#endif

end subroutine md_calculate_forces

end module force_wrap
