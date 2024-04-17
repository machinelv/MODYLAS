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
!! \brief  Module and subroutines to integrate the equation of motion
!!         of barostat.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module containing subroutines to update barostat.
!! \author Yoshimichi Andoh
!<
!----------------------------------------------------------------------
module barostat
  use mpi_tool
  implicit none

contains
!----------------------------------------------------------------------
!>
!! \brief  Subroutine for barostat located at outermost
!! \author Yoshimichi Andoh
!<
!! [Ref.] M.E.Tuckerman et al., J.Phys.A:Math.Gen.,39,5629(2006).
!----------------------------------------------------------------------
  subroutine update_vboxg_outermost(dtin)  ! long
    use extended_system
    use thermostat
    use md_multiplestep
    use md_condition
    use md_const
    use unit_cell
    use kinetic_energy
#include "timing.h90"
    implicit none
    real(8) :: totke,kinT
    real(8) :: pisov
    real(8) :: gboxg(3,3)
    real(8) :: dtin,virSca
    real(8) :: wkvirlong
    include 'mpif.h'
    integer(4) :: ierr

    !     ^^^ temperature ^^^
    call k_energy_scaler(totke)
    kinT = 2.0d0*totke*rvkbolz*degree_of_freedom_inverse

    !     ^^^ virial (long) ^^^
    wkvirlong=virlong(1)+virlong(2)+virlong(3)
    virSca=0d0
    TIME_BARRIER(TMB_ALLREDUCE_VIR)
    TIME_START(TM_ALLREDUCE_VIR)
    call mpi_allreduce(wkvirlong,virSca,1, &
         &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    TIME_STOP(TM_ALLREDUCE_VIR)

    !     ^^^ piston force (long) ^^^
    pisov=( 2.d0*totke+virSca ) * onethird
    gboxg(1,1)= (pisov-Poutermost*cellvol)*rvkbolz+kinT
    gboxg(2,2)= gboxg(1,1)  ! set same value to avoid bugs by ifort
    gboxg(3,3)= gboxg(1,1)

    !     ^^^ update cell velocity (long) ^^^
    vboxg(1,1)=vboxg(1,1)+0.5d0*dtin*gboxg(1,1)*rbaromass
    vboxg(1,2)=0d0
    vboxg(1,3)=0d0
    vboxg(2,1)=0d0
    vboxg(2,2)=vboxg(1,1)  ! set same value to avoid bugs by ifort
    vboxg(2,3)=0d0
    vboxg(3,1)=0d0
    vboxg(3,2)=0d0
    vboxg(3,3)=vboxg(1,1)  ! set same value to avoid bugs by ifort

  end subroutine update_vboxg_outermost

!----------------------------------------------------------------------
!>
!! \brief  Subroutine for barostat located at outer
!! \author Yoshimichi Andoh
!<
!! [Ref.] M.E.Tuckerman et al., J.Phys.A:Math.Gen.,39,5629(2006).
  subroutine update_vboxg_outer(dtin)  ! Middle
!----------------------------------------------------------------------
    use extended_system
    use thermostat
    use md_multiplestep
    use md_const
    use unit_cell
#include "timing.h90"
    implicit none
    real(8) :: pisov
    real(8) :: gboxg(3,3)
    real(8) :: dtin,virSca
    real(8) :: wkvirmiddle
    include 'mpif.h'
    integer(4) :: ierr

    !     ^^^ virial (middle) ^^^
    wkvirmiddle=virmiddle(1)+virmiddle(2)+virmiddle(3)
    virSca=0d0
    TIME_BARRIER(TMB_ALLREDUCE_VIR)
    TIME_START(TM_ALLREDUCE_VIR)
    call mpi_allreduce(wkvirmiddle,virSca,1, &
         &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    TIME_STOP(TM_ALLREDUCE_VIR)

    !     ^^^ piston force (middle) ^^^
    pisov=(           +virSca ) * onethird
    gboxg(1,1)= (pisov-Pouter*cellvol)*rvkbolz
    gboxg(2,2)= gboxg(1,1)
    gboxg(3,3)= gboxg(1,1)

    !     ^^^ update cell velocity (middle) ^^^
    vboxg(1,1)=vboxg(1,1)+0.5d0*dtin*gboxg(1,1)*rbaromass
    vboxg(1,2)=0d0
    vboxg(1,3)=0d0
    vboxg(2,1)=0d0
    vboxg(2,2)=vboxg(1,1)
    vboxg(2,3)=0d0
    vboxg(3,1)=0d0
    vboxg(3,2)=0d0
    vboxg(3,3)=vboxg(1,1)

  end subroutine update_vboxg_outer
!----------------------------------------------------------------------
!>
!! \brief  Subroutine for barostat located at inner
!! \author Yoshimichi Andoh
!<
!! [Ref.] M.E.Tuckerman et al., J.Phys.A:Math.Gen.,39,5629(2006).
  subroutine update_vboxg_inner(dtin)  ! Short
!----------------------------------------------------------------------
    use extended_system
    use thermostat
    use md_multiplestep
    use md_const
    use atom_virial
    use unit_cell
    use pressure_mod
    implicit none
    real(8) :: pisov
    real(8) :: gboxg(3,3)
    real(8) :: dtin

    !     ^^^ virial (short) ^^^
    call calc_virialtensor_inner()
    virialSca=virialTen(1)+virialTen(2)+virialTen(3)

    !     ^^^ piston force (short) ^^^
    pisov=(           +virialSca ) * onethird
    gboxg(1,1)= (pisov-Pinner*cellvol)*rvkbolz
    gboxg(2,2)= gboxg(1,1)
    gboxg(3,3)= gboxg(1,1)

    !     ^^^ update cell velocity (short) ^^^
    vboxg(1,1)=vboxg(1,1)+0.5d0*dtin*gboxg(1,1)*rbaromass
    vboxg(1,2)=0d0
    vboxg(1,3)=0d0
    vboxg(2,1)=0d0
    vboxg(2,2)=vboxg(1,1)
    vboxg(2,3)=0d0
    vboxg(3,1)=0d0
    vboxg(3,2)=0d0
    vboxg(3,3)=vboxg(1,1)
  end subroutine update_vboxg_inner

!----------------------------------------------------------------------
!>
!! \brief  Subroutine for barostat located at outermost
!! \author Yoshimichi Andoh
!<
!! [Ref.] M.E.Tuckerman et al., J.Phys.A:Math.Gen.,39,5629(2006).
  subroutine update_vboxg_outermost_pr(dtin)  ! long
!----------------------------------------------------------------------
    use extended_system
    use thermostat
    use md_multiplestep
    use md_condition
    use md_const
    use unit_cell
    use kinetic_energy
    implicit none
    real(8) :: totke,kinT
    real(8) :: vtemp(6)
    real(8) :: gboxg(3,3)
    real(8) :: dtin,virSca
    real(8) :: wkvirlong(6),k_ene(6),vir(6)
    include 'mpif.h'
    integer(4) :: ierr
    !NI-->
    ! For NtT ensemble
    real(8) :: Trans_box(3,3), tmp_matrix(3,3)
    !<--NI

    !     ^^^ temperature ^^^
    call k_energy_tensor(totke,k_ene)
    kinT = 2.0d0*totke*rvkbolz*degree_of_freedom_inverse

    !     ^^^ virial (long) ^^^
    wkvirlong=virlong
    call mpi_allreduce(wkvirlong,vir,6, &
         &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

    !     ^^^ piston force (long) ^^^
    vtemp(1)  = 2d0*k_ene(1)+vir(1) ! xx
    vtemp(2)  = 2d0*k_ene(2)+vir(2) ! yy
    vtemp(3)  = 2d0*k_ene(3)+vir(3) ! zz
    vtemp(4)  = 2d0*k_ene(4)+vir(4) ! xy,yx
    vtemp(5)  = 2d0*k_ene(5)+vir(5) ! xz,zx
    vtemp(6)  = 2d0*k_ene(6)+vir(6) ! yz,zy
    !NI-->
    ! For NtT ensemble
    tmp_matrix  = MATMUL( box, sigmaS )
    Trans_box   = TRANSPOSE( box )
    StressTerm  = MATMUL( tmp_matrix, Trans_box )
    gboxg(1,1)= (vtemp(1)-(Poutermost*cellvol+StressTerm(1,1)))*rvkbolz+kinT
    gboxg(2,2)= (vtemp(2)-(Poutermost*cellvol+StressTerm(2,2)))*rvkbolz+kinT
    gboxg(3,3)= (vtemp(3)-(Poutermost*cellvol+StressTerm(3,3)))*rvkbolz+kinT
    gboxg(1,2)= (vtemp(4)-StressTerm(1,2))*rvkbolz
    gboxg(1,3)= (vtemp(5)-StressTerm(1,3))*rvkbolz
    gboxg(2,3)= (vtemp(6)-StressTerm(2,3))*rvkbolz
!!!!!!!!!!!!!!!!!!!!!!!!
    ! Older version; if input stress is zero (-> StressTerm=0),
    !         the above time-evolution is the same with following one.
    !      gboxg(1,1)= (vtemp(1)-Poutermost*cellvol)*rvkbolz+kinT
    !      gboxg(2,2)= (vtemp(2)-Poutermost*cellvol)*rvkbolz+kinT
    !      gboxg(3,3)= (vtemp(3)-Poutermost*cellvol)*rvkbolz+kinT
    !      gboxg(1,2)= (vtemp(4)                   )*rvkbolz
    !      gboxg(1,3)= (vtemp(5)                   )*rvkbolz
    !      gboxg(2,3)= (vtemp(6)                   )*rvkbolz
!!!!!!!!!!!!!!!!!!!!!!!!
    !<--NI

    !     ^^^ update cell velocity (long) ^^^
    vboxg(1,1)=vboxg(1,1)+0.5d0*dtin*gboxg(1,1)*rbaromass
    vboxg(1,2)=vboxg(1,2)+0.5d0*dtin*gboxg(1,2)*rbaromass
    vboxg(1,3)=vboxg(1,3)+0.5d0*dtin*gboxg(1,3)*rbaromass
    vboxg(2,1)=0d0
    vboxg(2,2)=vboxg(2,2)+0.5d0*dtin*gboxg(2,2)*rbaromass
    vboxg(2,3)=vboxg(2,3)+0.5d0*dtin*gboxg(2,3)*rbaromass
    vboxg(3,1)=0d0
    vboxg(3,2)=0d0
    vboxg(3,3)=vboxg(3,3)+0.5d0*dtin*gboxg(3,3)*rbaromass
  end subroutine update_vboxg_outermost_pr

!----------------------------------------------------------------------
!>
!! \brief  Subroutine for barostat located at outer
!! \author Yoshimichi Andoh
!<
!! [Ref.] M.E.Tuckerman et al., J.Phys.A:Math.Gen.,39,5629(2006).
  subroutine update_vboxg_outer_pr(dtin)  ! Middle
!----------------------------------------------------------------------
    use extended_system
    use thermostat
    use md_multiplestep
    use md_const
    use unit_cell
    implicit none
    real(8) :: vtemp(6)
    real(8) :: gboxg(3,3)
    real(8) :: dtin,virSca
    real(8) :: wkvirmiddle(6),vir(6)
    include 'mpif.h'
    integer(4) :: ierr

    !     ^^^ virial (middle) ^^^
    wkvirmiddle=virmiddle
    call mpi_allreduce(wkvirmiddle,vir,6, &
         &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

    !     ^^^ piston force (middle) ^^^
    vtemp(1)  =             +vir(1) ! xx
    vtemp(2)  =             +vir(2) ! yy
    vtemp(3)  =             +vir(3) ! zz
    vtemp(4)  =             +vir(4) ! xy,yx
    vtemp(5)  =             +vir(5) ! xz,zx
    vtemp(6)  =             +vir(6) ! yz,zy
    gboxg(1,1)= (vtemp(1)-Pouter*cellvol)*rvkbolz
    gboxg(2,2)= (vtemp(2)-Pouter*cellvol)*rvkbolz
    gboxg(3,3)= (vtemp(3)-Pouter*cellvol)*rvkbolz
    gboxg(1,2)= (vtemp(4)               )*rvkbolz
    gboxg(1,3)= (vtemp(5)               )*rvkbolz
    gboxg(2,3)= (vtemp(6)               )*rvkbolz

    !     ^^^ update cell velocity (middle) ^^^
    vboxg(1,1)=vboxg(1,1)+0.5d0*dtin*gboxg(1,1)*rbaromass
    vboxg(1,2)=vboxg(1,2)+0.5d0*dtin*gboxg(1,2)*rbaromass
    vboxg(1,3)=vboxg(1,3)+0.5d0*dtin*gboxg(1,3)*rbaromass
    vboxg(2,1)=0d0
    vboxg(2,2)=vboxg(2,2)+0.5d0*dtin*gboxg(2,2)*rbaromass
    vboxg(2,3)=vboxg(2,3)+0.5d0*dtin*gboxg(2,3)*rbaromass
    vboxg(3,1)=0d0
    vboxg(3,2)=0d0
    vboxg(3,3)=vboxg(3,3)+0.5d0*dtin*gboxg(3,3)*rbaromass
  end subroutine update_vboxg_outer_pr
!----------------------------------------------------------------------
!>
!! \brief  Subroutine for barostat located at inner
!! \author Yoshimichi Andoh
!<
!! [Ref.] M.E.Tuckerman et al., J.Phys.A:Math.Gen.,39,5629(2006).
  subroutine update_vboxg_inner_pr(dtin)  ! Short
!----------------------------------------------------------------------
    use extended_system
    use thermostat
    use md_multiplestep
    use md_const
    use atom_virial
    use unit_cell
    use pressure_mod
    implicit none
    real(8) :: vtemp(6)
    real(8) :: gboxg(3,3)
    real(8) :: dtin

    !     ^^^ virial (short) ^^^
    call calc_virialtensor_inner()

    !     ^^^ piston force (short) ^^^
    vtemp(1)  =         +virialTen(1) ! xx
    vtemp(2)  =         +virialTen(2) ! yy
    vtemp(3)  =         +virialTen(3) ! zz
    vtemp(4)  =         +virialTen(4) ! xy,yx
    vtemp(5)  =         +virialTen(5) ! xz,zx
    vtemp(6)  =         +virialTen(6) ! yz,zy
    gboxg(1,1)= (vtemp(1)-Pinner*cellvol)*rvkbolz
    gboxg(2,2)= (vtemp(2)-Pinner*cellvol)*rvkbolz
    gboxg(3,3)= (vtemp(3)-Pinner*cellvol)*rvkbolz
    gboxg(1,2)= (vtemp(4)               )*rvkbolz
    gboxg(1,3)= (vtemp(5)               )*rvkbolz
    gboxg(2,3)= (vtemp(6)               )*rvkbolz

    !     ^^^ update cell velocity (short) ^^^
    vboxg(1,1)=vboxg(1,1)+0.5d0*dtin*gboxg(1,1)*rbaromass
    vboxg(1,2)=vboxg(1,2)+0.5d0*dtin*gboxg(1,2)*rbaromass
    vboxg(1,3)=vboxg(1,3)+0.5d0*dtin*gboxg(1,3)*rbaromass
    vboxg(2,1)=0d0
    vboxg(2,2)=vboxg(2,2)+0.5d0*dtin*gboxg(2,2)*rbaromass
    vboxg(2,3)=vboxg(2,3)+0.5d0*dtin*gboxg(2,3)*rbaromass
    vboxg(3,1)=0d0
    vboxg(3,2)=0d0
    vboxg(3,3)=vboxg(3,3)+0.5d0*dtin*gboxg(3,3)*rbaromass
  end subroutine update_vboxg_inner_pr

!----------------------------------------------------------------------
!>
!! \brief  Subroutine for barostat located at outermost
!! \author Yoshimichi Andoh
!<
!! [Ref.] M.E.Tuckerman et al., J.Phys.A:Math.Gen.,39,5629(2006).
  subroutine update_vboxg_outermost_z(dtin)  ! long
!----------------------------------------------------------------------
    use extended_system
    use thermostat
    use md_multiplestep
    use md_condition
    use md_const
    use unit_cell
    use ensemble_numbers
    use md_condition
    use kinetic_energy
    implicit none
    real(8) :: totke,kinT
    !     real(8) :: pisov
    real(8) :: pisovxy,pisovz,k_ene(6)
    real(8) :: pisovx,pisovy
    real(8) :: gboxg(3,3)
    real(8) :: dtin,virSca
    real(8) :: wkvirlong(6),vir(6)
    include 'mpif.h'
    integer(4) :: ierr

    !     ^^^ temperature ^^^
    call k_energy_tensor(totke,k_ene)
    kinT = 2.0d0*totke*rvkbolz*degree_of_freedom_inverse

    !     ^^^ virial (long) ^^^
    wkvirlong=virlong
    call mpi_allreduce(wkvirlong,vir,3, &
         &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

    !     ^^^ piston force (long) ^^^
    if(      md_condition__ensemble == NPT_Z )then
       pisovxy=( 2.d0*k_ene(1)+vir(1)+2.d0*k_ene(2)+vir(2) ) / 2d0
       pisovz =  2.d0*k_ene(3)+vir(3)
       gboxg(1,1)= (pisovxy-Poutermost*cellvol)*rvkbolz+kinT
       gboxg(2,2)= gboxg(1,1)
       gboxg(3,3)= (pisovz -Poutermost*cellvol)*rvkbolz+kinT
    elseif(  md_condition__ensemble == NPXPYPZT )then
       pisovx =  2.d0*k_ene(1)+vir(1)
       pisovy =  2.d0*k_ene(2)+vir(2)
       pisovz =  2.d0*k_ene(3)+vir(3)
       gboxg(1,1)= (pisovx -Poutermost*cellvol)*rvkbolz+kinT
       gboxg(2,2)= (pisovy -Poutermost*cellvol)*rvkbolz+kinT
       gboxg(3,3)= (pisovz -Poutermost*cellvol)*rvkbolz+kinT
    elseif(  md_condition__ensemble == NPTLZ )then
       pisovxy=( 2.d0*k_ene(1)+vir(1)+2.d0*k_ene(2)+vir(2) ) / 2d0
       gboxg(1,1)= (pisovxy-Poutermost*cellvol)*rvkbolz+kinT
       gboxg(2,2)= gboxg(1,1)
       gboxg(3,3)= 0d0 ; vboxg(3,3)=0d0
    elseif(  md_condition__ensemble == NLXLYPZT )then
       gboxg(1,1)= 0d0 ; vboxg(1,1)=0d0
       gboxg(2,2)= 0d0 ; vboxg(2,2)=0d0
       pisovz =  2.d0*k_ene(3)+vir(3)
       gboxg(3,3)= (pisovz -Poutermost*cellvol)*rvkbolz+kinT
    elseif(  md_condition__ensemble == NLXLYLZT )then
       gboxg(1,1)= 0d0 ; vboxg(1,1)=0d0
       gboxg(2,2)= 0d0 ; vboxg(2,2)=0d0
       gboxg(3,3)= 0d0 ; vboxg(3,3)=0d0
    elseif(  md_condition__ensemble == NPTLZxy )then
       pisovx =( 2.d0*k_ene(1)+vir(1) )
       pisovy =( 2.d0*k_ene(2)+vir(2) )
       gboxg(1,1)= (pisovx -Poutermost*cellvol)*rvkbolz+kinT
       gboxg(2,2)= (pisovy -Poutermost*cellvol)*rvkbolz+kinT
       gboxg(3,3)= 0d0 ; vboxg(3,3)=0d0
    endif

    !     ^^^ update cell velocity (long) ^^^
    vboxg(1,1)=vboxg(1,1)+0.5d0*dtin*gboxg(1,1)*rbaromass
    vboxg(1,2)=0d0
    vboxg(1,3)=0d0
    vboxg(2,1)=0d0
    if( md_condition__ensemble == NPTLZxy .or. &
      & md_condition__ensemble == NPXPYPZT ) then
      vboxg(2,2)=vboxg(2,2)+0.5d0*dtin*gboxg(2,2)*rbaromass
    else
      vboxg(2,2)=vboxg(1,1)
    endif
    vboxg(2,3)=0d0
    vboxg(3,1)=0d0
    vboxg(3,2)=0d0
    vboxg(3,3)=vboxg(3,3)+0.5d0*dtin*gboxg(3,3)*rbaromass
  end subroutine update_vboxg_outermost_z

!----------------------------------------------------------------------
!>
!! \brief  Subroutine for barostat located at outer
!! \author Yoshimichi Andoh
!<
!! [Ref.] M.E.Tuckerman et al., J.Phys.A:Math.Gen.,39,5629(2006).
  subroutine update_vboxg_outer_z(dtin)  ! Middle
!----------------------------------------------------------------------
    use extended_system
    use thermostat
    use md_multiplestep
    use md_const
    use unit_cell
    use ensemble_numbers
    use md_condition
    implicit none
    real(8) :: pisovxy,pisovz
    real(8) :: pisovx,pisovy
    real(8) :: gboxg(3,3)
    real(8) :: dtin,virSca
    real(8) :: wkvirmiddle(6),vir(6)
    include 'mpif.h'
    integer(4) :: ierr

    !     ^^^ virial (middle) ^^^
    wkvirmiddle=virmiddle
    call mpi_allreduce(wkvirmiddle,vir,3, &
         &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

    !     ^^^ piston force (middle) ^^^
    if(      md_condition__ensemble == NPT_Z )then
       pisovxy = (+vir(1)+vir(2) ) / 2d0
       pisovz =               +vir(3)
       gboxg(1,1)= (pisovxy-Pouter*cellvol)*rvkbolz
       gboxg(2,2)= gboxg(1,1)
       gboxg(3,3)= (pisovz -Pouter*cellvol)*rvkbolz
    elseif(  md_condition__ensemble == NPXPYPZT )then
       pisovx =               +vir(1)
       pisovy =               +vir(2) 
       pisovz =               +vir(3)
       gboxg(1,1)= (pisovx -Pouter*cellvol)*rvkbolz
       gboxg(2,2)= (pisovy -Pouter*cellvol)*rvkbolz
       gboxg(3,3)= (pisovz -Pouter*cellvol)*rvkbolz
    elseif(  md_condition__ensemble == NPTLZ )then
       pisovxy=(+vir(1)+vir(2) ) / 2d0
       gboxg(1,1)= (pisovxy-Pouter*cellvol)*rvkbolz
       gboxg(2,2)= gboxg(1,1)
       gboxg(3,3)= 0d0 ; vboxg(3,3)=0d0
    elseif(  md_condition__ensemble == NLXLYPZT )then
       pisovz =               +vir(3)
       gboxg(1,1)= 0d0 ; vboxg(1,1)=0d0 !(pisovxy-Pouter*cellvol)*rvkbolz
       gboxg(2,2)= 0d0 ; vboxg(2,2)=0d0 !gboxg(1,1)
       gboxg(3,3)= (pisovz -Pouter*cellvol)*rvkbolz
    elseif(  md_condition__ensemble == NLXLYLZT )then
       gboxg(1,1)= 0d0 ; vboxg(1,1)=0d0 !(pisovxy-Pouter*cellvol)*rvkbolz
       gboxg(2,2)= 0d0 ; vboxg(2,2)=0d0 !gboxg(1,1)
       gboxg(3,3)= 0d0 ; vboxg(3,3)=0d0
    elseif(  md_condition__ensemble == NPTLZxy )then
       pisovx =(              +vir(1) )
       pisovy =(              +vir(2) )
       gboxg(1,1)= (pisovx -Pouter*cellvol)*rvkbolz
       gboxg(2,2)= (pisovy -Pouter*cellvol)*rvkbolz
       gboxg(3,3)= 0d0 ; vboxg(3,3)=0d0
    endif

    !     ^^^ update cell velocity (middle) ^^^
    vboxg(1,1)=vboxg(1,1)+0.5d0*dtin*gboxg(1,1)*rbaromass
    vboxg(1,2)=0d0
    vboxg(1,3)=0d0
    vboxg(2,1)=0d0
    if( md_condition__ensemble == NPTLZxy .or. &
      & md_condition__ensemble == NPXPYPZT ) then
      vboxg(2,2)=vboxg(2,2)+0.5d0*dtin*gboxg(2,2)*rbaromass
    else
      vboxg(2,2)=vboxg(1,1)
    endif
    vboxg(2,3)=0d0
    vboxg(3,1)=0d0
    vboxg(3,2)=0d0
    vboxg(3,3)=vboxg(3,3)+0.5d0*dtin*gboxg(3,3)*rbaromass

  end subroutine update_vboxg_outer_z
!----------------------------------------------------------------------
!>
!! \brief  Subroutine for barostat located at inner
!! \author Yoshimichi Andoh
!<
!! [Ref.] M.E.Tuckerman et al., J.Phys.A:Math.Gen.,39,5629(2006).
  subroutine update_vboxg_inner_z(dtin)  ! Short
!----------------------------------------------------------------------
    use extended_system
    use thermostat
    use md_multiplestep
    use md_const
    use atom_virial
    use unit_cell
    use md_condition
    use ensemble_numbers
    use pressure_mod
    implicit none
    real(8) :: pisovxy,pisovz
    real(8) :: pisovx,pisovy
    real(8) :: gboxg(3,3)
    real(8) :: dtin

    !     ^^^ virial (short) ^^^
    call calc_virialtensor_inner()

    !     ^^^ piston force (short) ^^^
    if(      md_condition__ensemble == NPT_Z )then
       pisovxy=( virialTen(1)   +virialTen(2) ) /2d0
       pisovz =             virialTen(3)
       gboxg(1,1)= (pisovxy-Pinner*cellvol)*rvkbolz
       gboxg(2,2)= gboxg(1,1)
       gboxg(3,3)= (pisovz -Pinner*cellvol)*rvkbolz
    elseif(  md_condition__ensemble == NPXPYPZT )then
       pisovx =             virialTen(1)
       pisovy =             virialTen(2)
       pisovz =             virialTen(3)
       gboxg(1,1)= (pisovx -Pinner*cellvol)*rvkbolz
       gboxg(2,2)= (pisovy -Pinner*cellvol)*rvkbolz
       gboxg(3,3)= (pisovz -Pinner*cellvol)*rvkbolz
    elseif(  md_condition__ensemble == NPTLZ )then
       pisovxy=( virialTen(1)+virialTen(2) ) /2d0
       gboxg(1,1)= (pisovxy-Pinner*cellvol)*rvkbolz
       gboxg(2,2)= gboxg(1,1)
       gboxg(3,3)= 0d0 ; vboxg(3,3)=0d0
    elseif(  md_condition__ensemble == NLXLYPZT )then
       pisovz =             virialTen(3)
       gboxg(1,1)= 0d0 ; vboxg(1,1)=0d0 !(pisovxy-Pinner*cellvol)*rvkbolz
       gboxg(2,2)= 0d0 ; vboxg(2,2)=0d0 !gboxg(1,1)
       gboxg(3,3)= (pisovz -Pinner*cellvol)*rvkbolz
    elseif(  md_condition__ensemble == NLXLYLZT )then
       gboxg(1,1)= 0d0 ; vboxg(1,1)=0d0 !(pisovxy-Pinner*cellvol)*rvkbolz
       gboxg(2,2)= 0d0 ; vboxg(2,2)=0d0 !gboxg(1,1)
       gboxg(3,3)= 0d0 ; vboxg(3,3)=0d0
    elseif(  md_condition__ensemble == NPTLZxy )then
       pisovx =(            virialTen(1) )
       pisovy =(            virialTen(2) )
       gboxg(1,1)= (pisovx -Pinner*cellvol)*rvkbolz
       gboxg(2,2)= (pisovy -Pinner*cellvol)*rvkbolz
       gboxg(3,3)= 0d0 ; vboxg(3,3)=0d0
   endif

   !     ^^^ update cell velocity (short) ^^^
   vboxg(1,1)=vboxg(1,1)+0.5d0*dtin*gboxg(1,1)*rbaromass
   vboxg(1,2)=0d0
   vboxg(1,3)=0d0
   vboxg(2,1)=0d0
   if( md_condition__ensemble == NPTLZxy .or. &
     & md_condition__ensemble == NPXPYPZT ) then
     vboxg(2,2)=vboxg(2,2)+0.5d0*dtin*gboxg(2,2)*rbaromass
   else
     vboxg(2,2)=vboxg(1,1)
   endif
   vboxg(2,3)=0d0
   vboxg(3,1)=0d0
   vboxg(3,2)=0d0
   vboxg(3,3)=vboxg(3,3)+0.5d0*dtin*gboxg(3,3)*rbaromass
 end subroutine update_vboxg_inner_z

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to enter zero into variables/arrays
!! \author Yoshimichi Andoh
!<
  subroutine zero_barostat
!----------------------------------------------------------------------
    use thermostat
    implicit none

    vssb=0d0
    rssb=0d0

    vboxg=0d0
  end subroutine zero_barostat

end module barostat

