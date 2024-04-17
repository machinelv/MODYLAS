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
!!         (EoM) for NVT ensemble.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module for NVE ensemble (isotropic).
!! \author Noriyuki Yoshii, Kensuke Iwahashi
!<
module nvt_mod
  use omp_lib
  use mpi_tool
  use pressure_mod
  implicit none

contains

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to integrate EoM for NVT ensemble
!! \author Yoshimichi Andoh
!<
  subroutine nvt_integrate()
!----------------------------------------------------------------------
    use atom_mass
    use thermostat
    use shake_rattle_roll
    use atom_virial
    use md_condition
    use md_const
    use forces
    use md_monitors
    use trajectory_mpi
    use param
    use md_multiplestep
    use ensemble_numbers
    use bond_breaking_event, only : detect_bond_breaking
!default: MTD
    use comm_direct2_dr, only : comm_direct_2_dr
    use comm_direct3_dr, only : comm_direct_3_dr
    use comm_bound_mod
    use subcell
    use boundary
    use unit_cell
    use kinetic_energy
    use thermostat
    use force_wrap, only : md_calculate_forces
#ifdef SEGSHAKE
    use center_of_mass
#endif
#ifdef TIP4
    use tip4p, only : update_msite_coordinate, &
   &                  distribute_force_on_msite
#endif
    implicit none
    integer(4) :: i0,ipar,k0,l0,l1,iconstraints
    integer(4) :: iam

    dthL=dth/maxMTm/maxMTl
    dtL =dt /maxMTm/maxMTl

    iMT=0

    !     ### NHC (XO) ###
    call nhchain()  ! always XO

    DO MTl=1,maxMTl   !! == Long-range force ==
       DO MTm=1,maxMTm   !! == Middle-range force ==

          iMT=iMT+1

#ifdef SEGSHAKE
          CenterOfMassAOld = CenterOfMassA
          CenterOfMassBOld = CenterOfMassB
#endif
          if(totnconst > 0) then
             call init_wkdvir(iMT)
!$omp parallel default(shared) &
!$omp& private(k0,i0,iam) &
!$omp& private(l0,l1,iconstraints)
!$      iam=omp_get_thread_num()
!$omp do
             do k0=1,nselfseg
                do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
                   xyzstr(:,i0) = wkxyz(:,i0)
                   wk_dfc(:,i0)=0.d0
                enddo ! i0
             enddo ! k0
!$omp end do
             wkdvir2(:,iam)=0d0
!$omp do
             do l0=1,l0max
                iconstraints=ibseL(l0)
                do l1=1,iconstraints
                   lambL(l1,l0) = 0d0
                enddo
             enddo
!$omp end do
!$omp end parallel
          endif

!$omp parallel default(shared) &
!$omp& private(k0,i0,ipar)
!$omp do
          do k0=1,nselfseg
             do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
                ipar=paranum(m2i(i0))
                wkv(1,i0)=wkv(1,i0)+wk_f(1,i0)*dthL*r_mass(ipar)
                wkv(2,i0)=wkv(2,i0)+wk_f(2,i0)*dthL*r_mass(ipar)
                wkv(3,i0)=wkv(3,i0)+wk_f(3,i0)*dthL*r_mass(ipar)
                wkxyz(1,i0) = wkxyz(1,i0) + wkv(1,i0)*dtL
                wkxyz(2,i0) = wkxyz(2,i0) + wkv(2,i0)*dtL
                wkxyz(3,i0) = wkxyz(3,i0) + wkv(3,i0)*dtL
             enddo ! i0
          enddo ! k0
!$omp end do
!$omp end parallel

          call update_wkvir()
          if (totnconst .gt. 0) then
             call shake_roll(dtL)
             call update_wkdvir()
          endif
#ifdef TIP4
          call update_msite_coordinate
#endif

          if(allow_bond_breaking>0.0d0)then
              call detect_bond_breaking()
          endif

          if(MTl==maxMTl.and.MTm==maxMTm)then
             call update_wsegc()
             call comm_bound()
             if(totnconst > 0) call update_shake_local()
#ifdef SEGSHAKE
             call list_own_i0  !! must be befor comm_dir3
#endif
          ENDIF

#ifdef SEGSHAKE
          call  calc_center_of_mass
          call  shake_com(dtL)
#endif

!default: MTD
          call comm_direct_3_dr(wkxyz)
          call apply_pbc()

          call md_calculate_forces()  ! output wk_f

          iMT=iMT+1

          if(totnconst > 0) then
             call init_wkdvir(iMT)
!$omp parallel default(shared) &
!$omp& private(k0,i0,iam) &
!$omp& private(l0,l1,iconstraints)
!$      iam=omp_get_thread_num()
!$omp do
             do k0=1,nselfseg
                do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
                   wk_dfc(1,i0)=0.d0
                   wk_dfc(2,i0)=0.d0
                   wk_dfc(3,i0)=0.d0
                enddo ! i0
             enddo ! k0
!$omp end do
             wkdvir2(:,iam)=0d0
!$omp do
             do l0=1,l0max
                iconstraints=ibseL(l0)
                do l1=1,iconstraints
                   lambL(l1,l0) = 0d0
                enddo
             enddo
!$omp end do
!$omp end parallel
          endif

!$omp parallel default(shared) &
!$omp& private(k0,i0,ipar)
!$omp do
          do k0=1,nselfseg
             do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
                ipar=paranum(m2i(i0))
                wkv(1,i0) = wkv(1,i0) + wk_f(1,i0)*dthL*r_mass(ipar)
                wkv(2,i0) = wkv(2,i0) + wk_f(2,i0)*dthL*r_mass(ipar)
                wkv(3,i0) = wkv(3,i0) + wk_f(3,i0)*dthL*r_mass(ipar)
             enddo ! i0
          enddo ! k0
!$omp end do
!$omp end parallel

          call update_wkvir()
          if (totnconst .gt. 0) then
             call rattle_roll(dtL)
             call update_wkdvir()
          endif

#ifdef SEGSHAKE
          call calc_center_of_mass
          call calc_center_of_velocity()
          call rattle_com(dtL)
#endif

       ENDDO  !  MT middle
    ENDDO  !  MT long

    !     ### NHC (XO) ###
    call nhchain()  !! always XO

    call remove_system_momentum_para
#ifdef SEGSHAKE
    call calc_MeanForce
#endif

    if(md_condition__ensemble == NVLXLYLZT) then
        cellx = cellx + deformx
        celly = celly + deformy
        cellz = cellz + deformz
        cellxh = cellx * 0.5
        cellyh = celly * 0.5
        cellzh = cellz * 0.5
    endif

  end subroutine nvt_integrate
!----------------------------------------------------------------------
!
!     calculate thermodynamic values
!
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate outputs (NVT ensemble)
!! \author Yoshimichi Andoh
!<
  subroutine calc_hamiltonian_nvt()
!----------------------------------------------------------------------
    use thermostat
    use atom_virial
    use md_const
    use md_monitors
    use lj_mod
    use md_condition
    use unit_cell
    use kinetic_energy
    implicit none
    real(8) :: totke,k_ene(6)
    real(8) :: Uc
    include 'mpif.h'
    integer(4) :: ierr

    !     ^^^ reduce potential energy ^^^
    if(nprocs.eq.1) then
       p_energy = wk_p_energy
    else
       p_energy=0d0
       call mpi_allreduce(wk_p_energy,p_energy,1, &
            &       mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    end if
    if(LJ_LRC)then  ! add LJ-LRC to potential
       Uc = corrector_potential / cellvol
       p_energy = p_energy + Uc
    endif

    !     ^^^ reduce virial ^^^
    call calc_virialtensor() ! with correction

    !     ^^^ reduce kinetic energy ^^^
    call k_energy_tensor(totke,k_ene)

    !     ^^^ temperature ^^^
    k_energy    = totke
    temperature = 2.0d0*totke*degree_of_freedom_inverse    *rvkbolz
    t_energy    = p_energy + k_energy
    hamiltonian = t_energy + st
    pressure=(2.0d0*totke  +virialSca)/(3d0*cellvol)
    ptensor =(2.0d0*k_ene+virialTen)/     cellvol

    sum_P=sum_P+pressure

    !     ^^^ pressure division ^^^
    P_part(1)=(           +vir_part(1))/(3d0*cellvol)
    P_part(2)=(           +vir_part(2))/(3d0*cellvol)
    P_part(3)=(2.0d0*totke  +vir_part(3))/(3d0*cellvol)
    sum_P_part(:)=sum_P_part(:)+P_part(:)

  end subroutine calc_hamiltonian_nvt

end module nvt_mod
