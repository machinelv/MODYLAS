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
!!         (EoM) for semi-anisotropic NPT ensemble.
!>
!----------------------------------------------------------------------
!>
!! \brief  Module for NPT ensemble (semi-anisotropic).
!! \author Yoshimichi Andoh, Noriyuki Yoshii, Kensuke Iwahashi
!<
!----------------------------------------------------------------------
module npt_z_mod
  use omp_lib
  use mpi_tool
  use pressure_mod
  implicit none

contains

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to integrate EoM for npt ensemble (semi-anisotropic)
!! \author Yoshimichi Andoh
!<
  subroutine npt_z_integrate()
!----------------------------------------------------------------------
    use atom_mass
    use shake_rattle_roll
    use atom_virial
    use md_condition
    use md_const
    use md_forces
    use md_monitors
    use trajectory_mpi
    use param
    use md_multiplestep
    use bond_breaking_event, only : detect_bond_breaking
!default: MTD
    use comm_direct2_dr, only : comm_direct_2_dr
    use comm_direct3_dr, only : comm_direct_3_dr
    use npt_a_mod, only : recover_fshort,recover_fmiddle
    use comm_bound_mod
    use subcell
    use boundary
    use unit_cell
    use kinetic_energy
    use thermostat
    use cell_shape
    use update
    use barostat
    use force_short
    use force_middle
    use force_long
#ifdef SEGSHAKE
    use center_of_mass
#endif
#ifdef TIP4
    use tip4p, only : update_msite_coordinate,  &
   &                  distribute_force_on_msite
#endif
    implicit none
    integer(4) :: i,j
    integer(4) :: i0,k0,l0,l1,iconstraints
    real(8) :: ubox(3,3),aa2(3)
    real(8) :: veigv(3,3),vtemp(3,3)
    real(8) :: boxstr(3,3),vboxgstr(3,3)
    integer(4) :: liter,maxliter
    include 'mpif.h'
    integer(4) :: iam

    dthL=dth/maxMTm/maxMTl
    dtL =dt /maxMTm/maxMTl

    call cell_convert1(cellx,celly,cellz,alpha,beta,gamma,box)

    iMT=0

    if(totnconst > 0)then
       maxliter=1000
    endif

!default: MTD
    if(mdstep==0)then
      call recover_fshort()
      call recover_fmiddle()
    endif

    !     ### NHCA (XO) ###
    call nhc_z_thermo()
    call update_vboxg_outermost_z(dtL*maxMTm*maxMTl)
    call update_velocities_outermost(dtL*maxMTm*maxMTl)

    DO MTl=1,maxMTl   !! == Long-range force ==
       DO MTm=1,maxMTm   !! == Middle-range force ==

          iMT=iMT+1

#ifdef SEGSHAKE
          CenterOfMassAOld = CenterOfMassA
          CenterOfMassBOld = CenterOfMassB
#endif
          if(totnconst > 0) then
             call init_wkdvir(iMT)
             liter=0
!$omp parallel default(shared) &
!$omp& private(k0,i0) &
!$omp& private(l0,l1,iconstraints)
!$omp do
             do l0=1,l0max
                iconstraints=ibseL(l0)
                do l1=1,iconstraints
                   lambL(l1,l0) = 0d0
                enddo ! l1
             enddo ! l0
!$omp end do
!$omp do
             do k0=1,nselfseg
                do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
                   wk_dfc(1,i0)=0.d0
                   wk_dfc(2,i0)=0.d0
                   wk_dfc(3,i0)=0.d0
                enddo ! i0
             enddo ! k0
!$omp end do
             !
             !       store values
             !
!$omp do
             do k0=1,nselfseg
                do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
                   xyzstr(:,i0) = wkxyz(:,i0)
                   vstr(1:3,i0) = wkv(1:3,i0)
                enddo ! i0
             enddo ! k0
!$omp end do
!$omp end parallel
             boxstr   = box
             vboxgstr = vboxg
             rssstr   = rss
             vssstr   = vss
             rssbstr  = rssb
             vssbstr  = vssb
          endif
10        continue
          if(totnconst > 0) then
             if(liter.gt.maxliter) then
                write(0,'(a)')  'ERROR: Too many iterations at shake_roll'
                call modylas_abort()
             endif
             liter=liter+1
             wk_roll_flg=0d0; roll_flg=0d0
!$omp parallel default(shared) &
!$omp& private(k0,i0,iam) &
!$omp& private(l0,l1,iconstraints)
!$      iam=omp_get_thread_num()
             wkdvir2(:,iam)=0d0
!$omp do
             do l0=1,l0max
                iconstraints=ibseL(l0)
                do l1=1,iconstraints
                   lambLstr(l1,l0) = lambL(l1,l0)
                   lambL(l1,l0) = 0d0
                enddo
             enddo
!$omp end do
             !
             !       restore values
             !
!$omp do
             do k0=1,nselfseg
                do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
                   wkxyz(:,i0) = xyzstr(:,i0)
                   wkv(1:3,i0) = vstr(1:3,i0)
                enddo ! i0
             enddo ! k0
!$omp end do
!$omp end parallel
             box   = boxstr
             vboxg = vboxgstr
             rss   = rssstr
             vss   = vssstr
             rssb  = rssbstr
             vssb  = vssbstr

             call cell_convert2(cellx,celly,cellz,alpha,beta,gamma,box)
             call cell_convert3(cellx,celly,cellz,alpha,beta,gamma,cellxh,cellyh, &
                  &  cellzh,cellvol,sinbeta,cosbeta, &
                  &  singamma,cosgamma,b_factor,g_factor)
          endif

          if(scaleM==1d0)then
             call update_vboxg_outer_z(dtL*maxMTm)
             call update_velocities_outer(dtL*maxMTm)
          endif

          call update_vboxg_inner_z(dtL)
          call update_velocities_inner()

          call update_coordinates(vtemp,veigv,aa2)

          ubox = matmul(vtemp,box)
          do i=1,3
             do j=1,3
                ubox(i,j)=ubox(i,j)*aa2(i)
             enddo
          enddo
          box = matmul(veigv,ubox)
          call cell_convert2(cellx,celly,cellz,alpha,beta,gamma,box)
          call cell_convert3(cellx,celly,cellz,alpha,beta,gamma,cellxh,cellyh, &
               &  cellzh,cellvol,sinbeta,cosbeta, &
               &  singamma,cosgamma,b_factor,g_factor)

          call update_wkvir()
          if (totnconst .gt. 0) then
             call shake_roll(dtL)
             call update_wkdvir()
             call check_roll_convergence()
             if(roll_flg >= 1d0) goto 10
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
          !!    ^^^ short ^^^
          call md_calculate_forces_short(fshort,virshort,eneshort)
          call recover_fshort()
          !!    ^^^ middle ^^^
          if(MTm==maxMTm)then
             scaleM=1d0
             call md_calculate_forces_middle(fmiddle,virmiddle,enemiddle)
!default: MTD
             call recover_fmiddle()
#ifdef TIP4
             call distribute_force_on_msite(fmiddle)
#endif
          else
             scaleM=0d0
          endif
          !!    ^^^ long ^^^
          if(MTl==maxMTl.and.MTm==maxMTm)then
             scaleL=1d0
             call md_calculate_forces_long(flong,virlong,enelong)
#ifdef TIP4
             call distribute_force_on_msite(flong)
#endif
          else
             scaleL=0d0
          endif

          !!    ^^^ sum short, middle, and long ^^^
          wk_p_energy=eneshort+enemiddle+enelong

          iMT=iMT+1

          if(totnconst > 0) then
             call init_wkdvir(iMT)
             liter=0
!$omp parallel default(shared) &
!$omp& private(k0,i0) &
!$omp& private(l0,l1,iconstraints)
!$omp do
             do l0=1,l0max
                iconstraints=ibseL(l0)
                do l1=1,iconstraints
                   lambL(l1,l0) = 0d0
                enddo ! l1
             enddo ! l0
!$omp end do
!$omp do
             do k0=1,nselfseg
                do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
                   wk_dfc(1,i0)=0.d0
                   wk_dfc(2,i0)=0.d0
                   wk_dfc(3,i0)=0.d0
                enddo ! i0
             enddo ! k0
!$omp end do
             !
             !       store values
             !
!$omp do
             do k0=1,nselfseg
                do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
                   vstr(1:3,i0) = wkv(1:3,i0)
                enddo ! i0
             enddo ! k0
!$omp end do
!$omp end parallel
             vboxgstr = vboxg
             rssstr   = rss
             vssstr   = vss
             rssbstr  = rssb
             vssbstr  = vssb
          endif
20        continue
          if(totnconst > 0) then
             if(liter.gt.maxliter) then
                write(0,'(a)')  'ERROR: Too many iterations at rattle_roll'
                call modylas_abort()
             endif
             liter=liter+1
             wk_roll_flg=0d0; roll_flg=0d0
!$omp parallel default(shared) &
!$omp& private(k0,i0,iam) &
!$omp& private(l0,l1,iconstraints)
!$      iam=omp_get_thread_num()
             wkdvir2(:,iam)=0d0
!$omp do
             do l0=1,l0max
                iconstraints=ibseL(l0)
                do l1=1,iconstraints
                   lambLstr(l1,l0) = lambL(l1,l0)
                   lambL(l1,l0) = 0d0
                enddo
             enddo
!$omp end do
             !
             !       restore values
             !
!$omp do
             do k0=1,nselfseg
                do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
                   wkv(1:3,i0) = vstr(1:3,i0)
                enddo ! i0
             enddo ! k0
!$omp end do
!$omp end parallel
             vboxg = vboxgstr
             rss   = rssstr
             vss   = vssstr
             rssb  = rssbstr
             vssb  = vssbstr
          endif

          call update_velocities_inner()
          call update_vboxg_inner_z(dtL)

          if(scaleM==1d0)then
             call update_velocities_outer(dtL*maxMTm)
             call update_vboxg_outer_z(dtL*maxMTm)
          endif

          if(scaleM*scaleL==1d0)then
             call update_velocities_outermost(dtL*maxMTm*maxMTl)
             call update_vboxg_outermost_z(dtL*maxMTm*maxMTl)
             call nhc_a_thermo()
          endif

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
          if (totnconst .gt. 0) then
             call check_roll_convergence()
             if(roll_flg >= 1d0) goto 20
          endif

       ENDDO  !  MT middle
    ENDDO  !  MT long

    call remove_system_momentum_para
#ifdef SEGSHAKE
    call calc_MeanForce
#endif

    cellx = cellx + deformx
    celly = celly + deformy
    cellz = cellz + deformz
    cellxh = cellx * 0.5
    cellyh = celly * 0.5
    cellzh = cellz * 0.5

  end subroutine npt_z_integrate
!----------------------------------------------------------------------
!
!     calculate thermodynamic values
!
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate outputs (NPT ensemble)
!! \author Yoshimichi Andoh
!<
  subroutine calc_hamiltonian_npt_z()
!----------------------------------------------------------------------
    use atom_virial
    use md_const
    use md_monitors
    use lj_mod
    use md_condition
    use unit_cell
    use kinetic_energy
    use thermostat
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
    temperature = 2.0d0*totke*degree_of_freedom_inverse*rvkbolz
    t_energy    = p_energy + k_energy
    hamiltonian = t_energy + st
    pressure=(2.0d0*totke+virialSca)/(3d0*cellvol)
    ptensor =(2.0d0*k_ene+virialTen)/     cellvol

    sum_P=sum_P+pressure

    !     ^^^ pressure division ^^^
    P_part(1)=(           +vir_part(1))/(3d0*cellvol)
    P_part(2)=(           +vir_part(2))/(3d0*cellvol)
    P_part(3)=(2.0d0*totke+vir_part(3))/(3d0*cellvol)
    sum_P_part(:)=sum_P_part(:)+P_part(:)
  end subroutine calc_hamiltonian_npt_z

end module npt_z_mod
