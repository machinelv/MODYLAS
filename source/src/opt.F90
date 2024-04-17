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
!! \brief  Module and subroutines to perform optimization of structure.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to perform structure optimization
!! \author  Atsushi Yamada
!<
!-------------------------------------------------------------------------
module opt_mod
  use omp_lib
  use mpi_tool
  implicit none

  real(8),allocatable :: opt_direction(:,:)
  real(8) :: step_length=0.2d-10
  real(8) :: up_step_length=1.2d0, down_step_length=0.6d0
  real(8) :: conv=1.0d-13

contains
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to set variables of /optimize/step_length.
!! \author  Atsushi Yamada
!<
  subroutine fmod_opt_condition__step_length(value)
    implicit none
    real(8), intent(in) :: value
    step_length = value
  end subroutine fmod_opt_condition__step_length
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to set variables of /optimize/up_rate.
!! \author  Atsushi Yamada
!<
  subroutine fmod_opt_condition__step_lurate(value)
    implicit none
    real(8), intent(in) :: value
    up_step_length = value
  end subroutine fmod_opt_condition__step_lurate
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to set variables of /optimize/down_rate.
!! \author  Atsushi Yamada
!<
  subroutine fmod_opt_condition__step_ldrate(value)
    implicit none
    real(8), intent(in) :: value
    down_step_length = value
  end subroutine fmod_opt_condition__step_ldrate
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to set variables of /optimize/convergence.
!! \author  Atsushi Yamada
!<
  subroutine fmod_opt_condition__convergence(value)
    implicit none
    real(8), intent(in) :: value
    conv = value
  end subroutine fmod_opt_condition__convergence
!-------------------------------------------------------------------------
!
!     steepest descent method
!
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to perform structure optimization by
!!         the steepest descent method.
!! \author  Atsushi Yamada
!<
  subroutine opt_integrate()
    use md_multiplestep
    use forces
    use md_monitors
    use md_condition
!default: MTD
    use comm_direct2_dr, only : comm_direct_2_dr
    use comm_direct3_dr, only : comm_direct_3_dr
    use comm_bound_mod
    use position_constrain
    use shake_rattle_roll
    use trajectory_mpi
    use atom_virial
    use subcell
    use boundary
    use force_wrap, only : md_calculate_forces
    use atom_mass
    use param
#ifdef SEGSHAKE
    use center_of_mass
#endif
#ifdef TIP4
    use tip4p, only : update_msite_coordinate, &
   &                  distribute_force_on_msite
#endif
    implicit none
    integer(4) :: i,i0,k0, iam,l0,l1,iconstraints
    real(8) :: mean_force

    call update_wsegc()
    call comm_bound()
    if(totnconst > 0) call update_shake_local()
#ifdef SEGSHAKE
    call list_own_i0  !! must be befor comm_dir3
    call calc_center_of_mass
    CenterOfMassAOld = CenterOfMassA
    CenterOfMassBOld = CenterOfMassB
#endif
!default: MTD
    call comm_direct_3_dr(wkxyz)
    call apply_pbc()

    maxMTl=1    ! opt always assumes single time step(STS)
    maxMTm=1    ! opt always assumes single time step(STS)
    MTl=maxMTl  ! dummy to realize STS
    MTm=maxMTm  ! dummy to realize STS
    call md_calculate_forces()  ! output wk_f

    ! -- calculate summation of root mean square of force --
    call opt_calculate_mean_force(mean_force)

    ! --  analyze --sub
    !     if(opt_analyze_flag) then
    !       ! -- calculate maximum force --
    !       call calculate_max_force(i_max_force, max_force)

    !       ! -- calculate minimum distance --
    !       call calculate_min_distance(iatom1, iatom2, min_distance)

    !       ! -- output analyzed informations --
    !       call record_opt_analyze_information(mean_force, max_force, &
    !    &                   i_max_force, min_distance, iatom1, iatom2 )
    !     endif


    ! -- check the convergence condition --sub
    !     if(mean_force .lt. conv) then
    !#ifdef MPIPARA
    !       if(myrank.eq.0) then
    !#endif
    !          write(f_run,9500)
    !          write(f_run,9510)
    !#ifdef MPIPARA
    !       endif
    !#endif
    !       (output optimized coordinate)
    !        itmp1 = trj_interval
    !        itmp2 = restart_interval
    !        trj_interval = 1
    !        restart_interval = 1
    !        write(f_trj,9600)  ! for ascii mode
    !        call record_current_trajectory
    !        call record_restart
    !        if (ascii_output)  call md_ascii_output ! output *.mdxyz (ascii)
    !        trj_interval = itmp1
    !        restart_interval = itmp2
    !        stop
    !     endif

    ! -- calculate direction normalized vector --
    call opt_calc_direction_vector(mean_force)

    ! -- shake --
    if(totnconst > 0) then
!$omp parallel default(shared) &
!$omp& private(k0,i0,iam) &
!$omp& private(l0,l1,iconstraints)
!$      iam=omp_get_thread_num()
!$omp do
       do k0=1,nselfseg
          do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
             xyzstr(:,i0) = wkxyz(:,i0)
             wk_dfc(1,i0)=0.d0
             wk_dfc(2,i0)=0.d0
             wk_dfc(3,i0)=0.d0
          enddo ! i0
       enddo ! k0
!$omp end do
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

    ! -- 1 step optimization --
!$omp parallel do default(shared) &
!$omp& private(k0,i0)
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
#ifdef TIP4
    if(mass(paranum(m2i(i0))) .lt. 0) cycle
#endif
          wkxyz(1,i0)=wkxyz(1,i0)+opt_direction(1,i0)*step_length
          wkxyz(2,i0)=wkxyz(2,i0)+opt_direction(2,i0)*step_length
          wkxyz(3,i0)=wkxyz(3,i0)+opt_direction(3,i0)*step_length
       enddo
    enddo
#ifdef TIP4
    call update_msite_coordinate
#endif

    ! -- shake --
    if (totnconst .gt. 0) then
       dt = 1.0d-15  ![fs]
       call shake_roll(dt)

       if(type_p_constrain==3)then  ! for fix option
!$omp parallel do default(shared) &
!$omp& private(i,k0,i0) &
!$omp& shared(one_or_zero_fix_atom,m2i)
          do k0=1,nselfseg
             do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
                i = m2i(i0)
#ifdef TIP4
    if(mass(paranum(i)) .lt. 0) cycle
#endif
                wkxyz(1,i0) = xyzstr(1,i0)*(1-one_or_zero_fix_atom(i)) &
                     &                   +      wkxyz(1,i0)*one_or_zero_fix_atom(i)
                wkxyz(2,i0) = xyzstr(2,i0)*(1-one_or_zero_fix_atom(i)) &
                     &                   +      wkxyz(2,i0)*one_or_zero_fix_atom(i)
                wkxyz(3,i0) = xyzstr(3,i0)*(1-one_or_zero_fix_atom(i)) &
                     &                   +      wkxyz(3,i0)*one_or_zero_fix_atom(i)
                wk_f(:,i0)  = wk_f(:,i0) * one_or_zero_fix_atom(i)
             enddo
          enddo
       endif
    endif
#ifdef SEGSHAKE
    call calc_center_of_mass
    dt = 1.0d-15  ![fs]
    call  shake_com(dt)
#endif
#ifdef TIP4
    call update_msite_coordinate
#endif
    ! -- resize of step length --sub
    !       if(mdstep .ge. 2) then
    !       if(p_energy .lt. potential_energy_old)then
    !          step_length = step_length * up_step_length
    !       else
    !          step_length = step_length * down_step_length
    !          if(step_length .lt. eps_step_length) then
!#ifdef MPIPARA
    !             if(myrank.eq.0) then
!#endif
    !             write(f_run,9500)
    !             write(f_run,9200)
    !             write(f_run,9210) step_length
    !             write(f_trj,9610)
!#ifdef MPIPARA
    !             endif
!#endif
    !             itmp = trj_interval
    !             trj_interval = 1
    !             call record_current_trajectory
    !             trj_interval = itmp
    !             stop
    !          endif
    !       endif
    !       endif
    !       potential_energy_old = p_energy

9200 format('ERROR: step length in optimization is too small')
9210 format('step length:  ',e18.10)
9500 format('*************************')
9510 format('Optimization is Converged')
9600 format('# [ Optimized coordinate ]')
9610 format('# [ Not Optimized coordinate ]')
  end subroutine opt_integrate

!-------------------------------------------------------------------------
!
!     calculate summation of root mean square of force
!
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to summarize root meas square of force.
!! \author  Atsushi Yamada
!<
  subroutine opt_calculate_mean_force(mean_force)
    use forces
    use trajectory_org
    use trajectory_mpi
    use position_constrain
    use atom_mass
    use param
#ifdef TIP4
    use tip4p, only : msite
#endif
    implicit none
    integer(4) :: i,i0,k0, n_freeatom, ipar
    real(8) :: mean_force, tmp
    real(8) :: wk_mean_force
!#ifdef MPIPARA
    integer(4) :: ierr
    include 'mpif.h'
!#endif
    if(type_p_constrain==3) then
       n_freeatom = n - npconstrain
    else
       n_freeatom = n
    endif
#ifdef TIP4
    n_freeatom = n_freeatom - msite
#endif
    wk_mean_force = 0d0

    IF(.not.allocated(one_or_zero_fix_atom))THEN
!$omp parallel do default(none) &
!$omp& private(k0,i0,tmp) &
!$omp& shared(lsegtop,lseg_natoms,nselfseg) &
!$omp& shared(wk_f) &
!$omp& reduction(+:wk_mean_force) &
!$omp& shared(m2i,mass,paranum) &
!$omp& private(i,ipar)
       do k0=1,nselfseg
          do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
#ifdef TIP4
       i=m2i(i0)
       ipar=paranum(i) 
       if(mass(ipar) .lt. 0) cycle
#endif
             tmp = wk_f(1,i0)**2 + wk_f(2,i0)**2 + wk_f(3,i0)**2
             wk_mean_force = wk_mean_force + tmp
          enddo ! i0
       enddo ! k0
    ELSE
!$omp parallel do default(none) &
!$omp& private(k0,i0,tmp) &
!$omp& shared(lsegtop,lseg_natoms,nselfseg) &
!$omp& shared(wk_f,mass,paranum) &
!$omp& reduction(+:wk_mean_force) &
!$omp& shared(one_or_zero_fix_atom,m2i) &
!$omp& private(i,ipar)
       do k0=1,nselfseg
          do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
             i = m2i(i0)
#ifdef TIP4
       ipar=paranum(i)
       if(mass(ipar) .lt. 0) cycle
#endif
             tmp = wk_f(1,i0)**2 + wk_f(2,i0)**2 + wk_f(3,i0)**2
             wk_mean_force = wk_mean_force + tmp * dble(one_or_zero_fix_atom(i))
          enddo ! i0
       enddo ! k0
    ENDIF

    !     ^^^ reduce over process ^^^
    mean_force = 0d0
    call mpi_allreduce(wk_mean_force,mean_force,1, &
         &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

    mean_force = mean_force/dble(n_freeatom)
    mean_force = dsqrt( mean_force )

  end subroutine opt_calculate_mean_force

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate direction unit vector for optimization.
!! \author  Atsushi Yamada
!<
  subroutine opt_calc_direction_vector(mean_force)
    use forces
    use trajectory_org
    use trajectory_mpi
    use position_constrain
    use atom_mass
    use param
#ifdef TIP4
    use tip4p, only : msite
#endif
    implicit none
    integer(4) :: i,i0,k0,n_freeatom,ipar
    real(8) :: coef, mean_force

    if(type_p_constrain==3) then
       n_freeatom = n - npconstrain
    else
       n_freeatom = n
    endif
#ifdef TIP4
    n_freeatom = n_freeatom - msite
#endif
    coef = mean_force *dsqrt(dble(n_freeatom))
    coef = 1.0d0/coef

    IF(.not.allocated(one_or_zero_fix_atom))THEN
!$omp parallel do default(shared) &
!$omp& private(i0,k0) &
!$omp& private(i,ipar)
       do k0=1,nselfseg
          do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
#ifdef TIP4
       i=m2i(i0)
       ipar=paranum(i)
       if(mass(ipar) .lt. 0) cycle
#endif
             opt_direction(1,i0) = coef*wk_f(1,i0)
             opt_direction(2,i0) = coef*wk_f(2,i0)
             opt_direction(3,i0) = coef*wk_f(3,i0)
          enddo ! i0
       enddo ! k0
    ELSE
!$omp parallel do default(shared) &
!$omp& private(i0,k0) &
!$omp& shared(one_or_zero_fix_atom) &
!$omp& private(i,ipar)
       do k0=1,nselfseg
          do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
             i = m2i(i0)
#ifdef TIP4
       ipar=paranum(i)
       if(mass(ipar) .lt. 0) cycle
#endif
             opt_direction(1,i0) = coef*wk_f(1,i0) * dble(one_or_zero_fix_atom(i))
             opt_direction(2,i0) = coef*wk_f(2,i0) * dble(one_or_zero_fix_atom(i))
             opt_direction(3,i0) = coef*wk_f(3,i0) * dble(one_or_zero_fix_atom(i))
          enddo ! i0
       enddo ! k0
    ENDIF
  end subroutine opt_calc_direction_vector

!----------------------------------------------------------------------
!
!     calculate thermodynamic values
!
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate outputs (opt)
!! \author Atsushi Yamada
!<
  subroutine calc_potential_opt()
    use thermostat
    use atom_virial
    use md_const
    use md_monitors
    use lj_mod
    use md_condition
    use trajectory_mpi
    use unit_cell
    use subcell
    implicit none
    integer(8) :: i0,k0
    real(8) :: totke=0d0
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

    !!     ^^^ reduce virial ^^^
    !      call calc_scalervirial() ! with correction

    !     ^^^ initialize velocity ^^^
!$omp parallel do default(shared) &
!$omp& private(i0,k0)
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          wkv(1:3,i0)=0d0
       enddo
    enddo

    !     ^^^ temperature ^^^
    k_energy    = totke
    temperature = 2.0d0*totke*degree_of_freedom_inverse*rvkbolz
    t_energy    = p_energy + k_energy
    hamiltonian = t_energy !+ st
    pressure=(2.0d0*totke+virialSca)/(3d0*cellvol)
  end subroutine calc_potential_opt

end module opt_mod
