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
!! \brief  Modules to store global variables.
!<
!----------------------------------------------------------------------
!>
!! \brief  Modules to store force arrays.
!! \author Tatsuya Sakashita
!<
module forces
  real(8),allocatable :: f(:,:)
  real(8),allocatable :: wk_f(:,:)
  real(8),allocatable :: w3_f(:,:,:)
end module forces

!>
!! \brief  Modules to store variables for ensemble.
!! \author Tatsuya Sakashita
!<
module extended_system
  implicit none
  real(8) :: systemp=3.0d2,syspres=101325.00000
  real(8) :: tautherm=1.0d-12,taubaro=1.0d-12,taubaroW=1.0d-12
  real(8) :: fts0=1.d12,fts1=1.d12,fbs=1.d10,fbs0=1.d12,fbs1=1.d12
  real(8) :: sysvol0=0d0
  real(8) :: cellx0=0d0, celly0=0d0, cellz0=0d0
  real(8) :: alpha0=0d0, beta0=0d0, gamma0=0d0

contains

  subroutine fmod_systemp(value)
    implicit none
    real(8), intent(in) :: value
    systemp = value
  end subroutine fmod_systemp
!----------------------------------------------------------------------
  subroutine fmod_sysvolume(value)
    implicit none
    real(8), intent(in) :: value
    sysvol0 = value
  end subroutine fmod_sysvolume
!----------------------------------------------------------------------
  subroutine fmod_syscellx(value)
    use unit_cell
    implicit none
    real(8), intent(in) :: value
    cellx0 = value
  endsubroutine fmod_syscellx
!----------------------------------------------------------------------
  subroutine fmod_syscelly(value)
    use unit_cell
    implicit none
    real(8), intent(in) :: value
    celly0 = value
  endsubroutine fmod_syscelly
!----------------------------------------------------------------------
  subroutine fmod_syscellz(value)
    use unit_cell
    implicit none
    real(8), intent(in) :: value
    cellz0 = value
  endsubroutine fmod_syscellz
!----------------------------------------------------------------------
  subroutine fmod_sysalpha(value)
    use unit_cell
    implicit none
    real(8), intent(in) :: value
    alpha0 = value
  endsubroutine fmod_sysalpha
!----------------------------------------------------------------------
  subroutine fmod_sysbeta(value)
    use unit_cell
    implicit none
    real(8), intent(in) :: value
    beta0 = value
  endsubroutine fmod_sysbeta
!----------------------------------------------------------------------
  subroutine fmod_sysgamma(value)
    use unit_cell
    implicit none
    real(8), intent(in) :: value
    gamma0 = value
  endsubroutine fmod_sysgamma
!<--NI
!----------------------------------------------------------------------
  subroutine fmod_tauqtherm(value)
    implicit none
    real(8), intent(in) :: value
    tautherm = value
    fts0 = 1.0d0 / tautherm
    fts1 = fts0
  end subroutine fmod_tauqtherm
!----------------------------------------------------------------------
  subroutine fmod_tauqbaro(value)
    implicit none
    real(8), intent(in) :: value
    taubaro = value
    fbs0 = 1.0d0 / taubaro
    fbs1 = fbs0
  end subroutine fmod_tauqbaro
!----------------------------------------------------------------------
  subroutine fmod_tauwpres(value)
    implicit none
    real(8), intent(in) :: value
    taubaroW = value
    fbs = 1.0d0 / taubaroW
  end subroutine fmod_tauwpres

end module extended_system
!----------------------------------------------------------------------
!>
!! \brief  Modules to store variables for multiple time steps (RESPA).
!! \author Tatsuya Sakashita, Yoshimichi Andoh
!<
module md_multiplestep
  implicit none
  integer(4) :: MTs, MTm, MTl
  integer(4) :: maxMTm=1
  integer(4) :: maxMTl=1
  real(8),allocatable :: fshort(:,:),fmiddle(:,:),flong(:,:)
  real(8) :: virshort(6),virmiddle(6),virlong(6)
  real(8) :: eneshort,enemiddle,enelong
  real(8) :: dthL
  real(8) :: dtL
  integer(4) :: XORESPA=1
  integer(4) :: XIRESPA=0
  real(8) :: scaleM=1d0 !! YA: Never change this initial value
  real(8) :: scaleL=1d0 !! YA: Never change this initial value
  integer(4) :: iMT
  real(8) :: Pinner=1.01325d+5,Pouter=0d0,Poutermost=0d0

contains

  subroutine fmod_set_maxMTm(ivalue)
    use mpi_tool
    implicit none
    integer(4), intent(in) :: ivalue
    maxMTm = ivalue
  end subroutine fmod_set_maxMTm
!----------------------------------------------------------------------
  subroutine fmod_set_maxMTl(ivalue)
    use mpi_tool
    implicit none
    integer(4), intent(in) :: ivalue
    maxMTl = ivalue
  end subroutine fmod_set_maxMTl
!----------------------------------------------------------------------
  subroutine fmod_set_Pinner(value)
    use mpi_tool
    implicit none
    real(8), intent(in) :: value
    Pinner = value
  end subroutine fmod_set_Pinner
!----------------------------------------------------------------------
  subroutine fmod_set_Pouter(value)
    use extended_system
    use mpi_tool
    implicit none
    real(8), intent(in) :: value

    Pouter = value
  end subroutine fmod_set_Pouter
!----------------------------------------------------------------------
  subroutine fmod_set_Poutermost(value)
    use extended_system, only : syspres
    use mpi_tool
    implicit none
    real(8), intent(in) :: value
    real(8) :: rel
    Poutermost = value
    if(myrank==0)then
       write(*,*) 'Pref_outermost=', Poutermost
    endif
    rel=abs( 1d0-(Pinner+Pouter+Poutermost)/syspres )
    if(rel.ge.1d-3)then  !! three significant figures
       if(myrank==0)then
          write(0,*)'ERROR :: Pref_inner + _outer + _outermost /= syspres'
          write(0,*)Pinner,'+',Pouter,'+',Poutermost,'=', &
               &            Pinner+Pouter+Poutermost
          write(0,*) ',which /=', syspres
       endif
       call modylas_abort()
    endif
  end subroutine fmod_set_Poutermost

  subroutine fmod_syspres(value)
    use extended_system
    implicit none
    real(8), intent(in) :: value
    syspres = value
    Pinner= value
    Pouter=0d0
    Poutermost=0d0
  end subroutine fmod_syspres

end module md_multiplestep

!>
!! \brief  Modules to store variables for LJ cut-off.
!! \author Tatsuya Sakashita
!<
module cutoff_radius
  implicit none
  real(8) :: cutrad=0.0d0, cutrad2=0.0d0

contains

  subroutine fmod_cutrad(value)
    implicit none
    real(8), intent(in) :: value
    cutrad = value
  end subroutine fmod_cutrad

  subroutine fmod_cutrad2(value)
    implicit none
    real(8), intent(in) :: value
    cutrad2 = value
  end subroutine fmod_cutrad2

end module cutoff_radius

!----------------------------------------------------------------------
!>
!! \brief  Modules to store variables for P2P.
!! \author Tatsuya Sakashita, Jiachao Zhang
!<
module md_forces
#if defined(PRECISION_P2P_MIX) || defined(PRECISION_P2P_SP)
  real(4),allocatable :: chgv_table(:,:)
  real(4),allocatable :: epsilon_sqrt_table(:,:)
  real(4),allocatable :: R_half_table(:,:)
#else
  real(8),allocatable :: chgv_table(:,:)
  real(8),allocatable :: epsilon_sqrt_table(:,:)
  real(8),allocatable :: R_half_table(:,:)
#endif
end module md_forces
!----------------------------------------------------------------------
!>
!! \brief  Modules to store variables for maxwell distribution.
!! \author Tatsuya Sakashita
!<
module maxwell_distribution
  implicit none
  real(8) :: maxwell_temperature=0.0d0
  logical :: reset_maxwell=.false.  ! default

contains

  subroutine fmod_maxwell_temperature(value)
    implicit none
    real(8), intent(in) :: value
    maxwell_temperature = value
  end subroutine fmod_maxwell_temperature

  subroutine fmod_reset_maxwell(ivalue)
    use extended_system
    use mpi_tool
    implicit none
    integer(4), intent(in) :: ivalue
    if (ivalue == 1) then
       reset_maxwell = .true.
       maxwell_temperature = systemp
       if(myrank==0)then
          write(*,*) 'Velocities is reset to Maxwell distribution'
          write(*,*) 'with temperature=',maxwell_temperature
       endif
    endif
  end subroutine fmod_reset_maxwell

end module maxwell_distribution

!>
!! \brief  Modules to store variables for .mdmntr file.
!! \author Tatsuya Sakashita
!<
module md_monitors
  real(8) :: hamiltonian
  real(8) :: p_energy, wk_p_energy
  real(8) :: k_energy
  real(8) :: t_energy
  real(8) :: temperature
  real(8) :: pressure
  real(8) :: ptensor(6)
  real(8) :: P_part(3)=0d0      !! 3 means (short,middle,long)
  real(8) :: sum_P=0d0          !! 3 means (short,middle,long)
  real(8) :: sum_P_part(3)=0d0  !! 3 means (short,middle,long)
  !NI-->
  real(8) :: sum_box(3,3)=0.d0, dummy(1:6)=0.d0
  !<--NI
end module md_monitors

!>
!! \brief  Modules to store variables for non-charge neutral system.
!! \author Ryo Urano
!<
module nonneutral
    implicit none
    real(8):: qsum,q_squared_sum
    logical:: bNonNeutrality=.false. !this does not work due to init_fmm_ewald
    logical:: bNonNeutrality_tot=.false.
    integer:: FMM_interval=1
    real(8),parameter :: zero_threshold=1.0e-7
end module nonneutral
