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
!! \brief  Module to store LJ parameter information.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to store LJ parameter information.
!! \author Kensuke Iwahashi, Kazu Fujimoto, Yoshimichi Andoh
!<
module lj_mod
  implicit none
  real(8),allocatable :: epsilon_sqrt(:), R_half(:)
#ifdef CHARMMFSW
  real(8) :: Ron=8d-10,Roff=12d-10 !! default of CHARMM36 is 8,12A
  real(8) :: Ron2,Roff2,Ron3,Roff3,Ron6,Roff6
  real(8) :: b1,b2,b3,c6a,c6b,c6c,c6d
  real(8) :: iron2,iroff3,iroff6,ib3 !,R3,R6,R12
  real(8) :: c12a,c12b,dV2_12,dV2_6
  real(8) :: k_f12,k_f6,Sr
#endif

  real(8) :: corrector_potential=0.0d0, corrector_virial=0.0d0

contains

!>
!! \brief  Subroutine to set sqrt of epsilon.
!! \author Kazushi FUJIMOTO
!<
  subroutine fmod_set_mol_ljs__epsilon_sqrt(value)
    implicit none
    real(8), intent(in) :: value
    integer(4),save :: i=0

    i = i + 1
    epsilon_sqrt(i) = value

  end subroutine fmod_set_mol_ljs__epsilon_sqrt
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to set a half of R.
!! \author Kazushi FUJIMOTO
!<
  subroutine fmod_set_mol_ljs__r_half(value)
    implicit none
    real(8), intent(in) :: value
    integer(4),save :: i=0

    i = i + 1
    R_half(i) = value
  end subroutine fmod_set_mol_ljs__r_half
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to read epsilon_sqrt.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_epsilon_sqrt
    use device_numbers, only : f_mdff
    use param
    use mpi_tool
    implicit none
    include 'mpif.h'
    integer(4) :: ierr

    if(myrank==0) read(f_mdff) epsilon_sqrt
    call MPI_Bcast(epsilon_sqrt, npara, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  end subroutine read_epsilon_sqrt
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write epsilon_sqrt in mdff.bin File.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_epsilon_sqrt
    use device_numbers, only : f_mdff
    implicit none
    write(f_mdff) epsilon_sqrt
  end subroutine write_epsilon_sqrt
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write epsilon_sqrt.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_epsilon_sqrt
    implicit none
    write(*,*) '[write_epsilon_sqrt]'
    write(*,*) epsilon_sqrt
  end subroutine write_memory_epsilon_sqrt
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to read R_half.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_R_half
    use device_numbers, only : f_mdff
    use param
    use mpi_tool
    implicit none
    include 'mpif.h'
    integer(4) :: ierr
    if(myrank==0) read(f_mdff) R_half
    call MPI_Bcast(R_half, npara, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  end subroutine read_R_half
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write R_half in mdff.bin File.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_R_half
    use device_numbers, only : f_mdff
    implicit none
    write(f_mdff) R_half
  end subroutine write_R_half
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write R_half.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_R_half
    implicit none
    write(*,*) '[write_R_half]'
    write(*,*) R_half
  end subroutine write_memory_R_half
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate LJ correction terms
!! \author Yoshimichi Andoh
!<
  subroutine md_calc_corrector_constant()
!-----------------------------------------------------------------------
    use trajectory_org
    use md_forces
    use md_const
    use param
    use cutoff_radius
    use unit_cell
    use mpi_tool
    use md_condition
    implicit none
    integer(4) :: i, j, ipar
    real(8) :: i_rc3, eps, R3, R_rc3
    real(8) :: Uc, Wc
    integer(4) :: maxtype
    real(8) :: nn
    integer(4),allocatable :: type_n(:)
    real(8),allocatable :: type_eps(:),type_r(:)

    allocate( type_n(npara),type_eps(npara),type_r(npara) )
    maxtype=npara

    type_n = 0
    do i=1, n
       ipar=paranum(i)
       type_n(ipar)   = type_n(ipar)+1
       type_eps(ipar) = epsilon_sqrt(ipar)
       type_r(ipar)   = R_half(ipar)
    enddo

    i_rc3 = 1.0d0 / (cutrad2*cutrad)
    Uc = 0.0d0
    Wc = 0.0d0
    do i=1, maxtype
       !
       !       Part I: The same kind of particle
       !
       eps = type_eps(i) * type_eps(i)
       R3 = type_r(i) + type_r(i)
       R3 = R3 * R3 * R3
       R_rc3 = R3 * i_rc3
       nn = dble(type_n(i))*dble((type_n(i)-1))*0.5d0
       Uc = Uc + eps * R3 * R_rc3 * (R_rc3 * R_rc3 - 6.0d0) * nn
       Wc = Wc + eps * R3 * R_rc3 * (R_rc3 * R_rc3 - 3.0d0) * nn
       !
       !       Part II: Different kinds of particle
       !
       do j=i+1, maxtype
          eps = type_eps(i) * type_eps(j)
          R3 = type_r(i) + type_r(j)
          R3 = R3 * R3 * R3
          R_rc3 = R3 * i_rc3
          nn = dble(type_n(i))*dble(type_n(j))
          Uc = Uc + eps * R3 * R_rc3 * (R_rc3 * R_rc3 - 6.0d0) * nn
          Wc = Wc + eps * R3 * R_rc3 * (R_rc3 * R_rc3 - 3.0d0) * nn
       enddo
    enddo
    corrector_potential = Uc* 4.0d0* PI / 9.0d0
    corrector_virial    = Wc*16.0d0* PI / 9.0d0
    deallocate( type_n,type_eps,type_r ) !! closed here
#ifdef CHARMMFSW
    if(myrank==0) write(*,*) '-DCHARMMFSW is set, thus LJ_LRC ommited'
    LJ_LRC=.false.
#endif
  end subroutine md_calc_corrector_constant
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to initialize parameters for CHARMM force switching function
!! \author Yoshimichi Andoh
!<
!-----------------------------------------------------------------------
  subroutine init_charmmfsw
#ifdef CHARMMFSW
    implicit none

    Ron2=Ron**2; Roff2=Roff**2
    Ron3=Ron**3; Roff3=Roff**3
    Ron6=Ron**6; Roff6=Roff**6
    b1=Ron3*Roff3 ; b2=Roff3-Ron3
    b3=Roff2-Ron2 ; ib3=1d0/b3**3
    c6a=    1d0/b1
    c6b=  Roff3/b2
    c6c=   -2d0/b2
    c6d=    1d0/b2/Roff3
    iRon2 =1d0/Ron2
    iRoff3=1d0/Roff3
    iRoff6=1d0/Roff6
    c12a=    1d0/b1**2
    c12b=Roff6/(Roff6-Ron6)
#endif
  end subroutine 

end module lj_mod
