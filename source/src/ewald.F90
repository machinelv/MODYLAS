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
!! \brief  Module and subroutines which relate to the Ewald method.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to calculate long-range part of the Ewald method.
!! \author Kensuke Iwahashi, Noriyuki Yoshii, Ryo Urano
!<
module ewald_mod
  use mpi_tool
  implicit none
  real(8) :: potential_constant
  real(8) :: pot_charged
  logical:: bQcalc1
!Ewald method
  integer(4) :: max_h2, nhvectors
  integer(4),allocatable :: hvectors_ix(:),hvectors_iy(:), &
 &                          hvectors_iz(:)
  real(8),allocatable :: cos_stocks(:), sin_stocks(:)
contains

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to allocate arrays.
!! \author Kensuke Iwahashi, Noriyuki Yohsii
!<
  subroutine fmod_alloc_ewald_arraies()
    use trajectory_mpi, only : nadirect
    implicit none
!  allocate memory for cos and sin buffer
    call fmod_alloc_md_ewald__cos_stocks(nadirect)
    call fmod_alloc_md_ewald__sin_stocks(nadirect)
  end subroutine fmod_alloc_ewald_arraies
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to set parameter.
!! \author Kensuke Iwahashi, Noriyuki Yohsii
!<
  subroutine fmod_md_ewald__max_h2(ivalue)
  implicit none
  integer(4) :: ivalue
  max_h2 = ivalue
  end subroutine fmod_md_ewald__max_h2
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to set parameter.
!! \author Kensuke Iwahashi, Noriyuki Yohsii
!<
  subroutine fmod_md_ewald__howmany_hvectors(ivalue)
  implicit none
  integer(4) :: ivalue
  nhvectors = ivalue
  end subroutine fmod_md_ewald__howmany_hvectors
!----------------------------------------------------------------------
  subroutine fmod_alloc_hvectors_ix(ivalue)
  implicit none
  integer(4) :: ivalue
  allocate(hvectors_ix(ivalue))
  end subroutine fmod_alloc_hvectors_ix
!----------------------------------------------------------------------
  subroutine fmod_alloc_hvectors_iy(ivalue)
  implicit none
  integer(4) :: ivalue
  allocate(hvectors_iy(ivalue))
  end subroutine fmod_alloc_hvectors_iy
!----------------------------------------------------------------------
  subroutine fmod_alloc_hvectors_iz(ivalue)
  implicit none
  integer(4) :: ivalue
  allocate(hvectors_iz(ivalue))
  end subroutine fmod_alloc_hvectors_iz
!----------------------------------------------------------------------
  subroutine fmod_alloc_md_ewald__cos_stocks(ivalue)
  implicit none
  integer(4) :: ivalue
  allocate(cos_stocks(ivalue))
  end subroutine fmod_alloc_md_ewald__cos_stocks
!----------------------------------------------------------------------
  subroutine fmod_alloc_md_ewald__sin_stocks(ivalue)
  implicit none
  integer(4) :: ivalue
  allocate(sin_stocks(ivalue))
  end subroutine fmod_alloc_md_ewald__sin_stocks
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to initialize the Ewald method.
!! \author Kensuke Iwahashi, Noriyuki Yohsii
!<
  subroutine init_md_ewald()

!  prepare a set of h-vectors
   call md_prepare_hvectors_ewald()     !! note:: half-hvector
!  calculate constant part of potential energy in Ewalds method
   call md_calculate_ewald_e_constant()
!
  end subroutine init_md_ewald
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to prepare h-vectors.
!! \author Kensuke Iwahashi, Noriyuki Yohsii
!<
  subroutine md_prepare_hvectors_ewald()
!
!     make a set of hvectors which is used for calculation of
!     Ewalds long-range force
!     hvector is a vector whose all elements are integer
!     and its norm is less than or equal to dsqrt(md_ewald__max_h2)
!     to avoid time consuming process, we execute calculation on
!     a half of a set of hvectors because dcos(-x) = dcos(x)
!     so, the set of hvectors consists of such hvectors as:
!    
!     (1) hx > 0
!      hx = 1, 2, ... , MAX
!      hy = -MAX, -MAX + 1, ... , MAX
!      hz = -MAX, -MAX + 1, ... , MAX
!    
!     (2) hx == 0 && hy > 0
!      hx = 0
!      hx = 1, 2, ..., MAX
!      hz = -MAX, -MAX + 1, ... , MAX
!    
!     (3) hx == 0 && hy == 0 && hz > 0
!      hx = hy = 0
!      hz = 1, 2, ... , MAX */
!    
      implicit none
      integer(4) :: m, i, j, k, max_h, Nh
      max_h = int(dsqrt(dble(max_h2)))
!
!     count how many hvectors in a set
!     and allocate memory */
      Nh = 0
      do i=1, max_h
        do j=-max_h, max_h
          do k=-max_h, max_h
            if (i * i + j * j + k * k <= max_h2) then
              Nh = Nh + 1
            endif
          enddo
        enddo
      enddo
      do j=1, max_h
        do k=-max_h, max_h
          if(j * j + k * k <= max_h2) then
            Nh = Nh + 1
          endif
        enddo
      enddo
      Nh = Nh + max_h
      call fmod_md_ewald__howmany_hvectors(Nh)
      call fmod_alloc_hvectors_ix(Nh)
      call fmod_alloc_hvectors_iy(Nh)
      call fmod_alloc_hvectors_iz(Nh)

!     register hvectors
      m = 0
      do k=-max_h, max_h
        do j=-max_h, max_h
          do i=1, max_h
            if (i * i + j * j + k * k <= max_h2) then
              m = m + 1
              hvectors_ix(m) = i
              hvectors_iy(m) = j
              hvectors_iz(m) = k
            endif
          enddo
        enddo
      enddo
      do k=-max_h, max_h
        do j=1, max_h
          if (j * j + k * k <= max_h2) then
            m = m + 1
            hvectors_ix(m) = 0
            hvectors_iy(m) = j
            hvectors_iz(m) = k
          endif
        enddo
      enddo
      do k=1, max_h
        m = m + 1
        hvectors_ix(m) = 0
        hvectors_iy(m) = 0
        hvectors_iz(m) = k
      enddo
  end subroutine md_prepare_hvectors_ewald
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate 3rd term of the Ewald method.
!! \author Kensuke Iwahashi, Noriyuki Yohsii, Ryo Urano
!<
  subroutine md_calculate_ewald_e_constant()
    use trajectory_org
    use trajectory_mpi
    use md_const
    use coulomb_mod
    use ewald_variables
    use param
    use unit_cell
    use nonneutral
    implicit none
    include 'mpif.h'
    integer(4) :: i, ipar
    
    real(8) :: qsum_bg
    qsum=0.0d0
    qsum_bg=0.0d0

    do i=1,n
       ipar=paranum(i)
       qsum=qsum+chgv(ipar)
    enddo

    
    !### potconst ###!
    potential_constant=0d0
!$omp parallel default(shared) &
!$omp& private(i,ipar) &
!$omp& reduction(+:potential_constant)
!$omp do
    do i=1,n
       ipar=paranum(i)
       potential_constant=potential_constant+chgv(ipar)**2
    enddo
!$omp end do
!$omp end parallel


! write(*,*) "cellvol",cellvol,"qsum",qsum

    if(qsum .lt. 0.0000000001 .and.   qsum .gt. -0.000000001) then
       pot_charged=0.0d0
    else
       if (myrank.eq. 0) write(*,*) "***********************"
       if (myrank.eq. 0) write(*,*) "system is non-charge neutral"
       if (myrank.eq. 0) write(*,*) "***********************"
       pot_charged=-md_QQ_4PiE/ ewald_alpha / ewald_alpha/cellvol /2*qsum*qsum     
    endif

!! self energy
    potential_constant = -potential_constant * md_QQ_4PiE * ewald_alpha * r_PI_sqrt

  end subroutine md_calculate_ewald_e_constant
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate 3rd term of the Ewald method.
!! \author Kensuke Iwahashi, Noriyuki Yohsii, Ryo Urano
!<
  subroutine md_add_coulomb_ewald_const()
    use md_monitors
    use md_const
    use atom_virial
    implicit none
    include 'mpif.h'
    !     adding constant part of potential
    if (myrank == 0) then
       wk_p_energy = wk_p_energy + potential_constant

       wk_p_energy = wk_p_energy +pot_charged
       
#ifdef DEBUGFCE
#ifdef KCAL
       write(*,*) 'Pot(const)=', potential_constant*kJ_mol/4.184d0,'[kcal/mol]'
       write(*,*) 'Pot(chrgd)=', pot_charged*kJ_mol/4.184d0,'[kcal/mol]'
#else
       write(*,*) 'Pot(const)=', potential_constant*kJ_mol,'[kJ/mol]'
       write(*,*) 'Pot(chrgd)=', pot_charged*kJ_mol,'[kJ/mol]'
#endif
#endif
       
    endif
  end subroutine md_add_coulomb_ewald_const
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate 2nd term of the Ewald method.
!! \author Kensuke Iwahashi, Noriyuki Yohsii
!<
  subroutine md_add_ewald_lattice()
!
!    Ewald methods long-range forces (modified by yandoh)
!    Detail of method is described in Mol. Phys., 50, 1055 (1983).
!
  use math_const, only : PI
  use physics_const, only : md_QQ_4PiE
  use ewald_variables, only : ewald_alpha
  use unit_cell, only : cellx,celly,cellz, &
     &                  alpha,beta,gamma,cellvol
  use openmp_tool
  use trajectory_mpi, only : wkxyz
  use forces, only : w3_f
  use coulomb_mod, only : chgv
  use param, only : paranum
  use subcell, only : m2i, nselfseg, lsegtop, lseg_natoms
  use md_monitors, only : wk_p_energy
  use atom_virial
  use cell_shape, only : cell_convert1
  implicit none
  include 'mpif.h'
  integer(4) :: m
  integer(4) :: i0,k0
  integer(4) :: ipar,ierr
  real(8) :: hL2, exp_h, Scos, Ssin, ph, coef
  real(8) :: wk_Scos, wk_Ssin, wk_Sx(2),Sx(2)
  real(8) :: m_pi2_a2, m_pi_bai, Cf, Cu
  real(8) :: coefkakb, QSS, rcellvol
  real(8) :: kvec(3)
  real(8) :: Scostmp, Ssintmp, Cff_o, Cff_h
  real(8) :: ec, eewa
  real(8) :: box(3,3)
  real(8) :: vpab(3),vpbc(3),vpca(3)
  real(8) :: wk_v11,wk_v22,wk_v33
  real(8) :: wk_v21,wk_v31,wk_v32

  call cell_convert1(cellx,celly,cellz,alpha,beta,gamma,box)

! coefficient of reciprocal vectors
  vpbc(1) = box(2,2)*box(3,3)  ! br2*cr3  (recip11) 
  vpbc(2) =-box(3,3)*box(1,2)  !-cr3*br1  (recip21)
  vpbc(3) = box(1,2)*box(2,3)-box(1,3)*box(2,2) !br1*cr2-cr1*br2 (recip31)
  vpca(1) = 0d0
  vpca(2) = box(3,3)*box(1,1)  ! cr3*ar1  (recip22)
  vpca(3) =-box(1,1)*box(2,3)  ! -ar1*cr2 (recip32)
  vpab(1) = 0d0
  vpab(2) = 0d0
  vpab(3) = box(1,1)*box(2,2)  ! ar1*br2  (recip33)
  rcellvol = 1.0d0/cellvol

  m_pi2_a2 = + PI * PI / (ewald_alpha * ewald_alpha)
  m_pi_bai = 2d0*PI
  Cf = 4.0d0 * md_QQ_4PiE * rcellvol
  Cu =         md_QQ_4PiE * rcellvol / PI

  ec =0d0
  wk_v11=0d0; wk_v22=0d0; wk_v33=0d0
  wk_v21=0d0; wk_v31=0d0; wk_v32=0d0 

! summation over all h-vector
  DO m=1,nhvectors

! making kvector
    kvec(1)=vpbc(1)*dble(hvectors_ix(m))
    kvec(2)=vpbc(2)*dble(hvectors_ix(m))  &
 &         +vpca(2)*dble(hvectors_iy(m))
    kvec(3)=vpbc(3)*dble(hvectors_ix(m))  &
 &         +vpca(3)*dble(hvectors_iy(m))  &
 &         +vpab(3)*dble(hvectors_iz(m))
!
    kvec(1)=kvec(1)*rcellvol
    kvec(2)=kvec(2)*rcellvol
    kvec(3)=kvec(3)*rcellvol
    hL2 = kvec(1)*kvec(1) + kvec(2)*kvec(2) + kvec(3)*kvec(3)
    exp_h = dexp(-m_pi2_a2 * hL2) / hL2   ! Q(k)

! stocking partial value of triangle function
    wk_Scos = 0.0d0
    wk_Ssin = 0.0d0
!$omp parallel default(none) &
!$omp& private(i0,k0,ipar,ph) &
!$omp& shared(m_pi_bai,kvec,wkxyz,chgv,paranum,m2i) &
!$omp& shared(nselfseg,lsegtop,lseg_natoms) &
!$omp& shared(cos_stocks,sin_stocks) &
!$omp& reduction(+:wk_Scos,wk_Ssin)
!$omp do
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
            ph = m_pi_bai*( &
    &            kvec(1)*wkxyz(1,i0) &
    &          + kvec(2)*wkxyz(2,i0) &
    &          + kvec(3)*wkxyz(3,i0) )
            ipar=paranum(m2i(i0))
            cos_stocks(i0) = chgv(ipar) * dcos(ph)
            sin_stocks(i0) = chgv(ipar) * dsin(ph)
            wk_Scos = wk_Scos + cos_stocks(i0)
            wk_Ssin = wk_Ssin + sin_stocks(i0)
       enddo ! i0
    enddo ! k0
!$omp end do
!$omp end parallel

    wk_Sx(1)=wk_Scos
    wk_Sx(2)=wk_Ssin
! sum over all MPI process
    call mpi_allreduce(wk_Sx,Sx,2,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    Scos=Sx(1)
    Ssin=Sx(2)

    Scostmp=exp_h*Scos
    Ssintmp=exp_h*Ssin

! ^^^ potential, sum over k-vector ^^^
    QSS=Scostmp*Scos+Ssintmp*Ssin !QSS=exp_h*(Scos*Scos+Ssin*Ssin)  ! QSS depends k-vector ###
    QSS=Cu*QSS
    ec = ec + QSS

!   ^^^ virial, sum over k-vector ^^^
    coefkakb= (2d0/hL2)*( 1d0+m_pi2_a2*hL2 )*QSS
!   coefkakb=+2d0*( 1d0+m_pi2_a2*hL2 )*(QSS/hL2)
    wk_v11=wk_v11 - coefkakb*kvec(1)*kvec(1) + QSS 
    wk_v22=wk_v22 - coefkakb*kvec(2)*kvec(2) + QSS
    wk_v33=wk_v33 - coefkakb*kvec(3)*kvec(3) + QSS
    wk_v21=wk_v21 - coefkakb*kvec(1)*kvec(2)
    wk_v31=wk_v31 - coefkakb*kvec(1)*kvec(3)
    wk_v32=wk_v32 - coefkakb*kvec(2)*kvec(3)
!

!   ^^^ force, sum over k-vector ^^^
!$omp parallel default(none) &
!$omp& private(i0,k0,coef) &
!$omp& shared(nselfseg,lsegtop,lseg_natoms) &
!$omp& shared(Scostmp,sin_stocks,Ssintmp,cos_stocks) &
!$omp& shared(kvec,w3_f,Cf) 
!$omp do
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
         coef = Scostmp * sin_stocks(i0) - Ssintmp * cos_stocks(i0)
         w3_f(1,i0,0)   = w3_f(1,i0,0)   + kvec(1)*coef * Cf
         w3_f(2,i0,0)   = w3_f(2,i0,0)   + kvec(2)*coef * Cf
         w3_f(3,i0,0)   = w3_f(3,i0,0)   + kvec(3)*coef * Cf
       enddo ! i0
    enddo ! k0
!$omp end do
!$omp end parallel

  ENDDO ! h-vectors


  if(myrank==0)then
    wk_p_energy = wk_p_energy + ec
    wk_vir2(1,0)=wk_vir2(1,0)+wk_v11
    wk_vir2(2,0)=wk_vir2(2,0)+wk_v22
    wk_vir2(3,0)=wk_vir2(3,0)+wk_v33
    wk_vir2(4,0)=wk_vir2(4,0)+wk_v21
    wk_vir2(5,0)=wk_vir2(5,0)+wk_v31
    wk_vir2(6,0)=wk_vir2(6,0)+wk_v32
  endif

#ifdef DEBUGFCE
!    values in ec was already the sum of all nodes. 
!    please refer to line 381 (as of 2022.3.4) "call mpi_allreduce(wk_Sx,Sx,2,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)"
!    energy ec is computed from Sx which have been "mpi_allreduced" in line 381. 
!    call mpi_allreduce(ec,eewa,1, &
!         &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
#ifdef KCAL
    if(myrank==0) write(*,*) 'Pot(Ewald)=',ec*kJ_mol/4.184d0,'[kcal/mol]'
#else
    if(myrank==0) write(*,*) 'Pot(Ewald)=',ec*kJ_mol,'[kJ/mol]'
#endif
#endif

  end subroutine md_add_ewald_lattice

end module ewald_mod
