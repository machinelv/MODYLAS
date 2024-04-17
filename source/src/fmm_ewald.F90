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
!! \brief  Module and subroutines to calculate lattice sum of multipole (FMM ewald).
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to calculate FMM ewald.
!! \author Noriyuki Yoshii, Kensuke Iwahashi, Ryo Urano
!<
module fmm_ewald
  use mpi_tool
  use fmm_l_m_index
  use nonneutral , only :  FMM_interval, bNonNeutrality, &
                         & bNonNeutrality_tot
  implicit none
  complex(8), private, allocatable :: shew1(:,:), shew2(:,:)
  complex(8),allocatable,dimension(:,:) :: shew
  complex(8),allocatable,dimension(:,:) :: wk_we
  ! logical:: bNonNeutrality=.true.
  logical::bQcalc1= .true. 
#ifdef DEBUGFCE
  real(8),parameter :: kJ_mol=6.02214129d+23*1.0d-3
#endif

contains
!----------------------------------------------------------------------
!>
!! \brief Subroutine for initializing the Ewald sum of FMM calculation
!! \author Noriyuki Yoshii, Kensuke Iwahashi
!<
  subroutine init_fmm_ewald(nmax, lm_length)
!----------------------------------------------------------------------
    use unit_cell, only : cellx, celly, cellz
    use math_const, only : PI
    use spherical_harmonics, only : cart2angle, lngamma, incmpgamm
    use regular_singular_solid_harmonics
    use nonneutral
    implicit none
    integer(4), intent(in) :: nmax, lm_length
    complex(8) :: wk_we(0:2*nmax,-2*nmax:2*nmax)
    integer(4) :: nhmax,nrmax
    real(8) :: wkappa,nhcut,nrcut,wkappa_square
    real(8) :: nhcut_square,nrcut_square,vol,nh_square,nr_square
    integer(4) :: i,j,k,m,n
    real(8) :: xta,yta,zta,rad,the,csthe,phi
    real(8) :: fvh
    complex(8) :: fi
    real(8) :: s_plus_re, s_plus_im, s_minus_re, s_minus_im
    integer(4) :: m1, m2
    real(8), parameter :: PI_square = PI**2

    !******* set Ewald parameter
    nhmax=20
    nrmax=20
    wkappa=2.d0
    nhcut=nhmax
    nrcut=nrmax
    !     vol=1.d0
    vol=celly*cellz/cellx**2
    !******* set Ewald parameter
    wkappa_square = wkappa**2
    nhcut_square = nhcut**2
    nrcut_square = nrcut**2

    ! allocate(wk_we(((2*nmax+1))*(2*nmax+1)))
    allocate( shew1(lm_length,lm_length) )
    allocate( shew2(lm_length,lm_length) )

    wk_we = dcmplx(0.d0,0.d0)

 call Smat_1(wkappa,vol,nmax,nhmax,wk_we &
   &  ) 

!!! call calc_qtot()  !calc system charge  !! moved to init_fmm

 if(bNonNeutrality) wk_we(0,0) = wk_we(0,0) -2.0d0*wkappa/sqrt(pi) !! charged system term

 call Smat_2(wkappa,vol,nmax,nrmax,wk_we,shew1,shew2 &
    &       ) 

  end subroutine init_fmm_ewald


!----------------------------------------------------------------------
subroutine Smat_1(wkappa,vol,nmax,nhmax,wk_we &
             &   )
!----------------------------------------------------------------------
    use unit_cell, only : cellx, celly, cellz
    use math_const, only : PI
    use spherical_harmonics, only : cart2angle, lngamma, incmpgamm
    use regular_singular_solid_harmonics
    implicit none 
    integer,intent(in)::nmax
    integer,intent(in)::nhmax
    real(8),intent(in) :: wkappa,vol
    complex(8),intent(inout) :: wk_we(0:2*nmax,-2*nmax:2*nmax)

    real(8) :: nhcut,nrcut,wkappa_square
    real(8) :: nhcut_square,nrcut_square,nh_square,nr_square
    integer(4) :: i,j,k,m,n
    real(8) :: xta,yta,zta,rad,the,csthe,phi
    real(8) :: fvh
    complex(8) :: fi
    real(8) :: s_plus_re, s_plus_im, s_minus_re, s_minus_im
    integer(4) :: m1, m2
    real(8), parameter :: PI_square = PI**2

    !******* set Ewald parameter
    nhcut=nhmax
    !******* set Ewald parameter
    wkappa_square = wkappa**2
    nhcut_square = nhcut**2

    do i=-nhmax,nhmax
       do j=-nhmax,nhmax
          do k=-nhmax,nhmax
             nh_square=i*i+j*j+k*k
             xta=dble(i)
             yta=dble(j)/(celly/cellx)
             zta=dble(k)/(cellz/cellx)
             if(xta**2+yta**2+zta**2.gt.nhcut_square) cycle
             if(xta**2+yta**2+zta**2.eq.0.d0) cycle
             call cart2angle(xta,yta,zta,rad,the,csthe,phi)
             do n=0,2*nmax,FMM_interval
                fi=dcmplx(0.d0,1.d0)**n*PI**(n-0.5d0)
                fvh=dexp(-PI_square*rad**2/wkappa_square-lngamma(n+0.5d0))*rad**(2*n-1)/vol
                do m=-n,n,FMM_interval
                   wk_we(n,m) = wk_we(n,m) + fi*fvh*singular_harmonics(n, m, rad, csthe, phi)
                enddo
             enddo
          enddo
       enddo
    enddo

  end subroutine Smat_1

!----------------------------------------------------------------------
subroutine  Smat_2(wkappa,vol,nmax,nrmax,wk_we,Smat_shew1,Smat_shew2 &
               &  ) 
!----------------------------------------------------------------------
    use unit_cell, only : cellx, celly, cellz
    use math_const, only : PI
    use spherical_harmonics, only : cart2angle, lngamma, incmpgamm
    use regular_singular_solid_harmonics
    implicit none 
    integer,intent(in)::nmax
    integer,intent(in)::nrmax
    real(8),intent(in) :: wkappa,vol
    complex(8),intent(inout) :: wk_we(0:2*nmax,-2*nmax:2*nmax)
    complex(8), intent(inout) :: Smat_shew1(:,:), Smat_shew2(:,:)

    real(8) :: nhcut,nrcut,wkappa_square
    real(8) :: nhcut_square,nrcut_square,nh_square,nr_square
    integer(4) :: i,j,k,m,n
    real(8) :: xta,yta,zta,rad,the,csthe,phi
    real(8) :: fvh
    complex(8) :: fi
    real(8) :: s_plus_re, s_plus_im, s_minus_re, s_minus_im
    integer(4) :: m1, m2
    real(8), parameter :: PI_square = PI**2

    nrcut=nrmax
    wkappa_square = wkappa**2
    nrcut_square = nrcut**2


    do i=-nrmax,nrmax
       do j=-nrmax,nrmax
          do k=-nrmax,nrmax
             nr_square=i*i+j*j+k*k
             if(nr_square.gt.nrcut_square) cycle
             if(iabs(i).le.2.and.iabs(j).le.2.and.iabs(k).le.2)then
                if(nr_square.eq.0.d0) cycle
                xta=dble(i)
                yta=dble(j)*(celly/cellx)
                zta=dble(k)*(cellz/cellx)
                call cart2angle(xta,yta,zta,rad,the,csthe,phi)
                if(rad.ne.0d0)then
                   do n=0,2*nmax,FMM_interval
                      fvh=incmpgamm(n+0.5d0,wkappa_square*rad**2)
                      do m=-n,n
                         wk_we(n,m) = wk_we(n,m) - fvh * singular_harmonics(n, m, rad, csthe, phi)
                      enddo
                   enddo
                endif
             else
                xta=dble(i)
                yta=dble(j)*(celly/cellx)
                zta=dble(k)*(cellz/cellx)
                call cart2angle(xta,yta,zta,rad,the,csthe,phi)
                if(rad.ne.0d0)then
                   do n=0,2*nmax,FMM_interval
                      fvh=1d0-incmpgamm(n+0.5d0,wkappa_square*rad**2)
                      do m=-n,n
                         wk_we(n,m) = wk_we(n,m) + fvh * singular_harmonics(n, m, rad, csthe, phi)
                      enddo
                   enddo
                endif
             endif
          enddo
       enddo
    enddo

    ! multipole to local
    do j=0,nmax
       do k=0,j
          m1 = translate_l_m_to_1dim(j, k)
          do n=0,nmax
             m2 = translate_l_m_to_1dim(n, 0)
             Smat_shew1(m2,m1) = wk_we(j+n, -0-k)
             Smat_shew2(m2,m1) = dcmplx( 0d0, 0d0 )
             do m=1,n
                s_plus_re = (-1)**(m+k) * dreal(wk_we(j+n,m+k))
                s_plus_im = (-1)**(m+k) * (- dimag(wk_we(j+n,m+k)))
                s_minus_re = (-1)**m * dreal(wk_we(j+n,m-k))
                s_minus_im = (-1)**m * dimag(wk_we(j+n,m-k))
                m2 = translate_l_m_to_1dim(n, m)
                if (k /= 0) then
                   Smat_shew1(m2,m1) = dcmplx( s_plus_re + s_minus_re, s_plus_im + s_minus_im)
                   Smat_shew2(m2,m1) = dcmplx(-s_plus_im + s_minus_im, s_plus_re - s_minus_re)
                else
                   Smat_shew1(m2,m1) = dcmplx( s_plus_re + s_minus_re, 0d0)
                   Smat_shew2(m2,m1) = dcmplx(-s_plus_im + s_minus_im, 0d0)
                endif
             enddo
          enddo
       enddo
    enddo

end subroutine Smat_2


!----------------------------------------------------------------------
!>
!! \brief Subroutine which calculate the Ewald sum of FMM calculation
!! \author Noriyuki Yoshii, Kensuke Iwahashi
!<
  subroutine run_fmm_ewald(lm_length,wm,wl &
                      &   )
!----------------------------------------------------------------------
!   use TI_module
    implicit none
    integer(4), intent(in) :: lm_length
    complex(8), intent(in)  :: wm(lm_length)
    complex(8), intent(out) :: wl(lm_length)
    integer(4) :: m1,m2
    wl=dcmplx(0.d0,0.d0)
    do m1=1,lm_length
       do m2=1,lm_length
            wl(m1) = wl(m1) + dreal(wm(m2))*shew1(m2,m1) + dimag(wm(m2))*shew2(m2,m1)             
       enddo
    enddo

  end subroutine run_fmm_ewald

!----------------------------------------------------------------------
  subroutine background_charge_non_neutrality
!----------------------------------------------------------------------
    use omp_lib
    use param
    use coulomb_mod
    use subcell
    use unit_cell
    use trajectory_mpi
    use md_const
    use md_monitors
    use nonneutral
    implicit none
    integer(4) :: i,i0,k0, nomp
    real(8) :: coeff,coeff0

    real(8) :: E_self,Ebg_self2,Ebg_self
    real(8) :: wkappa,wkappa_akma

      integer(4) :: ierr
      include 'mpif.h' 

      wkappa=2.d0

      nomp = 1
!$    nomp = omp_get_max_threads()

!!    if(bQcalc1) call calc_qtot()  !calc system charge  -> moved to init_fmm

      if(myrank==0)then

!! at least cellx is the dimension of angstrom cellxh is the half of cellx
       wkappa_akma=wkappa/cellx

! charged system term
       coeff = 1/(8d0*cellvol*md_VACUUM_DIELECTRIC_CONSTANT*(wkappa_akma*wkappa_akma))

! self-energy
       coeff0=wkappa_akma/(4*(PI*sqrt(PI))*md_VACUUM_DIELECTRIC_CONSTANT) 

!!! self energy 
    Ebg_self= - (q_squared_sum*q_squared_sum)*coeff0

!! charged system term 
    Ebg_self2=-(qsum*qsum)*coeff

     wk_p_energy = wk_p_energy +Ebg_self2
     wk_p_energy = wk_p_energy +Ebg_self
    
#ifdef DEBUGFCE
#ifdef KCAL
     write(*,*) 'Pot(charged_system) =', Ebg_self2*kJ_mol/4.184d0,'[kcal/mol]'
     write(*,*) 'Pot(self) =', Ebg_self*kJ_mol/4.184d0,'[kcal/mol]' !
#else
     write(*,*) 'Pot(charged_system) =', Ebg_self2*kJ_mol,'[kJ/mol]'
     write(*,*) 'Pot(self) =', Ebg_self*kJ_mol,'[kJ/mol]' !
#endif
#endif
    endif
  end subroutine background_charge_non_neutrality

!>
!! \brief Subroutine to calculate total of parial charges.
!! \author Ryo Urano
!<
  subroutine calc_qtot()
      use omp_lib
      use param
      use coulomb_mod
      use subcell
      use md_const
      use nonneutral,only : qsum,q_squared_sum
      use fmm_parameters, only : nmax
      use segments
      implicit none

      integer(4) :: i, k, nomp
      real(8) :: temp,temp_self !,temp_self_s_sq
!     real(8) :: temp_s, temp_w,temp_bg1,temp_bg2,temp_bg
      real(8),parameter :: eps_numerical=0.0000001
      real(8),parameter :: zero_threshold=1.0e-7

      integer(4) :: ierr
      include 'mpif.h' 

      nomp = 1
      temp   = 0.0d0
      temp_self=0.0d0
      
!$omp parallel default(shared) &
!$omp& private(k,i) &
!$omp&  reduction(+:temp,temp_self)
!$omp do
    do k = 1, nsegments   ! nsegment is total number of segment in system
      do i = segtop(k)+1,segtop(k)+seg_natoms(k)
          temp=temp+ chgv(paranum(i))
          temp_self=temp_self+ chgv(paranum(i))**2
      end do ! i
    end do ! k
!$omp end do
!$omp end parallel
    qsum=temp
    q_squared_sum=temp_self

!!$omp do
!     do k0=1,nselfseg
!       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
!         i = m2i(i0)

!         temp=temp+ chgv(paranum(i))
!         temp_self=temp_self+ chgv(paranum(i))**2

!       end do ! i0
!     end do ! k0

!      call mpi_allreduce(temp,qsum,1, &
!    &    mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
 !! self-energy 
!     call mpi_allreduce(temp_self,q_squared_sum,1, &
!    &    mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

      if( IsRealZero(qsum, zero_threshold)) then
         if (myrank.eq. 0) write(*,*) "system is charge neutral"
         bNonNeutrality=.false.
         FMM_interval=2  !! summation interval on FMM
      else
         if (myrank.eq. 0) write(*,*) "***********************"
         if (myrank.eq. 0) write(*,*) "system is non-charge neutral"
         if (myrank.eq. 0) write(*,*) "Qtot =",qsum
         if (myrank.eq. 0) write(*,*) "***********************"
         bNonNeutrality=.true.
         bNonNeutrality_tot=.true.
         FMM_interval=1  !! summation interval on FMM
      endif

       qsum=qsum*md_ELEMENTARY_CHARGE
       q_squared_sum=q_squared_sum*md_ELEMENTARY_CHARGE*md_ELEMENTARY_CHARGE

      end subroutine calc_qtot

!----------------------------------------------------------------------
!>
!! \brief Subroutine for initializing the Ewald sum of FMM calculation on charged system
!! \author Ryo Urano
  !<
  subroutine init_fmm_ewald_charged(nmax, lm_length)
!----------------------------------------------------------------------
    use unit_cell, only : cellx, celly, cellz
    use math_const, only : PI
    use spherical_harmonics, only : cart2angle, lngamma, incmpgamm
    use regular_singular_solid_harmonics
    use nonneutral
    implicit none
    integer(4), intent(in) :: nmax, lm_length
    complex(8) :: wk_we(0:2*nmax,-2*nmax:2*nmax)
    integer(4) :: nhmax,nrmax
    real(8) :: wkappa,nhcut,nrcut,wkappa_square
    real(8) :: nhcut_square,nrcut_square,vol,nh_square,nr_square
    integer(4) :: i,j,k,m,n
    real(8) :: xta,yta,zta,rad,the,csthe,phi
    real(8) :: fvh
    complex(8) :: fi
    real(8) :: s_plus_re, s_plus_im, s_minus_re, s_minus_im
    integer(4) :: m1, m2
    real(8), parameter :: PI_square = PI**2

    !******* set Ewald parameter
    nhmax=20
    nrmax=20
    wkappa=2.d0
    nhcut=nhmax
    nrcut=nrmax
    !     vol=1.d0
    vol=celly*cellz/cellx**2
    !******* set Ewald parameter
    wkappa_square = wkappa**2
    nhcut_square = nhcut**2
    nrcut_square = nrcut**2

    wk_we = dcmplx(0.d0,0.d0)

    call Smat_1(wkappa,vol,nmax,nhmax,wk_we &
    & ) 

if(bNonNeutrality_tot) wk_we(0,0) = wk_we(0,0) -2.0d0*wkappa/sqrt(pi) !! charged system term

 call Smat_2(wkappa,vol,nmax,nrmax,wk_we,shew1,shew2 &
            ) 

!bQcalc1=.false. !! one call should be sufficient

  end subroutine init_fmm_ewald_charged

! emigrated from TI_module
        pure logical  function IsRealZero(realzero, limit0)
          real(8), intent (in)  ::realzero,limit0
          if(realzero <limit0 .and. realzero > -limit0) then
             IsRealZero  = .TRUE.
          else
             IsRealZero = .FALSE.
          endif
          
        end function IsRealZero

end module fmm_ewald

