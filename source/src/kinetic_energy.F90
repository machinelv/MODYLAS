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
!! \brief Module and subroutine to calculate the kinetic energy
!<
!----------------------------------------------------------------------
!>
!! \brief Module to calculate the kinetic energy
!! \author Noriyuki Yoshii
!<
module kinetic_energy
  use omp_lib
  use subcell
  use mpi_tool
  use domain, only : lxdiv, lydiv, lzdiv
  use dr_cntl, only : nbd
  implicit none
  real(8),allocatable :: wk_k_ene(:,:)

contains

!>
!! \brief  Subroutine to allocate array.
!! \author Kensuke Iwahashi
!<
  subroutine fmod_alloc_kinetic_energy
    use omp_lib
    implicit none
    integer(4) :: nomp=1
!$  nomp = omp_get_max_threads()
    allocate(wk_k_ene(6,0:nomp-1))
  end subroutine fmod_alloc_kinetic_energy
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to generate initial velocity if needed.
!! \author Kensuke Iwahashi
!<
  subroutine generate_velocity
    use atom_mass
    use trajectory_org
    use trajectory_org
    use md_condition
    use md_const
    use maxwell_distribution
    use param
    use random_seed
#ifdef TIP4
    use tip4p
#endif
    implicit none
    real(8),allocatable :: randbuf(:)
    real(8) :: w1,w2,M,vG__x,vG__y,vG__z,coef,T
    real(8) :: rand0, rand1, genrand_real2
    integer(4) :: i,nrand
    !     generate maxwell distribution of velocity of atom,
    !     in which posibility variable sqrt(m / kB / T) * v(x, y or z)
    !     obeys to normal distribution with
    !      expectation = 0
    !      standard deviation = 1
    !
    !     after random generation of velocity, they are adjusted for
    !     total momentum of the sytem to be 0 and temperature to be
    !     md_generic__maxwell temperature
    !
    !    handling special case: maxwell_temperature == 0
    if (maxwell_temperature < 1.0d-50) then
       v = 0.0d0
       return
    endif
    !     generate uniform random numbers
    if (mod(n,2) == 0) then
       nrand = 3*n
    else
       nrand = 3*n+1
    endif
    allocate(randbuf(nrand))
    call init_genrand(randomseed)
    !     conver uniform random numbers to normal random numbers
    do i=1,nrand-1,2
       rand0 = genrand_real2()
       rand1 = genrand_real2()
       w1=sqrt(-2.0d0*log(rand0))*cos(2.0d0*PI*rand1)
       w2=sqrt(-2.0d0*log(rand1))*cos(2.0d0*PI*rand0)
       randbuf(i+0) = w1
       randbuf(i+1) = w2
    enddo
    !     add velocities
    do i=1,n
#ifdef TIP4
       if(mass(paranum(i)) .lt. 0)then
         v(:,i)=0d0
       endif
       cycle
#endif
       coef = sqrt(md_BOLTZMANN * maxwell_temperature * r_mass(paranum(i)))
       v(1,i) = randbuf(i*3-2) * coef
       v(2,i) = randbuf(i*3-1) * coef
       v(3,i) = randbuf(i*3-0) * coef
    enddo
    deallocate(randbuf)
    !     velocity of center of mass
    M = 0.0d0
    vG__x = 0.0d0;  vG__y = 0.0d0;  vG__z = 0.0d0
    do i=1,n
#ifdef TIP4
       if(mass(paranum(i)) .lt. 0) cycle
#endif
       if (mass(paranum(i)) .lt. 1.0e+10) then
          M = M + mass(paranum(i))
          vG__x = vG__x + v(1,i) * mass(paranum(i))
          vG__y = vG__y + v(2,i) * mass(paranum(i))
          vG__z = vG__z + v(3,i) * mass(paranum(i))
       endif
    enddo
    vG__x = vG__x / M
    vG__y = vG__y / M
    vG__z = vG__z / M
    !     subtracting momentum of the system
    do i=1,n
#ifdef TIP4
       if(mass(paranum(i)) .lt. 0) cycle
#endif
       if (mass(paranum(i)) .lt. 1.0e+10) then
          v(1,i) = v(1,i) - vG__x
          v(2,i) = v(2,i) - vG__y
          v(3,i) = v(3,i) - vG__z
       endif
    enddo
    !     velocity scaling
    T = 0.0d0
    do i=1,n
#ifdef TIP4
       if(mass(paranum(i)) .lt. 0) cycle
#endif
       if (mass(paranum(i)) .lt. 1.0e+10) then
          T = T + mass(paranum(i))*(v(1,i)*v(1,i)+v(2,i)*v(2,i)+v(3,i)*v(3,i))
       endif
    enddo
    T = T * degree_of_freedom_inverse  * rvkbolz
    v = v * sqrt(maxwell_temperature / T)
  end subroutine generate_velocity

!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to initilize velocities.
!! \author Kensuke Iwahashi
!<
  subroutine init_md_velocity
    use maxwell_distribution
    use trajectory_org
    implicit none
    !     generate velocity of atom obeying to Maxwell distribution
    if (maxwell_temperature > 0.0d0 .or. reset_maxwell) then
       call generate_velocity
       call remove_system_momentum
    else
       call remove_system_momentum
    endif
  end subroutine init_md_velocity
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to remove system momentum.
!! \author Kensuke Iwahashi
!<
  subroutine remove_system_momentum
    use atom_mass
    use trajectory_org
    use param
    implicit none
    real(8) :: M, vG__x, vG__y, vG__z
    integer(4) :: i
    !     velocity of center of mass
    M = 0.0d0
    vG__x = 0.0d0
    vG__y = 0.0d0
    vG__z = 0.0d0
    do i=1, n
#ifdef TIP4
       if(mass(paranum(i)) .lt. 0) cycle
#endif
       M = M + mass(paranum(i))
       vG__x = vG__x + v(1,i) * mass(paranum(i))
       vG__y = vG__y + v(2,i) * mass(paranum(i))
       vG__z = vG__z + v(3,i) * mass(paranum(i))
    enddo
    vG__x = vG__x / M
    vG__y = vG__y / M
    vG__z = vG__z / M
    !     subtracting momentum of the system
    do i=1, n
#ifdef TIP4
       if(mass(paranum(i)) .lt. 0)then
         v(:,i)=0d0  ! set zero for M-site, explicitly
         cycle
       endif
#endif
       v(1,i) = v(1,i) - vG__x
       v(2,i) = v(2,i) - vG__y
       v(3,i) = v(3,i) - vG__z
    enddo
    totalmass=M
  end subroutine remove_system_momentum
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to remove system momentum. (mpi version)
!! \author Yoshimichi Andoh
!<
  subroutine remove_system_momentum_para
!-----------------------------------------------------------------------
    use atom_mass
    use param
    use trajectory_mpi, only : wkv
    use mpi_tool
    implicit none
    real(8) :: wk_vGx,wk_vGy,wk_vGz
    real(8) :: wk_vG(3),vG(3)
    integer(4) :: i0,k0,ipar
    include 'mpif.h'
    integer(4) :: ierr

    !     velocity of center of mass
    wk_vGx = 0.0d0
    wk_vGy = 0.0d0
    wk_vGz = 0.0d0
!$omp parallel do default(none) &
!$omp& private(i0,k0,ipar) &
!$omp& shared(nselfseg,lsegtop,lseg_natoms,paranum) &
!$omp& shared(m2i,mass,wkv) &
!$omp& reduction(+:wk_vGx,wk_vGy,wk_vGz)
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          ipar=paranum(m2i(i0))
#ifdef TIP4
          if(mass(ipar) < 0) cycle
#endif
          wk_vGx=wk_vGx+wkv(1,i0)*mass(ipar)
          wk_vGy=wk_vGy+wkv(2,i0)*mass(ipar)
          wk_vGz=wk_vGz+wkv(3,i0)*mass(ipar)
       enddo ! i0
    enddo ! k0

    wk_vG(1)=wk_vGx
    wk_vG(2)=wk_vGy
    wk_vG(3)=wk_vGz

    call mpi_allreduce(wk_vG,vG,3, mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

    vG = vG / totalmass

    !     subtracting momentum of the system
!$omp parallel do default(none) &
!$omp& private(k0,i0,ipar) &
!$omp& shared(nselfseg,lsegtop,lseg_natoms) &
!$omp& shared(wkv,vG,mass,paranum,m2i)
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
#ifdef TIP4
          ipar=paranum(m2i(i0))
          if(mass(ipar) < 0) cycle
#endif
          wkv(1,i0) = wkv(1,i0) - vG(1)
          wkv(2,i0) = wkv(2,i0) - vG(2)
          wkv(3,i0) = wkv(3,i0) - vG(3)
       enddo
    enddo
  end subroutine remove_system_momentum_para
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to apply velocity scaling
!! \author Yoshimichi Andoh (original by Masugata)
!<
!----------------------------------------------------------------------
  subroutine apply_velocity_scaling
    use atom_mass
    use trajectory_org
    use trajectory_mpi
    use extended_system
    use md_condition
    use md_const
    use param
    use shake_rattle_roll, only : l0max, ibseL
    implicit none
    include 'mpif.h'
    integer(4) :: ierr
    integer(4) :: ii, ipar
    integer(4) :: i0,icx0,icy0,icz0,icxyz0
    real(8) :: sc,tot,wkt
    integer(4) :: degree_of_freedom_local
    integer(4) :: nshake_local
                                
    tot = 0.0d0                 
    wkt = 0.0d0
    
    degree_of_freedom_local=0
!$omp parallel do default(shared) &
!$omp& private(ii,i0,icx0,icy0,icz0,icxyz0,ipar) &
!$omp& reduction(+:wkt,degree_of_freedom_local)
    do ii=1,lxdiv*lydiv*lzdiv
!default: MTD
       icz0=mod(ii-1,lzdiv)     +nbd
       icy0=mod(ii-1,lzdiv*lydiv)
       icy0=icy0/lzdiv          +nbd
       icx0=(ii-1)/(lzdiv*lydiv)+nbd
       do i0=tag(icz0,icy0,icx0), &
            &        tag(icz0,icy0,icx0)+na_per_cell(icz0,icy0,icx0)-1
          ipar=paranum(m2i(i0))
#ifdef TIP4
       if(mass(ipar) .lt. 0) cycle
#endif
          wkt = wkt + mass(ipar)*( &
               &      wkv(1,i0)*wkv(1,i0)+wkv(2,i0)*wkv(2,i0)+wkv(3,i0)*wkv(3,i0))
          degree_of_freedom_local = degree_of_freedom_local + 3
       enddo ! i0
    enddo ! ii
    if(velocity_scaling_region_wise) then
       tot = wkt;
       ! no communication necessary 
    else
       call mpi_allreduce(wkt,tot,1, mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    endif

    if(velocity_scaling_region_wise) then
       nshake_local = sum(ibseL(1:l0max)) 
       !!for debug only
       !  call mpi_allreduce(nshake_local,nshake_total,1, mpi_integer, mpi_sum,mpi_comm_world,ierr)
       !  if (myrank ==  0) write(*,*) nshake_total
       !!end debug

       degree_of_freedom_local = degree_of_freedom_local - nshake_local
       tot=(tot / degree_of_freedom_local) * rvkbolz
       sc  = sqrt(systemp / tot)
    else
       tot = tot * degree_of_freedom_inverse  * rvkbolz
       sc  = sqrt(systemp / tot)
    endif

!$omp parallel do default(shared) &
!$omp& private(ii,i0,icx0,icy0,icz0,icxyz0,ipar)
    do ii=1,lxdiv*lydiv*lzdiv
!default: MTD
       icz0=mod(ii-1,lzdiv)     +nbd
       icy0=mod(ii-1,lzdiv*lydiv)
       icy0=icy0/lzdiv          +nbd
       icx0=(ii-1)/(lzdiv*lydiv)+nbd
       do i0=tag(icz0,icy0,icx0), &
            &        tag(icz0,icy0,icx0)+na_per_cell(icz0,icy0,icx0)-1
#ifdef TIP4
       ipar=paranum(m2i(i0))
       if(mass(ipar) .lt. 0) cycle
#endif
          wkv(1,i0) = wkv(1,i0)*sc
          wkv(2,i0) = wkv(2,i0)*sc
          wkv(3,i0) = wkv(3,i0)*sc
       end do ! i0
    end do ! ii
  end subroutine apply_velocity_scaling
!-----------------------------------------------------------------------
!>
!! \brief Subroutine which calculate the kinetic energy
!! \author Noriyuki Yoshii
!<
  subroutine k_energy_scaler(k_ene_sum)
    use atom_mass
    use trajectory_org
    use trajectory_mpi
    use md_const
    use param
#include "timing.h90"
    implicit none
    include 'mpif.h'
    integer(4) :: ierr
    integer(4) :: ii,i0,ipar,nomp,iam
    integer(4) :: i,j,k
    integer(4) :: icx0,icy0,icz0,icxyz0
    real(8) :: k_ene_sum
    real(8) :: wk_ksum
    real(8) :: wmass

    nomp = 1
!$  nomp = omp_get_max_threads()
    iam = 0

!$omp parallel default(none) &
!$omp& private(iam,ii,i0,j,k,ipar,wmass) &
!$omp& private(icx0,icy0,icz0,icxyz0) &
!$omp& shared(wkv,mass,wk_k_ene,m2i,paranum) &
!$omp& shared(tag,na_per_cell,lxdiv,lydiv,lzdiv)
!$  iam = omp_get_thread_num()
    wk_k_ene(:,iam) = 0.0d0
!$omp do
    do ii=1,lxdiv*lydiv*lzdiv
!default: MTD
       icz0=mod(ii-1,lzdiv)     +nbd
       icy0=mod(ii-1,lzdiv*lydiv)
       icy0=icy0/lzdiv          +nbd
       icx0=(ii-1)/(lzdiv*lydiv)+nbd
       do i0=tag(icz0,icy0,icx0), &
            &      tag(icz0,icy0,icx0)+na_per_cell(icz0,icy0,icx0)-1
          ipar=paranum(m2i(i0))
#ifdef TIP4
       if(mass(ipar) .lt. 0) cycle
#endif
          wmass=mass(ipar)
          wk_k_ene(1,iam)=wk_k_ene(1,iam)+wmass*wkv(1,i0)*wkv(1,i0)
          wk_k_ene(2,iam)=wk_k_ene(2,iam)+wmass*wkv(2,i0)*wkv(2,i0)
          wk_k_ene(3,iam)=wk_k_ene(3,iam)+wmass*wkv(3,i0)*wkv(3,i0)
       enddo ! i0
    enddo ! ii
!$omp end do
!$omp end parallel

    wk_ksum=0d0
    do iam=0,nomp-1
       do i=1,3
          wk_ksum=wk_ksum+0.5d0*wk_k_ene(i,iam)
       enddo
    enddo

    TIME_BARRIER(TMB_ALLREDUCE_THERMO)
    TIME_START(TM_ALLREDUCE_THERMO)
    k_ene_sum=0d0
    call mpi_allreduce(wk_ksum,k_ene_sum,1, &
         &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    TIME_STOP(TM_ALLREDUCE_THERMO)
  end subroutine k_energy_scaler
!-----------------------------------------------------------------------
!>
!! \brief Subroutine which calculate the kinetic energy of rigid segment
!! \author Noriyuki Yoshii
!<
  subroutine k_energy_scaler_seg(ktrans,krot)
    use atom_mass
    use trajectory_org
    use trajectory_mpi
    use md_const
    use param
    use subcell, only : lsegtop,lseg_natoms, nselfseg
    use matrix_inverse_mod
#include "timing.h90"
    implicit none
    include 'mpif.h'
    integer(4) :: ierr
    integer(4) :: k0,i0,ipar,nomp,iam,k
    integer(4) :: icx0,icy0,icz0,icxyz0
    real(8) :: k_ene_sum, wk_krot
    real(8) :: wk_ksum(2), ksum(2)
    real(8) :: ktrans, krot
    real(8) :: wmass, dvel(3)
    real(8) :: segmom(3),segMa, segtemp, rottemp, tottemp
    real(8) :: segvel(3),angmom(3),dxyz(3),vangl(3)
    real(8) :: itensor(3,3),itinv(3,3)

    nomp = 1
!$  nomp = omp_get_max_threads()
    iam = 0

    wk_krot=0d0

!$omp parallel default(none) &
!$omp& private(iam,k0,i0,ipar,wmass) &
!$omp& private(segmom,segMa,dvel,krot) &
!$omp& private(segvel,angmom,itensor,itinv,dxyz,vangl) &
!$omp& shared(wkv,mass,wk_k_ene,m2i,paranum) &
!$omp& shared(wkxyz,wseg_cx,wseg_cy,wseg_cz) &
!$omp& shared(tag,na_per_cell,lxdiv,lydiv,lzdiv) &
!$omp& shared(lsegtop,lseg_natoms,nselfseg) &
!$omp& reduction(+:wk_krot)
!$  iam = omp_get_thread_num()
    wk_k_ene(:,iam) = 0.0d0
!$omp do
    do k0=1,nselfseg
!
! Translational degree of freedom
!
      segMa=0d0
      segmom(:)=0d0
      do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
        ipar=paranum(m2i(i0))
#ifdef TIP4
        if(mass(ipar) .lt. 0) cycle
#endif
        segMa=segMa+mass(ipar)
        segmom(:)=segmom(:)+mass(ipar)* wkv(:,i0)
      enddo ! i0
      wk_k_ene(1,iam)=wk_k_ene(1,iam)+segmom(1)**2/segMa  ! P*P/M
      wk_k_ene(2,iam)=wk_k_ene(2,iam)+segmom(2)**2/segMa  ! P*P/M
      wk_k_ene(3,iam)=wk_k_ene(3,iam)+segmom(3)**2/segMa  ! P*P/M
!
! Rotational degree of freedom
!
      if(    lseg_natoms(k0) .ge. 3)then

        segvel=segmom/segMa
        angmom=0d0   ! angular momentumn, L
        itensor=0d0  ! inertia tensor,    I
        do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          ipar=paranum(m2i(i0))
          wmass=mass(ipar)
          dxyz(1)=wkxyz(1,i0)-wseg_cx(k0) ! relative coordinate
          dxyz(2)=wkxyz(2,i0)-wseg_cy(k0) ! relative coordinate
          dxyz(3)=wkxyz(3,i0)-wseg_cz(k0) ! relative coordinate
          dvel(1)=wkv(1,i0)  -segvel(1)   ! relative velocity 
          dvel(2)=wkv(2,i0)  -segvel(2)   ! relative velocity 
          dvel(3)=wkv(3,i0)  -segvel(3)   ! relative velocity 
          angmom(1)=angmom(1)  &
                   +wmass*( dxyz(2)*dvel(3)-dxyz(3)*dvel(2) )
          angmom(2)=angmom(2)  &
                   +wmass*( dxyz(3)*dvel(1)-dxyz(1)*dvel(3) )
          angmom(3)=angmom(3)  &
                   +wmass*( dxyz(1)*dvel(2)-dxyz(2)*dvel(1) )
          itensor(1,1)=itensor(1,1)+wmass*( dxyz(2)**2+dxyz(3)**2 )
          itensor(2,2)=itensor(2,2)+wmass*( dxyz(3)**2+dxyz(1)**2 )
          itensor(3,3)=itensor(3,3)+wmass*( dxyz(1)**2+dxyz(2)**2 )
          itensor(1,2)=itensor(1,2)-wmass*( dxyz(1)*   dxyz(2)    )
          itensor(1,3)=itensor(1,3)-wmass*( dxyz(1)*   dxyz(3)    )
          itensor(2,3)=itensor(2,3)-wmass*( dxyz(2)*   dxyz(3)    )
        enddo
        itensor(2,1)=itensor(1,2)
        itensor(3,1)=itensor(1,3)
        itensor(3,2)=itensor(2,3)
        call matrix_inverse(itensor,itinv)  ! inverse of I
        vangl=matmul(itinv,angmom)          ! omega
        krot=dot_product(vangl,angmom)      ! omega*L
        wk_krot=wk_krot+krot

      elseif(lseg_natoms(k0) .eq. 2)then

        i0=lsegtop(k0)
        wmass=mass(paranum(m2i(i0  )))  &
     &       *mass(paranum(m2i(i0+1)))
        dvel(1)=wkv(1,i0) - wkv(1,i0+1)   ! relative velocity 
        dvel(2)=wkv(2,i0) - wkv(2,i0+1)   ! relative velocity 
        dvel(3)=wkv(3,i0) - wkv(3,i0+1)   ! relative velocity 

        krot=wmass/segMa*(dvel(1)**2+dvel(2)**2+dvel(3)**2)
        wk_krot=wk_krot+krot

      elseif(lseg_natoms(k0) .eq. 1)then

        krot=0d0                !! no-rotation
        wk_krot=wk_krot+krot

      endif
    enddo ! k0
!$omp end do nowait
!$omp end parallel

!
! Translational kinetic energy
!
    wk_ksum(1)=0d0
    do iam=0,nomp-1
       do k=1,3
          wk_ksum(1)=wk_ksum(1)+wk_k_ene(k,iam)
       enddo
    enddo
    wk_ksum(1)=0.5d0*wk_ksum(1)
!
! Rotational kinetic energy
!
    wk_ksum(2)=0.5d0*wk_krot
!
! Reduction over MPI processes
!
    TIME_BARRIER(TMB_ALLREDUCE_THERMO)
    TIME_START(TM_ALLREDUCE_THERMO)
    call mpi_allreduce(wk_ksum,ksum,2,  &
         &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    TIME_STOP(TM_ALLREDUCE_THERMO)
    ktrans=ksum(1)
    krot  =ksum(2)

  end subroutine k_energy_scaler_seg
!----------------------------------------------------------------------
!>
!! \brief Subroutine which calculate the kinetic energy tensor
!! \author Yoshimichi Andoh, Noriyuki Yoshii
!<
!----------------------------------------------------------------------
  subroutine k_energy_tensor(k_ene_sum,k_ene)
    use atom_mass
    use trajectory_org
    use trajectory_mpi
    use md_const
    use param
    implicit none
    include 'mpif.h'
    integer(4) :: ierr
    real(8) :: k_ene_sum, k_ene(6)
    integer(4) :: ii,i0,ipar,nomp,iam
    integer(4) :: i,j,k
    integer(4) :: icx0,icy0,icz0,icxyz0
    real(8) :: wmass
    real(8) :: wk_ksum(6)

    nomp = 1
    !$    nomp = omp_get_max_threads()
    iam = 0

!$omp parallel default(none) &
!$omp& private(iam,ii,i0,j,k,ipar,wmass) &
!$omp& private(icx0,icy0,icz0,icxyz0) &
!$omp& shared(wkv,mass,wk_k_ene,m2i,paranum) &
!$omp& shared(tag,na_per_cell,lxdiv,lydiv,lzdiv)
!$  iam = omp_get_thread_num()
    wk_k_ene(:,iam) = 0.0d0
!$omp do
    do ii=1,lxdiv*lydiv*lzdiv
!default: MTD
       icz0=mod(ii-1,lzdiv)     +nbd
       icy0=mod(ii-1,lzdiv*lydiv)
       icy0=icy0/lzdiv          +nbd
       icx0=(ii-1)/(lzdiv*lydiv)+nbd
       do i0=tag(icz0,icy0,icx0), &
            &      tag(icz0,icy0,icx0)+na_per_cell(icz0,icy0,icx0)-1
          ipar=paranum(m2i(i0))
#ifdef TIP4
       if(mass(ipar) .lt. 0) cycle
#endif
          wmass=mass(ipar)
          wk_k_ene(1,iam)=wk_k_ene(1,iam)+wmass*wkv(1,i0)*wkv(1,i0)
          wk_k_ene(2,iam)=wk_k_ene(2,iam)+wmass*wkv(2,i0)*wkv(2,i0)
          wk_k_ene(3,iam)=wk_k_ene(3,iam)+wmass*wkv(3,i0)*wkv(3,i0)
          wk_k_ene(4,iam)=wk_k_ene(4,iam)+wmass*wkv(2,i0)*wkv(1,i0)
          wk_k_ene(5,iam)=wk_k_ene(5,iam)+wmass*wkv(3,i0)*wkv(1,i0)
          wk_k_ene(6,iam)=wk_k_ene(6,iam)+wmass*wkv(3,i0)*wkv(2,i0)
       enddo ! i0
    enddo ! ii
!$omp end do
!$omp end parallel

    wk_ksum(:)=0d0
    do iam=0,nomp-1
       wk_ksum(:)=wk_ksum(:)+0.5d0*wk_k_ene(:,iam)
    enddo

    k_ene=0d0
    call mpi_allreduce(wk_ksum,k_ene,6, &
         &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    k_ene_sum=k_ene(1)+k_ene(2)+k_ene(3)
  end subroutine k_energy_tensor
!----------------------------------------------------------------------
!>
!! \brief Subroutine which calculate the barostat's kinetic energy
!! \author Noriyuki Yoshii
!----------------------------------------------------------------------
  subroutine k_energy_baro(vg,bmass,k_ene_baro)
    implicit none
    real(8) :: vg(3,3),tvg(3,3),prod(3,3)
    real(8) :: trc,k_ene_baro,bmass
    integer(4) :: i,j

    do i=1,3
       do j=1,3
          tvg(i,j)=vg(j,i)
       enddo
    enddo

    prod = matmul(tvg,vg)

    trc=prod(1,1)+prod(2,2)+prod(3,3)
    k_ene_baro = bmass*trc
  end subroutine k_energy_baro

end module kinetic_energy

