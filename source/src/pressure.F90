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
!! \brief  Module and subroutines to calculate atomic virial.
!<
!! [Ref] Y Andoh, N Yoshii, A Yamada, S Okazaki,
!!       J. Comput. Chem., 38, 704-713 (2017).
!----------------------------------------------------------------------
!>
!! \brief  Subroutines to calculate atomic virial
!! \author Yoshimichi Andoh
!<
module pressure_mod

#include "timing.h90"
  use omp_lib
  use mpi_tool
  implicit none

contains

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate atomic virial to be outputed
!! \author Yoshimichi Andoh
!<
  subroutine calc_scalervirial()
!----------------------------------------------------------------------
    use shake_rattle_roll
    use atom_virial
    use md_condition
    use lj_mod
    use md_multiplestep
    use unit_cell
    implicit none
    integer(4) :: i,k,j,l,ibase
    real(8) :: Wc, tmp
    real(8) :: tmpv(3)
    real(8),allocatable :: fS(:,:),fR(:,:)
    real(8),allocatable :: gS(:),gR(:)
    real(8) :: ffS(3),ffR(3)
    real(8) :: ggS,ggR
    include 'mpif.h'
    integer(4) :: ierr
    integer(4) :: nomp

    nomp = 1
!$  nomp = omp_get_max_threads()

    allocate( fS(3,maxMTl),fR(3,maxMTl) )
    allocate( gS(maxMTl),gR(maxMTl) )

!### virial from potential ###
    DO iMT=1,maxMTl

       ibase=2*maxMTm*(iMT-1) !! baseline for wkvirSca
       !=== f^S(0) ===!
       tmpv=0d0
       do k=1,maxMTm
          do i=1,k
             tmpv(1)=tmpv(1)+wkvirScaS(2*(i-1)+1+ibase)
             tmpv(2)=tmpv(2)+wkvirScaM(2*(i-1)+1+ibase)
             tmpv(3)=tmpv(3)+wkvirScaL(2*(i-1)+1+ibase)
          enddo
       enddo
       do k=1,maxMTm-1
          do i=1,k
             tmpv(1)=tmpv(1)+wkvirScaS(2*(i-1)+2+ibase)
             tmpv(2)=tmpv(2)+wkvirScaM(2*(i-1)+2+ibase)
             tmpv(3)=tmpv(3)+wkvirScaL(2*(i-1)+2+ibase)
          enddo
       enddo
       fS(:,iMT)=tmpv(:)/dble(maxMTm)**2
       !=== f^R(dt) ===!
       tmpv(:)=0d0
       do i=1,maxMTm
          tmpv(1) = tmpv(1) + wkvirScaS(2*(i-1)+1+ibase) + wkvirScaS(2*(i-1)+2+ibase)
          tmpv(2) = tmpv(2) + wkvirScaM(2*(i-1)+1+ibase) + wkvirScaM(2*(i-1)+2+ibase)
          tmpv(3) = tmpv(3) + wkvirScaL(2*(i-1)+1+ibase) + wkvirScaL(2*(i-1)+2+ibase)
       enddo
       fR(:,iMT)=tmpv(:)/dble(maxMTm)
       fR(:,iMT)=fR(:,iMT)-fS(:,iMT)

    ENDDO

    !=== F^S(0) ===!
    tmpv(:)=0d0
    do L=1,maxMTl
       do j=1,L
          tmpv(:)=tmpv(:)+fS(:,j)
       enddo
    enddo
    do L=1,maxMTl-1
       do j=1,L
          tmpv(:)=tmpv(:)+fR(:,j)
       enddo
    enddo
    ffS(:)=tmpv(:)/dble(maxMTl)**2
    !=== G^R(Dt) ===!
    tmpv(:)=0d0
    do j=1,maxMTl
       tmpv(:)=tmpv(:)+fS(:,j)+fR(:,j)
    enddo
    ffR(:)=tmpv(:)/dble(maxMTl)
    ffR(:)=ffR(:)-ffS(:)

    wkvir_part(1)=ffR(1)
    wkvir_part(2)=ffR(2)
    wkvir_part(3)=ffR(3)

    !### virial from constraint ###
    IF(totnconst > 0)THEN

       DO iMT=1,maxMTl

          ibase=2*maxMTm*(iMT-1) !! baseline for wkdvirSca
          !=== g^S(0) ===!
          tmp=0d0
          do k=1,maxMTm
             do i=1,k
                tmp=tmp+wkdvirSca(2*(i-1)+1+ibase)
             enddo
          enddo
          do k=1,maxMTm-1
             do i=1,k
                tmp=tmp+wkdvirSca(2*(i-1)+2+ibase)
             enddo
          enddo
          gS(iMT)=tmp/dble(maxMTm)**2
          !=== g^R(dt) ===!
          tmp=0d0
          do i=1,maxMTm
             tmp = tmp + wkdvirSca(2*(i-1)+1+ibase) + wkdvirSca(2*(i-1)+2+ibase)
          enddo
          gR(iMT)=tmp/dble(maxMTm)
          gR(iMT)=gR(iMT)-gS(iMT)

       ENDDO

       !=== G^S(0) ===!
       tmp=0d0
       do L=1,maxMTl
          do j=1,L
             tmp=tmp+gS(j)
          enddo
       enddo
       do L=1,maxMTl-1
          do j=1,L
             tmp=tmp+gR(j)
          enddo
       enddo
       ggS=tmp/dble(maxMTl)**2
       !=== G^R(Dt) ===!
       tmp=0d0
       do j=1,maxMTl
          tmp=tmp+gS(j)+gR(j)
       enddo
       ggR=tmp/dble(maxMTl)
       ggR=ggR-ggS

       wkvir_part(1)=wkvir_part(1)+ggR !! add RATTLE virial

    ENDIF
90  continue

    !### reduce virial over process ###
    vir_part(:)=0d0
    call mpi_allreduce(wkvir_part,vir_part,3, &
         &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

    if(LJ_LRC)then  ! add LJ-LRC to pressure
       Wc = corrector_virial / cellvol
       vir_part(1)=vir_part(1)+3d0*Wc
    endif

    virialSca=vir_part(1)+vir_part(2)+vir_part(3)  !! outputted in mdmntr

    !### summation partial virial ###
    sum_vir_part(:)=sum_vir_part(:)+vir_part(:)
    icount_vir=icount_vir+1  !! counter of virial sumation

    deallocate( fS,fR )
    deallocate( gS,gR )
  end subroutine calc_scalervirial
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate atomic virial at inner loop
!! \author Yoshimichi Andoh
!<
  subroutine calc_scalervirial_inner() ! used in npt_a_integrate,only
!----------------------------------------------------------------------
    use shake_rattle_roll
    use atom_virial
    use md_condition
    use lj_mod
    use md_multiplestep
    use mpi_tool
    use unit_cell
    implicit none
    integer(4) :: i
    real(8) :: Wc
    include 'mpif.h'
    integer(4) :: ierr
    integer(4) :: nomp

    nomp = 1
!$  nomp = omp_get_max_threads()

    wkvirialSca=0d0

    !### virial from potential (short) ###
    do i=1,3
       wkvirialSca=wkvirialSca+virshort(i)
    enddo

    !### virial from constraint (short) ###
    if(totnconst > 0)then
       wkvirialSca=wkvirialSca+wkdvirSca(iMT)  !! rattle only
    endif

    !### reduce virial over process ###
    virialSca=0d0
    TIME_BARRIER(TMB_ALLREDUCE_SCALERVIRIAL)
    TIME_START(TM_ALLREDUCE_SCALERVIRIAL)
    call mpi_allreduce(wkvirialSca,virialSca,1, &
         &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    TIME_STOP(TM_ALLREDUCE_SCALERVIRIAL)

    if(LJ_LRC)then  ! add LJ-LRC to pressure (short)
       Wc = corrector_virial / cellvol
       virialSca=virialSca+3d0*Wc
    endif

  end subroutine calc_scalervirial_inner
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate atomic virial
!! \author Yoshimichi Andoh
!<
  subroutine update_wkvirSca()
!----------------------------------------------------------------------
    use atom_virial
    use md_multiplestep
    implicit none
    integer(4)::i

    !### store virial from potential ###
    wkvirScaS(iMT)=0d0
    wkvirScaM(iMT)=0d0
    wkvirScaL(iMT)=0d0
    do i=1,3 ! x,y,z
       wkvirScaS(iMT) = wkvirScaS(iMT) + virshort(i)
       wkvirScaM(iMT) = wkvirScaM(iMT) + virmiddle(i)*scaleM*maxMTm
       wkvirScaL(iMT) = wkvirScaL(iMT) + virlong(i)  *scaleL*maxMTm*maxMTl
    enddo
  end subroutine update_wkvirSca
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate atomic virial by constraint force
!! \author Yoshimichi Andoh
!<
  subroutine update_wkdvirSca()
!----------------------------------------------------------------------
    use shake_rattle_roll
    use atom_virial
    use md_multiplestep
    implicit none
    integer(4) :: iam,nomp

    nomp=1
!$  nomp=omp_get_max_threads()

    do iam = 0,nomp-1
       wkdvirSca(iMT) = wkdvirSca(iMT) + wkdvir2(1,iam)+wkdvir2(2,iam)+wkdvir2(3,iam)
    enddo
  end subroutine update_wkdvirSca
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate atomic virial tensor to be outputed
!! \author Yoshimichi Andoh
!<
  subroutine calc_virialtensor()
!----------------------------------------------------------------------
![Ref] Y Andoh, N Yoshii, A Yamada, S Okazaki,
!      J. Comput. Chem., 38, 704-713 (2017).
    use shake_rattle_roll
    use atom_virial
    use md_condition
    use md_multiplestep
    use unit_cell
    use lj_mod, only : corrector_virial
    implicit none
    integer(4) :: i,k,j,l,ibase
    real(8) :: Wc, tmp(6)
    real(8) :: tmpv(6,3)
    real(8),allocatable :: fS(:,:,:),fR(:,:,:)
    real(8),allocatable :: gS(:,:),gR(:,:)
    real(8) :: ffS(6,3),ffR(6,3)
    real(8) :: ggS(6),ggR(6)
    include 'mpif.h'
    integer(4) :: ierr
    integer(4) :: nomp

    nomp = 1
!$  nomp = omp_get_max_threads()

    allocate( fS(6,3,maxMTl),fR(6,3,maxMTl) )
    allocate( gS(6,maxMTl),gR(6,maxMTl) )

    !### virial from potential ###
    DO iMT=1,maxMTl

       ibase=2*maxMTm*(iMT-1) !! baseline for wkvirTen
       !=== f^S(0) ===!
       tmpv=0d0
       do k=1,maxMTm
          do i=1,k
             tmpv(:,1)=tmpv(:,1)+wkvirTenS(:,2*(i-1)+1+ibase)
             tmpv(:,2)=tmpv(:,2)+wkvirTenM(:,2*(i-1)+1+ibase)
             tmpv(:,3)=tmpv(:,3)+wkvirTenL(:,2*(i-1)+1+ibase)
          enddo
       enddo
       do k=1,maxMTm-1
          do i=1,k
             tmpv(:,1)=tmpv(:,1)+wkvirTenS(:,2*(i-1)+2+ibase)
             tmpv(:,2)=tmpv(:,2)+wkvirTenM(:,2*(i-1)+2+ibase)
             tmpv(:,3)=tmpv(:,3)+wkvirTenL(:,2*(i-1)+2+ibase)
          enddo
       enddo
       fS(1:6,1:3,iMT)=tmpv(1:6,1:3)/dble(maxMTm)**2
       !=== f^R(dt) ===!
       tmpv=0d0
       do i=1,maxMTm
          tmpv(:,1) = tmpv(:,1) +wkvirTenS(:,2*(i-1)+1+ibase) + wkvirTenS(:,2*(i-1)+2+ibase)
          tmpv(:,2) = tmpv(:,2) + wkvirTenM(:,2*(i-1)+1+ibase) + wkvirTenM(:,2*(i-1)+2+ibase)
          tmpv(:,3) = tmpv(:,3) + wkvirTenL(:,2*(i-1)+1+ibase) + wkvirTenL(:,2*(i-1)+2+ibase)
       enddo
       fR(1:6,1:3,iMT)=tmpv(1:6,1:3)/dble(maxMTm)
       fR(1:6,1:3,iMT)=fR(1:6,1:3,iMT)-fS(1:6,1:3,iMT)

    ENDDO

    !=== F^S(0) ===!
    tmpv=0d0
    do L=1,maxMTl
       do j=1,L
          tmpv(1:6,1:3)=tmpv(1:6,1:3)+fS(1:6,1:3,j)
       enddo
    enddo
    do L=1,maxMTl-1
       do j=1,L
          tmpv(1:6,1:3)=tmpv(1:6,1:3)+fR(1:6,1:3,j)
       enddo
    enddo
    ffS(1:6,1:3)=tmpv(1:6,1:3)/dble(maxMTl)**2
    !=== G^R(Dt) ===!
    tmpv=0d0
    do j=1,maxMTl
       tmpv(1:6,1:3)=tmpv(1:6,1:3)+fS(1:6,1:3,j)+fR(1:6,1:3,j)
    enddo
    ffR(1:6,1:3)=tmpv(1:6,1:3)/dble(maxMTl)
    ffR(1:6,1:3)=ffR(1:6,1:3)-ffS(1:6,1:3)

    wkvir_partTen(:,1)=ffR(:,1) ! short
    wkvir_partTen(:,2)=ffR(:,2) ! middle
    wkvir_partTen(:,3)=ffR(:,3) ! long

    !### virial from constraint ###
    IF(totnconst > 0)THEN

       DO iMT=1,maxMTl

          ibase=2*maxMTm*(iMT-1) !! baseline for wkdvirSca
          !=== g^S(0) ===!
          tmp=0d0
          do k=1,maxMTm
             do i=1,k
                tmp(:)=tmp(:)+wkdvirTen(:,2*(i-1)+1+ibase)
             enddo
          enddo
          do k=1,maxMTm-1
             do i=1,k
                tmp(:)=tmp(:)+wkdvirTen(:,2*(i-1)+2+ibase)
             enddo
          enddo
          gS(:,iMT)=tmp(:)/dble(maxMTm)**2
          !=== g^R(dt) ===!
          tmp=0d0
          do i=1,maxMTm
             tmp(:) = tmp(:) + wkdvirTen(:,2*(i-1)+1+ibase) + wkdvirTen(:,2*(i-1)+2+ibase)
          enddo
          gR(:,iMT)=tmp(:)/dble(maxMTm)
          gR(:,iMT)=gR(:,iMT)-gS(:,iMT)

       ENDDO

       !=== G^S(0) ===!
       tmp=0d0
       do L=1,maxMTl
          do j=1,L
             tmp(:)=tmp(:)+gS(:,j)
          enddo
       enddo
       do L=1,maxMTl-1
          do j=1,L
             tmp(:)=tmp(:)+gR(:,j)
          enddo
       enddo
       ggS(:)=tmp(:)/dble(maxMTl)**2
       !=== G^R(Dt) ===!
       tmp=0d0
       do j=1,maxMTl
          tmp(:)=tmp(:)+gS(:,j)+gR(:,j)
       enddo
       ggR(:)=tmp(:)/dble(maxMTl)
       ggR(:)=ggR(:)-ggS(:)

       wkvir_partTen(:,1)=wkvir_partTen(:,1)+ggR(:) !! add RATTLE virial

    ENDIF

    !### reduce virial over process ###
    vir_partTen=0d0
    call mpi_allreduce(wkvir_partTen,vir_partTen,6*3, &
         &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

    if(LJ_LRC)then  ! add LJ-LRC to pressure
       Wc = corrector_virial / cellvol
       vir_partTen(1,1)=vir_partTen(1,1)+1d0*Wc  ! Pxx
       vir_partTen(2,1)=vir_partTen(2,1)+1d0*Wc  ! Pyy
       vir_partTen(3,1)=vir_partTen(3,1)+1d0*Wc  ! Pzz
    endif

    virialTen(:)=vir_partTen(:,1)+vir_partTen(:,2)+vir_partTen(:,3)  !! outputted in mdmntr
    virialSca=virialTen(1)+virialTen(2)+virialTen(3)

    !### summation partial virial ###
    vir_part(1)=vir_partTen(1,1)+vir_partTen(2,1)+vir_partTen(3,1)
    vir_part(2)=vir_partTen(1,2)+vir_partTen(2,2)+vir_partTen(3,2)
    vir_part(3)=vir_partTen(1,3)+vir_partTen(2,3)+vir_partTen(3,3)
    sum_vir_part(:)=sum_vir_part(:)+vir_part(:)
    icount_vir=icount_vir+1  !! counter of virial sumation

    deallocate( fS,fR )
    deallocate( gS,gR )
  end subroutine calc_virialtensor
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate atomic virial tensor at inner loop
!! \author Yoshimichi Andoh
!<
  subroutine calc_virialtensor_inner() ! used in npt_a_integrate,only
!----------------------------------------------------------------------
    use shake_rattle_roll
    use atom_virial
    use md_condition
    use md_multiplestep
    use unit_cell
    use lj_mod, only : corrector_virial
    implicit none
    integer(4) :: i
    real(8) :: Wc
    include 'mpif.h'
    integer(4) :: ierr
    integer(4) :: nomp

    nomp = 1
!$  nomp = omp_get_max_threads()

    !### virial from potential (short) ###
    wkvirialTen(1:6)=virshort(1:6)

    !### virial from constraint (short) ###
    if(totnconst > 0)then
       wkvirialTen=wkvirialTen+wkdvirTen(:,iMT)  !! rattle only
    endif

    !### reduce virial over process ###
    virialTen=0d0
    call mpi_allreduce(wkvirialTen,virialTen,6, &
         &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

    if(LJ_LRC)then  ! add LJ-LRC to pressure (short)
       Wc = corrector_virial / cellvol
       virialTen(1)=virialTen(1)+1d0*Wc ! Pxx
       virialTen(2)=virialTen(2)+1d0*Wc ! Pyy
       virialTen(3)=virialTen(3)+1d0*Wc ! Pzz
    endif
  end subroutine calc_virialtensor_inner
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate atomic virial
!! \author Yoshimichi Andoh
!<
  subroutine update_wkvirTen()
!----------------------------------------------------------------------
    use atom_virial
    use md_multiplestep
    implicit none

    !### store virial from potential ###
    wkvirTenS(:,iMT)=virshort(:)
    wkvirTenM(:,iMT)=virmiddle(:)*scaleM*maxMTm
    wkvirTenL(:,iMT)=virlong(:)  *scaleL*maxMTm*maxMTl
  end subroutine update_wkvirTen
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate atomic virial by constraint force
!! \author Yoshimichi Andoh
!<
  subroutine update_wkdvirTen()
!----------------------------------------------------------------------
    use shake_rattle_roll
    use atom_virial
    use md_multiplestep
    implicit none
    integer(4) :: iam,nomp

    nomp=1
!$  nomp=omp_get_max_threads()

    do iam = 0,nomp-1
       wkdvirTen(:,iMT)=wkdvirTen(:,iMT)+wkdvir2(:,iam)
    enddo
  end subroutine update_wkdvirTen
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to make branch update_wkvir scaler or tensor
!! \author Yoshimichi Andoh
!<
  subroutine update_wkvir()
!----------------------------------------------------------------------
    implicit none
    call update_wkvirTen()
  end subroutine update_wkvir
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to make branch update_wkdvir scaler or tensor
!! \author Yoshimichi Andoh
!<
  subroutine update_wkdvir()
!----------------------------------------------------------------------
    implicit none
    call update_wkdvirTen()
  end subroutine update_wkdvir
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to make branch update_wkdvir scaler or tensor
!! \author Yoshimichi Andoh
!<
  subroutine init_wkdvir(iMT)
!----------------------------------------------------------------------
    use shake_rattle_roll
    implicit none
    integer(4)::iMT
    wkdvirTen(1:6,iMT)=0d0
  end subroutine init_wkdvir

end module pressure_mod
