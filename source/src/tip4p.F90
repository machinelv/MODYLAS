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
!! \brief Module and subroutines to treat TIP4P water model
!<
!----------------------------------------------------------------------
!>
!! \brief Module to treat TIP4P water model
!! \author Yoshimichi Andoh
!<
module tip4p

 implicit none
 real(8):: dOM=0.150d0  ! default is original tip4p
 real(8):: dOH=0.9572d0
 real(8):: thetaHOH=104.52d0 
 real(8) dOHM, sf_rH, sf_rO
 integer(4) msite
 integer(4) :: type_of_msiteposition=0   !! 0: modylas, 1: gromacs

 contains
!-----------------------------------------------------------------------
 subroutine fmod_set_tip4p_geometry(ivalue)
  use mpi_tool
  implicit none
  integer*4 ivalue
!
  if(    ivalue == 0)then   ! TIP4P original
    dOM=0.150d0
    dOH=0.9572d0
    thetaHOH=104.52d0 
  elseif(ivalue == 1)then   ! TIP4P/2005
    dOM=0.1546d0
    dOH=0.9572d0
    thetaHOH=104.52d0 
  else
    if(myrank==0)then
      write(0,*) 'ERROR in <tip4p> type= </tip4p>'
    endif
    call modylas_abort
  endif

 end subroutine fmod_set_tip4p_geometry
!-----------------------------------------------------------------------
 subroutine fmod_set_msite_position(ivalue)
  use mpi_tool
  implicit none
  integer*4 ivalue

  if(    ivalue == 0)then   ! modylas style
    type_of_msiteposition=0
  elseif(ivalue == 1)then   ! gromacs style
    type_of_msiteposition=1
  else
    if(myrank==0)then
      write(0,*) 'ERROR in <tip4p> msite= </tip4p>'
    endif
    call modylas_abort
  endif

 end subroutine fmod_set_msite_position
!-----------------------------------------------------------------------
 subroutine init_tip4p
  use atom_mass, only : mass
  use param, only : paranum
  use trajectory_org, only : xyz, n
  use mpi_tool, only : myrank
  implicit none
  integer(4) i,ipar
  
  thetaHOH=thetaHOH/180.d0*3.14159265358979d0

  dOHM=dOH*dcos(thetaHOH*0.5d0)
  sf_rH=dOM/(2d0*dOHM)
  sf_rO=(1.d0-dOM/dOHM)

  msite=0
  do i=1,n
    ipar=paranum(i)
    if(mass(ipar)<0) then
      msite = msite + 1
    endif
  enddo

!for debug
  if(myrank==0)then
    write(*,*) 'TIP4P parameters are:'
    write(*,*) 'dOM  =', dOM
    write(*,*) 'dOH  =', dOH
    write(*,*) 'theta=', thetaHOH
!   write(*,*) 'msite=', msite
  endif

!$omp parallel do default(shared) &
!$omp& private(i,ipar)
  do i=1,n
    ipar=paranum(i)
    if(mass(ipar)<0) then
  IF(    type_of_msiteposition==0)THEN  !! modylas
      !i  : M
      !i+1: O
      !i+2: H
      !i+3: H
      xyz(1,i) = xyz(1,i+1)            *sf_rO &
     &         +(xyz(1,i+2)+xyz(1,i+3))*sf_rH
      xyz(2,i) = xyz(2,i+1)            *sf_rO &
     &         +(xyz(2,i+2)+xyz(2,i+3))*sf_rH
      xyz(3,i) = xyz(3,i+1)            *sf_rO &
     &         +(xyz(3,i+2)+xyz(3,i+3))*sf_rH
  ELSEIF(type_of_msiteposition==1)THEN  !! gromacs
      !i-3: O
      !i-2: H
      !i-1: H
      !i  : M
      xyz(1,i) = xyz(1,i-3)            *sf_rO &
     &         +(xyz(1,i-2)+xyz(1,i-1))*sf_rH
      xyz(2,i) = xyz(2,i-3)            *sf_rO &
     &         +(xyz(2,i-2)+xyz(2,i-1))*sf_rH
      xyz(3,i) = xyz(3,i-3)            *sf_rO &
     &         +(xyz(3,i-2)+xyz(3,i-1))*sf_rH
  ENDIF
        endif
      enddo
!$omp end parallel do

 end subroutine init_tip4p
!-----------------------------------------------------------------------
 subroutine update_msite_coordinate()
  use atom_mass, only : mass
  use param, only : paranum
  use subcell, only : nselfseg, lsegtop, lseg_natoms, m2i
  use trajectory_mpi, only : wkxyz
  implicit none
  integer(4) i0,k0,ipar

!$omp parallel do default(shared) &
!$omp& private(k0,i0,ipar)
  do k0=1,nselfseg
!
  if(lseg_natoms(k0)/=4)cycle  !! only check segments with 4 atoms
!
  IF(    type_of_msiteposition==0)THEN  !! modylas
    i0 = lsegtop(k0)
    ipar=paranum(m2i(i0))
    if(mass(ipar) <  0.d0) then
      !i0:   M
      !i0+1: O
      !i0+2: H
      !i0+3: H
      wkxyz(1,i0) = wkxyz(1,i0+1)               *sf_rO &
     &            +(wkxyz(1,i0+2)+wkxyz(1,i0+3))*sf_rH
      wkxyz(2,i0) = wkxyz(2,i0+1)               *sf_rO &
     &            +(wkxyz(2,i0+2)+wkxyz(2,i0+3))*sf_rH
      wkxyz(3,i0) = wkxyz(3,i0+1)               *sf_rO &
     &            +(wkxyz(3,i0+2)+wkxyz(3,i0+3))*sf_rH
    endif
  ELSEIF(type_of_msiteposition==1)THEN  !! gromacs
    i0 = lsegtop(k0)+3
    ipar=paranum(m2i(i0))
    if(mass(ipar) <  0.d0) then
      !i0-3: O
      !i0-2: H
      !i0-1: H
      !i0  : M
      wkxyz(1,i0) = wkxyz(1,i0-3)               *sf_rO &
     &            +(wkxyz(1,i0-2)+wkxyz(1,i0-1))*sf_rH
      wkxyz(2,i0) = wkxyz(2,i0-3)               *sf_rO &
     &            +(wkxyz(2,i0-2)+wkxyz(2,i0-1))*sf_rH
      wkxyz(3,i0) = wkxyz(3,i0-3)               *sf_rO &
     &            +(wkxyz(3,i0-2)+wkxyz(3,i0-1))*sf_rH
    endif
  ENDIF
  enddo
!$omp end parallel do

 end subroutine update_msite_coordinate
!-----------------------------------------------------------------------
 subroutine distribute_force_on_msite(wk_ftmp)
  use atom_mass, only : mass
  use param, only : paranum
  use subcell, only : nselfseg, lsegtop, lseg_natoms, m2i
  use trajectory_mpi, only : nadirect
  implicit none
  integer(4) i0,k0,ipar
  real(8)::wk_ftmp(3,nadirect)

!$omp parallel do default(shared) &
!$omp& private(k0,i0,ipar)
  do k0=1,nselfseg
!
  if(lseg_natoms(k0)/=4)cycle  !! only check segments with 4 atoms
!
  IF(    type_of_msiteposition==0)THEN  !! modylas
     i0 = lsegtop(k0)         
     ipar=paranum(m2i(i0))
     if(mass(ipar) <  0.d0) then
      !i0:   M
      !i0+1: O
      !i0+2: H
      !i0+3: H
       wk_ftmp(:,i0+1)=wk_ftmp(:,i0+1)+wk_ftmp(:,i0)*sf_rO
       wk_ftmp(:,i0+2)=wk_ftmp(:,i0+2)+wk_ftmp(:,i0)*sf_rH
       wk_ftmp(:,i0+3)=wk_ftmp(:,i0+3)+wk_ftmp(:,i0)*sf_rH
     endif
  ELSEIF(type_of_msiteposition==1)THEN  !! gromacs
     i0 = lsegtop(k0)+3
     ipar=paranum(m2i(i0))
     if(mass(ipar) <  0.d0) then
      !i0-3: O
      !i0-2: H
      !i0-1: H
      !i0  : M
       wk_ftmp(:,i0-3)=wk_ftmp(:,i0-3)+wk_ftmp(:,i0)*sf_rO
       wk_ftmp(:,i0-2)=wk_ftmp(:,i0-2)+wk_ftmp(:,i0)*sf_rH
       wk_ftmp(:,i0-1)=wk_ftmp(:,i0-1)+wk_ftmp(:,i0)*sf_rH
     endif
  ENDIF
  enddo 
!$omp end parallel do

 end subroutine distribute_force_on_msite

end module tip4p
