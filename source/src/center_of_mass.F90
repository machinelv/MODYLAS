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
!! \brief Module and subroutines which apply SHAKE/RATTLE method to constrain
!!        the distance between two mass centers.
!<
!----------------------------------------------------------------------
!>
!! \brief Module which apply SHAKE/RATTLE method to constrain
!!        the distance between two mass centers.
!! \author Yoshimichi Andoh
!<
module center_of_mass
  implicit none
  real(8) :: CenterOfMassA(3)=0.0d0, CenterOfMassB(3)=0.0d0
  real(8) :: CenterOfMassVeloA(3)=0.0d0, CenterOfMassVeloB(3)=0.0d0
  real(8) :: CenterOfMassAOld(3)=0.0d0, CenterOfMassBOld(3)=0.0d0
  real(8) :: TotalMassA=0.0d0, TotalMassB=0.0d0, mu
  real(8) :: RTotalMassA=0.0d0, RTotalMassB=0.0d0

  integer(4)::i0ownA,i0ownB
  integer(4),allocatable::i0listA(:),i0listB(:)

contains

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to initialize variables.
!! \author Yoshimichi Andoh, originated by K.Fujimoto
!<
  subroutine  init_center_of_mass
    use center_of_mass_variables
    use segments
    use atom_mass
    use param
    use trajectory_org
    use mpi_tool
    use unit_cell
    implicit none
    integer(4) :: i, j, k, i0
    integer(4) :: ipar
    
    !変数の初期化
    TotalMassA = 0.0d0
    TotalMassB = 0.0d0
    CenterOfMassA=0d0
    CenterOfMassB=0d0
    i0ownA=0
    i0ownB=0
    
    !グループAの重心の計算
    do i = groupAtop, groupAend
       do j = segtop(i)+1, segtop(i)+seg_natoms(i)
          ipar=paranum(j)
          TotalMassA = TotalMassA + mass(ipar)
          CenterOfMassA(:) = CenterOfMassA(:) + mass(ipar)*xyz(:,j)
          i0ownA=i0ownA+1
       enddo
    enddo
    CenterOfMassA = CenterOfMassA / TotalMassA
      
    !グループBの重心の計算
    do i = groupBtop, groupBend
       do j = segtop(i)+1, segtop(i)+seg_natoms(i)
          ipar=paranum(j)
          TotalMassB = TotalMassB + mass(ipar)
          CenterOfMassB(:) = CenterOfMassB(:) + mass(ipar)*xyz(:,j)
          i0ownB=i0ownB+1
       enddo
    enddo
    CenterOfMassB = CenterOfMassB / TotalMassB
    
    RTotalMassA = 1.0d0/TotalMassA
    RTotalMassB = 1.0d0/TotalMassB
    mu = (TotalMassA+TotalMassB)/(TotalMassA*TotalMassB)
    
    allocate(i0listA(i0ownA),i0listB(i0ownB))
    
    call checkPBCdr(CenterOfMassA,CenterOFMassB)
    !debug
    !     write(*,*) "TotalMassA=", TotalMassA
    !     write(*,*) "TotalMassB=", TotalMassB
    !debug
  end subroutine init_center_of_mass
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to make a list of i0 atoms belonging to myrank
!! \author Yoshimichi Andoh
!<
  subroutine  list_own_i0
    use center_of_mass_variables
    use segments
    use trajectory_mpi
    use subcell
    implicit none
    integer(4)::i,j,i0
    
    i0ownA=0;i0listA=0
    do i = groupAtop, groupAend
       do j = segtop(i)+1, segtop(i)+seg_natoms(i)
          i0=i2m(j)
          if(i0.eq.-1) cycle
          i0ownA=i0ownA+1
          i0listA(i0ownA)=i0
       enddo
    enddo
    
    i0ownB=0;i0listB=0
    do i = groupBtop, groupBend
       do j = segtop(i)+1, segtop(i)+seg_natoms(i)
          i0=i2m(j)
          if(i0.eq.-1) cycle
          i0ownB=i0ownB+1
          i0listB(i0ownB)=i0
       enddo
    enddo
  end subroutine list_own_i0
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate mass centers of two segments
!! \author Yoshimichi Andoh, originated by K.Fujimoto
!<
  subroutine  calc_center_of_mass
    use center_of_mass_variables
    use atom_mass
    use param
    use trajectory_mpi
    use subcell
    use mpi_tool
    implicit none
    include 'mpif.h'
    integer(4) :: i, j, k, i0
    integer(4) :: ipar,icheck,ierr
    real(8) :: wk_COMA(3),wk_COMB(3)
    
    !グループAの重心の計算
    wk_COMA=0d0
    do i=1,i0ownA
       i0=i0listA(i)
       ipar=paranum(m2i(i0))
       wk_COMA(1) = wk_COMA(1) + mass(ipar)*wkxyz(1,i0)
       wk_COMA(2) = wk_COMA(2) + mass(ipar)*wkxyz(2,i0)
       wk_COMA(3) = wk_COMA(3) + mass(ipar)*wkxyz(3,i0)
    enddo
    
    !グループBの重心の計算
    wk_COMB=0d0
    do i=1,i0ownB
       i0=i0listB(i)
       ipar=paranum(m2i(i0))
       wk_COMB(1) = wk_COMB(1) + mass(ipar)*wkxyz(1,i0)
       wk_COMB(2) = wk_COMB(2) + mass(ipar)*wkxyz(2,i0)
       wk_COMB(3) = wk_COMB(3) + mass(ipar)*wkxyz(3,i0)
    enddo
    
    call MPI_Allreduce(wk_COMA, CenterOfMassA,3,MPI_REAL8, &
         &                MPI_SUM,MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(wk_COMB, CenterOfMassB,3,MPI_REAL8, &
         &                MPI_SUM,MPI_COMM_WORLD, ierr)
    
    do i = 1, 3
       CenterOfMassA(i) = CenterOfMassA(i) / TotalMassA
       CenterOfMassB(i) = CenterOfMassB(i) / TotalMassB
    enddo
    
    !     call checkPBCdr(CenterOfMassA,CenterOFMassB)
  end subroutine calc_center_of_mass
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate mass centes's velocities of two segments
!! \author Yoshimichi Andoh, originated by K.Fujimoto
!<
  subroutine  calc_center_of_velocity
    use center_of_mass_variables
    use atom_mass
    use param
    use trajectory_mpi
    use subcell
    use mpi_tool
    implicit none
    include 'mpif.h'
    integer(4) :: i, j, k, i0
    integer(4) :: ipar,icheck,ierr
    real(8) :: wk_VCOMA(3),wk_VCOMB(3)
    
    !グループAの重心の計算
    wk_VCOMA=0d0
    do i=1,i0ownA
       i0=i0listA(i)
       ipar=paranum(m2i(i0))
       wk_VCOMA(1) = wk_VCOMA(1) + mass(ipar)*wkv(1,i0)
       wk_VCOMA(2) = wk_VCOMA(2) + mass(ipar)*wkv(2,i0)
       wk_VCOMA(3) = wk_VCOMA(3) + mass(ipar)*wkv(3,i0)
    enddo
    
    !グループBの重心の計算
    do i=1,i0ownB
       i0=i0listB(i)
       ipar=paranum(m2i(i0))
       wk_VCOMB(1) = wk_VCOMB(1) + mass(ipar)*wkv(1,i0)
       wk_VCOMB(2) = wk_VCOMB(2) + mass(ipar)*wkv(2,i0)
       wk_VCOMB(3) = wk_VCOMB(3) + mass(ipar)*wkv(3,i0)
    enddo
    
    call MPI_Allreduce(wk_VCOMA,CenterOfMassVeloA,3,MPI_REAL8, &
         &                MPI_SUM,MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(wk_VCOMB,CenterOfMassVeloB,3,MPI_REAL8, &
         &                MPI_SUM,MPI_COMM_WORLD, ierr)
    do i = 1, 3
       CenterOfMassVeloA(i) = CenterOfMassVeloA(i) / TotalMassA
       CenterOfMassVeloB(i) = CenterOfMassVeloB(i) / TotalMassB
    enddo
    
    !     call checkPBCdr(CenterOfMassA,CenterOFMassB)  !! check at t+dt
  end subroutine calc_center_of_velocity
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate radial distance between two COMs,
!!         under periodic boundary condition.
!! \author Yoshimichi Andoh, originated by K.Fujimoto
!<
  subroutine checkPBCdr(COA,COB)
    use center_of_mass_variables
    use atom_mass
    use param
    use trajectory_org
    use boundary
    use mpi_tool
    implicit none
    real(8)::COA(3),COB(3)
    real(8)::dx,dy,dz,dr
    
    dx=COA(1)-COB(1)
    dy=COA(2)-COB(2)
    dz=COA(3)-COB(3)
    
    call PBCij(dx,dy,dz)
    
    dr=sqrt(dx**2+dy**2+dz**2)
    if(myrank==0)then
       write(*,'(1x,a44,f12.5,a4)') &
            &  'Distance between groupA_COM and groupB_COM =',&
            &  dr*1d+10,' [A]'
    endif
  end subroutine checkPBCdr

end module center_of_mass
