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
!! \brief  Module and subroutines to calculate angle potential.
!<
!>
!! \brief  Module to parse parameters of angle potential.
!! \author Kensuke Iwahashi
!<
!----------------------------------------------------------------------
module parse_angle
  use mpi_tool
  implicit none

  type yyparse_angleA
     integer(4) :: atom2, atom3
     real(8) :: Ktheta, theta0
     type(yyparse_angleA),pointer :: next
  end type yyparse_angleA

  type yyparse_angleB
     integer(4) :: atom1, atom3
     real(8) :: Ktheta, theta0
     type(yyparse_angleB),pointer :: next
  end type yyparse_angleB

  type(yyparse_angleA),pointer :: yyparse_angleA_top(:)
  type(yyparse_angleB),pointer :: yyparse_angleB_top(:)
  integer(4),allocatable :: yyparse_angleA_n(:)
  integer(4),allocatable :: yyparse_angleB_n(:)
  integer(4) :: yyparse_angleA_total=0
  integer(4) :: yyparse_angleB_total=0

contains
!----------------------------------------------------------------------
!>
!! \brief  Subroutines to read inputs relating to angle potential.
!! \author Kensuke Iwahashi
!<
  subroutine fmod_alloc_yyparse_angleA_top(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    integer(4) :: i
    type(yyparse_angleA),pointer :: top
    allocate(yyparse_angleA_top(ivalue))
    do i=1,ivalue
       allocate(top)
       top%atom2 = -1
       top%atom3 = -1
       nullify(top%next)
       yyparse_angleA_top(i) = top
    enddo
  end subroutine fmod_alloc_yyparse_angleA_top
!----------------------------------------------------------------------
!>
!! \brief  Subroutines to read inputs relating to angle potential.
!! \author Kensuke Iwahashi
!<
  subroutine fmod_alloc_yyparse_angleA_n(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(yyparse_angleA_n(ivalue))
    yyparse_angleA_n = 0
  end subroutine fmod_alloc_yyparse_angleA_n
!----------------------------------------------------------------------
!>
!! \brief  Subroutines to read inputs relating to angle potential.
!! \author Kensuke Iwahashi
!<
  subroutine fmod_alloc_yyparse_angleB_top(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    integer(4) :: i
    type(yyparse_angleB),pointer :: top
    allocate(yyparse_angleB_top(ivalue))
    do i=1,ivalue
       allocate(top)
       top%atom1 = -1
       top%atom3 = -1
       nullify(top%next)
       yyparse_angleB_top(i) = top
    enddo
  end subroutine fmod_alloc_yyparse_angleB_top
!----------------------------------------------------------------------
!>
!! \brief  Subroutines to read inputs relating to angle potential.
!! \author Kensuke Iwahashi
!<
  subroutine fmod_alloc_yyparse_angleB_n(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(yyparse_angleB_n(ivalue))
    yyparse_angleB_n = 0
  end subroutine fmod_alloc_yyparse_angleB_n
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to keep angle pair in liner list temporally.
!! \author Kensuke Iwahashi
!<
  subroutine add_to_angleA_list(atom1,atom2,atom3,Ktheta,theta0)
    use trajectory_org
    use param
    implicit none
    integer(4), intent(in) :: atom1,atom2,atom3
    real(8), intent(in) :: Ktheta,theta0
    type(yyparse_angleA),pointer :: new1,p,next
    if (atom1 > npara) then
       write(0,'(a,i0,a)') &
            &    'The number of angle is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
    endif
    if (atom2 > npara) then
       write(0,'(a,i0,a)') &
            &    'The number of angle is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
    endif
    if (atom3 > npara) then
       write(0,'(a,i0,a)') &
            &    'The number of angle is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
    endif
    allocate(new1)
    nullify(new1%next)
    new1%atom2 = atom2
    new1%atom3 = atom3
    new1%Ktheta = Ktheta
    new1%theta0 = theta0
    p => yyparse_angleA_top(atom1)
    !      insert to list, sorting from smaller
    do while (.true.)
       next => p%next
       if (.not. associated(next)) then
          p%next => new1
          yyparse_angleA_n(atom1) = yyparse_angleA_n(atom1) + 1
          yyparse_angleA_total = yyparse_angleA_total + 1
          exit
       endif
       !        doubled must be merged
       if (atom2 == next%atom2 .and. atom3 == next%atom3) then
          next%Ktheta = new1%Ktheta
          next%theta0 = new1%theta0
          deallocate(new1)
          exit
       else if (atom2 <= next%atom2 .and. atom3 < next%atom3) then
          new1%next => p%next
          p%next => new1
          yyparse_angleA_n(atom1) = yyparse_angleA_n(atom1) + 1
          yyparse_angleA_total = yyparse_angleA_total + 1
          exit
       endif
       p => next
    enddo
  end subroutine add_to_angleA_list
! -----------------------------------------------------------------------
!>
!! \brief  Subroutine to keep angle pair in liner list temporally.
!! \author Kensuke Iwahashi
!<
  subroutine add_to_angleB_list(atom2,atom1,atom3,Ktheta,theta0)
    use trajectory_org
    use param
    implicit none
    integer(4), intent(in) :: atom1,atom2,atom3
    real(8), intent(in) :: Ktheta,theta0
    type(yyparse_angleB),pointer :: new1,p,next
    if (atom1 > npara) then
       write(0,'(a,i0,a)') &
            &    'The number of angle is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
    endif
    if (atom2 > npara) then
       write(0,'(a,i0,a)') &
            &    'The number of angle is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
    endif
    if (atom3 > npara) then
       write(0,'(a,i0,a)') &
            &    'The number of angle is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
    endif
    allocate(new1)
    nullify(new1%next)
    new1%atom1 = atom1
    new1%atom3 = atom3
    new1%Ktheta = Ktheta
    new1%theta0 = theta0
    p => yyparse_angleB_top(atom2)
    !      insert to list, sorting from smaller
    do while (.true.)
       next => p%next
       if (.not. associated(next)) then
          p%next => new1
          yyparse_angleB_n(atom2) = yyparse_angleB_n(atom2) + 1
          yyparse_angleB_total = yyparse_angleB_total + 1
          exit
       endif
       !        doubled must be merged
       if (atom1 == next%atom1 .and. atom3 == next%atom3) then
          next%Ktheta = new1%Ktheta
          next%theta0 = new1%theta0
          deallocate(new1)
          exit
       else if (atom1 <= next%atom1 .and. atom3 < next%atom3) then
          new1%next => p%next
          p%next => new1
          yyparse_angleB_n(atom2) = yyparse_angleB_n(atom2) + 1
          yyparse_angleB_total = yyparse_angleB_total + 1
          exit
       endif
       p => next
    enddo
  end subroutine add_to_angleB_list

end module parse_angle
!----------------------------------------------------------------------
!>
!! \brief  Module to calculate angle potential.
!! \author Kensuke Iwahashi
!<
!----------------------------------------------------------------------
module angle
  use omp_lib
  use parse_angle
  use subcell
  implicit none
  integer(4) :: nangles=0
  integer(4) :: nangleAs=0, nangleBs=0
  integer(4), allocatable :: atom1(:), atom2(:), atom3(:)
  real(8), allocatable :: Ktheta(:), theta0(:)
  integer(4), allocatable :: atom2A(:,:), atom3A(:,:)
  real(8), allocatable :: KthetaA(:,:), theta0A(:,:)
  integer(4), allocatable :: atom1B(:,:), atom3B(:,:)
  real(8), allocatable :: KthetaB(:,:), theta0B(:,:)
  integer(4), allocatable :: topA(:), nA(:)
  integer(4), allocatable :: topB(:), nB(:)

contains

!  subroutine init_angles()
!  end subroutine init_angles

!  subroutine intra_angles()
!  end subroutine intra_angles

!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to check whether j-atoms are communicated correctly
!! \author Yoshimichi Andoh
!<
  subroutine check_atombound_angle()
!-----------------------------------------------------------------------
    use trajectory_mpi
    use param
    use md_condition
    use mpi_tool
    implicit none
    integer(4) :: k0,i0,i0a,i0b,i0c,j
    integer(4) :: iA,iB,ipar
    integer(4) :: icheck

    icheck=0

!$omp parallel default(none) &
!$omp& private(k0,i0,i0a,i0b,i0c,iA,iB,j,ipar) &
!$omp& shared(nselfseg,lsegtop,lseg_natoms,m2i,i2m) &
!$omp& shared(paranum,nA,nB,atom2A,atom3A,atom1B,atom3B) &
!$omp& reduction(+:icheck)
!$omp do
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          iA  = m2i(i0)
          ipar = paranum(iA)
          do j = 1, nA(ipar)
             i0b = i2m(atom2A(ipar,j)+(iA-ipar))
             i0c = i2m(atom3A(ipar,j)+(iA-ipar))
             if(i0b.eq.-1) icheck = icheck -1
             if(i0c.eq.-1) icheck = icheck -1
          enddo
          do j = 1, nB(ipar)
             i0a = i2m(atom1B(ipar,j)+(iA-ipar))
             i0c = i2m(atom3B(ipar,j)+(iA-ipar))
             if(i0a.eq.-1) icheck = icheck -1
             if(i0c.eq.-1) icheck = icheck -1
          enddo
       enddo
    enddo
!$omp end do
!$omp end parallel

    if(icheck.le.-1) then
       write(0,*) 'ERROR[add_charmm_angle_a] &
            &   There is a particle outside the area.'
       write(0,*) 'Myrank,mdstep=',myrank,mdstep
       call modylas_abort()
    end if

  end subroutine check_atombound_angle
!-----------------------------------------------------------------------
  subroutine add_to_mol_angle_list(ivalue1,atom1,atom2,atom3,Ktheta,theta0)
    use parse_angle
    use mol_info, only : para_start
    implicit none
    integer(4), intent(in) :: ivalue1
    integer(4), intent(inout) :: atom1, atom2, atom3
    real(8), intent(in) :: Ktheta, theta0

    atom1 = atom1 + para_start(ivalue1)
    atom2 = atom2 + para_start(ivalue1)
    atom3 = atom3 + para_start(ivalue1)
    call add_to_angleA_list(atom1, atom2, atom3, Ktheta, theta0)
    call add_to_angleA_list(atom3, atom2, atom1, Ktheta, theta0)
    call add_to_angleB_list(atom2, atom1, atom3, Ktheta, theta0)
  end subroutine add_to_mol_angle_list
!-----------------------------------------------------------------------
  subroutine set_md_charmm_angleA
    use parse_angle
    use param
    implicit none
    integer(4) :: i,j,k,maxi(1)
    type(yyparse_angleA),pointer :: p,freed

    nangleAs = yyparse_angleA_total
    if (yyparse_angleA_total == 0)  return
    maxi = maxloc(yyparse_angleA_n)
    maxi = yyparse_angleA_n(maxi(1))
    allocate(atom2A(npara, maxi(1)))
    allocate(atom3A(npara, maxi(1)))
    allocate(KthetaA(npara, maxi(1)))
    allocate(theta0A(npara, maxi(1)))
    allocate(nA(npara))
    atom2A = 0
    atom3A = 0
    KthetaA = 0.0d0
    theta0A = 0.0d0
    nA = 0
    !     copying
    do i=1, npara
       nA(i) = yyparse_angleA_n(i)
       p => yyparse_angleA_top(i)
       do j=1,yyparse_angleA_n(i)
          p => p%next
          atom2A(i,j) = p%atom2
          atom3A(i,j) = p%atom3
          KthetaA(i,j) = p%Ktheta
          theta0A(i,j) = p%theta0
       enddo
    enddo

    deallocate(yyparse_angleA_top, yyparse_angleA_n)
  end subroutine set_md_charmm_angleA
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to read angle parameter B.
!! \author Kazushi FUJIMOTO
!<
  subroutine set_md_charmm_angleB
    use parse_angle
    use param
    implicit none
    integer(4) :: i,j,k,maxi(1)
    type(yyparse_angleB),pointer :: p,freed
    nangleBs =  yyparse_angleB_total
    if (yyparse_angleB_total == 0)  return
    maxi = maxloc(yyparse_angleB_n)
    maxi = yyparse_angleB_n(maxi(1))
    allocate(atom1B(npara, maxi(1)))
    allocate(atom3B(npara, maxi(1)))
    allocate(KthetaB(npara, maxi(1)))
    allocate(theta0B(npara, maxi(1)))
    allocate(nB(npara))
    atom1B = -11111
    atom3B = 0
    KthetaB = 0.0d0
    theta0B = 0.0d0
    nB = 0
    !     copying
    do i=1, npara
       nB(i) = yyparse_angleB_n(i)
       p => yyparse_angleB_top(i)
       do j=1,yyparse_angleB_n(i)
          p => p%next
          atom1B(i,j) = p%atom1
          atom3B(i,j) = p%atom3
          KthetaB(i,j) = p%Ktheta
          theta0B(i,j) = p%theta0
       enddo
    enddo

    deallocate(yyparse_angleB_top, yyparse_angleB_n)
  end subroutine set_md_charmm_angleB

!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to read angle parameter A.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_md_charmm_angleA
    use trajectory_org
    use device_numbers, only : f_mdff
    use param
    use mpi_tool
    implicit none
    include 'mpif.h'
    integer(4) :: ierr
    integer(4) :: num

    if(myrank==0) then
       read(f_mdff) nangleAs
    endif
    call MPI_Bcast(nangleAs,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    if(nangleAs.eq.0) return

    if(myrank==0) read(f_mdff) num
    call MPI_Bcast(num,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)

    allocate(nA(npara))
    allocate(atom2A(npara, num))
    allocate(atom3A(npara, num))
    allocate(KthetaA(npara, num))
    allocate(theta0A(npara, num))

    if(myrank==0) then
       !      read(f_mdff) topA
       read(f_mdff) nA
       read(f_mdff) atom2A
       read(f_mdff) atom3A
       read(f_mdff) KthetaA
       read(f_mdff) theta0A
    endif
    call MPI_Bcast(nA,npara,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom2A,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom3A,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(KthetaA,npara*num,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(theta0A,npara*num,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  end subroutine read_md_charmm_angleA
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write angle parameter A in mdff.bin File.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_md_charmm_angleA
    use device_numbers, only : f_mdff
    implicit none
    integer(4) :: num
    write(f_mdff) nangleAs
    if(nangleAs.eq.0) return
    num = size(atom2A(:,:), 2)
    write(f_mdff) num
    write(f_mdff) nA
    write(f_mdff) atom2A
    write(f_mdff) atom3A
    write(f_mdff) KthetaA
    write(f_mdff) theta0A
  end subroutine write_md_charmm_angleA
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write angle parameter A.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_md_charmm_angleA
    use device_numbers, only : f_mdff
    implicit none
    integer(4) :: num
    write(*,*) '[write_md_charmm_angleA]'
    write(*,*) nangleAs
    if(nangleAs.eq.0) return
    num = size(atom2A(:,:), 2)
    write(*,*) num
    write(*,*) nA
    write(*,*) atom2A
    write(*,*) atom3A
    write(*,*) KthetaA
    write(*,*) theta0A
  end subroutine write_memory_md_charmm_angleA
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to read angle parameter B.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_md_charmm_angleB
    use trajectory_org
    use device_numbers, only : f_mdff
    use param
    use mpi_tool
    implicit none
    include 'mpif.h'
    integer(4) :: ierr
    integer(4) :: num

    if(myrank==0) then
       read(f_mdff) nangleBs
    endif
    call MPI_Bcast(nangleBs,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    if(nangleBs.eq.0) return

    if(myrank==0) read(f_mdff) num
    call MPI_Bcast(num,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)

    allocate(nB(npara))
    allocate(atom1B(npara, num))
    allocate(atom3B(npara, num))
    allocate(KthetaB(npara, num))
    allocate(theta0B(npara, num))
    if(myrank==0) then
       !      read(f_mdff) topB
       read(f_mdff) nB
       read(f_mdff) atom1B
       read(f_mdff) atom3B
       read(f_mdff) KthetaB
       read(f_mdff) theta0B
    endif
    call MPI_Bcast(nB,npara,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom1B,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom3B,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(KthetaB,npara*num,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(theta0B,npara*num,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  end subroutine read_md_charmm_angleB
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write angle parameter B in mdff.bin File.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_md_charmm_angleB
    use device_numbers, only : f_mdff
    implicit none
    integer(4) :: num
    write(f_mdff) nangleBs
    if(nangleBs.eq.0) return
    num = size(atom1B(:,:), 2)
    write(f_mdff) num
    write(f_mdff) nB
    write(f_mdff) atom1B
    write(f_mdff) atom3B
    write(f_mdff) KthetaB
    write(f_mdff) theta0B
  end subroutine write_md_charmm_angleB
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write angle parameter B.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_md_charmm_angleB
    implicit none
    integer(4) :: num
    write(*,*) '[write_md_charmm_angleB]'
    write(*,*) nangleBs
    if(nangleBs.eq.0) return
    num = size(atom1B(:,:), 2)
    write(*,*) num
    write(*,*) nB
    write(*,*) atom1B
    write(*,*) atom3B
    write(*,*) KthetaB
    write(*,*) theta0B
  end subroutine write_memory_md_charmm_angleB
!-------------------------------------------------------------------------
  subroutine read_md_charmm_angle
    !fujimoto
    use device_numbers, only : f_mdff
    implicit none

    read(f_mdff) nangles
    if(nangles.eq.0) return
    allocate(atom1(nangles))
    allocate(atom2(nangles))
    allocate(atom3(nangles))
    allocate(Ktheta(nangles))
    allocate(theta0(nangles))
    read(f_mdff) atom1
    read(f_mdff) atom2
    read(f_mdff) atom3
    read(f_mdff) Ktheta
    read(f_mdff) theta0
  end subroutine read_md_charmm_angle
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate energy and force of angle.
!! \author Kensuke Iwahashi
!<
  subroutine add_charmm_angle_a()
!-----------------------------------------------------------------------
!
!     V = Ktheta(theta - theta0)^2
!
!          B
!         / \    theta = angle(A--B--C)
!        /   \      !! theta belongs to [0, pi]
!       A     C
!
!   /* theta => t */
    use trajectory_org
    use trajectory_mpi
    use atom_virial
    use forces
    use md_monitors
    use md_const
    use param
    use mpi_tool
#ifdef ONEPROC_AXIS
    use boundary, only : pbc_pair
#endif
    implicit none
    integer(4) :: j, i0a, i0b, i0c, ip,iAA, i0,ipar, k0
    real(8) :: t_t0, cost, ab, i_ab, cb, i_cb, fbuf
    real(8) :: dx1, dy1, dz1, dx2, dy2, dz2
    real(8) :: Fa__x, Fa__y, Fa__z
    real(8) :: Fb__x, Fb__y, Fb__z
    real(8) :: Fc__x, Fc__y, Fc__z
    real(8) :: vbuf1__x, vbuf1__y, vbuf1__z, vbuf2__x, vbuf2__y, vbuf2__z
    real(8) :: v11,v22,v33
    real(8) :: wk_v11, wk_v22, wk_v33
    real(8) :: wkvd(6), vd(6) ! debug
    real(8) :: v21,v31,v32
    real(8) :: wk_v21, wk_v31, wk_v32
    real(8) :: fi0(3)
    real(8) :: Uangle, oUa
    integer(4) :: iam
    include 'mpif.h'

    call check_atombound_angle()

    iam = 0
    Uangle = 0.0d0
    wk_v11=0d0;wk_v22=0d0;wk_v33=0d0
    wk_v21=0d0;wk_v31=0d0;wk_v32=0d0

!$omp parallel default(none) &
!$omp& private(v11,v22,v33) &
!$omp& reduction(+:wk_v11,wk_v22,wk_v33) &
!$omp& private(v21,v31,v32) &
!$omp& reduction(+:wk_v21,wk_v31,wk_v32) &
!$omp& private(k0,i0,ipar,iam) &
!$omp& private(j,ip,i0a,i0b,i0c,dx1,dy1,dz1,dx2,dy2,dz2,iAA) &
!$omp& private(ab,i_ab,cb,i_cb,cost,t_t0,fbuf) &
!$omp& private(vbuf1__x,vbuf1__y,vbuf1__z) &
!$omp& private(vbuf2__x,vbuf2__y,vbuf2__z) &
!$omp& private(Fa__x,Fa__y,Fa__z) &
!$omp& private(Fb__x,Fb__y,Fb__z) &
!$omp& private(Fc__x,Fc__y,Fc__z) &
!$omp& private(fi0) &
!$omp& shared(m2i,i2m,nA,atom2A,atom3A,topA,wkxyz) &
!$omp& shared(theta0A,KthetaA,w3_f,wk_vir2) &
!$omp& shared(nB,topB,atom1B,atom3B,theta0B,KthetaB,paranum) &
!$omp& shared(tag,na_per_cell) &
!$omp& shared(lsegtop,lseg_natoms,nselfseg) &
!$omp& reduction(+:Uangle)
!$  iam = omp_get_thread_num()
!$omp do
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          iAA  = m2i(i0)
          ipar = paranum(iAA)
          fi0(1:3)=0d0
!
          do j = 1, nB(ipar)
             i0a = i2m(atom1B(ipar,j)+(iAA-ipar))
             i0c = i2m(atom3B(ipar,j)+(iAA-ipar))
             dx1 = wkxyz(1,i0a) - wkxyz(1,i0)
             dy1 = wkxyz(2,i0a) - wkxyz(2,i0)
             dz1 = wkxyz(3,i0a) - wkxyz(3,i0)
             dx2 = wkxyz(1,i0c) - wkxyz(1,i0)
             dy2 = wkxyz(2,i0c) - wkxyz(2,i0)
             dz2 = wkxyz(3,i0c) - wkxyz(3,i0)
#ifdef ONEPROC_AXIS
             call pbc_pair(dx1,dy1,dz1)
             call pbc_pair(dx2,dy2,dz2)
#endif
             ab = sqrt(dx1*dx1+dy1*dy1+dz1*dz1)
             i_ab = 1.0d0 / ab
             cb = sqrt(dx2*dx2+dy2*dy2+dz2*dz2)
             i_cb = 1.0d0 / cb
             cost=(dx1*dx2+dy1*dy2+dz1*dz2)*i_ab*i_cb
             t_t0 = acos(cost) - theta0B(ipar,j)
             fbuf = -2.0d0 * KthetaB(ipar,j) * t_t0 / sqrt(1.0d0 - cost * cost)
             vbuf1__x = dx1 * cost * i_ab * i_ab
             vbuf1__y = dy1 * cost * i_ab * i_ab
             vbuf1__z = dz1 * cost * i_ab * i_ab
             vbuf2__x = dx2 * i_ab * i_cb
             vbuf2__y = dy2 * i_ab * i_cb
             vbuf2__z = dz2 * i_ab * i_cb
             Fa__x = (vbuf1__x - vbuf2__x) * fbuf
             Fa__y = (vbuf1__y - vbuf2__y) * fbuf
             Fa__z = (vbuf1__z - vbuf2__z) * fbuf
             vbuf1__x = dx2 * cost * i_cb * i_cb
             vbuf1__y = dy2 * cost * i_cb * i_cb
             vbuf1__z = dz2 * cost * i_cb * i_cb
             vbuf2__x = dx1 * i_ab * i_cb
             vbuf2__y = dy1 * i_ab * i_cb
             vbuf2__z = dz1 * i_ab * i_cb
             Fc__x = (vbuf1__x - vbuf2__x) * fbuf
             Fc__y = (vbuf1__y - vbuf2__y) * fbuf
             Fc__z = (vbuf1__z - vbuf2__z) * fbuf
             Fb__x = -Fa__x - Fc__x
             Fb__y = -Fa__y - Fc__y
             Fb__z = -Fa__z - Fc__z
             fi0(1)=fi0(1)+Fb__x
             fi0(2)=fi0(2)+Fb__y
             fi0(3)=fi0(3)+Fb__z
!default: MTD
             w3_f(1,i0a,iam) = w3_f(1,i0a,iam) + Fa__x
             w3_f(2,i0a,iam) = w3_f(2,i0a,iam) + Fa__y
             w3_f(3,i0a,iam) = w3_f(3,i0a,iam) + Fa__z
             w3_f(1,i0c,iam) = w3_f(1,i0c,iam) + Fc__x
             w3_f(2,i0c,iam) = w3_f(2,i0c,iam) + Fc__y
             w3_f(3,i0c,iam) = w3_f(3,i0c,iam) + Fc__z
             Uangle = Uangle + KthetaB(ipar,j) * t_t0 * t_t0
             v11 = dx1 * Fa__x + dx2 * Fc__x
             v22 = dy1 * Fa__y + dy2 * Fc__y
             v33 = dz1 * Fa__z + dz2 * Fc__z
             v21 = dy1 * Fa__x + dy2 * Fc__x
             v31 = dz1 * Fa__x + dz2 * Fc__x
             v32 = dz1 * Fa__y + dz2 * Fc__y
             wk_v11 = wk_v11 + v11
             wk_v22 = wk_v22 + v22
             wk_v33 = wk_v33 + v33
             wk_v21 = wk_v21 + v21
             wk_v31 = wk_v31 + v31
             wk_v32 = wk_v32 + v32
          enddo
          w3_f(1,i0,iam) = w3_f(1,i0,iam) + fi0(1)
          w3_f(2,i0,iam) = w3_f(2,i0,iam) + fi0(2)
          w3_f(3,i0,iam) = w3_f(3,i0,iam) + fi0(3)
       enddo ! i0
    enddo ! k0
!$omp end do
!$omp end parallel

    wk_p_energy = wk_p_energy + Uangle
    wk_vir2(1,0) = wk_vir2(1,0) + wk_v11
    wk_vir2(2,0) = wk_vir2(2,0) + wk_v22
    wk_vir2(3,0) = wk_vir2(3,0) + wk_v33
    wk_vir2(4,0) = wk_vir2(4,0) + wk_v21
    wk_vir2(5,0) = wk_vir2(5,0) + wk_v31
    wk_vir2(6,0) = wk_vir2(6,0) + wk_v32
#ifdef DEBUGFCE
    oUa=0d0
    call mpi_allreduce(Uangle,oUa,1, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,ipar)
#ifdef KCAL
    if(myrank==0)write(*,*)'Pot(angle)=',oUa *kJ_mol/4.184d0,'[kcal/mol]'
#else
    if(myrank==0)write(*,*)'Pot(angle)=',oUa *kJ_mol,'[kJ/mol]'
#endif
#endif

  end subroutine add_charmm_angle_a

end module angle

