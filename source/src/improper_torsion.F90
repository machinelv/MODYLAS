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
!! \brief  Module and subroutines to calculate improper torsion potential.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to parse parameters of improper torsion potential.
!! \author Kensuke Iwahashi
!<
module parse_improper_torsion
  use mpi_tool
  implicit none
!
  type yyparse_itorsionA
     integer(4) :: atom2, atom3, atom4
#if defined(GROEXT_OPLS) || defined(GAFF)
!gro-converter only
     real(8) :: Kd, nd, delta
!gro-conterter only
!note: gaff also uses these parameters: E=K*[1+d*cos(n*phi)],
!      where K=Kd, d=nd, n=delta
#else
     real(8) :: Kit, psi0
#endif
     type(yyparse_itorsionA),pointer :: next
  end type yyparse_itorsionA
  type(yyparse_itorsionA),pointer :: yyparse_itorsionA_top(:)
  integer(4),allocatable :: yyparse_itorsionA_n(:)
  integer(4) :: yyparse_itorsionA_total=0
!
  type yyparse_itorsionB
     integer(4) :: atom1, atom3, atom4
#if defined(GROEXT_OPLS) || defined(GAFF)
!gro-converter only
     real(8) :: Kd, nd, delta
!gro-conterter only
#else
     real(8) :: Kit, psi0
#endif
     type(yyparse_itorsionB),pointer :: next
  end type yyparse_itorsionB
  type(yyparse_itorsionB),pointer :: yyparse_itorsionB_top(:)
  integer(4),allocatable :: yyparse_itorsionB_n(:)
  integer(4) :: yyparse_itorsionB_total=0

contains

  subroutine fmod_alloc_yyparse_itorsA_top(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    integer(4) :: i
    type(yyparse_itorsionA),pointer :: top
    allocate(yyparse_itorsionA_top(ivalue))
    do i=1,ivalue
       allocate(top)
       top%atom2 = -1
       top%atom3 = -1
       top%atom4 = -1
       nullify(top%next)
       yyparse_itorsionA_top(i) = top
    enddo
  end subroutine fmod_alloc_yyparse_itorsA_top
!----------------------------------------------------------------------
  subroutine fmod_alloc_yyparse_itorsA_n(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(yyparse_itorsionA_n(ivalue))
    yyparse_itorsionA_n = 0
  end subroutine fmod_alloc_yyparse_itorsA_n
!----------------------------------------------------------------------
  subroutine fmod_alloc_yyparse_itorsB_top(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    integer(4) :: i
    type(yyparse_itorsionB),pointer :: top
    allocate(yyparse_itorsionB_top(ivalue))
    do i=1,ivalue
       allocate(top)
       top%atom1 = -1
       top%atom3 = -1
       top%atom4 = -1
       nullify(top%next)
       yyparse_itorsionB_top(i) = top
    enddo
  end subroutine fmod_alloc_yyparse_itorsB_top
!----------------------------------------------------------------------
  subroutine fmod_alloc_yyparse_itorsB_n(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(yyparse_itorsionB_n(ivalue))
    yyparse_itorsionB_n = 0
  end subroutine fmod_alloc_yyparse_itorsB_n
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to keep improper pair in liner list temporally.
!! \author Kensuke Iwahashi
!<
  subroutine add_to_itorsionA_list(atom1,atom2,atom3,atom4, &
#if defined(GROEXT_OPLS) || defined(GAFF)
    &                              Kd,nd,delta)
#else
    &                              Kit,psi0)
#endif
    use trajectory_org
    use param
    implicit none
    integer(4), intent(in) :: atom1,atom2,atom3,atom4
#if defined(GROEXT_OPLS) || defined(GAFF)
    real(8), intent(in) :: Kd,nd,delta
#else
    real(8) :: Kit,psi0
#endif
    type(yyparse_itorsionA),pointer :: new1,p,next
    if (atom1 > npara) then
       write(0,'(a,i0,a)') &
            &    'ERROR: The number of itorsion is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom2 > npara) then
       write(0,'(a,i0,a)') &
            &    'ERROR: The number of itorsion is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom3 > npara) then
       write(0,'(a,i0,a)') &
            &    'ERROR: The number of itorsion is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom4 > npara) then
       write(0,'(a,i0,a)') &
            &    'ERROR: The number of itorsion is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    allocate(new1)
    nullify(new1%next)
    new1%atom2 = atom2
    new1%atom3 = atom3
    new1%atom4 = atom4
#if defined(GROEXT_OPLS) || defined(GAFF)
    new1%Kd=Kd 
    new1%nd=nd
    new1%delta=delta
#else
    new1%Kit = Kit
    new1%psi0 = psi0
#endif
    p => yyparse_itorsionA_top(atom1)
    !     insert to list, sorting from smaller
    do while (.true.)
       next => p%next
       if (.not. associated(next)) then
          p%next => new1
          yyparse_itorsionA_n(atom1) = yyparse_itorsionA_n(atom1) + 1
          yyparse_itorsionA_total = yyparse_itorsionA_total + 1
          exit
       endif
       if (atom2 <= next%atom2 .and. atom3 <= next%atom3 .and. atom4 <= next%atom4) then
          new1%next => p%next
          p%next => new1
          yyparse_itorsionA_n(atom1) = yyparse_itorsionA_n(atom1) + 1
          yyparse_itorsionA_total = yyparse_itorsionA_total + 1
          exit
       endif
       p => next
    enddo
  end subroutine add_to_itorsionA_list
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to keep improper pair in liner list temporally.
!! \author Kensuke Iwahashi
!<
  subroutine add_to_itorsionB_list(atom2,atom1,atom3,atom4, &
#if defined(GROEXT_OPLS) || defined(GAFF)
    &                              Kd,nd,delta)
#else
    &                              Kit,psi0)
#endif
    use trajectory_org
    use param
    implicit none
    integer(4), intent(in) :: atom1,atom2,atom3,atom4
#if defined(GROEXT_OPLS) || defined(GAFF)
    real(8), intent(in) :: Kd,nd,delta
#else
    real(8), intent(in) :: Kit,psi0
#endif
    type(yyparse_itorsionB),pointer :: new1,p,next
    if (atom1 > npara) then
       write(0,'(a,i0,a)') &
            &    'ERROR: The number of itorsion is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom2 > npara) then
       write(0,'(a,i0,a)') &
            &    'ERROR: The number of itorsion is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom3 > npara) then
       write(0,'(a,i0,a)') &
            &    'ERROR: The number of itorsion is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom4 > npara) then
       write(0,'(a,i0,a)') &
            &    'ERROR: The number of itorsion is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    allocate(new1)
    nullify(new1%next)
    new1%atom1 = atom1
    new1%atom3 = atom3
    new1%atom4 = atom4
#if defined(GROEXT_OPLS) || defined(GAFF)
    new1%Kd=Kd 
    new1%nd=nd
    new1%delta=delta
#else
    new1%Kit = Kit
    new1%psi0 = psi0
#endif
    p => yyparse_itorsionB_top(atom2)
    !     insert to list, sorting from smaller
    do while (.true.)
       next => p%next
       if (.not. associated(next)) then
          p%next => new1
          yyparse_itorsionB_n(atom2) = yyparse_itorsionB_n(atom2) + 1
          yyparse_itorsionB_total = yyparse_itorsionB_total + 1
          exit
       endif
       if (atom1 <= next%atom1 .and. atom3 <= next%atom3 .and. atom4 < next%atom4) then
          new1%next => p%next
          p%next => new1
          yyparse_itorsionB_n(atom2) = yyparse_itorsionB_n(atom2) + 1
          yyparse_itorsionB_total = yyparse_itorsionB_total + 1
          exit
       endif
       p => next
    enddo
  end subroutine add_to_itorsionB_list

!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to keep improper pair in liner list temporally.
!!         used with gromacs-converter.
!! \author Yoshimichi ANDOH
!<
! subroutine add_to_itorsionA_list_gro(atom1,atom2,atom3,atom4, &
!   &                                  Kd,nd,delta)
! end subroutine add_to_itorsionA_list_gro

!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to keep improper pair in liner list temporally.
!!         used with gromacs-converter.
!! \author Yoshimichi ANDOH
!<
! subroutine add_to_itorsionB_list_gro(atom1,atom2,atom3,atom4, &
!   &                                  Kd,nd,delta)
! end subroutine add_to_itorsionB_list_gro

!-----------------------------------------------------------------------
  subroutine add_to_mol_itorsion_list &
#if defined(GROEXT_OPLS) || defined(GAFF)
    &              (ivalue1, atom1, atom2, atom3, atom4, Kd, nd, delta)
#else
    &              (ivalue1, atom1, atom2, atom3, atom4, Kit, psi0)
#endif
    use mol_info, only : para_start
    implicit none
    integer(4), intent(in) :: ivalue1
    integer(4), intent(inout) :: atom1, atom2, atom3, atom4
#if defined(GROEXT_OPLS) || defined(GAFF)
    real(8), intent(in) :: Kd, nd, delta
#else
    real(8), intent(in) :: Kit, psi0
#endif

    atom1 = atom1 + para_start(ivalue1)
    atom2 = atom2 + para_start(ivalue1)
    atom3 = atom3 + para_start(ivalue1)
    atom4 = atom4 + para_start(ivalue1)

#if defined(GROEXT_OPLS) || defined(GAFF)
    call add_to_itorsionA_list(atom1, atom2, atom3, atom4, Kd, nd, delta)
    call add_to_itorsionA_list(atom4, atom3, atom2, atom1, Kd, nd, delta)
    call add_to_itorsionB_list(atom2, atom1, atom3, atom4, Kd, nd, delta)
    call add_to_itorsionB_list(atom3, atom4, atom2, atom1, Kd, nd, delta)
#else
    call add_to_itorsionA_list(atom1, atom2, atom3, atom4, Kit, psi0)
    call add_to_itorsionA_list(atom4, atom3, atom2, atom1, Kit, psi0)
    call add_to_itorsionB_list(atom2, atom1, atom3, atom4, Kit, psi0)
    call add_to_itorsionB_list(atom3, atom4, atom2, atom1, Kit, psi0)
#endif

  end subroutine add_to_mol_itorsion_list

!-----------------------------------------------------------------------
! subroutine add_to_mol_itorsion_list_gro &
!   &              (ivalue1, atom1, atom2, atom3, atom4,   &
!   &                                  Kd,nd,delta)
! end subroutine add_to_mol_itorsion_list_gro

end module parse_improper_torsion
!----------------------------------------------------------------------
!>
!! \brief  Module to calculate improper torsion potential.
!! \author Kensuke Iwahashi
!<
!----------------------------------------------------------------------
module improper_torsion
  use omp_lib
  use parse_improper_torsion
  use subcell
  implicit none
  integer(4) :: nitorsions=0
  integer(4) :: nitorsionAs=0, nitorsionBs=0
  integer(4), allocatable :: atom1(:), atom2(:), atom3(:), atom4(:)
  integer(4), allocatable :: atom2A(:,:), atom3A(:,:), atom4A(:,:)
  integer(4), allocatable :: atom1B(:,:), atom3B(:,:), atom4B(:,:)
  integer(4), allocatable :: topA(:), nA(:)
  integer(4), allocatable :: topB(:), nB(:)
#if defined(GROEXT_OPLS) || defined(GAFF)
!gro-converter only
  real(8), allocatable :: Kd(:), nd(:), delta(:)
  real(8), allocatable :: KdA(:,:), ndA(:,:), deltaA(:,:)
  real(8), allocatable :: KdB(:,:), ndB(:,:), deltaB(:,:)
!gro-converter only
#else
  real(8), allocatable :: Kit(:), psi0(:)
  real(8), allocatable :: KitA(:,:), psi0A(:,:)
  real(8), allocatable :: KitB(:,:), psi0B(:,:)
#endif

contains

  !subroutine init_impropers
  !end subroutine init_impropers

  !subroutine intra_impropers
  !end subroutine intra_impropers

!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to check whether j-,k-,l-atoms are communicated correctly
!! \author Yoshimichi Andoh
!<
  subroutine check_atombound_itorsion()
!-----------------------------------------------------------------------
    use trajectory_mpi
    use param
    use md_condition
    use mpi_tool
    implicit none
    integer(4) :: k0,i0,j,i0a,i0b,i0c,i0d
    integer(4) :: iA,ipar
    integer(4) :: icheck

    icheck=0

    if(allow_bond_breaking<=0.0d0) then

!$omp parallel default(none) &
!$omp& private(k0,i0,i0a,i0b,i0c,i0d,j,iA,ipar) &
!$omp& shared(nselfseg,lsegtop,lseg_natoms,m2i,i2m) &
!$omp& shared(paranum,nA,nB) &
!$omp& shared(atom1B,atom3B,atom4B,atom2A,atom3A,atom4A) &
!$omp& reduction(+:icheck)
!$omp do
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          iA   = m2i(i0)
          ipar = paranum(iA)
          do j = 1, nB(ipar)
             i0a = i2m(atom1B(ipar,j)+(iA-ipar))
             i0c = i2m(atom3B(ipar,j)+(iA-ipar))
             i0d = i2m(atom4B(ipar,j)+(iA-ipar))
             if(i0a.eq.-1) icheck = icheck -1
             if(i0c.eq.-1) icheck = icheck -1
             if(i0d.eq.-1) icheck = icheck -1
          enddo
       enddo
    enddo
!$omp end do
!$omp end parallel

  else

!$omp parallel default(none) &
!$omp& private(k0,i0,i0a,i0b,i0c,i0d,j,iA,ipar) &
!$omp& shared(nselfseg,lsegtop,lseg_natoms,m2i,i2m) &
!$omp& shared(paranum,nA,nB) &
#if defined(GROEXT_OPLS) || defined(GAFF)
!$omp& shared(KdA, KdB) &
#else
!$omp& shared(KitA, KitB) &
#endif
!$omp& shared(atom1B,atom3B,atom4B,atom2A,atom3A,atom4A) &
!$omp& reduction(+:icheck)
!$omp do
    do k0=1,nselfseg
      do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
         iA   = m2i(i0)
         ipar = paranum(iA)
         do j = 1, nB(ipar)
            i0a = i2m(atom1B(ipar,j)+(iA-ipar))
            i0c = i2m(atom3B(ipar,j)+(iA-ipar))
            i0d = i2m(atom4B(ipar,j)+(iA-ipar))
#if defined(GROEXT_OPLS) || defined(GAFF)
            if(KdB(ipar,j)>=0.0d0 .and. i0a.eq.-1) icheck = icheck -1
            if(KdB(ipar,j)>=0.0d0 .and. i0c.eq.-1) icheck = icheck -1
            if(KdB(ipar,j)>=0.0d0 .and. i0d.eq.-1) icheck = icheck -1
#else
            if(KitB(ipar,j)>=0.0d0 .and. i0a.eq.-1) icheck = icheck -1
            if(KitB(ipar,j)>=0.0d0 .and. i0c.eq.-1) icheck = icheck -1
            if(KitB(ipar,j)>=0.0d0 .and. i0d.eq.-1) icheck = icheck -1
#endif
         enddo
      enddo
   enddo
!$omp end do
!$omp end parallel

  endif

    if(icheck.le.-1) then
       write(0,*) 'ERROR[add_charmm_itorsion_a] &
            &   There is a particle outside the area.'
       write(0,*) 'Myrank,mdstep=',myrank,mdstep
       call modylas_abort()
    end if
  end subroutine check_atombound_itorsion
!-----------------------------------------------------------------------
  subroutine set_md_charmm_itorsionA
    use param
    implicit none
    integer(4) :: i,j,k,maxi(1)
    type(yyparse_itorsionA),pointer :: p,freed
    nitorsionAs= yyparse_itorsionA_total
    if (yyparse_itorsionA_total == 0)  return
    maxi = maxloc(yyparse_itorsionA_n)
    maxi = yyparse_itorsionA_n(maxi(1))
    allocate(atom2A(npara,maxi(1)))
    allocate(atom3A(npara,maxi(1)))
    allocate(atom4A(npara,maxi(1)))
#if defined(GROEXT_OPLS) || defined(GAFF)
    allocate(KdA(npara,maxi(1)))
    allocate(ndA(npara,maxi(1)))
    allocate(deltaA(npara,maxi(1)))
#else
    allocate(KitA(npara,maxi(1)))
    allocate(psi0A(npara,maxi(1)))
#endif
    allocate(nA(npara))
    !     copying
    do i=1, npara
       nA(i) = yyparse_itorsionA_n(i)
       p => yyparse_itorsionA_top(i)
       do j=1,yyparse_itorsionA_n(i)
          p => p%next
          atom2A(i,j) = p%atom2
          atom3A(i,j) = p%atom3
          atom4A(i,j) = p%atom4
#if defined(GROEXT_OPLS) || defined(GAFF)
          KdA(i,j) = p%Kd
          ndA(i,j) = p%nd
          deltaA(i,j) = p%delta
#else
          KitA(i,j)  = p%Kit
          psi0A(i,j) = p%psi0
#endif
       enddo
    enddo

    deallocate(yyparse_itorsionA_top, yyparse_itorsionA_n)
  end subroutine set_md_charmm_itorsionA
!-----------------------------------------------------------------------
! subroutine set_md_charmm_itorsionA_gro
! end subroutine set_md_charmm_itorsionA_gro
!-----------------------------------------------------------------------
  subroutine set_md_charmm_itorsionB
    use param
    implicit none
    integer(4) :: i,j,k,maxi(1)
    type(yyparse_itorsionB),pointer :: p,freed
    nitorsionBs= yyparse_itorsionB_total
    if (yyparse_itorsionB_total == 0)  return
    maxi = maxloc(yyparse_itorsionB_n)
    maxi = yyparse_itorsionB_n(maxi(1))
    allocate(atom1B(npara,maxi(1)))
    allocate(atom3B(npara,maxi(1)))
    allocate(atom4B(npara,maxi(1)))
#if defined(GROEXT_OPLS) || defined(GAFF)
    allocate(KdB(npara,maxi(1)))
    allocate(ndB(npara,maxi(1)))
    allocate(deltaB(npara,maxi(1)))
#else
    allocate(KitB(npara,maxi(1)))
    allocate(psi0B(npara,maxi(1)))
#endif
    allocate(nB(npara))
    !     copying
    do i=1, npara
       nB(i) = yyparse_itorsionB_n(i)
       p => yyparse_itorsionB_top(i)
       do j=1,yyparse_itorsionB_n(i)
          p => p%next
          atom1B(i,j) = p%atom1
          atom3B(i,j) = p%atom3
          atom4B(i,j) = p%atom4
#if defined(GROEXT_OPLS) || defined(GAFF)
          KdB(i,j) = p%Kd
          ndB(i,j) = p%nd
          deltaB(i,j) = p%delta
#else
          KitB(i,j)  = p%Kit
          psi0B(i,j) = p%psi0
#endif
       enddo
    enddo

    deallocate(yyparse_itorsionB_top, yyparse_itorsionB_n)
  end subroutine set_md_charmm_itorsionB
!-----------------------------------------------------------------------
! subroutine set_md_charmm_itorsionB_gro
! end subroutine set_md_charmm_itorsionB_gro
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to read itorsion parameter A.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_md_charmm_itorsionA
    use trajectory_org
    use device_numbers, only : f_mdff
    use param
    use mpi_tool
    implicit none
    include 'mpif.h'
    integer(4) :: ierr
    integer(4) num

    if(myrank==0) then
       read(f_mdff) nitorsionAs
    endif
    call MPI_Bcast(nitorsionAs,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    if(nitorsionAs.eq.0) return

    if(myrank==0) read(f_mdff) num
    call MPI_Bcast(num,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)

    allocate(nA(npara))
    allocate(atom2A(npara,num))
    allocate(atom3A(npara,num))
    allocate(atom4A(npara,num))
#if defined(GROEXT_OPLS) || defined(GAFF)
    allocate(KdA(npara,num))
    allocate(ndA(npara,num))
    allocate(deltaA(npara,num))
#else
    allocate(KitA(npara,num))
    allocate(psi0A(npara,num))
#endif

    if(myrank==0) then
       !      read(f_mdff) topA
       read(f_mdff) nA
       read(f_mdff) atom2A
       read(f_mdff) atom3A
       read(f_mdff) atom4A
#if defined(GROEXT_OPLS) || defined(GAFF)
       read(f_mdff) KdA
       read(f_mdff) ndA
       read(f_mdff) deltaA
#else
       read(f_mdff) KitA
       read(f_mdff) psi0A
#endif
    endif
    call MPI_Bcast(nA,npara,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom2A,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom3A,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom4A,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
#if defined(GROEXT_OPLS) || defined(GAFF)
    call MPI_Bcast(KdA,npara*num,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(ndA,npara*num,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(deltaA,npara*num,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#else
    call MPI_Bcast(KitA,npara*num,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(psi0A,npara*num,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#endif
  end subroutine read_md_charmm_itorsionA
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write itorsion parameter A in mdff.bin File.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_md_charmm_itorsionA
    use device_numbers, only : f_mdff
    implicit none
    integer(4) num
    write(f_mdff) nitorsionAs
    if(nitorsionAs.eq.0) return
    num = size(atom2A(:,:), 2)
    write(f_mdff) num
    write(f_mdff) nA
    write(f_mdff) atom2A
    write(f_mdff) atom3A
    write(f_mdff) atom4A
#if defined(GROEXT_OPLS) || defined(GAFF)
    write(f_mdff) KdA
    write(f_mdff) ndA
    write(f_mdff) deltaA
#else
    write(f_mdff) KitA
    write(f_mdff) psi0A
#endif
  end subroutine write_md_charmm_itorsionA
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write itorsion parameter A.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_md_charmm_itorsionA
    implicit none
    integer(4) num
    write(*,*) '[write_md_charmm_itorsionA]'
    write(*,*) nitorsionAs
    if(nitorsionAs.eq.0) return
    num = size(atom2A(:,:), 2)
    write(*,*) num
    write(*,*) nA
    write(*,*) atom2A
    write(*,*) atom3A
    write(*,*) atom4A
#if defined(GROEXT_OPLS) || defined(GAFF)
    write(*,*) KdA
    write(*,*) ndA
    write(*,*) deltaA
#else
    write(*,*) KitA
    write(*,*) psi0A
#endif
  end subroutine write_memory_md_charmm_itorsionA
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write itorsion parameter B.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_md_charmm_itorsionB
    use trajectory_org
    use device_numbers, only : f_mdff
    use param
    use mpi_tool
    implicit none
    include 'mpif.h'
    integer(4) :: ierr
    integer(4) :: num

    if(myrank==0) then
       read(f_mdff) nitorsionBs
    endif
    call MPI_Bcast(nitorsionBs,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    if(nitorsionBs.eq.0) return

    if(myrank==0) read(f_mdff) num
    call MPI_Bcast(num,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)

    allocate(nB(npara))
    allocate(atom1B(npara,num))
    allocate(atom3B(npara,num))
    allocate(atom4B(npara,num))
#if defined(GROEXT_OPLS) || defined(GAFF)
    allocate(KdB(npara,num))
    allocate(ndB(npara,num))
    allocate(deltaB(npara,num))
#else
    allocate(KitB(npara,num))
    allocate(psi0B(npara,num))
#endif

    if(myrank==0) then
       !      read(f_mdff) topB
       read(f_mdff) nB
       read(f_mdff) atom1B
       read(f_mdff) atom3B
       read(f_mdff) atom4B
#if defined(GROEXT_OPLS) || defined(GAFF)
       read(f_mdff) KdB
       read(f_mdff) ndB
       read(f_mdff) deltaB
#else
       read(f_mdff) KitB
       read(f_mdff) psi0B
#endif
    endif
    call MPI_Bcast(nB,npara,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom1B,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom3B,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom4B,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
#if defined(GROEXT_OPLS) || defined(GAFF)
    call MPI_Bcast(KdB,npara*num,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(ndB,npara*num,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(deltaB,npara*num,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#else
    call MPI_Bcast(KitB,npara*num,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(psi0B,npara*num,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#endif
  end subroutine read_md_charmm_itorsionB
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write itorsion parameter B in mdff.bin File.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_md_charmm_itorsionB
    use device_numbers, only : f_mdff
    implicit none
    integer(4) :: num
    write(f_mdff) nitorsionBs
    if(nitorsionBs.eq.0) return
    num = size(atom1B(:,:), 2)
    write(f_mdff) num
    write(f_mdff) nB
    write(f_mdff) atom1B
    write(f_mdff) atom3B
    write(f_mdff) atom4B
#if defined(GROEXT_OPLS) || defined(GAFF)
    write(f_mdff) KdB
    write(f_mdff) ndB
    write(f_mdff) deltaB
#else
    write(f_mdff) KitB
    write(f_mdff) psi0B
#endif
  end subroutine write_md_charmm_itorsionB
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write itorsion parameter B.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_md_charmm_itorsionB
    implicit none
    integer(4) :: num
    write(*,*) '[write_md_charmm_itorsionB]'
    write(*,*) nitorsionBs
    if(nitorsionBs.eq.0) return
    num = size(atom1B(:,:), 2)
    write(*,*) num
    write(*,*) nB
    write(*,*) atom1B
    write(*,*) atom3B
    write(*,*) atom4B
#if defined(GROEXT_OPLS) || defined(GAFF)
    write(*,*) KdB
    write(*,*) ndB
    write(*,*) deltaB
#else
    write(*,*) KitB
    write(*,*) psi0B
#endif
  end subroutine write_memory_md_charmm_itorsionB
!-------------------------------------------------------------------------
  subroutine read_md_charmm_itorsion
    !fujimoto
    use device_numbers, only : f_mdff
    implicit none

    read(f_mdff) nitorsions
    if(nitorsions.eq.0) return
    allocate(atom1(nitorsions))
    allocate(atom2(nitorsions))
    allocate(atom3(nitorsions))
    allocate(atom4(nitorsions))
#if defined(GROEXT_OPLS) || defined(GAFF)
    allocate(Kd(nitorsions))
    allocate(nd(nitorsions))
    allocate(delta(nitorsions))
#else
    allocate(Kit(nitorsions))
    allocate(psi0(nitorsions))
#endif
    read(f_mdff) atom1
    read(f_mdff) atom2
    read(f_mdff) atom3
    read(f_mdff) atom4
#if defined(GROEXT_OPLS) || defined(GAFF)
    read(f_mdff) KdA
    read(f_mdff) ndA
    read(f_mdff) deltaA
#else
    read(f_mdff) Kit
    read(f_mdff) psi0
#endif
  end subroutine read_md_charmm_itorsion
!-----------------------------------------------------------------------
!>
!! \brief  Wrapper subroutine of improper_a
!! \author Zhiye Tang
!<
  subroutine add_charmm_itorsion_a()
!-----------------------------------------------------------------------
   use md_condition, only : allow_bond_breaking

   if(allow_bond_breaking<=0.0d0) then
       call add_charmm_itorsion_a_no_bond_breaking()
   else
       call add_charmm_itorsion_a_bond_breaking()
   endif

  end subroutine add_charmm_itorsion_a
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate energy and force of improper torsion.
!! \author Kensuke Iwahashi
!<
  subroutine add_charmm_itorsion_a_no_bond_breaking()
!
!      V = Kit(psi - psi0)^2
!
!         B
!         |     psi = angle between plane ABC and plane BCD
!         A      !! psi(A-B-C-D) = dihedral(A-B-C-D)
!        / \
!       C   D
!
    use trajectory_org
    use trajectory_mpi
    use atom_virial
    use forces
    use md_monitors
    use math_const, only : PI
    use param
    use mpi_tool
    use md_condition, only : allow_bond_breaking
#ifdef ONEPROC_AXIS
    use boundary, only : pbc_pair
#endif
    implicit none
    integer(4) :: j, ip, ia, ib, i0a,i0b,i0c,i0d, i0, ipar, k0
    real(8) :: a, b, s, s_s0, coef, tmp_angle
    real(8) :: tmp__x, tmp__y, tmp__z, i_norm
    real(8) :: tmp1__x, tmp1__y, tmp1__z, tmp2__x, tmp2__y, tmp2__z
    real(8) :: vBA__x, vBA__y, vBA__z, vBC__x, vBC__y, vBC__z
    real(8) :: vCD__x, vCD__y, vCD__z, vCA__x, vCA__y, vCA__z
    real(8) :: v1__x, v1__y, v1__z, v2__x, v2__y, v2__z
    real(8) :: v3__x, v3__y, v3__z
    real(8) :: Fa__x, Fa__y, Fa__z, Fb__x, Fb__y, Fb__z
    real(8) :: Fc__x, Fc__y, Fc__z, Fd__x, Fd__y, Fd__z
    real(8) :: sfact=1.0d0  ! scaling factor to potE and virial (default)
    real(8) :: wkvd(6), vd(6) ! debug
    real(8) :: v11,v22,v33
    real(8) :: wk_v11, wk_v22, wk_v33
    real(8) :: v21,v31,v32
    real(8) :: wk_v21, wk_v31, wk_v32
    real(8) :: wvr(6),vr(6) ! debug
real(8) :: ss1,ss2,ss3,r1,r2,r3,c0,c1,c2,s1,s2,s12,c
real(8) :: df,a11,a22,a33,a12,a13,a23,sx2,sy2,sz2
real(8) :: vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb4x,vb4y,vb4z

    real(8) :: Uitorsion, oUi
    real(8) :: fi0(3)
    integer(4) :: iam
    include 'mpif.h'

    call check_atombound_itorsion()

    iam = 0
    Uitorsion = 0.0d0
    wk_v11=0d0;wk_v22=0d0;wk_v33=0d0
    wk_v21=0d0;wk_v31=0d0;wk_v32=0d0

!$omp parallel default(none) &
!$omp& private(v11,v22,v33,v21,v31,v32) &
!$omp& reduction(+:wk_v11,wk_v22,wk_v33,wk_v21,wk_v31,wk_v32) &
!$omp& private(k0,i0,ipar,iam,j) &
!$omp& private(fi0,tmp_angle) &
!$omp& private(ia,ip,ib,i0a,i0b,i0c,i0d,vBA__x,vBA__y,vBA__z) &
!$omp& private(vBC__x,vBC__y,vBC__z,vCD__x,vCD__y,vCD__z) &
!$omp& private(vCA__x,vCA__y,vCA__z,v1__x,v1__y,v1__z) &
!$omp& private(v2__x,v2__y,v2__z,a,b,s,tmp__x,tmp__y,tmp__z) &
!$omp& private(tmp1__x,tmp1__y,tmp1__z,s_s0,coef) &
!$omp& private(v3__x,v3__y,v3__z,i_norm) &
!$omp& private(Fa__x,Fa__y,Fa__z,tmp2__x,tmp2__y,tmp2__z) &
!$omp& private(Fb__x,Fb__y,Fb__z,Fc__x,Fc__y,Fc__z) &
!$omp& private(Fd__x,Fd__y,Fd__z) &
!$omp& shared(m2i,i2m,nA,topA,atom2A,atom3A,atom4A) &
!$omp& shared(wkxyz,w3_f,wk_vir2) &
#if defined(GROEXT_OPLS) || defined(GAFF)
!$omp& shared(KdA,ndA,deltaA,kdB,ndB,deltaB) &
#else
!$omp& shared(psi0A,KitA,psi0B,KitB) &
#endif
!$omp& private(ss1,ss2,ss3,r1,r2,r3,c0,c1,c2,s1,s2,s12,c) &
!$omp& private(df,a11,a22,a33,a12,a13,a23,sx2,sy2,sz2) &
!$omp& private(vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb4x,vb4y,vb4z) &
!$omp& shared(nB,topB,atom1B,atom3B,atom4B,paranum) &
!$omp& shared(tag,na_per_cell) &
!$omp& shared(lsegtop,lseg_natoms,nselfseg) &
!$omp& shared(sfact) &
!$omp& reduction(+:Uitorsion)
!$  iam = omp_get_thread_num()
!$omp do
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          ia   = m2i(i0)
          ipar = paranum(ia)
          fi0(1:3) = 0d0

          do j = 1, nB(ipar)
!default: MTD
          if(ia < atom3B(ipar,j)+(ia-ipar))cycle
!         if(atom1B(ipar,j)+(ia-ipar) < atom4B(ipar,j)+(ia-ipar))cycle
             i0a = i2m(atom1B(ipar,j)+(ia-ipar))
             i0c = i2m(atom3B(ipar,j)+(ia-ipar))
             i0d = i2m(atom4B(ipar,j)+(ia-ipar))

             vb1x = wkxyz(1,i0a) - wkxyz(1,i0 )
             vb1y = wkxyz(2,i0a) - wkxyz(2,i0 )
             vb1z = wkxyz(3,i0a) - wkxyz(3,i0 )
             vb2x = wkxyz(1,i0c) - wkxyz(1,i0 )
             vb2y = wkxyz(2,i0c) - wkxyz(2,i0 )
             vb2z = wkxyz(3,i0c) - wkxyz(3,i0 )
             vb3x = wkxyz(1,i0d) - wkxyz(1,i0c)
             vb3y = wkxyz(2,i0d) - wkxyz(2,i0c)
             vb3z = wkxyz(3,i0d) - wkxyz(3,i0c)
             vb4x = wkxyz(1,i0a) - wkxyz(1,i0c)
             vb4y = wkxyz(2,i0a) - wkxyz(2,i0c)
             vb4z = wkxyz(3,i0a) - wkxyz(3,i0c)

#ifdef ONEPROC_AXIS
             call pbc_pair(vb1x,vb1y,vb1z)
             call pbc_pair(vb2x,vb2y,vb2z)
             call pbc_pair(vb3x,vb3y,vb3z)
             call pbc_pair(vb4x,vb4y,vb4z)
#endif
            ss1 = 1.0d0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z)
            ss2 = 1.0d0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z)
            ss3 = 1.0d0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z)

            r1 = sqrt(ss1);
            r2 = sqrt(ss2);
            r3 = sqrt(ss3);

            c0 = (vb1x * vb3x + vb1y * vb3y + vb1z * vb3z) * r1 * r3;
            c1 = (vb1x * vb2x + vb1y * vb2y + vb1z * vb2z) * r1 * r2;
            c2 = -(vb3x * vb2x + vb3y * vb2y + vb3z * vb2z) * r3 * r2;

            s1 = 1.0d0 - c1*c1;
            s1 = max(s1, 0.001d0)
            s1 = 1.0d0 / s1

            s2 = 1.0d0 - c2*c2;
            s2 = max(s2, 0.001d0)
            s2 = 1.0d0 / s2

            s12 = sqrt(s1*s2);
            c = (c1*c2 + c0) * s12;
            c=min(+1d0,c)
            c=max(-1d0,c)

            s = sqrt(1.0d0 - c*c)
            s = max(s, 0.001d0)

#ifdef GROEXT_OPLS
             !         force coefficient
             tmp_angle = ndB(ipar,j)*acos(c) - deltaB(ipar,j)
             coef = KdB(ipar,j) * ndB(ipar,j) * sin(tmp_angle)
#elif defined(GAFF)
             !         force coefficient
             tmp_angle = deltaB(ipar,j)*acos(c)
             coef = KdB(ipar,j) * ndB(ipar,j) * deltaB(ipar,j) &
     &              * sin(tmp_angle)
!note: gaff: Fcoef=+K*d*n*sin(n*phi)], where K=Kd, d=nd, n=delta
#else
             !         force coefficient
             s_s0 = acos(c) - psi0B(ipar,j)
             coef = -2.0d0 * KitB(ipar,j) * s_s0
#endif
            coef = coef / s;
            c = c * coef;
            s12 = s12 * coef;
            a11 = c*ss1*s1;
            a12 = -r1*r2*(c1*c*s1 + c2*s12);
            a13 = -r1*r3*s12;
            a22 = -ss2 * (2.0*c0*s12 - c*(s1+s2));
            a33 = c*ss3*s2;
            a23 = r2*r3*(c2*c*s2 + c1*s12);
            sx2  = a22*vb2x + a23*vb3x + a12*vb1x;
            sy2  = a22*vb2y + a23*vb3y + a12*vb1y;
            sz2  = a22*vb2z + a23*vb3z + a12*vb1z;

!default: MTD
            Fa__x = a12*vb2x + a13*vb3x + a11*vb1x;
            Fa__y = a12*vb2y + a13*vb3y + a11*vb1y;
            Fa__z = a12*vb2z + a13*vb3z + a11*vb1z;
            Fb__x = -sx2 - Fa__x;
            Fb__y = -sy2 - Fa__y;
            Fb__z = -sz2 - Fa__z;

             !         calculating & adding forces
             fi0(1) = fi0(1) + Fb__x
             fi0(2) = fi0(2) + Fb__y
             fi0(3) = fi0(3) + Fb__z
             !         Potential and virial are calculated in Loop A.
!default: MTD
            Fd__x = a23*vb2x + a33*vb3x + a13*vb1x;
            Fd__y = a23*vb2y + a33*vb3y + a13*vb1y;
            Fd__z = a23*vb2z + a33*vb3z + a13*vb1z;

            Fc__x = sx2 - Fd__x;
            Fc__y = sy2 - Fd__y;
            Fc__z = sz2 - Fd__z;

             v11 = (vb4x*Fa__x - vb2x*Fb__x + vb3x*Fd__x)
             v22 = (vb4y*Fa__y - vb2y*Fb__y + vb3y*Fd__y)
             v33 = (vb4z*Fa__z - vb2z*Fb__z + vb3z*Fd__z)
             v21 = (vb4y*Fa__x - vb2y*Fb__x + vb3y*Fd__x)
             v31 = (vb4z*Fa__x - vb2z*Fb__x + vb3z*Fd__x)
             v32 = (vb4z*Fa__y - vb2z*Fb__y + vb3z*Fd__y)
             wk_v11 = wk_v11 + v11
             wk_v22 = wk_v22 + v22
             wk_v33 = wk_v33 + v33
             wk_v21 = wk_v21 + v21
             wk_v31 = wk_v31 + v31
             wk_v32 = wk_v32 + v32
             w3_f(1,i0a,iam) = w3_f(1,i0a,iam) + Fa__x
             w3_f(2,i0a,iam) = w3_f(2,i0a,iam) + Fa__y
             w3_f(3,i0a,iam) = w3_f(3,i0a,iam) + Fa__z
             w3_f(1,i0c,iam) = w3_f(1,i0c,iam) + Fc__x
             w3_f(2,i0c,iam) = w3_f(2,i0c,iam) + Fc__y
             w3_f(3,i0c,iam) = w3_f(3,i0c,iam) + Fc__z
             w3_f(1,i0d,iam) = w3_f(1,i0d,iam) + Fd__x
             w3_f(2,i0d,iam) = w3_f(2,i0d,iam) + Fd__y
             w3_f(3,i0d,iam) = w3_f(3,i0d,iam) + Fd__z
#ifdef GROEXT_OPLS
             Uitorsion = Uitorsion + KdB(ipar,j)*(1d0+cos(tmp_angle))
#elif defined(GAFF)
             Uitorsion = Uitorsion + KdB(ipar,j)*(1d0+ndB(ipar,j)*cos(tmp_angle))
!note: gaff: E=K*[1+d*cos(n*phi)], where K=Kd, d=nd, n=delta
#else
             Uitorsion = Uitorsion + KitB(ipar,j) * s_s0 * s_s0
#endif
          enddo
          w3_f(1,i0,iam) = w3_f(1,i0,iam) + fi0(1)
          w3_f(2,i0,iam) = w3_f(2,i0,iam) + fi0(2)
          w3_f(3,i0,iam) = w3_f(3,i0,iam) + fi0(3)
       enddo ! i0
    enddo ! k0
!$omp end do
!$omp end parallel

    wk_p_energy = wk_p_energy + Uitorsion*sfact
    wk_vir2(1,0) = wk_vir2(1,0) + wk_v11 *sfact
    wk_vir2(2,0) = wk_vir2(2,0) + wk_v22 *sfact
    wk_vir2(3,0) = wk_vir2(3,0) + wk_v33 *sfact
    wk_vir2(4,0) = wk_vir2(4,0) + wk_v21 *sfact
    wk_vir2(5,0) = wk_vir2(5,0) + wk_v31 *sfact
    wk_vir2(6,0) = wk_vir2(6,0) + wk_v32 *sfact
#ifdef DEBUGFCE
    oUi=0d0
    call mpi_allreduce(Uitorsion*sfact,oUi,1, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,ipar)
#ifdef KCAL
    if(myrank==0)write(*,*)'Pot(itors)=',oUi *kJ_mol/4.184d0,'[kcal/mol]'
#else
    if(myrank==0)write(*,*)'Pot(itors)=',oUi *kJ_mol,'[kJ/mol]'
#endif
#endif
  end subroutine add_charmm_itorsion_a_no_bond_breaking

!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate energy and force of improper torsion.
!! \author Kensuke Iwahashi
!<
  subroutine add_charmm_itorsion_a_bond_breaking()
!
!      V = Kit(psi - psi0)^2
!
!         B
!         |     psi = angle between plane ABC and plane BCD
!         A      !! psi(A-B-C-D) = dihedral(A-B-C-D)
!        / \
!       C   D
!
    use trajectory_org
    use trajectory_mpi
    use atom_virial
    use forces
    use md_monitors
    use math_const, only : PI
    use param
    use mpi_tool
    use md_condition, only : allow_bond_breaking
#ifdef ONEPROC_AXIS
    use boundary, only : pbc_pair
#endif
    implicit none
    integer(4) :: j, ip, ia, ib, i0a,i0b,i0c,i0d, i0, ipar, k0
    real(8) :: a, b, s, s_s0, coef, tmp_angle
    real(8) :: tmp__x, tmp__y, tmp__z, i_norm
    real(8) :: tmp1__x, tmp1__y, tmp1__z, tmp2__x, tmp2__y, tmp2__z
    real(8) :: vBA__x, vBA__y, vBA__z, vBC__x, vBC__y, vBC__z
    real(8) :: vCD__x, vCD__y, vCD__z, vCA__x, vCA__y, vCA__z
    real(8) :: v1__x, v1__y, v1__z, v2__x, v2__y, v2__z
    real(8) :: v3__x, v3__y, v3__z
    real(8) :: Fa__x, Fa__y, Fa__z, Fb__x, Fb__y, Fb__z
    real(8) :: Fc__x, Fc__y, Fc__z, Fd__x, Fd__y, Fd__z
    real(8) :: sfact=1.0d0  ! scaling factor to potE and virial (default)
    real(8) :: wkvd(6), vd(6) ! debug
    real(8) :: v11,v22,v33
    real(8) :: wk_v11, wk_v22, wk_v33
    real(8) :: v21,v31,v32
    real(8) :: wk_v21, wk_v31, wk_v32
    real(8) :: wvr(6),vr(6) ! debug
real(8) :: ss1,ss2,ss3,r1,r2,r3,c0,c1,c2,s1,s2,s12,c
real(8) :: df,a11,a22,a33,a12,a13,a23,sx2,sy2,sz2
real(8) :: vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb4x,vb4y,vb4z

    real(8) :: Uitorsion, oUi
    real(8) :: fi0(3)
    integer(4) :: iam
    include 'mpif.h'

    call check_atombound_itorsion()

    iam = 0
    Uitorsion = 0.0d0
    wk_v11=0d0;wk_v22=0d0;wk_v33=0d0
    wk_v21=0d0;wk_v31=0d0;wk_v32=0d0

!$omp parallel default(none) &
!$omp& private(v11,v22,v33,v21,v31,v32) &
!$omp& reduction(+:wk_v11,wk_v22,wk_v33,wk_v21,wk_v31,wk_v32) &
!$omp& private(k0,i0,ipar,iam,j) &
!$omp& private(fi0,tmp_angle) &
!$omp& private(ia,ip,ib,i0a,i0b,i0c,i0d,vBA__x,vBA__y,vBA__z) &
!$omp& private(vBC__x,vBC__y,vBC__z,vCD__x,vCD__y,vCD__z) &
!$omp& private(vCA__x,vCA__y,vCA__z,v1__x,v1__y,v1__z) &
!$omp& private(v2__x,v2__y,v2__z,a,b,s,tmp__x,tmp__y,tmp__z) &
!$omp& private(tmp1__x,tmp1__y,tmp1__z,s_s0,coef) &
!$omp& private(v3__x,v3__y,v3__z,i_norm) &
!$omp& private(Fa__x,Fa__y,Fa__z,tmp2__x,tmp2__y,tmp2__z) &
!$omp& private(Fb__x,Fb__y,Fb__z,Fc__x,Fc__y,Fc__z) &
!$omp& private(Fd__x,Fd__y,Fd__z) &
!$omp& shared(m2i,i2m,nA,topA,atom2A,atom3A,atom4A) &
!$omp& shared(wkxyz,w3_f,wk_vir2) &
#if defined(GROEXT_OPLS) || defined(GAFF)
!$omp& shared(KdA,ndA,deltaA,kdB,ndB,deltaB) &
#else
!$omp& shared(psi0A,KitA,psi0B,KitB) &
#endif
!$omp& private(ss1,ss2,ss3,r1,r2,r3,c0,c1,c2,s1,s2,s12,c) &
!$omp& private(df,a11,a22,a33,a12,a13,a23,sx2,sy2,sz2) &
!$omp& private(vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb4x,vb4y,vb4z) &
!$omp& shared(nB,topB,atom1B,atom3B,atom4B,paranum) &
!$omp& shared(tag,na_per_cell) &
!$omp& shared(lsegtop,lseg_natoms,nselfseg) &
!$omp& shared(sfact) &
!$omp& reduction(+:Uitorsion)
!$  iam = omp_get_thread_num()
!$omp do
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          ia   = m2i(i0)
          ipar = paranum(ia)
          fi0(1:3) = 0d0

          do j = 1, nB(ipar)
!default: MTD
          if(ia < atom3B(ipar,j)+(ia-ipar))cycle
#if defined(GROEXT_OPLS) || defined(GAFF)
            if(KdB(ipar, j)<0.0d0) then
              !write(*,*) "improper removed", atom1B(ipar,j)+(ia-ipar), ia, &
              !  & atom3B(ipar,j)+(ia-ipar),atom4B(ipar,j)+(ia-ipar)
              cycle
            endif
#else
            if(KitB(ipar, j)<0.0d0) then
              !write(*,*) "improper removed", atom1B(ipar,j)+(ia-ipar), ia, &
              !  & atom3B(ipar,j)+(ia-ipar),atom4B(ipar,j)+(ia-ipar)
              cycle
            endif
#endif

             i0a = i2m(atom1B(ipar,j)+(ia-ipar))
             i0c = i2m(atom3B(ipar,j)+(ia-ipar))
             i0d = i2m(atom4B(ipar,j)+(ia-ipar))

             vb1x = wkxyz(1,i0a) - wkxyz(1,i0 )
             vb1y = wkxyz(2,i0a) - wkxyz(2,i0 )
             vb1z = wkxyz(3,i0a) - wkxyz(3,i0 )
             vb2x = wkxyz(1,i0c) - wkxyz(1,i0 )
             vb2y = wkxyz(2,i0c) - wkxyz(2,i0 )
             vb2z = wkxyz(3,i0c) - wkxyz(3,i0 )
             vb3x = wkxyz(1,i0d) - wkxyz(1,i0c)
             vb3y = wkxyz(2,i0d) - wkxyz(2,i0c)
             vb3z = wkxyz(3,i0d) - wkxyz(3,i0c)
             vb4x = wkxyz(1,i0a) - wkxyz(1,i0c)
             vb4y = wkxyz(2,i0a) - wkxyz(2,i0c)
             vb4z = wkxyz(3,i0a) - wkxyz(3,i0c)

#ifdef ONEPROC_AXIS
             call pbc_pair(vb1x,vb1y,vb1z)
             call pbc_pair(vb2x,vb2y,vb2z)
             call pbc_pair(vb3x,vb3y,vb3z)
             call pbc_pair(vb4x,vb4y,vb4z)
#endif
            ss1 = 1.0d0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z)
            ss2 = 1.0d0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z)
            ss3 = 1.0d0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z)

            r1 = sqrt(ss1);
            r2 = sqrt(ss2);
            r3 = sqrt(ss3);

            c0 = (vb1x * vb3x + vb1y * vb3y + vb1z * vb3z) * r1 * r3;
            c1 = (vb1x * vb2x + vb1y * vb2y + vb1z * vb2z) * r1 * r2;
            c2 = -(vb3x * vb2x + vb3y * vb2y + vb3z * vb2z) * r3 * r2;

            s1 = 1.0d0 - c1*c1;
            s1 = max(s1, 0.001d0)
            s1 = 1.0d0 / s1

            s2 = 1.0d0 - c2*c2;
            s2 = max(s2, 0.001d0)
            s2 = 1.0d0 / s2

            s12 = sqrt(s1*s2);
            c = (c1*c2 + c0) * s12;
            c=min(+1d0,c)
            c=max(-1d0,c)

            s = sqrt(1.0d0 - c*c)
            s = max(s, 0.001d0)

#ifdef GROEXT_OPLS
             !         force coefficient
             tmp_angle = ndB(ipar,j)*acos(c) - deltaB(ipar,j)
             coef = KdB(ipar,j) * ndB(ipar,j) * sin(tmp_angle)
#elif defined(GAFF)
             !         force coefficient
             tmp_angle = deltaB(ipar,j)*acos(c)
             coef = KdB(ipar,j) * ndB(ipar,j) * deltaB(ipar,j) &
     &              * sin(tmp_angle)
!note: gaff: Fcoef=+K*d*n*sin(n*phi)], where K=Kd, d=nd, n=delta
#else
             !         force coefficient
             s_s0 = acos(c) - psi0B(ipar,j)
             coef = -2.0d0 * KitB(ipar,j) * s_s0
#endif
            coef = coef / s;
            c = c * coef;
            s12 = s12 * coef;
            a11 = c*ss1*s1;
            a12 = -r1*r2*(c1*c*s1 + c2*s12);
            a13 = -r1*r3*s12;
            a22 = -ss2 * (2.0*c0*s12 - c*(s1+s2));
            a33 = c*ss3*s2;
            a23 = r2*r3*(c2*c*s2 + c1*s12);
            sx2  = a22*vb2x + a23*vb3x + a12*vb1x;
            sy2  = a22*vb2y + a23*vb3y + a12*vb1y;
            sz2  = a22*vb2z + a23*vb3z + a12*vb1z;

!default: MTD
            Fa__x = a12*vb2x + a13*vb3x + a11*vb1x;
            Fa__y = a12*vb2y + a13*vb3y + a11*vb1y;
            Fa__z = a12*vb2z + a13*vb3z + a11*vb1z;
            Fb__x = -sx2 - Fa__x;
            Fb__y = -sy2 - Fa__y;
            Fb__z = -sz2 - Fa__z;

             !         calculating & adding forces
             fi0(1) = fi0(1) + Fb__x
             fi0(2) = fi0(2) + Fb__y
             fi0(3) = fi0(3) + Fb__z
             !         Potential and virial are calculated in Loop A.
!default: MTD
            Fd__x = a23*vb2x + a33*vb3x + a13*vb1x;
            Fd__y = a23*vb2y + a33*vb3y + a13*vb1y;
            Fd__z = a23*vb2z + a33*vb3z + a13*vb1z;

            Fc__x = sx2 - Fd__x;
            Fc__y = sy2 - Fd__y;
            Fc__z = sz2 - Fd__z;

             v11 = (vb4x*Fa__x - vb2x*Fb__x + vb3x*Fd__x)
             v22 = (vb4y*Fa__y - vb2y*Fb__y + vb3y*Fd__y)
             v33 = (vb4z*Fa__z - vb2z*Fb__z + vb3z*Fd__z)
             v21 = (vb4y*Fa__x - vb2y*Fb__x + vb3y*Fd__x)
             v31 = (vb4z*Fa__x - vb2z*Fb__x + vb3z*Fd__x)
             v32 = (vb4z*Fa__y - vb2z*Fb__y + vb3z*Fd__y)
             wk_v11 = wk_v11 + v11
             wk_v22 = wk_v22 + v22
             wk_v33 = wk_v33 + v33
             wk_v21 = wk_v21 + v21
             wk_v31 = wk_v31 + v31
             wk_v32 = wk_v32 + v32
             w3_f(1,i0a,iam) = w3_f(1,i0a,iam) + Fa__x
             w3_f(2,i0a,iam) = w3_f(2,i0a,iam) + Fa__y
             w3_f(3,i0a,iam) = w3_f(3,i0a,iam) + Fa__z
             w3_f(1,i0c,iam) = w3_f(1,i0c,iam) + Fc__x
             w3_f(2,i0c,iam) = w3_f(2,i0c,iam) + Fc__y
             w3_f(3,i0c,iam) = w3_f(3,i0c,iam) + Fc__z
             w3_f(1,i0d,iam) = w3_f(1,i0d,iam) + Fd__x
             w3_f(2,i0d,iam) = w3_f(2,i0d,iam) + Fd__y
             w3_f(3,i0d,iam) = w3_f(3,i0d,iam) + Fd__z
#ifdef GROEXT_OPLS
             Uitorsion = Uitorsion + KdB(ipar,j)*(1d0+cos(tmp_angle))
#elif defined(GAFF)
             Uitorsion = Uitorsion + KdB(ipar,j)*(1d0+ndB(ipar,j)*cos(tmp_angle))
!note: gaff: E=K*[1+d*cos(n*phi)], where K=Kd, d=nd, n=delta
#else
             Uitorsion = Uitorsion + KitB(ipar,j) * s_s0 * s_s0
#endif
          enddo
          w3_f(1,i0,iam) = w3_f(1,i0,iam) + fi0(1)
          w3_f(2,i0,iam) = w3_f(2,i0,iam) + fi0(2)
          w3_f(3,i0,iam) = w3_f(3,i0,iam) + fi0(3)
       enddo ! i0
    enddo ! k0
!$omp end do
!$omp end parallel

    wk_p_energy = wk_p_energy + Uitorsion*sfact
    wk_vir2(1,0) = wk_vir2(1,0) + wk_v11 *sfact
    wk_vir2(2,0) = wk_vir2(2,0) + wk_v22 *sfact
    wk_vir2(3,0) = wk_vir2(3,0) + wk_v33 *sfact
    wk_vir2(4,0) = wk_vir2(4,0) + wk_v21 *sfact
    wk_vir2(5,0) = wk_vir2(5,0) + wk_v31 *sfact
    wk_vir2(6,0) = wk_vir2(6,0) + wk_v32 *sfact
#ifdef DEBUGFCE
    oUi=0d0
    call mpi_allreduce(Uitorsion*sfact,oUi,1, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,ipar)
#ifdef KCAL
    if(myrank==0)write(*,*)'Pot(itors)=',oUi *kJ_mol/4.184d0,'[kcal/mol]'
#else
    if(myrank==0)write(*,*)'Pot(itors)=',oUi *kJ_mol,'[kJ/mol]'
#endif
#endif
  end subroutine add_charmm_itorsion_a_bond_breaking

end module improper_torsion
