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
!! \brief  Module and subroutines to calculate morse-type bond potential.
!<
!>
!! \brief  Module to parse parameters of morse-type bond potential.
!! \author Kazushi Fujimoto
!<
!----------------------------------------------------------------------
module parse_bond_morse
    implicit none

    type yyparse_bond_morse
        integer(4) :: atom2
        real(8) :: Doh, aoh, Roh0
        type(yyparse_bond_morse),pointer :: next
    end type yyparse_bond_morse

    type(yyparse_bond_morse),pointer :: yyparse_bond_morse_top(:)
    integer(4),allocatable :: yyparse_bond_morse_n(:)
    integer(4) :: yyparse_bond_morse_total=0

    contains
!----------------------------------------------------------------------
    subroutine fmod_alloc_yyparse_bond_morse_top(ivalue)
        implicit none
        integer(4), intent(in) :: ivalue
        integer(4) :: i
        type(yyparse_bond_morse),pointer :: top
    allocate(yyparse_bond_morse_top(ivalue))
    do i=1,ivalue
       allocate(top)
       top%atom2 = -1
       nullify(top%next)
       yyparse_bond_morse_top(i) = top
    enddo
  end subroutine fmod_alloc_yyparse_bond_morse_top
!----------------------------------------------------------------------
  subroutine fmod_alloc_yyparse_bond_morse_n(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(yyparse_bond_morse_n(ivalue))
    yyparse_bond_morse_n = 0
  end subroutine fmod_alloc_yyparse_bond_morse_n
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to keep bond pair in liner list temporally.
!<
  subroutine add_to_bond_morse_list(atom1,atom2,Doh, aoh, Roh0)
    use trajectory_org
    use param
    use mpi_tool
    implicit none
    integer(4), intent(in) :: atom1,atom2
    real(8), intent(in) :: Doh, aoh, Roh0
    type(yyparse_bond_morse),pointer :: new1,new2,p,next
    if (atom1 > npara) then
       write(0,'(a, a,i0,a)')  'ERROR: '// &
            &    'The number of bond morse is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom2 > n) then
       write(0,'(a, a,i0,a)')  'ERROR: '// &
            &    'The number of bond morse is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    allocate(new1)
    nullify(new1%next)
    new1%atom2 = atom2
    new1%Doh = Doh
    new1%aoh = aoh
    new1%Roh0 = Roh0
    p => yyparse_bond_morse_top(atom1)
    !     insert to list, sorting from smaller
    do while (.true.)
       next => p%next
       if (.not. associated(next)) then
          p%next => new1
          yyparse_bond_morse_n(atom1) = yyparse_bond_morse_n(atom1) + 1
          yyparse_bond_morse_total = yyparse_bond_morse_total + 1
          exit
       endif
       !       doubled must be merged
       if (atom2 == next%atom2) then
          next%Doh = new1%Doh
          next%aoh = new1%aoh
          next%Roh0 = new1%Roh0
          deallocate(new1)
          exit
       else if (atom2 < next%atom2) then
          new1%next => p%next
          p%next => new1
          yyparse_bond_morse_n(atom1) = yyparse_bond_morse_n(atom1) + 1
          yyparse_bond_morse_total = yyparse_bond_morse_total + 1
          exit
       endif
       p => next
    enddo
  end subroutine add_to_bond_morse_list
!-----------------------------------------------------------------------
  subroutine add_to_mol_bond_morse_list(ivalue1, atom1, atom2, Doh, aoh, Roh0)
    use mol_info, only : para_start
    implicit none
    integer(4), intent(in) :: ivalue1
    integer(4), intent(inout) :: atom1, atom2
    real(8), intent(in) :: Doh, aoh, Roh0

    atom1 = atom1 + para_start(ivalue1)
    atom2 = atom2 + para_start(ivalue1)
    call add_to_bond_morse_list(atom1, atom2, Doh, aoh, Roh0)
    call add_to_bond_morse_list(atom2, atom1, Doh, aoh, Roh0)
  end subroutine add_to_mol_bond_morse_list

end module parse_bond_morse
!----------------------------------------------------------------------
!>
!! \brief  Module to calculate morse-type bond potential.
!! \author Kazushi Fujimoto
!<
!----------------------------------------------------------------------
module bond_morse
  use omp_lib
  use parse_bond_morse
  use subcell
  implicit none
  integer(4) :: nbondmorses=0
  integer(4) :: nbondmorse_as=0
  integer(4), allocatable :: atom1(:), atom2(:)
  real(8), allocatable :: Doh(:), aoh(:), Roh0(:)
  integer(4), allocatable :: atom2A(:,:)
  real(8), allocatable :: DohA(:,:), aohA(:,:), Roh0A(:,:)
  integer(4), allocatable :: topA(:), nA(:)

contains

!  subroutine init_bonds
!  end subroutine init_bonds

!  subroutine intra_bonds
!  end subroutine intra_bonds

!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to check whether j-atoms are communicated correctly
!! \author Yoshimichi Andoh
!<
  subroutine check_atombound_bond_morse()
!-----------------------------------------------------------------------
    use trajectory_mpi
    use param
    use md_condition
    use mpi_tool
    implicit none
    integer(4) :: k0,i0,j
    integer(4) :: iA,iB,ipar
    integer(4) :: icheck

    icheck=0

!$omp parallel default(none) &
!$omp& private(k0,i0,j,iA,iB,ipar) &
!$omp& shared(nselfseg,lsegtop,lseg_natoms) &
!$omp& shared(m2i,paranum,nA,atom2A,i2m) &
!$omp& reduction(+:icheck)
!$omp do
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          iA   = m2i(i0)
          ipar = paranum(iA)
          do iB = 1, nA(ipar)
             j = i2m(atom2A(ipar,iB)+(iA-ipar))
             if(j.eq.-1) icheck = icheck -1
          enddo
       enddo
    enddo
!$omp end do
!$omp end parallel

    if(icheck.le.-1) then
       write(*,*) 'Err[add_bond_morse] &
            &   There is a particle outside the area.'
       write(*,*) 'Myrank,mdstep=',myrank,mdstep
       call modylas_abort()
    end if
  end subroutine check_atombound_bond_morse
!-----------------------------------------------------------------------
  subroutine set_md_bond_morse
    use param
    implicit none
    integer(4) :: i,j,k, maxi(1)
    type(yyparse_bond_morse),pointer :: p,freed
    nbondmorse_as= yyparse_bond_morse_total
    if (yyparse_bond_morse_total == 0)  return
    maxi = maxloc(yyparse_bond_morse_n)
    maxi = yyparse_bond_morse_n(maxi)
    allocate(atom2A(npara, maxi(1)))
    allocate(DohA(npara, maxi(1)))
    allocate(aohA(npara, maxi(1)))
    allocate(Roh0A(npara, maxi(1)))
    allocate(nA(npara))
    atom2A = 0
    DohA(:,:) = 0.0d0
    aohA(:,:) = 0.0d0
    Roh0A(:,:) = 0.0d0
    nA = 0
    !     copying special pair entries

    do i=1, npara
       nA(i) = yyparse_bond_morse_n(i)
       p => yyparse_bond_morse_top(i)
       do j=1,yyparse_bond_morse_n(i)
          p => p%next
          atom2A(i,j) = p%atom2
          DohA(i,j) = p%Doh
          aohA(i,j) = p%aoh
          Roh0A(i,j) = p%Roh0
       enddo
    enddo

    deallocate(yyparse_bond_morse_top, yyparse_bond_morse_n)
  end subroutine set_md_bond_morse
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to read bond parameter.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_md_bond_morse
    use trajectory_org
    use device_numbers, only : f_mdff
    use param
    use mpi_tool
    implicit none
    include 'mpif.h'
    integer(4) :: ierr
    integer(4) :: num

    if(myrank==0 )read(f_mdff) nbondmorse_as
    call MPI_Bcast(nbondmorse_as,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    if(nbondmorse_as.eq.0) return

    if(myrank==0) read(f_mdff) num
    call MPI_Bcast(num,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)

    allocate(nA(npara))
    allocate(atom2A(npara,num))
    allocate(DohA(npara,num))
    allocate(aohA(npara,num))
    allocate(Roh0A(npara,num))

    if(myrank==0) then
       read(f_mdff) nA
       read(f_mdff) atom2A
       read(f_mdff) DohA
       read(f_mdff) aohA
       read(f_mdff) Roh0A
    endif
    call MPI_Bcast(nA,npara,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom2A,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(DohA,npara*num,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(aohA,npara*num,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(Roh0A,npara*num,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  end subroutine read_md_bond_morse
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write bond parameter in mdff.bin File.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_md_bond_morse
    use device_numbers, only : f_mdff
    implicit none
    integer(4) :: num
    write(f_mdff) nbondmorse_as
    if(nbondmorse_as.eq.0) return
    num = size(atom2A(:,:), 2)
    write(f_mdff) num
    write(f_mdff) nA
    write(f_mdff) atom2A
    write(f_mdff) DohA
    write(f_mdff) aohA
    write(f_mdff) Roh0A
  end subroutine write_md_bond_morse
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write bond parameter.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_md_bond_morse
    implicit none
    integer(4) :: num
    write(*,*) '[write_md_charmm_a_bond]'
    write(*,*) nbondmorse_as
    if(nbondmorse_as.eq.0) return
    num = size(atom2A(:,:), 2)
    write(*,*) num
    write(*,*) nA
    write(*,*) atom2A
    write(*,*) DohA
    write(*,*) aohA
    write(*,*) Roh0A
  end subroutine write_memory_md_bond_morse
!-------------------------------------------------------------------------
!!!!  subroutine read_md_bond_morse
!!!!    !fujimoto
!!!!    use device_numbers, only : f_mdff
!!!!    implicit none
!!!!
!!!!    read(f_mdff) nbondmorses
!!!!    if(nbondmorses.eq.0) return
!!!!    allocate(atom1(nbondmorses))
!!!!    allocate(atom2(nbondmorses))
!!!!    allocate(Doh(nbondmorses))
!!!!    allocate(aoh(nbondmorses))
!!!!    allocate(Roh0(nbondmorses))
!!!!    read(f_mdff) atom1
!!!!    read(f_mdff) atom2
!!!!    read(f_mdff) Doh
!!!!    read(f_mdff) aoh
!!!!    read(f_mdff) Roh0
!!!!  end subroutine read_md_bond_morse
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate energy and force of bond.
!! \author Kazushi Fujimoto
!<
  subroutine add_bond_morse()
!----------------------------------------------------------------------
!
!      V = Kb(b - b0)^2
!
!      A ----- B
!          b
!
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
    integer(4) :: iA, iB, i0, j0, ipar, k0
    real(8) :: aoh, Doh, b, b_b0, Fb__x, Fb__y, Fb__z, drAB__x, drAB__y, drAB__z
    real(8) :: v11,v22,v33
    real(8) :: wk_v11, wk_v22, wk_v33, Ubondmorse
    real(8) :: sfact=1.0d0  ! scaling factor to potE and virial (default)
    real(8) :: wkvd(6), vd(6) ! debug
    real(8) :: v21,v31,v32
    real(8) :: wk_v21, wk_v31, wk_v32
    real(8) :: wvr(6),vr(6) ! debug
    real(8) :: fi0(3), oUb
    integer(4) :: iam
    include 'mpif.h'

    call check_atombound_bond_morse()

    iam = 0
    Ubondmorse = 0.0d0
    wk_v11=0d0;wk_v22=0d0;wk_v33=0d0
    wk_v21=0d0;wk_v31=0d0;wk_v32=0d0

!$omp parallel default(none) &
!$omp& private(k0,i0,j0,ipar) &
!$omp& private(iam,iA,iB,drAB__x,drAB__y,drAB__z) &
!$omp& private(b,b_b0,Fb__x,Fb__y,Fb__z) &
!$omp& private(Doh, aoh) &
!$omp& private(v11,v22,v33,fi0) &
!$omp& shared(m2i,i2m,nA,atom2A,topA,wkxyz) &
!$omp& shared(DohA, aohA, Roh0A, w3_f,wk_vir2,paranum) &
!$omp& shared(tag,na_per_cell) &
!$omp& shared(lsegtop,lseg_natoms,nselfseg) &
!$omp& shared(sfact) &
!$omp& private(v21,v31,v32) &
!$omp& reduction(+:wk_v21,wk_v31,wk_v32) &
!$omp& reduction(+:Ubondmorse) &
!$omp& reduction(+:wk_v11,wk_v22,wk_v33)
!$    iam = omp_get_thread_num()
!$omp do
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          iA   = m2i(i0)
          ipar = paranum(iA)
          fi0(1:3)=0d0
          do iB = 1, nA(ipar)
!default: MTD
          if(iA < atom2A(ipar,iB)+(iA-ipar))cycle ! avoid double count
             j0 = i2m(atom2A(ipar,iB)+(iA-ipar))
             drAB__x = wkxyz(1,i0) - wkxyz(1,j0)
             drAB__y = wkxyz(2,i0) - wkxyz(2,j0)
             drAB__z = wkxyz(3,i0) - wkxyz(3,j0)
#ifdef ONEPROC_AXIS
             call pbc_pair(drAB__x,drAB__y,drAB__z)
#endif
             b = sqrt(drAB__x*drAB__x+drAB__y*drAB__y+drAB__z*drAB__z)
             b_b0 = b - Roh0A(ipar,iB)

             Doh = DohA(ipar, iB)
             aoh = aohA(ipar, iB)

             Fb__x =   drAB__x * 2.0d0 * aoh * Doh * exp(-aoh * b_b0) * &
             &         (exp(-aoh * b_b0) - 1) / b
             Fb__y =   drAB__y * 2.0d0 * aoh * Doh * exp(-aoh * b_b0) * &
             &         (exp(-aoh * b_b0) - 1) / b
             Fb__z =   drAB__z * 2.0d0 * aoh * Doh * exp(-aoh * b_b0) * &
             &         (exp(-aoh * b_b0) - 1) / b

             fi0(1) = fi0(1) + Fb__x
             fi0(2) = fi0(2) + Fb__y
             fi0(3) = fi0(3) + Fb__z
!default: MTD
             w3_f(1,j0,iam) = w3_f(1,j0,iam) - Fb__x
             w3_f(2,j0,iam) = w3_f(2,j0,iam) - Fb__y
             w3_f(3,j0,iam) = w3_f(3,j0,iam) - Fb__z
             Ubondmorse = Ubondmorse + Doh * (1 - exp(-aoh * b_b0))**2
             v11 = drAB__x * Fb__x * sfact !0.5d0  !virial11
             v22 = drAB__y * Fb__y * sfact !0.5d0
             v33 = drAB__z * Fb__z * sfact !0.5d0
             wk_v11 = wk_v11 + v11
             wk_v22 = wk_v22 + v22
             wk_v33 = wk_v33 + v33
             v21 = drAB__y * Fb__x * sfact !0.5d0
             v31 = drAB__z * Fb__x * sfact !0.5d0
             v32 = drAB__z * Fb__y * sfact !0.5d0
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

    wk_p_energy = wk_p_energy + Ubondmorse * sfact 
    wk_vir2(1,0) = wk_vir2(1,0) + wk_v11
    wk_vir2(2,0) = wk_vir2(2,0) + wk_v22
    wk_vir2(3,0) = wk_vir2(3,0) + wk_v33
    wk_vir2(4,0) = wk_vir2(4,0) + wk_v21
    wk_vir2(5,0) = wk_vir2(5,0) + wk_v31
    wk_vir2(6,0) = wk_vir2(6,0) + wk_v32
#ifdef DEBUGFCE
    oUb=0d0
    call mpi_allreduce(Ubondmorse * sfact, oUb, 1, &
                       mpi_double_precision, mpi_sum, mpi_comm_world, ipar)
    if(myrank==0) write(*,*) 'Pot(bondmorse) =', oUb *kJ_mol,'[kJ/mol]'
!debug molp
    wkvd(1)=wk_v11; wkvd(2)=wk_v22; wkvd(3)=wk_v33
    wkvd(4)=wk_v21; wkvd(5)=wk_v31; wkvd(6)=wk_v32
    call mpi_allreduce(wkvd,vd,6, &
    &                  mpi_double_precision, mpi_sum, mpi_comm_world, ipar)
    if(myrank==0) write(*,*) 'vir(bondmorse) =', vd(1:6), sum(vd(1:3))
#endif
  end subroutine add_bond_morse

end module bond_morse

