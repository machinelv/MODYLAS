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
!! \brief  Module and subroutines to calculate bond potential.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to parse parameters of bond potential.
!! \author Kensuke Iwahashi
!<
module parse_bond
  implicit none

  type yyparse_bond
     integer(4) :: atom2
     real(8) :: Kb, b0
     type(yyparse_bond),pointer :: next
  end type yyparse_bond
  type(yyparse_bond),pointer :: yyparse_bond_top(:)
  integer(4),allocatable :: yyparse_bond_n(:)
  integer(4) :: yyparse_bond_total=0

contains
!----------------------------------------------------------------------
!>
!! \brief  Subroutines to read inputs relating to bond potential.
!! \author Kensuke Iwahashi
!<
  subroutine fmod_alloc_ibond(ivalue)
    use shake_rattle_roll, only : ibond
    integer(4), intent(in) :: ivalue
    allocate(ibond(2,ivalue))
  end subroutine fmod_alloc_ibond
!----------------------------------------------------------------------
!>
!! \brief  Subroutines to read inputs relating to bond potential.
!! \author Kensuke Iwahashi
!<
  subroutine fmod_alloc_yyparse_bond_top(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    integer(4) :: i
    type(yyparse_bond),pointer :: top
    allocate(yyparse_bond_top(ivalue))
    do i=1,ivalue
       allocate(top)
       top%atom2 = -1
       nullify(top%next)
       yyparse_bond_top(i) = top
    enddo
  end subroutine fmod_alloc_yyparse_bond_top
!----------------------------------------------------------------------
!>
!! \brief  Subroutines to read inputs relating to bond potential.
!! \author Kensuke Iwahashi
!<
  subroutine fmod_alloc_yyparse_bond_n(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(yyparse_bond_n(ivalue))
    yyparse_bond_n = 0
  end subroutine fmod_alloc_yyparse_bond_n
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to keep bond pair in liner list temporally.
!! \author Kensuke Iwahashi
!<
  subroutine add_to_bond_list(atom1,atom2,Kb,b0)
    use trajectory_org
    use param
    use mpi_tool
    implicit none
    integer(4), intent(in) :: atom1,atom2
    real(8), intent(in) :: Kb,b0
    type(yyparse_bond),pointer :: new1,new2,p,next
    if (atom1 > npara) then
       write(0,'(a, a,i0,a)')  'ERROR: '// &
            &    'The number of bond is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom2 > n) then
       write(0,'(a, a,i0,a)')  'ERROR: '// &
            &    'The number of bond is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    allocate(new1)
    nullify(new1%next)
    new1%atom2 = atom2
    new1%Kb = Kb
    new1%b0 = b0
    p => yyparse_bond_top(atom1)
    !     insert to list, sorting from smaller
    do while (.true.)
       next => p%next
       if (.not. associated(next)) then
          p%next => new1
          yyparse_bond_n(atom1) = yyparse_bond_n(atom1) + 1
          yyparse_bond_total = yyparse_bond_total + 1
          exit
       endif
       !       doubled must be merged
       if (atom2 == next%atom2) then
          next%Kb = new1%Kb
          next%b0 = new1%b0
          deallocate(new1)
          exit
       else if (atom2 < next%atom2) then
          new1%next => p%next
          p%next => new1
          yyparse_bond_n(atom1) = yyparse_bond_n(atom1) + 1
          yyparse_bond_total = yyparse_bond_total + 1
          exit
       endif
       p => next
    enddo
  end subroutine add_to_bond_list
!-----------------------------------------------------------------------
!>
!! \brief  Subroutines to read inputs relating to bond potential.
!! \author Kensuke Iwahashi
!<
  subroutine add_to_mol_bond_list(ivalue1, atom1, atom2, kb, b0)
    use mol_info, only : para_start
    implicit none
    integer(4), intent(in) :: ivalue1
    integer(4), intent(inout) :: atom1, atom2
    real(8), intent(in) :: kb, b0

    atom1 = atom1 + para_start(ivalue1)
    atom2 = atom2 + para_start(ivalue1)
    call add_to_bond_list(atom1, atom2, kb, b0)
    call add_to_bond_list(atom2, atom1, kb, b0)
  end subroutine add_to_mol_bond_list

end module parse_bond
!----------------------------------------------------------------------
!>
!! \brief  Module to calculate bond potential.
!! \author Kensuke Iwahashi
!<
!----------------------------------------------------------------------
module bond
  use omp_lib
  use parse_bond
  use subcell
  implicit none
  integer(4) :: nbonds=0
  integer(4) :: nbond_as=0
  integer(4), allocatable :: atom1(:), atom2(:)
  real(8), allocatable :: Kb(:), b0(:)
  integer(4), allocatable :: atom2A(:,:)
  real(8), allocatable :: KbA(:,:), b0A(:,:)
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
  subroutine check_atombound_bond()
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
!$omp& shared(wkxyz) &
!$omp& reduction(+:icheck)
!$omp do
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          iA   = m2i(i0)
          ipar = paranum(iA)
          do iB = 1, nA(ipar)
             j = i2m(atom2A(ipar,iB)+(iA-ipar))
              if(j.eq.-1) then
                  icheck = icheck -1
                  write(6,*) "wth"
                  write(6,*) "bond is not in local domain",iA,atom2A(ipar,iB)+(iA-ipar)
                  write(6,*) "coor 0 ", wkxyz(1,i0), wkxyz(2,i0), wkxyz(3,i0)
                  flush(6)
              endif
          enddo
       enddo
    enddo
!$omp end do
!$omp end parallel

    if(icheck.le.-1) then
       write(*,*) 'Err[add_charmm_bond_a] &
            &   There is a particle outside the area.'
       write(*,*) 'Myrank,mdstep=',myrank,mdstep
       call modylas_abort()
    end if

  end subroutine check_atombound_bond
!-----------------------------------------------------------------------
  subroutine set_md_charmm_a_bond
    use param
    implicit none
    integer(4) :: i,j,k, maxi(1)
    type(yyparse_bond),pointer :: p,freed
    nbond_as= yyparse_bond_total
    if (yyparse_bond_total == 0)  return
    maxi = maxloc(yyparse_bond_n)
    maxi = yyparse_bond_n(maxi)
    allocate(atom2A(npara, maxi(1)))
    allocate(KbA(npara, maxi(1)))
    allocate(b0A(npara, maxi(1)))
    allocate(nA(npara))
    atom2A = 0
    KbA = 0.0d0
    b0A = 0.0d0
    nA = 0
    !     copying special pair entries

    do i=1, npara
       nA(i) = yyparse_bond_n(i)
       p => yyparse_bond_top(i)
       do j=1,yyparse_bond_n(i)
          p => p%next
          atom2A(i,j) = p%atom2
          KbA(i,j) = p%Kb
          b0A(i,j) = p%b0
       enddo
    enddo

    deallocate(yyparse_bond_top, yyparse_bond_n)
  end subroutine set_md_charmm_a_bond
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to read bond parameter.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_md_charmm_a_bond
    use trajectory_org
    use device_numbers, only : f_mdff
    use param
    use mpi_tool
    implicit none
    include 'mpif.h'
    integer(4) :: ierr
    integer(4) :: num

    if(myrank==0 )read(f_mdff) nbond_as
    call MPI_Bcast(nbond_as,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    if(nbond_as.eq.0) return

    if(myrank==0) read(f_mdff) num
    call MPI_Bcast(num,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)

    allocate(nA(npara))
    allocate(atom2A(npara,num))
    allocate(KbA(npara,num))
    allocate(b0A(npara,num))

    if(myrank==0) then
       read(f_mdff) nA
       read(f_mdff) atom2A
       read(f_mdff) KbA
       read(f_mdff) b0A
    endif
    call MPI_Bcast(nA,npara,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom2A,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(KbA,npara*num,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(b0A,npara*num,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  end subroutine read_md_charmm_a_bond
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write bond parameter in mdff.bin File.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_md_charmm_a_bond
    use device_numbers, only : f_mdff
    implicit none
    integer(4) :: num
    write(f_mdff) nbond_as
    if(nbond_as.eq.0) return
    num = size(atom2A(:,:), 2)
    write(f_mdff) num
    write(f_mdff) nA
    write(f_mdff) atom2A
    write(f_mdff) KbA
    write(f_mdff) b0A
  end subroutine write_md_charmm_a_bond
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write bond parameter.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_md_charmm_a_bond
    implicit none
    integer(4) :: num
    write(*,*) '[write_md_charmm_a_bond]'
    write(*,*) nbond_as
    if(nbond_as.eq.0) return
    num = size(atom2A(:,:), 2)
    write(*,*) num
    write(*,*) nA
    write(*,*) atom2A
    write(*,*) KbA
    write(*,*) b0A
  end subroutine write_memory_md_charmm_a_bond
!-------------------------------------------------------------------------
  subroutine read_md_charmm_bond
    !fujimoto
    use device_numbers, only : f_mdff
    implicit none

    read(f_mdff) nbonds
    if(nbonds.eq.0) return
    allocate(atom1(nbonds))
    allocate(atom2(nbonds))
    allocate(Kb(nbonds))
    allocate(b0(nbonds))
    read(f_mdff) atom1
    read(f_mdff) atom2
    read(f_mdff) Kb
    read(f_mdff) b0
  end subroutine read_md_charmm_bond
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate energy and force of bond.
!! \author Kensuke Iwahashi
!<
  subroutine add_charmm_bond_a()
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
    real(8) :: b, b_b0, Fb__x, Fb__y, Fb__z, drAB__x, drAB__y, drAB__z
    real(8) :: v11,v22,v33
    real(8) :: wk_v11, wk_v22, wk_v33, Ubond
    real(8) :: sfact=1.0d0  ! scaling factor to potE and virial (default)
    real(8) :: wkvd(6), vd(6) ! debug
    real(8) :: v21,v31,v32
    real(8) :: wk_v21, wk_v31, wk_v32
    real(8) :: wvr(6),vr(6) ! debug
    real(8) :: fi0(3), oUb
    integer(4) :: iam
    include 'mpif.h'

    call check_atombound_bond()

    iam = 0
    Ubond = 0.0d0
    wk_v11=0d0;wk_v22=0d0;wk_v33=0d0
    wk_v21=0d0;wk_v31=0d0;wk_v32=0d0

!$omp parallel default(none) &
!$omp& private(k0,i0,j0,ipar) &
!$omp& private(iam,iA,iB,drAB__x,drAB__y,drAB__z) &
!$omp& private(b,b_b0,Fb__x,Fb__y,Fb__z) &
!$omp& private(v11,v22,v33,fi0) &
!$omp& shared(m2i,i2m,nA,atom2A,topA,wkxyz) &
!$omp& shared(b0A,KbA,w3_f,wk_vir2,paranum) &
!$omp& shared(tag,na_per_cell) &
!$omp& shared(lsegtop,lseg_natoms,nselfseg) &
!$omp& shared(sfact) &
!$omp& private(v21,v31,v32) &
!$omp& reduction(+:wk_v21,wk_v31,wk_v32) &
!$omp& reduction(+:Ubond) &
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
             b_b0 = b - b0A(ipar,iB)
             Fb__x = - drAB__x * 2.0d0 * KbA(ipar,iB) * b_b0 / b
             Fb__y = - drAB__y * 2.0d0 * KbA(ipar,iB) * b_b0 / b
             Fb__z = - drAB__z * 2.0d0 * KbA(ipar,iB) * b_b0 / b
             fi0(1)=fi0(1)+Fb__x
             fi0(2)=fi0(2)+Fb__y
             fi0(3)=fi0(3)+Fb__z
!default: MTD
             w3_f(1,j0,iam) = w3_f(1,j0,iam) - Fb__x
             w3_f(2,j0,iam) = w3_f(2,j0,iam) - Fb__y
             w3_f(3,j0,iam) = w3_f(3,j0,iam) - Fb__z
             Ubond = Ubond + KbA(ipar,iB) * b_b0 * b_b0
             v11 = drAB__x * Fb__x * sfact !0.5d0  !virial11
             v22 = drAB__y * Fb__y * sfact !0.5d0
             v33 = drAB__z * Fb__z * sfact !0.5d0
             wk_v11=wk_v11+v11
             wk_v22=wk_v22+v22
             wk_v33=wk_v33+v33
             v21 = drAB__y * Fb__x * sfact !0.5d0
             v31 = drAB__z * Fb__x * sfact !0.5d0
             v32 = drAB__z * Fb__y * sfact !0.5d0
             wk_v21=wk_v21+v21
             wk_v31=wk_v31+v31
             wk_v32=wk_v32+v32
          enddo
          w3_f(1,i0,iam) = w3_f(1,i0,iam) + fi0(1)
          w3_f(2,i0,iam) = w3_f(2,i0,iam) + fi0(2)
          w3_f(3,i0,iam) = w3_f(3,i0,iam) + fi0(3)
       enddo ! i0
    enddo ! k0
!$omp end do
!$omp end parallel

    wk_p_energy = wk_p_energy + Ubond * sfact 
    wk_vir2(1,0) = wk_vir2(1,0) + wk_v11
    wk_vir2(2,0) = wk_vir2(2,0) + wk_v22
    wk_vir2(3,0) = wk_vir2(3,0) + wk_v33
    wk_vir2(4,0) = wk_vir2(4,0) + wk_v21
    wk_vir2(5,0) = wk_vir2(5,0) + wk_v31
    wk_vir2(6,0) = wk_vir2(6,0) + wk_v32

#ifdef DEBUGFCE
    oUb=0d0
    call mpi_allreduce(Ubond*sfact,oUb,1, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,ipar)
#ifdef KCAL
    if(myrank==0) write(*,*) 'Pot(bond) =', oUb *kJ_mol/4.184d0,'[kcal/mol]'
#else
    if(myrank==0) write(*,*) 'Pot(bond) =', oUb *kJ_mol,'[kJ/mol]'
#endif
#endif

  end subroutine add_charmm_bond_a

end module bond
