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
!! \brief  Module and subroutines to calculate morse-type bond potential. (part 2)
!<
!>
!! \brief  Module to parse parameters of morse-type bond potential. (part 2)
!! \author Kazushi Fujimoto
!<
!----------------------------------------------------------------------
module parse_bond_morse2
    implicit none

    type yyparse_bond_morse2
        integer(4) :: atom2
         real(8) :: D, am, bm, b0 
        type(yyparse_bond_morse2),pointer :: next
    end type yyparse_bond_morse2

    type(yyparse_bond_morse2),pointer :: yyparse_bond_morse2_top(:)
    integer(4),allocatable :: yyparse_bond_morse2_n(:)
    integer(4) :: yyparse_bond_morse2_total=0

    contains
!----------------------------------------------------------------------
    subroutine fmod_alloc_yyparse_bond_morse2_top(ivalue)
        implicit none
        integer(4), intent(in) :: ivalue
        integer(4) :: i
        type(yyparse_bond_morse2),pointer :: top
    allocate(yyparse_bond_morse2_top(ivalue))
    do i=1,ivalue
       allocate(top)
       top%atom2 = -1
       nullify(top%next)
       yyparse_bond_morse2_top(i) = top
    enddo
  end subroutine fmod_alloc_yyparse_bond_morse2_top
!----------------------------------------------------------------------
  subroutine fmod_alloc_yyparse_bond_morse2_n(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(yyparse_bond_morse2_n(ivalue))
    yyparse_bond_morse2_n = 0
  end subroutine fmod_alloc_yyparse_bond_morse2_n
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to keep bond pair in liner list temporally.
!<
  subroutine add_to_bond_morse2_list&
             (atom1,atom2,D, am, bm, b0)
    use trajectory_org
    use param
    use mpi_tool
    implicit none
    integer(4), intent(in) :: atom1,atom2
    real(8), intent(in) :: D, am, bm, b0
    type(yyparse_bond_morse2),pointer :: new1,new2,p,next
    if (atom1 > npara) then
       write(0,'(a, a,i0,a)')  'ERROR: '// &
            &    'The number of bond morse2 is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom2 > n) then
       write(0,'(a, a,i0,a)')  'ERROR: '// &
            &    'The number of bond morse2 is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    allocate(new1)
    nullify(new1%next)
    new1%atom2 = atom2
    new1%D = D
    new1%am = am
    new1%bm = bm
    new1%b0 = b0
    p => yyparse_bond_morse2_top(atom1)
    !     insert to list, sorting from smaller
    do while (.true.)
       next => p%next
       if (.not. associated(next)) then
          p%next => new1
          yyparse_bond_morse2_n(atom1) = yyparse_bond_morse2_n(atom1) + 1
          yyparse_bond_morse2_total = yyparse_bond_morse2_total + 1
          exit
       endif
       !       doubled must be merged
       if (atom2 == next%atom2) then
          next%D  = new1%D
          next%am = new1%am
          next%bm = new1%bm
          next%b0 = new1%b0
          deallocate(new1)
          exit
       else if (atom2 < next%atom2) then
          new1%next => p%next
          p%next => new1
          yyparse_bond_morse2_n(atom1) = yyparse_bond_morse2_n(atom1) + 1
          yyparse_bond_morse2_total = yyparse_bond_morse2_total + 1
          exit
       endif
       p => next
    enddo
  end subroutine add_to_bond_morse2_list
!-----------------------------------------------------------------------
  subroutine add_to_mol_bond_morse2_list &
             (ivalue1, atom1, atom2, D, am, bm, b0)
    use mol_info, only : para_start
    implicit none
    integer(4), intent(in) :: ivalue1
    integer(4), intent(inout) :: atom1, atom2
    real(8), intent(in) :: D, am, bm, b0

    atom1 = atom1 + para_start(ivalue1)
    atom2 = atom2 + para_start(ivalue1)
    call add_to_bond_morse2_list(atom1, atom2, D, am, bm, b0)
    call add_to_bond_morse2_list(atom2, atom1, D, am, bm, b0)
  end subroutine add_to_mol_bond_morse2_list

end module parse_bond_morse2
!----------------------------------------------------------------------
!>
!! \brief  Module to calculate morse-type bond potential. (part 2)
!! \author Kazushi Fujimoto
!<
!----------------------------------------------------------------------
module bond_morse2
  use omp_lib
  use parse_bond_morse2
  use subcell
  implicit none
  integer(4) :: nbondmorse2s=0
  integer(4) :: nbondmorse2_as=0
  integer(4), allocatable :: atom1(:), atom2(:)
  real(8), allocatable :: D(:), am(:), bm(:), bm0(:)
  integer(4), allocatable :: atom2A(:,:)
  real(8), allocatable :: DA(:,:), amA(:,:), bmA(:,:), b0A(:,:)
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
  subroutine check_atombound_bond_morse2()
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
    
    if(allow_bond_breaking<=0.0d0) then

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

    else

!$omp parallel default(none) &
!$omp& private(k0,i0,j,iA,iB,ipar) &
!$omp& shared(nselfseg,lsegtop,lseg_natoms) &
!$omp& shared(m2i,paranum,nA,atom2A,i2m) &
!$omp& reduction(+:icheck) &
!$omp& shared(DA)
!$omp do
      do k0=1,nselfseg
        do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
           iA   = m2i(i0)
           ipar = paranum(iA)
           do iB = 1, nA(ipar)
              j = i2m(atom2A(ipar,iB)+(iA-ipar))
              if(DA(ipar,iB)>=0.0d0 .and. j.eq.-1) icheck = icheck -1
           enddo
        enddo
     enddo
 !$omp end do
 !$omp end parallel
     endif

    if(icheck.le.-1) then
       write(*,*) 'Err[add_bond_morse2] &
            &   There is a particle outside the area.'
       write(*,*) 'Myrank,mdstep=',myrank,mdstep
       call modylas_abort()
    end if
  end subroutine check_atombound_bond_morse2
!-----------------------------------------------------------------------
  subroutine set_md_bond_morse2
    use param
    implicit none
    integer(4) :: i,j,k, maxi(1)
    type(yyparse_bond_morse2),pointer :: p,freed
    nbondmorse2_as= yyparse_bond_morse2_total
    if (yyparse_bond_morse2_total == 0)  return
    maxi = maxloc(yyparse_bond_morse2_n)
    maxi = yyparse_bond_morse2_n(maxi)
    allocate(atom2A(npara, maxi(1)))
    allocate(DA(npara, maxi(1)))
    allocate(amA(npara, maxi(1)))
    allocate(bmA(npara, maxi(1)))
    allocate(b0A(npara, maxi(1)))
    allocate(nA(npara))
    atom2A = 0
    DA(:,:) = 0.0d0
    amA(:,:) = 0.0d0
    bmA(:,:) = 0.0d0
    b0A(:,:) = 0.0d0
    nA = 0
    !     copying special pair entries

    do i=1, npara
       nA(i) = yyparse_bond_morse2_n(i)
       p => yyparse_bond_morse2_top(i)
       do j=1,yyparse_bond_morse2_n(i)
          p => p%next
          atom2A(i,j) = p%atom2
          DA(i,j)     = p%D
          amA(i,j)    = p%am
          bmA(i,j)    = p%bm
          b0A(i,j)    = p%b0
       enddo
    enddo

    deallocate(yyparse_bond_morse2_top, yyparse_bond_morse2_n)
  end subroutine set_md_bond_morse2
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to read bond parameter.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_md_bond_morse2
    use trajectory_org
    use device_numbers, only : f_mdff
    use param
    use mpi_tool
    implicit none
    include 'mpif.h'
    integer(4) :: ierr
    integer(4) :: num

    if(myrank==0 )read(f_mdff) nbondmorse2_as
    call MPI_Bcast(nbondmorse2_as,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    if(nbondmorse2_as.eq.0) return

    if(myrank==0) read(f_mdff) num
    call MPI_Bcast(num,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)

    allocate(nA(npara))
    allocate(atom2A(npara,num))
    allocate(DA(npara,num))
    allocate(amA(npara,num))
    allocate(bmA(npara,num))
    allocate(b0A(npara,num))

    if(myrank==0) then
       read(f_mdff) nA
       read(f_mdff) atom2A
       read(f_mdff) DA
       read(f_mdff) amA
       read(f_mdff) bmA
       read(f_mdff) b0A
    endif
    call MPI_Bcast(nA,npara,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom2A,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(DA,npara*num,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(amA,npara*num,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(bmA,npara*num,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(b0A,npara*num,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  end subroutine read_md_bond_morse2
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write bond parameter in mdff.bin File.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_md_bond_morse2
    use device_numbers, only : f_mdff
    implicit none
    integer(4) :: num
    write(f_mdff) nbondmorse2_as
    if(nbondmorse2_as.eq.0) return
    num = size(atom2A(:,:), 2)
    write(f_mdff) num
    write(f_mdff) nA
    write(f_mdff) atom2A
    write(f_mdff) DA
    write(f_mdff) amA
    write(f_mdff) bmA
    write(f_mdff) b0A
  end subroutine write_md_bond_morse2
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write bond parameter.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_md_bond_morse2
    implicit none
    integer(4) :: num
    write(*,*) '[write_md_bond_morse2]'
    write(*,*) nbondmorse2_as
    if(nbondmorse2_as.eq.0) return
    num = size(atom2A(:,:), 2)
    write(*,*) num
    write(*,*) nA
    write(*,*) atom2A
    write(*,*) DA
    write(*,*) amA
    write(*,*) bmA
    write(*,*) b0A
  end subroutine write_memory_md_bond_morse2
!-----------------------------------------------------------------------
!>
!! \brief  Wrapper subroutine of bond_morse2
!! \author Zhiye Tang
!<
  subroutine add_bond_morse2()
   use md_condition, only : allow_bond_breaking
   if(allow_bond_breaking<=0.0d0) then
       call add_bond_morse2_no_bond_breaking()
   else
       call add_bond_morse2_bond_breaking()
   endif
  end subroutine add_bond_morse2
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate energy and force of bond.
!! \author Kazushi Fujimoto, Zhiye Tang
!<
  subroutine add_bond_morse2_no_bond_breaking()
!----------------------------------------------------------------------
!
!      V = D*[1-exp{-a(b-b0)-b(b-b0)**2}]**2
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
    use md_condition, only : allow_bond_breaking
#ifdef ONEPROC_AXIS
    use boundary, only : pbc_pair
#endif
    implicit none
    integer(4) :: iA, iB, i0, j0, ipar, k0
    real(8) :: b, b_b0, Fb__x, Fb__y, Fb__z, drAB__x, drAB__y, drAB__z
    real(8) :: v11,v22,v33
    real(8) :: wD, wam, wbm, wb0, wEm, wFm
    real(8) :: wk_v11, wk_v22, wk_v33, Ubondmorse
    real(8) :: sfact=1.0d0  ! scaling factor to potE and virial (default)
    real(8) :: wkvd(6), vd(6) ! debug
    real(8) :: v21,v31,v32
    real(8) :: wk_v21, wk_v31, wk_v32
    real(8) :: wvr(6),vr(6) ! debug
    real(8) :: fi0(3), oUb
    integer(4) :: iam
    include 'mpif.h'

    call check_atombound_bond_morse2

    iam = 0
    Ubondmorse = 0.0d0
    wk_v11=0d0;wk_v22=0d0;wk_v33=0d0
    wk_v21=0d0;wk_v31=0d0;wk_v32=0d0

!$omp parallel default(none) &
!$omp& private(k0,i0,j0,ipar) &
!$omp& private(iam,iA,iB,drAB__x,drAB__y,drAB__z) &
!$omp& private(b,b_b0,Fb__x,Fb__y,Fb__z) &
!$omp& private(v11,v22,v33,fi0) &
!$omp& shared(m2i,i2m,nA,atom2A,topA,wkxyz) &
!$omp& shared(w3_f,wk_vir2,paranum) &
!$omp& private(wD,wam,wbm,wb0,wEm,wFm) &
!$omp& shared(DA,amA,bmA,b0A) &
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
            call morse_function(DA(ipar, iB), amA(ipar, iB), bmA(ipar, iB), b0A(ipar, iB),b,wEm,wFm)
            Fb__x = - wFm*drAB__x / b
            Fb__y = - wFm*drAB__y / b
            Fb__z = - wFm*drAB__z / b


             fi0(1)=fi0(1)+Fb__x
             fi0(2)=fi0(2)+Fb__y
             fi0(3)=fi0(3)+Fb__z
!default: MTD
             w3_f(1,j0,iam) = w3_f(1,j0,iam) - Fb__x
             w3_f(2,j0,iam) = w3_f(2,j0,iam) - Fb__y
             w3_f(3,j0,iam) = w3_f(3,j0,iam) - Fb__z
            Ubondmorse = Ubondmorse + wEm
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

    wk_p_energy = wk_p_energy + Ubondmorse * sfact 
    wk_vir2(1,0) = wk_vir2(1,0) + wk_v11
    wk_vir2(2,0) = wk_vir2(2,0) + wk_v22
    wk_vir2(3,0) = wk_vir2(3,0) + wk_v33
    wk_vir2(4,0) = wk_vir2(4,0) + wk_v21
    wk_vir2(5,0) = wk_vir2(5,0) + wk_v31
    wk_vir2(6,0) = wk_vir2(6,0) + wk_v32
#ifdef DEBUGFCE
    oUb=0d0
    call mpi_allreduce(Ubondmorse*sfact,oUb,1, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,ipar)
    if(myrank==0) write(*,*) 'Pot(bond morse2) =', oUb *kJ_mol,'[kJ/mol]'
!debug molp
    wkvd(1)=wk_v11; wkvd(2)=wk_v22; wkvd(3)=wk_v33
    wkvd(4)=wk_v21; wkvd(5)=wk_v31; wkvd(6)=wk_v32
    call mpi_allreduce(wkvd,vd,6, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,ipar)
    if(myrank==0) write(*,*) 'vir(bond) =', vd(1:6), sum(vd(1:3))
#endif
  end subroutine add_bond_morse2_no_bond_breaking
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate energy and force of bond.
!! \author Kazushi Fujimoto, Zhiye Tang
!<
  subroutine add_bond_morse2_bond_breaking()
!----------------------------------------------------------------------
!
!      V = D*[1-exp{-a(b-b0)-b(b-b0)**2}]**2
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
    use md_condition, only : allow_bond_breaking
#ifdef ONEPROC_AXIS
    use boundary, only : pbc_pair
#endif
    implicit none
    integer(4) :: iA, iB, i0, j0, ipar, k0
    real(8) :: b, b_b0, Fb__x, Fb__y, Fb__z, drAB__x, drAB__y, drAB__z
    real(8) :: v11,v22,v33
    real(8) :: wD, wam, wbm, wb0, wEm, wFm
    real(8) :: wk_v11, wk_v22, wk_v33, Ubondmorse
    real(8) :: sfact=1.0d0  ! scaling factor to potE and virial (default)
    real(8) :: wkvd(6), vd(6) ! debug
    real(8) :: v21,v31,v32
    real(8) :: wk_v21, wk_v31, wk_v32
    real(8) :: wvr(6),vr(6) ! debug
    real(8) :: fi0(3), oUb
    integer(4) :: iam
    include 'mpif.h'

    call check_atombound_bond_morse2

    iam = 0
    Ubondmorse = 0.0d0
    wk_v11=0d0;wk_v22=0d0;wk_v33=0d0
    wk_v21=0d0;wk_v31=0d0;wk_v32=0d0

!$omp parallel default(none) &
!$omp& private(k0,i0,j0,ipar) &
!$omp& private(iam,iA,iB,drAB__x,drAB__y,drAB__z) &
!$omp& private(b,b_b0,Fb__x,Fb__y,Fb__z) &
!$omp& private(v11,v22,v33,fi0) &
!$omp& shared(m2i,i2m,nA,atom2A,topA,wkxyz) &
!$omp& shared(w3_f,wk_vir2,paranum) &
!$omp& private(wD,wam,wbm,wb0,wEm,wFm) &
!$omp& shared(DA,amA,bmA,b0A) &
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

            if(DA(ipar,iB)<0.0d0) then
              !write(*,*) "removed bond_morse2", iA, atom2A(ipar,iB)+(iA-ipar)
              cycle ! removed entry
            endif

            j0 = i2m(atom2A(ipar,iB)+(iA-ipar))
            drAB__x = wkxyz(1,i0) - wkxyz(1,j0)
            drAB__y = wkxyz(2,i0) - wkxyz(2,j0)
            drAB__z = wkxyz(3,i0) - wkxyz(3,j0)
#ifdef ONEPROC_AXIS
            call pbc_pair(drAB__x,drAB__y,drAB__z)
#endif
            b = sqrt(drAB__x*drAB__x+drAB__y*drAB__y+drAB__z*drAB__z)
            call morse_function(DA(ipar, iB), amA(ipar, iB), bmA(ipar, iB), b0A(ipar, iB),b,wEm,wFm)
            Fb__x = - wFm*drAB__x / b
            Fb__y = - wFm*drAB__y / b
            Fb__z = - wFm*drAB__z / b


             fi0(1)=fi0(1)+Fb__x
             fi0(2)=fi0(2)+Fb__y
             fi0(3)=fi0(3)+Fb__z
!default: MTD
             w3_f(1,j0,iam) = w3_f(1,j0,iam) - Fb__x
             w3_f(2,j0,iam) = w3_f(2,j0,iam) - Fb__y
             w3_f(3,j0,iam) = w3_f(3,j0,iam) - Fb__z
            Ubondmorse = Ubondmorse + wEm
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

    wk_p_energy = wk_p_energy + Ubondmorse * sfact 
    wk_vir2(1,0) = wk_vir2(1,0) + wk_v11
    wk_vir2(2,0) = wk_vir2(2,0) + wk_v22
    wk_vir2(3,0) = wk_vir2(3,0) + wk_v33
    wk_vir2(4,0) = wk_vir2(4,0) + wk_v21
    wk_vir2(5,0) = wk_vir2(5,0) + wk_v31
    wk_vir2(6,0) = wk_vir2(6,0) + wk_v32
#ifdef DEBUGFCE
    oUb=0d0
    call mpi_allreduce(Ubondmorse*sfact,oUb,1, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,ipar)
    if(myrank==0) write(*,*) 'Pot(bond morse2) =', oUb *kJ_mol,'[kJ/mol]'
!debug molp
    wkvd(1)=wk_v11; wkvd(2)=wk_v22; wkvd(3)=wk_v33
    wkvd(4)=wk_v21; wkvd(5)=wk_v31; wkvd(6)=wk_v32
    call mpi_allreduce(wkvd,vd,6, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,ipar)
    if(myrank==0) write(*,*) 'vir(bond) =', vd(1:6), sum(vd(1:3))
#endif
  end subroutine add_bond_morse2_bond_breaking
!-----------------------------------------------------------------------
    subroutine morse_function(D,am,bm,b0,b,Em,Fm)
    implicit none 
!
!      V = D*[1-expr{-am*(b-b0)-bm*(b-b0)**2}]**2

!      Fm = dV/db
!
!      A ----- B
!          b
!
      real(8) :: b_b0, Ex
      real(8),intent(in) :: D, am, bm, b0, b
      real(8),intent(out) :: Em, Fm

      b_b0 = b - b0
      Ex = exp(-am*b_b0-bm*b_b0**2)
      Em = D*(1.0d0-Ex)**2
      Fm = 2.0d0*D*Ex*(1.0d0-Ex)*(am+2.0d0*bm*b_b0)

    end  subroutine morse_function


end module bond_morse2

