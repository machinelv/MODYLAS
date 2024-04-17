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
!! \brief  Module and subroutines to calculate UB potential.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to parse parameters of UB potential.
!! \author Kensuke Iwahashi
!<
module parse_UB
  implicit none
  type yyparse_ub
     integer(4) :: atom3
     real(8) :: Kub, s0
     type(yyparse_ub),pointer :: next
  end type yyparse_ub

  type(yyparse_ub),pointer :: yyparse_ub_top(:)
  integer(4),allocatable :: yyparse_ub_n(:)
  integer(4) :: yyparse_ub_total=0

contains

  subroutine fmod_alloc_yyparse_ub_top(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    integer(4) :: i
    type(yyparse_ub),pointer :: top
    allocate(yyparse_ub_top(ivalue))
    do i=1,ivalue
       allocate(top)
       top%atom3 = -1
       nullify(top%next)
       yyparse_ub_top(i) = top
    enddo
  end subroutine fmod_alloc_yyparse_ub_top
!----------------------------------------------------------------------
  subroutine fmod_alloc_yyparse_ub_n(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(yyparse_ub_n(ivalue))
    yyparse_ub_n = 0
  end subroutine fmod_alloc_yyparse_ub_n
! -----------------------------------------------------------------------
!>
!! \brief  Subroutine to keep UB pair in liner list temporally.
!! \author Kensuke Iwahashi
!<
  subroutine add_to_ub_list(atom1,atom3,Kub,s0)
    use trajectory_org
    use param
    use mpi_tool
    implicit none
    integer(4), intent(in) :: atom1,atom3
    real(8), intent(in) :: Kub,s0
    type(yyparse_ub),pointer :: new1,p,next
    if (atom1 > npara) then
       write(0,'(a,i0,a)') &
     &    'ERROR: The number of ub is out of bounds.  '// &
     &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom3 > npara) then
       write(0,'(a,i0,a)') &
     &    'ERROR: The number of ub is out of bounds.  '// &
     &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    allocate(new1)
    nullify(new1%next)
    new1%atom3 = atom3
    new1%Kub = Kub
    new1%s0 = s0
    p => yyparse_ub_top(atom1)
    !     insert to list, sorting from smaller
    do while (.true.)
       next => p%next
       if (.not. associated(next)) then
          p%next => new1
          yyparse_ub_n(atom1) = yyparse_ub_n(atom1) + 1
          yyparse_ub_total = yyparse_ub_total + 1
          exit
       endif
       !       doubled must be merged
       if (atom3 == next%atom3) then
          next%Kub = new1%Kub
          next%s0 = new1%s0
          deallocate(new1)
          exit
       else if (atom3 < next%atom3) then
          new1%next => p%next
          p%next => new1
          yyparse_ub_n(atom1) = yyparse_ub_n(atom1) + 1
          yyparse_ub_total = yyparse_ub_total + 1
          exit
       endif
       p => next
    enddo
  end subroutine add_to_ub_list
!-----------------------------------------------------------------------
  subroutine add_to_mol_ub_list(ivalue1, atom1, atom2, atom3, Kub, s0)
    use mol_info, only : para_start
    implicit none
    integer(4), intent(in) :: ivalue1
    integer(4), intent(inout) :: atom1, atom2, atom3
    real(8), intent(in) :: Kub, s0

    atom1 = atom1 + para_start(ivalue1)
    atom2 = atom2 + para_start(ivalue1)
    atom3 = atom3 + para_start(ivalue1)

    call add_to_ub_list(atom1, atom3, Kub, s0)
    call add_to_ub_list(atom3, atom1, Kub, s0)
  end subroutine add_to_mol_ub_list

end module parse_UB
!----------------------------------------------------------------------
!>
!! \brief  Module to calculate UB potential.
!! \author Kensuke Iwahashi
!<
!----------------------------------------------------------------------
module UB
  use omp_lib
  use parse_UB
  use subcell
  implicit none
  integer(4) :: nubs=0
  integer(4) :: nub_as=0
  integer(4), allocatable :: atom1(:), atom2(:), atom3(:)
  real(8), allocatable :: Kub(:), s0(:)
  integer(4), allocatable :: atom3A(:,:)
  real(8), allocatable :: KubA(:,:), s0A(:,:)
  integer(4), allocatable :: topA(:), nA(:)

contains

  !subroutine init_UB
  !end subroutine init_UB

  !subroutine intra_UBs
  !end subroutine intra_UBs

!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to check whether j-atoms are communicated correctly
!! \author Yoshimichi Andoh
!<
  subroutine check_atombound_UB()
!-----------------------------------------------------------------------
    use trajectory_mpi
    use param
    use md_condition
    use mpi_tool
    implicit none
    integer(4) :: k0,i0,j,i0c
    integer(4) :: iA,ipar
    integer(4) :: icheck

    icheck=0

!$omp parallel default(none) &
!$omp& private(k0,i0,i0c,j,iA,ipar) &
!$omp& shared(nselfseg,lsegtop,lseg_natoms) &
!$omp& shared(m2i,paranum,nA,atom3A,i2m) &
!$omp& reduction(+:icheck)
!$omp do
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          iA   = m2i(i0)
          ipar = paranum(iA)
          do j = 1, nA(ipar)
             i0c = i2m(atom3A(ipar,j)+(iA-ipar))
             if(i0c.eq.-1) icheck = icheck -1
          enddo
       enddo
    enddo
!$omp end do
!$omp end parallel

    if(icheck.le.-1) then
       write(0,*) 'ERROR[add_charmm_ub_a] &
            &   There is a particle outside the area.'
       write(0,*) 'Myrank,mdstep=',myrank,mdstep
       call modylas_abort()
    end if
  end subroutine check_atombound_UB
!-----------------------------------------------------------------------
  subroutine set_md_charmm_a_ub
    use param
    implicit none
    integer(4) :: i,j,k,maxi(1)
    type(yyparse_ub),pointer :: p,freed

    nub_as = yyparse_ub_total
    if (yyparse_ub_total == 0)  return
    maxi = maxloc(yyparse_ub_n)
    maxi =  yyparse_ub_n(maxi(1))
    allocate(atom3A(npara,maxi(1)))
    allocate(KubA(npara,maxi(1)))
    allocate(s0A(npara,maxi(1)))
    allocate(nA(npara))
    !     copying special pair entries

    do i=1, npara
       nA(i) = yyparse_ub_n(i)
       p => yyparse_ub_top(i)
       do j=1,yyparse_ub_n(i)
          p => p%next
          atom3A(i,j) = p%atom3
          KubA(i,j) = p%Kub
          s0A(i,j) = p%s0
       enddo
    enddo

    deallocate(yyparse_ub_top, yyparse_ub_n)
  end subroutine set_md_charmm_a_ub

!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to read UB parameter.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_md_charmm_a_ub
    use trajectory_org
    use device_numbers, only : f_mdff
    use param
    use mpi_tool
    implicit none
    include 'mpif.h'
    integer(4) :: ierr
    integer(4) :: num

    if(myrank==0) then
       read(f_mdff) nub_as
    endif
    call MPI_Bcast(nub_as,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    if(nub_as.eq.0) return

    if(myrank==0) read(f_mdff) num
    call MPI_Bcast(num,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)

    allocate(nA(npara))
    allocate(atom3A(npara,num))
    allocate(KubA(npara,num))
    allocate(s0A(npara,num))

    if(myrank==0) then
       !      read(f_mdff) topA
       read(f_mdff) nA
       read(f_mdff) atom3A
       read(f_mdff) KubA
       read(f_mdff) s0A
    endif
    call MPI_Bcast(nA,npara,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom3A,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(KubA,npara*num,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(s0A,npara*num,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  end subroutine read_md_charmm_a_ub
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write UB parameter in mdff.bin File.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_md_charmm_a_ub
    use device_numbers, only : f_mdff
    implicit none
    integer(4) :: num
    write(f_mdff) nub_as
    if(nub_as.eq.0) return
    num = size(atom3A(:,:), 2)
    write(f_mdff) num
    write(f_mdff) nA
    write(f_mdff) atom3A
    write(f_mdff) KubA
    write(f_mdff) s0A
  end subroutine write_md_charmm_a_ub
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write UB parameter.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_md_charmm_a_ub
    implicit none
    integer(4) :: num
    write(*,*) '[write_md_charmm_a_ub]'
    write(*,*) nub_as
    if(nub_as.eq.0) return
    num = size(atom3A(:,:), 2)
    write(*,*) num
    write(*,*) nA
    write(*,*) atom3A
    write(*,*) KubA
    write(*,*) s0A
  end subroutine write_memory_md_charmm_a_ub
!-------------------------------------------------------------------------
  subroutine read_md_charmm_ub
    !fujimoto
    use device_numbers, only : f_mdff
    implicit none

    read(f_mdff) nubs
    if(nubs.eq.0) return
    allocate(atom1(nubs))
    allocate(atom2(nubs))
    allocate(atom3(nubs))
    allocate(Kub(nubs))
    allocate(s0(nubs))
    read(f_mdff) atom1
    read(f_mdff) atom2
    read(f_mdff) atom3
    read(f_mdff) Kub
    read(f_mdff) s0
  end subroutine read_md_charmm_ub
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate energy and force of Urey-Bradley.
!! \author Kensuke Iwahashi
!<
  subroutine add_charmm_ub_a()
!-----------------------------------------------------------------------
!
!      V = Kub(s - s0)^2
!
!           B
!          / \    s is length between A and B
!         /   \     !! there is a bond between A and B and
!        A=====!    !! between B and C but not between A and C
!           s
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
    integer(4) :: i, j, ip, ia, i0c, i0, ipar, k0
    real(8) :: s, s_s0, Uub, oUu
    real(8) :: drAC__x, drAC__y, drAC__z, Fc__x, Fc__y, Fc__z
    real(8) :: v11, v22, v33
    real(8) :: wk_v11,wk_v22,wk_v33
    real(8) :: sfact=1.0d0  ! scaling factor to potE and virial (default)
    real(8) :: wkvd(6), vd(6) ! debug
    real(8) :: v21,v31,v32
    real(8) :: wk_v21, wk_v31, wk_v32
    real(8) :: wvr(6),vr(6) ! debug
    real(8) :: fi0(3)
    integer(4) :: iam
    include 'mpif.h'

    call check_atombound_UB()

    iam = 0
    Uub = 0.0d0
    wk_v11=0d0;wk_v22=0d0;wk_v33=0d0
    wk_v21=0d0;wk_v31=0d0;wk_v32=0d0

!$omp parallel default(none) &
!$omp& private(k0,i0,ipar,i,iam) &
!$omp& private(ia,i0c,ip,j,drAC__x,drAC__y,drAC__z) &
!$omp& private(s,s_s0,Fc__x,Fc__y,Fc__z) &
!$omp& private(v11,v22,v33) &
!$omp& private(fi0) &
!$omp& shared(m2i,i2m,nA,atom3A,topA,wkxyz) &
!$omp& shared(s0A,KubA,w3_f,wk_vir2,paranum) &
!$omp& shared(tag,na_per_cell) &
!$omp& shared(lsegtop,lseg_natoms,nselfseg) &
!$omp& shared(sfact) &
!$omp& private(v21,v31,v32) &
!$omp& reduction(+:wk_v21,wk_v31,wk_v32) &
!$omp& reduction(+:Uub) &
!$omp& reduction(+:wk_v11,wk_v22,wk_v33)
!$    iam = omp_get_thread_num()
!$omp do
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          ia   = m2i(i0)
          ipar = paranum(ia)
          fi0(1:3)=0d0
          do j = 1, nA(ipar)
!default: MTD
          if(ia < atom3A(ipar,j)+(ia-ipar))cycle ! avoid double count
             i0c = i2m(atom3A(ipar,j)+(ia-ipar))
             drAC__x = wkxyz(1,i0 ) - wkxyz(1,i0c)
             drAC__y = wkxyz(2,i0 ) - wkxyz(2,i0c)
             drAC__z = wkxyz(3,i0 ) - wkxyz(3,i0c)
#ifdef ONEPROC_AXIS
             call pbc_pair(drAC__x,drAC__y,drAC__z)
#endif
             s = sqrt(drAC__x*drAC__x+drAC__y*drAC__y+drAC__z*drAC__z)
             s_s0 = s - s0A(ipar,j)
             Fc__x = -drAC__x * 2.0d0 * KubA(ipar,j) * s_s0 / s
             Fc__y = -drAC__y * 2.0d0 * KubA(ipar,j) * s_s0 / s
             Fc__z = -drAC__z * 2.0d0 * KubA(ipar,j) * s_s0 / s
             fi0(1)=fi0(1)+Fc__x
             fi0(2)=fi0(2)+Fc__y
             fi0(3)=fi0(3)+Fc__z
!default: MTD
             w3_f(1,i0c,iam) = w3_f(1,i0c,iam) - Fc__x
             w3_f(2,i0c,iam) = w3_f(2,i0c,iam) - Fc__y
             w3_f(3,i0c,iam) = w3_f(3,i0c,iam) - Fc__z
             Uub = Uub + KubA(ipar,j) * s_s0 * s_s0
             v11 = drAC__x * Fc__x * sfact 
             v22 = drAC__y * Fc__y * sfact
             v33 = drAC__z * Fc__z * sfact
             wk_v11=wk_v11+v11
             wk_v22=wk_v22+v22
             wk_v33=wk_v33+v33
             v21 = drAC__y * Fc__x * sfact
             v31 = drAC__z * Fc__x * sfact
             v32 = drAC__z * Fc__y * sfact
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

    wk_p_energy = wk_p_energy + Uub * sfact
    wk_vir2(1,0) = wk_vir2(1,0) + wk_v11
    wk_vir2(2,0) = wk_vir2(2,0) + wk_v22
    wk_vir2(3,0) = wk_vir2(3,0) + wk_v33
    wk_vir2(4,0) = wk_vir2(4,0) + wk_v21
    wk_vir2(5,0) = wk_vir2(5,0) + wk_v31
    wk_vir2(6,0) = wk_vir2(6,0) + wk_v32
#ifdef DEBUGFCE
    oUu=0d0
    call mpi_allreduce(Uub*sfact,oUu,1, &
         &     mpi_double_precision,mpi_sum,mpi_comm_world,ipar)
#ifdef KCAL
    if(myrank==0) write(*,*) 'Pot(UB   )=', oUu *kJ_mol/4.184d0,'[kcal/mol]'
#else
    if(myrank==0) write(*,*) 'Pot(UB   )=', oUu *kJ_mol,'[kJ/mol]'
#endif
#endif
  end subroutine add_charmm_ub_a

end module UB
