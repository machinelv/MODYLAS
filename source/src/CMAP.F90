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
!! \brief  Module and subroutines to CMAP potential.
!<
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to parse parameters of CMAP potential.
!! \author Kazushi FUJIMOTO
!<
module parse_CMAP
  use mpi_tool
  implicit none
  type yyparse_cmapA
     integer(4) :: atom2,atom3,atom4,atom5
     real(8) :: KOC
     type(yyparse_cmapA),pointer :: next
  end type yyparse_cmapA
  type(yyparse_cmapA),pointer :: yyparse_cmap_topA(:)
  integer(4),allocatable :: yyparse_cmap_nA(:)
  integer(4) :: yyparse_cmap_totalA=0
  type yyparse_cmapB
     integer(4) :: atom1,atom3, atom4,atom5
     real(8) :: KOC
     type(yyparse_cmapB),pointer :: next
  end type yyparse_cmapB
  type(yyparse_cmapB),pointer :: yyparse_cmap_topB(:)
  integer(4),allocatable :: yyparse_cmap_nB(:)
  integer(4) :: yyparse_cmap_totalB=0
  type yyparse_cmapC
     integer(4) :: atom1,atom2,atom4,atom5
     real(8) :: KOC
     type(yyparse_cmapC),pointer :: next
  end type yyparse_cmapC
  type(yyparse_cmapC),pointer :: yyparse_cmap_topC(:)
  integer(4),allocatable :: yyparse_cmap_nC(:)
  integer(4) :: yyparse_cmap_totalC=0
  type yyparse_cmapD
     integer(4) :: atom1,atom2,atom3,atom5
     real(8) :: KOC
     type(yyparse_cmapD),pointer :: next
  end type yyparse_cmapD
  type(yyparse_cmapD),pointer :: yyparse_cmap_topD(:)
  integer(4),allocatable :: yyparse_cmap_nD(:)
  integer(4) :: yyparse_cmap_totalD=0
  type yyparse_cmapE
     integer(4) :: atom1,atom2,atom3,atom4
     real(8) :: KOC
     type(yyparse_cmapE),pointer :: next
  end type yyparse_cmapE
  type(yyparse_cmapE),pointer :: yyparse_cmap_topE(:)
  integer(4),allocatable :: yyparse_cmap_nE(:)
  integer(4) :: yyparse_cmap_totalE=0

contains
!----------------------------------------------------------------------
!>
!! \brief  Subroutines to read inputs relating to CMAP potential.
!! \author Kazushi FUJIMOTO
!<
  subroutine add_to_cmap_listA(atom1,atom2,atom3,atom4,atom5,KOC)
    !fujimoto
    use trajectory_org
    use param
    implicit none
    integer(4), intent(in) :: atom1,atom2,atom3,atom4,atom5,KOC
    type(yyparse_cmapA),pointer :: new1,p,next
    if (atom1 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of cmap is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom2 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
             &    'The number of cmap is out of bounds.  '// &
             &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom3 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of cmap is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom4 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of cmap is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom5 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of cmap is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif

    allocate(new1)
    nullify(new1%next)
    new1%atom2 = atom2
    new1%atom3 = atom3
    new1%atom4 = atom4
    new1%atom5 = atom5
    new1%KOC = KOC
    p => yyparse_cmap_topA(atom1)
    !     insert to list, sorting from smaller
    do while (.true.)
       next => p%next
       if (.not. associated(next)) then
          p%next => new1
          yyparse_cmap_nA(atom1) = yyparse_cmap_nA(atom1) + 1
          yyparse_cmap_totalA = yyparse_cmap_totalA + 1
          exit
       endif
       if (atom2 <= next%atom2 .and. atom3 <= next%atom3 .and. &
            &      atom4 <= next%atom4 .and. atom5 <= next%atom5) then
          new1%next => p%next
          p%next => new1
          yyparse_cmap_nA(atom1) = yyparse_cmap_nA(atom1) + 1
          yyparse_cmap_totalA = yyparse_cmap_totalA + 1
          exit
       endif
       p => next
    enddo
  end subroutine add_to_cmap_listA
!-----------------------------------------------------------------------
!>
!! \brief  Subroutines to read inputs relating to CMAP potential.
!! \author Kazushi FUJIMOTO
!<
  subroutine add_to_cmap_listB(atom1,atom2,atom3,atom4,atom5,KOC)
    !fujimoto
    use trajectory_org
    use param
    implicit none
    integer(4), intent(in) :: atom1,atom2,atom3,atom4,atom5,KOC
    type(yyparse_cmapB),pointer :: new1,p,next
    if (atom1 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of cmap is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom2 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of cmap is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom3 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of cmap is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom4 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of cmap is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom5 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of cmap is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif

    allocate(new1)
    nullify(new1%next)
    new1%atom1 = atom1
    new1%atom3 = atom3
    new1%atom4 = atom4
    new1%atom5 = atom5
    new1%KOC = KOC
    p => yyparse_cmap_topB(atom2)
    !     insert to list, sorting from smaller
    do while (.true.)
       next => p%next
       if (.not. associated(next)) then
          p%next => new1
          yyparse_cmap_nB(atom2) = yyparse_cmap_nB(atom2) + 1
          yyparse_cmap_totalB = yyparse_cmap_totalB + 1
          exit
       endif
       if (atom1 <= next%atom1 .and. atom3 <= next%atom3 .and. &
            &      atom4 <= next%atom4 .and. atom5 <= next%atom5) then
          new1%next => p%next
          p%next => new1
          yyparse_cmap_nB(atom2) = yyparse_cmap_nB(atom2) + 1
          yyparse_cmap_totalB = yyparse_cmap_totalB + 1
          exit
       endif
       p => next
    enddo
  end subroutine add_to_cmap_listB
!-----------------------------------------------------------------------
!>
!! \brief  Subroutines to read inputs relating to CMAP potential.
!! \author Kazushi FUJIMOTO
!<
  subroutine add_to_cmap_listC(atom1,atom2,atom3,atom4,atom5,KOC)
    !fujimoto
    use trajectory_org
    use param
    implicit none
    integer(4), intent(in) :: atom1,atom2,atom3,atom4,atom5,KOC
    type(yyparse_cmapC),pointer :: new1,p,next
    if (atom1 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of cmap is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom2 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of cmap is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom3 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of cmap is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom4 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of cmap is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom5 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of cmap is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif

    allocate(new1)
    nullify(new1%next)
    new1%atom1 = atom1
    new1%atom2 = atom2
    new1%atom4 = atom4
    new1%atom5 = atom5
    new1%KOC = KOC
    p => yyparse_cmap_topC(atom3)
    !     insert to list, sorting from smaller
    do while (.true.)
       next => p%next
       if (.not. associated(next)) then
          p%next => new1
          yyparse_cmap_nC(atom3) = yyparse_cmap_nC(atom3) + 1
          yyparse_cmap_totalC = yyparse_cmap_totalC + 1
          exit
       endif
       if (atom1 <= next%atom1 .and. atom2 <= next%atom2 .and. &
            &      atom4 <= next%atom4 .and. atom5 <= next%atom5) then
          new1%next => p%next
          p%next => new1
          yyparse_cmap_nC(atom3) = yyparse_cmap_nC(atom3) + 1
          yyparse_cmap_totalC = yyparse_cmap_totalC + 1
          exit
       endif
       p => next
    enddo
  end subroutine add_to_cmap_listC
!-----------------------------------------------------------------------
!>
!! \brief  Subroutines to read inputs relating to CMAP potential.
!! \author Kazushi FUJIMOTO
!<
  subroutine add_to_cmap_listD(atom1,atom2,atom3,atom4,atom5,KOC)
    !fujimoto
    use trajectory_org
    use param
    implicit none
    integer(4), intent(in) :: atom1,atom2,atom3,atom4,atom5,KOC
    type(yyparse_cmapD),pointer :: new1,p,next
    if (atom1 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of cmap is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom2 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of cmap is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom3 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of cmap is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom4 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of cmap is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom5 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of cmap is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif

    allocate(new1)
    nullify(new1%next)
    new1%atom1 = atom1
    new1%atom2 = atom2
    new1%atom3 = atom3
    new1%atom5 = atom5
    new1%KOC = KOC
    p => yyparse_cmap_topD(atom4)
    !     insert to list, sorting from smaller
    do while (.true.)
       next => p%next
       if (.not. associated(next)) then
          p%next => new1
          yyparse_cmap_nD(atom4) = yyparse_cmap_nD(atom4) + 1
          yyparse_cmap_totalD = yyparse_cmap_totalD + 1
          exit
       endif
       if (atom1 <= next%atom1 .and. atom2 <= next%atom2 .and. &
            &      atom3 <= next%atom3 .and. atom5 <= next%atom5) then
          new1%next => p%next
          p%next => new1
          yyparse_cmap_nD(atom4) = yyparse_cmap_nD(atom4) + 1
          yyparse_cmap_totalD = yyparse_cmap_totalD + 1
          exit
       endif
       p => next
    enddo
  end subroutine add_to_cmap_listD
!-----------------------------------------------------------------------
!>
!! \brief  Subroutines to read inputs relating to CMAP potential.
!! \author Kazushi FUJIMOTO
!<
  subroutine add_to_cmap_listE(atom1,atom2,atom3,atom4,atom5,KOC)
    !fujimoto
    use trajectory_org
    use param
    implicit none
    integer(4), intent(in) :: atom1,atom2,atom3,atom4,atom5,KOC
    type(yyparse_cmapE),pointer :: new1,p,next
    if (atom1 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of cmap is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom2 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of cmap is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom3 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of cmap is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom4 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of cmap is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom5 > npara) then
       write(0,'(a,i0,a)')  'ERROR: '// &
            &    'The number of cmap is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif

    allocate(new1)
    nullify(new1%next)
    new1%atom1 = atom1
    new1%atom2 = atom2
    new1%atom3 = atom3
    new1%atom4 = atom4
    new1%KOC = KOC
    p => yyparse_cmap_topE(atom5)
    !     insert to list, sorting from smaller
    do while (.true.)
       next => p%next
       if (.not. associated(next)) then
          p%next => new1
          yyparse_cmap_nE(atom5) = yyparse_cmap_nE(atom5) + 1
          yyparse_cmap_totalE = yyparse_cmap_totalE + 1
          exit
       endif
       if (atom1 <= next%atom1 .and. atom2 <= next%atom2 .and. &
            &      atom3 <= next%atom3 .and. atom4 <= next%atom4) then
          new1%next => p%next
          p%next => new1
          yyparse_cmap_nE(atom5) = yyparse_cmap_nE(atom5) + 1
          yyparse_cmap_totalE = yyparse_cmap_totalE + 1
          exit
       endif
       p => next
    enddo
  end subroutine add_to_cmap_listE

end module parse_CMAP
!----------------------------------------------------------------------
!>
!! \brief  Module to calculate CMAP potential.
!! \author Kazushi FUJIMOTO
!<
!----------------------------------------------------------------------
module CMAP
  use omp_lib
  use parse_CMAP
  use subcell
  use mpi_tool
  implicit none
  !map(i,j): charmmのCMAP。
  !          phi, psiともに、-180から180°まで,15°置きの25点のmap。
  !          i:psi j:phi
  !          単位は[J]. CMAPのオリジナルは[kcal/mol]
  !dmap(3,25,25): map(25,25)の偏微分値
  !      1: psiの1階の偏微分
  !      2: phiの1階の偏微分
  !      3: psi, phiの2階の偏微分
  !注: CMAPの計算にはラジアンを用いない。
  real(8) :: map(25,25,6),dmap(3,25,25,6)
  real(8),parameter :: PI=3.1415926535897932384626433832795029d0
  real(8),parameter :: r_PI=3.1830988618379D-01
  real(8), parameter :: rad2deg= 5.7295779513082D+01
  integer(4) :: CMAPVersion=36

  ! variables
  integer(4),allocatable :: atom2A(:,:),atom3A(:,:),&
       &                    atom4A(:,:),atom5A(:,:)
  integer(4),allocatable :: nA(:), KindOfCmapA(:,:)

  integer(4),allocatable :: atom1B(:,:),atom3B(:,:),&
       &                          atom4B(:,:),atom5B(:,:)
  integer(4),allocatable :: nB(:), KindOfCmapB(:,:)

  integer(4),allocatable :: atom1C(:,:),atom2C(:,:),&
       &                          atom4C(:,:),atom5C(:,:)
  integer(4),allocatable :: nC(:), KindOfCmapC(:,:)

  integer(4),allocatable :: atom1D(:,:),atom2D(:,:),&
       &                          atom3D(:,:),atom5D(:,:)
  integer(4),allocatable :: nD(:), KindOfCmapD(:,:)

  integer(4),allocatable :: atom1E(:,:),atom2E(:,:),&
       &                          atom3E(:,:),atom4E(:,:)
  integer(4),allocatable :: nE(:), KindOfCmapE(:,:)

  integer(4) :: ncmap=0  ! default must be 0

contains

!>
!! \brief  Subroutine to check whether j-,k-,l-atoms are communicated correctly
!! \author Yoshimichi Andoh
!<
  subroutine check_atombound_CMAP()
!-----------------------------------------------------------------------
    use trajectory_mpi
    use param
    use md_condition
    implicit none
    integer(4) :: k0,i0,j,i0a,i0b,i0c,i0d,i0e
    integer(4) :: iA,ipar
    integer(4) :: icheck

    icheck=0

!$omp parallel default(none) &
!$omp& private(k0,i0,i0a,i0b,i0c,i0d,i0e,j,iA,ipar) &
!$omp& shared(nselfseg,lsegtop,lseg_natoms,m2i,i2m) &
!$omp& shared(paranum,nA,nB,nC,nD,nE) &
!$omp& shared(atom1B,atom3B,atom4B,atom5B) &
!$omp& shared(atom1C,atom2C,atom4C,atom5C) &
!$omp& shared(atom1D,atom2D,atom3D,atom5D) &
!$omp& shared(atom1E,atom2E,atom3E,atom4E) &
!$omp& reduction(+:icheck)
!$omp do
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          iA   = m2i(i0)
          ipar = paranum(iA)
          ! C group: CA
          do j = 1, nC(ipar)
             i0b = i2m(atom1C(ipar,j)+(iA-ipar)) !C1
             i0c = i2m(atom2C(ipar,j)+(iA-ipar)) !N1
             i0d = i2m(atom4C(ipar,j)+(iA-ipar)) !C2
             i0e = i2m(atom5C(ipar,j)+(iA-ipar)) !N2
             if(i0b.eq.-1) icheck = icheck -1
             if(i0c.eq.-1) icheck = icheck -1
             if(i0d.eq.-1) icheck = icheck -1
             if(i0e.eq.-1) icheck = icheck -1
          enddo
       enddo
    enddo
!$omp end do
!$omp end parallel

    if(icheck.le.-1) then
       write(*,*) 'Err[add_charmm_CMAP_a] &
            &   There is a particle outside the area.'
       write(*,*) 'Myrank,mdstep=',myrank,mdstep
       call modylas_abort()
    end if
  end subroutine check_atombound_CMAP
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate dihedral angles.
!! \author Kazushi FUJIMOTO
!<
  subroutine cal_dihedral(r, theta)
    use md_const
!!!   use md_const
    !This program written by K.FUJIMOTO.
    !r is the coordinates of 4 atoms.
    !r is need for calculate dihedral.
    !alpah, beta, gamm are cros producte of r1 and r2, r2 and r3,
    !r3 and r4, respectively.
    !theta is dihedral and that unit is degree.
    !
    !
    !     1__                __ 4
    !     |\                  /|
    ! alpha \                / gamma
    !        \              /
    !         ============>
    !        2   beta      3

    implicit none
    !INPUT DATA
    !r(4,3) : the position of four atoms
    real(8), intent(in) :: r(3,4)
    !OUTPUT DATA
    !theta: dihedral angles, -180 <= theta <= 180
    real(8), intent(out) :: theta
    real(8) :: alpha(3), beta(3), gamm(3), n1(3), n2(3), cpn1n2(3)
    real(8) :: inner, scalar1, scalar2,cos_theta
    real(8) :: inner_beta_cpn1n2
    integer(4) :: i
    !calculate three vectors, alpha, beta, gamm
    do i = 1, 3
       alpha(i) = r(i,2) - r(i,1)
       beta(i)  = r(i,3) - r(i,2)
       gamm(i)  = r(i,4) - r(i,3)
    enddo

    !calculation of two cros products, n1, n2
    n1(1) = alpha(2)*beta(3) - alpha(3)*beta(2)
    n1(2) = alpha(3)*beta(1) - alpha(1)*beta(3)
    n1(3) = alpha(1)*beta(2) - alpha(2)*beta(1)

    n2(1) = beta(2)*gamm(3) - beta(3)*gamm(2)
    n2(2) = beta(3)*gamm(1) - beta(1)*gamm(3)
    n2(3) = beta(1)*gamm(2) - beta(2)*gamm(1)

    !calculation of scalars of n1 and n2
    scalar1 = dsqrt(n1(1)*n1(1) + n1(2)*n1(2) + n1(3)*n1(3))
    scalar2 = dsqrt(n2(1)*n2(1) + n2(2)*n2(2) + n2(3)*n2(3))
    !calculation of inner product of n1 and n2
    inner = n1(1)*n2(1) + n1(2)*n2(2) + n1(3)*n2(3)
    cos_theta = inner / (scalar1*scalar2)

    !calculation of cross product of n1 and n2
    cpn1n2(1) = n1(2)*n2(3) - n1(3)*n2(2)
    cpn1n2(2) = n1(3)*n2(1) - n1(1)*n2(3)
    cpn1n2(3) = n1(1)*n2(2) - n1(2)*n2(1)

    inner_beta_cpn1n2 = beta(1)*cpn1n2(1) + beta(2)*cpn1n2(2) + beta(3)*cpn1n2(3)

    if(inner_beta_cpn1n2 .ge. 0) then

       if (cos_theta > 1.0d0) then
          theta = 0.0d0
       elseif(cos_theta < -1.0d0) then
          theta = PI
       else
          theta = dacos(cos_theta)
       endif

       theta =  theta*180.0d0*r_PI
    else

       if (cos_theta > 1.0d0) then
          theta = 0.0d0
       elseif(cos_theta < -1.0d0) then
          theta = PI
       else
          theta = dacos(cos_theta)
       endif

       theta = -theta*180.0d0*r_PI
    endif
    return
  end subroutine cal_dihedral
!####################################################################
!>
!! \brief  Subroutine to calculate  coefficient of force
!! \brief  [Ref.] A. Blondel et al. J. Comput. Chem. 17, 1132-1141 (1996)
!! \author Kazushi FUJIMOTO
!<
  subroutine ForceCoeff(r,coeff)
    implicit none
    !INPUT
    real(8) :: r(3,4)
    !OUTPUT
    real(8) :: coeff(3,4)
    real(8) :: F(3), G(3), H(3)
    real(8) :: A(3), B(3) !A=FxG, B=HxG
    real(8) :: rA2, rB2, sqrtG2, rsqrtG2, FG, HG
    integer(4) :: i

    do i = 1, 3
       F(i) = r(i,1) - r(i,2)
       G(i) = r(i,2) - r(i,3)
       H(i) = r(i,4) - r(i,3)
    enddo
    A(1) = F(2)*G(3) - F(3)*G(2)
    A(2) = F(3)*G(1) - F(1)*G(3)
    A(3) = F(1)*G(2) - F(2)*G(1)
    B(1) = H(2)*G(3) - H(3)*G(2)
    B(2) = H(3)*G(1) - H(1)*G(3)
    B(3) = H(1)*G(2) - H(2)*G(1)

    rA2 = A(1)*A(1) + A(2)*A(2) + A(3)*A(3)
    rA2 = 1.0d0 / rA2
    rB2 = B(1)*B(1) + B(2)*B(2) + B(3)*B(3)
    rB2 = 1.0d0 / rB2
    sqrtG2 = dsqrt(G(1)*G(1) + G(2)*G(2) + G(3)*G(3))
    rsqrtG2 = 1.0d0 / sqrtG2
    FG = F(1)*G(1) + F(2)*G(2) + F(3)*G(3)
    HG = H(1)*G(1) + H(2)*G(2) + H(3)*G(3)

    !d /dr(i)
    coeff(1,1) = -sqrtG2*rA2*A(1)
    coeff(2,1) = -sqrtG2*rA2*A(2)
    coeff(3,1) = -sqrtG2*rA2*A(3)
    !d /dr(j)
    coeff(1,2) = (sqrtG2*rA2+FG*rA2*rsqrtG2)*A(1) - HG*rB2*rsqrtG2*B(1)
    coeff(2,2) = (sqrtG2*rA2+FG*rA2*rsqrtG2)*A(2) - HG*rB2*rsqrtG2*B(2)
    coeff(3,2) = (sqrtG2*rA2+FG*rA2*rsqrtG2)*A(3) - HG*rB2*rsqrtG2*B(3)
    !d /dr(k)
    coeff(1,3) = (HG*rB2*rsqrtG2-sqrtG2*rB2)*B(1) - FG*rA2*rsqrtG2*A(1)
    coeff(2,3) = (HG*rB2*rsqrtG2-sqrtG2*rB2)*B(2) - FG*rA2*rsqrtG2*A(2)
    coeff(3,3) = (HG*rB2*rsqrtG2-sqrtG2*rB2)*B(3) - FG*rA2*rsqrtG2*A(3)
    !d /dr(l)
    coeff(1,4) = sqrtG2*rB2*B(1)
    coeff(2,4) = sqrtG2*rB2*B(2)
    coeff(3,4) = sqrtG2*rB2*B(3)
  end subroutine ForceCoeff

!A. Bondel and M. Karplus, J. Comp. Chem., 17, 1132(1996)
!####################################################################
!>
!! \brief  Subroutine to calculated force and potential energy of CMAP
!! \author Kazushi FUJIMOTO
!<
  subroutine MainCMAP(r,kmap,F,Vcmap)
    ! Calculating potential energy and force of cmap
    implicit none
    !phi: C-N-CA-C      psi: N-CA-C-N
    !INPUT
    real(8) :: r(3,5)
    integer(4) :: kmap
    !OUTPUT
    real(8) :: Vcmap, F(3,5)
    real(8) :: phi, psi
    real(8) :: phil, phiu, psil, psiu
    real(8) :: Fpsi, Fphi
    real(8) :: y(4), y1(4), y2(4), y12(4)
    real(8) :: coeff(3,4)
    integer(4) :: tnphi, tnpsi
    integer(4) :: i, j

    do i =1, 5
       do j = 1, 3
          F(j,i) = 0.0d0
       enddo
    enddo
    call cal_dihedral(r(:,1:4), phi)
    call cal_dihedral(r(:,2:5), psi)

    !      if(phi.gt.165.0d0) then
    if(phi.gt.180.0d0) then
       phi = phi -360.0d0
    elseif(phi.le.-180.0d0) then
       phi = phi + 360.0d0
    endif
    !      if(psi.gt.165.0d0) then
    if(psi.gt.180.0d0) then
       psi = psi -360.0d0
    elseif(psi.le.-180.0d0) then
       psi = psi + 360.0d0
    endif
    tnphi = dint((phi+180.0d0)/15.0d0)
    tnpsi = dint((psi+180.0d0)/15.0d0)
    phil = -180.0d0 + 15.0d0 * tnphi
    phiu = phil + 15.0d0
    psil = -180.0d0 + 15.0d0 * tnpsi
    psiu = psil + 15.0d0
    tnphi = tnphi + 1
    tnpsi = tnpsi + 1
    !DEBUG WRITE
    !     write(*,*) tnphi, tnpsi
    !     write(*,*) "phi, psi", phi, psi
    !     write(*,*) phil, phiu
    !     write(*,*) psil, psiu

    y(1) = map(tnpsi,  tnphi  ,kmap)
    y(2) = map(tnpsi,  tnphi+1,kmap)
    y(3) = map(tnpsi+1,tnphi+1,kmap)
    y(4) = map(tnpsi+1,tnphi  ,kmap)
    y1(1) =  dmap(1,tnpsi,  tnphi  ,kmap)
    y1(2) =  dmap(1,tnpsi,  tnphi+1,kmap)
    y1(3) =  dmap(1,tnpsi+1,tnphi+1,kmap)
    y1(4) =  dmap(1,tnpsi+1,tnphi  ,kmap)
    y2(1) =  dmap(2,tnpsi,  tnphi  ,kmap)
    y2(2) =  dmap(2,tnpsi,  tnphi+1,kmap)
    y2(3) =  dmap(2,tnpsi+1,tnphi+1,kmap)
    y2(4) =  dmap(2,tnpsi+1,tnphi  ,kmap)
    y12(1) = dmap(3,tnpsi,  tnphi  ,kmap)
    y12(2) = dmap(3,tnpsi,  tnphi+1,kmap)
    y12(3) = dmap(3,tnpsi+1,tnphi+1,kmap)
    y12(4) = dmap(3,tnpsi+1,tnphi  ,kmap)
    !DEBUG WRITE
    !     write(*,'(4es20.10)') y(1), y(2), y(3), y(4)
    !     write(*,'(4es20.10)') y1(1), y1(2), y1(3), y1(4)
    !     write(*,'(4es20.10)') y2(1), y2(2), y2(3), y2(4)
    !     write(*,'(4es20.10)') y12(1), y12(2), y12(3), y12(4)

    call BicubicInterpolation(y,y1,y2,y12,phil,psil,phi,psi, Vcmap, Fphi, Fpsi)

    call ForceCoeff(r(:,1:4), coeff)
    !DEBUG WRITE
    !     write(*,'(3es23.12)') Vcmap, Fphi, Fpsi
    !     write(*,*) (coeff(i,1)*Fphi,i=1,3)
    !     write(*,*) (coeff(i,2)*Fphi,i=1,3)
    !     write(*,*) (coeff(i,3)*Fphi,i=1,3)
    !     write(*,*) (coeff(i,4)*Fphi,i=1,3)
    do i = 1, 4
       F(1,i) = F(1,i) - coeff(1,i)*Fphi
       F(2,i) = F(2,i) - coeff(2,i)*Fphi
       F(3,i) = F(3,i) - coeff(3,i)*Fphi
    enddo
    call ForceCoeff(r(:,2:5), coeff)
    !DEBUG WRITE
    !     write(*,*)
    !     write(*,*) (coeff(i,1)*Fpsi,i=1,3)
    !     write(*,*) (coeff(i,2)*Fpsi,i=1,3)
    !     write(*,*) (coeff(i,3)*Fpsi,i=1,3)
    !     write(*,*) (coeff(i,4)*Fpsi,i=1,3)
    !     stop
    do i = 2, 5
       F(1,i) = F(1,i) - coeff(1,i-1)*Fpsi
       F(2,i) = F(2,i) - coeff(2,i-1)*Fpsi
       F(3,i) = F(3,i) - coeff(3,i-1)*Fpsi
    enddo
  END subroutine MainCMAP
!#####################################################################
!>
!! \brief  Subroutine to prepare for calculating force and potential energy of CMAP
!! \author Kazushi FUJIMOTO
!<
  subroutine PreCMAP()
    !Calculating partial differential at the point of map.
    !25x25 Map is extended to 49x49 as follow.
    !
    !
    !   |***********||---------|*********||----------|
    ! -360         -180        0        180         360
    !                         phi
    !
    !|--------|: CMAP of -180 <= phi <  0
    !|********|: CMAP of    0 <  phi <= 180
    !psi is treated likewise.
    implicit none
    !TempMap(49,49): Extended 49x49 map。
    !d2TempMap(49,49) :: Second derived function of TempMap
    real(8) :: TempMap(49,49), d2TempMap(49,49)
    real(8) :: phi(49), psi(49)
    real(8) :: d1, d2, d12
    real(8) :: tmp(49), d1tmp(49), d2tmp(49)
    real(8) :: ytmp(49), y2tmp(49)
    integer(4) :: i, ii, j, i2, j2, k

    ! log
    !      if(myrank==0) then
    !        write(*,*) " You selected CMAP Version", CMAPVersion
    !      endif

    if(CMAPVersion==22) then
       call set_CMAP_22()
    else if(CMAPVersion==36) then
       call set_CMAP_36()
    else
       write(0,*) "ERROR: No such CMAP version number"
       call modylas_abort()
    endif

    !Unit of angle is degree instead of radian.
    do i = 1, 49
       phi(i) = -360.0d0 + 15.0d0 *(i-1)
       psi(i) = -360.0d0 + 15.0d0 *(i-1)
    enddo

    do ii = 1, 6
       do i = 1, 49
          if(i.le.12) then
             i2 = i + 12
          elseif(13.le.i .and. i.le.37) then
             i2 = i - 12
          elseif(38.le.i) then
             i2 = i - 36
          endif
          do j = 1, 49
             if(j.le.12) then
                j2 = j + 12
             elseif(13.le.j .and. j.le.37) then
                j2 = j - 12
             elseif(38.le.j) then
                j2 = j - 36
             endif
             TempMap(i,j) = map(i2,j2,ii)
          enddo
       enddo

       do i = 1, 49
          do j = 1, 49
             ytmp(j) = TempMap(j,i)
          enddo
          call SecondDifferentialBySpline(ytmp, 49, y2tmp)
          do j = 1, 49
             d2TempMap(j,i) = y2tmp(j)
          enddo
       enddo

       do i = 13, 37
          do j = 13, 37
             do k = 1, 49
                do i2 = 1, 49
                   ytmp(i2) = TempMap(i2,k)
                   y2tmp(i2) = d2TempMap(i2,k)
                enddo
                call SplineInterpolationAtLatticePoint(49, ytmp, y2tmp, j, tmp(k), d1tmp(k))
             enddo
             call SecondDifferentialBySpline(tmp, 49, d2tmp)
             call SplineInterpolationAtLatticePoint(49, tmp, d2tmp, i, d2, d1)
             call SecondDifferentialBySpline(d1tmp, 49, d2tmp)
             call SplineInterpolationAtLatticePoint(49, d1tmp, d2tmp, i, d2, d12)
             dmap(1,j-12,i-12,ii) = d1
             dmap(2,j-12,i-12,ii) = d2
             dmap(3,j-12,i-12,ii) = d12
          enddo
       enddo
    enddo !do loop of ii
  end subroutine PreCMAP
!####################################################################
!>
!! \brief  Subroutine to calculate second differential
!! \author Kazushi FUJIMOTO
!<
  subroutine SecondDifferentialBySpline(f, np, f2)
    implicit none
    integer(4), intent(in) :: np
    real(8), intent(in) :: f(np)
    real(8), intent(out) :: f2(np)
    real(8) :: h = 15.0d0
    real(8) :: delta(np), a(np), b(np), c(np)
    integer(4) :: i

    do i = 1, np-1
       delta(i) = (f(i+1)-f(i))/h
    enddo

    a(2) = 4.0d0*h
    b(2) = delta(2) - delta(1)
    do i = 3, np-1
       a(i) = 4.0d0*h - (h*h)/a(i-1)
       b(i) = delta(i) - delta(i-1) - b(i-1)*h/a(i-1)
    enddo

    c(1)  = 0.0d0
    c(np) = 0.0d0

    c(np-1) = b(np-1)/a(np-1)
    do i = np-2, 2, -1
       c(i) = (b(i) - h*c(i+1))/a(i)
    enddo

    do i = 1, np
       f2(i) = 6.0d0*c(i)
    enddo

  end subroutine SecondDifferentialBySpline
!####################################################################
!>
!! \brief  Subroutine to calculate spline interpolation at lattice point
!! \author Kazushi FUJIMOTO
!<
  subroutine SplineInterpolationAtLatticePoint(np, f, f2, pj, p, p1)
    implicit none
    integer(4), intent(in) :: np, pj
    real(8), intent(in) :: f(np)
    real(8), intent(in) :: f2(np)
    !       real(8), intent(in) :: x
    real(8), intent(out) :: p, p1
    real(8) :: h = 15.0d0
    real(8) :: delta(np), sigma(np)
    real(8) :: b(np), c(np), d(np)
    integer(4) :: i

    do i = 1, np-1
       delta(i) = (f(i+1)-f(i))/h
    enddo
    do i = 1, np
       sigma(i) = f2(i) / 6.0d0
    enddo

    do i = 1, np-1
       b(i) = delta(i) -h*(sigma(i+1)+2.0d0*sigma(i))
       c(i) = 3.0d0*sigma(i)
       d(i) = (sigma(i+1)-sigma(i))/h
    enddo

    p  = f(pj)
    p1 = b(pj)
  end subroutine SplineInterpolationAtLatticePoint
!####################################################################
!>
!! \brief  Subroutine to calculate bicubic interpolation at lattice point
!! \author Kazushi FUJIMOTO
!<
  subroutine BicubicInterpolation(y,y1,y2,y12,phil,psil,phi,psi, V, Fphi, Fpsi)
    implicit none
    !INPUT
    !y, y1, y2, y12
    !OUTPUT
    !V, Fphi, Fpsi

    real(8) :: phil, psil, phi, psi
    real(8) :: V, Fphi, Fpsi
    real(8) :: y(4), y1(4), y2(4), y12(4)
    integer(4) :: i, j
    real(8) :: t, u, c(4,4)

    y1(:) = y1(:)*15.0d0
    y2(:) = y2(:)*15.0d0
    y12(:) = y12(:)*15.0d0*15.0d0

    c(1,1) =  y(1)
    c(1,2) = y2(1)
    c(1,3) = -3.0d0*y(1)+3.0d0*y(4)-2.0d0*y2(1)-y2(4)
    c(1,4) =  2.0d0*y(1)-2.0d0*y(4)+      y2(1)+y2(4)

    c(2,1) =  y1(1)
    c(2,2) = y12(1)
    c(2,3) = -3.0d0*y1(1)+3.0d0*y1(4)-2.0d0*y12(1)-y12(4)
    c(2,4) =  2.0d0*y1(1)-2.0d0*y1(4)+      y12(1)+y12(4)

    c(3,1) = -3.0d0*y(1)+3.0d0*y(2)-2.0d0*y1(1)-y1(2)
    c(3,2) = -3.0d0*y2(1)+3.0d0*y2(2)-2.0d0*y12(1)-y12(2)
    c(3,3) =  9.0d0*  y(1)-9.0d0*  y(2)+9.0d0*  y(3)-9.0d0*  y(4) &
         &          +6.0d0* y1(1)+3.0d0* y1(2)-3.0d0* y1(3)-6.0d0* y1(4) &
         &          +6.0d0* y2(1)-6.0d0* y2(2)-3.0d0* y2(3)+3.0d0* y2(4) &
         &          +4.0d0*y12(1)+2.0d0*y12(2)+      y12(3)+2.0d0*y12(4)
    c(3,4) = -6.0d0*  y(1)+6.0d0*  y(2)-6.0d0*  y(3)+6.0d0*  y(4) &
         &          -4.0d0* y1(1)-2.0d0* y1(2)+2.0d0* y1(3)+4.0d0* y1(4) &
         &          -3.0d0* y2(1)+3.0d0* y2(2)+3.0d0* y2(3)-3.0d0* y2(4) &
         &          -2.0d0*y12(1)-      y12(2)-      y12(3)-2.0d0*y12(4)

    c(4,1) =  2.0d0* y(1)-2.0d0* y(2)+ y1(1)+ y1(2)
    c(4,2) =  2.0d0*y2(1)-2.0d0*y2(2)+y12(1)+y12(2)
    c(4,3) = -6.0d0*  y(1)+6.0d0*  y(2)-6.0d0*  y(3)+6.0d0*  y(4) &
         &          -3.0d0* y1(1)-3.0d0* y1(2)+3.0d0* y1(3)+3.0d0* y1(4) &
         &          -4.0d0* y2(1)+4.0d0* y2(2)+2.0d0* y2(3)-2.0d0* y2(4) &
         &          -2.0d0*y12(1)-2.0d0*y12(2)-      y12(3)-      y12(4)
    c(4,4) =  4.0d0*  y(1)-4.0d0*  y(2)+4.0d0*  y(3)-4.0d0*  y(4) &
         &          +2.0d0* y1(1)+2.0d0* y1(2)-2.0d0* y1(3)-2.0d0* y1(4) &
         &          +2.0d0* y2(1)-2.0d0* y2(2)-2.0d0* y2(3)+2.0d0* y2(4) &
         &          +      y12(1)+      y12(2)+      y12(3)+      y12(4)

    t=(phi-phil)/15.0d0
    u=(psi-psil)/15.0d0
    V =0.0d0
    Fphi=0.0d0
    Fpsi=0.0d0
    do i = 1, 4
       do j = 1, 4
          V = V + c(i,j)*t**(i-1)*u**(j-1)
    if(i/=1)then  !! avoid (1/t) with t=0  debug@20190722
          Fphi = Fphi + (i-1)*c(i,j)*t**(i-2)*u**(j-1)
    endif
    if(j/=1)then  !! avoid (1/u) with u=0  debug@20190722
          Fpsi = Fpsi + (j-1)*c(i,j)*t**(i-1)*u**(j-2)
    endif
       enddo
    enddo

    Fphi = Fphi/15.0d0*rad2deg
    Fpsi = Fpsi/15.0d0*rad2deg
  end subroutine BicubicInterpolation
!====================================================================

  subroutine fmod_set_CMAPVersion(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    CMAPVersion=ivalue
  end subroutine fmod_set_CMAPVersion

  subroutine fmod_alloc_yyparse_CMAP_topA(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    integer(4) :: i
    type(yyparse_cmapA),pointer :: top
    allocate(yyparse_cmap_topA(ivalue))
    do i=1,ivalue
       allocate(top)
       top%atom2 = -1
       top%atom3 = -1
       top%atom4 = -1
       top%atom5 = -1
       nullify(top%next)
       yyparse_cmap_topA(i) = top
    enddo
  end subroutine fmod_alloc_yyparse_CMAP_topA
!----------------------------------------------------------------------
  subroutine fmod_alloc_yyparse_cmap_nA(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(yyparse_cmap_nA(ivalue))
    yyparse_cmap_nA = 0
  end subroutine fmod_alloc_yyparse_cmap_nA
!----------------------------------------------------------------------
  subroutine fmod_alloc_yyparse_CMAP_topB(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    integer(4) :: i
    type(yyparse_cmapB),pointer :: top
    allocate(yyparse_cmap_topB(ivalue))
    do i=1,ivalue
       allocate(top)
       top%atom1 = -1
       top%atom3 = -1
       top%atom4 = -1
       top%atom5 = -1
       nullify(top%next)
       yyparse_cmap_topB(i) = top
    enddo
  end subroutine fmod_alloc_yyparse_CMAP_topB
!----------------------------------------------------------------------
  subroutine fmod_alloc_yyparse_cmap_nB(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(yyparse_cmap_nB(ivalue))
    yyparse_cmap_nB = 0
  end subroutine fmod_alloc_yyparse_cmap_nB
!----------------------------------------------------------------------
  subroutine fmod_alloc_yyparse_CMAP_topC(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    integer(4) :: i
    type(yyparse_cmapC),pointer :: top
    allocate(yyparse_cmap_topC(ivalue))
    do i=1,ivalue
       allocate(top)
       top%atom1 = -1
       top%atom2 = -1
       top%atom4 = -1
       top%atom5 = -1
       nullify(top%next)
       yyparse_cmap_topC(i) = top
    enddo
  end subroutine fmod_alloc_yyparse_CMAP_topC
!----------------------------------------------------------------------
  subroutine fmod_alloc_yyparse_cmap_nC(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(yyparse_cmap_nC(ivalue))
    yyparse_cmap_nC = 0
  end subroutine fmod_alloc_yyparse_cmap_nC
!----------------------------------------------------------------------
  subroutine fmod_alloc_yyparse_CMAP_topD(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    integer(4) :: i
    type(yyparse_cmapD),pointer :: top
    allocate(yyparse_cmap_topD(ivalue))
    do i=1,ivalue
       allocate(top)
       top%atom1 = -1
       top%atom2 = -1
       top%atom3 = -1
       top%atom5 = -1
       nullify(top%next)
       yyparse_cmap_topD(i) = top
    enddo
  end subroutine fmod_alloc_yyparse_CMAP_topD
!----------------------------------------------------------------------
  subroutine fmod_alloc_yyparse_cmap_nD(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(yyparse_cmap_nD(ivalue))
    yyparse_cmap_nD = 0
  end subroutine fmod_alloc_yyparse_cmap_nD
!----------------------------------------------------------------------
  subroutine fmod_alloc_yyparse_CMAP_topE(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    integer(4) :: i
    type(yyparse_cmapE),pointer :: top
    allocate(yyparse_cmap_topE(ivalue))
    do i=1,ivalue
       allocate(top)
       top%atom1 = -1
       top%atom2 = -1
       top%atom3 = -1
       top%atom4 = -1
       nullify(top%next)
       yyparse_cmap_topE(i) = top
    enddo
  end subroutine fmod_alloc_yyparse_CMAP_topE
!----------------------------------------------------------------------
  subroutine fmod_alloc_yyparse_cmap_nE(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(yyparse_cmap_nE(ivalue))
    yyparse_cmap_nE = 0
  end subroutine fmod_alloc_yyparse_cmap_nE
!-----------------------------------------------------------------------
  subroutine add_to_mol_cmap_list &
       &           (ivalue1, atom1, atom2, atom3, atom4,atom5, KOC)
    use mol_info, only : para_start
    implicit none
    integer(4), intent(in) :: ivalue1, KOC
    integer(4), intent(inout) :: atom1, atom2, atom3, atom4, atom5

    atom1 = atom1 + para_start(ivalue1)
    atom2 = atom2 + para_start(ivalue1)
    atom3 = atom3 + para_start(ivalue1)
    atom4 = atom4 + para_start(ivalue1)
    atom5 = atom5 + para_start(ivalue1)

    call add_to_cmap_listA(atom1, atom2, atom3, atom4,atom5, KOC)
    call add_to_cmap_listB(atom1, atom2, atom3, atom4,atom5, KOC)
    call add_to_cmap_listC(atom1, atom2, atom3, atom4,atom5, KOC)
    call add_to_cmap_listD(atom1, atom2, atom3, atom4,atom5, KOC)
    call add_to_cmap_listE(atom1, atom2, atom3, atom4,atom5, KOC)
  end subroutine add_to_mol_cmap_list
!----------------------------------------------------------------------
  subroutine set_md_charmm_CMAPA
    use param
    implicit none
    integer(4) :: i,j,k,maxi(1)
    type(yyparse_cmapA),pointer :: p,freed
    ncmap= yyparse_cmap_totalA
    if (yyparse_cmap_totalA == 0)  return
    maxi = maxloc(yyparse_cmap_nA)
    maxi = yyparse_cmap_nA(maxi(1))
    allocate(atom2A(npara,maxi(1)))
    allocate(atom3A(npara,maxi(1)))
    allocate(atom4A(npara,maxi(1)))
    allocate(atom5A(npara,maxi(1)))
    allocate(KindOfCmapA(npara,maxi(1)))
    allocate(nA(npara))
    !     copying
    !      k = 0
    do i=1, npara
       nA(i) = yyparse_cmap_nA(i)
       p => yyparse_cmap_topA(i)
       do j=1,yyparse_cmap_nA(i)
          p => p%next
          atom2A(i,j) = p%atom2
          atom3A(i,j) = p%atom3
          atom4A(i,j) = p%atom4
          atom5A(i,j) = p%atom5
          KindOfCmapA(i,j) = p%KOC
       enddo
    enddo
    deallocate(yyparse_cmap_topA, yyparse_cmap_nA)
  end subroutine set_md_charmm_CMAPA
!-----------------------------------------------------------------------
  subroutine set_md_charmm_CMAPB
    use param
    implicit none
    integer(4) :: i,j,k,maxi(1)
    type(yyparse_cmapB),pointer :: p,freed
    ncmap= yyparse_cmap_totalB
    if (yyparse_cmap_totalB == 0)  return
    maxi = maxloc(yyparse_cmap_nB)
    maxi = yyparse_cmap_nB(maxi(1))
    allocate(atom1B(npara,maxi(1)))
    allocate(atom3B(npara,maxi(1)))
    allocate(atom4B(npara,maxi(1)))
    allocate(atom5B(npara,maxi(1)))
    allocate(KindOfCmapB(npara,maxi(1)))
    allocate(nB(npara))
    !     copying
    do i=1, npara
       nB(i) = yyparse_cmap_nB(i)
       p => yyparse_cmap_topB(i)
       do j=1,yyparse_cmap_nB(i)
          p => p%next
          atom1B(i,j) = p%atom1
          atom3B(i,j) = p%atom3
          atom4B(i,j) = p%atom4
          atom5B(i,j) = p%atom5
          KindOfCmapB(i,j) = p%KOC
       enddo
    enddo

    deallocate(yyparse_cmap_topB, yyparse_cmap_nB)
  end subroutine set_md_charmm_CMAPB
  !-----------------------------------------------------------------------
  subroutine set_md_charmm_CMAPC
    use param
    implicit none
    integer(4) :: i,j,k,maxi(1)
    type(yyparse_cmapC),pointer :: p,freed
    ncmap= yyparse_cmap_totalC
    if (yyparse_cmap_totalC == 0)  return
    maxi = maxloc(yyparse_cmap_nC)
    maxi = yyparse_cmap_nC(maxi(1))
    allocate(atom1C(npara,maxi(1)))
    allocate(atom2C(npara,maxi(1)))
    allocate(atom4C(npara,maxi(1)))
    allocate(atom5C(npara,maxi(1)))
    allocate(KindOfCmapC(npara,maxi(1)))
    allocate(nC(npara))
    !     copying
    do i=1, npara
       nC(i) = yyparse_cmap_nC(i)
       p => yyparse_cmap_topC(i)
       do j=1,yyparse_cmap_nC(i)
          p => p%next
          atom1C(i,j) = p%atom1
          atom2C(i,j) = p%atom2
          atom4C(i,j) = p%atom4
          atom5C(i,j) = p%atom5
          KindOfCmapC(i,j) = p%KOC
       enddo
    enddo

    deallocate(yyparse_cmap_topC, yyparse_cmap_nC)
  end subroutine set_md_charmm_CMAPC
!-----------------------------------------------------------------------
  subroutine set_md_charmm_CMAPD
    use param
    implicit none
    integer(4) :: i,j,k,maxi(1)
    type(yyparse_cmapD),pointer :: p,freed
    ncmap= yyparse_cmap_totalD
    if (yyparse_cmap_totalD == 0)  return
    maxi = maxloc(yyparse_cmap_nD)
    maxi = yyparse_cmap_nD(maxi(1))
    allocate(atom1D(npara,maxi(1)))
    allocate(atom2D(npara,maxi(1)))
    allocate(atom3D(npara,maxi(1)))
    allocate(atom5D(npara,maxi(1)))
    allocate(KindOfCmapD(npara,maxi(1)))
    allocate(nD(npara))
    !     copying
    do i=1, npara
       nD(i) = yyparse_cmap_nD(i)
       p => yyparse_cmap_topD(i)
       do j=1,yyparse_cmap_nD(i)
          p => p%next
          atom1D(i,j) = p%atom1
          atom2D(i,j) = p%atom2
          atom3D(i,j) = p%atom3
          atom5D(i,j) = p%atom5
          KindOfCmapD(i,j) = p%KOC
       enddo
    enddo

    deallocate(yyparse_cmap_topD, yyparse_cmap_nD)
  end subroutine set_md_charmm_CMAPD
!-----------------------------------------------------------------------
  subroutine set_md_charmm_CMAPE
    use param
    implicit none
    integer(4) :: i,j,k,maxi(1)
    type(yyparse_cmapE),pointer :: p,freed
    ncmap= yyparse_cmap_totalE
    if (yyparse_cmap_totalE == 0)  return
    maxi = maxloc(yyparse_cmap_nE)
    maxi = yyparse_cmap_nE(maxi(1))
    allocate(atom1E(npara,maxi(1)))
    allocate(atom2E(npara,maxi(1)))
    allocate(atom3E(npara,maxi(1)))
    allocate(atom4E(npara,maxi(1)))
    allocate(KindOfCmapE(npara,maxi(1)))
    allocate(nE(npara))
    !     copying
    do i=1, npara
       nE(i) = yyparse_cmap_nE(i)
       p => yyparse_cmap_topE(i)
       do j=1,yyparse_cmap_nE(i)
          p => p%next
          atom1E(i,j) = p%atom1
          atom2E(i,j) = p%atom2
          atom3E(i,j) = p%atom3
          atom4E(i,j) = p%atom4
          KindOfCmapE(i,j) = p%KOC
       enddo
    enddo

    deallocate(yyparse_cmap_topE, yyparse_cmap_nE)
  end subroutine set_md_charmm_CMAPE
!------------------------------------------------------------------------
!>
!! \brief  Subroutine to read CMAP parameter A.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_md_charmm_CMAPA
    use trajectory_org
    use device_numbers, only : f_mdff
    use param
    implicit none
    include 'mpif.h'
    integer(4) :: ierr
    integer(4) :: num=0

    if(myrank==0) then
       read(f_mdff) ncmap
    endif
    call MPI_Bcast(ncmap,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    if(ncmap.eq.0) return

    if(myrank==0) read(f_mdff) num
    call MPI_Bcast(num,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)

    allocate(nA(npara))
    allocate(atom2A(npara,num))
    allocate(atom3A(npara,num))
    allocate(atom4A(npara,num))
    allocate(atom5A(npara,num))
    allocate(KindOfCmapA(npara,num))


    if(myrank==0) then
       read(f_mdff) nA
       read(f_mdff) atom2A
       read(f_mdff) atom3A
       read(f_mdff) atom4A
       read(f_mdff) atom5A
       read(f_mdff) KindOfCmapA
    endif
    call MPI_Bcast(nA,npara,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom2A,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom3A,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom4A,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom5A,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(KindOfCmapA,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
  end subroutine read_md_charmm_CMAPA
!------------------------------------------------------------------------
!>
!! \brief  Subroutine to wrte CMAP parameter A in mdff.bin File.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_md_charmm_CMAPA
    use device_numbers, only : f_mdff
    implicit none
    integer(4) :: num
    write(f_mdff) ncmap
    if(ncmap.eq.0) return
    num = size(atom2A(:,:), 2)
    write(f_mdff) num
    write(f_mdff) nA
    write(f_mdff) atom2A
    write(f_mdff) atom3A
    write(f_mdff) atom4A
    write(f_mdff) atom5A
    write(f_mdff) KindOfCmapA
  end subroutine write_md_charmm_CMAPA
!------------------------------------------------------------------------
!>
!! \brief  Subroutine to wrte CMAP parameter A.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_md_charmm_CMAPA
    implicit none
    integer(4) :: num
    write(*,*) '[write_md_charmm_CMAPA]'
    write(*,*) ncmap
    if(ncmap.eq.0) return
    num = size(atom2A(:,:), 2)
    write(*,*) num
    write(*,*) nA
    write(*,*) atom2A
    write(*,*) atom3A
    write(*,*) atom4A
    write(*,*) atom5A
    write(*,*) KindOfCmapA
  end subroutine write_memory_md_charmm_CMAPA
!------------------------------------------------------------------------
!>
!! \brief  Subroutine to read CMAP parameter B.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_md_charmm_CMAPB
    use trajectory_org
    use device_numbers, only : f_mdff
    use param
    implicit none
    include 'mpif.h'
    integer(4) :: ierr
    integer(4) :: num=0

    if(myrank==0) then
       read(f_mdff) ncmap
    endif
    call MPI_Bcast(ncmap,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    if(ncmap.eq.0) return

    if(myrank==0) read(f_mdff) num
    call MPI_Bcast(num,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)

    allocate(nB(npara))
    allocate(atom1B(npara,num))
    allocate(atom3B(npara,num))
    allocate(atom4B(npara,num))
    allocate(atom5B(npara,num))
    allocate(KindOfCmapB(npara,num))

    if(myrank==0) then
       read(f_mdff) nB
       read(f_mdff) atom1B
       read(f_mdff) atom3B
       read(f_mdff) atom4B
       read(f_mdff) atom5B
       read(f_mdff) KindOfCmapB
    endif
    call MPI_Bcast(nB,npara,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom1B,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom3B,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom4B,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom5B,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(KindOfCmapB,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
  end subroutine read_md_charmm_CMAPB
!------------------------------------------------------------------------
!>
!! \brief  Subroutine to wrte CMAP parameter B in mdff.bin File.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_md_charmm_CMAPB
    use device_numbers, only : f_mdff
    implicit none
    integer(4) :: num
    write(f_mdff) ncmap
    if(ncmap.eq.0) return
    num = size(atom1B(:,:), 2)
    write(f_mdff) num
    write(f_mdff) nB
    write(f_mdff) atom1B
    write(f_mdff) atom3B
    write(f_mdff) atom4B
    write(f_mdff) atom5B
    write(f_mdff) KindOfCmapB
  end subroutine write_md_charmm_CMAPB
!------------------------------------------------------------------------
!>
!! \brief  Subroutine to wrte CMAP parameter B.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_md_charmm_CMAPB
    implicit none
    integer(4) :: num
    write(*,*) '[write_md_charmm_CMAPB]'
    write(*,*) ncmap
    if(ncmap.eq.0) return
    num = size(atom1B(:,:), 2)
    write(*,*) num
    write(*,*) nB
    write(*,*) atom1B
    write(*,*) atom3B
    write(*,*) atom4B
    write(*,*) atom5B
    write(*,*) KindOfCmapB
  end subroutine write_memory_md_charmm_CMAPB
!------------------------------------------------------------------------
!>
!! \brief  Subroutine to read CMAP parameter C.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_md_charmm_CMAPC
    use trajectory_org
    use device_numbers, only : f_mdff
    use param
    implicit none
    include 'mpif.h'
    integer(4) :: ierr
    integer(4) :: num=0

    if(myrank==0) then
       read(f_mdff) ncmap
    endif
    call MPI_Bcast(ncmap,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    if(ncmap.eq.0) return

    if(myrank==0) read(f_mdff) num
    call MPI_Bcast(num,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)

    allocate(nC(npara))
    allocate(atom1C(npara,num))
    allocate(atom2C(npara,num))
    allocate(atom4C(npara,num))
    allocate(atom5C(npara,num))
    allocate(KindOfCmapC(npara,num))

    if(myrank==0) then
       read(f_mdff) nC
       read(f_mdff) atom1C
       read(f_mdff) atom2C
       read(f_mdff) atom4C
       read(f_mdff) atom5C
       read(f_mdff) KindOfCmapC
    endif
    call MPI_Bcast(nC,npara,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom1C,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom2C,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom4C,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom5C,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(KindOfCmapC,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
  end subroutine read_md_charmm_CMAPC
!------------------------------------------------------------------------
!>
!! \brief  Subroutine to wrte CMAP parameter C in mdff.bin File.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_md_charmm_CMAPC
    use device_numbers, only : f_mdff
    implicit none
    integer(4) :: num
    write(f_mdff) ncmap
    if(ncmap.eq.0) return
    num = size(atom1C(:,:), 2)
    write(f_mdff) num
    write(f_mdff) nC
    write(f_mdff) atom1C
    write(f_mdff) atom2C
    write(f_mdff) atom4C
    write(f_mdff) atom5C
    write(f_mdff) KindOfCmapC
  end subroutine write_md_charmm_CMAPC
!------------------------------------------------------------------------
!>
!! \brief  Subroutine to wrte CMAP parameter C.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_md_charmm_CMAPC
    implicit none
    integer(4) :: num
    write(*,*) '[write_md_charmm_CMAPC]'
    write(*,*) ncmap
    if(ncmap.eq.0) return
    num = size(atom1C(:,:), 2)
    write(*,*) num
    write(*,*) nC
    write(*,*) atom1C
    write(*,*) atom2C
    write(*,*) atom4C
    write(*,*) atom5C
    write(*,*) KindOfCmapC
  end subroutine write_memory_md_charmm_CMAPC
!------------------------------------------------------------------------
!>
!! \brief  Subroutine to read CMAP parameter D.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_md_charmm_CMAPD
    use trajectory_org
    use device_numbers, only : f_mdff
    use param
    implicit none
    include 'mpif.h'
    integer(4) :: ierr
    integer(4) :: num=0

    if(myrank==0) then
       read(f_mdff) ncmap
    endif
    call MPI_Bcast(ncmap,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    if(ncmap.eq.0) return

    if(myrank==0) read(f_mdff) num
    call MPI_Bcast(num,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)

    allocate(nD(npara))
    allocate(atom1D(npara,num))
    allocate(atom2D(npara,num))
    allocate(atom3D(npara,num))
    allocate(atom5D(npara,num))
    allocate(KindOfCmapD(npara,num))


    if(myrank==0) then
       read(f_mdff) nD
       read(f_mdff) atom1D
       read(f_mdff) atom2D
       read(f_mdff) atom3D
       read(f_mdff) atom5D
       read(f_mdff) KindOfCmapD
    endif
    call MPI_Bcast(nD,npara,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom1D,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom2D,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom3D,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom5D,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(KindOfCmapD,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
  end subroutine read_md_charmm_CMAPD
!------------------------------------------------------------------------
!>
!! \brief  Subroutine to wrte CMAP parameter D in mdff.bin File.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_md_charmm_CMAPD
    use device_numbers, only : f_mdff
    implicit none
    integer(4) :: num
    write(f_mdff) ncmap
    if(ncmap.eq.0) return
    num = size(atom1D(:,:), 2)
    write(f_mdff) num
    write(f_mdff) nD
    write(f_mdff) atom1D
    write(f_mdff) atom2D
    write(f_mdff) atom3D
    write(f_mdff) atom5D
    write(f_mdff) KindOfCmapD
  end subroutine write_md_charmm_CMAPD
!------------------------------------------------------------------------
!>
!! \brief  Subroutine to wrte CMAP parameter D.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_md_charmm_CMAPD
    implicit none
    integer(4) :: num
    write(*,*) '[write_md_charmm_CMAPD]'
    write(*,*) ncmap
    if(ncmap.eq.0) return
    num = size(atom1D(:,:), 2)
    write(*,*) num
    write(*,*) nD
    write(*,*) atom1D
    write(*,*) atom2D
    write(*,*) atom3D
    write(*,*) atom5D
    write(*,*) KindOfCmapD
  end subroutine write_memory_md_charmm_CMAPD
!------------------------------------------------------------------------
!>
!! \brief  Subroutine to read CMAP parameter E.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_md_charmm_CMAPE
    use trajectory_org
    use device_numbers, only : f_mdff
    use param
    implicit none
    include 'mpif.h'
    integer(4) :: ierr
    integer(4) :: num=0

    if(myrank==0) then
       read(f_mdff) ncmap
    endif
    call MPI_Bcast(ncmap,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    if(ncmap.eq.0) return

    if(myrank==0) read(f_mdff) num
    call MPI_Bcast(num,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)

    allocate(nE(npara))
    allocate(atom1E(npara,num))
    allocate(atom2E(npara,num))
    allocate(atom3E(npara,num))
    allocate(atom4E(npara,num))
    allocate(KindOfCmapE(npara,num))

    if(myrank==0) then
       read(f_mdff) nE
       read(f_mdff) atom1E
       read(f_mdff) atom2E
       read(f_mdff) atom3E
       read(f_mdff) atom4E
       read(f_mdff) KindOfCmapE
    endif
    call MPI_Bcast(nE,npara,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom1E,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom2E,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom3E,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(atom4E,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(KindOfCmapE,npara*num,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
  end subroutine read_md_charmm_CMAPE
!------------------------------------------------------------------------
!>
!! \brief  Subroutine to wrte CMAP parameter E in mdff.bin File.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_md_charmm_CMAPE
    use device_numbers, only : f_mdff
    implicit none
    integer(4) :: num
    write(f_mdff) ncmap
    if(ncmap.eq.0) return
    num = size(atom1E(:,:), 2)
    write(f_mdff) num
    write(f_mdff) nE
    write(f_mdff) atom1E
    write(f_mdff) atom2E
    write(f_mdff) atom3E
    write(f_mdff) atom4E
    write(f_mdff) KindOfCmapE
  end subroutine write_md_charmm_CMAPE
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to wrte CMAP parameter E.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_md_charmm_CMAPE
    implicit none
    integer(4) :: num
    write(*,*) '[write_md_charmm_CMAPE]'
    write(*,*) ncmap
    if(ncmap.eq.0) return
    num = size(atom1E(:,:), 2)
    write(*,*) num
    write(*,*) nE
    write(*,*) atom1E
    write(*,*) atom2E
    write(*,*) atom3E
    write(*,*) atom4E
    write(*,*) KindOfCmapE
  end subroutine write_memory_md_charmm_CMAPE
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to read CMAP Version.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_CMAPVersion
    use md_condition
    use device_numbers, only : f_mdff
    use mpi_tool
    implicit none
    include 'mpif.h'
    integer(4) :: ierr

    if(myrank==0) read(f_mdff,end=100) CMAPVersion
    call MPI_Bcast(CMAPVersion,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
    return

100 continue

    call MPI_Bcast(CMAPVersion,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
  end subroutine read_CMAPVersion
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write CMAP Version in mdff.bin File.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_CMAPVersion
    use device_numbers, only : f_mdff
    implicit none
    write(f_mdff) CMAPVersion
  end subroutine write_CMAPVersion
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write CMAP Version.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_CMAPVersion
    implicit none
    write(*,*) '[write_CMAPVersion]'
    write(*,*) CMAPVersion
  end subroutine write_memory_CMAPVersion

!-----------------------------------------------------------------------
  subroutine add_charmm_CMAP_a()
    !fujimoto
    !
    !
    !    *
    !    *             C2----N2
    !    *            /
    !    *    N1----CA
    !    *   /
    !    * C1
    !    *
    !    *phi=C1-N1-CA-C2
    !    *psi=N1-CA-C2-N2
    !
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
    integer(4) :: i, j,k,ip,ia, ib, i0a, i0b, i0c, i0d,i0e,i0,ipar, k0
    real(8) :: UCMAP,oUc
    real(8) :: rcmap(3,5), Fcmap(3,5), Vcmap
    real(8) :: dx, dy, dz
    integer(4) :: iam
    real(8) :: wkvd(6), vd(6) ! debug
    real(8) :: v11,v22,v33
    real(8) :: wk_v11,wk_v22,wk_v33
    real(8) :: v21,v31,v32
    real(8) :: wk_v21, wk_v31, wk_v32
    real(8) :: wvr(6),vr(6) ! debug
    include 'mpif.h'

    iam = 0
    UCMAP = 0.0d0
    wk_v11=0d0;wk_v22=0d0;wk_v33=0d0
    wk_v21=0d0;wk_v31=0d0;wk_v32=0d0

!$omp parallel default(none) &
!$omp& private(k0,i,iam,j,i0,ipar,k) &
!$omp& private(ip,ia,ib,i0a,i0b,i0c,i0d,i0e) &
!$omp& private(rcmap,Fcmap,Vcmap,dx,dy,dz) &
!$omp& shared(paranum) &
!$omp& shared(m2i,i2m) &
!$omp& shared(nA,atom2A,atom3A,atom4A,atom5A,KindOfCmapA) &
!$omp& shared(nB,atom1B,atom3B,atom4B,atom5B,KindOfCmapB) &
!$omp& shared(nC,atom1C,atom2C,atom4C,atom5C,KindOfCmapC) &
!$omp& shared(nD,atom1D,atom2D,atom3D,atom5D,KindOfCmapD) &
!$omp& shared(nE,atom1E,atom2E,atom3E,atom4E,KindOfCmapE) &
!$omp& shared(wkxyz,w3_f,wk_vir2) &
!$omp& shared(lsegtop,lseg_natoms,nselfseg) &
!$omp& private(v11,v22,v33,v21,v31,v32) &
!$omp& reduction(+:wk_v11,wk_v22,wk_v33) &
!$omp& reduction(+:wk_v21,wk_v31,wk_v32) &
!$omp& reduction(+:UCMAP)
!$  iam = omp_get_thread_num()
!$omp do
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          ia   = m2i(i0)
          ipar = paranum(ia)

          ! C group: CA
          do j = 1, nC(ipar)
             i0b = i2m(atom1C(ipar,j)+(ia-ipar)) !C1
             i0c = i2m(atom2C(ipar,j)+(ia-ipar)) !N1
             i0d = i2m(atom4C(ipar,j)+(ia-ipar)) !C2
             i0e = i2m(atom5C(ipar,j)+(ia-ipar)) !N2

             !C1
             rcmap(1,1) = wkxyz(1,i0b)
             rcmap(2,1) = wkxyz(2,i0b)
             rcmap(3,1) = wkxyz(3,i0b)

             !N1
             rcmap(1,2) = wkxyz(1,i0c)
             rcmap(2,2) = wkxyz(2,i0c)
             rcmap(3,2) = wkxyz(3,i0c)


             !CA
             rcmap(1,3) = wkxyz(1,i0 )
             rcmap(2,3) = wkxyz(2,i0 )
             rcmap(3,3) = wkxyz(3,i0 )

             !C2
             rcmap(1,4) = wkxyz(1,i0d)
             rcmap(2,4) = wkxyz(2,i0d)
             rcmap(3,4) = wkxyz(3,i0d)

             !N2
             rcmap(1,5) = wkxyz(1,i0e)
             rcmap(2,5) = wkxyz(2,i0e)
             rcmap(3,5) = wkxyz(3,i0e)

             do k = 1, 5
                dx = rcmap(1,k) - rcmap(1,3)
                dy = rcmap(2,k) - rcmap(2,3)
                dz = rcmap(3,k) - rcmap(3,3)
#ifdef ONEPROC_AXIS
                call pbc_pair(dx,dy,dz)
#endif
                rcmap(1,k) = dx + rcmap(1,3)
                rcmap(2,k) = dy + rcmap(2,3)
                rcmap(3,k) = dz + rcmap(3,3)
             enddo

             call MainCMAP(rcmap,KindOfCmapC(ipar,j),Fcmap,Vcmap)

             w3_f(1,i0 ,iam) = w3_f(1,i0 ,iam) + Fcmap(1,3)
             w3_f(2,i0 ,iam) = w3_f(2,i0 ,iam) + Fcmap(2,3)
             w3_f(3,i0 ,iam) = w3_f(3,i0 ,iam) + Fcmap(3,3)
!default: MTD
             w3_f(1,i0b,iam) = w3_f(1,i0b,iam) + Fcmap(1,1)
             w3_f(2,i0b,iam) = w3_f(2,i0b,iam) + Fcmap(2,1)
             w3_f(3,i0b,iam) = w3_f(3,i0b,iam) + Fcmap(3,1)
             w3_f(1,i0c,iam) = w3_f(1,i0c,iam) + Fcmap(1,2)
             w3_f(2,i0c,iam) = w3_f(2,i0c,iam) + Fcmap(2,2)
             w3_f(3,i0c,iam) = w3_f(3,i0c,iam) + Fcmap(3,2)
             w3_f(1,i0d,iam) = w3_f(1,i0d,iam) + Fcmap(1,4)
             w3_f(2,i0d,iam) = w3_f(2,i0d,iam) + Fcmap(2,4)
             w3_f(3,i0d,iam) = w3_f(3,i0d,iam) + Fcmap(3,4)
             w3_f(1,i0e,iam) = w3_f(1,i0e,iam) + Fcmap(1,5)
             w3_f(2,i0e,iam) = w3_f(2,i0e,iam) + Fcmap(2,5)
             w3_f(3,i0e,iam) = w3_f(3,i0e,iam) + Fcmap(3,5)
             v11 =       rcmap(1,1)*Fcmap(1,1) + rcmap(1,2)*Fcmap(1,2) &
                  &              + rcmap(1,3)*Fcmap(1,3) + rcmap(1,4)*Fcmap(1,4) &
                  &              + rcmap(1,5)*Fcmap(1,5)
             v22 =       rcmap(2,1)*Fcmap(2,1) + rcmap(2,2)*Fcmap(2,2) &
                  &              + rcmap(2,3)*Fcmap(2,3) + rcmap(2,4)*Fcmap(2,4) &
                  &              + rcmap(2,5)*Fcmap(2,5)
             v33 =       rcmap(3,1)*Fcmap(3,1) + rcmap(3,2)*Fcmap(3,2) &
                  &              + rcmap(3,3)*Fcmap(3,3) + rcmap(3,4)*Fcmap(3,4) &
                  &              + rcmap(3,5)*Fcmap(3,5)
             v21 =       rcmap(2,1)*Fcmap(1,1) + rcmap(2,2)*Fcmap(1,2) &
                  &              + rcmap(2,3)*Fcmap(1,3) + rcmap(2,4)*Fcmap(1,4) &
                  &              + rcmap(2,5)*Fcmap(1,5)
             v31 =       rcmap(3,1)*Fcmap(1,1) + rcmap(3,2)*Fcmap(1,2) &
                  &              + rcmap(3,3)*Fcmap(1,3) + rcmap(3,4)*Fcmap(1,4) &
                  &              + rcmap(3,5)*Fcmap(1,5)
             v32 =       rcmap(3,1)*Fcmap(2,1) + rcmap(3,2)*Fcmap(2,2) &
                  &              + rcmap(3,3)*Fcmap(2,3) + rcmap(3,4)*Fcmap(2,4) &
                  &              + rcmap(3,5)*Fcmap(2,5)
             wk_v11=wk_v11+v11
             wk_v22=wk_v22+v22
             wk_v33=wk_v33+v33
             wk_v21=wk_v21+v21
             wk_v31=wk_v31+v31
             wk_v32=wk_v32+v32
             UCMAP = UCMAP + Vcmap  ! 2020/3/16
          enddo ! j

       enddo ! i0
    enddo ! k0
!$omp end do
!$omp end parallel

    wk_p_energy = wk_p_energy + UCMAP
    wk_vir2(1,0) = wk_vir2(1,0) + wk_v11
    wk_vir2(2,0) = wk_vir2(2,0) + wk_v22
    wk_vir2(3,0) = wk_vir2(3,0) + wk_v33
    wk_vir2(4,0) = wk_vir2(4,0) + wk_v21
    wk_vir2(5,0) = wk_vir2(5,0) + wk_v31
    wk_vir2(6,0) = wk_vir2(6,0) + wk_v32
#ifdef DEBUGFCE
    oUc=0d0
    call mpi_allreduce(Ucmap,oUc,1, &
    &     mpi_double_precision,mpi_sum,mpi_comm_world,ipar)
#ifdef KCAL
    if(myrank==0)write(*,*)'Pot(CMAP )=',oUc *kJ_mol/4.184d0,'[kcal/mol]'
#else
    if(myrank==0)write(*,*)'Pot(CMAP )=',oUc *kJ_mol,'[kJ/mol]'
#endif
#endif
  end subroutine add_charmm_CMAP_a

end module CMAP
