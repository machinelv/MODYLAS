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
!! \brief Modules and subroutines which apply bond length constraints 
!!        by the PSHAKE/PRATTLE method.
!<
!----------------------------------------------------------------------
!>
!! \brief Modules to parse parameters for PSHAKE/PRATTLE method.
!! \author Hidekazu Kojima, and partially Yoshimichi Andoh
!<
module parse_shake

#include "timing.h90"
  implicit none
  type yyparse_shake
     integer(4) :: atom1,atom2
     real(8) :: r0
     type(yyparse_shake),pointer :: next
  end type yyparse_shake
  type(yyparse_shake),pointer :: yyparse_shake_top(:)
  integer(4), allocatable :: yyparse_shake_n(:)

contains

  !----------------------------------------------------------------------
  subroutine fmod_alloc_yyparse_shake(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    integer(4) :: i
    type(yyparse_shake),pointer :: top
    allocate(yyparse_shake_top(ivalue))
    allocate(yyparse_shake_n(ivalue))
    do i=1,ivalue
       allocate(top)
       top%atom1 = -1
       top%atom2 = -1
       nullify(top%next)
       yyparse_shake_top(i) = top
    enddo
  end subroutine fmod_alloc_yyparse_shake
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to keep shake data in liner list temporally.
!! \author Kensuke Iwahashi
!<
  subroutine add_to_shake_list(atom01,atom02,r0,line_number)
    use segments, only : seg_ready, segment_is_not_ready, irearrange
    use param, only : npara
    use mpi_tool
    implicit none
    integer(4) :: atom01, atom02
    integer(4) :: atom1,atom2,atom_small,atom_large,line_number
    real(8) :: r0
    type(yyparse_shake),pointer :: new1,new2,p,next
    if (seg_ready) then
       atom1 = irearrange(atom01 + 1)
       atom2 = irearrange(atom02 + 1)
    else
       call segment_is_not_ready
       return
    endif
    if (atom1 > npara) then
       write(0,'(a,i0)')  'ERROR (line ', line_number, ') : '
       write(0,'(a,i0,a)')  &
            &    'The number of shake.atom1 is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    if (atom2 > npara) then
       write(0,'(a,i0)')  'ERROR (line ', line_number, ') : '
       write(0,'(a,i0,a)')  &
            &    'The number of shake.atom2 is out of bounds.  '// &
            &    'It must be less than ', npara, '.'
       call modylas_abort()
    endif
    atom_small = min(atom1, atom2)
    atom_large = max(atom1, atom2)
    !
    allocate(new1)
    new1%atom1 = atom_small
    new1%atom2 = atom_large
    new1%r0 = r0
    nullify(new1%next)
    p => yyparse_shake_top(atom_small)
    !     insert to list, sorting from smaller
    do while (.true.)
       next => p%next
       if (.not. associated(next)) then
          p%next => new1
          exit
       endif
       !       doubled must be merged
       if (atom_large == next%atom2) then
          if (dabs(r0 - next%r0) < 1.0d-50) then
             deallocate(new1)
             return
          else
             write(0,'(a,i0)')  'ERROR (line ', line_number, ') : '
             write(0,'(a,i0,a)')  &
                  &        'Detect shake entry defined before: ', &
                  &        atom_small, '-', atom_large
             call modylas_abort()
          endif
       else if (atom_large < next%atom2) then
          new1%next => p%next
          p%next => new1
          exit
       endif
       p => next
    enddo
    !     yyparse_shake_n(atom_small) = yyparse_shake_n(atom_small) + 1
    !
    !     allocate(new2)
    !     new2%atom1 = atom_large
    !     new2%atom2 = atom_small
    !     new2%r0 = r0
    !     nullify(new2%next)
    !     p => yyparse_shake_top(atom_large)
    !     insert to list, sorting from smaller
    !     do while (.true.)
    !       next => p%next
    !       if (.not. associated(next)) then
    !         p%next => new2
    !         exit
    !       endif
    !       doubled must be merged
    !       if (atom_small == next%atom2) then
    !         if (dabs(r0 - next%r0) < 1.0d-50) then
    !           deallocate(new2)
    !           return
    !         else
    !           call abort_with_message_aiai( &
    !    &        'Detect shake entry defined before: ', &
    !    &        atom_large, '-', atom_small)
    !         endif
    !       else if (atom_small < next%atom2) then
    !         new2%next => p%next
    !         p%next => new2
    !         exit
    !       endif
    !       p => next
    !     enddo
    !yyparse_shake_n(atom_large) = yyparse_shake_n(atom_large) + 1
    !
    !yyparse_shake_total = yyparse_shake_total + 1
  end subroutine add_to_shake_list

end module parse_shake
!>
!! \brief Module to store arrays used in PSHAKE/PRATTLE method.
!! \author Hidekazu Kojima
!<
module pshake_init
  integer(4):: rngrp_ps
  real(8),allocatable::r_init_ps(:,:,:)
  integer(4),allocatable:: n_cnst_ps(:)
end module pshake_init
!>
!! \brief Module to apply PSHAKE/PRATTLE method.
!! \author Hidekazu Kojima
!<
module shake_rattle_roll
  use omp_lib
  use mpi_tool
  implicit none
  real(8),allocatable :: xyzstr(:,:)
  real(8),allocatable :: vstr(:,:)
  real(8),allocatable :: dij(:),lambstr(:),lamb(:)
  integer(4),allocatable :: shkij(:,:)
  integer(4) :: totnconst=0
  integer(4) :: totnconstL,maxibn,l0max
  integer(4) :: totnconstL_init=250
  integer(4),allocatable :: ibseL(:),shkijL(:,:,:)
  real(8),allocatable :: rmassL(:,:,:)
  real(8),allocatable :: lambL(:,:),lambLstr(:,:)
  real(8),allocatable :: dij2L(:,:)
  integer(4),allocatable :: igshakeL(:)
  integer(4), allocatable :: ShakeGroupLeader(:)
  integer(4), allocatable :: nconstraints(:)
  integer(4), allocatable :: atom1S(:,:), atom2S(:,:)
  real(8), allocatable :: slength(:,:)
  real(8),allocatable :: wk_dfc(:,:)
  real(8),allocatable :: wkdvir2(:,:)
  real(8),allocatable :: wkdvirSca(:) !scaler virial by constraint
  real(8),allocatable :: wkdvirTen(:,:) !virial tensor by constraint
  integer(4),allocatable :: ishake(:,:),kshake(:),ibond(:,:)
  real(8) :: shake_tolerance=1.0d-10,rattle_tolerance=1.0d-26
  real(8) :: roll_tolerance=1d-3
  integer(4) :: maxshakecycle=1000,ngshake=0
  real(8) :: wk_roll_flg, roll_flg

! integer(4),parameter:: n_cnst_max=20
  integer(4),parameter:: n_cnst_max=10
! integer(4),parameter:: n_cnst_max=4
  integer(4),parameter:: n_atom_max=20
  integer(4):: n_type_ps
  integer(4),allocatable:: type_psL(:)
  integer(4),allocatable:: gn2ips(:)
  integer(4),allocatable:: couple_ps(:,:,:), n_atom_ps(:)
  real(8),allocatable:: a_0(:,:,:), a_0_sym(:,:,:)
  real(8),allocatable:: rdij2_ps(:,:), mass_ps(:,:)

  real(8) :: sgamma, rgamma, drCOM,MForceA,MForceB,MF_AB
  real(8) :: sMF_AB=0d0

contains

  subroutine fmod_alloc_smp_wk
    use omp_lib
    use md_multiplestep
    implicit none
    integer(4) :: nomp = 1
!$  nomp = omp_get_max_threads()
    allocate(wkdvir2(6,0:nomp-1)) ! for SHAKE
    allocate(wkdvirSca(2*maxMTm*maxMTl)) ! for SHAKE
    wkdvir2=0d0
    wkdvirSca=0d0
    allocate(wkdvirTen(6,2*maxMTm*maxMTl)) ! for SHAKE
    wkdvirTen=0d0
  end subroutine fmod_alloc_smp_wk
!----------------------------------------------------------------------
  subroutine fmod_set_totnconstL(ivalue)
    use mpi_tool
    implicit none
    integer(4), intent(in) :: ivalue
    totnconstL_init = ivalue
    if(myrank==0) write(*,*) 'WARNING: totnconstL is set to ', ivalue
  end subroutine fmod_set_totnconstL
!----------------------------------------------------------------------
  subroutine shake_roll(dtin)
    use trajectory_mpi
    implicit none
    real(8) :: shktol
    integer(4) :: iterchk, iam
    integer(4) :: i_ps,l0,l1,iconstraints
    real(8) :: dtin,rdtin,rdtbin,rdt2in

    iam=0

    rdtin =1d0/dtin
    rdtbin=rdtin+rdtin
    rdt2in=rdtin*rdtin
    !
    iterchk=0
    shktol = shake_tolerance

!$omp parallel default(shared) &
!$omp& private(iam,l0,l1,i_ps,iconstraints)
!$    iam=omp_get_thread_num()
!$omp do
    do l0=1,l0max
       i_ps = type_psL(l0)
       !        if(i_ps == 0)then
       !          call  shake_roll_main(l0, rdtin, rdtbin, iterchk,shktol, iam)
       !        else

       call pshake_roll_main(l0, rdtin, rdtbin, iterchk,shktol,i_ps, iam)
       !        endif
    end do
!$omp end do
!$omp do
    do l0=1,l0max
       iconstraints=ibseL(l0)
       do l1=1,iconstraints
          lambL(l1,l0)=lambL(l1,l0)*rdt2in
       enddo
    enddo
!$omp end do
!$omp end parallel
    if(iterchk.ne.0) then
       write(0,*) 'ERROR: SHAKE/ROLL iteration no convergence.'
       call modylas_abort()
    end if

    return
  end subroutine shake_roll
  !----------------------------------------------------------------------
  !
  !     prattle, rattle/roll
  !
  !----------------------------------------------------------------------
  subroutine rattle_roll(dtin)
    implicit none
    real(8) :: shktol
    integer(4) :: iterchk, iam
    integer(4) :: i_ps,l0,l1,iconstraints
    real(8) :: dtin,rdtin,rdtbin

    iam = 0
    !
    rdtin =1d0/dtin
    rdtbin=rdtin+rdtin
    !
    iterchk=0
    shktol = shake_tolerance

!$omp parallel default(shared) &
!$omp& private(iam,l0,l1,i_ps,iconstraints)
!$    iam = omp_get_thread_num()
!$omp do
    do l0=1,l0max
       i_ps = type_psL(l0)
       !        if(i_ps == 0)then
       !          call  rattle_roll_main(l0,dtin,rdtbin,iterchk,shktol,iam)
       !        else
       call prattle_roll_main(l0,dtin,rdtbin,iterchk,shktol,i_ps,iam)
       !        endif
    end do
!$omp end do
!$omp do
    do l0=1,l0max
       iconstraints=ibseL(l0)
       do l1=1,iconstraints
          lambL(l1,l0)=lambL(l1,l0)*rdtin
       enddo
    enddo
!$omp end do
!$omp end parallel
    if(iterchk.ne.0) then
       write(0,*) 'ERROR: RATTLE/ROLL iteration no convergence.'
       call modylas_abort()
    end if

  end subroutine rattle_roll
!----------------------------------------------------------------------
!
!     shake/roll
!
!----------------------------------------------------------------------
  subroutine shake_roll_main(l0, rdtin, rdtbin, iterchk, shktol, iam)
    use trajectory_mpi
    implicit none
    integer(4) :: ida,idb
    integer(4) :: icycle
    real(8) :: dx1,dy1,dz1
    real(8) :: dx0,dy0,dz0
    real(8) :: shktol,spro
    real(8) :: rij2,sabun
    real(8) :: tmp,gmi,gmj
    real(8) :: deltax,deltay,deltaz
    real(8) :: rdtin,rdtbin
    integer(4) :: l0,l1,iconstraints
    integer(4) :: iflgloop
    integer(4) :: iterchk,iam

    iconstraints=ibseL(l0)

    do icycle=1,maxshakecycle
       iflgloop=0
       do l1 = 1, iconstraints
          ida = shkijL(1,l1,l0)
          idb = shkijL(2,l1,l0)
          dx1 = wkxyz(1,idb)-wkxyz(1,ida) ! L
          dy1 = wkxyz(2,idb)-wkxyz(2,ida)
          dz1 = wkxyz(3,idb)-wkxyz(3,ida)
          rij2 = dx1**2+dy1**2+dz1**2
          sabun = dij2L(l1,l0)-rij2 ! L
          if (abs(sabun)/dij2L(l1,l0) .gt. shktol) then  ! L
             dx0 = xyzstr(1,idb)-xyzstr(1,ida) ! L
             dy0 = xyzstr(2,idb)-xyzstr(2,ida)
             dz0 = xyzstr(3,idb)-xyzstr(3,ida)
             spro = dx1*dx0+dy1*dy0+dz1*dz0
             gmi = rmassL(1,l1,l0) ! L
             gmj = rmassL(2,l1,l0) ! L
             tmp = -sabun/(2d0*(gmi+gmj)*spro)
             lambL(l1,l0)=lambL(l1,l0)+tmp ! L
             deltax = dx0*tmp
             deltay = dy0*tmp
             deltaz = dz0*tmp
             wkxyz(1,ida) = wkxyz(1,ida)+deltax*gmi ! L
             wkxyz(2,ida) = wkxyz(2,ida)+deltay*gmi
             wkxyz(3,ida) = wkxyz(3,ida)+deltaz*gmi
             wkxyz(1,idb) = wkxyz(1,idb)-deltax*gmj ! L
             wkxyz(2,idb) = wkxyz(2,idb)-deltay*gmj
             wkxyz(3,idb) = wkxyz(3,idb)-deltaz*gmj
             deltax = deltax*rdtin
             deltay = deltay*rdtin
             deltaz = deltaz*rdtin
             wkv(1,ida) = wkv(1,ida)+deltax*gmi ! L
             wkv(2,ida) = wkv(2,ida)+deltay*gmi
             wkv(3,ida) = wkv(3,ida)+deltaz*gmi
             wkv(1,idb) = wkv(1,idb)-deltax*gmj ! L
             wkv(2,idb) = wkv(2,idb)-deltay*gmj
             wkv(3,idb) = wkv(3,idb)-deltaz*gmj
             deltax = deltax*rdtbin
             deltay = deltay*rdtbin
             deltaz = deltaz*rdtbin
             wk_dfc(1,ida) = wk_dfc(1,ida)+deltax ! L
             wk_dfc(2,ida) = wk_dfc(2,ida)+deltay
             wk_dfc(3,ida) = wk_dfc(3,ida)+deltaz
             wk_dfc(1,idb) = wk_dfc(1,idb)-deltax ! L
             wk_dfc(2,idb) = wk_dfc(2,idb)-deltay
             wk_dfc(3,idb) = wk_dfc(3,idb)-deltaz
             if(ida.lt.idb)then ! avoid double count
                wkdvir2(1,iam)=wkdvir2(1,iam)+dx0*deltax
                wkdvir2(2,iam)=wkdvir2(2,iam)+dy0*deltay
                wkdvir2(3,iam)=wkdvir2(3,iam)+dz0*deltaz
                wkdvir2(4,iam)=wkdvir2(4,iam)+dy0*deltax
                wkdvir2(5,iam)=wkdvir2(5,iam)+dz0*deltax
                wkdvir2(6,iam)=wkdvir2(6,iam)+dz0*deltay
             endif
          else
             iflgloop=iflgloop+1
          end if
       end do ! l1
       if(iflgloop.eq.iconstraints) then
          exit
       endif
    end do ! icycle

    if (icycle .ge. maxshakecycle) then
       iterchk = iterchk + 1
    end if
  end subroutine shake_roll_main
!----------------------------------------------------------------------
!
!     pshake/roll main
!
!----------------------------------------------------------------------
!>
!! \brief PSHAKE reaference: P. Gonnet, J. Comput. Phys. 220, 740 (2007).
!<
  subroutine pshake_roll_main(l0, rdtin, rdtbin, iterchk,shktol,i_ps, iam)
    use trajectory_mpi
    implicit none
    integer(4) :: i,ida,idb, j,l0,l1,iconstraints
    integer(4) :: icycle
    real(8) :: rdtin,rdtbin
    real(8) :: shktol,spro
    real(8) :: rij2
    real(8) :: gmi,gmj
    real(8) ::deltax(n_atom_max),deltay(n_atom_max),deltaz(n_atom_max)
    integer(4) :: iflgloop
    integer(4) :: iterchk, iam
    integer(4) :: i_ps, n_atom
    real(8) ::dx0_ps(n_cnst_max),dy0_ps(n_cnst_max),dz0_ps(n_cnst_max)
    real(8) ::dx1_ps(n_cnst_max),dy1_ps(n_cnst_max),dz1_ps(n_cnst_max)
    real(8) :: sabun_ps(n_cnst_max)
    real(8) :: tmp_ps(n_cnst_max)
    real(8) :: v_mat(3,n_cnst_max,n_atom_max)
    real(8) :: v_mat_old(3,n_cnst_max,n_atom_max)
    integer(4) :: atom_ps(n_cnst_max), l, k
    integer(4) :: couple_wk(2,n_cnst_max)
    real(8) :: mass_wk(n_atom_max)
    real(8) :: ma

    iconstraints=ibseL(l0)

    n_atom = n_atom_ps(i_ps)
    atom_ps(1) = shkijL(1,1,l0)
    v_mat_old(:,:iconstraints,:n_atom) = 0.0d0
    do l1=1, iconstraints
       couple_wk(1,l1) = couple_ps(1,l1,i_ps)
       couple_wk(2,l1) = couple_ps(2,l1,i_ps)
    enddo

    do i = 2, n_atom
       atom_ps(i) = atom_ps(1) + i - 1
    enddo

    do i=1, n_atom
       mass_wk(i) = mass_ps(i,i_ps)
    enddo

    do l1= 1, iconstraints !n_pshake = iconstraints
       ida = shkijL(1,l1,l0)
       idb = shkijL(2,l1,l0)
       gmi = rmassL(1,l1,l0)
       gmj = rmassL(2,l1,l0)
       dx0_ps(l1) = xyzstr(1,idb) - xyzstr(1,ida)
       dy0_ps(l1) = xyzstr(2,idb) - xyzstr(2,ida)
       dz0_ps(l1) = xyzstr(3,idb) - xyzstr(3,ida)
       v_mat_old(1,l1,couple_wk(1,l1)) = dx0_ps(l1)*gmi
       v_mat_old(2,l1,couple_wk(1,l1)) = dy0_ps(l1)*gmi
       v_mat_old(3,l1,couple_wk(1,l1)) = dz0_ps(l1)*gmi
       v_mat_old(1,l1,couple_wk(2,l1)) =-dx0_ps(l1)*gmj
       v_mat_old(2,l1,couple_wk(2,l1)) =-dy0_ps(l1)*gmj
       v_mat_old(3,l1,couple_wk(2,l1)) =-dz0_ps(l1)*gmj
    enddo

!ocl unroll(2)
    do j=1, n_atom
       v_mat(:,:iconstraints,j) = 0.0d0
       do i=1, iconstraints
          do l=1, iconstraints
             do k=1, 3
                v_mat(k,i,j) = v_mat(k,i,j) + v_mat_old(k,l,j) * a_0_sym(l,i,i_ps)
             enddo
          enddo
       enddo
    enddo

    do icycle=1,maxshakecycle
       iflgloop = 0
!ocl noswp
       do l1 = 1, iconstraints
          ida = shkijL(1,l1,l0)
          idb = shkijL(2,l1,l0)
          dx1_ps(l1) = wkxyz(1,idb) - wkxyz(1,ida)
          dy1_ps(l1) = wkxyz(2,idb) - wkxyz(2,ida)
          dz1_ps(l1) = wkxyz(3,idb) - wkxyz(3,ida)
          rij2 = dx1_ps(l1)*dx1_ps(l1) + dy1_ps(l1)*dy1_ps(l1) + dz1_ps(l1)*dz1_ps(l1)
          sabun_ps(l1) = dij2L(l1,l0) - rij2
          !judge of a constraint for pshake
!         if (dabs(sabun_ps(l1))/dij2L(l1,l0) .le. shktol) then !(relative val)**2
          if (dabs(sabun_ps(l1)) .le. dij2L(l1,l0)*shktol) then !(relative val)**2
             iflgloop=iflgloop+1
          end if
       enddo

       if(iflgloop.eq.iconstraints) then
          goto 123
       endif

       do l1 = 1, iconstraints
          ida = couple_wk(1,l1)
          idb = couple_wk(2,l1)
          spro = dx1_ps(l1)*(v_mat(1,l1,ida)-v_mat(1,l1,idb))+ &
               &          dy1_ps(l1)*(v_mat(2,l1,ida)-v_mat(2,l1,idb))+ &
               &          dz1_ps(l1)*(v_mat(3,l1,ida)-v_mat(3,l1,idb))
          tmp_ps(l1) = -sabun_ps(l1) * 0.5d0 / spro
          lambL(l1,l0)=lambL(l1,l0)+tmp_ps(l1)
       end do

       do i = 1, n_atom
          deltax(i) = 0.0d0
          deltay(i) = 0.0d0
          deltaz(i) = 0.0d0
          do l1=1, iconstraints
             deltax(i) = deltax(i) + v_mat(1,l1,i) * tmp_ps(l1)
             deltay(i) = deltay(i) + v_mat(2,l1,i) * tmp_ps(l1)
             deltaz(i) = deltaz(i) + v_mat(3,l1,i) * tmp_ps(l1)
          enddo
       enddo
!ocl norecurrence(wkxyz,wkv,wk_dfc)
       do i = 1, n_atom
          ida = atom_ps(i)
          ma = mass_wk(i)
          wkxyz(1,ida) = wkxyz(1,ida) + deltax(i)
          wkxyz(2,ida) = wkxyz(2,ida) + deltay(i)
          wkxyz(3,ida) = wkxyz(3,ida) + deltaz(i)
          deltax(i) = deltax(i) * rdtin
          deltay(i) = deltay(i) * rdtin
          deltaz(i) = deltaz(i) * rdtin
          wkv(1,ida) = wkv(1,ida) + deltax(i)
          wkv(2,ida) = wkv(2,ida) + deltay(i)
          wkv(3,ida) = wkv(3,ida) + deltaz(i)
          deltax(i) = deltax(i) * rdtbin * ma
          deltay(i) = deltay(i) * rdtbin * ma
          deltaz(i) = deltaz(i) * rdtbin * ma
          wk_dfc(1,ida) = wk_dfc(1,ida) + deltax(i)
          wk_dfc(2,ida) = wk_dfc(2,ida) + deltay(i)
          wk_dfc(3,ida) = wk_dfc(3,ida) + deltaz(i)
          wkdvir2(1,iam) = wkdvir2(1,iam) + xyzstr(1,ida)*deltax(i)
          wkdvir2(2,iam) = wkdvir2(2,iam) + xyzstr(2,ida)*deltay(i)
          wkdvir2(3,iam) = wkdvir2(3,iam) + xyzstr(3,ida)*deltaz(i)
          wkdvir2(4,iam) = wkdvir2(4,iam) + xyzstr(2,ida)*deltax(i)
          wkdvir2(5,iam) = wkdvir2(5,iam) + xyzstr(3,ida)*deltax(i)
          wkdvir2(6,iam) = wkdvir2(6,iam) + xyzstr(3,ida)*deltay(i)
       enddo

    end do
123 continue

    if (icycle .ge. maxshakecycle) then
       iterchk = iterchk + 1
    end if
 
  end subroutine pshake_roll_main
!----------------------------------------------------------------------
!
!     rattle/roll
!
!----------------------------------------------------------------------
  subroutine rattle_roll_main(l0,dtin,rdtbin,iterchk,shktol,iam)
    use trajectory_mpi
    implicit none
    integer(4) :: ida,idb
    integer(4) :: icycle
    real(8) :: dx1,dy1,dz1
    real(8) :: dvx,dvy,dvz
    real(8) :: shktol,spro
    real(8) :: gmi,gmj
    real(8) :: tmp,tmp2
    real(8) :: deltax,deltay,deltaz
    real(8) :: v11,v22,v33
    real(8) :: v21,v31,v32
    real(8) :: dtin, rdtbin
    integer(4) :: l0,l1,iconstraints
    integer(4) :: iflgloop
    integer(4) :: iterchk,iam

    iconstraints=ibseL(l0)

    do icycle=1,maxshakecycle
       iflgloop=0
       do l1=1,iconstraints
          ida = shkijL(1,l1,l0)
          idb = shkijL(2,l1,l0)
          dx1 = wkxyz(1,idb)-wkxyz(1,ida)
          dy1 = wkxyz(2,idb)-wkxyz(2,ida)
          dz1 = wkxyz(3,idb)-wkxyz(3,ida)
          dvx = wkv(1,idb)-wkv(1,ida)
          dvy = wkv(2,idb)-wkv(2,ida)
          dvz = wkv(3,idb)-wkv(3,ida)
          spro = dx1*dvx+dy1*dvy+dz1*dvz  ! unit: m*m/s
          tmp2 = -spro/dij2L(l1,l0)   !G   ! unit: 1/s
          tmp = tmp2*dtin            ! unit: none
          if (abs(tmp) .gt. shktol) then
             gmi = rmassL(1,l1,l0)
             gmj = rmassL(2,l1,l0)
             tmp = -tmp2/(gmi+gmj)
             lambL(l1,l0)=lambL(l1,l0)+tmp ! G
             deltax = dx1*tmp
             deltay = dy1*tmp
             deltaz = dz1*tmp
             wkv(1,ida) = wkv(1,ida)+deltax*gmi
             wkv(2,ida) = wkv(2,ida)+deltay*gmi
             wkv(3,ida) = wkv(3,ida)+deltaz*gmi
             wkv(1,idb) = wkv(1,idb)-deltax*gmj
             wkv(2,idb) = wkv(2,idb)-deltay*gmj
             wkv(3,idb) = wkv(3,idb)-deltaz*gmj
             deltax = deltax*rdtbin
             deltay = deltay*rdtbin
             deltaz = deltaz*rdtbin
             wk_dfc(1,ida)=wk_dfc(1,ida)+deltax
             wk_dfc(2,ida)=wk_dfc(2,ida)+deltay
             wk_dfc(3,ida)=wk_dfc(3,ida)+deltaz
             wk_dfc(1,idb)=wk_dfc(1,idb)-deltax
             wk_dfc(2,idb)=wk_dfc(2,idb)-deltay
             wk_dfc(3,idb)=wk_dfc(3,idb)-deltaz
             v11 = dx1*deltax
             v22 = dy1*deltay
             v33 = dz1*deltaz
             v21 = dy1*deltax
             v31 = dz1*deltax
             v32 = dz1*deltay
             if(ida.lt.idb)then  ! avoid double count
                wkdvir2(1,iam)=wkdvir2(1,iam)+v11
                wkdvir2(2,iam)=wkdvir2(2,iam)+v22
                wkdvir2(3,iam)=wkdvir2(3,iam)+v33
                wkdvir2(4,iam)=wkdvir2(4,iam)+v21
                wkdvir2(5,iam)=wkdvir2(5,iam)+v31
                wkdvir2(6,iam)=wkdvir2(6,iam)+v32
             endif
          else
             iflgloop=iflgloop+1
          end if
       end do ! l1
       if(iflgloop.eq.iconstraints) then
          exit
       endif
    end do ! icycle
    if(icycle.eq.maxshakecycle)then
       iterchk = iterchk + 1
    endif
  end subroutine rattle_roll_main
!----------------------------------------------------------------------
!
!     prattle/roll main
!
!----------------------------------------------------------------------
!>
!! \brief PRATTLE: The idea of PSHAKE method is applied to RATTLE.
!<
  subroutine prattle_roll_main(l0,dtin,rdtbin,iterchk,shktol,i_ps,iam)
    use trajectory_mpi
    implicit none
    real(8) :: dtin, rdtbin
    integer(4) :: i,ida,idb, j,l0,l1,iconstraints
    integer(4) :: icycle, iam
    real(8) :: dvx,dvy,dvz
    real(8) :: shktol,spro
    real(8) :: gmi,gmj
    real(8) :: tmp, tmp2
    real(8) ::deltax(n_atom_max),deltay(n_atom_max),deltaz(n_atom_max)
    real(8) :: v11,v22,v33
    real(8) :: v21,v31,v32
    integer(4) :: iflgloop
    integer(4) :: iterchk
    integer(4) :: i_ps, n_atom
    real(8) ::dx1_ps(n_cnst_max),dy1_ps(n_cnst_max),dz1_ps(n_cnst_max)
    real(8) :: tmp2_ps(n_cnst_max)
    real(8) :: v_mat(3,n_cnst_max,n_atom_max)
    real(8) :: v_mat_old(3,n_cnst_max,n_atom_max)
    integer(4) :: atom_ps(n_cnst_max), l, k
    integer(4) :: couple_wk(2,n_cnst_max)
    real(8) :: ma
    real(8) :: mass_wk(n_atom_max)

    iconstraints=ibseL(l0)

    n_atom = n_atom_ps(i_ps)
    atom_ps(1) = shkijL(1,1,l0)
    v_mat_old(:,:iconstraints,:n_atom) = 0.0d0

    do l1=1, iconstraints
       couple_wk(1,l1) = couple_ps(1,l1,i_ps)
       couple_wk(2,l1) = couple_ps(2,l1,i_ps)
    enddo

    do l1 = 2, n_atom
       atom_ps(l1) = atom_ps(1) + l1 - 1
    enddo

    do i=1, n_atom
       mass_wk(i) = mass_ps(i,i_ps)
    enddo

    do l1=1, iconstraints !n_pshake = iconstraints
       ida = shkijL(1,l1,l0)
       idb = shkijL(2,l1,l0)
       gmi = rmassL(1,l1,l0)
       gmj = rmassL(2,l1,l0)
       dx1_ps(l1) = wkxyz(1,idb) - wkxyz(1,ida)
       dy1_ps(l1) = wkxyz(2,idb) - wkxyz(2,ida)
       dz1_ps(l1) = wkxyz(3,idb) - wkxyz(3,ida)
       v_mat_old(1,l1,couple_wk(1,l1)) = dx1_ps(l1)*gmi
       v_mat_old(2,l1,couple_wk(1,l1)) = dy1_ps(l1)*gmi
       v_mat_old(3,l1,couple_wk(1,l1)) = dz1_ps(l1)*gmi
       v_mat_old(1,l1,couple_wk(2,l1)) =-dx1_ps(l1)*gmj
       v_mat_old(2,l1,couple_wk(2,l1)) =-dy1_ps(l1)*gmj
       v_mat_old(3,l1,couple_wk(2,l1)) =-dz1_ps(l1)*gmj
    enddo

!ocl unroll(2)
    do j=1, n_atom
       v_mat(:,:iconstraints,j) = 0.0d0
       do i=1, iconstraints
          do l=1, iconstraints
             do k=1, 3
                v_mat(k,i,j) = v_mat(k,i,j) + v_mat_old(k,l,j) * a_0_sym(l,i,i_ps)
             enddo
          enddo
       enddo
    enddo

    !
    !       apply the rattle algorithm to correct the velocities
    !
    do icycle=1,maxshakecycle
       iflgloop=0
!ocl noswp
       do l1 = 1, iconstraints
          ida = shkijL(1,l1,l0)
          idb = shkijL(2,l1,l0)
          dvx = wkv(1,idb) - wkv(1,ida)
          dvy = wkv(2,idb) - wkv(2,ida)
          dvz = wkv(3,idb) - wkv(3,ida)
          spro = dx1_ps(l1)*dvx + dy1_ps(l1)*dvy + dz1_ps(l1)*dvz
          tmp2_ps(l1) = spro * rdij2_ps(l1,i_ps) ! unit: 1/s
          lambL(l1,l0)=lambL(l1,l0)+tmp2_ps(l1)
          tmp = spro*dtin              ! unit: none
          if (abs(tmp) .le. dij2L(l1,l0)*shktol) then ! absolute value
             iflgloop=iflgloop+1
          endif
       enddo

       if(iflgloop.eq.iconstraints) then
          goto 123
       endif

       do i = 1, n_atom
          deltax(i) = 0.0d0
          deltay(i) = 0.0d0
          deltaz(i) = 0.0d0
          do l1=1, iconstraints
             deltax(i) = deltax(i) + v_mat(1,l1,i) * tmp2_ps(l1)
             deltay(i) = deltay(i) + v_mat(2,l1,i) * tmp2_ps(l1)
             deltaz(i) = deltaz(i) + v_mat(3,l1,i) * tmp2_ps(l1)
          enddo
       enddo
!ocl norecurrence(wkv,wk_dfc)
       do i = 1, n_atom
          ida = atom_ps(i)
          ma = mass_wk(i)
          wkv(1,ida) =  wkv(1,ida) + deltax(i)
          wkv(2,ida) =  wkv(2,ida) + deltay(i)
          wkv(3,ida) =  wkv(3,ida) + deltaz(i)
          deltax(i) = deltax(i) * rdtbin * ma
          deltay(i) = deltay(i) * rdtbin * ma
          deltaz(i) = deltaz(i) * rdtbin * ma
          wk_dfc(1,ida) = wk_dfc(1,ida) + deltax(i)
          wk_dfc(2,ida) = wk_dfc(2,ida) + deltay(i)
          wk_dfc(3,ida) = wk_dfc(3,ida) + deltaz(i)
          v11 = wkxyz(1,ida) * deltax(i)
          v22 = wkxyz(2,ida) * deltay(i)
          v33 = wkxyz(3,ida) * deltaz(i)
          wkdvir2(1,iam) = wkdvir2(1,iam) + v11
          wkdvir2(2,iam) = wkdvir2(2,iam) + v22
          wkdvir2(3,iam) = wkdvir2(3,iam) + v33
          v21 = wkxyz(2,ida) * deltax(i)
          v31 = wkxyz(3,ida) * deltax(i)
          v32 = wkxyz(3,ida) * deltay(i)
          wkdvir2(4,iam) = wkdvir2(4,iam) + v21
          wkdvir2(5,iam) = wkdvir2(5,iam) + v31
          wkdvir2(6,iam) = wkdvir2(6,iam) + v32
       enddo

    end do
123 continue
    if(icycle.eq.maxshakecycle)then
       iterchk = iterchk + 1
    endif

  end subroutine prattle_roll_main
!----------------------------------------------------------------------
!
!     pshake initialize1    2012/10/23
!
!----------------------------------------------------------------------
  subroutine pshake_initialize1
    use atom_mass
    use pshake_init
    use param
    use trajectory_mpi
    use trajectory_org
    implicit none
    integer(4):: num, GroupNumber, cnst_temp, atom_temp
    integer(4):: ia, ib, ia_max, ia_min
    integer(4):: i, j, k
    integer(4),allocatable:: n_cnst_temp(:), n_atom_temp(:)
    integer(4),allocatable:: couple_ps_temp(:,:,:)
    real(8),allocatable:: r_init_temp(:,:,:), mass_ps_temp(:,:)
    integer(4):: count_type_ps, i_atom, i_cnst, i_type
    integer(4):: nl_cnst_max, nl_atom_max
    logical,allocatable:: ps_flag(:)
    include 'mpif.h'


    n_type_ps = 0
    if(totnconst==0) return  !whole this subroutine

    if(myrank==0)then
       !       write(*,*)
       !       write(*,*)'##### PSHAKE info. #####'
    endif

    allocate(gn2ips(rngrp_ps))
    gn2ips(:) = 0

    allocate(ps_flag(rngrp_ps))
    ps_flag(:) = .true.

    allocate(n_cnst_temp(rngrp_ps))
    allocate(n_atom_temp(rngrp_ps))
    allocate(couple_ps_temp(2,n_cnst_max,rngrp_ps))
    allocate(r_init_temp(3,n_atom_max,rngrp_ps))
    allocate(mass_ps_temp(n_atom_max,rngrp_ps))
    n_cnst_temp(:) = 0
    n_atom_temp(:) = 0
    couple_ps_temp(:,:,:) = 0
    r_init_temp(:,:,:) = 0.0d0
    mass_ps_temp(:,:) = 0.0d0

    nl_cnst_max = 0
    nl_atom_max = 0
    count_type_ps = 0
    do i=1, n
       num = paranum(i)
       GroupNumber = ShakeGroupLeader(num)
       if(GroupNumber == 0)cycle

       if(ps_flag(GroupNumber))then
          ps_flag(GroupNumber) = .false.
          j = count_type_ps + 1

          !get the number of constraints in the shake group
          cnst_temp = nconstraints(GroupNumber)
          !not rigid
          if(cnst_temp .gt. n_cnst_max)then
             if(myrank==0)then
                write(0,*) 'ERROR: cnst_temp .gt. n_cnst_max at GroupNumber=', GroupNumber
                write(0,*) 'n_cnst_max=', n_cnst_max
                write(0,*) 'cnst_temp =', cnst_temp
                write(0,*)
             endif
             call modylas_abort() ! for update_shake_local
          endif
          !! less than 3 constraints
          !if(cnst_temp < 3)then
          !  goto 100
          !endif

          !get the number of atoms in the shake group
          ia_max = 0
          ia_min = 100000000
          do i_cnst= 1, cnst_temp
             ia=atom1S(GroupNumber,i_cnst)
             ib=atom2S(GroupNumber,i_cnst)
             if(ia .lt. ia_min)then
                ia_min = ia
             endif
             if(ib .gt. ia_max)then
                ia_max = ib
             endif
          enddo
          atom_temp = ia_max - ia_min + 1
          if(atom_temp .gt. n_atom_max)then
             if(myrank==0)then
                write(*,*)'WARN: atom_temp .gt. n_atom_max at GroupNumber=', GroupNumber
                write(*,*)'n_atom_max=', n_atom_max
                write(*,*)'atom_temp =', atom_temp
                write(*,*)
             endif
             goto 100
             !stop
          endif
          !!not rigid
          !if(atom_temp == 3 .and. cnst_temp /= 3)then
          !  goto 100
          !else if(atom_temp == 4 .and. cnst_temp /= 6)then
          !  goto 100
          !endif

          !temporary assignment for P-SHAKE
          !get the numbers of atom & constraint
          n_cnst_temp(j) = cnst_temp
          n_atom_temp(j) = atom_temp
          if(atom_temp > nl_atom_max)then
             nl_atom_max = atom_temp
          endif
          if(cnst_temp > nl_cnst_max)then
             nl_cnst_max = cnst_temp
          endif

          !get the constraint couples in the SHAKE group
          do i_cnst=1, n_cnst_temp(j)
             couple_ps_temp(1,i_cnst,j) = atom1S(GroupNumber,i_cnst) - ia_min + 1
             couple_ps_temp(2,i_cnst,j) = atom2S(GroupNumber,i_cnst) - ia_min + 1
          enddo

          !get mass of atoms in the SHAKE group
          do i_atom=1, n_atom_temp(j) !check
             ia = i + i_atom - 1
             !           write(*,*)ia, paranum(ia)
             mass_ps_temp(i_atom,j) = mass(paranum(ia))
             r_init_temp(1,i_atom,j) = xyz(1,ia)
             r_init_temp(2,i_atom,j) = xyz(2,ia)
             r_init_temp(3,i_atom,j) = xyz(3,ia)
          enddo

          count_type_ps = count_type_ps + 1
          gn2ips(GroupNumber) = count_type_ps
       endif

100    continue

    enddo

    n_type_ps = count_type_ps

    if(n_type_ps==0)then
       !       if(myrank==0)then
       !         write(*,*)'PSHAKE has NOT been applied.'
       !         write(*,*)
       !       endif
       return
    end if

    if(myrank==0)then
       !       write(*,*)'Pshake appllied:', n_type_ps, '/', rngrp_ps
       !       write(*,*)
    endif

    !assign the P-SHAKE informations

    allocate(n_cnst_ps(n_type_ps))
    allocate(n_atom_ps(n_type_ps))
    allocate(couple_ps(2,nl_cnst_max,n_type_ps))
    allocate(r_init_ps(3,nl_atom_max,n_type_ps))
    allocate(mass_ps(nl_atom_max,n_type_ps))
    n_cnst_ps(:)=0
    n_atom_ps(:)=0
    couple_ps(:,:,:)=0
    r_init_ps(:,:,:)=0.0d0
    mass_ps(:,:)=0.0d0

    allocate(a_0(nl_cnst_max,nl_cnst_max,n_type_ps))
    allocate(a_0_sym(nl_cnst_max,nl_cnst_max,n_type_ps))
    allocate(rdij2_ps(nl_cnst_max,n_type_ps))
    a_0(:,:,:)=0.0d0
    a_0_sym(:,:,:)=0.0d0
    rdij2_ps(:,:)=0.0d0

    do i_type=1, n_type_ps
       n_cnst_ps(i_type) = n_cnst_temp(i_type)
       n_atom_ps(i_type) = n_atom_temp(i_type)
       do i_cnst=1, n_cnst_temp(i_type)
          couple_ps(1,i_cnst,i_type) = couple_ps_temp(1,i_cnst,i_type)
          couple_ps(2,i_cnst,i_type) = couple_ps_temp(2,i_cnst,i_type)
       enddo
       do i_atom=1, n_atom_temp(i_type)
          do k=1, 3
             r_init_ps(k,i_atom,i_type) = r_init_temp(k,i_atom,i_type)
          enddo
       enddo
       do i_atom=1, n_atom_temp(i_type)
          mass_ps(i_atom,i_type) = mass_ps_temp(i_atom,i_type)
       enddo
    enddo

    deallocate(ps_flag)
    deallocate(n_cnst_temp)
    deallocate(n_atom_temp)
    deallocate(couple_ps_temp)
    deallocate(r_init_temp)
    deallocate(mass_ps_temp)
  end subroutine pshake_initialize1
!----------------------------------------------------------------------
!
!     finish p-shake 1
!
!----------------------------------------------------------------------
  subroutine pshake_finish1
    use pshake_init
    implicit none

    if(n_type_ps==0) return

    deallocate(n_cnst_ps)
    deallocate(r_init_ps)
  end subroutine pshake_finish1
!----------------------------------------------------------------------
!
!     finish p-shake 2
!
!----------------------------------------------------------------------
  subroutine pshake_finish2
    implicit none

    if(n_type_ps==0) return

    deallocate(couple_ps)
    deallocate(n_atom_ps)
    deallocate(a_0)
    deallocate(a_0_sym)
    deallocate(type_psL)
    deallocate(rdij2_ps)
    deallocate(gn2ips)
    deallocate(mass_ps)
  end subroutine pshake_finish2
!----------------------------------------------------------------------
!
!     initialize p-shake 2
!
!----------------------------------------------------------------------
!>
!! \brief Subroutine which pre-calculate matrix of constraint group for decoupling the constraint equations.
!<
  subroutine pshake_initialize2
    use pshake_init
    implicit none
    real(8),allocatable:: mal(:,:), mrv(:,:)
    integer(4),parameter:: n_piter_max=10
    integer(4):: k, i, j, l, n_piter, i_cnst, i_atom
    integer(4):: ia, ib, ja, jb, i_count, j_count
    real(8):: dr(3,n_cnst_max), dv(3,n_cnst_max,n_cnst_max)
    real(8):: dsdr(3,n_cnst_max,n_atom_max)
    real(8):: r_mat(n_cnst_max-1,n_cnst_max-1), r_clm(n_cnst_max-1)
    real(8):: r_mat_all(n_cnst_max,n_cnst_max)
    real(8):: r_mat_inv(n_cnst_max-1,n_cnst_max-1)
    real(8):: v_pshake(3,n_cnst_max,n_atom_max)
    real(8):: v_pshake_new(3,n_cnst_max,n_atom_max)
    real(8):: v_mat_old(3,n_cnst_max,n_atom_max)
    real(8):: v_mat(3,n_cnst_max,n_atom_max)
    real(8):: a_mat(n_cnst_max,n_cnst_max)
    real(8):: a_0_new(n_cnst_max,n_cnst_max)
    real(8):: a_0_wk(n_cnst_max,n_cnst_max)
    real(8):: gmi, gmj, r_mass_ps(n_atom_max)
    integer(4):: i_ps
    integer(4):: n_cnst_wk, n_atom_wk
    integer(4):: couple_wk(2,n_cnst_max)
    real(8):: r_init_wk(3,n_atom_max)
    real(8):: mass_wk(n_atom_max), len_cnst_wk(n_cnst_max)
    real(8):: temp

    do i_ps=1, n_type_ps  !whole subroutine

       !initialize
       couple_wk(:,:) = 0
       r_init_wk(:,:) = 0.0d0
       mass_wk(:) = 0.0d0
       len_cnst_wk(:) = 0.0d0

       n_cnst_wk = n_cnst_ps(i_ps)
       n_atom_wk = n_atom_ps(i_ps)
       !write(*,*)n_cnst_wk, n_atom_wk
       do i_cnst=1, n_cnst_wk
          couple_wk(1,i_cnst) = couple_ps(1,i_cnst,i_ps)
          couple_wk(2,i_cnst) = couple_ps(2,i_cnst,i_ps)
          !       write(*,*)'couple:',couple_wk(1,i_cnst), couple_wk(2,i_cnst)
       enddo
       do i_atom=1, n_atom_wk
          do k=1, 3
             r_init_wk(k,i_atom) = r_init_ps(k,i_atom,i_ps)
          enddo
          !write(*,*)'r:',(r_init_wk(i,i_atom),i=1,3)
       enddo
       do i_atom=1, n_atom_wk
          mass_wk(i_atom) = mass_ps(i_atom,i_ps)
          !write(*,*)'mass:',mass_wk(i_atom)
       enddo

       !check
       do i=1, n_atom_wk
          !mass_pshake(i) = mass_pshake(i) * md_ATOMIC_MASS_UNIT
          r_mass_ps(i) = 1.0d0 / mass_wk(i)
       enddo

       if(n_cnst_wk > 1)then !multi constraints in shake group

          allocate(mal(n_cnst_wk-1,(n_cnst_wk-1)*2))
          allocate(mrv(n_cnst_wk-1,n_cnst_wk-1))
          mal(:,:) = 0.0d0
          mrv(:,:) = 0.0d0

          !check
          !dt2 = dt*dt

          !calc eq. 6
          do i=1, n_cnst_wk
             ia = couple_wk(1,i)
             ib = couple_wk(2,i)
             !write(*,*)ia, ib
             do k=1, 3
                dr(k,i) = r_init_wk(k,ia) - r_init_wk(k,ib)
             enddo
             !temp = dr(1,i)*dr(1,i) + dr(2,i)*dr(2,i) + dr(3,i)*dr(3,i)
             !temp = dsqrt(temp)
             !write(*,*)temp, dr(1,i), dr(2,i), dr(3,i)
          enddo
          !write(*,*)

          dsdr(:,:,:) = 0.0d0
          do i=1, n_cnst_wk
             ia = couple_wk(1,i)
             ib = couple_wk(2,i)
             do k=1, 3
                dsdr(k,i,ia) = 2.0d0*dr(k,i)
                dsdr(k,i,ib) = -1.0d0*dsdr(k,i,ia)
             enddo
          enddo

          do j=1, n_atom_wk
             do i=1, n_cnst_wk
                do k=1, 3
                   v_pshake(k,i,j) = dsdr(k,i,j)*r_mass_ps(j)
                   !v_pshake(k,i,j) = dsdr(k,i,j)*dt2*r_mass_pshake(j)
                   !!There is (dt)^2 at every matrix element.
                enddo
                !write(*,*)v_pshake(1,i,j), v_pshake(2,i,j), v_pshake(3,i,j)
             enddo
          enddo

          !initialize A0
          a_0_wk(:,:) = 0.0d0
          do i=1, n_cnst_wk
             a_0_wk(i,i) = 1.0d0
          enddo

          do n_piter=1, n_piter_max

             !calc eq. 7
             dv(:,:,:) = 0.0d0
             do j=1, n_cnst_wk  !Yoko
                ja = couple_wk(1,j)
                jb = couple_wk(2,j)
                do i=1, n_cnst_wk  !tate
                   do k=1, 3
                      dv(k,i,j) = v_pshake(k,i,ja) - v_pshake(k,i,jb)
                   enddo
                   !write(*,*)i, j,dv(1,i,j), dv(2,i,j), dv(3,i,j)
                enddo
             enddo
             !write(*,*)

             !calc r matrix (Nc x Nc)
             do j=1, n_cnst_wk  !tate
                do i=1, n_cnst_wk  !Yoko
                   r_mat_all(i,j) = dv(1,i,j)*dr(1,j) + dv(2,i,j)*dr(2,j) + dv(3,i,j)*dr(3,j)
                enddo
                !write(*,*)i, r_mat_all(i,1), r_mat_all(i,2)
             enddo

             a_mat(:,:) = 0.0d0
             do i=1, n_cnst_wk
                a_mat(i,i) = 1.0d0
             enddo

          !iteration for A0
             do k=1, n_cnst_wk

                !calc r matrix (Nc-1 x Nc-1)
                j_count = 1
                do j=1, n_cnst_wk !Yoko
                   if(j == k)cycle
                   i_count = 1
                   do i=1, n_cnst_wk !tate
                      if(i == k)cycle
                      r_mat(i_count,j_count) = r_mat_all(i,j)
                      !write(*,*)'a', i_count, j_count, r_mat(i_count,j_count)
                      i_count = i_count + 1
                   enddo
                   !write(*,*)r_mat(i_count,1), r_mat(i_count,2)
                   j_count = j_count + 1
                enddo
             !write(*,*)

                !calc r column (Nc-1 x 1) in eq. 29
                j_count = 1
                do j=1, n_cnst_wk
                   if(j == k)cycle
                   r_clm(j_count) = -r_mat_all(k,j)
                   !write(*,*)'b', j_count, r_clm(j_count)
                   j_count = j_count + 1
                enddo

                !set for matrix inversino
                do j=1, n_cnst_wk-1
                   do i=1, n_cnst_wk-1
                      mal(i,j) = r_mat(i,j)
                   enddo
                enddo
                do i=1, n_cnst_wk-1
                   mal(i,i+n_cnst_wk-1) = 1.0d0
                enddo

                !matrix inversion:
                !mal --> mrv
                call inv_matrix(n_cnst_wk, mal, mrv)

                !assignment of the inversion matrix
                do j=1, n_cnst_wk-1
                   do i=1, n_cnst_wk-1
                      r_mat_inv(i,j) = mrv(i,j)
                   enddo
                enddo

                !calculate the elements in A
                i_count=1
                do i=1, n_cnst_wk ! tate for A
                   if(i == k)then
                      a_mat(k,i) = 1.0d0
                      cycle
                   endif

                   do j=1, n_cnst_wk-1
                      a_mat(k,i) = a_mat(k,i) + r_mat_inv(j,i_count)*r_clm(j)
                      !write(*,*)i, j, r_mat_inv(j,i_count), a_mat(k,i)
                   enddo

                   i_count = i_count + 1
                enddo

             enddo

             !"V <- VA"
             v_pshake_new(:,:,:) = 0.0d0
             a_0_new(:,:) = 0.0d0
             do j=1, n_atom_wk  !tate
                do i=1, n_cnst_wk !Yoko
                   do l=1, n_cnst_wk
                      do k=1, 3
                         v_pshake_new(k,i,j) = v_pshake_new(k,i,j) + v_pshake(k,l,j)*a_mat(i,l)
                      enddo
                   enddo
                enddo
             enddo

             do j=1, n_atom_wk
                do i=1, n_cnst_wk
                   do k=1, 3
                      v_pshake(k,i,j) = v_pshake_new(k,i,j)
                   enddo
                enddo
             enddo

             !"(A_0) <- (A_0)A"
             do j=1, n_cnst_wk
                do i=1, n_cnst_wk
                   do l=1, n_cnst_wk
                      a_0_new(i,j) = a_0_new(i,j) + a_0_wk(l,j)*a_mat(i,l)
                   enddo
                enddo
             enddo

             do j=1, n_cnst_wk
                do i=1, n_cnst_wk
                   a_0_wk(i,j) = a_0_new(i,j)
                enddo
             enddo


          enddo

          !set transposed matrix of a_0 for cost cut
          do j=1, n_cnst_wk
             do i=1, n_cnst_wk
                a_0(i,j,i_ps) = a_0_wk(i,j)
                a_0_sym(j,i,i_ps) = a_0_wk(i,j)
             enddo
          enddo

          deallocate(mal)
          deallocate(mrv)

       else ! 1 constraints in shake group

          do j=1, n_cnst_wk
             do i=1, n_cnst_wk
                a_0(i,j,i_ps) = 1.0d0
                a_0_sym(j,i,i_ps) = 1.0d0
             enddo
          enddo

       endif

       v_mat_old(:,:,:)=0.0d0
       do i_cnst=1, n_cnst_wk
          ia = couple_wk(1,i_cnst)
          ib = couple_wk(2,i_cnst)
          do k=1, 3
             dr(k,i_cnst) = r_init_wk(k,ia) - r_init_wk(k,ib)
          enddo
          gmi= r_mass_ps(ia)
          gmj= r_mass_ps(ib)
          v_mat_old(1,i_cnst,ia) = dr(1,i_cnst)*gmi
          v_mat_old(2,i_cnst,ia) = dr(2,i_cnst)*gmi
          v_mat_old(3,i_cnst,ia) = dr(3,i_cnst)*gmi
          v_mat_old(1,i_cnst,ib) =-dr(1,i_cnst)*gmj
          v_mat_old(2,i_cnst,ib) =-dr(2,i_cnst)*gmj
          v_mat_old(3,i_cnst,ib) =-dr(3,i_cnst)*gmj
       enddo

       v_mat(:,:,:) = 0.0d0
       do i=1, n_cnst_wk
          do j=1, n_atom_wk
             do l=1, n_cnst_wk
                do k=1, 3
                   v_mat(k,i,j) = v_mat(k,i,j) + v_mat_old(k,l,j) * a_0_sym(l,i,i_ps)
                enddo
             enddo
          enddo
       enddo

       do i_cnst=1, n_cnst_wk
          ia = couple_wk(1,i_cnst)
          ib = couple_wk(2,i_cnst)
          temp = (v_mat(1,i_cnst,ia)-v_mat(1,i_cnst,ib))*dr(1,i_cnst) + &
               &         (v_mat(2,i_cnst,ia)-v_mat(2,i_cnst,ib))*dr(2,i_cnst) + &
               &         (v_mat(3,i_cnst,ia)-v_mat(3,i_cnst,ib))*dr(3,i_cnst)
          rdij2_ps(i_cnst,i_ps) = 1.0d0/temp
       enddo

    enddo

  end subroutine pshake_initialize2

!----------------------------------------------------------------------
  subroutine inv_matrix(n_cnst, mal, mrv)
!----------------------------------------------------------------------
    implicit none
    integer(4), intent(in) :: n_cnst
    real(8), intent(inout) :: mal(:,:)
    real(8), intent(out) :: mrv(:,:)
    integer(4) :: i, j, k, n
    real(8) :: value
    real(8) :: tmp, aki

    ! for getting scaling parameter

    n = n_cnst-1

    value = 0.0d0
    do j=1, n
       do i=1, n
          value = value + dabs(mal(i,j))
       enddo
    enddo

    value = value / dble(n*n)

    ! for making upper triangular matrix
    do i=1, n-1

       ! exchange when diagonal element is zero
       if(dabs(mal(i,i)/value) < 1.0d-6)then
          do k=i+1, n

             if(dabs(mal(k,i)/value) > 1.0d-6)then
                do j=1, n*2
                   tmp = mal(k,j)
                   mal(k,j) = mal(i,j)
                   mal(i,j) = tmp
                enddo
                exit
             endif

          enddo
       endif

       ! for making diagonal element 1
       tmp = 1.0d0 / mal(i,i)
       do j=i, n*2
          mal(i,j) = mal(i,j) * tmp
       enddo

       ! for making element under diagonal 0
       do k=i+1, n

          if(dabs(mal(k,i)/value) > 1.0d-6)then
             aki = mal(k,i)
             do j=i, n*2
                mal(k,j) = mal(k,j) - mal(i,j)*aki
             enddo
          endif

       enddo

       !          do t1=1, n
       !            write(*,'(4es12.4)')(mal(t1,t2),t2=1,n*2)
       !          enddo

    enddo

    ! for making diagonal element 1 at 'n' th line
    tmp = 1.0d0 / mal(n,n)
    do j=n, n*2
       mal(i,j) = mal(i,j) * tmp
    enddo

    ! for making diagonal matrix
    do i=n, 2, -1
       do k=1, i-1
          aki = mal(k,i)
          do j=i, n*2
             mal(k,j) = mal(k,j) - mal(i,j)*aki
          enddo
       enddo
    enddo

    do j=1, n
       do i=1, n
          mrv(i,j) = mal(i,j+n)
       enddo
    enddo

  end subroutine inv_matrix
!----------------------------------------------------------------------
!>
!! \brief Subroutine to initialize SHAKE arrays (localized)
!! \author Yoshimichi Andoh
!<
  subroutine init_shake_local()
!----------------------------------------------------------------------
    use domain, only : lxdiv, lydiv, lzdiv, ncellx, ncelly, ncellz
    !pshake
    implicit none

    !     maxibn=5
    !     maxibn=10
    totnconstL=max( int(totnconst/(ncellx*ncelly*ncellz)*lxdiv*lydiv*lzdiv*2.5d0), &
                    totnconstL_init)
    allocate(ibseL(totnconstL))
    allocate(shkijL(2,n_cnst_max,totnconstL))
    allocate(rmassL(2,n_cnst_max,totnconstL))
    allocate(lambL(n_cnst_max,totnconstL))
    allocate(dij2L(n_cnst_max,totnconstL))
    allocate(igshakeL(totnconstL))
    allocate(lambLstr(n_cnst_max,totnconstL))
    allocate(type_psL(totnconstL))

  end subroutine init_shake_local
!----------------------------------------------------------------------
!>
!! \brief Subroutine to update SHAKE arrays (localized)
!! \author Yoshimichi Andoh
!<
  subroutine update_shake_local()
!----------------------------------------------------------------------
    use atom_mass
    use subcell, only : nselfseg, lsegtop, lseg_natoms, i2m, m2i
    use param
    use trajectory_mpi
    implicit none
    include 'mpif.h'
    integer(4) :: k0,i0,l0,l1
    integer(4) :: i,ibn
    integer(4) :: ia,ib,ida,idb
    integer(4) :: num, GroupNumber, dn1, dn2

    !=== copy constants to local arrays ===!
    l0=0
    kshake=0
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          if(kshake(i0).ne.0) cycle
          i = m2i(i0) ! G (global)
          num=paranum(i)
          if(ShakeGroupLeader(num)==0) cycle
          GroupNumber=ShakeGroupLeader(num)
          l0=l0+1
          if(l0.ge.totnconstL)then
             write(0,*) 'ERROR: l0 >= totnconstL'
             write(0,*) l0,totnconstL
             call modylas_abort()
          endif
          ibseL(l0)=nconstraints(GroupNumber)
          type_psL(l0)=gn2ips(GroupNumber) !pshake

          l1=0
          do ibn = 1, nconstraints(GroupNumber)
             l1=l1+1
             if(l1.ge.n_cnst_max) then
                write(0,*) 'ERROR: l1 >= n_cnst_max'
                call modylas_abort()
             endif
             !         l1=l1+1 ; if(l1.ge.maxibn) stop'ERROR: l1.ge.maxibn'
             dn1=abs(atom1S(GroupNumber,ibn)-num)
             dn2=abs(atom2S(GroupNumber,ibn)-num)
             ia=i + dn1
             ib=i + dn2
             ida = i2m(ia)  ! L (local)
             idb = i2m(ib)
             kshake(ida)=+1
             kshake(idb)=+1
             shkijL(1,l1,l0) = ida
             shkijL(2,l1,l0) = idb
             dij2L(l1,l0) = slength(GroupNumber,ibn)**2
             rmassL(1,l1,l0) = r_mass(paranum(ia)) ! G
             rmassL(2,l1,l0) = r_mass(paranum(ib)) ! G
          enddo ! ibn
       enddo ! i0
    enddo ! ii
    l0max=l0

  end subroutine update_shake_local

!----------------------------------------------------------------------
  subroutine fmod_alloc_vstr(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(vstr(3,ivalue))
  end subroutine fmod_alloc_vstr
!----------------------------------------------------------------------
  subroutine fmod_alloc_dij(ivalue)
    implicit none
    integer(4) :: ivalue
    allocate(dij(ivalue))
  end subroutine fmod_alloc_dij
!----------------------------------------------------------------------
  subroutine fmod_alloc_lambstr(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(lambstr(ivalue))
  end subroutine fmod_alloc_lambstr
!----------------------------------------------------------------------
  subroutine fmod_alloc_lamb(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(lamb(ivalue))
    lamb = 0.0d0
  end subroutine fmod_alloc_lamb
!----------------------------------------------------------------------
  subroutine fmod_alloc_shkij(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(shkij(2,ivalue))
  end subroutine fmod_alloc_shkij
!----------------------------------------------------------------------
  subroutine fmod_totnconst(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    totnconst = ivalue
  end subroutine fmod_totnconst
!----------------------------------------------------------------------
  subroutine fmod_alloc_ishake(ivalue)
    integer(4), intent(in) :: ivalue
    allocate(ishake(2,ivalue))
  end subroutine fmod_alloc_ishake
!----------------------------------------------------------------------
  subroutine fmod_alloc_kshake(ivalue)
    integer(4), intent(in) :: ivalue
    allocate(kshake(ivalue))
  end subroutine fmod_alloc_kshake
!----------------------------------------------------------------------
  subroutine fmod_md_shake__shake_tolerance(value)
    implicit none
    real(8), intent(in) :: value
    shake_tolerance = value
  end subroutine fmod_md_shake__shake_tolerance
!----------------------------------------------------------------------
  subroutine fmod_md_shake__rattle_tolerance(value)
    implicit none
    real(8), intent(in) :: value
    rattle_tolerance = value
  end subroutine fmod_md_shake__rattle_tolerance
!----------------------------------------------------------------------
  subroutine fmod_shake_max_iteration(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    maxshakecycle = ivalue
  end subroutine fmod_shake_max_iteration
!----------------------------------------------------------------------
  subroutine fmod_shake_ngshake(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    ngshake = ivalue
  end subroutine fmod_shake_ngshake
!----------------------------------------------------------------------
  subroutine set_shake()
    use parse_shake
    use param
    use pshake_init
    implicit none
    integer(4):: i, j, k
    integer(4):: atom1, atom2
    integer(4) :: ngrp
    type(yyparse_shake), pointer :: p1, p2, tempp1
!#ifdef MPIPARA
    include 'mpif.h'
    integer(4) :: ierr
!#endif /*MPI*/


    allocate(tempp1)
    do i = 1, npara
       !if(myrank==0) write(*,*) "i=", i
       p1 => yyparse_shake_top(i)%next
       do while(.true.)
          if(.not.associated(p1)) exit
          atom1 = p1%atom1
          atom2 = p1%atom2
          p2 => yyparse_shake_top(atom2)
          !call mpi_barrier(MPI_COMM_WORLD,ierr)
          if(.not.associated(p2%next)) then
            !write(*,*) "OK2?"
             !call mpi_barrier(MPI_COMM_WORLD,ierr)
             p1 => p1%next
             cycle
          endif
          p2 => p2%next
          tempp1%next => p1%next
          p1%next => p2
          do while(.true.)
             if(.not. associated(p2%next)) exit
             p2 => p2%next
          enddo
          p2%next=>tempp1%next
          nullify(yyparse_shake_top(atom2)%next)
          p1 => p1%next
       enddo
    enddo


    allocate(ShakeGroupLeader(npara))
    ngrp = 0
    do i = 1, npara
       if(.not.associated(yyparse_shake_top(i)%next)) then
          ShakeGroupLeader(i) = 0
       else
          ngrp = ngrp + 1
          ShakeGroupLeader(i) = ngrp
       endif
    enddo
    rngrp_ps = ngrp !pshake
    !     if(myrank==0) then
    !       write(*,*) "The number of SHAKE GOURP:", ngrp
    !     endif
    call mpi_barrier(MPI_COMM_WORLD,ierr)

    allocate(nconstraints(ngrp))
    allocate(atom1S(ngrp,10), atom2S(ngrp,10))
    allocate(slength(ngrp,10))
    atom1S(:,:) = -1
    atom2S(:,:) = -1
    slength(:,:) = -1.0d0
    ngrp = 0
    do i = 1, npara
       p1 => yyparse_shake_top(i)%next
       if(.not.associated(yyparse_shake_top(i)%next)) cycle
       ngrp = ngrp + 1
       j = 0
       do while(.true.)
          if(.not.associated(p1)) exit
          j = j + 1
          atom1S(ngrp,j) = p1%atom1
          atom2S(ngrp,j) = p1%atom2
          slength(ngrp,j)= p1%r0
          p1 => p1%next
       enddo
       nconstraints(ngrp) = j
    enddo

    deallocate(yyparse_shake_top)

    call mpi_barrier(MPI_COMM_WORLD,ierr)

  end subroutine set_shake
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to check ROLL convergence
!! \author Yoshimichi Andoh
!<
  subroutine check_roll_convergence()
!----------------------------------------------------------------------
    use md_forces
    use atom_virial
    implicit none

    real(8) :: rlamb_thr
    integer(4) :: ierr
    include 'mpif.h'
    integer(4) :: l0,l1,iconstraints

    TIME_START(TM_ROLL_CONV)
    rlamb_thr=1d0/roll_tolerance

!$omp parallel default(shared) &
!$omp& private(l0,l1,iconstraints) &
!$omp& reduction(+:wk_roll_flg)
!$omp do
    do l0=1,l0max
       iconstraints=ibseL(l0)
       do l1=1,iconstraints
          wk_roll_flg=wk_roll_flg+AINT(rlamb_thr*abs(lambLstr(l1,l0)-lambL(l1,l0)))
       enddo
    enddo
!$omp end do
!$omp end parallel
    TIME_STOP(TM_ROLL_CONV)

    TIME_BARRIER(TMB_ALLREDUCE_ROLL)
    TIME_START(TM_ALLREDUCE_ROLL)
    call mpi_allreduce(wk_roll_flg,roll_flg,1, &
         &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    TIME_STOP(TM_ALLREDUCE_ROLL)
  end subroutine check_roll_convergence

!>
!! \brief  Subroutine to apply SHAKE to two mass centers
!! \author Yoshimichi Andoh, originated by K.Fujimoto,H.Kojima
!<
  subroutine shake_com(dtin)
    use boundary
    use center_of_mass
    use center_of_mass_variables
    use trajectory_mpi
    implicit none
    integer(4) :: i, j, i0
    real(8) dx, dy, dz
    real(8) dxo, dyo, dzo
    real(8) :: COMterm(3)
    real(8) :: deltaR, deltaRold
    real(8) :: A, B, C
    real(8) :: tmp, dtin, rdtin !hk
    include 'mpif.h'

    if(constrain_COM)then
       continue
    else
       return
    endif

    if(change_distCOM)then !.true.
       KratCOM  = KratCOM + d_COMR
       KratCOM2 = KratCOM * KratCOM
       !       write(*,*) "KratCOM is updated", KratCOM, d_COMR
    endif

#ifdef ZSHAKE
    !     dxo = 0.0d0
    !     dyo = 0.0d0
    !     dzo = CenterOfMassAOld(3) - CenterOfMassBOld(3)
    !     call PBCij(dxo, dyo, dzo)
    !     dx = 0.0d0
    !     dy = 0.0d0
    !     dz = CenterOfMassA(3) - CenterOfMassB(3)
    !     call PBCij(dxo, dyo, dzo)
    write(0,*) 'ERROR: ZSHAKE not supported'
    call modylas_abort()
#else
    dxo = CenterOfMassAOld(1) - CenterOfMassBOld(1)
    dyo = CenterOfMassAOld(2) - CenterOfMassBOld(2)
    dzo = CenterOfMassAOld(3) - CenterOfMassBOld(3)
    call PBCij(dxo, dyo, dzo)
    dx = CenterOfMassA(1) - CenterOfMassB(1)
    dy = CenterOfMassA(2) - CenterOfMassB(2)
    dz = CenterOfMassA(3) - CenterOfMassB(3)
    call PBCij(dx, dy, dz)
#endif
    deltaR    = dx**2 + dy**2 + dz**2
    deltaRold = dxo**2 + dyo**2 + dzo**2

    A = mu**2 * deltaRold
    B = 2.0d0*mu*(dx*dxo + dy*dyo + dz*dzo)
    C = deltaR - KratCOM2
    sgamma = (-B + dsqrt(B**2-4.0d0*A*C))/(2.0d0*A)

    COMterm(1) = sgamma * dxo
    COMterm(2) = sgamma * dyo
    COMterm(3) = sgamma * dzo

    !position
    !集団Aの原子の位置更新
    do i=1,i0ownA
       i0=i0listA(i)
       wkxyz(1,i0) = wkxyz(1,i0) + COMterm(1)*RTotalMassA
       wkxyz(2,i0) = wkxyz(2,i0) + COMterm(2)*RTotalMassA
       wkxyz(3,i0) = wkxyz(3,i0) + COMterm(3)*RTotalMassA
    enddo

    !集団Bの原子の位置更新
    do i=1,i0ownB
       i0=i0listB(i)
       wkxyz(1,i0) = wkxyz(1,i0) - COMterm(1)*RTotalMassB
       wkxyz(2,i0) = wkxyz(2,i0) - COMterm(2)*RTotalMassB
       wkxyz(3,i0) = wkxyz(3,i0) - COMterm(3)*RTotalMassB
    enddo

    !velocity
    rdtin = 1.0d0/dtin
    COMterm(1) = COMterm(1) * rdtin
    COMterm(2) = COMterm(2) * rdtin
    COMterm(3) = COMterm(3) * rdtin

    !集団Aの原子の速度
    do i=1,i0ownA
       i0=i0listA(i)
       wkv(1,i0) = wkv(1,i0) + COMterm(1)*RTotalMassA
       wkv(2,i0) = wkv(2,i0) + COMterm(2)*RTotalMassA
       wkv(3,i0) = wkv(3,i0) + COMterm(3)*RTotalMassA
    enddo

    !集団Bの原子の速度更新
    do i=1,i0ownB
       i0=i0listB(i)
       wkv(1,i0) = wkv(1,i0) - COMterm(1)*RTotalMassB
       wkv(2,i0) = wkv(2,i0) - COMterm(2)*RTotalMassB
       wkv(3,i0) = wkv(3,i0) - COMterm(3)*RTotalMassB
    enddo
  end subroutine shake_com
!====================================================================
!>
!! \brief  Subroutine to apply RATTLE to two mass centers
!! \author Yoshimichi Andoh, originated by K.Fujimoto,H.Kojima
!<
  subroutine rattle_com(dtin)
    use boundary
    use center_of_mass
    use center_of_mass_variables
    use trajectory_mpi
    implicit none
    integer(4) :: i, j, i0
    real(8) :: dx, dy, dz
    real(8) :: dxo, dyo, dzo
    real(8) :: dvx, dvy, dvz
    real(8) :: dvdr, dr0dr, drdr
    real(8) :: A, B, term(3)
    real(8) :: dtin, rdtin
    real(8) :: tmp2 !hk
    real(8) :: rchk_old, rchk_new !hk
    include 'mpif.h'

    if(constrain_COM)then
       continue
    else
       return
    endif

    rdtin = 1.0d0/dtin
#ifdef ZSHAKE
    !     dx = 0.0d0
    !     dy = 0.0d0
    !     dz = CenterOfMassA(3) - CenterOfMassB(3)
    !     call PBCij(dx, dy, dz)
    !     dvx = 0.0d0
    !     dvy = 0.0d0
    !     dvz = CenterOfMassVeloA(3) - CenterOfMassVeloB(3)
#else
    dx = CenterOfMassA(1) - CenterOfMassB(1)
    dy = CenterOfMassA(2) - CenterOfMassB(2)
    dz = CenterOfMassA(3) - CenterOfMassB(3)
    call PBCij(dx, dy, dz)
    dvx = CenterOfMassVeloA(1) - CenterOfMassVeloB(1)
    dvy = CenterOfMassVeloA(2) - CenterOfMassVeloB(2)
    dvz = CenterOfMassVeloA(3) - CenterOfMassVeloB(3)
#endif

    dvdr = dx*dvx + dy*dvy + dz*dvz
    drdr  = dx*dx  + dy*dy  + dz*dz

    rgamma = -dvdr/(drdr*mu)

    term(1) = rgamma*dx
    term(2) = rgamma*dy
    term(3) = rgamma*dz

    !集団Aの原子の速度
    do i=1,i0ownA
       i0=i0listA(i)
       wkv(1,i0) = wkv(1,i0) + term(1)*RTotalMassA
       wkv(2,i0) = wkv(2,i0) + term(2)*RTotalMassA
       wkv(3,i0) = wkv(3,i0) + term(3)*RTotalMassA
    enddo

    !集団Bの原子の速度更新
    do i=1,i0ownB
       i0=i0listB(i)
       wkv(1,i0) = wkv(1,i0) - term(1)*RTotalMassB
       wkv(2,i0) = wkv(2,i0) - term(2)*RTotalMassB
       wkv(3,i0) = wkv(3,i0) - term(3)*RTotalMassB
    enddo
  end subroutine rattle_com
!--------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate mean force
!! \author Yoshimichi Andoh, originated by K.Fujimoto
!<
  subroutine calc_MeanForce()
    use device_numbers
    use md_forces
    use center_of_mass_variables
    use center_of_mass
    use trajectory_mpi
    use md_multiplestep
    use boundary
    implicit none
    integer(4) :: i, j, i0, ierr
    real(8) :: dx, dy, dz,dr
    real(8) :: wk_MFA(3),wk_MFB(3)
    real(8) :: MeanForceA(3), MeanForceB(3)
    !     real(8) :: ForceA, ForceB
    !     integer(4), save :: ic=0
    real(8) :: ftot(3) !! total force at t=t+dt
    include 'mpif.h'

    MeanForceA(:) = 0.0d0
    MeanForceB(:) = 0.0d0

    !集団Aの原子の力の足し込み
    wk_MFA=0d0
    do i=1,i0ownA
       i0=i0listA(i)
       ftot(:)=fshort(:,i0)+fmiddle(:,i0)+flong(:,i0)
       wk_MFA(1)=wk_MFA(1)+ftot(1)
       wk_MFA(2)=wk_MFA(2)+ftot(2)
       wk_MFA(3)=wk_MFA(3)+ftot(3)
    enddo

    !集団Bの原子の力の足し込み
    wk_MFB=0d0
    do i=1,i0ownB
       i0=i0listB(i)
       ftot(:)=fshort(:,i0)+fmiddle(:,i0)+flong(:,i0)
       wk_MFB(1)=wk_MFB(1)+ftot(1)
       wk_MFB(2)=wk_MFB(2)+ftot(2)
       wk_MFB(3)=wk_MFB(3)+ftot(3)
    enddo

    call MPI_Allreduce(wk_MFA,MeanForceA,3,MPI_REAL8, MPI_SUM,MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(wk_MFB,MeanForceB,3,MPI_REAL8, MPI_SUM,MPI_COMM_WORLD, ierr)

#ifdef ZSHAKE
    !     dx = 0.0d0
    !     dy = 0.0d0
    !     dz = CenterofMassA(3) - CenterofMassB(3)
    !     call PBCij(dx, dy, dz)
    !     dr = dsqrt(dx*dx + dy*dy + dz*dz)
    !     dx = dx/dr
    !     dy = dy/dr
    !     dz = dz/dr
#else
    dx = CenterOfMassA(1) - CenterOfMassB(1)
    dy = CenterOfMassA(2) - CenterOfMassB(2)
    dz = CenterOfMassA(3) - CenterOfMassB(3)
    call PBCij(dx, dy, dz)
    dr = dsqrt(dx*dx + dy*dy + dz*dz)
    dx = dx/dr
    dy = dy/dr
    dz = dz/dr
#endif

    MForceA = dx*MeanForceA(1) + dy*MeanForceA(2) + dz*MeanForceA(3)
    MForceB = dx*MeanForceB(1) + dy*MeanForceB(2) + dz*MeanForceB(3)
    drCOM = dr

    MF_AB = (TotalMassB)/(TotalMassA+TotalMassB)*MForceA &
         &       -(TotalMassA)/(TotalMassA+TotalMassB)*MForceB

    !     aveMFA=aveMFA+MForceA
    !     aveMFB=aveMFB+MForceB
    sMF_AB=sMF_AB+MF_AB
  end subroutine calc_MeanForce

end module shake_rattle_roll
