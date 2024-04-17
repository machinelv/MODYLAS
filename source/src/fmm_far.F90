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
!! \brief Module and subroutines to calculate far part of FMM.
!<
!----------------------------------------------------------------------
!>
!! \brief Module to calculate far part of FMM.
!! \author Noriyuki Yoshii, Tatsuya Sakashita, Yoshimichi Andoh, 
!!         Jiachao Zhang, Shin-ichi Ichikawa, Kensuke Iwahashi
!<
module fmm_far
  use precision
  use omp_lib
  use mpi_tool
  use fmm_parameters, only : nmax
  use md_condition, only : mdstep
#include "timing.h90"
  implicit none
  integer(kind=wip),allocatable::nbd_zm(:),nbd_zp(:)
  integer(kind=wip),allocatable::nbd_ym(:),nbd_yp(:)
  integer(kind=wip),allocatable::nbd_xm(:),nbd_xp(:)
  integer(kind=wip),allocatable::mcells_x(:),mcells_y(:),mcells_z(:)
  integer(kind=wip),allocatable::nscllx(:),nsclly(:),nscllz(:)
  integer(kind=wip),allocatable::nscdvx(:),nscdvy(:),nscdvz(:)
  integer(kind=wip),allocatable::wl_size_x(:),wl_size_y(:),wl_size_z(:)
  integer(kind=wip),allocatable::wm_size_x(:),wm_size_y(:),wm_size_z(:)
  integer(kind=wip),allocatable,dimension(:,:,:) :: lddir
  integer(kind=wip),allocatable,dimension(:)     :: nload, nchunk
  ! boundary size of inter-process superposition multipole moment array: wl_local.
  ! for cell number of power of 2 and 3 factors extention.
  integer(4), allocatable :: mbd_x_lvl(:)
  integer(4), allocatable :: mbd_y_lvl(:)
  integer(4), allocatable :: mbd_z_lvl(:)
  integer(4) :: nsbsize
  integer(4) mbd_ux, mbd_uy, mbd_uz  !! mbd for upper level L2L
  integer(4) mbd_lx, mbd_ly, mbd_lz  !! mbd for lower level L2L
  ! basic parameters for power of 2 or power of 3 number cells.
  integer(kind=wip), parameter :: max_ncref = 8   ! maximum referencing cell-range for M2L translation.
  integer(kind=wip), parameter :: ncdir = 2       ! cell-range of direct force calculation.
  integer(kind=wip),allocatable,dimension(:)     :: npowx, npowy, npowz
  complex(kind=wmp),allocatable,dimension(:,:,:,:,:,:) :: shmm1,shmm2,shll1,shll2  
  ! to store precomputed translation operators for M2M and L2L
  real(kind=wvp),allocatable,dimension(:,:,:,:) :: xtall,ytall,ztall  
  ! to store precomputed dx,dy,dz for L2L
  real(kind=wmp), allocatable,dimension(:,:,:,:,:) :: shml 
  ! to store precomputed translation operators for M2L
  real(kind=wvap),allocatable,dimension(:,:,:,:,:,:) :: wl_omp
  real(kind=wmp), allocatable,dimension(:,:,:,:,:,:,:)  :: shml1, shml2 
  real(kind=wmp),allocatable,dimension(:,:,:,:,:,:)  :: sh_rot1, sh_rot2, sh_inv_rot1, sh_inv_rot2  
  ! to store precomputed Wigner D-matrix for direct and inverse rotation
  real(kind=wvap),    allocatable :: wm_omp(:,:,:,:,:,:) ! for reductive summation in P2M
  integer(kind=wip),allocatable :: nlx(:),nly(:),nlz(:)
  integer(kind=wip) :: lgflg=1 ! 0,1,2,3, when switch local2global
  integer(kind=wip), protected :: lm_length
  real(kind=dp), protected :: dxcell, dycell, dzcell  ! length of one subcell
  real(kind=dp), private :: x_center_offset, y_center_offset, z_center_offset  ! for P2M, L2P
  integer(4), allocatable :: ind_m2m(:,:)

#ifdef MOMENT
  ! relative cell address list and cell-pair list for simultaneous multiple cell-pair callculation in M2L.
  integer(kind=wip),allocatable,dimension(:,:,:) :: lddir8p, lddir6p, lddir4p, lddir3p, lddir2p, lddir1p
  integer(kind=wip),allocatable,dimension(:)     :: nload8p, nload6p, nload4p, nload3p, nload2p, nload1p
  integer(kind=wip),allocatable,dimension(:,:,:) :: cp8_list0, cp6_list0, cp4_list0, cp3_list0,cp2_list0, cp1_list0
  integer(kind=wip),allocatable,dimension(:,:,:) :: cp8_list1, cp6_list1, cp4_list1, cp3_list1,cp2_list1, cp1_list1
  integer(kind=wip),allocatable,dimension(:,:,:) :: cp8_irxyz, cp6_irxyz, cp4_irxyz, cp3_irxyz, cp2_irxyz, cp1_irxyz
  integer(kind=wip),allocatable,dimension(:)     :: cp8_count, cp6_count, cp4_count, cp3_count,cp2_count, cp1_count
#endif /** MOMENT **/

  type :: wm_type
     real(kind=wvp), allocatable    :: array(:,:,:,:,:)
  end type wm_type
  type(wm_type), allocatable :: wm(:)

  type :: wl_type
     real(kind=wvap), allocatable    :: array(:,:,:,:,:)
  end type wl_type
!
  type(wl_type), allocatable :: wl(:)

#ifdef DEBUGFCE
  real(8),parameter :: kJ_mol=6.02214129d+23*1.0d-3
#endif

contains

!---------------------------------------------------------------------
!>
!! \brief Subroutine to set lgflg (to switich upper/lower level).
!! \author Yoshimichi Andoh, Shin-ichi Ichikawa
!<
  subroutine fmod_set_lgflg(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    lgflg = ivalue
  end subroutine fmod_set_lgflg
!---------------------------------------------------------------------
!>
!! \brief Subroutine to set nlevel. 
!! \author Yoshimichi Andoh
!<
  subroutine fmod_set_fmm_nlevel(ivalue)
    use domain, only : nlevel, nlevel_input
    implicit none
    integer(kind=wip), intent(in) :: ivalue

    nlevel=ivalue
    nlevel_input=.true.
  end subroutine fmod_set_fmm_nlevel
!---------------------------------------------------------------------
!>
!! \brief Subroutine to set nlevel_input. 
!! \author Yoshimichi Andoh
!<
  subroutine fmod_set_fmm_nlevel_23(ivalue)
    use domain, only : nlevel, nlevel_input
    use mpi_tool ! for debug
    implicit none
    integer(kind=wip), intent(in) :: ivalue
    integer(kind=wip) ntmp

    ntmp=ivalue  ! ncell

    nlevel=0
    do while(mod(ntmp,2)==0)
      ntmp=ntmp/2
      nlevel=nlevel+1
    enddo
    do while(mod(ntmp,3)==0)
      ntmp=ntmp/3
      nlevel=nlevel+1
    enddo
    nlevel_input=.true.
  end subroutine fmod_set_fmm_nlevel_23
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to allocate arrays for multipoles
!! \author Yoshimichi Andoh
!<
  subroutine fmod_alloc_multipole
!----------------------------------------------------------------------
    use mpi_3d_grid
    use domain, only : ncellx, ncelly, ncellz, nlevel
    use fmm_subcell_index
    use openmp_tool, only : nomp
    use comm_fmm_mod, only : &
   &        sclist_xps,sclist_xms,sclist_xpr,sclist_xmr, &
   &        sclist_yps,sclist_yms,sclist_ypr,sclist_ymr, &
   &        sclist_zps,sclist_zms,sclist_zpr,sclist_zmr, &
   &        init_comm_fmm_local_super !, &
    use dr_cntl, only : mbd, mbdp1, mtbd, mtbdp2, &
   &                    get_halo_number
    implicit none
    integer(kind=wip) :: mcell_size_x,mcell_size_y,mcell_size_z
    integer(kind=wip) :: nscellx,nscelly,nscellz
    integer(kind=wip) :: npowz,npowy,npowx
    integer(kind=wip) :: iczg0,icyg0,icxg0
    integer(kind=wip) :: iczg1,icyg1,icxg1
    integer(kind=wip) :: mbd_z, mbd_y, mbd_x
    integer(kind=wip) :: nbsizet
    integer(kind=wip) :: nsczdiv, nscydiv, nscxdiv
! for M2L_lower_level_dr working array wl_omp.
! for M2L_lower_level working array wl_omp.
    integer(kind=wip) :: mbdx_max, mbdy_max, mbdz_max
    integer(kind=wip) :: ilevel,il

    lm_length = (nmax+1)*(nmax+2)/2

    call check_nlevels_along_axes()

!detecting input error 
    if(nlevel==0)then
      if(myrank==0)then
write(0,*) &
& 'ERROR: nlevel=0 in your input. It must be greater than 3.          ', &
& '       If ncell contains 3-powers and/or ncell is not uniform alog ', &
& '       each axis, nlevel must be described explicitly in <fmm> tag.'
      endif
      call modylas_abort
    endif

    if(nlevel <= 1)then
      if(myrank==0) write(0,*) &
   & 'ERROR: nlevel is too small. Now, it is ', nlevel, &
   & ', but it must be ≥ 2, or equivaneltly, ncell must be ≥ 4.'
      call modylas_abort
    endif

!checking ULswitch value. note lgflg = ULswitch in .mddef file
    if(lgflg < 0)then
      if(myrank==0) write(0,*) &
   & 'ERROR: ULswitch must be ≥ 0. ', &
   & 'Your selected ULswitch value in .mddef file is:', lgflg
      call modylas_abort
    elseif(lgflg .ge. nlevel)then
      if(myrank==0) write(0,*) &
   & 'ERROR: ULswitch must be < nlevel. ', &
   & 'Your selected nlevel value is:', nlevel, &
   & '. So, ULswitch must be ≤', nlevel-1
      call modylas_abort
    endif

!allocate arrays with (0:nlevel)
    allocate(nlx(0:nlevel),nly(0:nlevel),nlz(0:nlevel))
    allocate(nbd_xm(0:nlevel),nbd_xp(0:nlevel))
    allocate(nbd_ym(0:nlevel),nbd_yp(0:nlevel))
    allocate(nbd_zm(0:nlevel),nbd_zp(0:nlevel))
    allocate(mcells_x(0:nlevel),mcells_y(0:nlevel))
    allocate(mcells_z(0:nlevel))
    allocate(nscllx(0:nlevel),nsclly(0:nlevel),nscllz(0:nlevel))
    allocate(nscdvx(0:nlevel),nscdvy(0:nlevel),nscdvz(0:nlevel))
    allocate(wl_size_x(0:nlevel),wl_size_y(0:nlevel),wl_size_z(0:nlevel))
    allocate(wm_size_x(0:nlevel),wm_size_y(0:nlevel),wm_size_z(0:nlevel))
    ! structured array !
    allocate(wm(0:nlevel))
    allocate(wl(0:nlevel))

    allocate(sclist_xps(0:nlevel))
    allocate(sclist_xms(0:nlevel))
    allocate(sclist_xpr(0:nlevel))
    allocate(sclist_xmr(0:nlevel))
    allocate(sclist_yps(0:nlevel))
    allocate(sclist_yms(0:nlevel))
    allocate(sclist_ypr(0:nlevel))
    allocate(sclist_ymr(0:nlevel))
    allocate(sclist_zps(0:nlevel))
    allocate(sclist_zms(0:nlevel))
    allocate(sclist_zpr(0:nlevel))
    allocate(sclist_zmr(0:nlevel))
! for p2p3 DR method.
    allocate(mbd_x_lvl(0:nlevel))
    allocate(mbd_y_lvl(0:nlevel))
    allocate(mbd_z_lvl(0:nlevel))
    
    !### set nlx ###!
    call cmerge_number(nlevel, ncellx, npx, nlx)
    nlx(nlevel) = 1
    !### set nly ###!
    call cmerge_number(nlevel, ncelly, npy, nly)
    nly(nlevel) = 1
    !### set nlz ###!
    call cmerge_number(nlevel, ncellz, npz, nlz)
    nlz(nlevel) = 1

    mcell_size_z=1
    mcell_size_y=1
    mcell_size_x=1
    mbdx_max  = 0
    mbdy_max  = 0
    mbdz_max  = 0

    do ilevel=0,nlevel
    !### set cell parameters ##
       mcells_x(ilevel)=mcell_size_x
       mcells_y(ilevel)=mcell_size_y
       mcells_z(ilevel)=mcell_size_z

       wl_size_z(ilevel) = get_cell_size(ncellz, mcell_size_z, npz)  ! NOT same as lzdiv
       wl_size_y(ilevel) = get_cell_size(ncelly, mcell_size_y, npy)  ! NOT same as lzdiv
       wl_size_x(ilevel) = get_cell_size(ncellx, mcell_size_x, npx)  ! NOT same as lzdiv

       nsczdiv = wl_size_z(ilevel)
       nscydiv = wl_size_y(ilevel)
       nscxdiv = wl_size_x(ilevel)

       nscdvz(ilevel)=nsczdiv
       nscdvy(ilevel)=nscydiv
       nscdvx(ilevel)=nscxdiv

    !### allocate wm WITH halo ###
       nscellz= (ncellz-1)/mcell_size_z+ 1
       nscelly= (ncelly-1)/mcell_size_y+ 1
       nscellx= (ncellx-1)/mcell_size_x+ 1
       nscllz(ilevel)=nscellz
       nsclly(ilevel)=nscelly
       nscllx(ilevel)=nscellx
!      ^^^ z-axis ^^^
       call get_halo_number(ipz, &
     &           npz,nlz(ilevel),nscellz,  &
     &           nbd_zm(ilevel),nbd_zp(ilevel),mbd_z_lvl(ilevel))

!      ^^^ y-axis ^^^
       call get_halo_number(ipy, &
     &           npy,nly(ilevel),nscelly, &
     &           nbd_ym(ilevel),nbd_yp(ilevel),mbd_y_lvl(ilevel))

!      ^^^ x-axis ^^^
       call get_halo_number(ipx, &
     &           npx,nlx(ilevel),nscellx, &
     &           nbd_xm(ilevel),nbd_xp(ilevel),mbd_x_lvl(ilevel))

       if(mbdx_max < mbd_x_lvl(ilevel)) mbdx_max = mbd_x_lvl(ilevel)
       if(mbdy_max < mbd_y_lvl(ilevel)) mbdy_max = mbd_y_lvl(ilevel)
       if(mbdz_max < mbd_z_lvl(ilevel)) mbdz_max = mbd_z_lvl(ilevel)

! constants used for allocation
       mbd_z=mbd_z_lvl(ilevel)
       mbd_y=mbd_y_lvl(ilevel)
       mbd_x=mbd_x_lvl(ilevel)

       if(ilevel.le.lgflg) then
          wm_size_z(ilevel) = nbd_zm(ilevel)+nsczdiv+nbd_zp(ilevel)
          wm_size_y(ilevel) = nbd_ym(ilevel)+nscydiv+nbd_yp(ilevel)
          wm_size_x(ilevel) = nbd_xm(ilevel)+nscxdiv+nbd_xp(ilevel)
       else
          wm_size_z(ilevel) = nscellz
          wm_size_y(ilevel) = nscelly
          wm_size_x(ilevel) = nscellx
       endif

       allocate( wm(ilevel)%array(lm_length, 2, wm_size_z(ilevel), &
    &                                           wm_size_y(ilevel), &
    &                                           wm_size_x(ilevel)) )

       if(ilevel.le.lgflg) then
      !local
          allocate(wl(ilevel)%array(lm_length, 2, mbd_z+nsczdiv+mbd_z, &
               &                                  mbd_y+nscydiv+mbd_y, &
               &                                  mbd_x+nscxdiv+mbd_x))

       else
      !global
          allocate(wl(ilevel)%array(lm_length, 2, nsczdiv, nscydiv, nscxdiv))
       endif

    !### allocate arrays for lowever-level communication ###
       if(ilevel.le.lgflg) then
          nbsizet = (mbd_y*2+nscydiv)*(mbd_z*2+nsczdiv)*mbd_x
          allocate(sclist_xps(ilevel)%iarray(3,0:nbsizet) )
          allocate(sclist_xms(ilevel)%iarray(3,0:nbsizet) )
          allocate(sclist_xpr(ilevel)%iarray(3,0:nbsizet) )
          allocate(sclist_xmr(ilevel)%iarray(3,0:nbsizet) )
          nbsizet = (mbd_z*2+nsczdiv)*(mbd_x*2+nscxdiv)*mbd_y
          allocate(sclist_yps(ilevel)%iarray(3,0:nbsizet) )
          allocate(sclist_yms(ilevel)%iarray(3,0:nbsizet) )
          allocate(sclist_ypr(ilevel)%iarray(3,0:nbsizet) )
          allocate(sclist_ymr(ilevel)%iarray(3,0:nbsizet) )
          nbsizet = (mbd_x*2+nscxdiv)*(mbd_y*2+nscydiv)*mbd_z
          allocate(sclist_zps(ilevel)%iarray(3,0:nbsizet) )
          allocate(sclist_zms(ilevel)%iarray(3,0:nbsizet) )
          allocate(sclist_zpr(ilevel)%iarray(3,0:nbsizet) )
          allocate(sclist_zmr(ilevel)%iarray(3,0:nbsizet) )
          call init_comm_fmm_local_super( &
       &      sclist_zps(ilevel)%iarray, sclist_zms(ilevel)%iarray, &
       &      sclist_zpr(ilevel)%iarray, sclist_zmr(ilevel)%iarray, &
       &      sclist_yps(ilevel)%iarray, sclist_yms(ilevel)%iarray, &
       &      sclist_ypr(ilevel)%iarray, sclist_ymr(ilevel)%iarray, &
       &      sclist_xps(ilevel)%iarray, sclist_xms(ilevel)%iarray, &
       &      sclist_xpr(ilevel)%iarray, sclist_xmr(ilevel)%iarray, &
       &      nsczdiv, nscydiv, nscxdiv, &
       &      mbd_z, mbd_y, mbd_x)
       endif

       !=== update for next level ===!
       mcell_size_x = mcell_size_x * nlx(ilevel)
       mcell_size_y = mcell_size_y * nly(ilevel)
       mcell_size_z = mcell_size_z * nlz(ilevel)
    ENDDO ! ilevel

!add for DR method.
    nsbsize = max((mbd_x_lvl(0)*2+nscdvx(0))*(mbd_y_lvl(0)*2+nscdvy(0)), &
          &       (mbd_y_lvl(0)*2+nscdvy(0))*(mbd_z_lvl(0)*2+nscdvz(0))  )
    nsbsize = max((mbd_z_lvl(0)*2+nscdvz(0))*(mbd_x_lvl(0)*2+nscdvx(0)), &
          &       nsbsize  )
    nsbsize = max(mbd_x_lvl(0),mbd_y_lvl(0),mbd_z_lvl(0))*2 * nsbsize

!   ! For reductive summation in P2M
    allocate( wm_omp(lm_length, 2, wm_size_z(0), wm_size_y(0), wm_size_x(0), 0:nomp-1) )
    allocate( wl_omp( &
         &                  lm_length,2, &
         &                  1-mbdz_max:wl_size_z(0)+mbdz_max, &
         &                  1-mbdy_max:wl_size_y(0)+mbdy_max, &
         &                  1-mbdx_max:wl_size_x(0)+mbdx_max, 0:nomp-1) )

  end subroutine fmod_alloc_multipole
!---------------------------------------------------------------------
!>
!! \brief  Subroutine to check correspondence of nlevels along each axis 
!! \author Yoshimichi Andoh
!<
  subroutine check_nlevels_along_axes()
!---------------------------------------------------------------------
    use domain, only : nlevel, ncellx, ncelly, ncellz, &
   &                   nlevel_input, ncellx_input, ncelly_input, ncellz_input
    implicit none
    integer(kind=wip) :: nlevelx, nlevely, nlevelz
    integer(kind=wip) :: ntmp

      !^^ x ^^!
      nlevelx=0
      ntmp=ncellx
      do while(mod(ntmp,2)==0)
        ntmp=ntmp/2
        nlevelx=nlevelx+1
      enddo
      do while(mod(ntmp,3)==0)
        ntmp=ntmp/3
        nlevelx=nlevelx+1
      enddo

      !^^ y ^^!
      nlevely=0
      ntmp=ncelly
      do while(mod(ntmp,2)==0)
        ntmp=ntmp/2
        nlevely=nlevely+1
      enddo
      do while(mod(ntmp,3)==0)
        ntmp=ntmp/3
        nlevely=nlevely+1
      enddo

      !^^ z ^^!
      nlevelz=0
      ntmp=ncellz
      do while(mod(ntmp,2)==0)
        ntmp=ntmp/2
        nlevelz=nlevelz+1
      enddo
      do while(mod(ntmp,3)==0)
        ntmp=ntmp/3
        nlevelz=nlevelz+1
      enddo
!debug
!     if(myrank==0) write(*,*) 'nlevelxyz=',nlevelx,nlevely,nlevelz
!debug
      if( nlevelx == nlevely .and. nlevely == nlevelz )then
        continue
      else
        if(myrank==0)then
          write(0,*) 'ERROR: nlevel along x, y, z axes does not match.'
          write(0,*) '  ncellx, nlevelx=', ncellx, nlevelx
          write(0,*) '  ncelly, nlevely=', ncelly, nlevely
          write(0,*) '  ncellz, nlevelz=', ncellz, nlevelz
        endif
        call modylas_abort
      endif

      nlevel=nlevelx

  end subroutine check_nlevels_along_axes
!---------------------------------------------------------------------
!>
!! \brief  Subroutine to check factors in ncell and nprocs
!! \author Shin-ichi Ichikawa
!<
!---------------------------------------------------------------------
  subroutine cmerge_number(nlevel, ncellxyz, npxyz, nlxyz)
    integer(kind=wip) :: ncellxyz, npxyz, nlevel
    integer(kind=wip) :: nlxyz(0:nlevel)
    integer(kind=wip) :: pp, qp, nc, mc, maxcell, maxproc, il, nerror
!
!.... Cell Merge Number Definition by Arbitrary Order.             ....
!....                                                              ....
!.... cell merge number definition formula:                        ....
!....  cells = 2^n*3^m,  procs = 2^p*3^q                           ....
!....      +-------+-----------+-----------+------------+          ....
!....      |       |    n > p  |   n = p   |   n < p    |          ....
!....      +-------+-----------+-----------+------------+          ....
!....      | m > q |  (n - 1)  |   m - 1   |  input     |          ....
!....      |       | or m - 1  |           |  error     |          ....
!....      +-------+-----------+-----------+------------+          ....
!....      | m = q |   n - 1   | n - 1 >= 0: n - 1      |          ....
!....      +-------+-----------+ n - 1 < 0 :            |          ....
!....      | m < q |  input    |  m - 1 >= 0: m - 1     |          ....
!....      |       |  error    |  m - 1 < 0 : level err |          ....
!....      +-------+-----------+-----------+------------+          ....
!....      ( ) means primary choice.                               ....

    nc = 0
    maxcell = ncellxyz
    do while(mod(maxcell,2).eq.0)
       maxcell = maxcell/2
       nc = nc + 1
    enddo
    mc = 0
    do while(mod(maxcell,3).eq.0)
       maxcell = maxcell/3
       mc = mc + 1
    enddo
    pp = 0
    maxproc = npxyz
    do while(mod(maxproc,2).eq.0)
       maxproc = maxproc/2
       pp = pp + 1
    enddo
    qp = 0
    do while(mod(maxproc,3).eq.0)
       maxproc = maxproc/3
       qp = qp + 1
    enddo
    if(maxcell /= 1 .or. maxproc /= 1) then
      write(0,*) 'ERROR: number of cells or number of processes ', &
          &      'have factor other than 2 or 3. ncellxyz: ',ncellxyz, &
          &      ' npxyz: ',npxyz
      call modylas_abort
    endif

!    write(myrank+NBFU,*) "+++ cmerge_number: nc,mc,pp,qp: ",nc,mc,pp,qp

    nerror = 0
    do il = 0, nlevel-1
      if(nc > pp) then
        if(mc >= qp) then
            nc = nc - 1
            nlxyz(il) = 2
        else
            nerror = 2 ; exit  ! reason 2 input cnum or pnum error.
        end if
      else if(nc == pp) then
        if(mc > qp) then
            mc = mc - 1
            nlxyz(il) = 3
        else
          if(nc - 1 >= 0) then
            nc = nc - 1
            nlxyz(il) = 2
          else if(mc - 1 >= 0) then
            mc = mc - 1
            nlxyz(il) = 3
          else
            nerror = 1 ; exit  ! reason 1 nlevel error.
          endif
        end if    ! mc <= qp.
      else   ! nc < pp.
        if(mc > qp) then
            nerror = 2 ; exit  ! reason 2 input cnum or pnum error.
        else
          if(nc-1 >= 0) then
            nc = nc - 1
            nlxyz(il) = 2
          else if(mc-1 >= 0) then
            mc = mc - 1
            nlxyz(il) = 3
          else
            nerror = 1 ; exit   ! reason 1 nlevel error.
          endif
        end if    ! mc <= qp.
      end if    ! nc < pp.

!      write(myrank+NBFU,*) "*** cmerge_number: il,nlxyz(il): ",il,nlxyz(il)

    end do

    if(nerror == 1) then
      write(0,*) 'ERROR: number of levels: nlevel may be wrong because ', &
          &      'cell number is already one ', &
          &      'even when the level is less than or equal to nlevel-1.'
      call modylas_abort
    endif
    if(nerror == 2) then
      write(0,*) 'ERROR: input cell number can not be devided evenly by process number.'
      call modylas_abort
    endif

  end subroutine cmerge_number

!---------------------------------------------------------------------
!>
!! \brief Subroutine for initialization of FMM calculation
!! \author Noriyuki Yoshii, Kensuke Iwahashi
!<
  subroutine init_fmm()
!---------------------------------------------------------------------
    use fmm_ewald, only : init_fmm_ewald, init_fmm_ewald_charged, &
             &            calc_qtot
    use fmm_l_m_index, only : init_m2m_indices
    use domain, only : nlevel ! debug
    use mpi_tool
    use ewald_mod, only : bQcalc1
    implicit none
    integer(kind=wip)::l,ilevel

    call generate_relative_cell_list_for_M2L()

#ifdef MOMENT
    ! add for relative cell list with assorted cell-pairs.
    call generate_assorted_relative_cell_list_for_M2L()
    ! add for cell-pair list with assorted cell-pairs.
    call generate_assorted_cellpair_list_for_M2L()
#endif
    
    call initialize_translation_operators()

    call init_m2m_indices(nmax, ind_m2m)

    call calc_qtot()  !calc system charge

    call init_fmm_ewald(nmax, lm_length)

    if(bQcalc1) call init_fmm_ewald_charged(nmax, lm_length) ! only once call

  end subroutine init_fmm

!---------------------------------------------------------------------
!>
!! \brief Subroutine which sotres multipole moments in all subcells of all levels
!! \author Noriyuki Yoshii, Kensuke Iwahashi, Tatsuya Sakashita
!<
!---------------------------------------------------------------------
  subroutine calc_mm ! P2M, M2M
!---------------------------------------------------------------------
    use domain, only : nlevel
    use comm_fmm_mod, only : comm_fmm_local_top, comm_fmm_local_multi_dr

    implicit none
    integer(kind=wip) :: il, il_higher
    integer(kind=wip) :: nomp
    integer(kind=wip) nbsize

    nomp = 1
!$  nomp = omp_get_max_threads()
    
    TIME_START(TM_P2M)
    call P2M(wm(0)%array) ! calculate_particle2multipole
    TIME_STOP(TM_P2M)

    DO il=0,nlevel-1

       if(il.le.lgflg)then
          TIME_BARRIER(TMB_FMM_COMM_LOCAL)
          TIME_START(TM_FMM_COMM_LOCAL)
          nbsize = min( max(nbd_zm(il),nbd_zp(il)), nscdvz(il) ) ! debugged
          call comm_fmm_local_multi_dr(lm_length, &
               &       nlz(il), nly(il), nlx(il), &
               &       wm(il)%array,nbsize,nscdvz(il),nscdvy(il),nscdvx(il), &
               &       nscllz(il),nsclly(il),nscllx(il), &
               &       nbd_zm(il), nbd_ym(il), nbd_xm(il), &
               &       nbd_zp(il), nbd_yp(il), nbd_xp(il), il)
          TIME_STOP(TM_FMM_COMM_LOCAL)
       else
          TIME_BARRIER(TMB_FMM_COMM_GLOBAL)
          TIME_START(TM_FMM_COMM_GLOBAL)
          call comm_fmm_local_top(lm_length, &
               &       wm(il)%array, &
               &       wl_size_z(il),wl_size_y(il),wl_size_x(il), &
               &       nscllz(il),nsclly(il),nscllx(il), &
               &       mcells_z(il),mcells_y(il),mcells_x(il))
          TIME_STOP(TM_FMM_COMM_GLOBAL)
       endif
       TIME_START(TM_M2M)
       il_higher = il + 1
       call M2M(il,il_higher, wm(il)%array,wm(il_higher)%array,nomp)
       TIME_STOP(TM_M2M)
       
    ENDDO ! il
  end subroutine calc_mm
!---------------------------------------------------------------------
!>
!! \brief Subroutine which calculate multipole moment (Particle to Multipole)
!! \author Noriyuki Yoshii, Kensuke Iwahashi
!<
!---------------------------------------------------------------------
  subroutine P2M(wm_lowest) ! calculate_particle2multipole
!---------------------------------------------------------------------
    use trajectory_mpi, only : wkxyz
    use coulomb_mod, only : chgv
    use domain, only : ncellx, ncelly, ncellz, lxdiv, lydiv, lzdiv, ixmin, iymin, izmin
    use subcell, only : tag, na_per_cell, m2i
    use param, only : paranum
    use unit_cell, only : cellx, celly, cellz
    use regular_solid_harmonics_cartesian
    use dr_cntl, only : nbd

    implicit none
    real(kind=wvp), intent(out)    :: wm_lowest(lm_length  , 2, &
   &                                     wm_size_z(0), wm_size_y(0), wm_size_x(0))
    real(kind=dp) :: q0
    integer(kind=wip) :: icx0, icy0, icz0
    integer(kind=wip) :: i,ii
    integer(kind=wip) :: m1
    real(kind=dp) :: xta, yta, zta
    complex(kind=dp) :: R_array(lm_length)
    integer(kind=wip) :: i0,ipar,iam,nomp
    integer(kind=wip) :: icx00,icy00,icz00
    integer(kind=wip) :: nbound_zm, nbound_ym, nbound_xm
    real(kind=dp) :: factorial_inverse_array(0:2*nmax)
    real(kind=dp) :: minus_double_factorial_array(0:nmax)
    real(kind=dp) :: x_center, y_center, z_center

    nomp = 1
    iam = 0
!$  nomp = omp_get_max_threads()

    nbound_zm=nbd_zm(0)
    nbound_ym=nbd_ym(0)
    nbound_xm=nbd_xm(0)
    
    call prepare_factorial_arrays(nmax, factorial_inverse_array, minus_double_factorial_array)

!$omp parallel default(none) &
!$omp& private(iam) &
!$omp& private(ii) &
!$omp& private(icx0,icy0,icz0,icx00,icy00,icz00) &
!$omp& private(x_center, y_center, z_center) &
!$omp& private(i) &
!$omp& private(q0,xta,yta,zta) &
!$omp& private(m1) &
!$omp& private(i0,ipar) &
!$omp& private(R_array) &
!$omp& shared(dxcell, dycell, dzcell, x_center_offset, y_center_offset, z_center_offset) &
!$omp& shared(factorial_inverse_array, minus_double_factorial_array) &
!$omp& shared(nomp,nmax,lm_length) &
!$omp& shared(wkxyz,chgv) &
!$omp& shared(paranum) &
!$omp& shared(tag,na_per_cell,m2i) &
!$omp& shared(lxdiv,lydiv,lzdiv) &
!$omp& shared(nbound_xm,nbound_ym,nbound_zm) &
!$omp& shared(wm_lowest,wm_omp)
!$  iam = omp_get_thread_num()
    do ii=1,lxdiv*lydiv*lzdiv
       icz0=mod(ii-1,lzdiv)     +nbd 
       icy0=mod(ii-1,lzdiv*lydiv)
       icy0=icy0/lzdiv          +nbd 
       icx0=(ii-1)/(lzdiv*lydiv)+nbd 
       icx00=icx0-nbd+1+nbound_xm ! local -> local FMM
       icy00=icy0-nbd+1+nbound_ym ! local -> local FMM
       icz00=icz0-nbd+1+nbound_zm ! local -> local FMM
       x_center = icx0 * dxcell + x_center_offset
       y_center = icy0 * dycell + y_center_offset
       z_center = icz0 * dzcell + z_center_offset
       wm_omp(:,1,icz00,icy00,icx00,iam) = Cwvzero ! wvap
       wm_omp(:,2,icz00,icy00,icx00,iam) = Cwvzero ! wvap
!$omp do
       do i0=tag(icz0,icy0,icx0), &
     &       tag(icz0,icy0,icx0)+na_per_cell(icz0,icy0,icx0)-1
          i = m2i(i0)
          ipar = paranum(i)
          q0 = chgv(ipar)

#if defined(PRECISION_M2L_MIX) || defined(PRECISION_M2L_SP)
            !In order to use Single-Precision, Reduce units to avoid overflow. 
            xta = x_center - wkxyz(1,i0) * 1d10
            yta = y_center - wkxyz(2,i0) * 1d10
            zta = z_center - wkxyz(3,i0) * 1d10
#else
            !By default
            xta = x_center - wkxyz(1,i0)
            yta = y_center - wkxyz(2,i0)
            zta = z_center - wkxyz(3,i0)
#endif

          call calculate_regular_harmonics_1dim_array(nmax, &
               & factorial_inverse_array, minus_double_factorial_array, &
               & xta, yta, zta, R_array)

          !*** calculate m-matrices
          do m1 = 1, lm_length
             wm_omp(m1,1,icz00,icy00,icx00,iam) = wm_omp(m1,1,icz00,icy00,icx00,iam) &
                  & + q0 * real( R_array(m1) )
             wm_omp(m1,2,icz00,icy00,icx00,iam) = wm_omp(m1,2,icz00,icy00,icx00,iam) &
                  & + q0 * imag( R_array(m1) )
          enddo ! m1
       enddo ! i0
!$omp end do !!!!nowait <= nowait cause bug
    enddo ! ii
!$omp do
    do ii=1,lxdiv*lydiv*lzdiv
       icz0=mod(ii-1,lzdiv)     +nbd 
       icy0=mod(ii-1,lzdiv*lydiv)
       icy0=icy0/lzdiv          +nbd 
       icx0=(ii-1)/(lzdiv*lydiv)+nbd 
       icx00=icx0-nbd+1+nbound_xm ! local -> local FMM
       icy00=icy0-nbd+1+nbound_ym ! local -> local FMM
       icz00=icz0-nbd+1+nbound_zm ! local -> local FMM
       wm_lowest(:,1,icz00,icy00,icx00) = Cvzero ! wvp
       wm_lowest(:,2,icz00,icy00,icx00) = Cvzero ! wvp
       do iam=0,nomp-1
          wm_lowest(:,1,icz00,icy00,icx00) = &
        & wm_lowest(:,1,icz00,icy00,icx00) + wm_omp(:,1,icz00,icy00,icx00,iam)
          wm_lowest(:,2,icz00,icy00,icx00) = &
        & wm_lowest(:,2,icz00,icy00,icx00) + wm_omp(:,2,icz00,icy00,icx00,iam)
       enddo ! iam
    enddo ! ii
!$omp end do nowait
!$omp end parallel

#ifdef DEBUG_MTDFMM
    do ii=1,lxdiv*lydiv*lzdiv
       icz0=mod(ii-1,lzdiv)     +nbd 
       icy0=mod(ii-1,lzdiv*lydiv)
       icy0=icy0/lzdiv          +nbd 
       icx0=(ii-1)/(lzdiv*lydiv)+nbd 
       icx00=icx0-nbd+1+nbound_xm ! local -> local FMM
       icy00=icy0-nbd+1+nbound_ym ! local -> local FMM
       icz00=icz0-nbd+1+nbound_zm ! local -> local FMM
       do m1=1,lm_length
          write(myrank+10000,'(2es22.15)') wm_lowest(m1,1:2,icz00,icy00,icx00)
       enddo
    enddo ! ii
    call flush(myrank+10000)
!   call mpiend ! P2M
#endif /**DEBUG_MTDFMM**/
  end subroutine P2M
!---------------------------------------------------------------------
!>
!! \brief Subroutine which calculate M2M transformation
!! \author Noriyuki Yoshii, Kensuke Iwahashi
!<
  subroutine M2M(ilevel1,ilevel2, wm_preM2M,wm_postM2M,nomp)
!---------------------------------------------------------------------
    use mpi_3d_grid
    use domain, only : ncellx, ncelly, ncellz
    use fmm_l_m_index
    use dr_cntl, only : mbd, mtbd
    implicit none
    integer(kind=wip), intent(in) :: ilevel1, ilevel2
    real(kind=wvp), intent(in) :: wm_preM2M(lm_length,2,wm_size_z(ilevel1), &
   &                                                    wm_size_y(ilevel1), &
   &                                                    wm_size_x(ilevel1))
    real(kind=wvp), intent(out)::wm_postM2M(lm_length,2,wm_size_z(ilevel2), &
   &                                                    wm_size_y(ilevel2), &
   &                                                    wm_size_x(ilevel2))
    real(kind=wvp) :: postM2M_omp(lm_length,2,wm_size_z(ilevel2), &
   &                                          wm_size_y(ilevel2), &
   &                                          wm_size_x(ilevel2), 0:nomp-1)
    integer(kind=wip) :: icxyz, j, m1, m2
    integer(kind=wip) :: nsczdiv,nscydiv,nscxdiv
    integer(kind=wip) :: x_pre_shift, y_pre_shift, z_pre_shift
    integer(kind=wip) :: x_post_shift, y_post_shift, z_post_shift
    integer(kind=wip) :: nbound_zm, nbound_ym, nbound_xm
    integer(kind=wip) :: iam,nomp
    integer(kind=wip) :: icx0,icy0,icz0
    integer(kind=wip) :: icx_pre, icy_pre, icz_pre
    integer(kind=wip) :: icx_post,icy_post,icz_post
    integer(kind=wip) :: mx,my,mz
    iam = 0

    !NOTE: ilevel2 -> level n
    !      ilevel1 -> level n-1
    !
    !### parameters for level n ###
    nsczdiv = wl_size_z(ilevel2)
    nscydiv = wl_size_y(ilevel2)
    nscxdiv = wl_size_x(ilevel2)

    !### parameters for level n-1 ###
    nbound_xm=nbd_xm(ilevel1)
    nbound_ym=nbd_ym(ilevel1)
    nbound_zm=nbd_zm(ilevel1)

    ! process >= number of supercell
    if(    nlz(ilevel1)==2)then
       if(ncellz / mcells_z(ilevel1) / npz .le.1) nbound_zm=mbd  !! YA modified
    elseif(nlz(ilevel1)==3)then
       if(ncellz / mcells_z(ilevel1) / npz .le.1) nbound_zm=mtbd
    endif
    if(    nly(ilevel1)==2)then
       if(ncelly / mcells_y(ilevel1) / npy .le.1) nbound_ym=mbd  !! YA modified
    elseif(nly(ilevel1)==3)then
       if(ncelly / mcells_y(ilevel1) / npy .le.1) nbound_ym=mtbd
    endif
    if(    nlx(ilevel1)==2)then
       if(ncellx / mcells_x(ilevel1) / npx .le.1) nbound_xm=mbd  !! YA modified
    elseif(nlx(ilevel1)==3)then
       if(ncellx / mcells_x(ilevel1) / npx .le.1) nbound_xm=mtbd
    endif

    if(ilevel1 < lgflg)then
       z_post_shift = nbd_zm(ilevel2) ! local (nl) -> local FMM (nl)
       y_post_shift = nbd_ym(ilevel2) ! local (nl) -> local FMM (nl)
       x_post_shift = nbd_xm(ilevel2) ! local (nl) -> local FMM (nl)
    else
       z_post_shift = get_izmin_in_fmm(ilevel2) - 1 ! local (nl) -> global (nl)
       y_post_shift = get_iymin_in_fmm(ilevel2) - 1 ! local (nl) -> global (nl)
       x_post_shift = get_ixmin_in_fmm(ilevel2) - 1 ! local (nl) -> global (nl)
    endif

    if(ilevel1 <= lgflg)then
       z_pre_shift = -(nlz(ilevel1)-1)+nbound_zm ! local (nl) -> local FMM (nl-1)
       y_pre_shift = -(nly(ilevel1)-1)+nbound_ym ! local (nl) -> local FMM (nl-1)
       x_pre_shift = -(nlx(ilevel1)-1)+nbound_xm ! local (nl) -> local FMM (nl-1)
    else
       z_pre_shift = -(nlz(ilevel1)-1) + (get_izmin_in_fmm(ilevel2) - 1) * nlz(ilevel1) ! local (nl) -> global (nl-1)
       y_pre_shift = -(nly(ilevel1)-1) + (get_iymin_in_fmm(ilevel2) - 1) * nly(ilevel1) ! local (nl) -> global (nl-1)
       x_pre_shift = -(nlx(ilevel1)-1) + (get_ixmin_in_fmm(ilevel2) - 1) * nlx(ilevel1) ! local (nl) -> global (nl-1)
    endif
    
!$omp parallel default(none) &
!$omp& shared(nmax) &
!$omp& private(iam) &
!$omp& private(icxyz,m1,m2,j,mx,my,mz) &
!$omp& private(icx0,icy0,icz0) &
!$omp& private(icx_pre,icy_pre,icz_pre,icx_post,icy_post,icz_post) &
!$omp& shared(nlx,nly,nlz) &
!$omp& shared(x_pre_shift, y_pre_shift, z_pre_shift) &
!$omp& shared(x_post_shift, y_post_shift, z_post_shift) &
!$omp& shared(myrank) &
!$omp& shared(ilevel1,ind_m2m,len_m2m,shmm1,shmm2) &
!$omp& shared(nsczdiv,nscydiv,nscxdiv) &
!$omp& shared(wm_postM2M,wm_preM2M) &
!$omp& shared(postM2M_omp,nomp)
!$  iam = omp_get_thread_num()
    do icxyz=1,nscxdiv*nscydiv*nsczdiv
       icz0=mod(icxyz-1,nsczdiv)       +1
       icy0=mod(icxyz-1,nsczdiv*nscydiv)
       icy0=icy0/nsczdiv               +1
       icx0=(icxyz-1)/(nsczdiv*nscydiv)+1

       icz_post = icz0 + z_post_shift
       icy_post = icy0 + y_post_shift
       icx_post = icx0 + x_post_shift
       icz_pre = icz0 * nlz(ilevel1) + z_pre_shift
       icy_pre = icy0 * nly(ilevel1) + y_pre_shift
       icx_pre = icx0 * nlx(ilevel1) + x_pre_shift
       postM2M_omp(:,1,icz_post,icy_post,icx_post,iam)=Cvzero ! wvp
       postM2M_omp(:,2,icz_post,icy_post,icx_post,iam)=Cvzero ! wvp
       do mx=0,nlx(ilevel1)-1
          do my=0,nly(ilevel1)-1
             do mz=0,nlz(ilevel1)-1
!$omp do
                do j = 1, len_m2m
                   m1 = ind_m2m(1,j)
                   m2 = ind_m2m(2,j)
                   postM2M_omp(m1,1,icz_post,icy_post,icx_post,iam) = postM2M_omp(m1,1,icz_post,icy_post,icx_post,iam) &
                        & + wm_preM2M(m2,1,icz_pre+mz,icy_pre+my,icx_pre+mx) * real(shmm1(m2,m1,mx,my,mz,ilevel1)) &
                        & + wm_preM2M(m2,2,icz_pre+mz,icy_pre+my,icx_pre+mx) * real(shmm2(m2,m1,mx,my,mz,ilevel1))
                   postM2M_omp(m1,2,icz_post,icy_post,icx_post,iam) = postM2M_omp(m1,2,icz_post,icy_post,icx_post,iam) &
                        & + wm_preM2M(m2,1,icz_pre+mz,icy_pre+my,icx_pre+mx) * imag(shmm1(m2,m1,mx,my,mz,ilevel1)) &
                        & + wm_preM2M(m2,2,icz_pre+mz,icy_pre+my,icx_pre+mx) * imag(shmm2(m2,m1,mx,my,mz,ilevel1))
                enddo
!$omp end do
             enddo
          enddo
       enddo
    enddo ! icxyz
!$omp do
    do icxyz=1,nscxdiv*nscydiv*nsczdiv
       icz0=mod(icxyz-1,nsczdiv)       +1
       icy0=mod(icxyz-1,nsczdiv*nscydiv)
       icy0=icy0/nsczdiv               +1
       icx0=(icxyz-1)/(nsczdiv*nscydiv)+1
       icz_post = icz0 + z_post_shift
       icy_post = icy0 + y_post_shift
       icx_post = icx0 + x_post_shift
       wm_postM2M(:,1,icz_post,icy_post,icx_post)=Cvzero ! wvp
       wm_postM2M(:,2,icz_post,icy_post,icx_post)=Cvzero ! wvp
       do iam=0,nomp-1
          wm_postM2M(:,1:2,icz_post,icy_post,icx_post) = &
          wm_postM2M(:,1:2,icz_post,icy_post,icx_post) &
               & + postM2M_omp(:,1:2,icz_post,icy_post,icx_post,iam)
       enddo
    enddo
!$omp end do
!$omp end parallel

#ifdef DEBUG_MTDFMM
    do icxyz=1,nscxdiv*nscydiv*nsczdiv
       icz0=mod(icxyz-1,nsczdiv)       +1
       icy0=mod(icxyz-1,nsczdiv*nscydiv)
       icy0=icy0/nsczdiv               +1
       icx0=(icxyz-1)/(nsczdiv*nscydiv)+1
       icz_post = icz0 + z_post_shift
       icy_post = icy0 + y_post_shift
       icx_post = icx0 + x_post_shift
       do m1=1,lm_length
       write(myrank+ilevel1*1000+20000,'(2es22.15)') wm_postM2M(m1,1:2,icz_post,icy_post,icx_post)
       enddo ! m1
    enddo
    call flush(myrank+ilevel1*1000+20000)
#endif /**DEBUG_MTDFMM**/
  end subroutine M2M
!---------------------------------------------------------------------
!>
!! \brief Subroutine which calculate lattice sum for multipoles
!! \author Noriyuki Yoshii
!<
  subroutine run_scaled_fmm_ewald(wm_highest, wl_highest)
!---------------------------------------------------------------------
    use unit_cell, only : cellx
    use fmm_ewald, only : run_fmm_ewald
    use fmm_l_m_index
    implicit none
    real(kind=wvp),intent(in)     :: wm_highest(lm_length,2)  ! for reductive summation in P2M
    real(kind=wvap),intent(inout) :: wl_highest(lm_length,2)  ! for reductive summation in P2M
    complex(8) :: winput(lm_length)  !Use original array
    complex(8) :: woutput(lm_length) !Use original array
    integer(4) :: j, k, m1

    !*** Scaling multipole moment by unit cell length at highest level
    do j=0,nmax
       do k=0,j
          m1 = translate_l_m_to_1dim(j, k)
!
#if defined(PRECISION_M2L_MIX) || defined(PRECISION_M2L_SP)
            winput(m1) = dcmplx( wm_highest(m1,1), wm_highest(m1,2) ) / (cellx*1d10)**j
#else /**not SP**/
            winput(m1) = dcmplx( wm_highest(m1,1), wm_highest(m1,2) ) / cellx**j
#endif /** PRECISION_M2L_X **/
!
       enddo
    enddo

    call run_fmm_ewald(lm_length, winput, woutput)

    !*** Scaling local expansion coefficient by unit cell length at highest level
    do j=0,nmax
       do k=0,j
          m1 = translate_l_m_to_1dim(j, k)
#if defined(PRECISION_M2L_MIX) || defined(PRECISION_M2L_SP)
            wl_highest(m1,1) = real(woutput(m1)) / (cellx*1d10)**(j+1)  !! note real -> sp/db precision
            wl_highest(m1,2) = imag(woutput(m1)) / (cellx*1d10)**(j+1)  !! note imag -> sp/db precision
#else
            wl_highest(m1,1) = dreal(woutput(m1)) / cellx**(j+1)  !! note dreal -> dble precision
            wl_highest(m1,2) = dimag(woutput(m1)) / cellx**(j+1)  !! note dimag -> dble precision
#endif /** PRECISION_M2L_X **/
       enddo
    enddo
  end subroutine run_scaled_fmm_ewald

!---------------------------------------------------------------------
!>
!! \brief Subroutine which calculate force and potential energy from
!!        local expansion coefficients
!! \author Noriyuki Yoshii, Kensuke Iwahashi
!<
!---------------------------------------------------------------------
  subroutine energy_fmm()  !! FMM Ewald, M2L, L2L
!---------------------------------------------------------------------
    use fmm_ewald
    use domain, only : nlevel
    use fmm_ewald
    use nonneutral 
    use comm_fmm_mod, only : &
   &        sclist_xps,sclist_xms,sclist_xpr,sclist_xmr, &
   &        sclist_yps,sclist_yms,sclist_ypr,sclist_ymr, &
   &        sclist_zps,sclist_zms,sclist_zpr,sclist_zmr
    use comm_fmm_mod, only : comm_fmm_local_super
    integer(kind=wip) :: nomp
    integer(kind=wip) :: ilevel
    integer(kind=wip) :: il_spec     ! /* FOR SIMD_BPREC TEST ONLY. */
    integer(kind=wip) :: itr         ! /* FOR SIMD_BPREC TEST ONLY. */
    integer(kind=wip) :: itrs        ! /* FOR SIMD_BPREC TEST ONLY. */

    nomp = 1
!$  nomp = omp_get_max_threads()

    do ilevel=0, nlevel
       wl(ilevel)%array=Cwvzero
    enddo

#ifdef DEBUGFCE
       bNonNeutrality=.true.
#endif
    TIME_START(TM_BG_ENE)
    if(bNonNeutrality)  call background_charge_non_neutrality
!   if(bQcalc1) call init_fmm_ewald_charged(nmax, lm_length) ! only once call
    TIME_STOP(TM_BG_ENE)

    !******* multipole to local translation
    TIME_START(TM_FMM_EWALD)
    call run_scaled_fmm_ewald( &
   &                      wm(nlevel)%array(:,:,1,1,1), wl(nlevel)%array(:,:,1,1,1))
    TIME_STOP(TM_FMM_EWALD)

#ifdef DEBUG_MTDFMM
    write(myrank+30000,*) wl(nlevel)%array(:,1:2,1,1,1)
    call flush(myrank+30000)
#endif /**DEBUG_MTDFMM**/

    DO ilevel = nlevel-1, 0, -1
       if(ilevel .le. lgflg) then
          TIME_START(TM_M2L)
          TIMED_START(TM_M2L_LEVEL_OFFSET+ilevel)  !! TIMED lines closed, because gfortran can't compile them.

          call M2L_lower_level_dr( &
               &     ilevel,wm(ilevel)%array,wl(ilevel)%array, &
               &     nscdvz(ilevel),nscdvy(ilevel),nscdvx(ilevel), &
               &     nscllz(ilevel),nsclly(ilevel),nscllx(ilevel), &
               &     mbd_z_lvl(ilevel), mbd_y_lvl(ilevel), mbd_x_lvl(ilevel), &
               &     nomp)

          TIMED_STOP(TM_M2L_LEVEL_OFFSET+ilevel)
          TIME_STOP(TM_M2L)
          TIME_BARRIER(TMB_COMM_FMM_SUPER)
          TIME_START(TM_COMM_FMM_SUPER)
          call comm_fmm_local_super( &
               &     ilevel,lm_length,nsbsize, wl(ilevel)%array, &
               &     sclist_zps(ilevel)%iarray, sclist_zms(ilevel)%iarray, &
               &     sclist_zpr(ilevel)%iarray, sclist_zmr(ilevel)%iarray, &
               &     sclist_yps(ilevel)%iarray, sclist_yms(ilevel)%iarray, &
               &     sclist_ypr(ilevel)%iarray, sclist_ymr(ilevel)%iarray, &
               &     sclist_xps(ilevel)%iarray, sclist_xms(ilevel)%iarray, &
               &     sclist_xpr(ilevel)%iarray, sclist_xmr(ilevel)%iarray, &
               &     nscdvz(ilevel),nscdvy(ilevel),nscdvx(ilevel), &
               &     mcells_z(ilevel),mcells_y(ilevel),mcells_x(ilevel), &
               &     mbd_z_lvl(ilevel), mbd_y_lvl(ilevel), mbd_x_lvl(ilevel))
          TIME_STOP(TM_COMM_FMM_SUPER)
               
       else

          TIME_START(TM_M2L)
          TIMED_START(TM_M2L_LEVEL_OFFSET+ilevel)
          call M2L_upper_level(ilevel, wm(ilevel)%array, wl(ilevel)%array, nomp)
          TIMED_STOP(TM_M2L_LEVEL_OFFSET+ilevel)
          TIME_STOP(TM_M2L)
 
       endif
    ENDDO ! ilevel

    TIME_START(TM_L2L)
    DO ilevel = nlevel, 1, -1
       if(    ilevel-1.gt.lgflg)then
         mbd_uz=0; mbd_uy=0; mbd_ux=0 
         mbd_lz=0; mbd_ly=0; mbd_lx=0
       elseif(ilevel-1.eq.lgflg)then
         mbd_uz=0; mbd_uy=0; mbd_ux=0
         mbd_lz=mbd_z_lvl(ilevel-1)
         mbd_ly=mbd_y_lvl(ilevel-1)
         mbd_lx=mbd_x_lvl(ilevel-1)
       else
         mbd_uz=mbd_z_lvl(ilevel)
         mbd_uy=mbd_y_lvl(ilevel)
         mbd_ux=mbd_x_lvl(ilevel)
         mbd_lz=mbd_z_lvl(ilevel-1)
         mbd_ly=mbd_y_lvl(ilevel-1)
         mbd_lx=mbd_x_lvl(ilevel-1)
       endif
       call L2L(ilevel,wl(ilevel)%array,wl(ilevel-1)%array)
    ENDDO ! ilevel
    TIME_STOP(TM_L2L)

    TIME_START(TM_L2P)
    call L2P(wl(0)%array) ! calculate_fmm_potential_force
    TIME_STOP(TM_L2P)

#ifdef DEBUG_MTDFMM
call mpiend()  
#endif
  end subroutine energy_fmm
!---------------------------------------------------------------------
!>
!! \brief Subroutine which calculate L2L transformation
!! \author Noriyuki Yoshii, Kensuke Iwahashi
!<
  subroutine L2L(ilevel,wl,wl_lower)
!---------------------------------------------------------------------
    use domain, only : ixmin, iymin, izmin, ixmax, iymax, izmax, lxdiv, lydiv, lzdiv
    use fmm_l_m_index
    implicit none
    integer(4), intent(in) :: ilevel
    real(kind=wvap), intent(in)    :: wl(lm_length,2,          &
                                1-mbd_uz:wl_size_z(ilevel)+mbd_uz, &
                                1-mbd_uy:wl_size_y(ilevel)+mbd_uy, &
                                1-mbd_ux:wl_size_x(ilevel)+mbd_ux)
    real(kind=wvap), intent(inout) :: wl_lower(lm_length,2,          &
                                1-mbd_lz:wl_size_z(ilevel-1)+mbd_lz, &
                                1-mbd_ly:wl_size_y(ilevel-1)+mbd_ly, &
                                1-mbd_lx:wl_size_x(ilevel-1)+mbd_lx)
    integer(kind=wip) :: izst,iyst,ixst,izen,iyen,ixen
    integer(kind=wip) :: izsts,iysts,ixsts
    integer(kind=wip) :: ixyzsize,iysize,izsize
    integer(kind=wip) :: iall
    integer(kind=wip) :: icx0, icy0, icz0
    integer(kind=wip) :: icx00, icy00, icz00
    integer(kind=wip) :: icx, icy, icz
    integer(kind=wip) :: j, k, l, m
    integer(kind=wip) :: jxp, jyp, jzp, ic, jc, kc
    integer(kind=wip) :: m1, m2
    
    izsts =(izmin-1)/mcells_z(ilevel)+1
    iysts =(iymin-1)/mcells_y(ilevel)+1
    ixsts =(ixmin-1)/mcells_x(ilevel)+1
    
    iysize = (lydiv - 1) / mcells_y(ilevel) + 1
    izsize = (lzdiv - 1) / mcells_z(ilevel) + 1
    ixyzsize= ((lxdiv - 1) / mcells_x(ilevel) + 1) *iysize*izsize

    izst=(izmin-1)/mcells_z(ilevel-1)+1
    iyst=(iymin-1)/mcells_y(ilevel-1)+1
    ixst=(ixmin-1)/mcells_x(ilevel-1)+1
    izen=(izmax-1)/mcells_z(ilevel-1)+1
    iyen=(iymax-1)/mcells_y(ilevel-1)+1
    ixen=(ixmax-1)/mcells_x(ilevel-1)+1
    
!$omp parallel default(none) &
!$omp& private(iall,icz0,icy0,icx0,icz00,icy00,icx00) &
!$omp& private(icx,icy,icz) &
!$omp& private(jxp,jyp,jzp,ic,jc,kc,m1,m2) &
!$omp& private(j,k,l,m) &
!$omp& shared(ixyzsize,izsize,iysize,izsts,ixsts,iysts,izst,iyst,ixst) &
!$omp& shared(nmax,ilevel,nlx,nly,nlz,ixen,iyen,izen) &
!$omp& shared(wl_lower,wl,shll1,shll2)
!$omp do
    do iall=0,ixyzsize-1
       icz0 = mod(iall,izsize)+izsts
       icy0 = iall/izsize
       icx0 = icy0/iysize+ixsts
       icy0 = mod(icy0,iysize)+iysts
       
       icx00=icx0-ixsts+1 ! global -> local FMM L2L
       icy00=icy0-iysts+1 ! global -> local FMM L2L
       icz00=icz0-izsts+1 ! global -> local FMM L2L
       
       do jxp = 0, nlx(ilevel-1)-1
          do jyp = 0, nly(ilevel-1)-1
             do jzp = 0, nlz(ilevel-1)-1
                ic = icx0*nlx(ilevel-1)-(nlx(ilevel-1)-1) + jxp
                jc = icy0*nly(ilevel-1)-(nly(ilevel-1)-1) + jyp
                kc = icz0*nlz(ilevel-1)-(nlz(ilevel-1)-1) + jzp
                
                if(ic.ge.ixst.and.ic.le.ixen.and. &
                     &       jc.ge.iyst.and.jc.le.iyen.and. &
                     &       kc.ge.izst.and.kc.le.izen &
                     &      )then
                   icx=ic-ixst+1 ! global (ilevel-1) -> local FMM L2L (ilevel-1)
                   icy=jc-iyst+1 ! global (ilevel-1) -> local FMM L2L (ilevel-1)
                   icz=kc-izst+1 ! global (ilevel-1) -> local FMM L2L (ilevel-1)

                   do j=0,nmax; do k=0,j
                      m1 = translate_l_m_to_1dim(j, k)
                      do l=j,nmax; do m=max(0,-l+j+k), l-j+k
                         m2 = translate_l_m_to_1dim(l, m)
                         wl_lower(m1,1,icz,icy,icx) = wl_lower(m1,1,icz,icy,icx) &
                              & + wl(m2,1,icz00,icy00,icx00)*real(shll1(m2,m1,jxp,jyp,jzp,ilevel)) &
                              & + wl(m2,2,icz00,icy00,icx00)*real(shll2(m2,m1,jxp,jyp,jzp,ilevel))
                         wl_lower(m1,2,icz,icy,icx) = wl_lower(m1,2,icz,icy,icx) &
                              & + wl(m2,1,icz00,icy00,icx00)*imag(shll1(m2,m1,jxp,jyp,jzp,ilevel)) &
                              & + wl(m2,2,icz00,icy00,icx00)*imag(shll2(m2,m1,jxp,jyp,jzp,ilevel))
                      enddo; enddo
                   enddo; enddo
                endif
                
             enddo ! jzp = 0, 1
          enddo ! jyp = 0, 1
       enddo ! jxp = 0, 1
    enddo ! iall=0,ixyzsize-1
!$omp end do nowait
!$omp end parallel

#ifdef DEBUG_MTDFMM
    do iall=0,ixyzsize-1
       icz0 = mod(iall,izsize)+izsts
       icy0 = iall/izsize
       icx0 = icy0/iysize+ixsts
       icy0 = mod(icy0,iysize)+iysts
       do jxp = 0, nlx(ilevel-1)-1
          do jyp = 0, nly(ilevel-1)-1
             do jzp = 0, nlz(ilevel-1)-1
                ic = icx0*nlx(ilevel-1)-(nlx(ilevel-1)-1) + jxp
                jc = icy0*nly(ilevel-1)-(nly(ilevel-1)-1) + jyp
                kc = icz0*nlz(ilevel-1)-(nlz(ilevel-1)-1) + jzp
                if(ic.ge.ixst.and.ic.le.ixen.and. &
                     &       jc.ge.iyst.and.jc.le.iyen.and. &
                     &       kc.ge.izst.and.kc.le.izen &
                     &      )then
                   icx=ic-ixst+1 ! global (ilevel-1) -> local FMM L2L (ilevel-1)
                   icy=jc-iyst+1 ! global (ilevel-1) -> local FMM L2L (ilevel-1)
                   icz=kc-izst+1 ! global (ilevel-1) -> local FMM L2L (ilevel-1)
                   write(myrank+ilevel*1000+60000,111) ic,jc,kc
                   do j=0,nmax; do k=0,j
                      m1 = translate_l_m_to_1dim(j, k)
                   write(myrank+ilevel*1000+60000,  *) wl_lower(m1,1:2,icz,icy,icx)
                   enddo; enddo
111 format(3i5)
222 format(2d23.15)
                endif
             enddo ! jzp = 0, 1
          enddo ! jyp = 0, 1
       enddo ! jxp = 0, 1
    enddo
    call flush(myrank+ilevel*1000+60000)
#endif /**DEBUG_MTDFMM**/

  end subroutine L2L
!---------------------------------------------------------------------
!>
!! \brief Subroutine which calculate L2P transformation
!! \author Noriyuki Yoshii, Kensuke Iwahashi
!<
  subroutine L2P(wl_lowest) ! calculate_fmm_potential_force
!---------------------------------------------------------------------
    use trajectory_mpi, only : wkxyz
    use coulomb_mod, only : chgv
    use domain, only : ncellx, ncelly, ncellz, lxdiv, lydiv, lzdiv, ixmin, iymin, izmin
    use subcell, only : tag, na_per_cell, m2i
    use forces, only : w3_f
    use md_monitors, only : wk_p_energy
    use md_const
    use param, only : paranum
    use atom_virial
    use unit_cell, only : cellx, celly, cellz
    use regular_solid_harmonics_cartesian
    use dr_cntl, only : nbd
    implicit none
    real(kind=wvap), intent(in) &
                           :: wl_lowest(lm_length,2, &
   &                          1-mbd_lz:lzdiv+mbd_lz,   &
   &                          1-mbd_ly:lydiv+mbd_ly,   &
   &                          1-mbd_lx:lxdiv+mbd_lx)
    real(kind=dp) :: fxfmm, fyfmm, fzfmm
    integer(kind=wip) :: i,ii,i0,ipar
    integer(kind=wip) :: icx0, icy0, icz0
    integer(kind=wip) :: icx00, icy00, icz00
    integer(kind=wip) :: iam,nomp
    integer(kind=wip) :: m1
    real(dp) :: pot1, qta, xta, yta, zta
    real(dp) :: pot_nocoeff
    real(dp) :: R_re(lm_length), R_im(lm_length)
    complex(dp) :: f_coeff(3,lm_length)
    real(dp) :: l_re, l_im
    real(dp) :: factorial_array(0:2*nmax)
    real(dp) :: minus_double_factorial_array(0:nmax)
    real(dp) :: x_center, y_center, z_center

    !debug
    real(dp)::pot1sum
    real(dp),allocatable :: pot1omp(:)

    integer(kind=wip)::ierr
    include 'mpif.h'

    nomp = 1
    iam = 0
!$  nomp = omp_get_max_threads()

    pot1 = 0.0d0

    allocate(pot1omp(0:nomp-1))
    pot1omp=0d0

    pot_nocoeff = 0d0
    
    call prepare_factorial_arrays(nmax, factorial_array, minus_double_factorial_array)

!$omp parallel default(none) &
!$omp& private(ii,iam,i,i0,ipar) &
!$omp& private(icx0,icy0,icz0) &
!$omp& private(icx00,icy00,icz00) &
!$omp& private(x_center, y_center, z_center) &
!$omp& private(qta,xta,yta,zta) &
!$omp& private(m1,fxfmm,fyfmm,fzfmm) &
!$omp& private(l_re, l_im) &
!$omp& private(R_re,R_im,f_coeff) &
!$omp& private(pot_nocoeff) &
!$omp& shared(dxcell, dycell, dzcell, x_center_offset, y_center_offset, z_center_offset) &
!$omp& shared(factorial_array, minus_double_factorial_array) &
!$omp& shared(nmax,lm_length,chgv,paranum,wl_lowest) &
!$omp& shared(wkxyz) &
!$omp& shared(w3_f) &
!$omp& shared(tag,na_per_cell,m2i) &
!$omp& shared(lxdiv,lydiv,lzdiv,myrank) &
!$omp& shared(pot1omp)
!$  iam = omp_get_thread_num()
    do ii=1,lxdiv*lydiv*lzdiv
       icz0=mod(ii-1,lzdiv)     +nbd 
       icy0=mod(ii-1,lzdiv*lydiv)
       icy0=icy0/lzdiv          +nbd 
       icx0=(ii-1)/(lzdiv*lydiv)+nbd 
       icx00=icx0-nbd+1             ! local -> local FMM
       icy00=icy0-nbd+1             ! local -> local FMM
       icz00=icz0-nbd+1             ! local -> local FMM
       x_center = icx0 * dxcell + x_center_offset
       y_center = icy0 * dycell + y_center_offset
       z_center = icz0 * dzcell + z_center_offset
!$omp do
       do i0=tag(icz0,icy0,icx0), &
     &       tag(icz0,icy0,icx0)+na_per_cell(icz0,icy0,icx0)-1
          i=m2i(i0)
          ipar= paranum(i)
          qta = md_QQ_4PiE*chgv(ipar)
          
          !zjc
#if defined(PRECISION_M2L_MIX) || defined(PRECISION_M2L_SP)
            xta = wkxyz(1,i0) * 1d10 - x_center
            yta = wkxyz(2,i0) * 1d10 - y_center
            zta = wkxyz(3,i0) * 1d10 - z_center
#else
            xta = wkxyz(1,i0) - x_center
            yta = wkxyz(2,i0) - y_center
            zta = wkxyz(3,i0) - z_center
#endif
          call calculate_regular_harmonics_force_1dim_array( &
               & nmax, factorial_array, minus_double_factorial_array, &
               & xta, yta, zta, R_re, R_im, f_coeff)
          !*** calculate multipole potential
          pot_nocoeff = 0d0
          do m1=1,lm_length
             pot_nocoeff = pot_nocoeff + &
                  & wl_lowest(m1,1,icz00,icy00,icx00) * R_re(m1) &
                  - wl_lowest(m1,2,icz00,icy00,icx00) * R_im(m1)
          enddo ! for m1
          pot1omp(iam) = pot1omp(iam) + qta * pot_nocoeff

          !******* calculate multipole force
          fxfmm=0.d0
          fyfmm=0.d0
          fzfmm=0.d0
          do m1=2,lm_length  ! start index 2 means (l,m)=(1,0)
             l_re = wl_lowest(m1,1,icz00,icy00,icx00)
             l_im = wl_lowest(m1,2,icz00,icy00,icx00)

             fxfmm = fxfmm + l_re * real(f_coeff(1,m1)) + l_im * imag(f_coeff(1,m1))
             fyfmm = fyfmm + l_re * real(f_coeff(2,m1)) + l_im * imag(f_coeff(2,m1))
             fzfmm = fzfmm + l_re * real(f_coeff(3,m1)) + l_im * imag(f_coeff(3,m1))
          enddo ! m1

#if defined(PRECISION_M2L_MIX) || defined(PRECISION_M2L_SP)
          fxfmm = fxfmm * 1d20
          fyfmm = fyfmm * 1d20
          fzfmm = fzfmm * 1d20
#endif

          w3_f(1,i0,iam) = w3_f(1,i0,iam)+       qta * fxfmm
          w3_f(2,i0,iam) = w3_f(2,i0,iam)+       qta * fyfmm
          w3_f(3,i0,iam) = w3_f(3,i0,iam)+ 2d0 * qta * fzfmm

       enddo ! i0
!$omp end do nowait
    enddo ! ii
!$omp end parallel
    pot1=0d0
    do iam=0,nomp-1
       pot1=pot1+pot1omp(iam)
    enddo

    deallocate(pot1omp)

#if defined(PRECISION_M2L_MIX) || defined(PRECISION_M2L_SP)
    pot1 = pot1 * 1d10
#endif

    wk_p_energy = wk_p_energy + pot1

    wk_vir2(1,0)=wk_vir2(1,0)+pot1/3d0
    wk_vir2(2,0)=wk_vir2(2,0)+pot1/3d0
    wk_vir2(3,0)=wk_vir2(3,0)+pot1/3d0

#ifdef DEBUGFCE
    call mpi_allreduce(pot1,pot1sum,1, &
         &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
#ifdef KCAL
    if(myrank==0) write(*,*) 'Pot(FMM  )=', pot1sum*kJ_mol/4.184d0,'[kcal/mol]'
#else
    if(myrank==0) write(*,*) 'Pot(FMM  )=', pot1sum*kJ_mol,'[kJ/mol]'
#endif
#endif /**DEBUGFCE**/

#ifdef DEBUG_MTDFMM
    call mpi_allreduce(pot1,pot1sum,1, &
         &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    if(myrank==0) write(*,*) 'Pot(FMM  )=', pot1sum/md_QQ_4PiE
#endif /**DEBUG_MTDFMM**/

  end subroutine L2P
!--------------------------------------------------------------------- 
! M2L codes: M2L_lower_level_dr/M2L_upper_level
!---------------------------------------------------------------------
! There are following branches in M2L codes for performance optimization
!
!      #ifdef MOMEMNT
!      #else
!      #endif
!
!---------------------------------------------------------------------
!>
!! \brief  Subroutine to execute M2L (with local arrays, lower level) 
!!         with real and imaginary parts explicitly separeted.
!! \author Noriyuki Yoshii, Yoshimichi Andoh, Shin-chi Ichikawa
!<
  subroutine M2L_lower_level_dr(ilevel,wm,wl, &
       &                        nsczdiv,nscydiv,nscxdiv, &
       &                        nscellz,nscelly,nscellx, &
       &                        mbd_z,  mbd_y,  mbd_x, &
       &                        nomp)
!---------------------------------------------------------------------
! Main routine of M2L_lower_level_dr, handling rotation or not, 
! nmax=4 or general, multiple cell pairs or single cell pair.
!
    implicit none
    integer(kind=wip), intent(in) :: ilevel
    real(kind=wvp), intent(in) :: wm(lm_length,2, &
         &              1-nbd_zm(ilevel):nsczdiv+nbd_zp(ilevel), &
         &              1-nbd_ym(ilevel):nscydiv+nbd_yp(ilevel), &
         &              1-nbd_xm(ilevel):nscxdiv+nbd_xp(ilevel))
    real(kind=wvap), intent(inout) :: wl(lm_length,2, &
         &              1-mbd_z         :nsczdiv+mbd_z, &
         &              1-mbd_y         :nscydiv+mbd_y, &
         &              1-mbd_x         :nscxdiv+mbd_x      )
    integer(kind=wip), intent(in) :: nsczdiv,nscydiv,nscxdiv
    integer(kind=wip), intent(in) :: nscellz, nscelly, nscellx
    integer(kind=wip), intent(in) :: mbd_z, mbd_y, mbd_x
    integer(kind=wip) :: nomp

#ifdef MOMENT
    if(nmax==4) then    !/*---- for Matrix Order 15 ----*/

      call M2L_lower_level_dr_mpairs_nmax4(ilevel,wm,wl, &
       &                        nsczdiv,nscydiv,nscxdiv, &
       &                        mbd_z,  mbd_y,  mbd_x, &
       &                        nomp)

    else                !/*---- for General Value of NMAX ----*/

     call M2L_lower_level_dr_mpairs(ilevel,wm,wl, &
       &                        nsczdiv,nscydiv,nscxdiv, &
       &                        mbd_z,  mbd_y,  mbd_x, &
       &                        nomp)

    endif               !/*---- for General Value of NMAX ----*/
#else   /* Not MOMENT */
    if(nmax==4) then    !/*---- for Matrix Order 15 ----*/
! Sakashita's code
      call M2L_lower_level_dr_nmax4(ilevel,wm,wl, &
       &                        nsczdiv,nscydiv,nscxdiv, &
       &                        nscellz,nscelly,nscellx, &
       &                        mbd_z,  mbd_y,  mbd_x, &
       &                        nomp)

    else                !/*---- for General Value of NMAX ----*/

      call M2L_lower_level_dr_general(ilevel,wm,wl, &
       &                        nsczdiv,nscydiv,nscxdiv, &
       &                        nscellz,nscelly,nscellx, &
       &                        mbd_z,  mbd_y,  mbd_x, &
       &                        nomp)

    endif
#endif  /* Not MOMENT */

  end subroutine M2L_lower_level_dr

#ifdef MOMENT
!---------------------------------------------------------------------
!>
!! \brief  Subroutine to execute M2L (with local arrays, lower level)
!! \author Noriyuki Yoshii, Yoshimichi Andoh, Shin-chi Ichikawa
!<
  subroutine M2L_lower_level_dr_mpairs_nmax4(ilevel,wm,wl, &
       &                        nsczdiv,nscydiv,nscxdiv, &
       &                        mbd_z,  mbd_y,  mbd_x, &
       &                        nomp)
!---------------------------------------------------------------------
! symmetry and rotation code.
! pipelined multiple cell pairs and unrolled for nmax=4.
! 
    use mpi_3d_grid, only : npx, npy, npz, ipx, ipy, ipz
    use domain, only : nlevel
    implicit none
    integer(kind=wip), intent(in) :: ilevel
    real(kind=wvp), intent(in) :: wm(lm_length,2, &
         &              1-nbd_zm(ilevel):nsczdiv+nbd_zp(ilevel), &
         &              1-nbd_ym(ilevel):nscydiv+nbd_yp(ilevel), &
         &              1-nbd_xm(ilevel):nscxdiv+nbd_xp(ilevel))
    real(kind=wvap), intent(inout) :: wl(lm_length,2, &
         &              1-mbd_z         :nsczdiv+mbd_z, &
         &              1-mbd_y         :nscydiv+mbd_y, &
         &              1-mbd_x         :nscxdiv+mbd_x      )
    integer(kind=wip), intent(in) :: nsczdiv,nscydiv,nscxdiv
    integer(kind=wip), intent(in) :: mbd_z, mbd_y, mbd_x
    real(kind=wvap) :: wl_omp(lm_length+5, 2, &
         &              1-mbd_z         :nsczdiv+2+mbd_z, &
         &              1-mbd_y         :nscydiv  +mbd_y, &
         &              1-mbd_x         :nscxdiv  +mbd_x, 0:nomp-1)
    integer(kind=wip) :: icx0, icy0, icz0
    integer(kind=wip) :: icx1, icy1, icz1
    integer(kind=wip) :: nomp
    integer(kind=wip) :: nt, iam
    integer(kind=wip) :: icxyz
    integer(kind=wip) :: mbd_z2, mbd_y2, mbd_x2
    integer(kind=wip) :: lcbx, lcby, lcbz
    integer(kind=wip) :: irx, iry, irz
    integer(kind=wip) :: ncall, ncblk
    real(kind=wvp) :: wwm_north(lm_length,2,8), wwl_north(lm_length,2,8)
!add for simultaneous multiple cell pairs.
    integer(kind=wip) :: icp
    integer(kind=wip) :: icx0_1,icy0_1,icz0_1,icx1_1,icy1_1,icz1_1
    integer(kind=wip) :: icx0_2,icy0_2,icz0_2,icx1_2,icy1_2,icz1_2
    integer(kind=wip) :: icx0_3,icy0_3,icz0_3,icx1_3,icy1_3,icz1_3
    integer(kind=wip) :: icx0_4,icy0_4,icz0_4,icx1_4,icy1_4,icz1_4
    integer(kind=wip) :: icx0_5,icy0_5,icz0_5,icx1_5,icy1_5,icz1_5
    integer(kind=wip) :: icx0_6,icy0_6,icz0_6,icx1_6,icy1_6,icz1_6
    integer(kind=wip) :: icx0_7,icy0_7,icz0_7,icx1_7,icy1_7,icz1_7
    integer(kind=wip) :: icx0_8,icy0_8,icz0_8,icx1_8,icy1_8,icz1_8
    real(kind=wvp) :: shml_n0, shml_n1, shml_n2, shml_n3, shml_n4
    real(kind=wvp) :: shrot1_q0_r, shrot1_q0_i, shrot2_q0_r, shrot2_q0_i
    real(kind=wvp) :: shrot1_q1_r, shrot1_q1_i, shrot2_q1_r, shrot2_q1_i
    real(kind=wvp) :: shrot1_q2_r, shrot1_q2_i, shrot2_q2_r, shrot2_q2_i
    real(kind=wvp) :: shrot1_q3_r, shrot1_q3_i, shrot2_q3_r, shrot2_q3_i
    real(kind=wvp) :: shrot1_q4_r, shrot1_q4_i, shrot2_q4_r, shrot2_q4_i
    real(kind=wvp) :: shinvrot1_q0_r, shinvrot1_q0_i, shinvrot2_q0_r, shinvrot2_q0_i
    real(kind=wvp) :: shinvrot1_q1_r, shinvrot1_q1_i, shinvrot2_q1_r, shinvrot2_q1_i
    real(kind=wvp) :: shinvrot1_q2_r, shinvrot1_q2_i, shinvrot2_q2_r, shinvrot2_q2_i
    real(kind=wvp) :: shinvrot1_q3_r, shinvrot1_q3_i, shinvrot2_q3_r, shinvrot2_q3_i
    real(kind=wvp) :: shinvrot1_q4_r, shinvrot1_q4_i, shinvrot2_q4_r, shinvrot2_q4_i
    real(kind=wvp) :: wm_q0_1r, wm_q0_1i, wm_q0_2r, wm_q0_2i, wm_q0_3r, wm_q0_3i, wm_q0_4r, wm_q0_4i
    real(kind=wvp) :: wm_q1_1r, wm_q1_1i, wm_q1_2r, wm_q1_2i, wm_q1_3r, wm_q1_3i, wm_q1_4r, wm_q1_4i
    real(kind=wvp) :: wm_q2_1r, wm_q2_1i, wm_q2_2r, wm_q2_2i, wm_q2_3r, wm_q2_3i, wm_q2_4r, wm_q2_4i
    real(kind=wvp) :: wm_q3_1r, wm_q3_1i, wm_q3_2r, wm_q3_2i, wm_q3_3r, wm_q3_3i, wm_q3_4r, wm_q3_4i
    real(kind=wvp) :: wm_q4_1r, wm_q4_1i, wm_q4_2r, wm_q4_2i, wm_q4_3r, wm_q4_3i, wm_q4_4r, wm_q4_4i
    real(kind=wvp) :: wm_q0_5r, wm_q0_5i, wm_q0_6r, wm_q0_6i, wm_q0_7r, wm_q0_7i, wm_q0_8r, wm_q0_8i
    real(kind=wvp) :: wm_q1_5r, wm_q1_5i, wm_q1_6r, wm_q1_6i, wm_q1_7r, wm_q1_7i, wm_q1_8r, wm_q1_8i
    real(kind=wvp) :: wm_q2_5r, wm_q2_5i, wm_q2_6r, wm_q2_6i, wm_q2_7r, wm_q2_7i, wm_q2_8r, wm_q2_8i
    real(kind=wvp) :: wm_q3_5r, wm_q3_5i, wm_q3_6r, wm_q3_6i, wm_q3_7r, wm_q3_7i, wm_q3_8r, wm_q3_8i
    real(kind=wvp) :: wm_q4_5r, wm_q4_5i, wm_q4_6r, wm_q4_6i, wm_q4_7r, wm_q4_7i, wm_q4_8r, wm_q4_8i
    real(kind=wvp) :: wlnorth_q0_1r, wlnorth_q0_1i, wlnorth_q0_2r, wlnorth_q0_2i, &
   &                  wlnorth_q0_3r, wlnorth_q0_3i, wlnorth_q0_4r, wlnorth_q0_4i
    real(kind=wvp) :: wlnorth_q1_1r, wlnorth_q1_1i, wlnorth_q1_2r, wlnorth_q1_2i, &
   &                  wlnorth_q1_3r, wlnorth_q1_3i, wlnorth_q1_4r, wlnorth_q1_4i
    real(kind=wvp) :: wlnorth_q2_1r, wlnorth_q2_1i, wlnorth_q2_2r, wlnorth_q2_2i, &
   &                  wlnorth_q2_3r, wlnorth_q2_3i, wlnorth_q2_4r, wlnorth_q2_4i
    real(kind=wvp) :: wlnorth_q3_1r, wlnorth_q3_1i, wlnorth_q3_2r, wlnorth_q3_2i, &
   &                  wlnorth_q3_3r, wlnorth_q3_3i, wlnorth_q3_4r, wlnorth_q3_4i
    real(kind=wvp) :: wlnorth_q4_1r, wlnorth_q4_1i, wlnorth_q4_2r, wlnorth_q4_2i, &
   &                  wlnorth_q4_3r, wlnorth_q4_3i, wlnorth_q4_4r, wlnorth_q4_4i
    real(kind=wvp) :: wlnorth_q0_5r, wlnorth_q0_5i, wlnorth_q0_6r, wlnorth_q0_6i, &
   &                  wlnorth_q0_7r, wlnorth_q0_7i, wlnorth_q0_8r, wlnorth_q0_8i
    real(kind=wvp) :: wlnorth_q1_5r, wlnorth_q1_5i, wlnorth_q1_6r, wlnorth_q1_6i, &
   &                  wlnorth_q1_7r, wlnorth_q1_7i, wlnorth_q1_8r, wlnorth_q1_8i
    real(kind=wvp) :: wlnorth_q2_5r, wlnorth_q2_5i, wlnorth_q2_6r, wlnorth_q2_6i, &
   &                  wlnorth_q2_7r, wlnorth_q2_7i, wlnorth_q2_8r, wlnorth_q2_8i
    real(kind=wvp) :: wlnorth_q3_5r, wlnorth_q3_5i, wlnorth_q3_6r, wlnorth_q3_6i, &
   &                  wlnorth_q3_7r, wlnorth_q3_7i, wlnorth_q3_8r, wlnorth_q3_8i
    real(kind=wvp) :: wlnorth_q4_5r, wlnorth_q4_5i, wlnorth_q4_6r, wlnorth_q4_6i, &
   &                  wlnorth_q4_7r, wlnorth_q4_7i, wlnorth_q4_8r, wlnorth_q4_8i
    real(kind=wvp) :: tmp_1r, tmp_1i, tmp_2r, tmp_2i, tmp_3r, tmp_3i, tmp_4r, tmp_4i
    real(kind=wvp) :: tmp_5r, tmp_5i, tmp_6r, tmp_6i, tmp_7r, tmp_7i, tmp_8r, tmp_8i
    real(kind=wvp) :: wm_north(lm_length,2), wl_north(lm_length,2)
    real(kind=wvp) :: tmp_shml, tmp_ce_r, tmp_ce_i
    integer(kind=wip) :: k,p, m1

    iam = 0

    wl_omp = Cdzero

    ! global address offset of a cell prior to starting cell.
    mbd_z2 = 2*mbd_z
    mbd_y2 = 2*mbd_y
    mbd_x2 = 2*mbd_x

!$omp parallel default(none) &
!$omp& private(iam) &
!$omp& private(icx0,icy0,icz0,icx1,icy1,icz1) &
!$omp& private(irx,iry,irz) &
!$omp& private(wwm_north,wwl_north) &
!$omp& private(icp) &
!$omp& private(icx0_1,icy0_1,icz0_1,icx1_1,icy1_1,icz1_1) &
!$omp& private(icx0_2,icy0_2,icz0_2,icx1_2,icy1_2,icz1_2) &
!$omp& private(icx0_3,icy0_3,icz0_3,icx1_3,icy1_3,icz1_3) &
!$omp& private(icx0_4,icy0_4,icz0_4,icx1_4,icy1_4,icz1_4) &
!$omp& private(icx0_5,icy0_5,icz0_5,icx1_5,icy1_5,icz1_5) &
!$omp& private(icx0_6,icy0_6,icz0_6,icx1_6,icy1_6,icz1_6) &
!$omp& private(icx0_7,icy0_7,icz0_7,icx1_7,icy1_7,icz1_7) &
!$omp& private(icx0_8,icy0_8,icz0_8,icx1_8,icy1_8,icz1_8) &
!$omp& private(shml_n0, shml_n1, shml_n2, shml_n3, shml_n4) &
!$omp& private(shrot1_q0_r, shrot1_q0_i, shrot2_q0_r, shrot2_q0_i) &
!$omp& private(shrot1_q1_r, shrot1_q1_i, shrot2_q1_r, shrot2_q1_i) &
!$omp& private(shrot1_q2_r, shrot1_q2_i, shrot2_q2_r, shrot2_q2_i) &
!$omp& private(shrot1_q3_r, shrot1_q3_i, shrot2_q3_r, shrot2_q3_i) &
!$omp& private(shrot1_q4_r, shrot1_q4_i, shrot2_q4_r, shrot2_q4_i) &
!$omp& private(shinvrot1_q0_r, shinvrot1_q0_i, shinvrot2_q0_r, shinvrot2_q0_i) &
!$omp& private(shinvrot1_q1_r, shinvrot1_q1_i, shinvrot2_q1_r, shinvrot2_q1_i) &
!$omp& private(shinvrot1_q2_r, shinvrot1_q2_i, shinvrot2_q2_r, shinvrot2_q2_i) &
!$omp& private(shinvrot1_q3_r, shinvrot1_q3_i, shinvrot2_q3_r, shinvrot2_q3_i) &
!$omp& private(shinvrot1_q4_r, shinvrot1_q4_i, shinvrot2_q4_r, shinvrot2_q4_i) &
!$omp& private(wm_q0_1r, wm_q0_1i, wm_q0_2r, wm_q0_2i, wm_q0_3r, wm_q0_3i, wm_q0_4r, wm_q0_4i) &
!$omp& private(wm_q1_1r, wm_q1_1i, wm_q1_2r, wm_q1_2i, wm_q1_3r, wm_q1_3i, wm_q1_4r, wm_q1_4i) &
!$omp& private(wm_q2_1r, wm_q2_1i, wm_q2_2r, wm_q2_2i, wm_q2_3r, wm_q2_3i, wm_q2_4r, wm_q2_4i) &
!$omp& private(wm_q3_1r, wm_q3_1i, wm_q3_2r, wm_q3_2i, wm_q3_3r, wm_q3_3i, wm_q3_4r, wm_q3_4i) &
!$omp& private(wm_q4_1r, wm_q4_1i, wm_q4_2r, wm_q4_2i, wm_q4_3r, wm_q4_3i, wm_q4_4r, wm_q4_4i) &
!$omp& private(wm_q0_5r, wm_q0_5i, wm_q0_6r, wm_q0_6i, wm_q0_7r, wm_q0_7i, wm_q0_8r, wm_q0_8i) &
!$omp& private(wm_q1_5r, wm_q1_5i, wm_q1_6r, wm_q1_6i, wm_q1_7r, wm_q1_7i, wm_q1_8r, wm_q1_8i) &
!$omp& private(wm_q2_5r, wm_q2_5i, wm_q2_6r, wm_q2_6i, wm_q2_7r, wm_q2_7i, wm_q2_8r, wm_q2_8i) &
!$omp& private(wm_q3_5r, wm_q3_5i, wm_q3_6r, wm_q3_6i, wm_q3_7r, wm_q3_7i, wm_q3_8r, wm_q3_8i) &
!$omp& private(wm_q4_5r, wm_q4_5i, wm_q4_6r, wm_q4_6i, wm_q4_7r, wm_q4_7i, wm_q4_8r, wm_q4_8i) &
!$omp& private(wlnorth_q0_1r, wlnorth_q0_1i, wlnorth_q0_2r, wlnorth_q0_2i) &
!$omp& private(wlnorth_q0_3r, wlnorth_q0_3i, wlnorth_q0_4r, wlnorth_q0_4i) &
!$omp& private(wlnorth_q1_1r, wlnorth_q1_1i, wlnorth_q1_2r, wlnorth_q1_2i) &
!$omp& private(wlnorth_q1_3r, wlnorth_q1_3i, wlnorth_q1_4r, wlnorth_q1_4i) &
!$omp& private(wlnorth_q2_1r, wlnorth_q2_1i, wlnorth_q2_2r, wlnorth_q2_2i) &
!$omp& private(wlnorth_q2_3r, wlnorth_q2_3i, wlnorth_q2_4r, wlnorth_q2_4i) &
!$omp& private(wlnorth_q3_1r, wlnorth_q3_1i, wlnorth_q3_2r, wlnorth_q3_2i) &
!$omp& private(wlnorth_q3_3r, wlnorth_q3_3i, wlnorth_q3_4r, wlnorth_q3_4i) &
!$omp& private(wlnorth_q4_1r, wlnorth_q4_1i, wlnorth_q4_2r, wlnorth_q4_2i) &
!$omp& private(wlnorth_q4_3r, wlnorth_q4_3i, wlnorth_q4_4r, wlnorth_q4_4i) &
!$omp& private(wlnorth_q0_5r, wlnorth_q0_5i, wlnorth_q0_6r, wlnorth_q0_6i) &
!$omp& private(wlnorth_q0_7r, wlnorth_q0_7i, wlnorth_q0_8r, wlnorth_q0_8i) &
!$omp& private(wlnorth_q1_5r, wlnorth_q1_5i, wlnorth_q1_6r, wlnorth_q1_6i) &
!$omp& private(wlnorth_q1_7r, wlnorth_q1_7i, wlnorth_q1_8r, wlnorth_q1_8i) &
!$omp& private(wlnorth_q2_5r, wlnorth_q2_5i, wlnorth_q2_6r, wlnorth_q2_6i) &
!$omp& private(wlnorth_q2_7r, wlnorth_q2_7i, wlnorth_q2_8r, wlnorth_q2_8i) &
!$omp& private(wlnorth_q3_5r, wlnorth_q3_5i, wlnorth_q3_6r, wlnorth_q3_6i) &
!$omp& private(wlnorth_q3_7r, wlnorth_q3_7i, wlnorth_q3_8r, wlnorth_q3_8i) &
!$omp& private(wlnorth_q4_5r, wlnorth_q4_5i, wlnorth_q4_6r, wlnorth_q4_6i) &
!$omp& private(wlnorth_q4_7r, wlnorth_q4_7i, wlnorth_q4_8r, wlnorth_q4_8i) &
!$omp& private(tmp_1r, tmp_1i, tmp_2r, tmp_2i, tmp_3r, tmp_3i, tmp_4r, tmp_4i) &
!$omp& private(tmp_5r, tmp_5i, tmp_6r, tmp_6i, tmp_7r, tmp_7i, tmp_8r, tmp_8i) &
!$omp& private(wm_north,wl_north) &
!$omp& private(k,p, m1) &
!$omp& private(tmp_shml, tmp_ce_r, tmp_ce_i) &
!$omp& shared(shml, sh_rot1, sh_rot2, sh_inv_rot1, sh_inv_rot2) &
!$omp& shared(ilevel,lm_length) &
!$omp& shared(wl_omp,wm) &
!$omp& shared(cp8_count,cp6_count,cp4_count,cp3_count,cp2_count,cp1_count) &
!$omp& shared(cp8_list0,cp6_list0,cp4_list0,cp3_list0,cp2_list0,cp1_list0) &
!$omp& shared(cp8_list1,cp6_list1,cp4_list1,cp3_list1,cp2_list1,cp1_list1)
!$  iam = omp_get_thread_num()

    if(cp8_count(ilevel) > 0) then    !! /*------- pairs multiple of 8 for each relative cell address -------*/

!$omp do
    do icp = 1, cp8_count(ilevel), 8
      irx = cp8_list1(1,icp,ilevel) - cp8_list0(1,icp,ilevel)  ! x-dir relative cell address of transformation matrix.
      iry = cp8_list1(2,icp,ilevel) - cp8_list0(2,icp,ilevel)  ! y-dir relative cell address of transformation matrix.
      irz = cp8_list1(3,icp,ilevel) - cp8_list0(3,icp,ilevel)  ! z-dir relative cell address of transformation matrix.

    !*** multipole to local translation
      icx0_1 = cp8_list0(1,icp  ,ilevel) ; icy0_1 = cp8_list0(2,icp  ,ilevel) ; icz0_1 = cp8_list0(3,icp  ,ilevel)
      icx1_1 = cp8_list1(1,icp  ,ilevel) ; icy1_1 = cp8_list1(2,icp  ,ilevel) ; icz1_1 = cp8_list1(3,icp  ,ilevel)
      icx0_2 = cp8_list0(1,icp+1,ilevel) ; icy0_2 = cp8_list0(2,icp+1,ilevel) ; icz0_2 = cp8_list0(3,icp+1,ilevel)
      icx1_2 = cp8_list1(1,icp+1,ilevel) ; icy1_2 = cp8_list1(2,icp+1,ilevel) ; icz1_2 = cp8_list1(3,icp+1,ilevel)
      icx0_3 = cp8_list0(1,icp+2,ilevel) ; icy0_3 = cp8_list0(2,icp+2,ilevel) ; icz0_3 = cp8_list0(3,icp+2,ilevel)
      icx1_3 = cp8_list1(1,icp+2,ilevel) ; icy1_3 = cp8_list1(2,icp+2,ilevel) ; icz1_3 = cp8_list1(3,icp+2,ilevel)
      icx0_4 = cp8_list0(1,icp+3,ilevel) ; icy0_4 = cp8_list0(2,icp+3,ilevel) ; icz0_4 = cp8_list0(3,icp+3,ilevel)
      icx1_4 = cp8_list1(1,icp+3,ilevel) ; icy1_4 = cp8_list1(2,icp+3,ilevel) ; icz1_4 = cp8_list1(3,icp+3,ilevel)
      icx0_5 = cp8_list0(1,icp+4,ilevel) ; icy0_5 = cp8_list0(2,icp+4,ilevel) ; icz0_5 = cp8_list0(3,icp+4,ilevel)
      icx1_5 = cp8_list1(1,icp+4,ilevel) ; icy1_5 = cp8_list1(2,icp+4,ilevel) ; icz1_5 = cp8_list1(3,icp+4,ilevel)
      icx0_6 = cp8_list0(1,icp+5,ilevel) ; icy0_6 = cp8_list0(2,icp+5,ilevel) ; icz0_6 = cp8_list0(3,icp+5,ilevel)
      icx1_6 = cp8_list1(1,icp+5,ilevel) ; icy1_6 = cp8_list1(2,icp+5,ilevel) ; icz1_6 = cp8_list1(3,icp+5,ilevel)
      icx0_7 = cp8_list0(1,icp+6,ilevel) ; icy0_7 = cp8_list0(2,icp+6,ilevel) ; icz0_7 = cp8_list0(3,icp+6,ilevel)
      icx1_7 = cp8_list1(1,icp+6,ilevel) ; icy1_7 = cp8_list1(2,icp+6,ilevel) ; icz1_7 = cp8_list1(3,icp+6,ilevel)
      icx0_8 = cp8_list0(1,icp+7,ilevel) ; icy0_8 = cp8_list0(2,icp+7,ilevel) ; icz0_8 = cp8_list0(3,icp+7,ilevel)
      icx1_8 = cp8_list1(1,icp+7,ilevel) ; icy1_8 = cp8_list1(2,icp+7,ilevel) ; icz1_8 = cp8_list1(3,icp+7,ilevel)

      include 'tuned_M2L/fmm_far_m2l_rotation_mpair8_nmax4_exsimd-reduct_unroll.h90'
    end do  ! cp8_count.
!$omp end do

    endif    ! cp8_count>0.

    if(cp6_count(ilevel) > 0) then    !! /*------- pairs multiple of 6 for each relative cell address -------*/

!$omp do
    do icp = 1, cp6_count(ilevel), 6
      irx = cp6_list1(1,icp,ilevel) - cp6_list0(1,icp,ilevel)  ! x-dir relative cell address of transformation matrix.
      iry = cp6_list1(2,icp,ilevel) - cp6_list0(2,icp,ilevel)  ! y-dir relative cell address of transformation matrix.
      irz = cp6_list1(3,icp,ilevel) - cp6_list0(3,icp,ilevel)  ! z-dir relative cell address of transformation matrix.

    !*** multipole to local translation
      icx0_1 = cp6_list0(1,icp  ,ilevel) ; icy0_1 = cp6_list0(2,icp  ,ilevel) ; icz0_1 = cp6_list0(3,icp  ,ilevel)
      icx1_1 = cp6_list1(1,icp  ,ilevel) ; icy1_1 = cp6_list1(2,icp  ,ilevel) ; icz1_1 = cp6_list1(3,icp  ,ilevel)
      icx0_2 = cp6_list0(1,icp+1,ilevel) ; icy0_2 = cp6_list0(2,icp+1,ilevel) ; icz0_2 = cp6_list0(3,icp+1,ilevel)
      icx1_2 = cp6_list1(1,icp+1,ilevel) ; icy1_2 = cp6_list1(2,icp+1,ilevel) ; icz1_2 = cp6_list1(3,icp+1,ilevel)
      icx0_3 = cp6_list0(1,icp+2,ilevel) ; icy0_3 = cp6_list0(2,icp+2,ilevel) ; icz0_3 = cp6_list0(3,icp+2,ilevel)
      icx1_3 = cp6_list1(1,icp+2,ilevel) ; icy1_3 = cp6_list1(2,icp+2,ilevel) ; icz1_3 = cp6_list1(3,icp+2,ilevel)
      icx0_4 = cp6_list0(1,icp+3,ilevel) ; icy0_4 = cp6_list0(2,icp+3,ilevel) ; icz0_4 = cp6_list0(3,icp+3,ilevel)
      icx1_4 = cp6_list1(1,icp+3,ilevel) ; icy1_4 = cp6_list1(2,icp+3,ilevel) ; icz1_4 = cp6_list1(3,icp+3,ilevel)
      icx0_5 = cp6_list0(1,icp+4,ilevel) ; icy0_5 = cp6_list0(2,icp+4,ilevel) ; icz0_5 = cp6_list0(3,icp+4,ilevel)
      icx1_5 = cp6_list1(1,icp+4,ilevel) ; icy1_5 = cp6_list1(2,icp+4,ilevel) ; icz1_5 = cp6_list1(3,icp+4,ilevel)
      icx0_6 = cp6_list0(1,icp+5,ilevel) ; icy0_6 = cp6_list0(2,icp+5,ilevel) ; icz0_6 = cp6_list0(3,icp+5,ilevel)
      icx1_6 = cp6_list1(1,icp+5,ilevel) ; icy1_6 = cp6_list1(2,icp+5,ilevel) ; icz1_6 = cp6_list1(3,icp+5,ilevel)

      include 'tuned_M2L/fmm_far_m2l_rotation_mpair6_nmax4_exsimd-reduct_unroll.h90'
    end do  ! cp6_count.
!$omp end do
    endif    ! cp6_count>0.

    if(cp4_count(ilevel) > 0) then    !! /*------- pairs multiple of 4 for each relative cell address -------*/

!$omp do
    do icp = 1, cp4_count(ilevel), 4
      irx = cp4_list1(1,icp,ilevel) - cp4_list0(1,icp,ilevel)  ! x-dir relative cell address of transformation matrix.
      iry = cp4_list1(2,icp,ilevel) - cp4_list0(2,icp,ilevel)  ! y-dir relative cell address of transformation matrix.
      irz = cp4_list1(3,icp,ilevel) - cp4_list0(3,icp,ilevel)  ! z-dir relative cell address of transformation matrix.

    !*** multipole to local translation
      icx0_1 = cp4_list0(1,icp  ,ilevel) ; icy0_1 = cp4_list0(2,icp  ,ilevel) ; icz0_1 = cp4_list0(3,icp  ,ilevel)
      icx1_1 = cp4_list1(1,icp  ,ilevel) ; icy1_1 = cp4_list1(2,icp  ,ilevel) ; icz1_1 = cp4_list1(3,icp  ,ilevel)
      icx0_2 = cp4_list0(1,icp+1,ilevel) ; icy0_2 = cp4_list0(2,icp+1,ilevel) ; icz0_2 = cp4_list0(3,icp+1,ilevel)
      icx1_2 = cp4_list1(1,icp+1,ilevel) ; icy1_2 = cp4_list1(2,icp+1,ilevel) ; icz1_2 = cp4_list1(3,icp+1,ilevel)
      icx0_3 = cp4_list0(1,icp+2,ilevel) ; icy0_3 = cp4_list0(2,icp+2,ilevel) ; icz0_3 = cp4_list0(3,icp+2,ilevel)
      icx1_3 = cp4_list1(1,icp+2,ilevel) ; icy1_3 = cp4_list1(2,icp+2,ilevel) ; icz1_3 = cp4_list1(3,icp+2,ilevel)
      icx0_4 = cp4_list0(1,icp+3,ilevel) ; icy0_4 = cp4_list0(2,icp+3,ilevel) ; icz0_4 = cp4_list0(3,icp+3,ilevel)
      icx1_4 = cp4_list1(1,icp+3,ilevel) ; icy1_4 = cp4_list1(2,icp+3,ilevel) ; icz1_4 = cp4_list1(3,icp+3,ilevel)

      include 'tuned_M2L/fmm_far_m2l_rotation_mpair4_nmax4_exsimd-reduct_unroll.h90'
    end do  ! cp4_count.
!$omp end do
    endif    ! cp4_count>0.

    if(cp3_count(ilevel) > 0) then    !! /*------- pairs multiple of 3 for each relative cell address -------*/

!$omp do
    do icp = 1, cp3_count(ilevel), 3
      irx = cp3_list1(1,icp,ilevel) - cp3_list0(1,icp,ilevel)  ! x-dir relative cell address of transformation matrix.
      iry = cp3_list1(2,icp,ilevel) - cp3_list0(2,icp,ilevel)  ! y-dir relative cell address of transformation matrix.
      irz = cp3_list1(3,icp,ilevel) - cp3_list0(3,icp,ilevel)  ! z-dir relative cell address of transformation matrix.

    !*** multipole to local translation
      icx0_1 = cp3_list0(1,icp  ,ilevel) ; icy0_1 = cp3_list0(2,icp  ,ilevel) ; icz0_1 = cp3_list0(3,icp  ,ilevel)
      icx1_1 = cp3_list1(1,icp  ,ilevel) ; icy1_1 = cp3_list1(2,icp  ,ilevel) ; icz1_1 = cp3_list1(3,icp  ,ilevel)
      icx0_2 = cp3_list0(1,icp+1,ilevel) ; icy0_2 = cp3_list0(2,icp+1,ilevel) ; icz0_2 = cp3_list0(3,icp+1,ilevel)
      icx1_2 = cp3_list1(1,icp+1,ilevel) ; icy1_2 = cp3_list1(2,icp+1,ilevel) ; icz1_2 = cp3_list1(3,icp+1,ilevel)
      icx0_3 = cp3_list0(1,icp+2,ilevel) ; icy0_3 = cp3_list0(2,icp+2,ilevel) ; icz0_3 = cp3_list0(3,icp+2,ilevel)
      icx1_3 = cp3_list1(1,icp+2,ilevel) ; icy1_3 = cp3_list1(2,icp+2,ilevel) ; icz1_3 = cp3_list1(3,icp+2,ilevel)

      include 'tuned_M2L/fmm_far_m2l_rotation_mpair3_nmax4_exsimd-reduct_unroll.h90'
    end do  ! cp3_count.
!$omp end do
    endif    ! cp3_count>0.

    if(cp2_count(ilevel) > 0) then    !! /*------- pairs multiple of 2 for each relative cell address -------*/

!$omp do
    do icp = 1, cp2_count(ilevel), 2
      irx = cp2_list1(1,icp,ilevel) - cp2_list0(1,icp,ilevel)  ! x-dir relative cell address of transformation matrix.
      iry = cp2_list1(2,icp,ilevel) - cp2_list0(2,icp,ilevel)  ! y-dir relative cell address of transformation matrix.
      irz = cp2_list1(3,icp,ilevel) - cp2_list0(3,icp,ilevel)  ! z-dir relative cell address of transformation matrix.

    !*** multipole to local translation
      icx0_1 = cp2_list0(1,icp  ,ilevel) ; icy0_1 = cp2_list0(2,icp  ,ilevel) ; icz0_1 = cp2_list0(3,icp  ,ilevel)
      icx1_1 = cp2_list1(1,icp  ,ilevel) ; icy1_1 = cp2_list1(2,icp  ,ilevel) ; icz1_1 = cp2_list1(3,icp  ,ilevel)
      icx0_2 = cp2_list0(1,icp+1,ilevel) ; icy0_2 = cp2_list0(2,icp+1,ilevel) ; icz0_2 = cp2_list0(3,icp+1,ilevel)
      icx1_2 = cp2_list1(1,icp+1,ilevel) ; icy1_2 = cp2_list1(2,icp+1,ilevel) ; icz1_2 = cp2_list1(3,icp+1,ilevel)

      include 'tuned_M2L/fmm_far_m2l_rotation_mpair2_nmax4_exsimd-reduct_unroll.h90'
    end do  ! cp2_count.
!$omp end do
    endif    ! cp2_count>0.

    if(cp1_count(ilevel) > 0) then    !! /*------- 1 pair for each relative cell address -------*/

!$omp do
    do icp = 1, cp1_count(ilevel)
      irx = cp1_list1(1,icp,ilevel) - cp1_list0(1,icp,ilevel)  ! x-dir relative cell address of transformation matrix.
      iry = cp1_list1(2,icp,ilevel) - cp1_list0(2,icp,ilevel)  ! y-dir relative cell address of transformation matrix.
      irz = cp1_list1(3,icp,ilevel) - cp1_list0(3,icp,ilevel)  ! z-dir relative cell address of transformation matrix.

      !*** multipole to local translation
      icx0 = cp1_list0(1,icp,ilevel) ; icy0 = cp1_list0(2,icp,ilevel) ; icz0 = cp1_list0(3,icp,ilevel)
      icx1 = cp1_list1(1,icp,ilevel) ; icy1 = cp1_list1(2,icp,ilevel) ; icz1 = cp1_list1(3,icp,ilevel)

      include 'tuned_M2L/fmm_far_m2l_rotation_nmax4_exsimd-reduct_unroll.h90'
    end do  ! cp1_count.
!$omp end do
    endif    ! cp1_count>0.

!$omp end parallel


    lcbx = nscxdiv + mbd_x2                          ! dimension size of superposition array.
    lcby = nscydiv + mbd_y2
    lcbz = nsczdiv + mbd_z2
    ncall = lcbx*lcby*lcbz                           ! overall candidates of superposition cell.
    ncblk = (ncall - 1)/nomp + 1                     ! block size for thread parallel. 

!$omp parallel default(none) &
!$omp& private(nt,icx0,icy0,icz0,icxyz,m1,iam) &
!$omp& shared(lm_length,wl,wl_omp) &
!$omp& shared(nscxdiv,nscydiv,nsczdiv) &
!$omp& shared(mbd_x ,mbd_y ,mbd_z) &
!$omp& shared(mbd_x2,mbd_y2,mbd_z2) &
!$omp& shared(lcby,lcbz) &
!$omp& shared(ncblk,ncall,nomp)
!$  iam = omp_get_thread_num()
     do nt = 0, nomp-1
       do icxyz = ncblk*iam + 1, min(ncblk*(iam+1), ncall)
        icz0 = mod(icxyz - 1, lcbz) + 1 - mbd_z      ! 1-dimensional to 3-dimensional cell index.
        icy0 = mod(icxyz - 1, lcbz*lcby)
        icy0 = icy0/lcbz + 1 - mbd_y                 ! 1-dimensional to 3-dimensional cell index.
        icx0 = (icxyz - 1)/(lcbz*lcby) + 1 - mbd_x   ! 1-dimensional to 3-dimensional cell index.
        ! exclude vacant target cell.
        if(nscxdiv == 1      .and. 1 - mbd_x /= icx0 .and. icx0 /= nscxdiv + mbd_x) cycle
        if(nscxdiv <  mbd_x2 .and. 1         <= icx0 .and. icx0 <= nscxdiv        ) cycle
        if(nscydiv == 1      .and. 1 - mbd_y /= icy0 .and. icy0 /= nscydiv + mbd_y) cycle
        if(nscydiv <  mbd_y2 .and. 1         <= icy0 .and. icy0 <= nscydiv        ) cycle
        if(nsczdiv == 1      .and. 1 - mbd_z /= icz0 .and. icz0 /= nsczdiv + mbd_z) cycle
        if(nsczdiv <  mbd_z2 .and. 1         <= icz0 .and. icz0 <= nsczdiv        ) cycle

        do m1=1,lm_length
          wl(m1,1,icz0,icy0,icx0) = &
     &    wl(m1,1,icz0,icy0,icx0) + wl_omp(m1,1,icz0,icy0,icx0,nt)
          wl(m1,2,icz0,icy0,icx0) = &
     &    wl(m1,2,icz0,icy0,icx0) + wl_omp(m1,2,icz0,icy0,icx0,nt)
        enddo

      end do  ! icxyz.
    enddo ! nt.
!$omp end parallel

  end subroutine M2L_lower_level_dr_mpairs_nmax4
!---------------------------------------------------------------------
!>
!! \brief  Subroutine to execute M2L (with local arrays, lower level)
!! \author Noriyuki Yoshii, Yoshimichi Andoh, Shin-chi Ichikawa
!<
  subroutine M2L_lower_level_dr_mpairs(ilevel,wm,wl, &
       &                        nsczdiv,nscydiv,nscxdiv, &
       &                        mbd_z,  mbd_y,  mbd_x, &
       &                        nomp)
!---------------------------------------------------------------------
! symmetry and rotation code.
! pipelined multiple cell pairs and general value of nmax.
! 
    use mpi_3d_grid, only : npx, npy, npz, ipx, ipy, ipz
    use domain, only : nlevel
    implicit none
    integer(kind=wip), intent(in) :: ilevel
    real(kind=wvp), intent(in) :: wm(lm_length,2, &
         &              1-nbd_zm(ilevel):nsczdiv+nbd_zp(ilevel), &
         &              1-nbd_ym(ilevel):nscydiv+nbd_yp(ilevel), &
         &              1-nbd_xm(ilevel):nscxdiv+nbd_xp(ilevel))
    real(kind=wvap), intent(inout) :: wl(lm_length,2, &
         &              1-mbd_z         :nsczdiv+mbd_z, &
         &              1-mbd_y         :nscydiv+mbd_y, &
         &              1-mbd_x         :nscxdiv+mbd_x      )
    integer(kind=wip), intent(in) :: nsczdiv,nscydiv,nscxdiv
    integer(kind=wip), intent(in) :: mbd_z, mbd_y, mbd_x
! add padding(+5,,+2,,,).
    real(kind=wvap) :: wl_omp(lm_length+5, 2, &
         &              1-mbd_z         :nsczdiv+2+mbd_z, &
         &              1-mbd_y         :nscydiv  +mbd_y, &
         &              1-mbd_x         :nscxdiv  +mbd_x, 0:nomp-1)
    integer(kind=wip) :: icx0, icy0, icz0
    integer(kind=wip) :: icx1, icy1, icz1
    integer(kind=wip) :: nomp
    integer(kind=wip) :: nt, iam
    integer(kind=wip) :: icxyz
    integer(kind=wip) :: mbd_z2, mbd_y2, mbd_x2
    integer(kind=wip) :: lcbx, lcby, lcbz
    integer(kind=wip) :: irx, iry, irz
    integer(kind=wip) :: ncall, ncblk
    real(kind=wvp) :: wwm_north(lm_length,2,8), wwl_north(lm_length,2,8)
!add for simultaneous multiple cell pairs.
    integer(kind=wip) :: icp
    integer(kind=wip) :: icx0_1,icy0_1,icz0_1,icx1_1,icy1_1,icz1_1
    integer(kind=wip) :: icx0_2,icy0_2,icz0_2,icx1_2,icy1_2,icz1_2
    integer(kind=wip) :: icx0_3,icy0_3,icz0_3,icx1_3,icy1_3,icz1_3
    integer(kind=wip) :: icx0_4,icy0_4,icz0_4,icx1_4,icy1_4,icz1_4
    integer(kind=wip) :: icx0_5,icy0_5,icz0_5,icx1_5,icy1_5,icz1_5
    integer(kind=wip) :: icx0_6,icy0_6,icz0_6,icx1_6,icy1_6,icz1_6
    integer(kind=wip) :: icx0_7,icy0_7,icz0_7,icx1_7,icy1_7,icz1_7
    integer(kind=wip) :: icx0_8,icy0_8,icz0_8,icx1_8,icy1_8,icz1_8
    real(kind=wvp) :: shrot1_r, shrot2_r, shrot1_i, shrot2_i
    real(kind=wvp) :: shinvrot1_r, shinvrot2_r, shinvrot1_i, shinvrot2_i
    real(kind=wvp) :: wm_1r, wm_2r, wm_3r, wm_4r, wm_5r, wm_6r, wm_7r, wm_8r
    real(kind=wvp) :: wm_1i, wm_2i, wm_3i, wm_4i, wm_5i, wm_6i, wm_7i, wm_8i
    real(kind=wvp) :: wlnorth_1r,wlnorth_2r,wlnorth_3r,wlnorth_4r,wlnorth_5r,wlnorth_6r,wlnorth_7r,wlnorth_8r
    real(kind=wvp) :: wlnorth_1i,wlnorth_2i,wlnorth_3i,wlnorth_4i,wlnorth_5i,wlnorth_6i,wlnorth_7i,wlnorth_8i
    real(kind=wvp) :: wm_north(lm_length,2), wl_north(lm_length,2)
    real(kind=wvp) :: tmp_shml, tmp_ce_r, tmp_ce_i
    integer(kind=wip) :: j,k,n,p,q, m1,m2,m2_b, ind_1dim,ind_1dim_b

    iam = 0

    wl_omp = Cdzero

    ! global address offset of a cell prior to starting cell.
    mbd_z2 = 2*mbd_z
    mbd_y2 = 2*mbd_y
    mbd_x2 = 2*mbd_x

!$omp parallel default(none) &
!$omp& private(iam) &
!$omp& private(icx0,icy0,icz0,icx1,icy1,icz1) &
!$omp& private(irx,iry,irz) &
!$omp& private(wwm_north,wwl_north) &
!$omp& private(icp) &
!$omp& private(icx0_1,icy0_1,icz0_1,icx1_1,icy1_1,icz1_1) &
!$omp& private(icx0_2,icy0_2,icz0_2,icx1_2,icy1_2,icz1_2) &
!$omp& private(icx0_3,icy0_3,icz0_3,icx1_3,icy1_3,icz1_3) &
!$omp& private(icx0_4,icy0_4,icz0_4,icx1_4,icy1_4,icz1_4) &
!$omp& private(icx0_5,icy0_5,icz0_5,icx1_5,icy1_5,icz1_5) &
!$omp& private(icx0_6,icy0_6,icz0_6,icx1_6,icy1_6,icz1_6) &
!$omp& private(icx0_7,icy0_7,icz0_7,icx1_7,icy1_7,icz1_7) &
!$omp& private(icx0_8,icy0_8,icz0_8,icx1_8,icy1_8,icz1_8) &
!$omp& private(shrot1_r, shrot2_r, shrot1_i, shrot2_i) &
!$omp& private(shinvrot1_r, shinvrot2_r, shinvrot1_i, shinvrot2_i) &
!$omp& private(wm_1r, wm_2r, wm_3r, wm_4r, wm_5r, wm_6r, wm_7r, wm_8r) &
!$omp& private(wm_1i, wm_2i, wm_3i, wm_4i, wm_5i, wm_6i, wm_7i, wm_8i) &
!$omp& private(wlnorth_1r,wlnorth_2r,wlnorth_3r,wlnorth_4r,wlnorth_5r,wlnorth_6r,wlnorth_7r,wlnorth_8r) &
!$omp& private(wlnorth_1i,wlnorth_2i,wlnorth_3i,wlnorth_4i,wlnorth_5i,wlnorth_6i,wlnorth_7i,wlnorth_8i) &
!$omp& private(wm_north,wl_north) &
!$omp& private(j,k,n,p,q, m1,m2,m2_b, ind_1dim,ind_1dim_b) &
!$omp& private(tmp_shml, tmp_ce_r, tmp_ce_i) &
!$omp& shared(shml, sh_rot1, sh_rot2, sh_inv_rot1, sh_inv_rot2) &
!$omp& shared(nmax) &
!$omp& shared(ilevel,lm_length) &
!$omp& shared(wl_omp,wm) &
!$omp& shared(cp8_count,cp6_count,cp4_count,cp3_count,cp2_count,cp1_count) &
!$omp& shared(cp8_list0,cp6_list0,cp4_list0,cp3_list0,cp2_list0,cp1_list0) &
!$omp& shared(cp8_list1,cp6_list1,cp4_list1,cp3_list1,cp2_list1,cp1_list1)
!$  iam = omp_get_thread_num()

    if(cp8_count(ilevel) > 0) then    !! /*------- pairs multiple of 8 for each relative cell address -------*/

!$omp do
    do icp = 1, cp8_count(ilevel), 8
      irx = cp8_list1(1,icp,ilevel) - cp8_list0(1,icp,ilevel)  ! x-dir relative cell address of transformation matrix.
      iry = cp8_list1(2,icp,ilevel) - cp8_list0(2,icp,ilevel)  ! y-dir relative cell address of transformation matrix.
      irz = cp8_list1(3,icp,ilevel) - cp8_list0(3,icp,ilevel)  ! z-dir relative cell address of transformation matrix.

    !*** multipole to local translation
      icx0_1 = cp8_list0(1,icp  ,ilevel) ; icy0_1 = cp8_list0(2,icp  ,ilevel) ; icz0_1 = cp8_list0(3,icp  ,ilevel)
      icx1_1 = cp8_list1(1,icp  ,ilevel) ; icy1_1 = cp8_list1(2,icp  ,ilevel) ; icz1_1 = cp8_list1(3,icp  ,ilevel)
      icx0_2 = cp8_list0(1,icp+1,ilevel) ; icy0_2 = cp8_list0(2,icp+1,ilevel) ; icz0_2 = cp8_list0(3,icp+1,ilevel)
      icx1_2 = cp8_list1(1,icp+1,ilevel) ; icy1_2 = cp8_list1(2,icp+1,ilevel) ; icz1_2 = cp8_list1(3,icp+1,ilevel)
      icx0_3 = cp8_list0(1,icp+2,ilevel) ; icy0_3 = cp8_list0(2,icp+2,ilevel) ; icz0_3 = cp8_list0(3,icp+2,ilevel)
      icx1_3 = cp8_list1(1,icp+2,ilevel) ; icy1_3 = cp8_list1(2,icp+2,ilevel) ; icz1_3 = cp8_list1(3,icp+2,ilevel)
      icx0_4 = cp8_list0(1,icp+3,ilevel) ; icy0_4 = cp8_list0(2,icp+3,ilevel) ; icz0_4 = cp8_list0(3,icp+3,ilevel)
      icx1_4 = cp8_list1(1,icp+3,ilevel) ; icy1_4 = cp8_list1(2,icp+3,ilevel) ; icz1_4 = cp8_list1(3,icp+3,ilevel)
      icx0_5 = cp8_list0(1,icp+4,ilevel) ; icy0_5 = cp8_list0(2,icp+4,ilevel) ; icz0_5 = cp8_list0(3,icp+4,ilevel)
      icx1_5 = cp8_list1(1,icp+4,ilevel) ; icy1_5 = cp8_list1(2,icp+4,ilevel) ; icz1_5 = cp8_list1(3,icp+4,ilevel)
      icx0_6 = cp8_list0(1,icp+5,ilevel) ; icy0_6 = cp8_list0(2,icp+5,ilevel) ; icz0_6 = cp8_list0(3,icp+5,ilevel)
      icx1_6 = cp8_list1(1,icp+5,ilevel) ; icy1_6 = cp8_list1(2,icp+5,ilevel) ; icz1_6 = cp8_list1(3,icp+5,ilevel)
      icx0_7 = cp8_list0(1,icp+6,ilevel) ; icy0_7 = cp8_list0(2,icp+6,ilevel) ; icz0_7 = cp8_list0(3,icp+6,ilevel)
      icx1_7 = cp8_list1(1,icp+6,ilevel) ; icy1_7 = cp8_list1(2,icp+6,ilevel) ; icz1_7 = cp8_list1(3,icp+6,ilevel)
      icx0_8 = cp8_list0(1,icp+7,ilevel) ; icy0_8 = cp8_list0(2,icp+7,ilevel) ; icz0_8 = cp8_list0(3,icp+7,ilevel)
      icx1_8 = cp8_list1(1,icp+7,ilevel) ; icy1_8 = cp8_list1(2,icp+7,ilevel) ; icz1_8 = cp8_list1(3,icp+7,ilevel)

      include 'tuned_M2L/fmm_far_m2l_rotation_mpair8_exsimd-reduct.h90'
    end do  ! cp8_count.
!$omp end do

    endif    ! cp8_count>0.

    if(cp6_count(ilevel) > 0) then    !! /*------- pairs multiple of 6 for each relative cell address -------*/

!$omp do
    do icp = 1, cp6_count(ilevel), 6
      irx = cp6_list1(1,icp,ilevel) - cp6_list0(1,icp,ilevel)  ! x-dir relative cell address of transformation matrix.
      iry = cp6_list1(2,icp,ilevel) - cp6_list0(2,icp,ilevel)  ! y-dir relative cell address of transformation matrix.
      irz = cp6_list1(3,icp,ilevel) - cp6_list0(3,icp,ilevel)  ! z-dir relative cell address of transformation matrix.

    !*** multipole to local translation
      icx0_1 = cp6_list0(1,icp  ,ilevel) ; icy0_1 = cp6_list0(2,icp  ,ilevel) ; icz0_1 = cp6_list0(3,icp  ,ilevel)
      icx1_1 = cp6_list1(1,icp  ,ilevel) ; icy1_1 = cp6_list1(2,icp  ,ilevel) ; icz1_1 = cp6_list1(3,icp  ,ilevel)
      icx0_2 = cp6_list0(1,icp+1,ilevel) ; icy0_2 = cp6_list0(2,icp+1,ilevel) ; icz0_2 = cp6_list0(3,icp+1,ilevel)
      icx1_2 = cp6_list1(1,icp+1,ilevel) ; icy1_2 = cp6_list1(2,icp+1,ilevel) ; icz1_2 = cp6_list1(3,icp+1,ilevel)
      icx0_3 = cp6_list0(1,icp+2,ilevel) ; icy0_3 = cp6_list0(2,icp+2,ilevel) ; icz0_3 = cp6_list0(3,icp+2,ilevel)
      icx1_3 = cp6_list1(1,icp+2,ilevel) ; icy1_3 = cp6_list1(2,icp+2,ilevel) ; icz1_3 = cp6_list1(3,icp+2,ilevel)
      icx0_4 = cp6_list0(1,icp+3,ilevel) ; icy0_4 = cp6_list0(2,icp+3,ilevel) ; icz0_4 = cp6_list0(3,icp+3,ilevel)
      icx1_4 = cp6_list1(1,icp+3,ilevel) ; icy1_4 = cp6_list1(2,icp+3,ilevel) ; icz1_4 = cp6_list1(3,icp+3,ilevel)
      icx0_5 = cp6_list0(1,icp+4,ilevel) ; icy0_5 = cp6_list0(2,icp+4,ilevel) ; icz0_5 = cp6_list0(3,icp+4,ilevel)
      icx1_5 = cp6_list1(1,icp+4,ilevel) ; icy1_5 = cp6_list1(2,icp+4,ilevel) ; icz1_5 = cp6_list1(3,icp+4,ilevel)
      icx0_6 = cp6_list0(1,icp+5,ilevel) ; icy0_6 = cp6_list0(2,icp+5,ilevel) ; icz0_6 = cp6_list0(3,icp+5,ilevel)
      icx1_6 = cp6_list1(1,icp+5,ilevel) ; icy1_6 = cp6_list1(2,icp+5,ilevel) ; icz1_6 = cp6_list1(3,icp+5,ilevel)

      include 'tuned_M2L/fmm_far_m2l_rotation_mpair6_exsimd-reduct.h90'
    end do  ! cp6_count.
!$omp end do
    endif    ! cp6_count>0.

    if(cp4_count(ilevel) > 0) then    !! /*------- pairs multiple of 4 for each relative cell address -------*/

!$omp do
    do icp = 1, cp4_count(ilevel), 4
      irx = cp4_list1(1,icp,ilevel) - cp4_list0(1,icp,ilevel)  ! x-dir relative cell address of transformation matrix.
      iry = cp4_list1(2,icp,ilevel) - cp4_list0(2,icp,ilevel)  ! y-dir relative cell address of transformation matrix.
      irz = cp4_list1(3,icp,ilevel) - cp4_list0(3,icp,ilevel)  ! z-dir relative cell address of transformation matrix.

    !*** multipole to local translation
      icx0_1 = cp4_list0(1,icp  ,ilevel) ; icy0_1 = cp4_list0(2,icp  ,ilevel) ; icz0_1 = cp4_list0(3,icp  ,ilevel)
      icx1_1 = cp4_list1(1,icp  ,ilevel) ; icy1_1 = cp4_list1(2,icp  ,ilevel) ; icz1_1 = cp4_list1(3,icp  ,ilevel)
      icx0_2 = cp4_list0(1,icp+1,ilevel) ; icy0_2 = cp4_list0(2,icp+1,ilevel) ; icz0_2 = cp4_list0(3,icp+1,ilevel)
      icx1_2 = cp4_list1(1,icp+1,ilevel) ; icy1_2 = cp4_list1(2,icp+1,ilevel) ; icz1_2 = cp4_list1(3,icp+1,ilevel)
      icx0_3 = cp4_list0(1,icp+2,ilevel) ; icy0_3 = cp4_list0(2,icp+2,ilevel) ; icz0_3 = cp4_list0(3,icp+2,ilevel)
      icx1_3 = cp4_list1(1,icp+2,ilevel) ; icy1_3 = cp4_list1(2,icp+2,ilevel) ; icz1_3 = cp4_list1(3,icp+2,ilevel)
      icx0_4 = cp4_list0(1,icp+3,ilevel) ; icy0_4 = cp4_list0(2,icp+3,ilevel) ; icz0_4 = cp4_list0(3,icp+3,ilevel)
      icx1_4 = cp4_list1(1,icp+3,ilevel) ; icy1_4 = cp4_list1(2,icp+3,ilevel) ; icz1_4 = cp4_list1(3,icp+3,ilevel)

      include 'tuned_M2L/fmm_far_m2l_rotation_mpair4_exsimd-reduct.h90'
    end do  ! cp4_count.
!$omp end do
    endif    ! cp4_count>0.

    if(cp3_count(ilevel) > 0) then    !! /*------- pairs multiple of 3 for each relative cell address -------*/

!$omp do
    do icp = 1, cp3_count(ilevel), 3
      irx = cp3_list1(1,icp,ilevel) - cp3_list0(1,icp,ilevel)  ! x-dir relative cell address of transformation matrix.
      iry = cp3_list1(2,icp,ilevel) - cp3_list0(2,icp,ilevel)  ! y-dir relative cell address of transformation matrix.
      irz = cp3_list1(3,icp,ilevel) - cp3_list0(3,icp,ilevel)  ! z-dir relative cell address of transformation matrix.

    !*** multipole to local translation
      icx0_1 = cp3_list0(1,icp  ,ilevel) ; icy0_1 = cp3_list0(2,icp  ,ilevel) ; icz0_1 = cp3_list0(3,icp  ,ilevel)
      icx1_1 = cp3_list1(1,icp  ,ilevel) ; icy1_1 = cp3_list1(2,icp  ,ilevel) ; icz1_1 = cp3_list1(3,icp  ,ilevel)
      icx0_2 = cp3_list0(1,icp+1,ilevel) ; icy0_2 = cp3_list0(2,icp+1,ilevel) ; icz0_2 = cp3_list0(3,icp+1,ilevel)
      icx1_2 = cp3_list1(1,icp+1,ilevel) ; icy1_2 = cp3_list1(2,icp+1,ilevel) ; icz1_2 = cp3_list1(3,icp+1,ilevel)
      icx0_3 = cp3_list0(1,icp+2,ilevel) ; icy0_3 = cp3_list0(2,icp+2,ilevel) ; icz0_3 = cp3_list0(3,icp+2,ilevel)
      icx1_3 = cp3_list1(1,icp+2,ilevel) ; icy1_3 = cp3_list1(2,icp+2,ilevel) ; icz1_3 = cp3_list1(3,icp+2,ilevel)

      include 'tuned_M2L/fmm_far_m2l_rotation_mpair3_exsimd-reduct.h90'
    end do  ! cp3_count.
!$omp end do
    endif    ! cp3_count>0.

    if(cp2_count(ilevel) > 0) then    !! /*------- pairs multiple of 2 for each relative cell address -------*/

!$omp do
    do icp = 1, cp2_count(ilevel), 2
      irx = cp2_list1(1,icp,ilevel) - cp2_list0(1,icp,ilevel)  ! x-dir relative cell address of transformation matrix.
      iry = cp2_list1(2,icp,ilevel) - cp2_list0(2,icp,ilevel)  ! y-dir relative cell address of transformation matrix.
      irz = cp2_list1(3,icp,ilevel) - cp2_list0(3,icp,ilevel)  ! z-dir relative cell address of transformation matrix.

    !*** multipole to local translation
      icx0_1 = cp2_list0(1,icp  ,ilevel) ; icy0_1 = cp2_list0(2,icp  ,ilevel) ; icz0_1 = cp2_list0(3,icp  ,ilevel)
      icx1_1 = cp2_list1(1,icp  ,ilevel) ; icy1_1 = cp2_list1(2,icp  ,ilevel) ; icz1_1 = cp2_list1(3,icp  ,ilevel)
      icx0_2 = cp2_list0(1,icp+1,ilevel) ; icy0_2 = cp2_list0(2,icp+1,ilevel) ; icz0_2 = cp2_list0(3,icp+1,ilevel)
      icx1_2 = cp2_list1(1,icp+1,ilevel) ; icy1_2 = cp2_list1(2,icp+1,ilevel) ; icz1_2 = cp2_list1(3,icp+1,ilevel)

      include 'tuned_M2L/fmm_far_m2l_rotation_mpair2_exsimd-reduct.h90'
    end do  ! cp2_count.
!$omp end do
    endif    ! cp2_count>0.

    if(cp1_count(ilevel) > 0) then    !! /*------- 1 pair for each relative cell address -------*/

!$omp do
    do icp = 1, cp1_count(ilevel)
      irx = cp1_list1(1,icp,ilevel) - cp1_list0(1,icp,ilevel)  ! x-dir relative cell address of transformation matrix.
      iry = cp1_list1(2,icp,ilevel) - cp1_list0(2,icp,ilevel)  ! y-dir relative cell address of transformation matrix.
      irz = cp1_list1(3,icp,ilevel) - cp1_list0(3,icp,ilevel)  ! z-dir relative cell address of transformation matrix.

      !*** multipole to local translation
      icx0 = cp1_list0(1,icp,ilevel) ; icy0 = cp1_list0(2,icp,ilevel) ; icz0 = cp1_list0(3,icp,ilevel)
      icx1 = cp1_list1(1,icp,ilevel) ; icy1 = cp1_list1(2,icp,ilevel) ; icz1 = cp1_list1(3,icp,ilevel)

      include 'tuned_M2L/fmm_far_m2l_rotation_exsimd-reduct.h90'
    end do  ! cp1_count.
!$omp end do
    endif    ! cp1_count>0.

!$omp end parallel


    lcbx = nscxdiv + mbd_x2                          ! dimension size of superposition array.
    lcby = nscydiv + mbd_y2
    lcbz = nsczdiv + mbd_z2
    ncall = lcbx*lcby*lcbz                           ! overall candidates of superposition cell.
    ncblk = (ncall - 1)/nomp + 1                     ! block size for thread parallel. 

!$omp parallel default(none) &
!$omp& private(nt,icx0,icy0,icz0,icxyz,m1,iam) &
!$omp& shared(lm_length,wl,wl_omp) &
!$omp& shared(nscxdiv,nscydiv,nsczdiv) &
!$omp& shared(mbd_x ,mbd_y ,mbd_z) &
!$omp& shared(mbd_x2,mbd_y2,mbd_z2) &
!$omp& shared(lcby,lcbz) &
!$omp& shared(ncblk,ncall,nomp)
!$  iam = omp_get_thread_num()
     do nt = 0, nomp-1
       do icxyz = ncblk*iam + 1, min(ncblk*(iam+1), ncall)
        icz0 = mod(icxyz - 1, lcbz) + 1 - mbd_z      ! 1-dimensional to 3-dimensional cell index.
        icy0 = mod(icxyz - 1, lcbz*lcby)
        icy0 = icy0/lcbz + 1 - mbd_y                 ! 1-dimensional to 3-dimensional cell index.
        icx0 = (icxyz - 1)/(lcbz*lcby) + 1 - mbd_x   ! 1-dimensional to 3-dimensional cell index.
        ! exclude vacant target cell.
        if(nscxdiv == 1      .and. 1 - mbd_x /= icx0 .and. icx0 /= nscxdiv + mbd_x) cycle
        if(nscxdiv <  mbd_x2 .and. 1         <= icx0 .and. icx0 <= nscxdiv        ) cycle
        if(nscydiv == 1      .and. 1 - mbd_y /= icy0 .and. icy0 /= nscydiv + mbd_y) cycle
        if(nscydiv <  mbd_y2 .and. 1         <= icy0 .and. icy0 <= nscydiv        ) cycle
        if(nsczdiv == 1      .and. 1 - mbd_z /= icz0 .and. icz0 /= nsczdiv + mbd_z) cycle
        if(nsczdiv <  mbd_z2 .and. 1         <= icz0 .and. icz0 <= nsczdiv        ) cycle

        do m1=1,lm_length
          wl(m1,1,icz0,icy0,icx0) = &
     &    wl(m1,1,icz0,icy0,icx0) + wl_omp(m1,1,icz0,icy0,icx0,nt)
          wl(m1,2,icz0,icy0,icx0) = &
     &    wl(m1,2,icz0,icy0,icx0) + wl_omp(m1,2,icz0,icy0,icx0,nt)
        enddo

      end do  ! icxyz.
    enddo ! nt.
!$omp end parallel

  end subroutine M2L_lower_level_dr_mpairs
!---------------------------------------------------------------------
#else   /* Not MOMENT */
!---------------------------------------------------------------------
!>
!! \brief  Subroutine to execute M2L (with local arrays, lower level)
!! \author Noriyuki Yoshii, Yoshimichi Andoh, Shin-chi Ichikawa
!<
  subroutine M2L_lower_level_dr_nmax4(ilevel,wm,wl, &
       &                        nsczdiv,nscydiv,nscxdiv, &
       &                        nscellz,nscelly,nscellx, &
       &                        mbd_z,  mbd_y,  mbd_x, &
       &                        nomp)
!---------------------------------------------------------------------
! single pair symmetry and rotation code. fully unrolled for nmax=4.
!
    use mpi_3d_grid, only : npx, npy, npz, ipx, ipy, ipz
    use domain, only : nlevel
    implicit none
    integer(kind=wip), intent(in) :: ilevel
    real(kind=wvp), intent(in) :: wm(lm_length,2, &
         &              1-nbd_zm(ilevel):nsczdiv+nbd_zp(ilevel), &
         &              1-nbd_ym(ilevel):nscydiv+nbd_yp(ilevel), &
         &              1-nbd_xm(ilevel):nscxdiv+nbd_xp(ilevel))
    real(kind=wvap), intent(inout) :: wl(lm_length,2, &
         &              1-mbd_z         :nsczdiv+mbd_z, &
         &              1-mbd_y         :nscydiv+mbd_y, &
         &              1-mbd_x         :nscxdiv+mbd_x      )
    integer(kind=wip), intent(in) :: nsczdiv,nscydiv,nscxdiv
    integer(kind=wip), intent(in) :: nscellx, nscelly, nscellz
    integer(kind=wip), intent(in) :: mbd_z, mbd_y, mbd_x
! add padding(+5,,+2,,,).
    real(kind=wvap) :: wl_omp(lm_length+5, 2, &
         &              1-mbd_z         :nsczdiv+2+mbd_z, &
         &              1-mbd_y         :nscydiv  +mbd_y, &
         &              1-mbd_x         :nscxdiv  +mbd_x, 0:nomp-1)
    integer(kind=wip) :: icx0, icy0, icz0
    integer(kind=wip) :: icx1, icy1, icz1
    integer(kind=wip) :: icx1s, icy1s, icz1s
    integer(kind=wip) :: icx1e, icy1e, icz1e
    integer(kind=wip) :: nomp
    integer(kind=wip) :: nt, iam
    integer(kind=wip) :: icxgo, icygo, iczgo
    integer(kind=wip) :: icxg,  icyg,  iczg
    integer(kind=wip) :: icxyz
    integer(kind=wip) :: mbd_z2, mbd_y2, mbd_x2
    integer(kind=wip) :: ngord, nrgord
    integer(kind=wip) :: lcbx, lcby, lcbz
    integer(kind=wip) :: irx, iry, irz, load
    integer(kind=wip) :: ncall, ncblk
    real(kind=wvp) :: wm_north(lm_length,2), wl_north(lm_length,2)
    real(kind=wvp) :: tmp_ce_r, tmp_ce_i
    integer(kind=wip) :: m1

    iam = 0

    wl_omp = Cdzero

    ! global address offset of a cell prior to starting cell.
    icxgo = (nscellx * ipx) / npx
    icygo = (nscelly * ipy) / npy
    iczgo = (nscellz * ipz) / npz
    mbd_z2 = 2*mbd_z
    mbd_y2 = 2*mbd_y
    mbd_x2 = 2*mbd_x

!$omp parallel default(none) &
!$omp& private(iam) &
!$omp& private(load) &
!$omp& private(icx0,icy0,icz0,icx1,icy1,icz1) &
!$omp& private(icx1s,icy1s,icz1s) &
!$omp& private(icx1e,icy1e,icz1e) &
!$omp& private(icxg,icyg,iczg) &
!$omp& private(irx,iry,irz) &
!$omp& private(ngord,nrgord) &
!$omp& private(wm_north,wl_north) &
!$omp& private(m1) &
!$omp& private(tmp_ce_r, tmp_ce_i) &
!$omp& shared(shml, sh_rot1, sh_rot2, sh_inv_rot1, sh_inv_rot2) &
!$omp& shared(ilevel,lm_length) &
!$omp& shared(wl_omp,wm) &
!$omp& shared(nload,lddir) &
!$omp& shared(nscellx,nscelly,nscellz) &
!$omp& shared(nscxdiv,nscydiv,nsczdiv) &
!$omp& shared(mbd_x ,mbd_y ,mbd_z) &
!$omp& shared(mbd_x2,mbd_y2,mbd_z2) &
!$omp& shared(icxgo,icygo,iczgo)
!$  iam = omp_get_thread_num()
!$omp do schedule(static,1)
    do load = 1, nload(ilevel)
      irx = lddir(1,load,ilevel)     ! x-dir relative cell address of transformation matrix.
      iry = lddir(2,load,ilevel)     ! y-dir relative cell address of transformation matrix.
      irz = lddir(3,load,ilevel)     ! z-dir relative cell address of transformation matrix.

      if(abs(irx).le.2 .and. abs(iry).le.2 .and. abs(irz).le.2) cycle  ! lddir does not include true case.

    do icx0 = 1-mbd_x, nscxdiv+mbd_x      ! target cell x-address.
      if(nscxdiv == 1      .and. 1 - mbd_x /= icx0 .and. icx0 /= nscxdiv + mbd_x) cycle
      if(nscxdiv <  mbd_x2 .and. 1         <= icx0 .and. icx0 <= nscxdiv        ) cycle
      icxg  = modulo(icxgo + icx0 - 1, nscellx) + 1  ! target cell global x-address.
      ngord  = mod(icxg - 1, mbd_x)                  ! sequence number within cell merge group minus 1.
      nrgord = mbd_x - ngord - 1                     ! reverse sequence number within cell merge group minus 1.
      icx1s = icx0 - mbd_x2 - ngord                  ! first cell of the cell range to be transformed.
      icx1e = icx0 + mbd_x2 + nrgord                 ! final cell of the cell range to be transformed.
      if(icx0 <= mbd_x          ) icx1s = icx0       ! plus-direction half-space transformation.
      if(icx0 >  nscxdiv - mbd_x) icx1e = icx0 - 1   ! minus-direction half-space transformation. inter-process exclusive sontrol.
      icx1 = icx0 + irx                              ! real transforming cell defined by relative cell address.
      if(icx1 < icx1s .or. icx1e < icx1) cycle

      do icy0 = 1-mbd_y, nscydiv+mbd_y    ! target cell y-address.
        if(nscydiv == 1     .and. 1 - mbd_y /= icy0 .and. icy0 /= nscydiv + mbd_y) cycle
        if(nscydiv < mbd_y2 .and. 1         <= icy0 .and. icy0 <= nscydiv        ) cycle
        icyg  = modulo(icygo + icy0 - 1, nscelly) + 1
        ngord  = mod(icyg - 1, mbd_y)
        nrgord = mbd_y - ngord - 1
        icy1s = icy0 - mbd_y2 - ngord
        icy1e = icy0 + mbd_y2 + nrgord
        if(icy0 <= mbd_y          ) icy1s = icy0
        if(icy0 >  nscydiv - mbd_y) icy1e = icy0 - 1
        icy1 = icy0 + iry
        if(icy1 < icy1s .or. icy1e < icy1) cycle

        do icz0 = 1-mbd_z, nsczdiv+mbd_z  ! target cell z-address.
          if(nsczdiv == 1     .and. 1 - mbd_z /= icz0 .and. icz0 /= nsczdiv + mbd_z) cycle
          if(nsczdiv < mbd_z2 .and. 1         <= icz0 .and. icz0 <= nsczdiv        ) cycle
          iczg  = modulo(iczgo + icz0 - 1, nscellz) + 1
          ngord  = mod(iczg - 1, mbd_z)
          nrgord = mbd_z - ngord - 1
          icz1s = icz0 - mbd_z2 - ngord
          icz1e = icz0 + mbd_z2 + nrgord
          if(icz0 <= mbd_z          ) icz1s = icz0
          if(icz0 >  nsczdiv - mbd_z) icz1e = icz0 - 1
          icz1 = icz0 + irz
          if(icz1 < icz1s .or. icz1e < icz1) cycle

          include 'tuned_M2L/fmm_far_m2l_rotation_nmax4_unrolled_fj.h90'

        end do  ! icz0.
      end do  ! icy0.
    end do  ! icx0.

    end do  ! load.
!$omp end do
!$omp end parallel

    lcbx = nscxdiv + mbd_x2                          ! dimension size of superposition array.
    lcby = nscydiv + mbd_y2
    lcbz = nsczdiv + mbd_z2
    ncall = lcbx*lcby*lcbz                           ! overall candidates of superposition cell.
    ncblk = (ncall - 1)/nomp + 1                     ! block size for thread parallel. 

!$omp parallel default(none) &
!$omp& private(nt,icx0,icy0,icz0,icxyz,m1,iam) &
!$omp& shared(lm_length,wl,wl_omp) &
!$omp& shared(nscxdiv,nscydiv,nsczdiv) &
!$omp& shared(mbd_x ,mbd_y ,mbd_z) &
!$omp& shared(mbd_x2,mbd_y2,mbd_z2) &
!$omp& shared(lcby,lcbz) &
!$omp& shared(ncblk,ncall,nomp)
!$  iam = omp_get_thread_num()
     do nt = 0, nomp-1
       do icxyz = ncblk*iam + 1, min(ncblk*(iam+1), ncall)
        icz0 = mod(icxyz - 1, lcbz) + 1 - mbd_z      ! 1-dimensional to 3-dimensional cell index.
        icy0 = mod(icxyz - 1, lcbz*lcby)
        icy0 = icy0/lcbz + 1 - mbd_y                 ! 1-dimensional to 3-dimensional cell index.
        icx0 = (icxyz - 1)/(lcbz*lcby) + 1 - mbd_x   ! 1-dimensional to 3-dimensional cell index.
        ! exclude vacant target cell.
        if(nscxdiv == 1      .and. 1 - mbd_x /= icx0 .and. icx0 /= nscxdiv + mbd_x) cycle
        if(nscxdiv <  mbd_x2 .and. 1         <= icx0 .and. icx0 <= nscxdiv        ) cycle
        if(nscydiv == 1      .and. 1 - mbd_y /= icy0 .and. icy0 /= nscydiv + mbd_y) cycle
        if(nscydiv <  mbd_y2 .and. 1         <= icy0 .and. icy0 <= nscydiv        ) cycle
        if(nsczdiv == 1      .and. 1 - mbd_z /= icz0 .and. icz0 /= nsczdiv + mbd_z) cycle
        if(nsczdiv <  mbd_z2 .and. 1         <= icz0 .and. icz0 <= nsczdiv        ) cycle

        do m1=1,lm_length
          wl(m1,1,icz0,icy0,icx0) = &
     &    wl(m1,1,icz0,icy0,icx0) + wl_omp(m1,1,icz0,icy0,icx0,nt)
          wl(m1,2,icz0,icy0,icx0) = &
     &    wl(m1,2,icz0,icy0,icx0) + wl_omp(m1,2,icz0,icy0,icx0,nt)
        enddo

      end do  ! icxyz.
    enddo ! nt.
!$omp end parallel

  end subroutine M2L_lower_level_dr_nmax4
!---------------------------------------------------------------------
!>
!! \brief  Subroutine to execute M2L (with local arrays, lower level)
!! \author Noriyuki Yoshii, Yoshimichi Andoh, Shin-chi Ichikawa
!<
  subroutine M2L_lower_level_dr_general(ilevel,wm,wl, &
       &                        nsczdiv,nscydiv,nscxdiv, &
       &                        nscellz,nscelly,nscellx, &
       &                        mbd_z,  mbd_y,  mbd_x, &
       &                        nomp)
!---------------------------------------------------------------------
! single pair symmetry and rotation code. general value of nmax.
!
    use mpi_3d_grid, only : npx, npy, npz, ipx, ipy, ipz
    use domain, only : nlevel
    implicit none
    integer(kind=wip), intent(in) :: ilevel
    real(kind=wvp), intent(in) :: wm(lm_length,2, &
         &              1-nbd_zm(ilevel):nsczdiv+nbd_zp(ilevel), &
         &              1-nbd_ym(ilevel):nscydiv+nbd_yp(ilevel), &
         &              1-nbd_xm(ilevel):nscxdiv+nbd_xp(ilevel))
    real(kind=wvap), intent(inout) :: wl(lm_length,2, &
         &              1-mbd_z         :nsczdiv+mbd_z, &
         &              1-mbd_y         :nscydiv+mbd_y, &
         &              1-mbd_x         :nscxdiv+mbd_x      )
    integer(kind=wip), intent(in) :: nsczdiv,nscydiv,nscxdiv
    integer(kind=wip), intent(in) :: nscellx, nscelly, nscellz
    integer(kind=wip), intent(in) :: mbd_z, mbd_y, mbd_x
! add padding(+5,,+2,,,).
    real(kind=wvap) :: wl_omp(lm_length+5, 2, &
         &              1-mbd_z         :nsczdiv+2+mbd_z, &
         &              1-mbd_y         :nscydiv  +mbd_y, &
         &              1-mbd_x         :nscxdiv  +mbd_x, 0:nomp-1)
    integer(kind=wip) :: icx0, icy0, icz0
    integer(kind=wip) :: icx1, icy1, icz1
    integer(kind=wip) :: icx1s, icy1s, icz1s
    integer(kind=wip) :: icx1e, icy1e, icz1e
    integer(kind=wip) :: nomp
    integer(kind=wip) :: nt, iam
    integer(kind=wip) :: icxgo, icygo, iczgo
    integer(kind=wip) :: icxg,  icyg,  iczg
    integer(kind=wip) :: icxyz
    integer(kind=wip) :: mbd_z2, mbd_y2, mbd_x2
    integer(kind=wip) :: ngord, nrgord
    integer(kind=wip) :: lcbx, lcby, lcbz
    integer(kind=wip) :: irx, iry, irz, load
    integer(kind=wip) :: ncall, ncblk
    real(kind=wvp) :: wm_north(lm_length,2), wl_north(lm_length,2)
    real(kind=wvp) :: tmp_shml, tmp_ce_r, tmp_ce_i
    integer(kind=wip) :: j,k,n,p,q, m1,m2,m2_b, ind_1dim,ind_1dim_b

    iam = 0

    wl_omp = Cdzero

    ! global address offset of a cell prior to starting cell.
    icxgo = (nscellx * ipx) / npx
    icygo = (nscelly * ipy) / npy
    iczgo = (nscellz * ipz) / npz
    mbd_z2 = 2*mbd_z
    mbd_y2 = 2*mbd_y
    mbd_x2 = 2*mbd_x

!$omp parallel default(none) &
!$omp& private(iam) &
!$omp& private(load) &
!$omp& private(icx0,icy0,icz0,icx1,icy1,icz1) &
!$omp& private(icx1s,icy1s,icz1s) &
!$omp& private(icx1e,icy1e,icz1e) &
!$omp& private(icxg,icyg,iczg) &
!$omp& private(irx,iry,irz) &
!$omp& private(ngord,nrgord) &
!$omp& private(wm_north,wl_north) &
!$omp& private(j,k,n,p,q, m1,m2,m2_b, ind_1dim,ind_1dim_b) &
!$omp& private(tmp_shml, tmp_ce_r, tmp_ce_i) &
!$omp& shared(shml, sh_rot1, sh_rot2, sh_inv_rot1, sh_inv_rot2) &
!$omp& shared(nmax) &
!$omp& shared(ilevel,lm_length) &
!$omp& shared(wl_omp,wm) &
!$omp& shared(nload,lddir) &
!$omp& shared(nscellx,nscelly,nscellz) &
!$omp& shared(nscxdiv,nscydiv,nsczdiv) &
!$omp& shared(mbd_x ,mbd_y ,mbd_z) &
!$omp& shared(mbd_x2,mbd_y2,mbd_z2) &
!$omp& shared(icxgo,icygo,iczgo)
!$  iam = omp_get_thread_num()
!$omp do schedule(static,1)
    do load = 1, nload(ilevel)
      irx = lddir(1,load,ilevel)     ! x-dir relative cell address of transformation matrix.
      iry = lddir(2,load,ilevel)     ! y-dir relative cell address of transformation matrix.
      irz = lddir(3,load,ilevel)     ! z-dir relative cell address of transformation matrix.

      if(abs(irx).le.2 .and. abs(iry).le.2 .and. abs(irz).le.2) cycle  ! lddir does not include true case.

    do icx0 = 1-mbd_x, nscxdiv+mbd_x      ! target cell x-address.
      if(nscxdiv == 1      .and. 1 - mbd_x /= icx0 .and. icx0 /= nscxdiv + mbd_x) cycle
      if(nscxdiv <  mbd_x2 .and. 1         <= icx0 .and. icx0 <= nscxdiv        ) cycle
      icxg  = modulo(icxgo + icx0 - 1, nscellx) + 1  ! target cell global x-address.
      ngord  = mod(icxg - 1, mbd_x)                  ! sequence number within cell merge group minus 1.
      nrgord = mbd_x - ngord - 1                     ! reverse sequence number within cell merge group minus 1.
      icx1s = icx0 - mbd_x2 - ngord                  ! first cell of the cell range to be transformed.
      icx1e = icx0 + mbd_x2 + nrgord                 ! final cell of the cell range to be transformed.
      if(icx0 <= mbd_x          ) icx1s = icx0       ! plus-direction half-space transformation.
      if(icx0 >  nscxdiv - mbd_x) icx1e = icx0 - 1   ! minus-direction half-space transformation. inter-process exclusive sontrol.
      icx1 = icx0 + irx                              ! real transforming cell defined by relative cell address.
      if(icx1 < icx1s .or. icx1e < icx1) cycle

      do icy0 = 1-mbd_y, nscydiv+mbd_y    ! target cell y-address.
        if(nscydiv == 1     .and. 1 - mbd_y /= icy0 .and. icy0 /= nscydiv + mbd_y) cycle
        if(nscydiv < mbd_y2 .and. 1         <= icy0 .and. icy0 <= nscydiv        ) cycle
        icyg  = modulo(icygo + icy0 - 1, nscelly) + 1
        ngord  = mod(icyg - 1, mbd_y)
        nrgord = mbd_y - ngord - 1
        icy1s = icy0 - mbd_y2 - ngord
        icy1e = icy0 + mbd_y2 + nrgord
        if(icy0 <= mbd_y          ) icy1s = icy0
        if(icy0 >  nscydiv - mbd_y) icy1e = icy0 - 1
        icy1 = icy0 + iry
        if(icy1 < icy1s .or. icy1e < icy1) cycle

        do icz0 = 1-mbd_z, nsczdiv+mbd_z  ! target cell z-address.
          if(nsczdiv == 1     .and. 1 - mbd_z /= icz0 .and. icz0 /= nsczdiv + mbd_z) cycle
          if(nsczdiv < mbd_z2 .and. 1         <= icz0 .and. icz0 <= nsczdiv        ) cycle
          iczg  = modulo(iczgo + icz0 - 1, nscellz) + 1
          ngord  = mod(iczg - 1, mbd_z)
          nrgord = mbd_z - ngord - 1
          icz1s = icz0 - mbd_z2 - ngord
          icz1e = icz0 + mbd_z2 + nrgord
          if(icz0 <= mbd_z          ) icz1s = icz0
          if(icz0 >  nsczdiv - mbd_z) icz1e = icz0 - 1
          icz1 = icz0 + irz
          if(icz1 < icz1s .or. icz1e < icz1) cycle

          include 'tuned_M2L/fmm_far_m2l_rotation_exsimd-reduct.h90'

        end do  ! icz0.
      end do  ! icy0.
    end do  ! icx0.

    end do  ! load.
!$omp end do
!$omp end parallel

    lcbx = nscxdiv + mbd_x2                          ! dimension size of superposition array.
    lcby = nscydiv + mbd_y2
    lcbz = nsczdiv + mbd_z2
    ncall = lcbx*lcby*lcbz                           ! overall candidates of superposition cell.
    ncblk = (ncall - 1)/nomp + 1                     ! block size for thread parallel. 

!$omp parallel default(none) &
!$omp& private(nt,icx0,icy0,icz0,icxyz,m1,iam) &
!$omp& shared(lm_length,wl,wl_omp) &
!$omp& shared(nscxdiv,nscydiv,nsczdiv) &
!$omp& shared(mbd_x ,mbd_y ,mbd_z) &
!$omp& shared(mbd_x2,mbd_y2,mbd_z2) &
!$omp& shared(lcby,lcbz) &
!$omp& shared(ncblk,ncall,nomp)
!$  iam = omp_get_thread_num()
     do nt = 0, nomp-1
       do icxyz = ncblk*iam + 1, min(ncblk*(iam+1), ncall)
        icz0 = mod(icxyz - 1, lcbz) + 1 - mbd_z      ! 1-dimensional to 3-dimensional cell index.
        icy0 = mod(icxyz - 1, lcbz*lcby)
        icy0 = icy0/lcbz + 1 - mbd_y                 ! 1-dimensional to 3-dimensional cell index.
        icx0 = (icxyz - 1)/(lcbz*lcby) + 1 - mbd_x   ! 1-dimensional to 3-dimensional cell index.
        ! exclude vacant target cell.
        if(nscxdiv == 1      .and. 1 - mbd_x /= icx0 .and. icx0 /= nscxdiv + mbd_x) cycle
        if(nscxdiv <  mbd_x2 .and. 1         <= icx0 .and. icx0 <= nscxdiv        ) cycle
        if(nscydiv == 1      .and. 1 - mbd_y /= icy0 .and. icy0 /= nscydiv + mbd_y) cycle
        if(nscydiv <  mbd_y2 .and. 1         <= icy0 .and. icy0 <= nscydiv        ) cycle
        if(nsczdiv == 1      .and. 1 - mbd_z /= icz0 .and. icz0 /= nsczdiv + mbd_z) cycle
        if(nsczdiv <  mbd_z2 .and. 1         <= icz0 .and. icz0 <= nsczdiv        ) cycle

        do m1=1,lm_length
          wl(m1,1,icz0,icy0,icx0) = &
     &    wl(m1,1,icz0,icy0,icx0) + wl_omp(m1,1,icz0,icy0,icx0,nt)
          wl(m1,2,icz0,icy0,icx0) = &
     &    wl(m1,2,icz0,icy0,icx0) + wl_omp(m1,2,icz0,icy0,icx0,nt)
        enddo

      end do  ! icxyz.
    enddo ! nt.
!$omp end parallel

  end subroutine M2L_lower_level_dr_general
!---------------------------------------------------------------------
#endif  /* Not MOMENT */
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!>
!! \brief  Subroutine to execute M2L (with global arrays, upper level)
!! \author Noriyuki Yoshii, Yoshimichi Andoh, Shin-chi Ichikawa
!<
  subroutine M2L_upper_level(ilevel, wm, wl, nomp)
!---------------------------------------------------------------------
! Main routine of M2L_upper_level, handling rotation or not, 
! nmax=4 or general, multiple cell pairs or single cell pair.
! M2L_upper_level is for global type array wm dedicated to communication.
!
    implicit none
    integer(kind=wip), intent(in)  :: ilevel
    real(kind=wmp), intent(in)     :: wm(lm_length, 2, &
  &                                      wm_size_z(ilevel), &
  &                                      wm_size_y(ilevel), &
  &                                      wm_size_x(ilevel))
    real(kind=wvap), intent(inout) :: wl(lm_length, 2, &
  &                                      wl_size_z(ilevel), &
  &                                      wl_size_y(ilevel), &
  &                                      wl_size_x(ilevel))
    integer(kind=wip) :: nomp

#ifdef MOMENT    /* multiple cell pairs */
    if(nmax==4) then    !/*---- for Matrix Order 15 ----*/

      call M2L_upper_level_mpairs_nmax4(ilevel, wm, wl, nomp)

    else                !/*---- for General Value of NMAX ----*/

      call M2L_upper_level_mpairs(ilevel, wm, wl, nomp)

    endif               !/*---- for General Value of NMAX ----*/
#else   /* Not MOMENT */
    if(nmax==4) then    !/*---- for Matrix Order 15 ----*/

      call M2L_upper_level_nmax4(ilevel, wm, wl, nomp)

    else                !/*---- for General Value of NMAX ----*/

      call M2L_upper_level_general(ilevel, wm, wl, nomp)

    endif
#endif  /* Not MOMENT */

  end subroutine M2L_upper_level

#ifdef MOMENT    /* multiple cell pairs */
!---------------------------------------------------------------------
!>
!! \brief  Subroutine to execute M2L (with global arrays, upper level)
!! \author Noriyuki Yoshii, Yoshimichi Andoh, Shin-chi Ichikawa
!<
  subroutine M2L_upper_level_mpairs_nmax4(ilevel, wm, wl, nomp)
!---------------------------------------------------------------------
! symmetry and rotation code.
! pipelined multiple cell pairs and unrolled for nmax=4.
!
    use mpi_3d_grid, only : npx, npy, npz, ipx, ipy, ipz
    implicit none
    integer(kind=wip), intent(in)  :: ilevel
    real(kind=wmp), intent(in)     :: wm(lm_length, 2, &
  &                                         wm_size_z(ilevel), &
  &                                         wm_size_y(ilevel), &
  &                                         wm_size_x(ilevel))
    real(kind=wvap), intent(inout) :: wl(lm_length, 2, &
  &                                         wl_size_z(ilevel), &
  &                                         wl_size_y(ilevel), &
  &                                         wl_size_x(ilevel))
    integer(kind=wip) :: iall
    integer(kind=wip) :: icx0, icy0, icz0
    integer(kind=wip) :: icx1, icy1, icz1
    integer(kind=wip) :: irx, iry, irz
    integer(kind=wip) :: nomp, iam, iiam
    integer(kind=wip) :: nsczdiv,nscydiv,nscxdiv
! add padding(+5,,+2,,,).
    real(kind=wvap) :: wl_omp(lm_length+5, 2, &
                                    wl_size_z(ilevel)+2, &
                                    wl_size_y(ilevel)  , &
                                    wl_size_x(ilevel)  , 0:nomp-1)

!add for simultaneous multiple cell pairs.
    real(kind=wvp) :: wwm_north(lm_length,2,8), wwl_north(lm_length,2,8)
    integer(kind=wip) :: icp
    integer(kind=wip) :: icx0_1,icy0_1,icz0_1,icx1_1,icy1_1,icz1_1
    integer(kind=wip) :: icx0_2,icy0_2,icz0_2,icx1_2,icy1_2,icz1_2
    integer(kind=wip) :: icx0_3,icy0_3,icz0_3,icx1_3,icy1_3,icz1_3
    integer(kind=wip) :: icx0_4,icy0_4,icz0_4,icx1_4,icy1_4,icz1_4
    integer(kind=wip) :: icx0_5,icy0_5,icz0_5,icx1_5,icy1_5,icz1_5
    integer(kind=wip) :: icx0_6,icy0_6,icz0_6,icx1_6,icy1_6,icz1_6
    integer(kind=wip) :: icx0_7,icy0_7,icz0_7,icx1_7,icy1_7,icz1_7
    integer(kind=wip) :: icx0_8,icy0_8,icz0_8,icx1_8,icy1_8,icz1_8
    real(kind=wvp) :: shml_n0, shml_n1, shml_n2, shml_n3, shml_n4
    real(kind=wvp) :: shrot1_q0_r, shrot1_q0_i, shrot2_q0_r, shrot2_q0_i
    real(kind=wvp) :: shrot1_q1_r, shrot1_q1_i, shrot2_q1_r, shrot2_q1_i
    real(kind=wvp) :: shrot1_q2_r, shrot1_q2_i, shrot2_q2_r, shrot2_q2_i
    real(kind=wvp) :: shrot1_q3_r, shrot1_q3_i, shrot2_q3_r, shrot2_q3_i
    real(kind=wvp) :: shrot1_q4_r, shrot1_q4_i, shrot2_q4_r, shrot2_q4_i
    real(kind=wvp) :: shinvrot1_q0_r, shinvrot1_q0_i, shinvrot2_q0_r, shinvrot2_q0_i
    real(kind=wvp) :: shinvrot1_q1_r, shinvrot1_q1_i, shinvrot2_q1_r, shinvrot2_q1_i
    real(kind=wvp) :: shinvrot1_q2_r, shinvrot1_q2_i, shinvrot2_q2_r, shinvrot2_q2_i
    real(kind=wvp) :: shinvrot1_q3_r, shinvrot1_q3_i, shinvrot2_q3_r, shinvrot2_q3_i
    real(kind=wvp) :: shinvrot1_q4_r, shinvrot1_q4_i, shinvrot2_q4_r, shinvrot2_q4_i
    real(kind=wvp) :: wm_q0_1r, wm_q0_1i, wm_q0_2r, wm_q0_2i, wm_q0_3r, wm_q0_3i, wm_q0_4r, wm_q0_4i
    real(kind=wvp) :: wm_q1_1r, wm_q1_1i, wm_q1_2r, wm_q1_2i, wm_q1_3r, wm_q1_3i, wm_q1_4r, wm_q1_4i
    real(kind=wvp) :: wm_q2_1r, wm_q2_1i, wm_q2_2r, wm_q2_2i, wm_q2_3r, wm_q2_3i, wm_q2_4r, wm_q2_4i
    real(kind=wvp) :: wm_q3_1r, wm_q3_1i, wm_q3_2r, wm_q3_2i, wm_q3_3r, wm_q3_3i, wm_q3_4r, wm_q3_4i
    real(kind=wvp) :: wm_q4_1r, wm_q4_1i, wm_q4_2r, wm_q4_2i, wm_q4_3r, wm_q4_3i, wm_q4_4r, wm_q4_4i
    real(kind=wvp) :: wm_q0_5r, wm_q0_5i, wm_q0_6r, wm_q0_6i, wm_q0_7r, wm_q0_7i, wm_q0_8r, wm_q0_8i
    real(kind=wvp) :: wm_q1_5r, wm_q1_5i, wm_q1_6r, wm_q1_6i, wm_q1_7r, wm_q1_7i, wm_q1_8r, wm_q1_8i
    real(kind=wvp) :: wm_q2_5r, wm_q2_5i, wm_q2_6r, wm_q2_6i, wm_q2_7r, wm_q2_7i, wm_q2_8r, wm_q2_8i
    real(kind=wvp) :: wm_q3_5r, wm_q3_5i, wm_q3_6r, wm_q3_6i, wm_q3_7r, wm_q3_7i, wm_q3_8r, wm_q3_8i
    real(kind=wvp) :: wm_q4_5r, wm_q4_5i, wm_q4_6r, wm_q4_6i, wm_q4_7r, wm_q4_7i, wm_q4_8r, wm_q4_8i
    real(kind=wvp) :: wlnorth_q0_1r, wlnorth_q0_1i, wlnorth_q0_2r, wlnorth_q0_2i, &
   &                  wlnorth_q0_3r, wlnorth_q0_3i, wlnorth_q0_4r, wlnorth_q0_4i
    real(kind=wvp) :: wlnorth_q1_1r, wlnorth_q1_1i, wlnorth_q1_2r, wlnorth_q1_2i, &
   &                  wlnorth_q1_3r, wlnorth_q1_3i, wlnorth_q1_4r, wlnorth_q1_4i
    real(kind=wvp) :: wlnorth_q2_1r, wlnorth_q2_1i, wlnorth_q2_2r, wlnorth_q2_2i, &
   &                  wlnorth_q2_3r, wlnorth_q2_3i, wlnorth_q2_4r, wlnorth_q2_4i
    real(kind=wvp) :: wlnorth_q3_1r, wlnorth_q3_1i, wlnorth_q3_2r, wlnorth_q3_2i, &
   &                  wlnorth_q3_3r, wlnorth_q3_3i, wlnorth_q3_4r, wlnorth_q3_4i
    real(kind=wvp) :: wlnorth_q4_1r, wlnorth_q4_1i, wlnorth_q4_2r, wlnorth_q4_2i, &
   &                  wlnorth_q4_3r, wlnorth_q4_3i, wlnorth_q4_4r, wlnorth_q4_4i
    real(kind=wvp) :: wlnorth_q0_5r, wlnorth_q0_5i, wlnorth_q0_6r, wlnorth_q0_6i, &
   &                  wlnorth_q0_7r, wlnorth_q0_7i, wlnorth_q0_8r, wlnorth_q0_8i
    real(kind=wvp) :: wlnorth_q1_5r, wlnorth_q1_5i, wlnorth_q1_6r, wlnorth_q1_6i, &
   &                  wlnorth_q1_7r, wlnorth_q1_7i, wlnorth_q1_8r, wlnorth_q1_8i
    real(kind=wvp) :: wlnorth_q2_5r, wlnorth_q2_5i, wlnorth_q2_6r, wlnorth_q2_6i, &
   &                  wlnorth_q2_7r, wlnorth_q2_7i, wlnorth_q2_8r, wlnorth_q2_8i
    real(kind=wvp) :: wlnorth_q3_5r, wlnorth_q3_5i, wlnorth_q3_6r, wlnorth_q3_6i, &
   &                  wlnorth_q3_7r, wlnorth_q3_7i, wlnorth_q3_8r, wlnorth_q3_8i
    real(kind=wvp) :: wlnorth_q4_5r, wlnorth_q4_5i, wlnorth_q4_6r, wlnorth_q4_6i, &
   &                  wlnorth_q4_7r, wlnorth_q4_7i, wlnorth_q4_8r, wlnorth_q4_8i
    real(kind=wvp) :: tmp_1r, tmp_1i, tmp_2r, tmp_2i, tmp_3r, tmp_3i, tmp_4r, tmp_4i
    real(kind=wvp) :: tmp_5r, tmp_5i, tmp_6r, tmp_6i, tmp_7r, tmp_7i, tmp_8r, tmp_8i
    real(kind=wvp) :: wm_north(lm_length,2), wl_north(lm_length,2)
    real(kind=wvp) :: tmp_shml, tmp_ce_r, tmp_ce_i
    integer(kind=wip) :: k,p, m1
    !debug
    !     integer(kind=wip)::numjc=0
    !debug

!#ifndef _OPENMP
    iam = 0
!#endif

    wl_omp = Cdzero

    nsczdiv = wl_size_z(ilevel);  nscydiv = wl_size_y(ilevel);  nscxdiv = wl_size_x(ilevel)

!$omp parallel default(none) &
!$omp& private(iam) &
!$omp& private(icx0,icy0,icz0,icx1,icy1,icz1) &
!$omp& private(irx,iry,irz) &
!$omp& private(wwm_north,wwl_north) &
!$omp& private(icp) &
!$omp& private(icx0_1,icy0_1,icz0_1,icx1_1,icy1_1,icz1_1) &
!$omp& private(icx0_2,icy0_2,icz0_2,icx1_2,icy1_2,icz1_2) &
!$omp& private(icx0_3,icy0_3,icz0_3,icx1_3,icy1_3,icz1_3) &
!$omp& private(icx0_4,icy0_4,icz0_4,icx1_4,icy1_4,icz1_4) &
!$omp& private(icx0_5,icy0_5,icz0_5,icx1_5,icy1_5,icz1_5) &
!$omp& private(icx0_6,icy0_6,icz0_6,icx1_6,icy1_6,icz1_6) &
!$omp& private(icx0_7,icy0_7,icz0_7,icx1_7,icy1_7,icz1_7) &
!$omp& private(icx0_8,icy0_8,icz0_8,icx1_8,icy1_8,icz1_8) &
!$omp& private(shml_n0, shml_n1, shml_n2, shml_n3, shml_n4) &
!$omp& private(shrot1_q0_r, shrot1_q0_i, shrot2_q0_r, shrot2_q0_i) &
!$omp& private(shrot1_q1_r, shrot1_q1_i, shrot2_q1_r, shrot2_q1_i) &
!$omp& private(shrot1_q2_r, shrot1_q2_i, shrot2_q2_r, shrot2_q2_i) &
!$omp& private(shrot1_q3_r, shrot1_q3_i, shrot2_q3_r, shrot2_q3_i) &
!$omp& private(shrot1_q4_r, shrot1_q4_i, shrot2_q4_r, shrot2_q4_i) &
!$omp& private(shinvrot1_q0_r, shinvrot1_q0_i, shinvrot2_q0_r, shinvrot2_q0_i) &
!$omp& private(shinvrot1_q1_r, shinvrot1_q1_i, shinvrot2_q1_r, shinvrot2_q1_i) &
!$omp& private(shinvrot1_q2_r, shinvrot1_q2_i, shinvrot2_q2_r, shinvrot2_q2_i) &
!$omp& private(shinvrot1_q3_r, shinvrot1_q3_i, shinvrot2_q3_r, shinvrot2_q3_i) &
!$omp& private(shinvrot1_q4_r, shinvrot1_q4_i, shinvrot2_q4_r, shinvrot2_q4_i) &
!$omp& private(wm_q0_1r, wm_q0_1i, wm_q0_2r, wm_q0_2i, wm_q0_3r, wm_q0_3i, wm_q0_4r, wm_q0_4i) &
!$omp& private(wm_q1_1r, wm_q1_1i, wm_q1_2r, wm_q1_2i, wm_q1_3r, wm_q1_3i, wm_q1_4r, wm_q1_4i) &
!$omp& private(wm_q2_1r, wm_q2_1i, wm_q2_2r, wm_q2_2i, wm_q2_3r, wm_q2_3i, wm_q2_4r, wm_q2_4i) &
!$omp& private(wm_q3_1r, wm_q3_1i, wm_q3_2r, wm_q3_2i, wm_q3_3r, wm_q3_3i, wm_q3_4r, wm_q3_4i) &
!$omp& private(wm_q4_1r, wm_q4_1i, wm_q4_2r, wm_q4_2i, wm_q4_3r, wm_q4_3i, wm_q4_4r, wm_q4_4i) &
!$omp& private(wm_q0_5r, wm_q0_5i, wm_q0_6r, wm_q0_6i, wm_q0_7r, wm_q0_7i, wm_q0_8r, wm_q0_8i) &
!$omp& private(wm_q1_5r, wm_q1_5i, wm_q1_6r, wm_q1_6i, wm_q1_7r, wm_q1_7i, wm_q1_8r, wm_q1_8i) &
!$omp& private(wm_q2_5r, wm_q2_5i, wm_q2_6r, wm_q2_6i, wm_q2_7r, wm_q2_7i, wm_q2_8r, wm_q2_8i) &
!$omp& private(wm_q3_5r, wm_q3_5i, wm_q3_6r, wm_q3_6i, wm_q3_7r, wm_q3_7i, wm_q3_8r, wm_q3_8i) &
!$omp& private(wm_q4_5r, wm_q4_5i, wm_q4_6r, wm_q4_6i, wm_q4_7r, wm_q4_7i, wm_q4_8r, wm_q4_8i) &
!$omp& private(wlnorth_q0_1r, wlnorth_q0_1i, wlnorth_q0_2r, wlnorth_q0_2i) &
!$omp& private(wlnorth_q0_3r, wlnorth_q0_3i, wlnorth_q0_4r, wlnorth_q0_4i) &
!$omp& private(wlnorth_q1_1r, wlnorth_q1_1i, wlnorth_q1_2r, wlnorth_q1_2i) &
!$omp& private(wlnorth_q1_3r, wlnorth_q1_3i, wlnorth_q1_4r, wlnorth_q1_4i) &
!$omp& private(wlnorth_q2_1r, wlnorth_q2_1i, wlnorth_q2_2r, wlnorth_q2_2i) &
!$omp& private(wlnorth_q2_3r, wlnorth_q2_3i, wlnorth_q2_4r, wlnorth_q2_4i) &
!$omp& private(wlnorth_q3_1r, wlnorth_q3_1i, wlnorth_q3_2r, wlnorth_q3_2i) &
!$omp& private(wlnorth_q3_3r, wlnorth_q3_3i, wlnorth_q3_4r, wlnorth_q3_4i) &
!$omp& private(wlnorth_q4_1r, wlnorth_q4_1i, wlnorth_q4_2r, wlnorth_q4_2i) &
!$omp& private(wlnorth_q4_3r, wlnorth_q4_3i, wlnorth_q4_4r, wlnorth_q4_4i) &
!$omp& private(wlnorth_q0_5r, wlnorth_q0_5i, wlnorth_q0_6r, wlnorth_q0_6i) &
!$omp& private(wlnorth_q0_7r, wlnorth_q0_7i, wlnorth_q0_8r, wlnorth_q0_8i) &
!$omp& private(wlnorth_q1_5r, wlnorth_q1_5i, wlnorth_q1_6r, wlnorth_q1_6i) &
!$omp& private(wlnorth_q1_7r, wlnorth_q1_7i, wlnorth_q1_8r, wlnorth_q1_8i) &
!$omp& private(wlnorth_q2_5r, wlnorth_q2_5i, wlnorth_q2_6r, wlnorth_q2_6i) &
!$omp& private(wlnorth_q2_7r, wlnorth_q2_7i, wlnorth_q2_8r, wlnorth_q2_8i) &
!$omp& private(wlnorth_q3_5r, wlnorth_q3_5i, wlnorth_q3_6r, wlnorth_q3_6i) &
!$omp& private(wlnorth_q3_7r, wlnorth_q3_7i, wlnorth_q3_8r, wlnorth_q3_8i) &
!$omp& private(wlnorth_q4_5r, wlnorth_q4_5i, wlnorth_q4_6r, wlnorth_q4_6i) &
!$omp& private(wlnorth_q4_7r, wlnorth_q4_7i, wlnorth_q4_8r, wlnorth_q4_8i) &
!$omp& private(tmp_1r, tmp_1i, tmp_2r, tmp_2i, tmp_3r, tmp_3i, tmp_4r, tmp_4i) &
!$omp& private(tmp_5r, tmp_5i, tmp_6r, tmp_6i, tmp_7r, tmp_7i, tmp_8r, tmp_8i) &
!$omp& private(wm_north,wl_north) &
!$omp& private(k,p, m1) &
!$omp& private(tmp_shml, tmp_ce_r, tmp_ce_i) &
!$omp& shared(shml, sh_rot1, sh_rot2, sh_inv_rot1, sh_inv_rot2) &
!$omp& shared(nmax) &
!debug
!!!$omp& reduction(+:numjc) &
!debug
!$omp& shared(ilevel) &
!$omp& shared(wl_omp,wm) &
!$omp& shared(cp8_count,cp6_count,cp4_count,cp3_count,cp2_count,cp1_count) &
!$omp& shared(cp8_list0,cp6_list0,cp4_list0,cp3_list0,cp2_list0,cp1_list0) &
!$omp& shared(cp8_list1,cp6_list1,cp4_list1,cp3_list1,cp2_list1,cp1_list1) &
!$omp& shared(cp8_irxyz,cp6_irxyz,cp4_irxyz,cp3_irxyz,cp2_irxyz,cp1_irxyz) &
!$omp& shared(lm_length) &
!$omp& shared(nchunk)
!$  iam = omp_get_thread_num()
!!$omp do schedule(static,nchunk(ilevel))

    if(cp8_count(ilevel) > 0) then    !! /*---- pairs multiple of 8 for each relative cell address ----*/

!$omp do
    do icp = 1, cp8_count(ilevel), 8
      irx = cp8_irxyz(1,icp,ilevel)  ! x-dir relative cell address of transformation matrix.
      iry = cp8_irxyz(2,icp,ilevel)  ! y-dir relative cell address of transformation matrix.
      irz = cp8_irxyz(3,icp,ilevel)  ! z-dir relative cell address of transformation matrix.

    !*** multipole to local translation
      icx0_1 = cp8_list0(1,icp  ,ilevel) ; icy0_1 = cp8_list0(2,icp  ,ilevel) ; icz0_1 = cp8_list0(3,icp  ,ilevel)
      icx1_1 = cp8_list1(1,icp  ,ilevel) ; icy1_1 = cp8_list1(2,icp  ,ilevel) ; icz1_1 = cp8_list1(3,icp  ,ilevel)
      icx0_2 = cp8_list0(1,icp+1,ilevel) ; icy0_2 = cp8_list0(2,icp+1,ilevel) ; icz0_2 = cp8_list0(3,icp+1,ilevel)
      icx1_2 = cp8_list1(1,icp+1,ilevel) ; icy1_2 = cp8_list1(2,icp+1,ilevel) ; icz1_2 = cp8_list1(3,icp+1,ilevel)
      icx0_3 = cp8_list0(1,icp+2,ilevel) ; icy0_3 = cp8_list0(2,icp+2,ilevel) ; icz0_3 = cp8_list0(3,icp+2,ilevel)
      icx1_3 = cp8_list1(1,icp+2,ilevel) ; icy1_3 = cp8_list1(2,icp+2,ilevel) ; icz1_3 = cp8_list1(3,icp+2,ilevel)
      icx0_4 = cp8_list0(1,icp+3,ilevel) ; icy0_4 = cp8_list0(2,icp+3,ilevel) ; icz0_4 = cp8_list0(3,icp+3,ilevel)
      icx1_4 = cp8_list1(1,icp+3,ilevel) ; icy1_4 = cp8_list1(2,icp+3,ilevel) ; icz1_4 = cp8_list1(3,icp+3,ilevel)
      icx0_5 = cp8_list0(1,icp+4,ilevel) ; icy0_5 = cp8_list0(2,icp+4,ilevel) ; icz0_5 = cp8_list0(3,icp+4,ilevel)
      icx1_5 = cp8_list1(1,icp+4,ilevel) ; icy1_5 = cp8_list1(2,icp+4,ilevel) ; icz1_5 = cp8_list1(3,icp+4,ilevel)
      icx0_6 = cp8_list0(1,icp+5,ilevel) ; icy0_6 = cp8_list0(2,icp+5,ilevel) ; icz0_6 = cp8_list0(3,icp+5,ilevel)
      icx1_6 = cp8_list1(1,icp+5,ilevel) ; icy1_6 = cp8_list1(2,icp+5,ilevel) ; icz1_6 = cp8_list1(3,icp+5,ilevel)
      icx0_7 = cp8_list0(1,icp+6,ilevel) ; icy0_7 = cp8_list0(2,icp+6,ilevel) ; icz0_7 = cp8_list0(3,icp+6,ilevel)
      icx1_7 = cp8_list1(1,icp+6,ilevel) ; icy1_7 = cp8_list1(2,icp+6,ilevel) ; icz1_7 = cp8_list1(3,icp+6,ilevel)
      icx0_8 = cp8_list0(1,icp+7,ilevel) ; icy0_8 = cp8_list0(2,icp+7,ilevel) ; icz0_8 = cp8_list0(3,icp+7,ilevel)
      icx1_8 = cp8_list1(1,icp+7,ilevel) ; icy1_8 = cp8_list1(2,icp+7,ilevel) ; icz1_8 = cp8_list1(3,icp+7,ilevel)

      include 'tuned_M2L/fmm_far_m2l_rotation_mpair8_nmax4_exsimd-reduct_unroll.h90'
    end do  ! cp8_count.
!$omp end do

    endif    ! cp8_count>0.

    if(cp6_count(ilevel) > 0) then    !! /*---- pairs multiple of 6 for each relative cell address ----*/

!$omp do
    do icp = 1, cp6_count(ilevel), 6
      irx = cp6_irxyz(1,icp,ilevel)  ! x-dir relative cell address of transformation matrix.
      iry = cp6_irxyz(2,icp,ilevel)  ! y-dir relative cell address of transformation matrix.
      irz = cp6_irxyz(3,icp,ilevel)  ! z-dir relative cell address of transformation matrix.

    !*** multipole to local translation
      icx0_1 = cp6_list0(1,icp  ,ilevel) ; icy0_1 = cp6_list0(2,icp  ,ilevel) ; icz0_1 = cp6_list0(3,icp  ,ilevel)
      icx1_1 = cp6_list1(1,icp  ,ilevel) ; icy1_1 = cp6_list1(2,icp  ,ilevel) ; icz1_1 = cp6_list1(3,icp  ,ilevel)
      icx0_2 = cp6_list0(1,icp+1,ilevel) ; icy0_2 = cp6_list0(2,icp+1,ilevel) ; icz0_2 = cp6_list0(3,icp+1,ilevel)
      icx1_2 = cp6_list1(1,icp+1,ilevel) ; icy1_2 = cp6_list1(2,icp+1,ilevel) ; icz1_2 = cp6_list1(3,icp+1,ilevel)
      icx0_3 = cp6_list0(1,icp+2,ilevel) ; icy0_3 = cp6_list0(2,icp+2,ilevel) ; icz0_3 = cp6_list0(3,icp+2,ilevel)
      icx1_3 = cp6_list1(1,icp+2,ilevel) ; icy1_3 = cp6_list1(2,icp+2,ilevel) ; icz1_3 = cp6_list1(3,icp+2,ilevel)
      icx0_4 = cp6_list0(1,icp+3,ilevel) ; icy0_4 = cp6_list0(2,icp+3,ilevel) ; icz0_4 = cp6_list0(3,icp+3,ilevel)
      icx1_4 = cp6_list1(1,icp+3,ilevel) ; icy1_4 = cp6_list1(2,icp+3,ilevel) ; icz1_4 = cp6_list1(3,icp+3,ilevel)
      icx0_5 = cp6_list0(1,icp+4,ilevel) ; icy0_5 = cp6_list0(2,icp+4,ilevel) ; icz0_5 = cp6_list0(3,icp+4,ilevel)
      icx1_5 = cp6_list1(1,icp+4,ilevel) ; icy1_5 = cp6_list1(2,icp+4,ilevel) ; icz1_5 = cp6_list1(3,icp+4,ilevel)
      icx0_6 = cp6_list0(1,icp+5,ilevel) ; icy0_6 = cp6_list0(2,icp+5,ilevel) ; icz0_6 = cp6_list0(3,icp+5,ilevel)
      icx1_6 = cp6_list1(1,icp+5,ilevel) ; icy1_6 = cp6_list1(2,icp+5,ilevel) ; icz1_6 = cp6_list1(3,icp+5,ilevel)

      include 'tuned_M2L/fmm_far_m2l_rotation_mpair6_nmax4_exsimd-reduct_unroll.h90'
    end do  ! cp6_count.
!$omp end do

    endif    ! cp6_count>0.

    if(cp4_count(ilevel) > 0) then    !! /*---- pairs multiple of 4 for each relative cell address ----*/

!$omp do
    do icp = 1, cp4_count(ilevel), 4
      irx = cp4_irxyz(1,icp,ilevel)  ! x-dir relative cell address of transformation matrix.
      iry = cp4_irxyz(2,icp,ilevel)  ! y-dir relative cell address of transformation matrix.
      irz = cp4_irxyz(3,icp,ilevel)  ! z-dir relative cell address of transformation matrix.

    !*** multipole to local translation
      icx0_1 = cp4_list0(1,icp  ,ilevel) ; icy0_1 = cp4_list0(2,icp  ,ilevel) ; icz0_1 = cp4_list0(3,icp  ,ilevel)
      icx1_1 = cp4_list1(1,icp  ,ilevel) ; icy1_1 = cp4_list1(2,icp  ,ilevel) ; icz1_1 = cp4_list1(3,icp  ,ilevel)
      icx0_2 = cp4_list0(1,icp+1,ilevel) ; icy0_2 = cp4_list0(2,icp+1,ilevel) ; icz0_2 = cp4_list0(3,icp+1,ilevel)
      icx1_2 = cp4_list1(1,icp+1,ilevel) ; icy1_2 = cp4_list1(2,icp+1,ilevel) ; icz1_2 = cp4_list1(3,icp+1,ilevel)
      icx0_3 = cp4_list0(1,icp+2,ilevel) ; icy0_3 = cp4_list0(2,icp+2,ilevel) ; icz0_3 = cp4_list0(3,icp+2,ilevel)
      icx1_3 = cp4_list1(1,icp+2,ilevel) ; icy1_3 = cp4_list1(2,icp+2,ilevel) ; icz1_3 = cp4_list1(3,icp+2,ilevel)
      icx0_4 = cp4_list0(1,icp+3,ilevel) ; icy0_4 = cp4_list0(2,icp+3,ilevel) ; icz0_4 = cp4_list0(3,icp+3,ilevel)
      icx1_4 = cp4_list1(1,icp+3,ilevel) ; icy1_4 = cp4_list1(2,icp+3,ilevel) ; icz1_4 = cp4_list1(3,icp+3,ilevel)

      include 'tuned_M2L/fmm_far_m2l_rotation_mpair4_nmax4_exsimd-reduct_unroll.h90'
    end do  ! cp4_count.
!$omp end do

    endif    ! cp4_count>0.

    if(cp3_count(ilevel) > 0) then    !! /*---- pairs multiple of 3 for each relative cell address ----*/

!$omp do
    do icp = 1, cp3_count(ilevel), 3
      irx = cp3_irxyz(1,icp,ilevel)  ! x-dir relative cell address of transformation matrix.
      iry = cp3_irxyz(2,icp,ilevel)  ! y-dir relative cell address of transformation matrix.
      irz = cp3_irxyz(3,icp,ilevel)  ! z-dir relative cell address of transformation matrix.

    !*** multipole to local translation
      icx0_1 = cp3_list0(1,icp  ,ilevel) ; icy0_1 = cp3_list0(2,icp  ,ilevel) ; icz0_1 = cp3_list0(3,icp  ,ilevel)
      icx1_1 = cp3_list1(1,icp  ,ilevel) ; icy1_1 = cp3_list1(2,icp  ,ilevel) ; icz1_1 = cp3_list1(3,icp  ,ilevel)
      icx0_2 = cp3_list0(1,icp+1,ilevel) ; icy0_2 = cp3_list0(2,icp+1,ilevel) ; icz0_2 = cp3_list0(3,icp+1,ilevel)
      icx1_2 = cp3_list1(1,icp+1,ilevel) ; icy1_2 = cp3_list1(2,icp+1,ilevel) ; icz1_2 = cp3_list1(3,icp+1,ilevel)
      icx0_3 = cp3_list0(1,icp+2,ilevel) ; icy0_3 = cp3_list0(2,icp+2,ilevel) ; icz0_3 = cp3_list0(3,icp+2,ilevel)
      icx1_3 = cp3_list1(1,icp+2,ilevel) ; icy1_3 = cp3_list1(2,icp+2,ilevel) ; icz1_3 = cp3_list1(3,icp+2,ilevel)

      include 'tuned_M2L/fmm_far_m2l_rotation_mpair3_nmax4_exsimd-reduct_unroll.h90'
    end do  ! cp3_count.
!$omp end do

    endif    ! cp3_count>0.

    if(cp2_count(ilevel) > 0) then    !! /*---- pairs multiple of 2 for each relative cell address ----*/

!$omp do
    do icp = 1, cp2_count(ilevel), 2
      irx = cp2_irxyz(1,icp,ilevel)  ! x-dir relative cell address of transformation matrix.
      iry = cp2_irxyz(2,icp,ilevel)  ! y-dir relative cell address of transformation matrix.
      irz = cp2_irxyz(3,icp,ilevel)  ! z-dir relative cell address of transformation matrix.

    !*** multipole to local translation
      icx0_1 = cp2_list0(1,icp  ,ilevel) ; icy0_1 = cp2_list0(2,icp  ,ilevel) ; icz0_1 = cp2_list0(3,icp  ,ilevel)
      icx1_1 = cp2_list1(1,icp  ,ilevel) ; icy1_1 = cp2_list1(2,icp  ,ilevel) ; icz1_1 = cp2_list1(3,icp  ,ilevel)
      icx0_2 = cp2_list0(1,icp+1,ilevel) ; icy0_2 = cp2_list0(2,icp+1,ilevel) ; icz0_2 = cp2_list0(3,icp+1,ilevel)
      icx1_2 = cp2_list1(1,icp+1,ilevel) ; icy1_2 = cp2_list1(2,icp+1,ilevel) ; icz1_2 = cp2_list1(3,icp+1,ilevel)

      include 'tuned_M2L/fmm_far_m2l_rotation_mpair2_nmax4_exsimd-reduct_unroll.h90'
    end do  ! cp2_count.
!$omp end do

    endif    ! cp2_count>0.

    if(cp1_count(ilevel) > 0) then    !! /*---- pairs multiple of 1 for each relative cell address ----*/

!$omp do
    do icp = 1, cp1_count(ilevel)
      irx = cp1_irxyz(1,icp,ilevel)  ! x-dir relative cell address of transformation matrix.
      iry = cp1_irxyz(2,icp,ilevel)  ! y-dir relative cell address of transformation matrix.
      irz = cp1_irxyz(3,icp,ilevel)  ! z-dir relative cell address of transformation matrix.

    !*** multipole to local translation
      icx0 = cp1_list0(1,icp  ,ilevel) ; icy0 = cp1_list0(2,icp  ,ilevel) ; icz0 = cp1_list0(3,icp  ,ilevel)
      icx1 = cp1_list1(1,icp  ,ilevel) ; icy1 = cp1_list1(2,icp  ,ilevel) ; icz1 = cp1_list1(3,icp  ,ilevel)

      include 'tuned_M2L/fmm_far_m2l_rotation_nmax4_exsimd-reduct_unroll.h90'
    end do  ! cp1_count.
!$omp end do

    endif    ! cp1_count>0.

!$omp end parallel


!$omp parallel default(none) &
!$omp& private(iall,icz0,icy0,icx0) &
!$omp& private(m1) &
!$omp& private(iiam) &
!$omp& shared(lm_length) &
!$omp& shared(nscxdiv,nscydiv,nsczdiv) &
!$omp& shared(wl,wl_omp) &
!$omp& shared(nomp)
!$omp do
    DO iall=0, nscxdiv * nscydiv * nsczdiv - 1
       icz0 = mod(iall,nsczdiv) + 1
       icy0 = iall/nsczdiv
       icx0 = icy0/nscydiv + 1
       icy0 = mod(icy0,nscydiv) + 1
       do iiam=0,nomp-1
!c !$omp do
       do m1=1,lm_length
!c          do iiam=0,nomp-1
             wl(m1,1,icz0,icy0,icx0) = wl(m1,1,icz0,icy0,icx0) + wl_omp(m1,1,icz0,icy0,icx0,iiam)
             wl(m1,2,icz0,icy0,icx0) = wl(m1,2,icz0,icy0,icx0) + wl_omp(m1,2,icz0,icy0,icx0,iiam)
!c          enddo
       enddo
       enddo
!c !$omp end do nowait
    ENDDO
!$omp end do
!$omp end parallel

#ifdef DEBUG_MTDFMM
    DO iall=0, nscxdiv * nscydiv * nsczdiv - 1
       icz0 = mod(iall,nsczdiv) + 1
       icy0 = iall/nsczdiv
       icx0 = icy0/nscydiv + 1
       icy0 = mod(icy0,nscydiv) + 1
       write(myrank+ilevel*1000+40000,*) wl(:,1,icz0,icy0,icx0)
       write(myrank+ilevel*1000+40000,*) wl(:,2,icz0,icy0,icx0)
    ENDDO
    call flush(myrank+ilevel*1000+40000)
#endif

  end subroutine M2L_upper_level_mpairs_nmax4
!---------------------------------------------------------------------
!>
!! \brief  Subroutine to execute M2L (with global arrays, upper level)
!! \author Noriyuki Yoshii, Yoshimichi Andoh, Shin-chi Ichikawa
!<
  subroutine M2L_upper_level_mpairs(ilevel, wm, wl, nomp)
!---------------------------------------------------------------------
! symmetry and rotation code.
! pipelined multiple cell pairs. for general value of nmax.
!
    implicit none
    integer(kind=wip), intent(in)  :: ilevel
    real(kind=wmp), intent(in)     :: wm(lm_length, 2, &
  &                                         wm_size_z(ilevel), &
  &                                         wm_size_y(ilevel), &
  &                                         wm_size_x(ilevel))
    real(kind=wvap), intent(inout) :: wl(lm_length, 2, &
  &                                         wl_size_z(ilevel), &
  &                                         wl_size_y(ilevel), &
  &                                         wl_size_x(ilevel))
    integer(kind=wip) :: iall
    integer(kind=wip) :: icx0, icy0, icz0
    integer(kind=wip) :: icx1, icy1, icz1
    integer(kind=wip) :: iam
    integer(kind=wip) :: irx, iry, irz
    integer(kind=wip) :: nomp,iiam
    integer(kind=wip) :: nsczdiv,nscydiv,nscxdiv
! add padding(+5,,+2,,,).
    real(kind=wvap) :: wl_omp(lm_length+5, 2, &
                                    wl_size_z(ilevel)+2, &
                                    wl_size_y(ilevel)  , &
                                    wl_size_x(ilevel)  , 0:nomp-1)

    real(kind=wvp) :: wwm_north(lm_length,2,8), wwl_north(lm_length,2,8)
    integer(kind=wip) :: icp
    integer(kind=wip) :: icx0_1,icy0_1,icz0_1,icx1_1,icy1_1,icz1_1
    integer(kind=wip) :: icx0_2,icy0_2,icz0_2,icx1_2,icy1_2,icz1_2
    integer(kind=wip) :: icx0_3,icy0_3,icz0_3,icx1_3,icy1_3,icz1_3
    integer(kind=wip) :: icx0_4,icy0_4,icz0_4,icx1_4,icy1_4,icz1_4
    integer(kind=wip) :: icx0_5,icy0_5,icz0_5,icx1_5,icy1_5,icz1_5
    integer(kind=wip) :: icx0_6,icy0_6,icz0_6,icx1_6,icy1_6,icz1_6
    integer(kind=wip) :: icx0_7,icy0_7,icz0_7,icx1_7,icy1_7,icz1_7
    integer(kind=wip) :: icx0_8,icy0_8,icz0_8,icx1_8,icy1_8,icz1_8
    real(kind=wvp) :: shrot1_r, shrot2_r, shrot1_i, shrot2_i
    real(kind=wvp) :: shinvrot1_r, shinvrot2_r, shinvrot1_i, shinvrot2_i
    real(kind=wvp) :: wm_1r, wm_2r, wm_3r, wm_4r, wm_5r, wm_6r, wm_7r, wm_8r
    real(kind=wvp) :: wm_1i, wm_2i, wm_3i, wm_4i, wm_5i, wm_6i, wm_7i, wm_8i
    real(kind=wvp) :: wlnorth_1r,wlnorth_2r,wlnorth_3r,wlnorth_4r,wlnorth_5r,wlnorth_6r,wlnorth_7r,wlnorth_8r
    real(kind=wvp) :: wlnorth_1i,wlnorth_2i,wlnorth_3i,wlnorth_4i,wlnorth_5i,wlnorth_6i,wlnorth_7i,wlnorth_8i
    real(kind=wvp) :: wm_north(lm_length,2), wl_north(lm_length,2)
    real(kind=wvp) :: tmp_shml, tmp_ce_r, tmp_ce_i
    integer(kind=wip) :: j,k,n,p,q, m1,m2,m2_b, ind_1dim,ind_1dim_b

!#ifndef _OPENMP
    iam = 0
!#endif

    wl_omp = Cdzero

    nsczdiv = wl_size_z(ilevel);  nscydiv = wl_size_y(ilevel);  nscxdiv = wl_size_x(ilevel)

!$omp parallel default(none) &
!$omp& private(iam) &
!$omp& private(icx0,icy0,icz0,icx1,icy1,icz1) &
!$omp& private(irx,iry,irz) &
!$omp& private(wwm_north,wwl_north) &
!$omp& private(icp) &
!$omp& private(icx0_1,icy0_1,icz0_1,icx1_1,icy1_1,icz1_1) &
!$omp& private(icx0_2,icy0_2,icz0_2,icx1_2,icy1_2,icz1_2) &
!$omp& private(icx0_3,icy0_3,icz0_3,icx1_3,icy1_3,icz1_3) &
!$omp& private(icx0_4,icy0_4,icz0_4,icx1_4,icy1_4,icz1_4) &
!$omp& private(icx0_5,icy0_5,icz0_5,icx1_5,icy1_5,icz1_5) &
!$omp& private(icx0_6,icy0_6,icz0_6,icx1_6,icy1_6,icz1_6) &
!$omp& private(icx0_7,icy0_7,icz0_7,icx1_7,icy1_7,icz1_7) &
!$omp& private(icx0_8,icy0_8,icz0_8,icx1_8,icy1_8,icz1_8) &
!$omp& private(shrot1_r, shrot2_r, shrot1_i, shrot2_i) &
!$omp& private(shinvrot1_r, shinvrot2_r, shinvrot1_i, shinvrot2_i) &
!$omp& private(wm_1r, wm_2r, wm_3r, wm_4r, wm_5r, wm_6r, wm_7r, wm_8r) &
!$omp& private(wm_1i, wm_2i, wm_3i, wm_4i, wm_5i, wm_6i, wm_7i, wm_8i) &
!$omp& private(wlnorth_1r,wlnorth_2r,wlnorth_3r,wlnorth_4r,wlnorth_5r,wlnorth_6r,wlnorth_7r,wlnorth_8r) &
!$omp& private(wlnorth_1i,wlnorth_2i,wlnorth_3i,wlnorth_4i,wlnorth_5i,wlnorth_6i,wlnorth_7i,wlnorth_8i) &
!$omp& private(wm_north,wl_north) &
!$omp& private(j,k,n,p,q, m1,m2,m2_b, ind_1dim,ind_1dim_b) &
!$omp& private(tmp_shml, tmp_ce_r, tmp_ce_i) &
!$omp& shared(shml, sh_rot1, sh_rot2, sh_inv_rot1, sh_inv_rot2) &
!$omp& shared(nmax) &
!debug
!!!$omp& reduction(+:numjc) &
!debug
!$omp& shared(ilevel) &
!$omp& shared(wl_omp,wm) &
!$omp& shared(cp8_count,cp6_count,cp4_count,cp3_count,cp2_count,cp1_count) &
!$omp& shared(cp8_list0,cp6_list0,cp4_list0,cp3_list0,cp2_list0,cp1_list0) &
!$omp& shared(cp8_list1,cp6_list1,cp4_list1,cp3_list1,cp2_list1,cp1_list1) &
!$omp& shared(cp8_irxyz,cp6_irxyz,cp4_irxyz,cp3_irxyz,cp2_irxyz,cp1_irxyz) &
!$omp& shared(nscxdiv,nscydiv,nsczdiv) &
!$omp& shared(lm_length) &
!$omp& shared(nchunk)
!$  iam = omp_get_thread_num()
!!$omp do schedule(static,nchunk(ilevel))

    if(cp8_count(ilevel) > 0) then    !! /*---- pairs multiple of 8 for each relative cell address ----*/

!$omp do
    do icp = 1, cp8_count(ilevel), 8
      irx = cp8_irxyz(1,icp,ilevel)  ! x-dir relative cell address of transformation matrix.
      iry = cp8_irxyz(2,icp,ilevel)  ! y-dir relative cell address of transformation matrix.
      irz = cp8_irxyz(3,icp,ilevel)  ! z-dir relative cell address of transformation matrix.

    !*** multipole to local translation
      icx0_1 = cp8_list0(1,icp  ,ilevel) ; icy0_1 = cp8_list0(2,icp  ,ilevel) ; icz0_1 = cp8_list0(3,icp  ,ilevel)
      icx1_1 = cp8_list1(1,icp  ,ilevel) ; icy1_1 = cp8_list1(2,icp  ,ilevel) ; icz1_1 = cp8_list1(3,icp  ,ilevel)
      icx0_2 = cp8_list0(1,icp+1,ilevel) ; icy0_2 = cp8_list0(2,icp+1,ilevel) ; icz0_2 = cp8_list0(3,icp+1,ilevel)
      icx1_2 = cp8_list1(1,icp+1,ilevel) ; icy1_2 = cp8_list1(2,icp+1,ilevel) ; icz1_2 = cp8_list1(3,icp+1,ilevel)
      icx0_3 = cp8_list0(1,icp+2,ilevel) ; icy0_3 = cp8_list0(2,icp+2,ilevel) ; icz0_3 = cp8_list0(3,icp+2,ilevel)
      icx1_3 = cp8_list1(1,icp+2,ilevel) ; icy1_3 = cp8_list1(2,icp+2,ilevel) ; icz1_3 = cp8_list1(3,icp+2,ilevel)
      icx0_4 = cp8_list0(1,icp+3,ilevel) ; icy0_4 = cp8_list0(2,icp+3,ilevel) ; icz0_4 = cp8_list0(3,icp+3,ilevel)
      icx1_4 = cp8_list1(1,icp+3,ilevel) ; icy1_4 = cp8_list1(2,icp+3,ilevel) ; icz1_4 = cp8_list1(3,icp+3,ilevel)
      icx0_5 = cp8_list0(1,icp+4,ilevel) ; icy0_5 = cp8_list0(2,icp+4,ilevel) ; icz0_5 = cp8_list0(3,icp+4,ilevel)
      icx1_5 = cp8_list1(1,icp+4,ilevel) ; icy1_5 = cp8_list1(2,icp+4,ilevel) ; icz1_5 = cp8_list1(3,icp+4,ilevel)
      icx0_6 = cp8_list0(1,icp+5,ilevel) ; icy0_6 = cp8_list0(2,icp+5,ilevel) ; icz0_6 = cp8_list0(3,icp+5,ilevel)
      icx1_6 = cp8_list1(1,icp+5,ilevel) ; icy1_6 = cp8_list1(2,icp+5,ilevel) ; icz1_6 = cp8_list1(3,icp+5,ilevel)
      icx0_7 = cp8_list0(1,icp+6,ilevel) ; icy0_7 = cp8_list0(2,icp+6,ilevel) ; icz0_7 = cp8_list0(3,icp+6,ilevel)
      icx1_7 = cp8_list1(1,icp+6,ilevel) ; icy1_7 = cp8_list1(2,icp+6,ilevel) ; icz1_7 = cp8_list1(3,icp+6,ilevel)
      icx0_8 = cp8_list0(1,icp+7,ilevel) ; icy0_8 = cp8_list0(2,icp+7,ilevel) ; icz0_8 = cp8_list0(3,icp+7,ilevel)
      icx1_8 = cp8_list1(1,icp+7,ilevel) ; icy1_8 = cp8_list1(2,icp+7,ilevel) ; icz1_8 = cp8_list1(3,icp+7,ilevel)

      include 'tuned_M2L/fmm_far_m2l_rotation_mpair8_exsimd-reduct.h90'
    end do  ! cp8_count.
!$omp end do

    endif    ! cp8_count>0.

    if(cp6_count(ilevel) > 0) then    !! /*---- pairs multiple of 6 for each relative cell address ----*/

!$omp do
    do icp = 1, cp6_count(ilevel), 6
      irx = cp6_irxyz(1,icp,ilevel)  ! x-dir relative cell address of transformation matrix.
      iry = cp6_irxyz(2,icp,ilevel)  ! y-dir relative cell address of transformation matrix.
      irz = cp6_irxyz(3,icp,ilevel)  ! z-dir relative cell address of transformation matrix.

    !*** multipole to local translation
      icx0_1 = cp6_list0(1,icp  ,ilevel) ; icy0_1 = cp6_list0(2,icp  ,ilevel) ; icz0_1 = cp6_list0(3,icp  ,ilevel)
      icx1_1 = cp6_list1(1,icp  ,ilevel) ; icy1_1 = cp6_list1(2,icp  ,ilevel) ; icz1_1 = cp6_list1(3,icp  ,ilevel)
      icx0_2 = cp6_list0(1,icp+1,ilevel) ; icy0_2 = cp6_list0(2,icp+1,ilevel) ; icz0_2 = cp6_list0(3,icp+1,ilevel)
      icx1_2 = cp6_list1(1,icp+1,ilevel) ; icy1_2 = cp6_list1(2,icp+1,ilevel) ; icz1_2 = cp6_list1(3,icp+1,ilevel)
      icx0_3 = cp6_list0(1,icp+2,ilevel) ; icy0_3 = cp6_list0(2,icp+2,ilevel) ; icz0_3 = cp6_list0(3,icp+2,ilevel)
      icx1_3 = cp6_list1(1,icp+2,ilevel) ; icy1_3 = cp6_list1(2,icp+2,ilevel) ; icz1_3 = cp6_list1(3,icp+2,ilevel)
      icx0_4 = cp6_list0(1,icp+3,ilevel) ; icy0_4 = cp6_list0(2,icp+3,ilevel) ; icz0_4 = cp6_list0(3,icp+3,ilevel)
      icx1_4 = cp6_list1(1,icp+3,ilevel) ; icy1_4 = cp6_list1(2,icp+3,ilevel) ; icz1_4 = cp6_list1(3,icp+3,ilevel)
      icx0_5 = cp6_list0(1,icp+4,ilevel) ; icy0_5 = cp6_list0(2,icp+4,ilevel) ; icz0_5 = cp6_list0(3,icp+4,ilevel)
      icx1_5 = cp6_list1(1,icp+4,ilevel) ; icy1_5 = cp6_list1(2,icp+4,ilevel) ; icz1_5 = cp6_list1(3,icp+4,ilevel)
      icx0_6 = cp6_list0(1,icp+5,ilevel) ; icy0_6 = cp6_list0(2,icp+5,ilevel) ; icz0_6 = cp6_list0(3,icp+5,ilevel)
      icx1_6 = cp6_list1(1,icp+5,ilevel) ; icy1_6 = cp6_list1(2,icp+5,ilevel) ; icz1_6 = cp6_list1(3,icp+5,ilevel)

      include 'tuned_M2L/fmm_far_m2l_rotation_mpair6_exsimd-reduct.h90'
    end do  ! cp6_count.
!$omp end do

    endif    ! cp6_count>0.

    if(cp4_count(ilevel) > 0) then    !! /*---- pairs multiple of 4 for each relative cell address ----*/

!$omp do
    do icp = 1, cp4_count(ilevel), 4
      irx = cp4_irxyz(1,icp,ilevel)  ! x-dir relative cell address of transformation matrix.
      iry = cp4_irxyz(2,icp,ilevel)  ! y-dir relative cell address of transformation matrix.
      irz = cp4_irxyz(3,icp,ilevel)  ! z-dir relative cell address of transformation matrix.

    !*** multipole to local translation
      icx0_1 = cp4_list0(1,icp  ,ilevel) ; icy0_1 = cp4_list0(2,icp  ,ilevel) ; icz0_1 = cp4_list0(3,icp  ,ilevel)
      icx1_1 = cp4_list1(1,icp  ,ilevel) ; icy1_1 = cp4_list1(2,icp  ,ilevel) ; icz1_1 = cp4_list1(3,icp  ,ilevel)
      icx0_2 = cp4_list0(1,icp+1,ilevel) ; icy0_2 = cp4_list0(2,icp+1,ilevel) ; icz0_2 = cp4_list0(3,icp+1,ilevel)
      icx1_2 = cp4_list1(1,icp+1,ilevel) ; icy1_2 = cp4_list1(2,icp+1,ilevel) ; icz1_2 = cp4_list1(3,icp+1,ilevel)
      icx0_3 = cp4_list0(1,icp+2,ilevel) ; icy0_3 = cp4_list0(2,icp+2,ilevel) ; icz0_3 = cp4_list0(3,icp+2,ilevel)
      icx1_3 = cp4_list1(1,icp+2,ilevel) ; icy1_3 = cp4_list1(2,icp+2,ilevel) ; icz1_3 = cp4_list1(3,icp+2,ilevel)
      icx0_4 = cp4_list0(1,icp+3,ilevel) ; icy0_4 = cp4_list0(2,icp+3,ilevel) ; icz0_4 = cp4_list0(3,icp+3,ilevel)
      icx1_4 = cp4_list1(1,icp+3,ilevel) ; icy1_4 = cp4_list1(2,icp+3,ilevel) ; icz1_4 = cp4_list1(3,icp+3,ilevel)

      include 'tuned_M2L/fmm_far_m2l_rotation_mpair4_exsimd-reduct.h90'
    end do  ! cp4_count.
!$omp end do

    endif    ! cp4_count>0.

    if(cp3_count(ilevel) > 0) then    !! /*---- pairs multiple of 3 for each relative cell address ----*/

!$omp do
    do icp = 1, cp3_count(ilevel), 3
      irx = cp3_irxyz(1,icp,ilevel)  ! x-dir relative cell address of transformation matrix.
      iry = cp3_irxyz(2,icp,ilevel)  ! y-dir relative cell address of transformation matrix.
      irz = cp3_irxyz(3,icp,ilevel)  ! z-dir relative cell address of transformation matrix.

    !*** multipole to local translation
      icx0_1 = cp3_list0(1,icp  ,ilevel) ; icy0_1 = cp3_list0(2,icp  ,ilevel) ; icz0_1 = cp3_list0(3,icp  ,ilevel)
      icx1_1 = cp3_list1(1,icp  ,ilevel) ; icy1_1 = cp3_list1(2,icp  ,ilevel) ; icz1_1 = cp3_list1(3,icp  ,ilevel)
      icx0_2 = cp3_list0(1,icp+1,ilevel) ; icy0_2 = cp3_list0(2,icp+1,ilevel) ; icz0_2 = cp3_list0(3,icp+1,ilevel)
      icx1_2 = cp3_list1(1,icp+1,ilevel) ; icy1_2 = cp3_list1(2,icp+1,ilevel) ; icz1_2 = cp3_list1(3,icp+1,ilevel)
      icx0_3 = cp3_list0(1,icp+2,ilevel) ; icy0_3 = cp3_list0(2,icp+2,ilevel) ; icz0_3 = cp3_list0(3,icp+2,ilevel)
      icx1_3 = cp3_list1(1,icp+2,ilevel) ; icy1_3 = cp3_list1(2,icp+2,ilevel) ; icz1_3 = cp3_list1(3,icp+2,ilevel)

      include 'tuned_M2L/fmm_far_m2l_rotation_mpair3_exsimd-reduct.h90'
    end do  ! cp3_count.
!$omp end do

    endif    ! cp3_count>0.

    if(cp2_count(ilevel) > 0) then    !! /*---- pairs multiple of 2 for each relative cell address ----*/

!$omp do
    do icp = 1, cp2_count(ilevel), 2
      irx = cp2_irxyz(1,icp,ilevel)  ! x-dir relative cell address of transformation matrix.
      iry = cp2_irxyz(2,icp,ilevel)  ! y-dir relative cell address of transformation matrix.
      irz = cp2_irxyz(3,icp,ilevel)  ! z-dir relative cell address of transformation matrix.

    !*** multipole to local translation
      icx0_1 = cp2_list0(1,icp  ,ilevel) ; icy0_1 = cp2_list0(2,icp  ,ilevel) ; icz0_1 = cp2_list0(3,icp  ,ilevel)
      icx1_1 = cp2_list1(1,icp  ,ilevel) ; icy1_1 = cp2_list1(2,icp  ,ilevel) ; icz1_1 = cp2_list1(3,icp  ,ilevel)
      icx0_2 = cp2_list0(1,icp+1,ilevel) ; icy0_2 = cp2_list0(2,icp+1,ilevel) ; icz0_2 = cp2_list0(3,icp+1,ilevel)
      icx1_2 = cp2_list1(1,icp+1,ilevel) ; icy1_2 = cp2_list1(2,icp+1,ilevel) ; icz1_2 = cp2_list1(3,icp+1,ilevel)

      include 'tuned_M2L/fmm_far_m2l_rotation_mpair2_exsimd-reduct.h90'
    end do  ! cp2_count.
!$omp end do

    endif    ! cp2_count>0.

    if(cp1_count(ilevel) > 0) then    !! /*---- pairs multiple of 1 for each relative cell address ----*/

!$omp do
    do icp = 1, cp1_count(ilevel)
      irx = cp1_irxyz(1,icp,ilevel)  ! x-dir relative cell address of transformation matrix.
      iry = cp1_irxyz(2,icp,ilevel)  ! y-dir relative cell address of transformation matrix.
      irz = cp1_irxyz(3,icp,ilevel)  ! z-dir relative cell address of transformation matrix.

    !*** multipole to local translation
      icx0 = cp1_list0(1,icp  ,ilevel) ; icy0 = cp1_list0(2,icp  ,ilevel) ; icz0 = cp1_list0(3,icp  ,ilevel)
      icx1 = cp1_list1(1,icp  ,ilevel) ; icy1 = cp1_list1(2,icp  ,ilevel) ; icz1 = cp1_list1(3,icp  ,ilevel)

      include 'tuned_M2L/fmm_far_m2l_rotation_exsimd-reduct.h90'
    end do  ! cp1_count.
!$omp end do

    endif    ! cp1_count>0.

!$omp end parallel


!$omp parallel default(none) &
!$omp& private(iall,icz0,icy0,icx0) &
!$omp& private(m1) &
!$omp& private(iiam) &
!$omp& shared(lm_length) &
!$omp& shared(nscxdiv,nscydiv,nsczdiv) &
!$omp& shared(wl,wl_omp) &
!$omp& shared(nomp)
!$omp do
    DO iall=0, nscxdiv * nscydiv * nsczdiv - 1
       icz0 = mod(iall,nsczdiv) + 1
       icy0 = iall/nsczdiv
       icx0 = icy0/nscydiv + 1
       icy0 = mod(icy0,nscydiv) + 1
       do iiam=0,nomp-1
!c !$omp do
       do m1=1,lm_length
!c          do iiam=0,nomp-1
             wl(m1,1,icz0,icy0,icx0) = wl(m1,1,icz0,icy0,icx0) + wl_omp(m1,1,icz0,icy0,icx0,iiam)
             wl(m1,2,icz0,icy0,icx0) = wl(m1,2,icz0,icy0,icx0) + wl_omp(m1,2,icz0,icy0,icx0,iiam)
!c          enddo
       enddo
       enddo
!c !$omp end do nowait
    ENDDO
!$omp end do
!$omp end parallel

#ifdef DEBUG_MTDFMM
    DO iall=0, nscxdiv * nscydiv * nsczdiv - 1
       icz0 = mod(iall,nsczdiv) + 1
       icy0 = iall/nsczdiv
       icx0 = icy0/nscydiv + 1
       icy0 = mod(icy0,nscydiv) + 1
       write(myrank+ilevel*1000+40000,*) wl(:,1,icz0,icy0,icx0)
       write(myrank+ilevel*1000+40000,*) wl(:,2,icz0,icy0,icx0)
    ENDDO
    call flush(myrank+ilevel*1000+40000)
#endif

  end subroutine M2L_upper_level_mpairs
!---------------------------------------------------------------------
#else   /* Not MOMENT */
!---------------------------------------------------------------------
!>
!! \brief  Subroutine to execute M2L (with global arrays, upper level)
!! \author Noriyuki Yoshii, Yoshimichi Andoh, Shin-chi Ichikawa
!<
  subroutine M2L_upper_level_nmax4(ilevel, wm, wl, nomp)
!---------------------------------------------------------------------
! symmetry and rotation code. fully unrolled for nmax=4.
!
    implicit none
    integer(kind=wip), intent(in)  :: ilevel
    real(kind=wmp), intent(in)     :: wm(lm_length, 2, &
  &                                         wm_size_z(ilevel), &
  &                                         wm_size_y(ilevel), &
  &                                         wm_size_x(ilevel))
    real(kind=wvap), intent(inout) :: wl(lm_length, 2, &
  &                                         wl_size_z(ilevel), &
  &                                         wl_size_y(ilevel), &
  &                                         wl_size_x(ilevel))
    integer(kind=wip) :: iall
    integer(kind=wip) :: icx0, icy0, icz0
    integer(kind=wip) :: icx1, icy1, icz1
    integer(kind=wip) :: irx, iry, irz
    integer(kind=wip) :: nomp,iam,iiam
    integer(kind=wip) :: nsczdiv,nscydiv,nscxdiv
! add padding(+5,,+2,,,).
    real(kind=wvap) :: wl_omp(lm_length+5, 2, &
                                    wl_size_z(ilevel)+2, &
                                    wl_size_y(ilevel)  , &
                                    wl_size_x(ilevel)  , 0:nomp-1)
    integer(kind=wip) :: nscellx, nscelly, nscellz
    integer(kind=wip) :: icxg0, icyg0, iczg0
    integer(kind=wip) :: load
    real(kind=wvp) :: wm_north(lm_length,2), wl_north(lm_length,2)
    real(kind=wvp) :: tmp_ce_r, tmp_ce_i
    integer(kind=wip) :: m1

!#ifndef _OPENMP
    iam = 0
!#endif

    wl_omp = Cdzero

    nsczdiv = wl_size_z(ilevel);  nscydiv = wl_size_y(ilevel);  nscxdiv = wl_size_x(ilevel)
    nscellz = nscllz(ilevel);  nscelly = nsclly(ilevel);  nscellx = nscllx(ilevel)
    
    ! global address of starting cell.
    icxg0 = get_ixmin_in_fmm(ilevel)
    icyg0 = get_iymin_in_fmm(ilevel)
    iczg0 = get_izmin_in_fmm(ilevel)

!$omp parallel default(none) &
!$omp& private(iam) &
!$omp& private(load) &
!$omp& private(icx0,icy0,icz0,icx1,icy1,icz1) &
!$omp& private(irx,iry,irz) &
!$omp& private(wm_north,wl_north) &
!$omp& private(tmp_ce_r, tmp_ce_i) &
!$omp& shared(shml, sh_rot1, sh_rot2, sh_inv_rot1, sh_inv_rot2) &
!$omp& shared(ilevel) &
!$omp& shared(wl_omp,wm) &
!$omp& shared(nload,lddir) &
!$omp& shared(nscellx,nscelly,nscellz) &
!$omp& shared(icxg0,icyg0,iczg0) &
!$omp& shared(nscxdiv,nscydiv,nsczdiv) &
!$omp& shared(nlx,nly,nlz) &
!$omp& shared(lm_length) &
!$omp& shared(nchunk)
!$  iam = omp_get_thread_num()
!$omp do schedule(static,1)
    DO load = 1, nload(ilevel)
       irx = lddir(1,load,ilevel)
       iry = lddir(2,load,ilevel)
       irz = lddir(3,load,ilevel)

       do icx0 = max(1, -2*nscellx-icxg0+2-abs(irx)), min(nscxdiv, 3*nscellx-icxg0+1+abs(irx))

          if (.not. is_valid_cell(icxg0, icx0, irx, nlx(ilevel)) ) cycle
          icx1 = modulo(icxg0-1 + icx0-1 + irx, nscellx) + 1

          do icy0 = max(1, -2*nscelly-icyg0+2-abs(iry)), min(nscydiv, 3*nscelly-icyg0+1+abs(iry))

             if (.not. is_valid_cell(icyg0, icy0, iry, nly(ilevel)) ) cycle
             icy1 = modulo(icyg0-1 + icy0-1 + iry, nscelly) + 1

             do icz0 = max(1, -2*nscellz-iczg0+2-abs(irz)), min(nsczdiv, 3*nscellz-iczg0+1+abs(irz))

                if (.not. is_valid_cell(iczg0, icz0, irz, nlz(ilevel)) ) cycle
                icz1 = modulo(iczg0-1 + icz0-1 + irz, nscellz) + 1

                !*** multipole to local translation
                include 'tuned_M2L/fmm_far_m2l_rotation_nmax4_unrolled_fj.h90'

             end do ! icz0
          end do ! icy0
       end do ! icx0

    enddo ! load
!$omp end do nowait
!$omp end parallel

!$omp parallel default(none) &
!$omp& private(iall,icz0,icy0,icx0) &
!$omp& private(m1) &
!$omp& private(iiam) &
!$omp& shared(lm_length) &
!$omp& shared(nscxdiv,nscydiv,nsczdiv) &
!$omp& shared(wl,wl_omp) &
!$omp& shared(nomp)
!$omp do
    DO iall=0, nscxdiv * nscydiv * nsczdiv - 1
       icz0 = mod(iall,nsczdiv) + 1
       icy0 = iall/nsczdiv
       icx0 = icy0/nscydiv + 1
       icy0 = mod(icy0,nscydiv) + 1
       do iiam=0,nomp-1
       do m1=1,lm_length
             wl(m1,1,icz0,icy0,icx0) = wl(m1,1,icz0,icy0,icx0) + wl_omp(m1,1,icz0,icy0,icx0,iiam)
             wl(m1,2,icz0,icy0,icx0) = wl(m1,2,icz0,icy0,icx0) + wl_omp(m1,2,icz0,icy0,icx0,iiam)
       enddo
       enddo
    ENDDO
!$omp end do
!$omp end parallel

#ifdef DEBUG_MTDFMM
    DO iall=0, nscxdiv * nscydiv * nsczdiv - 1
       icz0 = mod(iall,nsczdiv) + 1
       icy0 = iall/nsczdiv
       icx0 = icy0/nscydiv + 1
       icy0 = mod(icy0,nscydiv) + 1
       write(myrank+ilevel*1000+40000,*) wl(:,1,icz0,icy0,icx0)
       write(myrank+ilevel*1000+40000,*) wl(:,2,icz0,icy0,icx0)
    ENDDO
    call flush(myrank+ilevel*1000+40000)
#endif

  end subroutine M2L_upper_level_nmax4
!---------------------------------------------------------------------
!>
!! \brief  Subroutine to execute M2L (with global arrays, upper level)
!! \author Noriyuki Yoshii, Yoshimichi Andoh, Shin-chi Ichikawa
!<
  subroutine M2L_upper_level_general(ilevel, wm, wl, nomp)
!---------------------------------------------------------------------
! symmetry and rotation code. for general value of nmax.
!
    implicit none
    integer(kind=wip), intent(in)  :: ilevel
    real(kind=wmp), intent(in)     :: wm(lm_length, 2, &
  &                                         wm_size_z(ilevel), &
  &                                         wm_size_y(ilevel), &
  &                                         wm_size_x(ilevel))
    real(kind=wvap), intent(inout) :: wl(lm_length, 2, &
  &                                         wl_size_z(ilevel), &
  &                                         wl_size_y(ilevel), &
  &                                         wl_size_x(ilevel))
    integer(kind=wip) :: iall
    integer(kind=wip) :: icx0, icy0, icz0
    integer(kind=wip) :: icx1, icy1, icz1
    integer(kind=wip) :: irx, iry, irz
    integer(kind=wip) :: nomp,iam,iiam
    integer(kind=wip) :: nsczdiv,nscydiv,nscxdiv
! add padding(+5,,+2,,,).
    real(kind=wvap) :: wl_omp(lm_length+5, 2, &
                                    wl_size_z(ilevel)+2, &
                                    wl_size_y(ilevel)  , &
                                    wl_size_x(ilevel)  , 0:nomp-1)
    integer(kind=wip) :: nscellx, nscelly, nscellz
    integer(kind=wip) :: icxg0, icyg0, iczg0
    integer(kind=wip) :: load
    real(kind=wvp) :: wm_north(lm_length,2), wl_north(lm_length,2)
    real(kind=wvp) :: tmp_shml, tmp_ce_r, tmp_ce_i
    integer(kind=wip) :: j,k,n,p,q, m1,m2,m2_b, ind_1dim,ind_1dim_b

    iam = 0

    wl_omp = Cdzero

    nsczdiv = wl_size_z(ilevel);  nscydiv = wl_size_y(ilevel);  nscxdiv = wl_size_x(ilevel)
    nscellz = nscllz(ilevel);  nscelly = nsclly(ilevel);  nscellx = nscllx(ilevel)
    
    ! global address of starting cell.
    icxg0 = get_ixmin_in_fmm(ilevel)
    icyg0 = get_iymin_in_fmm(ilevel)
    iczg0 = get_izmin_in_fmm(ilevel)

!$omp parallel default(none) &
!$omp& private(iam) &
!$omp& private(load) &
!$omp& private(icx0,icy0,icz0,icx1,icy1,icz1) &
!$omp& private(irx,iry,irz) &
!$omp& private(wm_north,wl_north) &
!$omp& private(j,k,n,p,q, m1,m2,m2_b, ind_1dim,ind_1dim_b) &
!$omp& private(tmp_shml, tmp_ce_r, tmp_ce_i) &
!$omp& shared(shml, sh_rot1, sh_rot2, sh_inv_rot1, sh_inv_rot2) &
!$omp& shared(nmax) &
!$omp& shared(ilevel) &
!$omp& shared(wl_omp,wm) &
!$omp& shared(nload,lddir) &
!$omp& shared(nscellx,nscelly,nscellz) &
!$omp& shared(icxg0,icyg0,iczg0) &
!$omp& shared(nscxdiv,nscydiv,nsczdiv) &
!$omp& shared(nlx,nly,nlz) &
!$omp& shared(lm_length) &
!$omp& shared(nchunk)
!$  iam = omp_get_thread_num()
!$omp do schedule(static,1)
    DO load = 1, nload(ilevel)
       irx = lddir(1,load,ilevel)
       iry = lddir(2,load,ilevel)
       irz = lddir(3,load,ilevel)

       do icx0 = max(1, -2*nscellx-icxg0+2-abs(irx)), min(nscxdiv, 3*nscellx-icxg0+1+abs(irx))

          if (.not. is_valid_cell(icxg0, icx0, irx, nlx(ilevel)) ) cycle
          icx1 = modulo(icxg0-1 + icx0-1 + irx, nscellx) + 1

          do icy0 = max(1, -2*nscelly-icyg0+2-abs(iry)), min(nscydiv, 3*nscelly-icyg0+1+abs(iry))

             if (.not. is_valid_cell(icyg0, icy0, iry, nly(ilevel)) ) cycle
             icy1 = modulo(icyg0-1 + icy0-1 + iry, nscelly) + 1

             do icz0 = max(1, -2*nscellz-iczg0+2-abs(irz)), min(nsczdiv, 3*nscellz-iczg0+1+abs(irz))

                if (.not. is_valid_cell(iczg0, icz0, irz, nlz(ilevel)) ) cycle
                icz1 = modulo(iczg0-1 + icz0-1 + irz, nscellz) + 1

                !*** multipole to local translation
                include 'tuned_M2L/fmm_far_m2l_rotation_exsimd-reduct.h90'

             end do ! icz0
          end do ! icy0
       end do ! icx0

    enddo ! load
!$omp end do nowait
!$omp end parallel

!$omp parallel default(none) &
!$omp& private(iall,icz0,icy0,icx0) &
!$omp& private(m1) &
!$omp& private(iiam) &
!$omp& shared(lm_length) &
!$omp& shared(nscxdiv,nscydiv,nsczdiv) &
!$omp& shared(wl,wl_omp) &
!$omp& shared(nomp)
!$omp do
    DO iall=0, nscxdiv * nscydiv * nsczdiv - 1
       icz0 = mod(iall,nsczdiv) + 1
       icy0 = iall/nsczdiv
       icx0 = icy0/nscydiv + 1
       icy0 = mod(icy0,nscydiv) + 1
       do iiam=0,nomp-1
       do m1=1,lm_length
             wl(m1,1,icz0,icy0,icx0) = wl(m1,1,icz0,icy0,icx0) + wl_omp(m1,1,icz0,icy0,icx0,iiam)
             wl(m1,2,icz0,icy0,icx0) = wl(m1,2,icz0,icy0,icx0) + wl_omp(m1,2,icz0,icy0,icx0,iiam)
       enddo
       enddo
    ENDDO
!$omp end do
!$omp end parallel

#ifdef DEBUG_MTDFMM
    DO iall=0, nscxdiv * nscydiv * nsczdiv - 1
       icz0 = mod(iall,nsczdiv) + 1
       icy0 = iall/nsczdiv
       icx0 = icy0/nscydiv + 1
       icy0 = mod(icy0,nscydiv) + 1
       write(myrank+ilevel*1000+40000,*) wl(:,1,icz0,icy0,icx0)
       write(myrank+ilevel*1000+40000,*) wl(:,2,icz0,icy0,icx0)
    ENDDO
    call flush(myrank+ilevel*1000+40000)
#endif

  end subroutine M2L_upper_level_general
#endif  /* Not MOMENT */
!---------------------------------------------------------------------

!---------------------------------------------------------------------
! >> subroutines called in init_fmm
!---------------------------------------------------------------------
!>
!! \brief  Subroutine to generate list of relative cells in M2L.
!! \author Tatsuya Sakashita
!<
  subroutine generate_relative_cell_list_for_M2L()
!---------------------------------------------------------------------
  ! relative cell address list generation for any number of threads
  !... A multi-level cubic shells ordering with cyclic assignment is  ...
  !... employed to balance computational load of M2L translation among...
  !... threads is employed. Cyclic assignment causes frequent         ...
  !... replacement of j-cell moment.                                  ...
    use domain, only : nlevel
    use mpi_tool
    implicit none
    integer(kind=wip) :: ncrefx, ncrefy, ncrefz
    integer(kind=wip) :: max_cell_dist
    integer(kind=wip) :: nrx, nry, nrz
    integer(kind=wip) :: nex, ney, nez
    integer(kind=wip) :: ithr
    integer(kind=wip) :: ith,nscan, iht
    integer(kind=wip) :: ist,ien
    integer(kind=wip) :: isign,mst
    integer(kind=wip) :: ilevel
    integer(kind=wip) :: irx, iry, irz
    integer(kind=wip) :: i
    integer(kind=wip) :: ncref
    integer(kind=wip) :: nshell,ishl
    integer(kind=wip) :: nomp
    integer(kind=wip) :: ixst( (6+2)*max_ncref ), ixen( (6+2)*max_ncref )
    integer(kind=wip) :: iyst( (6+2)*max_ncref ), iyen( (6+2)*max_ncref )
    integer(kind=wip) :: izst( (6+2)*max_ncref ), izen( (6+2)*max_ncref )

    nomp = 1
!$  nomp = omp_get_max_threads()
    
    allocate( lddir(3, (2*max_ncref+1)**3 - (2*ncdir+1)**3, 0:nlevel-1) )      ! relative cell address list
    allocate( nload(0:nlevel-1) )
    allocate( nchunk(0:nlevel-1) )

    ! relative cell address list in an order to balance computational load.
    ! calculate reference range, shell definition and the list.

    DO ilevel = 0, nlevel-1    ! ---  level -loop.  ---

       ncrefx = 5
       if(nlx(ilevel) == 3) ncrefx = 8
       ncrefy = 5
       if(nly(ilevel) == 3) ncrefy = 8
       ncrefz = 5
       if(nlz(ilevel) == 3) ncrefz = 8
       
       max_cell_dist = ncrefx
       if(max_cell_dist < ncrefy) max_cell_dist = ncrefy
       if(max_cell_dist < ncrefz) max_cell_dist = ncrefz
       ncref = ncrefx
       if(ncref > ncrefy) ncref = ncrefy
       if(ncref > ncrefz) ncref = ncrefz

       ! shell definition within minimum reference range among direction.
       nshell = 0
       do i = 1, ncref
          ! -- 6 squares + 2 vertices. --
          ! z-x plane (y:negative)
          nshell = nshell + 1
          izst(nshell) = -i + 1 ; izen(nshell) =  i
          iyst(nshell) = -i     ; iyen(nshell) = -i
          ixst(nshell) = -i     ; ixen(nshell) =  i - 1
          ! z-x plane (y:positive)
          nshell = nshell + 1
          izst(nshell) = -i     ; izen(nshell) =  i - 1
          iyst(nshell) =  i     ; iyen(nshell) =  i
          ixst(nshell) = -i + 1 ; ixen(nshell) =  i
          ! y-z plane (x:negative)
          nshell = nshell + 1
          izst(nshell) = -i     ; izen(nshell) =  i - 1
          iyst(nshell) = -i + 1 ; iyen(nshell) =  i
          ixst(nshell) = -i     ; ixen(nshell) = -i
          ! y-z plane (x:positive)
          nshell = nshell + 1
          izst(nshell) = -i + 1 ; izen(nshell) =  i
          iyst(nshell) = -i     ; iyen(nshell) =  i - 1
          ixst(nshell) =  i     ; ixen(nshell) =  i
          ! x-y plane (z:negative)
          nshell = nshell + 1
          izst(nshell) = -i     ; izen(nshell) = -i
          iyst(nshell) = -i     ; iyen(nshell) =  i - 1
          ixst(nshell) = -i + 1 ; ixen(nshell) =  i
          ! x-y plane (z:positive)
          nshell = nshell + 1
          izst(nshell) =  i     ; izen(nshell) =  i
          iyst(nshell) = -i + 1 ; iyen(nshell) =  i
          ixst(nshell) = -i     ; ixen(nshell) =  i - 1
          ! vertex (-i,-i,-i)
          nshell = nshell + 1
          izst(nshell) = -i     ; izen(nshell) = -i
          iyst(nshell) = -i     ; iyen(nshell) = -i
          ixst(nshell) = -i     ; ixen(nshell) = -i
          ! vertex (i,i,i)
          nshell = nshell + 1
          izst(nshell) =  i     ; izen(nshell) =  i
          iyst(nshell) =  i     ; iyen(nshell) =  i
          ixst(nshell) =  i     ; ixen(nshell) =  i
       end do

       ! shell definition beyond the minimum reference range.
       if(ncref < max_cell_dist) then
          do i = ncref + 1, max_cell_dist
             if(nlx(ilevel) == 3) then
                nrx = i; nex = 1
             else
                nrx = 5; nex = 0
             endif
             if(nly(ilevel) == 3) then
                nry = i; ney = 1
             else
                nry = 5; ney = 0
             endif
             if(nlz(ilevel) == 3) then
                nrz = i; nez = 1
             else
                nrz = 5; nez = 0
             endif
             ! -- 6 squares + 2 vertices. --
             if(nly(ilevel) == 3) then
                ! z-x plane (y:negative)
                nshell = nshell + 1
                izst(nshell) = -nrz + nez; izen(nshell) = nrz
                iyst(nshell) = -i        ; iyen(nshell) = -i
                ixst(nshell) = -nrx      ; ixen(nshell) = nrx - nex
                ! z-x plane (y:positive)
                nshell = nshell + 1
                izst(nshell) = -nrz      ; izen(nshell) = nrz - nez
                iyst(nshell) = i         ; iyen(nshell) = i
                ixst(nshell) = -nrx + nex; ixen(nshell) = nrx
             endif
             if(nlx(ilevel) == 3) then
                ! y-z plane (x:negative)
                nshell = nshell + 1
                izst(nshell) = -nrz      ; izen(nshell) = nrz - nez
                iyst(nshell) = -nry + ney; iyen(nshell) = nry
                ixst(nshell) = -i        ; ixen(nshell) = -i
                ! y-z plane (x:positive)
                nshell = nshell + 1
                izst(nshell) = -nrz + nez; izen(nshell) = nrz
                iyst(nshell) = -nry      ; iyen(nshell) = nry - ney
                ixst(nshell) = i         ; ixen(nshell) = i
             endif
             if(nlz(ilevel) == 3) then
                ! x-y plane (z:negative)
                nshell = nshell + 1
                izst(nshell) = -i        ; izen(nshell) = -i
                iyst(nshell) = -nry      ; iyen(nshell) = nry - ney
                ixst(nshell) = -nrx + nex; ixen(nshell) = nrx
                ! x-y plane (z:positive)
                nshell = nshell + 1
                izst(nshell) = i         ; izen(nshell) = i
                iyst(nshell) = -nry + ney; iyen(nshell) = nry
                ixst(nshell) = -nrx      ; ixen(nshell) = nrx - nex
             endif
             if(nlx(ilevel) == 3) then
                if(nly(ilevel) == 3) then
                   if(nlz(ilevel) == 3) then
                      ! vertex (-i,-i,-i)
                      nshell = nshell + 1
                      izst(nshell) = -i; izen(nshell) = -i
                      iyst(nshell) = -i; iyen(nshell) = -i
                      ixst(nshell) = -i; ixen(nshell) = -i
                      ! vertex (i,i,i)
                      nshell = nshell + 1
                      izst(nshell) =  i; izen(nshell) =  i
                      iyst(nshell) =  i; iyen(nshell) =  i
                      ixst(nshell) =  i; ixen(nshell) =  i
                   endif
                endif
             endif
          end do     ! i

       endif     ! ncref < max_cell_dist.

       ! Relative cell address list ordered on multi-layer cubic shell.

       nload(ilevel) = 0

       ist = 1
       ien = nshell
       mst = 1
       nscan = 0

       do ishl = ist, ien, mst
          do irx = ixst(ishl), ixen(ishl)
             do iry = iyst(ishl), iyen(ishl)
                do irz = izst(ishl), izen(ishl)
                   iht = 0
                   if( irx < -2 .OR. irx > 2 ) iht = 1
                   if( iry < -2 .OR. iry > 2 ) iht = 1
                   if( irz < -2 .OR. irz > 2 ) iht = 1
                   if( iht == 1 ) then
                      nscan = nscan + 1
                      nload(ilevel) = nload(ilevel) + 1
                      lddir(1,nload(ilevel),ilevel) = irx
                      lddir(2,nload(ilevel),ilevel) = iry
                      lddir(3,nload(ilevel),ilevel) = irz
                   endif
                end do  ! irz
             end do  ! iry
          end do  ! irx

       end do  ! ishl

    END DO  ! nl

  end subroutine generate_relative_cell_list_for_M2L

#ifdef MOMENT
!---------------------------------------------------------------------
!>
!! \brief  Subroutine to assort list of relative cells in M2L.
!! \author Shi-ichi Ichikawa
!<
  ! generation of relative cell address list assorted by possible     ...
  ! number of cell pairs for each relative cell address.              ...
  ! the purpose is to execute particular code unrolled by cell-pairs  ...
  ! with unrolled matrix dimensions in M2L, and to avoid frequent     ...
  ! instruction access miss.                                          ...
  !... 
  subroutine generate_assorted_relative_cell_list_for_M2L()
    use mpi_3d_grid, only : npx, npy, npz, ipx, ipy, ipz
    use domain, only : nlevel
    use mpi_tool
    implicit none
    integer(kind=wip) :: ilevel
    integer(kind=wip) :: irx, iry, irz

    integer(kind=wip) :: nomp
  ! integer(kind=wip), parameter :: max_ncref = 8   ! maximum referencing cell-range for M2L translation
  ! integer(kind=wip), parameter :: ncdir = 2   ! cell range of direct force calculation

    integer(kind=wip) :: nscxdiv, nscydiv, nsczdiv
    integer(kind=wip) :: nscellx, nscelly, nscellz
    integer(kind=wip) :: icxg0, icyg0, iczg0
    integer(kind=wip) :: icxgo, icygo, iczgo
    integer(kind=wip) :: icxg,  icyg,  iczg
    integer(kind=wip) :: icx0, icy0, icz0, icx1, icy1, icz1
    integer(kind=wip) :: icx1s, icy1s, icz1s
    integer(kind=wip) :: icx1e, icy1e, icz1e

    integer(kind=wip) :: mbd_z, mbd_y, mbd_x
    integer(kind=wip) :: mbd_z2, mbd_y2, mbd_x2
    integer(kind=wip) :: ngord, nrgord
    integer(kind=wip) :: load
    integer(kind=wip) :: ncp

    allocate( lddir8p(3,(2*max_ncref+1)**3-(2*ncdir+1)**3,0:nlevel-1) )   ! relative cell address list for 8 pairs.
    allocate( lddir6p(3,(2*max_ncref+1)**3-(2*ncdir+1)**3,0:nlevel-1) )   ! relative cell address list for 6 pairs.
    allocate( lddir4p(3,(2*max_ncref+1)**3-(2*ncdir+1)**3,0:nlevel-1) )   ! relative cell address list for 4 pairs.
    allocate( lddir3p(3,(2*max_ncref+1)**3-(2*ncdir+1)**3,0:nlevel-1) )   ! relative cell address list for 3 pairs.
    allocate( lddir2p(3,(2*max_ncref+1)**3-(2*ncdir+1)**3,0:nlevel-1) )   ! relative cell address list for 2 pairs.
    allocate( lddir1p(3,(2*max_ncref+1)**3-(2*ncdir+1)**3,0:nlevel-1) )   ! relative cell address list for 1 pair.
    allocate( nload8p(0:nlevel-1) )
    allocate( nload6p(0:nlevel-1) )
    allocate( nload4p(0:nlevel-1) )
    allocate( nload3p(0:nlevel-1) )
    allocate( nload2p(0:nlevel-1) )
    allocate( nload1p(0:nlevel-1) )

    nload8p = 0; nload6p = 0; nload4p = 0; nload3p = 0; nload2p = 0; nload1p = 0

    do ilevel = 0, nlevel-1

   !----------------------------------------------------------------
    if(ilevel <= lgflg) then
   !----------------------------------------------------------------
    nsczdiv = wl_size_z(ilevel);  nscydiv = wl_size_y(ilevel);  nscxdiv = wl_size_x(ilevel)
    nscellz = nscllz(ilevel);  nscelly = nsclly(ilevel);  nscellx = nscllx(ilevel)

    ! global address offset of a cell prior to starting cell.
    icxgo = (nscellx * ipx) / npx
    icygo = (nscelly * ipy) / npy
    iczgo = (nscellz * ipz) / npz
    mbd_x = mbd_x_lvl(ilevel) ; mbd_x2 = 2*mbd_x
    mbd_y = mbd_y_lvl(ilevel) ; mbd_y2 = 2*mbd_y
    mbd_z = mbd_z_lvl(ilevel) ; mbd_z2 = 2*mbd_z

    do load = 1, nload(ilevel)
      irx = lddir(1,load,ilevel)     ! x-dir relative cell address of transformation matrix.
      iry = lddir(2,load,ilevel)     ! y-dir relative cell address of transformation matrix.
      irz = lddir(3,load,ilevel)     ! z-dir relative cell address of transformation matrix.
      ncp = 0                        ! cell pair count.

      if(abs(irx).le.2 .and. abs(iry).le.2 .and. abs(irz).le.2) cycle  ! lddir does not include true case.

    do icx0 = 1-mbd_x, nscxdiv+mbd_x      ! target cell x-address.
      if(nscxdiv == 1      .and. 1 - mbd_x /= icx0 .and. icx0 /= nscxdiv + mbd_x) cycle
      if(nscxdiv <  mbd_x2 .and. 1         <= icx0 .and. icx0 <= nscxdiv        ) cycle
      icxg  = modulo(icxgo + icx0 - 1, nscellx) + 1  ! target cell global x-address.
      ngord  = mod(icxg - 1, mbd_x)                  ! sequence number within cell merge group minus 1.
      nrgord = mbd_x - ngord - 1                     ! reverse sequence number within cell merge group minus 1.
      icx1s = icx0 - mbd_x2 - ngord                  ! first cell of the cell range to be transformed.
      icx1e = icx0 + mbd_x2 + nrgord                 ! final cell of the cell range to be transformed.
      if(icx0 <= mbd_x          ) icx1s = icx0       ! plus-direction half-space transformation.
      if(icx0 >  nscxdiv - mbd_x) icx1e = icx0 - 1   ! minus-direction half-space transformation. inter-process exclusive sontrol.
      icx1 = icx0 + irx                              ! real transforming cell defined by relative cell address.
      if(icx1 < icx1s .or. icx1e < icx1) cycle

      do icy0 = 1-mbd_y, nscydiv+mbd_y    ! target cell y-address.
        if(nscydiv == 1     .and. 1 - mbd_y /= icy0 .and. icy0 /= nscydiv + mbd_y) cycle
        if(nscydiv < mbd_y2 .and. 1         <= icy0 .and. icy0 <= nscydiv        ) cycle
        icyg  = modulo(icygo + icy0 - 1, nscelly) + 1
        ngord  = mod(icyg - 1, mbd_y)
        nrgord = mbd_y - ngord - 1
        icy1s = icy0 - mbd_y2 - ngord
        icy1e = icy0 + mbd_y2 + nrgord
        if(icy0 <= mbd_y          ) icy1s = icy0
        if(icy0 >  nscydiv - mbd_y) icy1e = icy0 - 1
        icy1 = icy0 + iry
        if(icy1 < icy1s .or. icy1e < icy1) cycle

        do icz0 = 1-mbd_z, nsczdiv+mbd_z  ! target cell z-address.
          if(nsczdiv == 1     .and. 1 - mbd_z /= icz0 .and. icz0 /= nsczdiv + mbd_z) cycle
          if(nsczdiv < mbd_z2 .and. 1         <= icz0 .and. icz0 <= nsczdiv        ) cycle
          iczg  = modulo(iczgo + icz0 - 1, nscellz) + 1
          ngord  = mod(iczg - 1, mbd_z)
          nrgord = mbd_z - ngord - 1
          icz1s = icz0 - mbd_z2 - ngord
          icz1e = icz0 + mbd_z2 + nrgord
          if(icz0 <= mbd_z          ) icz1s = icz0
          if(icz0 >  nsczdiv - mbd_z) icz1e = icz0 - 1
          icz1 = icz0 + irz
          if(icz1 < icz1s .or. icz1e < icz1) cycle

          ncp = ncp + 1

        end do  ! icz0.
      end do  ! icy0.
    end do  ! icx0.

    if(ncp > 0) then
      if(mod(ncp,8) == 0) then
         nload8p(ilevel) = nload8p(ilevel) + 1
         lddir8p(1,nload8p(ilevel),ilevel) = irx
         lddir8p(2,nload8p(ilevel),ilevel) = iry
         lddir8p(3,nload8p(ilevel),ilevel) = irz
      elseif(mod(ncp,6) == 0) then
         nload6p(ilevel) = nload6p(ilevel) + 1
         lddir6p(1,nload6p(ilevel),ilevel) = irx
         lddir6p(2,nload6p(ilevel),ilevel) = iry
         lddir6p(3,nload6p(ilevel),ilevel) = irz
      elseif(mod(ncp,4) == 0) then
         nload4p(ilevel) = nload4p(ilevel) + 1
         lddir4p(1,nload4p(ilevel),ilevel) = irx
         lddir4p(2,nload4p(ilevel),ilevel) = iry
         lddir4p(3,nload4p(ilevel),ilevel) = irz
      elseif(mod(ncp,3) == 0) then
         nload3p(ilevel) = nload3p(ilevel) + 1
         lddir3p(1,nload3p(ilevel),ilevel) = irx
         lddir3p(2,nload3p(ilevel),ilevel) = iry
         lddir3p(3,nload3p(ilevel),ilevel) = irz
      elseif(mod(ncp,2) == 0) then
         nload2p(ilevel) = nload2p(ilevel) + 1
         lddir2p(1,nload2p(ilevel),ilevel) = irx
         lddir2p(2,nload2p(ilevel),ilevel) = iry
         lddir2p(3,nload2p(ilevel),ilevel) = irz
      elseif(ncp == 1) then
         nload1p(ilevel) = nload1p(ilevel) + 1
         lddir1p(1,nload1p(ilevel),ilevel) = irx
         lddir1p(2,nload1p(ilevel),ilevel) = iry
         lddir1p(3,nload1p(ilevel),ilevel) = irz
      else
         ! error.
         write(0,*) "ERROR.1: In generate_assorted_relative_cell_list: ", &
              &     "cell pairs are not the multiple of 8, nor 6, nor 4, ", &
              &     "nor 3, nor 2 for each relative cell address."
         write(0,*) "ERROR.2: Level; ",ilevel," lgflg: ",lgflg, &
              &     "  Relative cell address irx,iry,irz; ",irx,iry,irz
         call modylas_abort
      endif
      !write(*,*) "$$$ Number of pairs for relative cell address (irx,iry,irz)=(",irx,",",iry,",",irz,") : ",ncp
    else
         ! no pair case.
    endif

    end do  ! load.

!    write(*,*) "+++ Level: ",ilevel," Number of rel.addr for pairs: nload8p; ",nload8p(ilevel), &
!      &        ", nload6p; ",nload6p(ilevel),", nload4p; ",nload4p(ilevel),", nload3p; ",nload3p(ilevel), &
!      &        ", nload2p; ",nload2p(ilevel),", nload1p; ",nload1p(ilevel)

   !----------------------------------------------------------------
    else    ! ilevel>lgflg.
   !----------------------------------------------------------------
    nsczdiv = wl_size_z(ilevel);  nscydiv = wl_size_y(ilevel);  nscxdiv = wl_size_x(ilevel)
    nscellz = nscllz(ilevel);  nscelly = nsclly(ilevel);  nscellx = nscllx(ilevel)
    
    ! global address of starting cell.
    icxg0 = get_ixmin_in_fmm(ilevel)
    icyg0 = get_iymin_in_fmm(ilevel)
    iczg0 = get_izmin_in_fmm(ilevel)

    DO load = 1, nload(ilevel)
       irx = lddir(1,load,ilevel)
       iry = lddir(2,load,ilevel)
       irz = lddir(3,load,ilevel)
       ncp = 0                        ! cell pair count.

       do icx0 = max(1, -2*nscellx-icxg0+2-abs(irx)), min(nscxdiv, 3*nscellx-icxg0+1+abs(irx))

          if (.not. is_valid_cell(icxg0, icx0, irx, nlx(ilevel)) ) cycle
          icx1 = modulo(icxg0-1 + icx0-1 + irx, nscellx) + 1

          do icy0 = max(1, -2*nscelly-icyg0+2-abs(iry)), min(nscydiv, 3*nscelly-icyg0+1+abs(iry))

             if (.not. is_valid_cell(icyg0, icy0, iry, nly(ilevel)) ) cycle
             icy1 = modulo(icyg0-1 + icy0-1 + iry, nscelly) + 1

             do icz0 = max(1, -2*nscellz-iczg0+2-abs(irz)), min(nsczdiv, 3*nscellz-iczg0+1+abs(irz))

                if (.not. is_valid_cell(iczg0, icz0, irz, nlz(ilevel)) ) cycle
                icz1 = modulo(iczg0-1 + icz0-1 + irz, nscellz) + 1

                ncp = ncp + 1

             end do ! icz0

          end do ! icy0

       end do ! icx0

       if(ncp > 0) then
         if(mod(ncp,8) == 0) then
            nload8p(ilevel) = nload8p(ilevel) + 1
            lddir8p(1,nload8p(ilevel),ilevel) = irx
            lddir8p(2,nload8p(ilevel),ilevel) = iry
            lddir8p(3,nload8p(ilevel),ilevel) = irz
         elseif(mod(ncp,6) == 0) then
            nload6p(ilevel) = nload6p(ilevel) + 1
            lddir6p(1,nload6p(ilevel),ilevel) = irx
            lddir6p(2,nload6p(ilevel),ilevel) = iry
            lddir6p(3,nload6p(ilevel),ilevel) = irz
         elseif(mod(ncp,4) == 0) then
            nload4p(ilevel) = nload4p(ilevel) + 1
            lddir4p(1,nload4p(ilevel),ilevel) = irx
            lddir4p(2,nload4p(ilevel),ilevel) = iry
            lddir4p(3,nload4p(ilevel),ilevel) = irz
         elseif(mod(ncp,3) == 0) then
            nload3p(ilevel) = nload3p(ilevel) + 1
            lddir3p(1,nload3p(ilevel),ilevel) = irx
            lddir3p(2,nload3p(ilevel),ilevel) = iry
            lddir3p(3,nload3p(ilevel),ilevel) = irz
         elseif(mod(ncp,2) == 0) then
            nload2p(ilevel) = nload2p(ilevel) + 1
            lddir2p(1,nload2p(ilevel),ilevel) = irx
            lddir2p(2,nload2p(ilevel),ilevel) = iry
            lddir2p(3,nload2p(ilevel),ilevel) = irz
         elseif(ncp == 1) then
            nload1p(ilevel) = nload1p(ilevel) + 1
            lddir1p(1,nload1p(ilevel),ilevel) = irx
            lddir1p(2,nload1p(ilevel),ilevel) = iry
            lddir1p(3,nload1p(ilevel),ilevel) = irz
         else
            ! error. there should not be this case.
            write(0,*) "ERROR.1: In generate_assorted_relative_cell_list: ", &
                 &     "cell pairs are not the multiple of 8, nor 6, nor 4, ", &
                 &     "nor 3, nor 2 for each relative cell address."
            write(0,*) "ERROR.2: Level; ",ilevel," lgflg: ",lgflg, &
                 &     "  Relative cell address irx,iry,irz; ",irx,iry,irz
            call modylas_abort
         endif
         !write(*,*) "$$$ Number of pairs for relative cell address (irx,iry,irz)=(",irx,",",iry,",",irz,") : ",ncp
       else
            ! no pair case.
       endif

    enddo ! load

!    write(*,*) "+++ Level: ",ilevel," Number of rel.addr for pairs: nload8p; ",nload8p(ilevel), &
!      &        ", nload6p; ",nload6p(ilevel),", nload4p; ",nload4p(ilevel),", nload3p; ",nload3p(ilevel), &
!      &        ", nload2p; ",nload2p(ilevel),", nload1p; ",nload1p(ilevel)

   !----------------------------------------------------------------
    endif    ! ilevel>lgflg.
   !----------------------------------------------------------------

    end do    ! ilevel.

  end subroutine generate_assorted_relative_cell_list_for_M2L

!---------------------------------------------------------------------
!>
!! \brief  Subroutine to assort list of relative cells in M2L.
!! \author Shi-ichi Ichikawa
!<
  ! generation of cell-pair address list assorted by possible number  ...
  ! of cell pairs for each relative cell address.                     ...
  ! the purpose is to execute particular code unrolled by cell-pairs  ...
  ! with unrolled matrix dimensions in M2L, and to avoid frequent     ...
  ! instruction access miss. also, is to avoid the cost of cell-pair  ...
  ! calculation.                                                      ...
  !... 
  subroutine generate_assorted_cellpair_list_for_M2L()
    use mpi_3d_grid, only : npx, npy, npz, ipx, ipy, ipz
    use domain, only : nlevel
    use mpi_tool
    implicit none
    integer(kind=wip) :: ilevel
    integer(kind=wip) :: irx, iry, irz

    integer(kind=wip) :: nomp=1
  ! integer(kind=wip), parameter :: max_ncref = 8   ! maximum referencing cell-range for M2L translation
  ! integer(kind=wip), parameter :: ncdir = 2   ! cell range of direct force calculation

    integer(kind=wip) :: nscxdiv, nscydiv, nsczdiv
    integer(kind=wip) :: nscellx, nscelly, nscellz
    integer(kind=wip) :: icxg0, icyg0, iczg0
    integer(kind=wip) :: icxgo, icygo, iczgo
    integer(kind=wip) :: icxg,  icyg,  iczg
    integer(kind=wip) :: icx0, icy0, icz0, icx1, icy1, icz1
    integer(kind=wip) :: icx1s, icy1s, icz1s
    integer(kind=wip) :: icx1e, icy1e, icz1e

    integer(kind=wip) :: mbd_z, mbd_y, mbd_x
    integer(kind=wip) :: mbd_z2, mbd_y2, mbd_x2
    integer(kind=wip) :: ngord, nrgord
    integer(kind=wip) :: load
    integer(kind=wip) :: ithr
    integer(kind=wip) :: ncp8, ncp6, ncp4, ncp3, ncp2, ncp1

    nsczdiv = wl_size_z(0);  nscydiv = wl_size_y(0);  nscxdiv = wl_size_x(0)
    allocate( cp8_list0(3, ((2*max_ncref+1)**3-(2*ncdir+1)**3)*(nscxdiv*nscydiv*nsczdiv) ,0:nlevel-1) )
    allocate( cp6_list0(3, ((2*max_ncref+1)**3-(2*ncdir+1)**3)*(nscxdiv*nscydiv*nsczdiv) ,0:nlevel-1) )
    allocate( cp4_list0(3, ((2*max_ncref+1)**3-(2*ncdir+1)**3)*(nscxdiv*nscydiv*nsczdiv) ,0:nlevel-1) )
    allocate( cp3_list0(3, ((2*max_ncref+1)**3-(2*ncdir+1)**3)*(nscxdiv*nscydiv*nsczdiv) ,0:nlevel-1) )
    allocate( cp2_list0(3, ((2*max_ncref+1)**3-(2*ncdir+1)**3)*(nscxdiv*nscydiv*nsczdiv) ,0:nlevel-1) )
    allocate( cp1_list0(3, ((2*max_ncref+1)**3-(2*ncdir+1)**3)*(nscxdiv*nscydiv*nsczdiv) ,0:nlevel-1) )
    allocate( cp8_list1(3, ((2*max_ncref+1)**3-(2*ncdir+1)**3)*(nscxdiv*nscydiv*nsczdiv) ,0:nlevel-1) )
    allocate( cp6_list1(3, ((2*max_ncref+1)**3-(2*ncdir+1)**3)*(nscxdiv*nscydiv*nsczdiv) ,0:nlevel-1) )
    allocate( cp4_list1(3, ((2*max_ncref+1)**3-(2*ncdir+1)**3)*(nscxdiv*nscydiv*nsczdiv) ,0:nlevel-1) )
    allocate( cp3_list1(3, ((2*max_ncref+1)**3-(2*ncdir+1)**3)*(nscxdiv*nscydiv*nsczdiv) ,0:nlevel-1) )
    allocate( cp2_list1(3, ((2*max_ncref+1)**3-(2*ncdir+1)**3)*(nscxdiv*nscydiv*nsczdiv) ,0:nlevel-1) )
    allocate( cp1_list1(3, ((2*max_ncref+1)**3-(2*ncdir+1)**3)*(nscxdiv*nscydiv*nsczdiv) ,0:nlevel-1) )
    allocate( cp8_irxyz(3, ((2*max_ncref+1)**3-(2*ncdir+1)**3)*(nscxdiv*nscydiv*nsczdiv) ,0:nlevel-1) )
    allocate( cp6_irxyz(3, ((2*max_ncref+1)**3-(2*ncdir+1)**3)*(nscxdiv*nscydiv*nsczdiv) ,0:nlevel-1) )
    allocate( cp4_irxyz(3, ((2*max_ncref+1)**3-(2*ncdir+1)**3)*(nscxdiv*nscydiv*nsczdiv) ,0:nlevel-1) )
    allocate( cp3_irxyz(3, ((2*max_ncref+1)**3-(2*ncdir+1)**3)*(nscxdiv*nscydiv*nsczdiv) ,0:nlevel-1) )
    allocate( cp2_irxyz(3, ((2*max_ncref+1)**3-(2*ncdir+1)**3)*(nscxdiv*nscydiv*nsczdiv) ,0:nlevel-1) )
    allocate( cp1_irxyz(3, ((2*max_ncref+1)**3-(2*ncdir+1)**3)*(nscxdiv*nscydiv*nsczdiv) ,0:nlevel-1) )

    allocate( cp8_count(0:nlevel-1) )
    allocate( cp6_count(0:nlevel-1) )
    allocate( cp4_count(0:nlevel-1) )
    allocate( cp3_count(0:nlevel-1) )
    allocate( cp2_count(0:nlevel-1) )
    allocate( cp1_count(0:nlevel-1) )

!$  nomp = omp_get_max_threads()

    do ilevel = 0, nlevel-1

   !----------------------------------------------------------------
    if(ilevel <= lgflg) then
   !----------------------------------------------------------------
    nsczdiv = wl_size_z(ilevel);  nscydiv = wl_size_y(ilevel);  nscxdiv = wl_size_x(ilevel)
    nscellz = nscllz(ilevel);  nscelly = nsclly(ilevel);  nscellx = nscllx(ilevel)

    ! global address offset of a cell prior to starting cell.
    icxgo = (nscellx * ipx) / npx
    icygo = (nscelly * ipy) / npy
    iczgo = (nscellz * ipz) / npz
    mbd_x = mbd_x_lvl(ilevel) ; mbd_x2 = 2*mbd_x
    mbd_y = mbd_y_lvl(ilevel) ; mbd_y2 = 2*mbd_y
    mbd_z = mbd_z_lvl(ilevel) ; mbd_z2 = 2*mbd_z

    ncp8 = 0

    do ithr = 0, nomp-1

      do load = 1, nload8p(ilevel)
        irx = lddir8p(1,load,ilevel)     ! x-dir relative cell address of transformation matrix.
        iry = lddir8p(2,load,ilevel)     ! y-dir relative cell address of transformation matrix.
        irz = lddir8p(3,load,ilevel)     ! z-dir relative cell address of transformation matrix.

        if(abs(irx).le.2 .and. abs(iry).le.2 .and. abs(irz).le.2) cycle  ! lddir does not include true case.

      do icx0 = 1-mbd_x, nscxdiv+mbd_x      ! target cell x-address.
        if(nscxdiv == 1      .and. 1 - mbd_x /= icx0 .and. icx0 /= nscxdiv + mbd_x) cycle
        if(nscxdiv <  mbd_x2 .and. 1         <= icx0 .and. icx0 <= nscxdiv        ) cycle
        icxg  = modulo(icxgo + icx0 - 1, nscellx) + 1  ! target cell global x-address.
        ngord  = mod(icxg - 1, mbd_x)                  ! sequence number within cell merge group minus 1.
        nrgord = mbd_x - ngord - 1                     ! reverse sequence number within cell merge group minus 1.
        icx1s = icx0 - mbd_x2 - ngord                  ! first cell of the cell range to be transformed.
        icx1e = icx0 + mbd_x2 + nrgord                 ! final cell of the cell range to be transformed.
        if(icx0 <= mbd_x          ) icx1s = icx0       ! plus-direction half-space transformation.
        if(icx0 >  nscxdiv - mbd_x) icx1e = icx0 - 1   ! minus-direction half-space transformation. inter-process exclusive sontrol.
        icx1 = icx0 + irx                              ! real transforming cell defined by relative cell address.
        if(icx1 < icx1s .or. icx1e < icx1) cycle
  
        do icy0 = 1-mbd_y, nscydiv+mbd_y    ! target cell y-address.
          if(nscydiv == 1     .and. 1 - mbd_y /= icy0 .and. icy0 /= nscydiv + mbd_y) cycle
          if(nscydiv < mbd_y2 .and. 1         <= icy0 .and. icy0 <= nscydiv        ) cycle
          icyg  = modulo(icygo + icy0 - 1, nscelly) + 1
          ngord  = mod(icyg - 1, mbd_y)
          nrgord = mbd_y - ngord - 1
          icy1s = icy0 - mbd_y2 - ngord
          icy1e = icy0 + mbd_y2 + nrgord
          if(icy0 <= mbd_y          ) icy1s = icy0
          if(icy0 >  nscydiv - mbd_y) icy1e = icy0 - 1
          icy1 = icy0 + iry
          if(icy1 < icy1s .or. icy1e < icy1) cycle
  
          do icz0 = 1-mbd_z, nsczdiv+mbd_z  ! target cell z-address.
            if(nsczdiv == 1     .and. 1 - mbd_z /= icz0 .and. icz0 /= nsczdiv + mbd_z) cycle
            if(nsczdiv < mbd_z2 .and. 1         <= icz0 .and. icz0 <= nsczdiv        ) cycle
            iczg  = modulo(iczgo + icz0 - 1, nscellz) + 1
            ngord  = mod(iczg - 1, mbd_z)
            nrgord = mbd_z - ngord - 1
            icz1s = icz0 - mbd_z2 - ngord
            icz1e = icz0 + mbd_z2 + nrgord
            if(icz0 <= mbd_z          ) icz1s = icz0
            if(icz0 >  nsczdiv - mbd_z) icz1e = icz0 - 1
            icz1 = icz0 + irz
            if(icz1 < icz1s .or. icz1e < icz1) cycle
  
            if(mod(load-1, nomp) == ithr) then
              ncp8 = ncp8 + 1
              cp8_list0(1,ncp8,ilevel) = icx0
              cp8_list0(2,ncp8,ilevel) = icy0
              cp8_list0(3,ncp8,ilevel) = icz0
              cp8_list1(1,ncp8,ilevel) = icx1
              cp8_list1(2,ncp8,ilevel) = icy1
              cp8_list1(3,ncp8,ilevel) = icz1
            endif
  
          end do  ! icz0.
        end do  ! icy0.
      end do  ! icx0.

      end do  ! load.

!      write(*,*) "+++ Level: ",ilevel," Number of cell-pairs: ",ncp8," for nload8p. accumulated thread; ",ithr

    end do  ! ithr
    cp8_count(ilevel) = ncp8

    ncp6 = 0

    do ithr = 0, nomp-1

      do load = 1, nload6p(ilevel)
        irx = lddir6p(1,load,ilevel)     ! x-dir relative cell address of transformation matrix.
        iry = lddir6p(2,load,ilevel)     ! y-dir relative cell address of transformation matrix.
        irz = lddir6p(3,load,ilevel)     ! z-dir relative cell address of transformation matrix.
  
!        if(abs(irx).le.2 .and. abs(iry).le.2 .and. abs(irz).le.2) cycle  ! lddir does not include true case.
  
      do icx0 = 1-mbd_x, nscxdiv+mbd_x      ! target cell x-address.
        if(nscxdiv == 1      .and. 1 - mbd_x /= icx0 .and. icx0 /= nscxdiv + mbd_x) cycle
        if(nscxdiv <  mbd_x2 .and. 1         <= icx0 .and. icx0 <= nscxdiv        ) cycle
        icxg  = modulo(icxgo + icx0 - 1, nscellx) + 1  ! target cell global x-address.
        ngord  = mod(icxg - 1, mbd_x)                  ! sequence number within cell merge group minus 1.
        nrgord = mbd_x - ngord - 1                     ! reverse sequence number within cell merge group minus 1.
        icx1s = icx0 - mbd_x2 - ngord                  ! first cell of the cell range to be transformed.
        icx1e = icx0 + mbd_x2 + nrgord                 ! final cell of the cell range to be transformed.
        if(icx0 <= mbd_x          ) icx1s = icx0       ! plus-direction half-space transformation.
        if(icx0 >  nscxdiv - mbd_x) icx1e = icx0 - 1   ! minus-direction half-space transformation. inter-process exclusive sontrol.
        icx1 = icx0 + irx                              ! real transforming cell defined by relative cell address.
        if(icx1 < icx1s .or. icx1e < icx1) cycle
  
        do icy0 = 1-mbd_y, nscydiv+mbd_y    ! target cell y-address.
          if(nscydiv == 1     .and. 1 - mbd_y /= icy0 .and. icy0 /= nscydiv + mbd_y) cycle
          if(nscydiv < mbd_y2 .and. 1         <= icy0 .and. icy0 <= nscydiv        ) cycle
          icyg  = modulo(icygo + icy0 - 1, nscelly) + 1
          ngord  = mod(icyg - 1, mbd_y)
          nrgord = mbd_y - ngord - 1
          icy1s = icy0 - mbd_y2 - ngord
          icy1e = icy0 + mbd_y2 + nrgord
          if(icy0 <= mbd_y          ) icy1s = icy0
          if(icy0 >  nscydiv - mbd_y) icy1e = icy0 - 1
          icy1 = icy0 + iry
          if(icy1 < icy1s .or. icy1e < icy1) cycle
  
          do icz0 = 1-mbd_z, nsczdiv+mbd_z  ! target cell z-address.
            if(nsczdiv == 1     .and. 1 - mbd_z /= icz0 .and. icz0 /= nsczdiv + mbd_z) cycle
            if(nsczdiv < mbd_z2 .and. 1         <= icz0 .and. icz0 <= nsczdiv        ) cycle
            iczg  = modulo(iczgo + icz0 - 1, nscellz) + 1
            ngord  = mod(iczg - 1, mbd_z)
            nrgord = mbd_z - ngord - 1
            icz1s = icz0 - mbd_z2 - ngord
            icz1e = icz0 + mbd_z2 + nrgord
            if(icz0 <= mbd_z          ) icz1s = icz0
            if(icz0 >  nsczdiv - mbd_z) icz1e = icz0 - 1
            icz1 = icz0 + irz
            if(icz1 < icz1s .or. icz1e < icz1) cycle
  
            if(mod(load-1, nomp) == ithr) then
              ncp6 = ncp6 + 1
              cp6_list0(1,ncp6,ilevel) = icx0
              cp6_list0(2,ncp6,ilevel) = icy0
              cp6_list0(3,ncp6,ilevel) = icz0
              cp6_list1(1,ncp6,ilevel) = icx1
              cp6_list1(2,ncp6,ilevel) = icy1
              cp6_list1(3,ncp6,ilevel) = icz1
            endif
  
          end do  ! icz0.
        end do  ! icy0.
      end do  ! icx0.
  
      end do  ! load.

!      write(*,*) "+++ Level: ",ilevel," Number of cell-pairs: ",ncp6," for nload6p, accumulated thread; ",ithr

    end do  ! ithr
    cp6_count(ilevel) = ncp6

    ncp4 = 0

    do ithr = 0, nomp-1

      do load = 1, nload4p(ilevel)
        irx = lddir4p(1,load,ilevel)     ! x-dir relative cell address of transformation matrix.
        iry = lddir4p(2,load,ilevel)     ! y-dir relative cell address of transformation matrix.
        irz = lddir4p(3,load,ilevel)     ! z-dir relative cell address of transformation matrix.
  
        if(abs(irx).le.2 .and. abs(iry).le.2 .and. abs(irz).le.2) cycle  ! lddir does not include true case.
  
      do icx0 = 1-mbd_x, nscxdiv+mbd_x      ! target cell x-address.
        if(nscxdiv == 1      .and. 1 - mbd_x /= icx0 .and. icx0 /= nscxdiv + mbd_x) cycle
        if(nscxdiv <  mbd_x2 .and. 1         <= icx0 .and. icx0 <= nscxdiv        ) cycle
        icxg  = modulo(icxgo + icx0 - 1, nscellx) + 1  ! target cell global x-address.
        ngord  = mod(icxg - 1, mbd_x)                  ! sequence number within cell merge group minus 1.
        nrgord = mbd_x - ngord - 1                     ! reverse sequence number within cell merge group minus 1.
        icx1s = icx0 - mbd_x2 - ngord                  ! first cell of the cell range to be transformed.
        icx1e = icx0 + mbd_x2 + nrgord                 ! final cell of the cell range to be transformed.
        if(icx0 <= mbd_x          ) icx1s = icx0       ! plus-direction half-space transformation.
        if(icx0 >  nscxdiv - mbd_x) icx1e = icx0 - 1   ! minus-direction half-space transformation. inter-process exclusive sontrol.
        icx1 = icx0 + irx                              ! real transforming cell defined by relative cell address.
        if(icx1 < icx1s .or. icx1e < icx1) cycle
  
        do icy0 = 1-mbd_y, nscydiv+mbd_y    ! target cell y-address.
          if(nscydiv == 1     .and. 1 - mbd_y /= icy0 .and. icy0 /= nscydiv + mbd_y) cycle
          if(nscydiv < mbd_y2 .and. 1         <= icy0 .and. icy0 <= nscydiv        ) cycle
          icyg  = modulo(icygo + icy0 - 1, nscelly) + 1
          ngord  = mod(icyg - 1, mbd_y)
          nrgord = mbd_y - ngord - 1
          icy1s = icy0 - mbd_y2 - ngord
          icy1e = icy0 + mbd_y2 + nrgord
          if(icy0 <= mbd_y          ) icy1s = icy0
          if(icy0 >  nscydiv - mbd_y) icy1e = icy0 - 1
          icy1 = icy0 + iry
          if(icy1 < icy1s .or. icy1e < icy1) cycle
  
          do icz0 = 1-mbd_z, nsczdiv+mbd_z  ! target cell z-address.
            if(nsczdiv == 1     .and. 1 - mbd_z /= icz0 .and. icz0 /= nsczdiv + mbd_z) cycle
            if(nsczdiv < mbd_z2 .and. 1         <= icz0 .and. icz0 <= nsczdiv        ) cycle
            iczg  = modulo(iczgo + icz0 - 1, nscellz) + 1
            ngord  = mod(iczg - 1, mbd_z)
            nrgord = mbd_z - ngord - 1
            icz1s = icz0 - mbd_z2 - ngord
            icz1e = icz0 + mbd_z2 + nrgord
            if(icz0 <= mbd_z          ) icz1s = icz0
            if(icz0 >  nsczdiv - mbd_z) icz1e = icz0 - 1
            icz1 = icz0 + irz
            if(icz1 < icz1s .or. icz1e < icz1) cycle
  
            if(mod(load-1, nomp) == ithr) then
              ncp4 = ncp4 + 1
              cp4_list0(1,ncp4,ilevel) = icx0
              cp4_list0(2,ncp4,ilevel) = icy0
              cp4_list0(3,ncp4,ilevel) = icz0
              cp4_list1(1,ncp4,ilevel) = icx1
              cp4_list1(2,ncp4,ilevel) = icy1
              cp4_list1(3,ncp4,ilevel) = icz1
            endif
  
          end do  ! icz0.
        end do  ! icy0.
      end do  ! icx0.
  
      end do  ! load.

!      write(*,*) "+++ Level: ",ilevel," Number of cell-pairs: ",ncp4," for nload4p, accumulated thread; ",ithr

    end do  ! ithr
    cp4_count(ilevel) = ncp4

    ncp3 = 0

    do ithr = 0, nomp-1

      do load = 1, nload3p(ilevel)
        irx = lddir3p(1,load,ilevel)     ! x-dir relative cell address of transformation matrix.
        iry = lddir3p(2,load,ilevel)     ! y-dir relative cell address of transformation matrix.
        irz = lddir3p(3,load,ilevel)     ! z-dir relative cell address of transformation matrix.
  
        if(abs(irx).le.2 .and. abs(iry).le.2 .and. abs(irz).le.2) cycle  ! lddir does not include true case.
  
      do icx0 = 1-mbd_x, nscxdiv+mbd_x      ! target cell x-address.
        if(nscxdiv == 1      .and. 1 - mbd_x /= icx0 .and. icx0 /= nscxdiv + mbd_x) cycle
        if(nscxdiv <  mbd_x2 .and. 1         <= icx0 .and. icx0 <= nscxdiv        ) cycle
        icxg  = modulo(icxgo + icx0 - 1, nscellx) + 1  ! target cell global x-address.
        ngord  = mod(icxg - 1, mbd_x)                  ! sequence number within cell merge group minus 1.
        nrgord = mbd_x - ngord - 1                     ! reverse sequence number within cell merge group minus 1.
        icx1s = icx0 - mbd_x2 - ngord                  ! first cell of the cell range to be transformed.
        icx1e = icx0 + mbd_x2 + nrgord                 ! final cell of the cell range to be transformed.
        if(icx0 <= mbd_x          ) icx1s = icx0       ! plus-direction half-space transformation.
        if(icx0 >  nscxdiv - mbd_x) icx1e = icx0 - 1   ! minus-direction half-space transformation. inter-process exclusive sontrol.
        icx1 = icx0 + irx                              ! real transforming cell defined by relative cell address.
        if(icx1 < icx1s .or. icx1e < icx1) cycle
  
        do icy0 = 1-mbd_y, nscydiv+mbd_y    ! target cell y-address.
          if(nscydiv == 1     .and. 1 - mbd_y /= icy0 .and. icy0 /= nscydiv + mbd_y) cycle
          if(nscydiv < mbd_y2 .and. 1         <= icy0 .and. icy0 <= nscydiv        ) cycle
          icyg  = modulo(icygo + icy0 - 1, nscelly) + 1
          ngord  = mod(icyg - 1, mbd_y)
          nrgord = mbd_y - ngord - 1
          icy1s = icy0 - mbd_y2 - ngord
          icy1e = icy0 + mbd_y2 + nrgord
          if(icy0 <= mbd_y          ) icy1s = icy0
          if(icy0 >  nscydiv - mbd_y) icy1e = icy0 - 1
          icy1 = icy0 + iry
          if(icy1 < icy1s .or. icy1e < icy1) cycle
  
          do icz0 = 1-mbd_z, nsczdiv+mbd_z  ! target cell z-address.
            if(nsczdiv == 1     .and. 1 - mbd_z /= icz0 .and. icz0 /= nsczdiv + mbd_z) cycle
            if(nsczdiv < mbd_z2 .and. 1         <= icz0 .and. icz0 <= nsczdiv        ) cycle
            iczg  = modulo(iczgo + icz0 - 1, nscellz) + 1
            ngord  = mod(iczg - 1, mbd_z)
            nrgord = mbd_z - ngord - 1
            icz1s = icz0 - mbd_z2 - ngord
            icz1e = icz0 + mbd_z2 + nrgord
            if(icz0 <= mbd_z          ) icz1s = icz0
            if(icz0 >  nsczdiv - mbd_z) icz1e = icz0 - 1
            icz1 = icz0 + irz
            if(icz1 < icz1s .or. icz1e < icz1) cycle
  
            if(mod(load-1, nomp) == ithr) then
              ncp3 = ncp3 + 1
              cp3_list0(1,ncp3,ilevel) = icx0
              cp3_list0(2,ncp3,ilevel) = icy0
              cp3_list0(3,ncp3,ilevel) = icz0
              cp3_list1(1,ncp3,ilevel) = icx1
              cp3_list1(2,ncp3,ilevel) = icy1
              cp3_list1(3,ncp3,ilevel) = icz1
            endif
  
          end do  ! icz0.
        end do  ! icy0.
      end do  ! icx0.
  
      end do  ! load.

!      write(*,*) "+++ Level: ",ilevel," Number of cell-pairs: ",ncp3," for nload3p accumulated thread; ",ithr

    end do  ! ithr
    cp3_count(ilevel) = ncp3

    ncp2 = 0

    do ithr = 0, nomp-1

      do load = 1, nload2p(ilevel)
        irx = lddir2p(1,load,ilevel)     ! x-dir relative cell address of transformation matrix.
        iry = lddir2p(2,load,ilevel)     ! y-dir relative cell address of transformation matrix.
        irz = lddir2p(3,load,ilevel)     ! z-dir relative cell address of transformation matrix.
  
        if(abs(irx).le.2 .and. abs(iry).le.2 .and. abs(irz).le.2) cycle  ! lddir does not include true case.
  
      do icx0 = 1-mbd_x, nscxdiv+mbd_x      ! target cell x-address.
        if(nscxdiv == 1      .and. 1 - mbd_x /= icx0 .and. icx0 /= nscxdiv + mbd_x) cycle
        if(nscxdiv <  mbd_x2 .and. 1         <= icx0 .and. icx0 <= nscxdiv        ) cycle
        icxg  = modulo(icxgo + icx0 - 1, nscellx) + 1  ! target cell global x-address.
        ngord  = mod(icxg - 1, mbd_x)                  ! sequence number within cell merge group minus 1.
        nrgord = mbd_x - ngord - 1                     ! reverse sequence number within cell merge group minus 1.
        icx1s = icx0 - mbd_x2 - ngord                  ! first cell of the cell range to be transformed.
        icx1e = icx0 + mbd_x2 + nrgord                 ! final cell of the cell range to be transformed.
        if(icx0 <= mbd_x          ) icx1s = icx0       ! plus-direction half-space transformation.
        if(icx0 >  nscxdiv - mbd_x) icx1e = icx0 - 1   ! minus-direction half-space transformation. inter-process exclusive sontrol.
        icx1 = icx0 + irx                              ! real transforming cell defined by relative cell address.
        if(icx1 < icx1s .or. icx1e < icx1) cycle
  
        do icy0 = 1-mbd_y, nscydiv+mbd_y    ! target cell y-address.
          if(nscydiv == 1     .and. 1 - mbd_y /= icy0 .and. icy0 /= nscydiv + mbd_y) cycle
          if(nscydiv < mbd_y2 .and. 1         <= icy0 .and. icy0 <= nscydiv        ) cycle
          icyg  = modulo(icygo + icy0 - 1, nscelly) + 1
          ngord  = mod(icyg - 1, mbd_y)
          nrgord = mbd_y - ngord - 1
          icy1s = icy0 - mbd_y2 - ngord
          icy1e = icy0 + mbd_y2 + nrgord
          if(icy0 <= mbd_y          ) icy1s = icy0
          if(icy0 >  nscydiv - mbd_y) icy1e = icy0 - 1
          icy1 = icy0 + iry
          if(icy1 < icy1s .or. icy1e < icy1) cycle
  
          do icz0 = 1-mbd_z, nsczdiv+mbd_z  ! target cell z-address.
            if(nsczdiv == 1     .and. 1 - mbd_z /= icz0 .and. icz0 /= nsczdiv + mbd_z) cycle
            if(nsczdiv < mbd_z2 .and. 1         <= icz0 .and. icz0 <= nsczdiv        ) cycle
            iczg  = modulo(iczgo + icz0 - 1, nscellz) + 1
            ngord  = mod(iczg - 1, mbd_z)
            nrgord = mbd_z - ngord - 1
            icz1s = icz0 - mbd_z2 - ngord
            icz1e = icz0 + mbd_z2 + nrgord
            if(icz0 <= mbd_z          ) icz1s = icz0
            if(icz0 >  nsczdiv - mbd_z) icz1e = icz0 - 1
            icz1 = icz0 + irz
            if(icz1 < icz1s .or. icz1e < icz1) cycle
  
            if(mod(load-1, nomp) == ithr) then
              ncp2 = ncp2 + 1
              cp2_list0(1,ncp2,ilevel) = icx0
              cp2_list0(2,ncp2,ilevel) = icy0
              cp2_list0(3,ncp2,ilevel) = icz0
              cp2_list1(1,ncp2,ilevel) = icx1
              cp2_list1(2,ncp2,ilevel) = icy1
              cp2_list1(3,ncp2,ilevel) = icz1
            endif
  
          end do  ! icz0.
        end do  ! icy0.
      end do  ! icx0.
  
      end do  ! load.

!      write(*,*) "+++ Level: ",ilevel," Number of cell-pairs: ",ncp2," for nload2p accumulated thread; ",ithr

    end do  ! ithr
    cp2_count(ilevel) = ncp2

    ncp1 = 0

    do ithr = 0, nomp-1

      do load = 1, nload1p(ilevel)
        irx = lddir1p(1,load,ilevel)     ! x-dir relative cell address of transformation matrix.
        iry = lddir1p(2,load,ilevel)     ! y-dir relative cell address of transformation matrix.
        irz = lddir1p(3,load,ilevel)     ! z-dir relative cell address of transformation matrix.
  
        if(abs(irx).le.2 .and. abs(iry).le.2 .and. abs(irz).le.2) cycle  ! lddir does not include true case.
  
      do icx0 = 1-mbd_x, nscxdiv+mbd_x      ! target cell x-address.
        if(nscxdiv == 1      .and. 1 - mbd_x /= icx0 .and. icx0 /= nscxdiv + mbd_x) cycle
        if(nscxdiv <  mbd_x2 .and. 1         <= icx0 .and. icx0 <= nscxdiv        ) cycle
        icxg  = modulo(icxgo + icx0 - 1, nscellx) + 1  ! target cell global x-address.
        ngord  = mod(icxg - 1, mbd_x)                  ! sequence number within cell merge group minus 1.
        nrgord = mbd_x - ngord - 1                     ! reverse sequence number within cell merge group minus 1.
        icx1s = icx0 - mbd_x2 - ngord                  ! first cell of the cell range to be transformed.
        icx1e = icx0 + mbd_x2 + nrgord                 ! final cell of the cell range to be transformed.
        if(icx0 <= mbd_x          ) icx1s = icx0       ! plus-direction half-space transformation.
        if(icx0 >  nscxdiv - mbd_x) icx1e = icx0 - 1   ! minus-direction half-space transformation. inter-process exclusive sontrol.
        icx1 = icx0 + irx                              ! real transforming cell defined by relative cell address.
        if(icx1 < icx1s .or. icx1e < icx1) cycle
  
        do icy0 = 1-mbd_y, nscydiv+mbd_y    ! target cell y-address.
          if(nscydiv == 1     .and. 1 - mbd_y /= icy0 .and. icy0 /= nscydiv + mbd_y) cycle
          if(nscydiv < mbd_y2 .and. 1         <= icy0 .and. icy0 <= nscydiv        ) cycle
          icyg  = modulo(icygo + icy0 - 1, nscelly) + 1
          ngord  = mod(icyg - 1, mbd_y)
          nrgord = mbd_y - ngord - 1
          icy1s = icy0 - mbd_y2 - ngord
          icy1e = icy0 + mbd_y2 + nrgord
          if(icy0 <= mbd_y          ) icy1s = icy0
          if(icy0 >  nscydiv - mbd_y) icy1e = icy0 - 1
          icy1 = icy0 + iry
          if(icy1 < icy1s .or. icy1e < icy1) cycle
  
          do icz0 = 1-mbd_z, nsczdiv+mbd_z  ! target cell z-address.
            if(nsczdiv == 1     .and. 1 - mbd_z /= icz0 .and. icz0 /= nsczdiv + mbd_z) cycle
            if(nsczdiv < mbd_z2 .and. 1         <= icz0 .and. icz0 <= nsczdiv        ) cycle
            iczg  = modulo(iczgo + icz0 - 1, nscellz) + 1
            ngord  = mod(iczg - 1, mbd_z)
            nrgord = mbd_z - ngord - 1
            icz1s = icz0 - mbd_z2 - ngord
            icz1e = icz0 + mbd_z2 + nrgord
            if(icz0 <= mbd_z          ) icz1s = icz0
            if(icz0 >  nsczdiv - mbd_z) icz1e = icz0 - 1
            icz1 = icz0 + irz
            if(icz1 < icz1s .or. icz1e < icz1) cycle
  
            if(mod(load-1, nomp) == ithr) then
              ncp1 = ncp1 + 1
              cp1_list0(1,ncp1,ilevel) = icx0
              cp1_list0(2,ncp1,ilevel) = icy0
              cp1_list0(3,ncp1,ilevel) = icz0
              cp1_list1(1,ncp1,ilevel) = icx1
              cp1_list1(2,ncp1,ilevel) = icy1
              cp1_list1(3,ncp1,ilevel) = icz1
            endif
  
          end do  ! icz0.
        end do  ! icy0.
      end do  ! icx0.
  
      end do  ! load.

!      write(*,*) "+++ Level: ",ilevel," Number of cell-pairs: ",ncp1," for nload1p accumulated thread; ",ithr

    end do  ! ithr
    cp1_count(ilevel) = ncp1

   !----------------------------------------------------------------
    else    ! ilevel>lgflg.
   !----------------------------------------------------------------
    nsczdiv = wl_size_z(ilevel);  nscydiv = wl_size_y(ilevel);  nscxdiv = wl_size_x(ilevel)
    nscellz = nscllz(ilevel);  nscelly = nsclly(ilevel);  nscellx = nscllx(ilevel)
    
    ! global address of starting cell.
    icxg0 = get_ixmin_in_fmm(ilevel)
    icyg0 = get_iymin_in_fmm(ilevel)
    iczg0 = get_izmin_in_fmm(ilevel)

    ncp8 = 0

    do ithr = 0, nomp-1

    DO load = 1, nload8p(ilevel)
       irx = lddir8p(1,load,ilevel)
       iry = lddir8p(2,load,ilevel)
       irz = lddir8p(3,load,ilevel)

       do icx0 = max(1, -2*nscellx-icxg0+2-abs(irx)), min(nscxdiv, 3*nscellx-icxg0+1+abs(irx))
          if (.not. is_valid_cell(icxg0, icx0, irx, nlx(ilevel)) ) cycle
          icx1 = modulo(icxg0-1 + icx0-1 + irx, nscellx) + 1

          do icy0 = max(1, -2*nscelly-icyg0+2-abs(iry)), min(nscydiv, 3*nscelly-icyg0+1+abs(iry))
             if (.not. is_valid_cell(icyg0, icy0, iry, nly(ilevel)) ) cycle
             icy1 = modulo(icyg0-1 + icy0-1 + iry, nscelly) + 1

             do icz0 = max(1, -2*nscellz-iczg0+2-abs(irz)), min(nsczdiv, 3*nscellz-iczg0+1+abs(irz))
                if (.not. is_valid_cell(iczg0, icz0, irz, nlz(ilevel)) ) cycle
                icz1 = modulo(iczg0-1 + icz0-1 + irz, nscellz) + 1

                if(mod(load-1, nomp) == ithr) then
                  ncp8 = ncp8 + 1
                  cp8_list0(1,ncp8,ilevel) = icx0
                  cp8_list0(2,ncp8,ilevel) = icy0
                  cp8_list0(3,ncp8,ilevel) = icz0
                  cp8_list1(1,ncp8,ilevel) = icx1
                  cp8_list1(2,ncp8,ilevel) = icy1
                  cp8_list1(3,ncp8,ilevel) = icz1
                  cp8_irxyz(1,ncp8,ilevel) = irx
                  cp8_irxyz(2,ncp8,ilevel) = iry
                  cp8_irxyz(3,ncp8,ilevel) = irz
                endif

             end do ! icz0
          end do ! icy0
       end do ! icx0

    enddo ! load

!       write(*,*) "+++ Level: ",ilevel," Number of cell-pairs: ",ncp8," for nload8p accumulated thread; ",ithr

    end do  ! ithr
    cp8_count(ilevel) = ncp8

    ncp6 = 0

    do ithr = 0, nomp-1

    DO load = 1, nload6p(ilevel)
       irx = lddir6p(1,load,ilevel)
       iry = lddir6p(2,load,ilevel)
       irz = lddir6p(3,load,ilevel)

       do icx0 = max(1, -2*nscellx-icxg0+2-abs(irx)), min(nscxdiv, 3*nscellx-icxg0+1+abs(irx))
          if (.not. is_valid_cell(icxg0, icx0, irx, nlx(ilevel)) ) cycle
          icx1 = modulo(icxg0-1 + icx0-1 + irx, nscellx) + 1

          do icy0 = max(1, -2*nscelly-icyg0+2-abs(iry)), min(nscydiv, 3*nscelly-icyg0+1+abs(iry))
             if (.not. is_valid_cell(icyg0, icy0, iry, nly(ilevel)) ) cycle
             icy1 = modulo(icyg0-1 + icy0-1 + iry, nscelly) + 1

             do icz0 = max(1, -2*nscellz-iczg0+2-abs(irz)), min(nsczdiv, 3*nscellz-iczg0+1+abs(irz))
                if (.not. is_valid_cell(iczg0, icz0, irz, nlz(ilevel)) ) cycle
                icz1 = modulo(iczg0-1 + icz0-1 + irz, nscellz) + 1

                if(mod(load-1, nomp) == ithr) then
                  ncp6 = ncp6 + 1
                  cp6_list0(1,ncp6,ilevel) = icx0
                  cp6_list0(2,ncp6,ilevel) = icy0
                  cp6_list0(3,ncp6,ilevel) = icz0
                  cp6_list1(1,ncp6,ilevel) = icx1
                  cp6_list1(2,ncp6,ilevel) = icy1
                  cp6_list1(3,ncp6,ilevel) = icz1
                  cp6_irxyz(1,ncp6,ilevel) = irx
                  cp6_irxyz(2,ncp6,ilevel) = iry
                  cp6_irxyz(3,ncp6,ilevel) = irz
                endif

             end do ! icz0
          end do ! icy0
       end do ! icx0

    enddo ! load

!       write(*,*) "+++ Level: ",ilevel," Number of cell-pairs: ",ncp6," for nload6p accumulated thread; ",ithr

    end do  ! ithr
    cp6_count(ilevel) = ncp6

    ncp4 = 0

    do ithr = 0, nomp-1

    DO load = 1, nload4p(ilevel)
       irx = lddir4p(1,load,ilevel)
       iry = lddir4p(2,load,ilevel)
       irz = lddir4p(3,load,ilevel)

       do icx0 = max(1, -2*nscellx-icxg0+2-abs(irx)), min(nscxdiv, 3*nscellx-icxg0+1+abs(irx))
          if (.not. is_valid_cell(icxg0, icx0, irx, nlx(ilevel)) ) cycle
          icx1 = modulo(icxg0-1 + icx0-1 + irx, nscellx) + 1

          do icy0 = max(1, -2*nscelly-icyg0+2-abs(iry)), min(nscydiv, 3*nscelly-icyg0+1+abs(iry))
             if (.not. is_valid_cell(icyg0, icy0, iry, nly(ilevel)) ) cycle
             icy1 = modulo(icyg0-1 + icy0-1 + iry, nscelly) + 1

             do icz0 = max(1, -2*nscellz-iczg0+2-abs(irz)), min(nsczdiv, 3*nscellz-iczg0+1+abs(irz))
                if (.not. is_valid_cell(iczg0, icz0, irz, nlz(ilevel)) ) cycle
                icz1 = modulo(iczg0-1 + icz0-1 + irz, nscellz) + 1

                if(mod(load-1, nomp) == ithr) then
                  ncp4 = ncp4 + 1
                  cp4_list0(1,ncp4,ilevel) = icx0
                  cp4_list0(2,ncp4,ilevel) = icy0
                  cp4_list0(3,ncp4,ilevel) = icz0
                  cp4_list1(1,ncp4,ilevel) = icx1
                  cp4_list1(2,ncp4,ilevel) = icy1
                  cp4_list1(3,ncp4,ilevel) = icz1
                  cp4_irxyz(1,ncp4,ilevel) = irx
                  cp4_irxyz(2,ncp4,ilevel) = iry
                  cp4_irxyz(3,ncp4,ilevel) = irz
                endif

             end do ! icz0
          end do ! icy0
       end do ! icx0

    enddo ! load

!       write(*,*) "+++ Level: ",ilevel," Number of cell-pairs: ",ncp4," for nload4p accumulated thread; ",ithr

    end do  ! ithr
    cp4_count(ilevel) = ncp4

    ncp3 = 0

    do ithr = 0, nomp-1

    DO load = 1, nload3p(ilevel)
       irx = lddir3p(1,load,ilevel)
       iry = lddir3p(2,load,ilevel)
       irz = lddir3p(3,load,ilevel)

       do icx0 = max(1, -2*nscellx-icxg0+2-abs(irx)), min(nscxdiv, 3*nscellx-icxg0+1+abs(irx))
          if (.not. is_valid_cell(icxg0, icx0, irx, nlx(ilevel)) ) cycle
          icx1 = modulo(icxg0-1 + icx0-1 + irx, nscellx) + 1

          do icy0 = max(1, -2*nscelly-icyg0+2-abs(iry)), min(nscydiv, 3*nscelly-icyg0+1+abs(iry))
             if (.not. is_valid_cell(icyg0, icy0, iry, nly(ilevel)) ) cycle
             icy1 = modulo(icyg0-1 + icy0-1 + iry, nscelly) + 1

             do icz0 = max(1, -2*nscellz-iczg0+2-abs(irz)), min(nsczdiv, 3*nscellz-iczg0+1+abs(irz))
                if (.not. is_valid_cell(iczg0, icz0, irz, nlz(ilevel)) ) cycle
                icz1 = modulo(iczg0-1 + icz0-1 + irz, nscellz) + 1

                if(mod(load-1, nomp) == ithr) then
                  ncp3 = ncp3 + 1
                  cp3_list0(1,ncp3,ilevel) = icx0
                  cp3_list0(2,ncp3,ilevel) = icy0
                  cp3_list0(3,ncp3,ilevel) = icz0
                  cp3_list1(1,ncp3,ilevel) = icx1
                  cp3_list1(2,ncp3,ilevel) = icy1
                  cp3_list1(3,ncp3,ilevel) = icz1
                  cp3_irxyz(1,ncp3,ilevel) = irx
                  cp3_irxyz(2,ncp3,ilevel) = iry
                  cp3_irxyz(3,ncp3,ilevel) = irz
                endif

             end do ! icz0
          end do ! icy0
       end do ! icx0

    enddo ! load

!       write(*,*) "+++ Level: ",ilevel," Number of cell-pairs: ",ncp3," for nload3p accumulated thread; ",ithr

    end do  ! ithr
    cp3_count(ilevel) = ncp3

    ncp2 = 0

    do ithr = 0, nomp-1

    DO load = 1, nload2p(ilevel)
       irx = lddir2p(1,load,ilevel)
       iry = lddir2p(2,load,ilevel)
       irz = lddir2p(3,load,ilevel)

       do icx0 = max(1, -2*nscellx-icxg0+2-abs(irx)), min(nscxdiv, 3*nscellx-icxg0+1+abs(irx))
          if (.not. is_valid_cell(icxg0, icx0, irx, nlx(ilevel)) ) cycle
          icx1 = modulo(icxg0-1 + icx0-1 + irx, nscellx) + 1

          do icy0 = max(1, -2*nscelly-icyg0+2-abs(iry)), min(nscydiv, 3*nscelly-icyg0+1+abs(iry))
             if (.not. is_valid_cell(icyg0, icy0, iry, nly(ilevel)) ) cycle
             icy1 = modulo(icyg0-1 + icy0-1 + iry, nscelly) + 1

             do icz0 = max(1, -2*nscellz-iczg0+2-abs(irz)), min(nsczdiv, 3*nscellz-iczg0+1+abs(irz))
                if (.not. is_valid_cell(iczg0, icz0, irz, nlz(ilevel)) ) cycle
                icz1 = modulo(iczg0-1 + icz0-1 + irz, nscellz) + 1

                if(mod(load-1, nomp) == ithr) then
                  ncp2 = ncp2 + 1
                  cp2_list0(1,ncp2,ilevel) = icx0
                  cp2_list0(2,ncp2,ilevel) = icy0
                  cp2_list0(3,ncp2,ilevel) = icz0
                  cp2_list1(1,ncp2,ilevel) = icx1
                  cp2_list1(2,ncp2,ilevel) = icy1
                  cp2_list1(3,ncp2,ilevel) = icz1
                  cp2_irxyz(1,ncp2,ilevel) = irx
                  cp2_irxyz(2,ncp2,ilevel) = iry
                  cp2_irxyz(3,ncp2,ilevel) = irz
                endif

             end do ! icz0
          end do ! icy0
       end do ! icx0

    enddo ! load

!       write(*,*) "+++ Level: ",ilevel," Number of cell-pairs: ",ncp2," for nload2p accumulated thread; ",ithr

    end do  ! ithr
    cp2_count(ilevel) = ncp2

    ncp1 = 0

    do ithr = 0, nomp-1

    DO load = 1, nload1p(ilevel)
       irx = lddir1p(1,load,ilevel)
       iry = lddir1p(2,load,ilevel)
       irz = lddir1p(3,load,ilevel)

       do icx0 = max(1, -2*nscellx-icxg0+2-abs(irx)), min(nscxdiv, 3*nscellx-icxg0+1+abs(irx))
          if (.not. is_valid_cell(icxg0, icx0, irx, nlx(ilevel)) ) cycle
          icx1 = modulo(icxg0-1 + icx0-1 + irx, nscellx) + 1

          do icy0 = max(1, -2*nscelly-icyg0+2-abs(iry)), min(nscydiv, 3*nscelly-icyg0+1+abs(iry))
             if (.not. is_valid_cell(icyg0, icy0, iry, nly(ilevel)) ) cycle
             icy1 = modulo(icyg0-1 + icy0-1 + iry, nscelly) + 1

             do icz0 = max(1, -2*nscellz-iczg0+2-abs(irz)), min(nsczdiv, 3*nscellz-iczg0+1+abs(irz))
                if (.not. is_valid_cell(iczg0, icz0, irz, nlz(ilevel)) ) cycle
                icz1 = modulo(iczg0-1 + icz0-1 + irz, nscellz) + 1

                if(mod(load-1, nomp) == ithr) then
                  ncp1 = ncp1 + 1
                  cp1_list0(1,ncp1,ilevel) = icx0
                  cp1_list0(2,ncp1,ilevel) = icy0
                  cp1_list0(3,ncp1,ilevel) = icz0
                  cp1_list1(1,ncp1,ilevel) = icx1
                  cp1_list1(2,ncp1,ilevel) = icy1
                  cp1_list1(3,ncp1,ilevel) = icz1
                  cp1_irxyz(1,ncp1,ilevel) = irx
                  cp1_irxyz(2,ncp1,ilevel) = iry
                  cp1_irxyz(3,ncp1,ilevel) = irz
                endif

             end do ! icz0
          end do ! icy0
       end do ! icx0

    enddo ! load

!       write(*,*) "+++ Level: ",ilevel," Number of cell-pairs: ",ncp1," for nload1p accumulated thread; ",ithr

    end do  ! ithr
    cp1_count(ilevel) = ncp1

   !----------------------------------------------------------------
    endif    ! ilevel>lgflg.
   !----------------------------------------------------------------

    end do    ! ilevel.

  end subroutine generate_assorted_cellpair_list_for_M2L
#endif /**MOMENT**/

!---------------------------------------------------------------------
!>
!! \brief  Subroutine to initialize translation operators in
!!         M2M, M2L, L2L, and L2P.
!! \author Tatsuya Sakashita, Noriyuki Yoshii
!<
  subroutine initialize_translation_operators()
!---------------------------------------------------------------------
    use domain, only : ncellx, ncelly, ncellz, nlevel, ixmin, iymin, izmin
    use unit_cell, only : cellx, celly, cellz
    use spherical_harmonics, only : cart2angle
    use regular_singular_solid_harmonics_north
    use regular_singular_rotation_1dim
    use regular_singular_solid_harmonics
    use fmm_l_m_index
    implicit none
    real(kind=dp) :: r_plus_re, r_plus_im, r_minus_re, r_minus_im
    real(kind=dp) :: s_plus_re, s_plus_im, s_minus_re, s_minus_im
    real(kind=dp) :: re_tmp, im_tmp
    real(kind=dp) :: xtb,ytb,ztb,rad,the,csthe,phi
    real(kind=dp) :: dx_supercell, dy_supercell, dz_supercell
    real(kind=dp) :: sum_rot1i, sum_rot1r, sum_rot2i, sum_rot2r   ! add for SIMD_BPREC_TEST.
    real(kind=dp) :: sum_shml, tmp                                ! add for SIMD_BPREC_TEST.
    integer(kind=wip) :: mnlx,mnly,mnlz
    integer(kind=wip) :: j, k, l, m, m1, m2
    integer(kind=wip) :: ix, iy, iz
    integer(kind=wip) :: il
    integer(kind=wip) :: ind_1dim, n   ! add for SIMD_BPREC_TEST.
!M2M
    allocate(shmm1(lm_length,lm_length,0:2,0:2,0:2,0:nlevel-1))
    allocate(shmm2(lm_length,lm_length,0:2,0:2,0:2,0:nlevel-1))
!M2L
    allocate(shml(0:2*nmax,-8:8,-8:8,-8:8,0:nlevel-1))
    allocate(shml1(lm_length,2,lm_length,-8:8,-8:8,-8:8,0:nlevel-1))
    allocate(shml2(lm_length,2,lm_length,-8:8,-8:8,-8:8,0:nlevel-1))

    allocate(sh_rot1((nmax+1)*(nmax+2)*(2*nmax+3)/6,2,-8:8,-8:8,-8:8,0:nlevel-1))
    allocate(sh_rot2((nmax+1)*(nmax+2)*(2*nmax+3)/6,2,-8:8,-8:8,-8:8,0:nlevel-1))
    allocate(sh_inv_rot1((nmax+1)*(nmax+2)*(2*nmax+3)/6,2,-8:8,-8:8,-8:8,0:nlevel-1))
    allocate(sh_inv_rot2((nmax+1)*(nmax+2)*(2*nmax+3)/6,2,-8:8,-8:8,-8:8,0:nlevel-1))

!L2L
    allocate(xtall(0:2,0:2,0:2,1:nlevel))
    allocate(ytall(0:2,0:2,0:2,1:nlevel))
    allocate(ztall(0:2,0:2,0:2,1:nlevel))
    allocate(shll1(lm_length,lm_length,0:2,0:2,0:2,1:nlevel))
    allocate(shll2(lm_length,lm_length,0:2,0:2,0:2,1:nlevel))

! for P2M, L2P
    dxcell = cellx / dble(ncellx)
    dycell = celly / dble(ncelly)
    dzcell = cellz / dble(ncellz)

#if defined(PRECISION_M2L_MIX) || defined(PRECISION_M2L_SP)
      dxcell = dxcell * 1d10
      dycell = dycell * 1d10
      dzcell = dzcell * 1d10
#endif

!
! for P2M, L2P
!
#if defined(PRECISION_M2L_MIX) || defined(PRECISION_M2L_SP)
      x_center_offset = (ixmin-1-0.5d0) * dxcell - 0.5d0 * cellx * 1d10
      y_center_offset = (iymin-1-0.5d0) * dycell - 0.5d0 * celly * 1d10
      z_center_offset = (izmin-1-0.5d0) * dzcell - 0.5d0 * cellz * 1d10
#else
      x_center_offset = (ixmin-1-0.5d0) * dxcell - 0.5d0 * cellx
      y_center_offset = (iymin-1-0.5d0) * dycell - 0.5d0 * celly
      z_center_offset = (izmin-1-0.5d0) * dzcell - 0.5d0 * cellz
#endif
    
!
! M2M
!
    do il=0,nlevel-1
       !debug
       !       if(myrank==0) write(*,*) 'mcell_size_x,y,z'
       !       if(myrank==0) write(*,*) mcell_size_x,mcell_size_y,mcell_size_z
       !       if(myrank==0) write(*,*) 'M2M checking ',il,'-th level'
       !debug

       dx_supercell = dxcell * mcells_x(il)
       dy_supercell = dycell * mcells_y(il)
       dz_supercell = dzcell * mcells_z(il)
       
       do ix=0,nlx(il)-1
          xtb = ( ix + 0.5d0*(-nlx(il)+1) ) * dx_supercell
          do iy=0,nly(il)-1
             ytb = ( iy + 0.5d0*(-nly(il)+1) ) * dy_supercell
             do iz=0,nlz(il)-1
                ztb = ( iz + 0.5d0*(-nlz(il)+1) ) * dz_supercell
                !       if(myrank==0) write(*,*) xtb,ytb,ztb
                call cart2angle(-xtb,-ytb,-ztb,rad,the,csthe,phi)
                do j=0,nmax
                   do k=0,j
                      m1 = translate_l_m_to_1dim(j, k)
                      do l=0,j
                         m2 = translate_l_m_to_1dim(l, 0)
                         call regular_harmonics_reals(j-l, k-0, rad, csthe, phi, re_tmp, im_tmp) ! When m = 0
                         shmm1(m2,m1,ix,iy,iz,il) = dcmplx(re_tmp, im_tmp)
                         shmm2(m2,m1,ix,iy,iz,il) = dcmplx(0d0, 0d0)
                         do m=max(1,-(j-l)+k), min(l,(j-l)+k)
                            m2 = translate_l_m_to_1dim(l, m)
                            call regular_harmonics_reals(j-l, k-m, rad, csthe, phi, r_plus_re, r_plus_im)
                            call regular_harmonics_pz_reals(j-l, k+m, rad, csthe, phi, re_tmp, im_tmp)
                            r_minus_re = (-1)**m * re_tmp
                            r_minus_im = (-1)**m * im_tmp
                            if (k /= 0) then
                               shmm1(m2,m1,ix,iy,iz,il) = dcmplx( r_plus_re + r_minus_re, r_plus_im + r_minus_im)
                               shmm2(m2,m1,ix,iy,iz,il) = dcmplx(-r_plus_im + r_minus_im, r_plus_re - r_minus_re)
                            else
                               shmm1(m2,m1,ix,iy,iz,il) = dcmplx( r_plus_re + r_minus_re, 0d0)
                               shmm2(m2,m1,ix,iy,iz,il) = dcmplx(-r_plus_im + r_minus_im, 0d0)
                            endif
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo ! il
!
! M2L
!
    do il=0,nlevel-1
       mnlx=2*nlx(il)+(nlx(il)-1)
       mnly=2*nly(il)+(nly(il)-1)
       mnlz=2*nlz(il)+(nlz(il)-1)

       dx_supercell = dxcell * mcells_x(il)
       dy_supercell = dycell * mcells_y(il)
       dz_supercell = dzcell * mcells_z(il)
       
       !debug
       !       if(myrank==0) write(*,*) 'M2L checking ',il,'-th level'
       !       if(myrank==0) write(*,*)
       !    &  dxcell*mcell_size_x,dycell*mcell_size_y,dzcell*mcell_size_z
       !debug
       do ix=-mnlx,mnlx
          xtb = ix * dx_supercell
          do iy=-mnly,mnly
             ytb = iy * dy_supercell
             do iz=-mnlz,mnlz
                ztb = iz * dz_supercell
                !       if(myrank==0) write(*,*) xtb,ytb,ztb
                call cart2angle(-xtb,-ytb,-ztb,rad,the,csthe,phi)
                if(rad.ne.0d0)then
                   do j=0,2*nmax
                      shml(j,iz,iy,ix,il) = singular_harmonics_north(j, rad)
                   enddo
                   call generate_rotation_matrix_1dim(the, phi, nmax, &
                        & sh_rot1(:,1,iz,iy,ix,il), sh_rot2(:,1,iz,iy,ix,il), &
                        & sh_inv_rot1(:,1,iz,iy,ix,il), sh_inv_rot2(:,1,iz,iy,ix,il) )

                endif
             enddo
          enddo
       enddo
    enddo ! il
!
! L2L
!
    do il=1,nlevel
       !debug
       !       if(myrank==0) write(*,*) 'L2L checking ',il,'-th level'
       !       if(myrank==0) write(*,*) 'nlxyz=', nlx(il-1),nly(il-1),nlz(il-1)
       !debug

       dx_supercell = dxcell * mcells_x(il)
       dy_supercell = dycell * mcells_y(il)
       dz_supercell = dzcell * mcells_z(il)
       
       do ix=0,nlx(il-1)-1
          xtb = get_shift_from_center(nlx(il-1), ix) * dx_supercell
          do iy=0,nly(il-1)-1
             ytb = get_shift_from_center(nly(il-1), iy) * dy_supercell
             do iz=0,nlz(il-1)-1
                ztb = get_shift_from_center(nlz(il-1), iz) * dz_supercell
!L2L
xtall(ix,iy,iz,il)=-xtb
ytall(ix,iy,iz,il)=-ytb
ztall(ix,iy,iz,il)=-ztb
                call cart2angle(-xtb,-ytb,-ztb,rad,the,csthe,phi)
                do j=0,nmax
                   do k=0,j
                      m1 = translate_l_m_to_1dim(j, k)
                      do l=j,nmax
                         call regular_harmonics_pz_reals(l-j, 0+k, rad, csthe, phi, re_tmp, im_tmp) ! When m=0, we need R(n-j, 0-k)
                         m2 = translate_l_m_to_1dim(l, 0)
                         shll1(m2,m1,ix,iy,iz,il) = dcmplx((-1)**k *re_tmp, (-1)**k *(-im_tmp) )
                         shll2(m2,m1,ix,iy,iz,il) = dcmplx(0d0, 0d0)
                         do m=1,l
                            call regular_harmonics_reals(l-j, m-k, rad, csthe, phi, r_plus_re, r_plus_im)
                            call regular_harmonics_pz_reals(l-j, m+k, rad, csthe, phi, re_tmp, im_tmp)
                            r_minus_re = (-1)**k * re_tmp
                            r_minus_im = (-1)**k * (- im_tmp)
                            m2 = translate_l_m_to_1dim(l, m)
                            if (k /= 0) then
                               shll1(m2,m1,ix,iy,iz,il) = dcmplx( r_plus_re + r_minus_re, r_plus_im + r_minus_im)
                               shll2(m2,m1,ix,iy,iz,il) = dcmplx(-r_plus_im + r_minus_im, r_plus_re - r_minus_re)
                            else
                               shll1(m2,m1,ix,iy,iz,il) = dcmplx( r_plus_re + r_minus_re, 0d0 )
                               shll2(m2,m1,ix,iy,iz,il) = dcmplx(-r_plus_im + r_minus_im, 0d0 )
                            endif
                         enddo ! m
                      enddo ! l
                   enddo ! k
                enddo ! j
             enddo ! iz
          enddo ! iy
       enddo ! ix
    enddo ! il

  end subroutine initialize_translation_operators
!---------------------------------------------------------------------
! << subroutines called in init_fmm
!---------------------------------------------------------------------

!---------------------------------------------------------------------
! Functions
!---------------------------------------------------------------------
  integer function get_ixmin_in_fmm(ilevel)
    use domain, only : ixmin
    implicit none
    integer, intent(in) :: ilevel
    
    get_ixmin_in_fmm = (ixmin-1) / mcells_x(ilevel) + 1
  end function get_ixmin_in_fmm
!---------------------------------------------------------------------
  integer function get_iymin_in_fmm(ilevel)
    use domain, only : iymin
    implicit none
    integer, intent(in) :: ilevel
    
    get_iymin_in_fmm = (iymin-1) / mcells_y(ilevel) + 1
  end function get_iymin_in_fmm
!---------------------------------------------------------------------
  integer function get_izmin_in_fmm(ilevel)
    use domain, only : izmin
    implicit none
    integer, intent(in) :: ilevel
    
    get_izmin_in_fmm = (izmin-1) / mcells_z(ilevel) + 1
  end function get_izmin_in_fmm
!---------------------------------------------------------------------
  logical function is_valid_cell(icxg0, icx0, ic, npow)
    implicit none
    integer, intent(in) :: icxg0, icx0, ic, npow
    integer :: ni
    integer :: i_mod_orig, i_mod_mine

    ni = icxg0 + icx0 - 1
    i_mod_orig = mod(ni-1,npow)   ! =1,2 or 1,2,3 when ic <= 0.  0-started index
    i_mod_mine = i_mod_orig
    if(ic > 0)  i_mod_mine = (npow-1) - i_mod_orig ! 0-started index

    is_valid_cell = iabs(ic) <= i_mod_mine + 2 * npow
  end function is_valid_cell
!---------------------------------------------------------------------
  ! difference of center of child cell and of my cell
  ! n : divided number
  ! i : cell index (i=0,...,n-1)
  real(8) pure elemental function get_shift_from_center(n, i)
    implicit none
    integer(4), intent(in) :: n, i
    get_shift_from_center = 0.5d0 - (0.5d0+dble(i)) / n
  end function get_shift_from_center

!---------------------------------------------------------------------
end module fmm_far
!---------------------------------------------------------------------
