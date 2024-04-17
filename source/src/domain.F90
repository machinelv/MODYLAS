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
!! \brief  Module and subroutines to set up domain decomposition. 
!<
!----------------------------------------------------------------------
!>
!! \brief  Module which relates to domain decomposition.
!! \author Yoshimichi Andoh
!<
!----------------------------------------------------------------------
module domain
  use omp_lib
  use mpi_tool
  implicit none

  integer(4), protected :: ixmin,iymin,izmin
  integer(4), protected :: ixmax,iymax,izmax

  integer(4), protected :: lxdiv,lydiv,lzdiv

  integer(4) :: ndirect=2, nlevel=0
  logical :: nlevel_input=.false.
  integer(4) :: ncellz=0,ncelly=0,ncellx=0
  logical :: ncellz_input=.false.
  logical :: ncelly_input=.false.
  logical :: ncellx_input=.false.
  integer(4) :: nimage=96
  integer(4) :: max_nsegments_per_cell
  integer(4) :: max_nsegments_per_cell_init=250

contains

  integer(4) function get_ixmin(ipx)
    use mpi_3d_grid, only : npx
    implicit none
    integer(4), intent(in) :: ipx
    integer(4) :: n_cell_per_px
    n_cell_per_px = ncellx / npx
    get_ixmin = n_cell_per_px * ipx + 1
  end function get_ixmin

  integer(4) function get_ixmax(ipx)
    use mpi_3d_grid, only : npx
    implicit none
    integer(4), intent(in) :: ipx
    integer(4) :: n_cell_per_px
    n_cell_per_px = ncellx / npx
    get_ixmax = n_cell_per_px * (ipx+1) -1 + 1
  end function get_ixmax

  integer(4) function get_iymin(ipy)
    use mpi_3d_grid, only : npy
    implicit none
    integer(4), intent(in) :: ipy
    integer(4) :: n_cell_per_py
    n_cell_per_py = ncelly / npy
    get_iymin = n_cell_per_py * ipy + 1
  end function get_iymin

  integer(4) function get_iymax(ipy)
    use mpi_3d_grid, only : npy
    implicit none
    integer(4), intent(in) :: ipy
    integer(4) :: n_cell_per_py
    n_cell_per_py = ncelly / npy
    get_iymax = n_cell_per_py * (ipy+1) -1 + 1
  end function get_iymax

  integer(4) function get_izmin(ipz)
    use mpi_3d_grid, only : npz
    implicit none
    integer(4), intent(in) :: ipz
    integer(4) :: n_cell_per_pz
    n_cell_per_pz = ncellz / npz
    get_izmin = n_cell_per_pz * ipz + 1
  end function get_izmin

  integer(4) function get_izmax(ipz)
    use mpi_3d_grid, only : npz
    implicit none
    integer(4), intent(in) :: ipz
    integer(4) :: n_cell_per_pz
    n_cell_per_pz = ncellz / npz
    get_izmax = n_cell_per_pz * (ipz+1) -1 + 1
  end function get_izmax

  integer(4) function get_my_ncell()
    implicit none
    get_my_ncell = lxdiv * lydiv * lzdiv
  end function get_my_ncell

  subroutine translate_cell_local_i_to_global_xyz(ic0, icx, icy, icz)
    implicit none
    integer(4), intent(in) :: ic0
    integer(4), intent(out) :: icx, icy, icz
    icx = mod(ic0-1,lxdiv) ! 0-based index
    icy = (ic0-1) / lxdiv
    icz = icy / lydiv ! 0-based index
    icy = mod(icy,lydiv) ! 0-based index

    icx = icx + ixmin ! 1-based index
    icy = icy + iymin ! 1-based index
    icz = icz + izmin ! 1-based index
  end subroutine translate_cell_local_i_to_global_xyz

  logical function is_my_cell(icx, icy, icz)
    use mpi_3d_grid, only : ipx, ipy, ipz
    implicit none
    integer(4), intent(in) :: icx, icy, icz
    is_my_cell = (ixmin <= icx) .and. (icx <= ixmax) .and. &
         (iymin <= icy) .and. (icy <= iymax) .and. &
         (izmin <= icz) .and. (icz <= izmax)
  end function is_my_cell

!----------------------------------------------------------------------
  subroutine fmod_set_fmm_ndirect(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    ndirect = ivalue
  end subroutine fmod_set_fmm_ndirect
!----------------------------------------------------------------------
  subroutine fmod_set_fmm_ncellx(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    ncellx = ivalue
    ncellx_input=.true.
  end subroutine fmod_set_fmm_ncellx
!----------------------------------------------------------------------
  subroutine fmod_set_fmm_ncelly(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    ncelly = ivalue
    ncelly_input=.true.
  end subroutine fmod_set_fmm_ncelly
!----------------------------------------------------------------------
  subroutine fmod_set_fmm_ncellz(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    ncellz = ivalue
    ncellz_input=.true.
  end subroutine fmod_set_fmm_ncellz
!----------------------------------------------------------------------
  subroutine fmod_set_ncell(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue

    call fmod_set_fmm_ncellx( ivalue )
    call fmod_set_fmm_ncelly( ivalue )
    call fmod_set_fmm_ncellz( ivalue )
  end subroutine fmod_set_ncell
!----------------------------------------------------------------------
  subroutine fmod_set_maxnsegmentspercell(ivalue)
    implicit none
    integer(4) :: ivalue
    max_nsegments_per_cell_init = ivalue
    if(myrank==0) write(*,*) &
      'WARNING: max_nsegments_per_cell_init is set to ', ivalue
  end subroutine fmod_set_maxnsegmentspercell
!----------------------------------------------------------------------
  subroutine fmod_set_maxsegments
    use segments, only : nsegments
    use margin_sizes
    implicit none
    integer(4) :: ivalue
    ivalue = ncellx*ncelly*ncellz
    max_nsegments_per_cell = max(int(nsegments/ivalue*margin_mnpc), &
        max_nsegments_per_cell_init)
  end subroutine fmod_set_maxsegments
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to check domain condition.
!! \author Yoshimichi Andoh
!<
!----------------------------------------------------------------------
  subroutine init_fmm_domain_div()
    use mpi_3d_grid, only : set_ip_x_y_z, npx, npy, npz, ipx, ipy, ipz, get_ip
    implicit none
    include 'mpif.h'
    integer(4) :: ierr

    !## checking mddef input
    if(nlevel_input .and. &
         & .not.(ncellx_input .and. ncelly_input .and. ncellz_input) )then
       if(myrank==0)then
          write(0,*) 'ERROR: nlevel is set, but ncellx,ncelly,ncellz ', &
               &               'are not set in *.mddef'
       endif
       call modylas_abort()
    endif
    !## checking mddef input
    if(.not. nlevel_input .and. &
         & .not.(ncellx_input .and. ncelly_input .and. ncellz_input) )then
       if(myrank==0)then
          write(0,*) 'ERROR: ncellx,ncelly,ncellz are set, but ', &
               &               'nlevel is not set in *.mddef'
       endif
       call modylas_abort()
    endif

    lxdiv = ncellx/npx
    lydiv = ncelly/npy
    lzdiv = ncellz/npz

    call mpi_barrier(mpi_comm_world,ierr)

    ixmin = get_ixmin(ipx)
    ixmax = get_ixmax(ipx)
    iymin = get_iymin(ipy)
    iymax = get_iymax(ipy)
    izmin = get_izmin(ipz)
    izmax = get_izmax(ipz)
  end subroutine init_fmm_domain_div

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to check subcell length, comparing with cutoff
!! \author Yoshimichi Andoh
!<
  subroutine check_cutofflength
!----------------------------------------------------------------------
    use cutoff_radius
    use unit_cell
    implicit none
    real(8) :: cellxd,cellyd,cellzd
    integer(4) :: ierr=0
    cellxd=cellx/dble(ncellx)
    cellyd=celly/dble(ncelly)
    cellzd=cellz/dble(ncellz)

    cellxd=2d0*cellxd
    cellyd=2d0*cellyd
    cellzd=2d0*cellzd
    if(cutrad.eq.0d0)then
       if(myrank==0)then
          write(0,*)'ERROR: LJ-cutoff length=0, which is unlikely situation.'
       endif
       ierr=1 !call mpiend
    endif
! TEMP COMMENT K. FUJIMOTO
    if(    cellxd.lt.cutrad)then  ! Cubic cell only
       if(myrank==0)then
          write(0,*)'ERROR: Length of 2*subcellx=', cellxd*1d+10, ' Aungstrom'
          write(0,*)'LJ-cutoff length=', cutrad*1d+10, ' Aungstrom'
          write(0,*)'This situation is unlikely, reduce ncellx value.'
       endif
       ierr=1 !call mpiend
    elseif(cellyd.lt.cutrad)then  ! Cubic cell only
       if(myrank==0)then
          write(0,*)'ERROR: Length of 2*subcelly=', cellyd*1d+10, ' Aungstrom'
          write(0,*)'LJ-cutoff length=', cutrad*1d+10, ' Aungstrom'
          write(0,*)'This situation is unlikely, reduce ncelly value.'
       endif
       ierr=1 !call mpiend
    elseif(cellzd.lt.cutrad)then  ! Cubic cell only
       if(myrank==0)then
          write(0,*)'ERROR: Length of 2*subcellz=', cellzd*1d+10, ' Aungstrom'
          write(0,*)'LJ-cutoff length=', cutrad*1d+10, ' Aungstrom'
          write(0,*)'This situation is unlikely, reduce ncellz value.'
       endif
       ierr=1 !call mpiend
    endif
    if(ierr.ne.0) call modylas_abort()
  end subroutine check_cutofflength

  logical function is_cell_for_me(icx, icy, icz)
    use margin_sizes
    implicit none
    integer(4), intent(in) :: icx, icy, icz
    integer(4) :: ndirect2,mrcsafe

    if (icx .ge. 2*ncellx .or. icx .le. -2*ncellx .or. &
         &        icy .ge. 2*ncelly .or. icy .le. -2*ncelly .or. &
         &        icz .ge. 2*ncellz .or. icz .le. -2*ncellz) then
       is_cell_for_me = .false.
       return
    endif

    ndirect2 = ndirect + ndcellmargin
    mrcsafe = min(nimage-1, ndirect2)

    is_cell_for_me = (icx >= (ixmin-mrcsafe)) .and. (icx <= (ixmax+mrcsafe)) .and. &
         (icy >= (iymin-mrcsafe)) .and. (icy <= (iymax+mrcsafe)) .and. &
         (icz >= (izmin-mrcsafe)) .and. (icz <= (izmax+mrcsafe))
  end function is_cell_for_me
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to check whether process number is valid or not
!! \author Yoshimichi Andoh
!<
  subroutine check_parallel_condition
!-----------------------------------------------------------------------
    use mpi_3d_grid, only : npx, npy, npz
    implicit none
    include 'mpif.h'
    integer(4) :: ierr

!### checm nprocs ###!
    if(nprocs .lt. 8)then
      if(myrank==0)then
         write(0,*) 'ERROR: nprocs must be greater or equal than 8.'
      endif
      call mpi_barrier(mpi_comm_world, ierr)
      call mpiend
    endif

!### checm nprocs ###!
    if(mod(nprocs,2).ne.0)then
       if(mod(nprocs,3).ne.0)then
          if(myrank==0)then
             write(0,*) 'ERROR: nprocs is not equal to 2 and 3 powers.'
          endif
          call mpi_barrier(mpi_comm_world, ierr)
          call mpiend
       endif
    endif

    if(npx==1 .or. npy==1 .or. npz==1)then
#ifndef ONEPROC_AXIS
       if(myrank==0)then
          write(0,*) 'ERROR: one MPI proc per axis (npx, npy, or npz=1)' &
               &  ,' should be accompanied with -DONEPROC_AXIS in Makefile.'
          write(0,*) 'nprocs              =', nprocs
          write(0,*) 'npx   ,npy   ,npz   =', npx,npy,npz
          write(0,*) 'ncellx,ncelly,ncellz=', ncellx,ncelly,ncellz
       endif
       call modylas_abort
#else /** defined ONEPROC_AXIS **/
#if defined(PRECISION_P2P_SP) || defined(PRECISION_P2P_MIX)
       if(myrank==0)then
          write(0,*) 'ERROR: modylas compiled with -DPRECISION_P2P_SP or -DPRECISION_P2P_MIX ', &
                     'does not run, when npx, npy, or npz = 1.'
          write(0,*) 'nprocs              =', nprocs
          write(0,*) 'npx   ,npy   ,npz   =', npx,npy,npz
       endif
       call modylas_abort
#endif
#endif
    endif

  end subroutine check_parallel_condition

end module domain
