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
!! \brief  Module and subroutines which relate to MPI communications 
!!         of multipoles in the MTD method.
!<
!! Ref: Y.Andoh, S.Ichikawa, T Sakashita, N Yoshii, S Okazaki,
!!      J. Comput. Chem., 42, 1073-1087 (2021).
!----------------------------------------------------------------------
!>
!! \brief  Module which relate to MPI communications of multipoles.
!! \author Yoshimichi Andoh, Shin-ichi Ichikawa
!<
!---------------------------------------------------------------------
module comm_fmm_mod
  use precision
  use mpi_tool
#ifdef PROFILE_COMM
  use profile_comm
#endif
  implicit none
  type :: local_superposition_cell_type
     integer(4), allocatable :: iarray(:,:)
  end type local_superposition_cell_type
  type(local_superposition_cell_type), allocatable :: sclist_xps(:)
  type(local_superposition_cell_type), allocatable :: sclist_xms(:)
  type(local_superposition_cell_type), allocatable :: sclist_xpr(:)
  type(local_superposition_cell_type), allocatable :: sclist_xmr(:)
  type(local_superposition_cell_type), allocatable :: sclist_yps(:)
  type(local_superposition_cell_type), allocatable :: sclist_yms(:)
  type(local_superposition_cell_type), allocatable :: sclist_ypr(:)
  type(local_superposition_cell_type), allocatable :: sclist_ymr(:)
  type(local_superposition_cell_type), allocatable :: sclist_zps(:)
  type(local_superposition_cell_type), allocatable :: sclist_zms(:)
  type(local_superposition_cell_type), allocatable :: sclist_zpr(:)
  type(local_superposition_cell_type), allocatable :: sclist_zmr(:)
#ifdef FJ_RDMA
  type :: address_for_llrdma_type
     integer(4), allocatable :: iaddr(:)
  end type address_for_llrdma_type
  type(address_for_llrdma_type), allocatable :: icbp0_test(:)
  type(address_for_llrdma_type), allocatable :: icbp1_test(:)
  type(address_for_llrdma_type), allocatable :: icbm0_test(:)
  type(address_for_llrdma_type), allocatable :: icbm1_test(:)
!Y+,Y-
  type(address_for_llrdma_type), allocatable :: icbp0_y(:)
  type(address_for_llrdma_type), allocatable :: icbp1_y(:)
  type(address_for_llrdma_type), allocatable :: icbm0_y(:)
  type(address_for_llrdma_type), allocatable :: icbm1_y(:)
!X+,X-
  type(address_for_llrdma_type), allocatable :: icbp0_x(:)
  type(address_for_llrdma_type), allocatable :: icbp1_x(:)
  type(address_for_llrdma_type), allocatable :: icbm0_x(:)
  type(address_for_llrdma_type), allocatable :: icbm1_x(:)
!
  integer(4) nitr_y, nitr_x
#endif /** FJ_RDMA **/

contains

!---------------------------------------------------------------------
!>
!! \brief  Subroutine to initialize arrays used in MPI communications.
!! \author Yoshimichi Andoh, Shin-ichi Ichikawa
!<
  subroutine init_comm_fmm_local_super( &
       &                      scl_zps, scl_zms, scl_zpr, scl_zmr, &
       &                      scl_yps, scl_yms, scl_ypr, scl_ymr, &
       &                      scl_xps, scl_xms, scl_xpr, scl_xmr, &
       &                      nsczdiv, nscydiv, nscxdiv, &
       &                      mbd_z, mbd_y, mbd_x)

!.... Initialization of cell list for thread in Superposition      ....
!.... communication of Multipole moment after M2L.                 ....
    use mpi_3d_grid

    implicit none

    integer(kind=wip), intent(in) :: nscxdiv, nscydiv, nsczdiv
    integer(kind=wip), intent(in) :: mbd_z, mbd_y, mbd_x

    integer(kind=wip) scl_xps(3,0:(mbd_y*2+nscydiv)*(mbd_z*2+nsczdiv)*mbd_x)
    integer(kind=wip) scl_xpr(3,0:(mbd_y*2+nscydiv)*(mbd_z*2+nsczdiv)*mbd_x)
    integer(kind=wip) scl_xms(3,0:(mbd_y*2+nscydiv)*(mbd_z*2+nsczdiv)*mbd_x)
    integer(kind=wip) scl_xmr(3,0:(mbd_y*2+nscydiv)*(mbd_z*2+nsczdiv)*mbd_x)
    integer(kind=wip) scl_yps(3,0:(mbd_z*2+nsczdiv)*(mbd_x*2+nscxdiv)*mbd_y)
    integer(kind=wip) scl_ypr(3,0:(mbd_z*2+nsczdiv)*(mbd_x*2+nscxdiv)*mbd_y)
    integer(kind=wip) scl_yms(3,0:(mbd_z*2+nsczdiv)*(mbd_x*2+nscxdiv)*mbd_y)
    integer(kind=wip) scl_ymr(3,0:(mbd_z*2+nsczdiv)*(mbd_x*2+nscxdiv)*mbd_y)
    integer(kind=wip) scl_zps(3,0:(mbd_x*2+nscxdiv)*(mbd_y*2+nscydiv)*mbd_z)
    integer(kind=wip) scl_zpr(3,0:(mbd_x*2+nscxdiv)*(mbd_y*2+nscydiv)*mbd_z)
    integer(kind=wip) scl_zms(3,0:(mbd_x*2+nscxdiv)*(mbd_y*2+nscydiv)*mbd_z)
    integer(kind=wip) scl_zmr(3,0:(mbd_x*2+nscxdiv)*(mbd_y*2+nscydiv)*mbd_z)

    integer(kind=wip) icz, icy, icx
    integer(kind=wip) iczp0 , iczp1
    integer(kind=wip) iczbp0, iczbp1
    integer(kind=wip) iczm0 , iczm1
    integer(kind=wip) iczbm0, iczbm1
    integer(kind=wip) icyp0 , icyp1
    integer(kind=wip) icybp0, icybp1
    integer(kind=wip) icym0 , icym1
    integer(kind=wip) icybm0, icybm1
    integer(kind=wip) icxp0 , icxp1
    integer(kind=wip) icxbp0, icxbp1
    integer(kind=wip) icxm0 , icxm1
    integer(kind=wip) icxbm0, icxbm1
    integer(kind=wip) ncc
    integer(kind=wip) nccps, nccpr
    integer(kind=wip) nccms, nccmr

! +-X send cell list.
    icxp0 = nscxdiv + 1
    icxp1 = nscxdiv + mbd_x
    icxbp0 = 1
    ! plus direction boundary in source rank is same as this rank for nscxdiv>1.
    icxbp1 = mbd_x
    if(nscxdiv == 1) then
      icxp0 = nscxdiv + mbd_x
      icxp1 = nscxdiv + mbd_x
      icxbp0 = 1
      icxbp1 = 1
    endif
    nccps = 0
    do icx = icxp0, icxp1
      do icy = 1 - mbd_y, nscydiv + mbd_y
        if(nscydiv == 1 .and. icy /= 1 - mbd_y .and. icy /= nscydiv + mbd_y) cycle
        if(nscydiv <= mbd_y .and. 1 <= icy .and. icy <= nscydiv) cycle
        do icz = 1 - mbd_z, nsczdiv + mbd_z
          if(nsczdiv == 1 .and. icz /= 1 - mbd_z .and. icz /= nsczdiv + mbd_z) cycle
          if(nsczdiv <= mbd_z .and. 1 <= icz .and. icz <= nsczdiv) cycle
          nccps = nccps + 1
          scl_xps(1,nccps) = icx
          scl_xps(2,nccps) = icy
          scl_xps(3,nccps) = icz
        end do
      end do
    end do
    scl_xps(1,0) = nccps

    icxm0 = 1 - mbd_x
    icxm1 = 0
    ! minus direction boundary in source rank is same as this rank for nscxdiv>1.
    icxbm0 = nscxdiv - mbd_x + 1
    icxbm1 = nscxdiv
    if(nscxdiv == 1) then 
      icxm0 = 1 - mbd_x
      icxm1 = 1 - mbd_x
      icxbm0 = nscxdiv
      icxbm1 = nscxdiv
    endif
    nccms = 0
    do icx = icxm0, icxm1
      do icy = 1 - mbd_y, nscydiv + mbd_y
        if(nscydiv == 1 .and. icy /= 1 - mbd_y .and. icy /= nscydiv + mbd_y) cycle
        if(nscydiv <= mbd_y .and. 1 <= icy .and. icy <= nscydiv) cycle
        do icz = 1 - mbd_z, nsczdiv + mbd_z
          if(nsczdiv == 1 .and. icz /= 1 - mbd_z .and. icz /= nsczdiv + mbd_z) cycle
          if(nsczdiv <= mbd_z .and. 1 <= icz .and. icz <= nsczdiv) cycle
          nccms = nccms + 1
          scl_xms(1,nccms) = icx
          scl_xms(2,nccms) = icy
          scl_xms(3,nccms) = icz
        end do
      end do
    end do
    scl_xms(1,0) = nccms
! +-X superposition cell list.
    nccpr = 0
    do icx = icxbp0, icxbp1
      do icy = 1 - mbd_y, nscydiv + mbd_y
        if(nscydiv == 1 .and. icy /= 1 - mbd_y .and. icy /= nscydiv + mbd_y) cycle
        if(nscydiv <= mbd_y .and. 1 <= icy .and. icy <= nscydiv) cycle
        do icz = 1 - mbd_z, nsczdiv + mbd_z
          if(nsczdiv == 1 .and. icz /= 1 - mbd_z .and. icz /= nsczdiv + mbd_z) cycle
          if(nsczdiv <= mbd_z .and. 1 <= icz .and. icz <= nsczdiv) cycle
          nccpr = nccpr + 1
          scl_xpr(1,nccpr) = icx
          scl_xpr(2,nccpr) = icy
          scl_xpr(3,nccpr) = icz
        end do
      end do
    end do
    scl_xpr(1,0) = nccpr
    nccmr = 0
    do icx = icxbm0, icxbm1
      do icy = 1 - mbd_y, nscydiv + mbd_y
        if(nscydiv == 1 .and. icy /= 1 - mbd_y .and. icy /= nscydiv + mbd_y) cycle
        if(nscydiv <= mbd_y .and. 1 <= icy .and. icy <= nscydiv) cycle
        do icz = 1 - mbd_z, nsczdiv + mbd_z
          if(nsczdiv == 1 .and. icz /= 1 - mbd_z .and. icz /= nsczdiv + mbd_z) cycle
          if(nsczdiv <= mbd_z .and. 1 <= icz .and. icz <= nsczdiv) cycle
          nccmr = nccmr + 1
          scl_xmr(1,nccmr) = icx
          scl_xmr(2,nccmr) = icy
          scl_xmr(3,nccmr) = icz
        end do
      end do
    end do
    scl_xmr(1,0) = nccmr

#ifdef DBG_COM_LS
    write(myrank+NBFU,*) "+X         : ", &
         &                    "icxp0,icxp1,icxbp0,icxbp1: ", &
         &                     icxp0,icxp1,icxbp0,icxbp1
    write(myrank+NBFU,*) "+X         : ", &
         &            "nccps, nccpr: ", &
         &             nccps, nccpr
    write(myrank+NBFU,*) "-X         : ", &
         &                    "icxm0,icxm1,icxbm0,icxbm1: ", &
         &                     icxm0,icxm1,icxbm0,icxbm1
    write(myrank+NBFU,*) "-X         : ", &
         &            "nccms, nccmr: ", &
         &             nccms, nccmr
#endif

! +-Y send cell list.
    icyp0 = nscydiv + 1
    icyp1 = nscydiv + mbd_y
    icybp0 = 1
    ! plus direction boundary in source rank is same as this rank for nscydiv>1.
    icybp1 = mbd_y
    if(nscydiv == 1) then
      icyp0 = nscydiv + mbd_y
      icyp1 = nscydiv + mbd_y
      icybp0 = 1
      icybp1 = 1
    endif
    nccps = 0
    do icx = 1, nscxdiv
      do icy = icyp0, icyp1
        do icz = 1 - mbd_z, nsczdiv + mbd_z
          if(nsczdiv == 1 .and. icz /= 1 - mbd_z .and. icz /= nsczdiv + mbd_z) cycle
          if(nsczdiv <= mbd_z .and. 1 <= icz .and. icz <= nsczdiv) cycle
          nccps = nccps + 1
          scl_yps(1,nccps) = icx
          scl_yps(2,nccps) = icy
          scl_yps(3,nccps) = icz
        end do
      end do
    end do
    scl_yps(1,0) = nccps

    icym0 = 1 - mbd_y
    icym1 = 0
    ! minus direction boundary in source rank is same as this rank for nscydiv>1.
    icybm0 = nscydiv - mbd_y + 1
    icybm1 = nscydiv
    if(nscydiv == 1) then 
      icym0 = 1 - mbd_y
      icym1 = 1 - mbd_y
      icybm0 = nscydiv
      icybm1 = nscydiv
    endif
    nccms = 0
    do icx = 1, nscxdiv
      do icy = icym0, icym1
        do icz = 1 - mbd_z, nsczdiv + mbd_z
          if(nsczdiv == 1 .and. icz /= 1 - mbd_z .and. icz /= nsczdiv + mbd_z) cycle
          if(nsczdiv <= mbd_z .and. 1 <= icz .and. icz <= nsczdiv) cycle
          nccms = nccms + 1
          scl_yms(1,nccms) = icx
          scl_yms(2,nccms) = icy
          scl_yms(3,nccms) = icz
        end do
      end do
    end do
    scl_yms(1,0) = nccms

! +-Y superposition cell list.
    nccpr = 0
    do icx = 1, nscxdiv
      do icy = icybp0, icybp1
        do icz = 1 - mbd_z, nsczdiv + mbd_z
          if(nsczdiv == 1 .and. icz /= 1 - mbd_z .and. icz /= nsczdiv + mbd_z) cycle
          if(nsczdiv <= mbd_z .and. 1 <= icz .and. icz <= nsczdiv) cycle
          nccpr = nccpr + 1
          scl_ypr(1,nccpr) = icx
          scl_ypr(2,nccpr) = icy
          scl_ypr(3,nccpr) = icz
        end do
      end do
    end do
    scl_ypr(1,0) = nccpr
    nccmr = 0
    do icx = 1, nscxdiv
      do icy = icybm0, icybm1
        do icz = 1 - mbd_z, nsczdiv + mbd_z
          if(nsczdiv == 1 .and. icz /= 1 - mbd_z .and. icz /= nsczdiv + mbd_z) cycle
          if(nsczdiv <= mbd_z .and. 1 <= icz .and. icz <= nsczdiv) cycle
          nccmr = nccmr + 1
          scl_ymr(1,nccmr) = icx
          scl_ymr(2,nccmr) = icy
          scl_ymr(3,nccmr) = icz
        end do
      end do
    end do
    scl_ymr(1,0) = nccmr

#ifdef DBG_COM_LS
    write(myrank+NBFU,*) "+Y         : ", &
         &                    "icyp0,icyp1,icybp0,icybp1: ", &
         &                     icyp0,icyp1,icybp0,icybp1
    write(myrank+NBFU,*) "+Y         : ", &
         &            "nccps, nccpr: ", &
         &             nccps, nccpr
    write(myrank+NBFU,*) "-Y         : ", &
         &                    "icym0,icym1,icybm0,icybm1: ", &
         &                     icym0,icym1,icybm0,icybm1
    write(myrank+NBFU,*) "-Y         : ", &
         &            "nccps, nccpr: ", &
         &             nccps, nccpr
#endif

! +-Z send cell list.
    iczp0 = nsczdiv + 1
    iczp1 = nsczdiv + mbd_z
    iczbp0 = 1
    ! plus direction boundary in source rank is same as this rank for nsczdiv>1.
    iczbp1 = mbd_z
    if(nsczdiv == 1) then
      iczp0 = nsczdiv + mbd_z
      iczp1 = nsczdiv + mbd_z
      iczbp0 = 1
      iczbp1 = 1
    endif
    nccps = 0
    do icx = 1, nscxdiv
      do icy = 1, nscydiv
        do icz = iczp0, iczp1
          nccps = nccps + 1
          scl_zps(1,nccps) = icx
          scl_zps(2,nccps) = icy
          scl_zps(3,nccps) = icz
        end do
      end do
    end do
    scl_zps(1,0) = nccps

    iczm0 = 1 - mbd_z
    iczm1 = 0
    ! minus direction boundary in source rank is same as this rank for nsczdiv>1.
    iczbm0 = nsczdiv - mbd_z + 1
    iczbm1 = nsczdiv
    if(nsczdiv == 1) then 
      iczm0 = 1 - mbd_z
      iczm1 = 1 - mbd_z
      iczbm0 = nsczdiv
      iczbm1 = nsczdiv
    endif
    nccms = 0
    do icx = 1, nscxdiv
      do icy = 1, nscydiv
        do icz = iczm0, iczm1
          nccms = nccms + 1
          scl_zms(1,nccms) = icx
          scl_zms(2,nccms) = icy
          scl_zms(3,nccms) = icz
        end do
      end do
    end do
    scl_zms(1,0) = nccms

! +-Z superposition cell list.
    nccpr = 0
    do icx = 1, nscxdiv
      do icy = 1, nscydiv
        do icz = iczbp0, iczbp1
          nccpr = nccpr + 1
          scl_zpr(1,nccpr) = icx
          scl_zpr(2,nccpr) = icy
          scl_zpr(3,nccpr) = icz
        end do
      end do
    end do
    scl_zpr(1,0) = nccpr

    nccmr = 0
    do icx = 1, nscxdiv
      do icy = 1, nscydiv
        do icz = iczbm0, iczbm1
          nccmr = nccmr + 1
          scl_zmr(1,nccmr) = icx
          scl_zmr(2,nccmr) = icy
          scl_zmr(3,nccmr) = icz
        end do
      end do
    end do
    scl_zmr(1,0) = nccmr

#ifdef DBG_COM_LS
    write(myrank+NBFU,*) "+Z         : ", &
         &            "iczp0,iczp1,iczbp0,iczbp1: ", &
         &             iczp0,iczp1,iczbp0,iczbp1
    write(myrank+NBFU,*) "+Z        : ", &
         &            "nccps, nccpr: ", &
         &             nccps, nccpr
    write(myrank+NBFU,*) "-Z        : ", &
         &            "iczm0,iczm1,iczbm0,iczbm1: ", &
         &             iczm0,iczm1,iczbm0,iczbm1
    write(myrank+NBFU,*) "-Z        : ", &
         &            "nccms, nccmr: ", &
         &             nccms, nccmr
#endif

  end subroutine init_comm_fmm_local_super
!---------------------------------------------------------------------
#ifdef FJ_RDMA
!>
!! \brief  Subroutine to initialize arrays used in RDMA communications.
!! \author Yoshimichi Andoh, Shin-ichi Ichikawa
!<
  subroutine init_comm_fmm_local_multi_dr( &
                      il,nlz,nly,nlx, &
                      nscellz,nscelly,nscellx,nsczdiv,nscydiv,nscxdiv)
    use domain, only : ncellx, ncelly, ncellz
    use dr_cntl, only : get_halo_number
    use mpi_tool, only : myrank
    use mpi_3d_grid, only : ipz,ipy,ipx,npz,npy,npx
    implicit none
    integer(4),intent(in) :: il
    integer(4),intent(in) :: nlz,nly,nlx
    integer(4),intent(in) :: nsczdiv,nscydiv,nscxdiv
    integer(4),intent(in) :: nscellz,nscelly,nscellx
    integer(4) nbd_m,nbd_p,idum
    integer(4) nb_destm,nb_destp,nitr
    integer(4) np_supercelly,np_supercellx
    integer(4) ipy_pdest,ipy_psrc,ipy_mdest,ipy_msrc
    integer(4) ipx_pdest,ipx_psrc,ipx_mdest,ipx_msrc
! for RDMA
    integer(4) nsy,nsx
! for RDMA

    ! LLmoment +Y
    np_supercelly = (npy - 1) / nscelly + 1
    ipy_pdest = ipz*npy*npx + mod(ipy + np_supercelly - 1/npy + npy, npy)*npx + ipx
    ipy_psrc = ipz*npy*npx + mod(ipy - np_supercelly + 1/npy + npy, npy)*npx + ipx
    ! LLmoment -Y
    ipy_mdest = ipz*npy*npx + mod(ipy - np_supercelly + 1/npy + npy, npy)*npx + ipx
    ipy_msrc = ipz*npy*npx + mod(ipy + np_supercelly - 1/npy + npy, npy)*npx + ipx

    ! LLmoment +X
    np_supercellx = (npx - 1) / nscellx + 1
    ipx_pdest = ipz*npy*npx + ipy*npx + mod(ipx + np_supercellx - 1/npx + npx, npx)
    ipx_psrc = ipz*npy*npx + ipy*npx + mod(ipx - np_supercellx + 1/npx + npx, npx)
    ! LLmoment -X
    ipx_mdest = ipz*npy*npx + ipy*npx + mod(ipx - np_supercellx + 1/npx + npx, npx)
    ipx_msrc = ipz*npy*npx + ipy*npx + mod(ipx + np_supercellx - 1/npx + npx, npx)

    if(il==0)then
      nsy=1
    else
      nsy=max(npy/nscelly,1)
    endif
!
! information about ipy+1 (Y+)
!
    call get_halo_number(ipy+nsy, &
     &              npy,nly,nscelly, &
     &              nbd_m,nbd_p,idum)
    call preprocess_comm(ipy+nsy,nbd_m,nbd_p, &
     &                   npy,nly,nscelly,nscydiv, &
     &                   nb_destm,nb_destp,nitr_y)    
    allocate(icbp0_y(il)%iaddr(nitr_y))
    allocate(icbp1_y(il)%iaddr(nitr_y))
    call get_icbp01(nbd_m,nbd_p,nitr_y, &
     &              nscydiv,nb_destp, &
     &              icbp0_y(il)%iaddr,icbp1_y(il)%iaddr)
!
! information about ipy-1 (Y-)
!
    call get_halo_number(ipy-nsy, &
     &              npy,nly,nscelly, &
     &              nbd_m,nbd_p,idum)
    call preprocess_comm(ipy-nsy,nbd_m,nbd_p, &
     &                   npy,nly,nscelly,nscydiv, &
     &                   nb_destm,nb_destp,nitr_y)    
    allocate(icbm0_y(il)%iaddr(nitr_y))
    allocate(icbm1_y(il)%iaddr(nitr_y))
    call get_icbm01(nbd_m,nbd_p,nitr_y, &
     &               nscydiv,nb_destm, &
     &               icbm0_y(il)%iaddr,icbm1_y(il)%iaddr)
    if(il==0)then
      nsx=1
    else
      nsx=max(npx/nscellx,1)
    endif
!
! information about ipx+1 (X+)
!
    call get_halo_number(ipx+nsx, &
     &              npx,nlx,nscellx, &
     &              nbd_m,nbd_p,idum)
    call preprocess_comm(ipx+nsx,nbd_m,nbd_p, &
     &                   npx,nlx,nscellx,nscxdiv, &
     &                   nb_destm,nb_destp,nitr_x)    
    allocate(icbp0_x(il)%iaddr(nitr_x))
    allocate(icbp1_x(il)%iaddr(nitr_x))
    call get_icbp01(nbd_m,nbd_p,nitr_x, &
     &               nscxdiv,nb_destp, &
     &               icbp0_x(il)%iaddr,icbp1_x(il)%iaddr)
!
! information about ipx-1 (X-)
!
    call get_halo_number(ipx-nsx, &
     &              npx,nlx,nscellx, &
     &              nbd_m,nbd_p,idum)
    call preprocess_comm(ipx-nsx,nbd_m,nbd_p, &
     &                   npx,nlx,nscellx,nscxdiv, &
     &                   nb_destm,nb_destp,nitr_x)    
    allocate(icbm0_x(il)%iaddr(nitr_x))
    allocate(icbm1_x(il)%iaddr(nitr_x))
    call get_icbm01(nbd_m,nbd_p,nitr_x, &
     &               nscxdiv,nb_destm, &
     &               icbm0_x(il)%iaddr,icbm1_x(il)%iaddr)
  end subroutine init_comm_fmm_local_multi_dr
#endif /** FJ_RDMA **/
!---------------------------------------------------------------------
!>
!! \brief  Subroutine to perform MPI communication to superpose
!!         partial local expansion coefficient after the M2L.
!! \author Yoshimichi Andoh, Shin-ichi Ichikawa
!<
  subroutine comm_fmm_local_super(nl,mylm0, nsbsize, wl, &
       &                      scl_zps, scl_zms, scl_zpr, scl_zmr, &
       &                      scl_yps, scl_yms, scl_ypr, scl_ymr, &
       &                      scl_xps, scl_xms, scl_xpr, scl_xmr, &
       &                      nsczdiv,nscydiv,nscxdiv, &
       &                      mcell_size_z,mcell_size_y,mcell_size_x, &
       &                      mbd_z, mbd_y, mbd_x)
!.... Superposition communication of Multipole moment after M2L.   ....
    use mpi_3d_grid
    use domain, only : ncellx, ncelly, ncellz
    use dr_cntl, only : mbd, mtbd
    use md_condition, only : mdstep
    implicit none
    integer(kind=wip) nl
    integer(kind=wip) mylm,mylm0
    integer(kind=wip) nsbsize
    integer(kind=wip) nscxdiv, nscydiv, nsczdiv
    integer(kind=wip) mcell_size_x, mcell_size_y, mcell_size_z
    integer(kind=wip), intent(in) :: mbd_z, mbd_y, mbd_x
    real(kind=wvap) wl(mylm0,2, 1 - mbd_z:nsczdiv + mbd_z, &
         &                      1 - mbd_y:nscydiv + mbd_y, &
         &                      1 - mbd_x:nscxdiv + mbd_x)
#ifndef FJ_RDMA
    real(kind=wvap) ccbufp (2*mylm0*nsbsize)
    real(kind=wvap) rccbufp(2*mylm0*nsbsize)
    real(kind=wvap) ccbufm (2*mylm0*nsbsize)
    real(kind=wvap) rccbufm(2*mylm0*nsbsize)
#endif
    integer(kind=wip) scl_xps(3,0:(mbd_y*2+nscydiv)*(mbd_z*2+nsczdiv)*mbd_x)
    integer(kind=wip) scl_xpr(3,0:(mbd_y*2+nscydiv)*(mbd_z*2+nsczdiv)*mbd_x)
    integer(kind=wip) scl_xms(3,0:(mbd_y*2+nscydiv)*(mbd_z*2+nsczdiv)*mbd_x)
    integer(kind=wip) scl_xmr(3,0:(mbd_y*2+nscydiv)*(mbd_z*2+nsczdiv)*mbd_x)
    integer(kind=wip) scl_yps(3,0:(mbd_z*2+nsczdiv)*(mbd_x*2+nscxdiv)*mbd_y)
    integer(kind=wip) scl_ypr(3,0:(mbd_z*2+nsczdiv)*(mbd_x*2+nscxdiv)*mbd_y)
    integer(kind=wip) scl_yms(3,0:(mbd_z*2+nsczdiv)*(mbd_x*2+nscxdiv)*mbd_y)
    integer(kind=wip) scl_ymr(3,0:(mbd_z*2+nsczdiv)*(mbd_x*2+nscxdiv)*mbd_y)
    integer(kind=wip) scl_zps(3,0:(mbd_x*2+nscxdiv)*(mbd_y*2+nscydiv)*mbd_z)
    integer(kind=wip) scl_zpr(3,0:(mbd_x*2+nscxdiv)*(mbd_y*2+nscydiv)*mbd_z)
    integer(kind=wip) scl_zms(3,0:(mbd_x*2+nscxdiv)*(mbd_y*2+nscydiv)*mbd_z)
    integer(kind=wip) scl_zmr(3,0:(mbd_x*2+nscxdiv)*(mbd_y*2+nscydiv)*mbd_z)

    integer(kind=wip) np_supercell, neig_dist
    integer(kind=wip) ipz_pdest, ipy_pdest, ipx_pdest
    integer(kind=wip) ipz_psrc, ipy_psrc, ipx_psrc
    integer(kind=wip) ipz_mdest, ipy_mdest, ipx_mdest
    integer(kind=wip) ipz_msrc, ipy_msrc, ipx_msrc

    integer(kind=wip) icz, icy, icx
    integer(kind=wip) icc
    integer(kind=wip) ncc
    integer(kind=wip) nccps, nccpr
    integer(kind=wip) nccms, nccmr
    integer(kind=wip) m

    integer(kind=wip) istatus(mpi_status_size, 4), ierr
    integer(kind=wip),dimension(4) :: irq
    integer(kind=wip) nrq

!debug
    integer(kind=wip) :: icx0, icy0, icz0
    integer(kind=wip) :: icxyz
    integer(kind=wip) :: lcbx, lcby, lcbz

#ifdef FJ_RDMA
  real(kind=wvap),pointer :: ccbufp (:)
  real(kind=wvap),pointer :: rccbufp(:)
  real(kind=wvap),pointer :: ccbufm (:)
  real(kind=wvap),pointer :: rccbufm(:)
  integer(8) :: heap_laddr
  integer(8), pointer :: heap_raddr(:)
  type(c_ptr) :: cptr
  integer(8)  :: ccbufp_offset
  integer(8)  :: rccbufp_offset
  integer(8)  :: ccbufm_offset
  integer(8)  :: rccbufm_offset
  integer(4) ibds  ! data size for one data-block
  !
#ifdef PRECISION_M2L_SP
  ibds=8   ! SP for wvap
#else /** PRECISION_M2L_MIX, or no preprocessor **/
  ibds=16  ! DP for wvap
#endif
  !
  cptr = rdma_allocate_coarray((mylm0*nsbsize)*ibds, ccbufp_offset)
  call c_f_pointer(cptr, fptr=ccbufp, shape=[mylm0*nsbsize])
  cptr = rdma_allocate_coarray((mylm0*nsbsize)*ibds, rccbufp_offset)
  call c_f_pointer(cptr, fptr=rccbufp, shape=[mylm0*nsbsize])
  cptr = rdma_allocate_coarray((mylm0*nsbsize)*ibds, ccbufm_offset)
  call c_f_pointer(cptr, fptr=ccbufm, shape=[mylm0*nsbsize])
  cptr = rdma_allocate_coarray((mylm0*nsbsize)*ibds, rccbufm_offset)
  call c_f_pointer(cptr, fptr=rccbufm, shape=[mylm0*nsbsize])
  !
  heap_laddr = rdma_get_heap_laddr()
  cptr       = rdma_get_heap_raddr()
  call c_f_pointer(cptr, fptr=heap_raddr, shape=[nprocs])
#endif /** FJ_RDMA **/

    mylm=mylm0*1

    lcbx = nscxdiv !+ mbd_x2
    lcby = nscydiv !+ mbd_y2
    lcbz = nsczdiv !+ mbd_z2

#ifdef COMLS_TIMER
    te = MPI_WTIME()
    tcomls_ent = tcomls_ent + te - ts
    ts = te
#endif

    ! ---- 3D rank order rule. ----
    !     ipx=mod(myrank, npx)
    !     ipy=mod((myrank - ipx) / npx, npy)
    !     ipz=mod((myrank - ipx - ipy*npx) / (npx*npy), npz)

#ifdef DBG_COMLS
    write(myrank+NBFU,*) "npz,npy,npx:",npz,npy,npx
    write(myrank+NBFU,*) "myrank,ipz,ipy,ipx:",myrank,ipz,ipy,ipx
#endif

    !----- common parameters for lower level moment communication. -----
    !      nsczdiv = (ncellz / mcell_size_z - 1) / npz + 1
    !      nscydiv = (ncelly / mcell_size_y - 1) / npy + 1
    !      nscxdiv = (ncellx / mcell_size_x - 1) / npx + 1

    !----- lower level moment superposition communication starts here. -----

    ! LSmoment +X
    np_supercell = (npx* mcell_size_x - 1) / ncellx + 1
    neig_dist = np_supercell
    if(nscxdiv == 1) neig_dist = mbd_x*np_supercell
    ipx_pdest = ipz*npy*npx + ipy*npx + modulo(ipx + neig_dist*(1 - 1/npx), npx)
    ipx_psrc  = ipz*npy*npx + ipy*npx + modulo(ipx - neig_dist*(1 - 1/npx), npx)
    ! LSmoment -X
    ipx_mdest = ipz*npy*npx + ipy*npx + modulo(ipx - neig_dist*(1 - 1/npx), npx)
    ipx_msrc  = ipz*npy*npx + ipy*npx + modulo(ipx + neig_dist*(1 - 1/npx), npx)

#ifdef DBG_COMLS
    write(myrank+NBFU,*) "+X: ipx_pdest,ipx_psrc,np_supercell: ", &
         &                   ipx_pdest,ipx_psrc,np_supercell
    write(myrank+NBFU,*) "+X: neig_dist;",neig_dist
    write(myrank+NBFU,*) "-X: ipx_mdest,ipx_msrc,np_supercell: ", &
         &                   ipx_mdest,ipx_msrc,np_supercell
    write(myrank+NBFU,*) "-X: neig_dist;",neig_dist
#endif

!$omp parallel default(none) &
!$omp& private(icc,icx,icy,icz,ncc,m) &
!$omp& shared(mylm,wl) &
!$omp& shared(scl_xps,ccbufp,scl_xms,ccbufm)
!$omp do
    do icc = 1, scl_xps(1,0)
      icx = scl_xps(1,icc)
      icy = scl_xps(2,icc)
      icz = scl_xps(3,icc)
      ncc = 2*(icc - 1)*mylm
      do m = 1, mylm
        ccbufp(ncc + m       ) = wl(m,1,icz,icy,icx)
        ccbufp(ncc + mylm + m) = wl(m,2,icz,icy,icx)
      end do
    end do
!$omp end do
!$omp do
    do icc = 1, scl_xms(1,0)
      icx = scl_xms(1,icc)
      icy = scl_xms(2,icc)
      icz = scl_xms(3,icc)
      ncc = 2*(icc - 1)*mylm
      do m = 1, mylm
        ccbufm(ncc + m       ) = wl(m,1,icz,icy,icx)
        ccbufm(ncc + mylm + m) = wl(m,2,icz,icy,icx)
      end do
    end do
!$omp end do
!$omp end parallel
    nccps = scl_xps(1,0)*mylm*2
    nccms = scl_xms(1,0)*mylm*2
    nccpr = nccps
    nccmr = nccms

#ifdef COMLS_TIMER
    te = MPI_WTIME()
    tcomls_lcl = tcomls_lcl + te - ts
    call mpi_barrier(mpi_comm_world,ierr)
    ts = MPI_WTIME()
    tcomls_pbar = tcomls_pbar + ts - te
#endif

#ifdef PROFILE_COMM
    call mpi_barrier(mpi_comm_world, ierr)
    call start_profile(comm_fmm_local_super_id)
#endif
#ifdef FJ_RDMA
    call rdma_post(ipx_psrc)
    if (ipx_psrc .ne. ipx_msrc) then
       call rdma_post(ipx_msrc)
       call rdma_wait(ipx_mdest)
    endif
    call rdma_wait(ipx_pdest)
    !
    if (ipx_pdest .ne. ipx_mdest) then
       call rdma_put_post(ipx_pdest, heap_raddr(ipx_pdest+1)+rccbufp_offset, heap_laddr+ccbufp_offset, nccps*ibds)
       call rdma_put_post(ipx_mdest, heap_raddr(ipx_mdest+1)+rccbufm_offset, heap_laddr+ccbufm_offset, nccms*ibds)
       call rdma_wait(ipx_psrc)
       call rdma_wait(ipx_msrc)
    else
       call rdma_put     (ipx_pdest, heap_raddr(ipx_pdest+1)+rccbufp_offset, heap_laddr+ccbufp_offset, nccps*ibds)
       call rdma_put_post(ipx_mdest, heap_raddr(ipx_mdest+1)+rccbufm_offset, heap_laddr+ccbufm_offset, nccms*ibds)
       call rdma_wait(ipx_msrc)
    endif
#else
    call mpi_irecv(rccbufp, nccpr, &
         &             MPI_PREC_MOMVEC_SUPER, ipx_psrc,    1, &
         &             mpi_comm_world, irq(1), ierr)
    call mpi_isend(ccbufp, nccps, &
         &             MPI_PREC_MOMVEC_SUPER, ipx_pdest,   1, &
         &             mpi_comm_world, irq(2), ierr)
    call mpi_irecv(rccbufm, nccmr, &
         &             MPI_PREC_MOMVEC_SUPER, ipx_msrc,    2, &
         &             mpi_comm_world, irq(3), ierr)
    call mpi_isend(ccbufm, nccms, &
         &             MPI_PREC_MOMVEC_SUPER, ipx_mdest,   2, &
         &             mpi_comm_world, irq(4), ierr)
    nrq = 4
    call mpi_waitall(nrq, irq, istatus, ierr)
#endif
#ifdef PROFILE_COMM
    call end_profile(comm_fmm_local_super_id)
#endif

#ifdef COMLS_TIMER
    te = MPI_WTIME()
    tcomls_com = tcomls_com + te - ts
    nccomls_com = nccomls_com + 1
#ifdef COMLS_TIMER_DETAIL
    write(myrank+NBFU,*) "Rank: ", myrank, &
         &                    " TcomLS-X: ", te - ts
#endif
    ts = te
#endif

!$omp parallel default(none) &
!$omp& private(icc,icx,icy,icz,ncc,m) &
!$omp& shared(mylm,wl) &
!$omp& shared(scl_xpr,rccbufp,scl_xmr,rccbufm)
!$omp do
    do icc = 1, scl_xpr(1,0)
      icx = scl_xpr(1,icc)
      icy = scl_xpr(2,icc)
      icz = scl_xpr(3,icc)
      ncc = 2*(icc - 1)*mylm
      do m = 1, mylm
        wl(m,1,icz,icy,icx) = wl(m,1,icz,icy,icx) + rccbufp(ncc + m       )
        wl(m,2,icz,icy,icx) = wl(m,2,icz,icy,icx) + rccbufp(ncc + mylm + m)
      end do
    end do
!$omp end do
!$omp do
    do icc = 1, scl_xmr(1,0)
      icx = scl_xmr(1,icc)
      icy = scl_xmr(2,icc)
      icz = scl_xmr(3,icc)
      ncc = 2*(icc - 1)*mylm
      do m = 1, mylm
        wl(m,1,icz,icy,icx) = wl(m,1,icz,icy,icx) + rccbufm(ncc + m       )
        wl(m,2,icz,icy,icx) = wl(m,2,icz,icy,icx) + rccbufm(ncc + mylm + m)
      end do
    end do
!$omp end do
!$omp end parallel

    ! LSmoment +Y
    np_supercell = (npy * mcell_size_y - 1) / ncelly + 1
    neig_dist = np_supercell
    if(nscydiv == 1) neig_dist = mbd_y*np_supercell
    ipy_pdest = ipz*npy*npx + modulo(ipy + neig_dist*(1 - 1/npy), npy)*npx + ipx
    ipy_psrc  = ipz*npy*npx + modulo(ipy - neig_dist*(1 - 1/npy), npy)*npx + ipx
    ! LSmoment -Y
    ipy_mdest = ipz*npy*npx + modulo(ipy - neig_dist*(1 - 1/npy), npy)*npx + ipx
    ipy_msrc  = ipz*npy*npx + modulo(ipy + neig_dist*(1 - 1/npy), npy)*npx + ipx

#ifdef DBG_COMLS
    write(myrank+NBFU,*) "+Y: ipy_pdest,ipy_psrc,np_supercell: ", &
         &                   ipy_pdest,ipy_psrc,np_supercell
    write(myrank+NBFU,*) "-Y: ipy_pdest,ipy_psrc,np_supercell: ", &
         &                   ipy_mdest,ipy_msrc,np_supercell
#endif

!$omp parallel default(none) &
!$omp& private(icc,icx,icy,icz,ncc,m) &
!$omp& shared(mylm,wl) &
!$omp& shared(scl_yps,ccbufp,scl_yms,ccbufm)
!$omp do
    do icc = 1, scl_yps(1,0)
      icx = scl_yps(1,icc)
      icy = scl_yps(2,icc)
      icz = scl_yps(3,icc)
      ncc = 2*(icc - 1)*mylm
      do m = 1, mylm
        ccbufp(ncc + m       ) = wl(m,1,icz,icy,icx)
        ccbufp(ncc + mylm + m) = wl(m,2,icz,icy,icx)
      end do
    end do
!$omp end do
!$omp do
    do icc = 1, scl_yms(1,0)
      icx = scl_yms(1,icc)
      icy = scl_yms(2,icc)
      icz = scl_yms(3,icc)
      ncc = 2*(icc - 1)*mylm
      do m = 1, mylm
        ccbufm(ncc + m       ) = wl(m,1,icz,icy,icx)
        ccbufm(ncc + mylm + m) = wl(m,2,icz,icy,icx)
      end do
    end do
!$omp end do
!$omp end parallel
    nccps = scl_yps(1,0)*mylm*2
    nccms = scl_yms(1,0)*mylm*2
    nccpr = nccps
    nccmr = nccms

#ifdef COMLS_TIMER
    te = MPI_WTIME()
    tcomls_lcl = tcomls_lcl + te - ts
    call mpi_barrier(mpi_comm_world,ierr)
    ts = MPI_WTIME()
    tcomls_pbar = tcomls_pbar + ts - te
#endif

#ifdef PROFILE_COMM
    call mpi_barrier(mpi_comm_world, ierr)
    call start_profile(comm_fmm_local_super_id)
#endif
#ifdef FJ_RDMA
    call rdma_post(ipy_psrc)
    if (ipy_psrc .ne. ipy_msrc) then
       call rdma_post(ipy_msrc)
       call rdma_wait(ipy_mdest)
    endif
    call rdma_wait(ipy_pdest)
    !
    if (ipy_pdest .ne. ipy_mdest) then
       call rdma_put_post(ipy_pdest, heap_raddr(ipy_pdest+1)+rccbufp_offset, heap_laddr+ccbufp_offset, nccps*ibds)
       call rdma_put_post(ipy_mdest, heap_raddr(ipy_mdest+1)+rccbufm_offset, heap_laddr+ccbufm_offset, nccms*ibds)
       call rdma_wait(ipy_psrc)
       call rdma_wait(ipy_msrc)
    else
       call rdma_put     (ipy_pdest, heap_raddr(ipy_pdest+1)+rccbufp_offset, heap_laddr+ccbufp_offset, nccps*ibds)
       call rdma_put_post(ipy_mdest, heap_raddr(ipy_mdest+1)+rccbufm_offset, heap_laddr+ccbufm_offset, nccms*ibds)
       call rdma_wait(ipy_msrc)
    endif
#else
    call mpi_irecv(rccbufp, nccpr, &
         &              MPI_PREC_MOMVEC_SUPER, ipy_psrc,   1, &
         &              mpi_comm_world, irq(1), ierr)
    call mpi_isend(ccbufp, nccps, &
         &              MPI_PREC_MOMVEC_SUPER, ipy_pdest,  1, &
         &              mpi_comm_world, irq(2), ierr)
    call mpi_irecv(rccbufm, nccmr, &
         &              MPI_PREC_MOMVEC_SUPER, ipy_msrc,   2, &
         &              mpi_comm_world, irq(3), ierr)
    call mpi_isend(ccbufm, nccms, &
         &              MPI_PREC_MOMVEC_SUPER, ipy_mdest,  2, &
         &              mpi_comm_world, irq(4), ierr)
    nrq = 4
    call mpi_waitall(nrq, irq, istatus, ierr)
#endif
#ifdef PROFILE_COMM
    call end_profile(comm_fmm_local_super_id)
#endif

#ifdef COMLS_TIMER
    te = MPI_WTIME()
    tcomls_com = tcomls_com + te - ts
    nccomls_com = nccomls_com + 1
#ifdef COMLS_TIMER_DETAIL
    write(myrank+NBFU,*) "Rank: ", myrank, &
         &                    " TcomLS-Y: ", te - ts
#endif
    ts = te
#endif

!$omp parallel default(none) &
!$omp& private(icc,icx,icy,icz,ncc,m) &
!$omp& shared(mylm,wl) &
!$omp& shared(scl_ypr,rccbufp,scl_ymr,rccbufm)
!$omp do
    do icc = 1, scl_ypr(1,0)
      icx = scl_ypr(1,icc)
      icy = scl_ypr(2,icc)
      icz = scl_ypr(3,icc)
      ncc = 2*(icc - 1)*mylm
      do m = 1, mylm
        wl(m,1,icz,icy,icx) = wl(m,1,icz,icy,icx) + rccbufp(ncc + m       )
        wl(m,2,icz,icy,icx) = wl(m,2,icz,icy,icx) + rccbufp(ncc + mylm + m)
      end do
    end do
!$omp end do
!$omp do
    do icc = 1, scl_ymr(1,0)
      icx = scl_ymr(1,icc)
      icy = scl_ymr(2,icc)
      icz = scl_ymr(3,icc)
      ncc = 2*(icc - 1)*mylm
      do m = 1, mylm
        wl(m,1,icz,icy,icx) = wl(m,1,icz,icy,icx) + rccbufm(ncc + m       )
        wl(m,2,icz,icy,icx) = wl(m,2,icz,icy,icx) + rccbufm(ncc + mylm + m)
      end do
    end do
!$omp end do
!$omp end parallel

    ! LSmoment +Z
    np_supercell = (npz * mcell_size_z - 1) / ncellz + 1
    neig_dist = np_supercell
    if(nsczdiv == 1) neig_dist = mbd_z*np_supercell
    ipz_pdest = modulo(ipz + neig_dist*(1 - 1/npz), npz)*npy*npx + ipy*npx + ipx
    ipz_psrc  = modulo(ipz - neig_dist*(1 - 1/npz), npz)*npy*npx + ipy*npx + ipx
    ! LSmoment -Z
    ipz_mdest = modulo(ipz - neig_dist*(1 - 1/npz), npz)*npy*npx + ipy*npx + ipx
    ipz_msrc  = modulo(ipz + neig_dist*(1 - 1/npz), npz)*npy*npx + ipy*npx + ipx

#ifdef DBG_COMLS
    write(myrank+NBFU,*) "+Z: ipz_pdest,ipz_psrc,np_supercell: ", &
         &                   ipz_pdest,ipz_psrc,np_supercell
    write(myrank+NBFU,*) "-Z: ipz_mdest,ipz_msrc,np_supercell: ", &
         &                   ipz_mdest,ipz_msrc,np_supercell
#endif

!$omp parallel default(none) &
!$omp& private(icc,icx,icy,icz,ncc,m) &
!$omp& shared(mylm,wl) &
!$omp& shared(scl_zps,ccbufp,scl_zms,ccbufm)
!$omp do
    do icc = 1, scl_zps(1,0)
      icx = scl_zps(1,icc)
      icy = scl_zps(2,icc)
      icz = scl_zps(3,icc)
      ncc = 2*(icc - 1)*mylm
      do m = 1, mylm
        ccbufp(ncc + m       ) = wl(m,1,icz,icy,icx)
        ccbufp(ncc + mylm + m) = wl(m,2,icz,icy,icx)
      end do
    end do
!$omp end do
!$omp do
    do icc = 1, scl_zms(1,0)
      icx = scl_zms(1,icc)
      icy = scl_zms(2,icc)
      icz = scl_zms(3,icc)
      ncc = 2*(icc - 1)*mylm
      do m = 1, mylm
        ccbufm(ncc + m       ) = wl(m,1,icz,icy,icx)
        ccbufm(ncc + mylm + m) = wl(m,2,icz,icy,icx)
      end do
    end do
!$omp end do
!$omp end parallel
    nccps = scl_zps(1,0)*mylm*2
    nccms = scl_zms(1,0)*mylm*2
    nccpr = nccps
    nccmr = nccms

#ifdef COMLS_TIMER
    te = MPI_WTIME()
    tcomls_lcl = tcomls_lcl + te - ts
    call mpi_barrier(mpi_comm_world,ierr)
    ts = MPI_WTIME()
    tcomls_pbar = tcomls_pbar + ts - te
#endif

#ifdef PROFILE_COMM
    call mpi_barrier(mpi_comm_world, ierr)
    call start_profile(comm_fmm_local_super_id)
#endif
#ifdef FJ_RDMA
    call rdma_post(ipz_psrc)
    if (ipz_psrc .ne. ipz_msrc) then
       call rdma_post(ipz_msrc)
       call rdma_wait(ipz_mdest)
    endif
    call rdma_wait(ipz_pdest)
    !
    if (ipz_pdest .ne. ipz_mdest) then
       call rdma_put_post(ipz_pdest, heap_raddr(ipz_pdest+1)+rccbufp_offset, heap_laddr+ccbufp_offset, nccps*ibds)
       call rdma_put_post(ipz_mdest, heap_raddr(ipz_mdest+1)+rccbufm_offset, heap_laddr+ccbufm_offset, nccms*ibds)
       call rdma_wait(ipz_psrc)
       call rdma_wait(ipz_msrc)
    else
       call rdma_put     (ipz_pdest, heap_raddr(ipz_pdest+1)+rccbufp_offset, heap_laddr+ccbufp_offset, nccps*ibds)
       call rdma_put_post(ipz_mdest, heap_raddr(ipz_mdest+1)+rccbufm_offset, heap_laddr+ccbufm_offset, nccms*ibds)
       call rdma_wait(ipz_msrc)
    endif
#else
    call mpi_irecv(rccbufp, nccpr, &
         &              MPI_PREC_MOMVEC_SUPER, ipz_psrc,   1, &
         &              mpi_comm_world, irq(1), ierr)
    call mpi_isend(ccbufp, nccps, &
         &              MPI_PREC_MOMVEC_SUPER, ipz_pdest,  1, &
         &              mpi_comm_world, irq(2), ierr)
    call mpi_irecv(rccbufm, nccmr, &
         &              MPI_PREC_MOMVEC_SUPER, ipz_msrc,   2, &
         &              mpi_comm_world, irq(3), ierr)
    call mpi_isend(ccbufm, nccms, &
         &              MPI_PREC_MOMVEC_SUPER, ipz_mdest,  2, &
         &              mpi_comm_world, irq(4), ierr)
    nrq = 4
    call mpi_waitall(nrq, irq, istatus, ierr)
#endif
#ifdef PROFILE_COMM
    call end_profile(comm_fmm_local_super_id)
#endif

#ifdef COMLS_TIMER
    te = MPI_WTIME()
    tcomls_com = tcomls_com + te - ts
    nccomls_com = nccomls_com + 1
#ifdef COMLS_TIMER_DETAIL
    write(myrank+NBFU,*) "Rank: ", myrank, &
         &                    " TcomLS-Z: ", te - ts
#endif
    ts = te
#endif

!$omp parallel default(none) &
!$omp& private(icc,icx,icy,icz,ncc,m) &
!$omp& shared(mylm,wl) &
!$omp& shared(scl_zpr,rccbufp,scl_zmr,rccbufm)
!$omp do
    do icc = 1, scl_zpr(1,0)
      icx = scl_zpr(1,icc)
      icy = scl_zpr(2,icc)
      icz = scl_zpr(3,icc)
      ncc = 2*(icc - 1)*mylm
      do m = 1, mylm
        wl(m,1,icz,icy,icx) = wl(m,1,icz,icy,icx) + rccbufp(ncc + m       )
        wl(m,2,icz,icy,icx) = wl(m,2,icz,icy,icx) + rccbufp(ncc + mylm + m)
      end do
    end do
!$omp end do
!$omp do
    do icc = 1, scl_zmr(1,0)
      icx = scl_zmr(1,icc)
      icy = scl_zmr(2,icc)
      icz = scl_zmr(3,icc)
      ncc = 2*(icc - 1)*mylm
      do m = 1, mylm
        wl(m,1,icz,icy,icx) = wl(m,1,icz,icy,icx) + rccbufm(ncc + m       )
        wl(m,2,icz,icy,icx) = wl(m,2,icz,icy,icx) + rccbufm(ncc + mylm + m)
      end do
    end do
!$omp end do
!$omp end parallel

#ifdef COMLS_TIMER
    te = MPI_WTIME()
    tcomls_lcl = tcomls_lcl + te - ts
    ts = te
#endif

#ifdef DEBUG_MTDFMM
    do icxyz = 1, lcbx*lcby*lcbz
        icz0 = mod(icxyz - 1, lcbz)    + 1 !- mbd_z
        icy0 = mod(icxyz - 1, lcbz*lcby)
        icy0 = icy0/lcbz               + 1 !- mbd_y
        icx0 = (icxyz - 1)/(lcbz*lcby) + 1 !- mbd_x
      do m=1,mylm
        write(myrank+nl*1000+50000,'(2es15.7)') wl(m,1:2,icz0,icy0,icx0)
      enddo ! m
    ENDDO
    call flush(myrank+nl*1000+50000)
#endif

#ifdef FAPP
              if(mdstep.ge.10) CALL fapp_stop("commfmm_low_back",1,1)
#endif
#ifdef FJ_RDMA
              call rdma_deallocate(4)
#endif
  end subroutine comm_fmm_local_super
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to perform MPI communication to gather 
!!         multipole expansion coefficient before the M2M. \n
!!         (grobal communications in the upper levels)
!! \author Yoshimichi Andoh, Shin-ichi Ichikawa
!<
  subroutine comm_fmm_local_top(mylm0, &
       &           wm,nsczdiv,nscydiv,nscxdiv, &
       &           nscellz,nscelly,nscellx, &
       &           mcell_size_z,mcell_size_y,mcell_size_x)
!---------------------------------------------------------------------
!.... Adjacent communication module for FMM upper level moment.    ....
!.... Parallel bucket relay method between processes, which the    ....
!.... relative position within super-cell is same, is employed.    ....
!.... Case of sigle cell for each process in each direction is one ....
!.... major case of parallel bucket relay. At the same time, case  ....
!.... of multiple cell for each process in each direction is also  ....
!.... assumed.                                                     ....
!.... The super-cell tends to be shared by multiple processes.     ....
!.... The moment calculation is assumed to be carried out          ....
!.... redundantly among the processes within the same super-cell.  ....
!.... The number of super-cell per process is one for shared case. ....
!.... For these cases, adjacent communication is established with  ....
!.... far process. This means communication for large super-cell   ....
!.... is inefficient for directly linked mesh network.             ....
!.... Moment calculation (M2L transformation) requires adjacent    ....
!.... 10x10x10 cube of cells. This means even more than whole cells....
!.... with periodic boundary condition are required for large      ....
!.... super-cell. In these cases, same copy of the cell is treated ....
!.... as different cell. Therefore, in contrast with lower level   ....
!.... moment, a strategy to gather no more than whole cells to     ....
!.... processes is superior in network traffic for upper level     ....
!.... moment. By this reason, the cell array, which has local index....
!.... in lower level moment, has global index for upperlevel moment....
!.... When number of cells includes factor 3 three cells are merged....
!.... to construct upper level cell. The level to apply this merge ....
!.... number three is assumed to be upper level, and it is assumed ....
!.... that no merge number two is applied to the level upper than  ....
!.... the the level to which merge number three is applied.        ....
!.... When the cell size is 2^n*3^m the number of processes assumed....
!.... are 1, 3, , 3^m, 2*3^m, 2^2*3^m, 2^3*3^m, , 2^n*3^m in each  ....
!.... dimension.                                                   ....
!
    use domain, only : ncellx, ncelly, ncellz
    use mpi_3d_grid
    implicit none
    integer(kind=wip) mylm,mylm0
    integer(kind=wip) nsczdiv, nscydiv, nscxdiv
    integer(kind=wip) nscellz, nscelly, nscellx
    real(kind=wvp) wm(mylm0,2,nscellz,nscelly,nscellx)
    real(kind=wvp) ccbuf(2*mylm0*nsczdiv*nscydiv*nscxdiv)
    real(kind=wvp) rccbufp(2*mylm0*nsczdiv*nscydiv*nscxdiv,2)
    real(kind=wvp) rccbufm(2*mylm0*nsczdiv*nscydiv*nscxdiv,2)
    integer mcell_size_z, mcell_size_y, mcell_size_x
    integer(kind=wip) np_supercell
    integer(kind=wip) m
    integer(kind=wip) ipz_pdest, ipy_pdest, ipx_pdest
    integer(kind=wip) ipz_psrc, ipy_psrc, ipx_psrc
    integer(kind=wip) ipz_mdest, ipy_mdest, ipx_mdest
    integer(kind=wip) ipz_msrc, ipy_msrc, ipx_msrc
    integer(kind=wip) icx, icy, icz
    integer(kind=wip) icz0, icz1
    integer(kind=wip) icy0, icy1
    integer(kind=wip) icx0, icx1
    integer(kind=wip) iczbp0, iczbp1
    integer(kind=wip) icybp0, icybp1
    integer(kind=wip) icxbp0, icxbp1
    integer(kind=wip) iczbm0, iczbm1
    integer(kind=wip) icybm0, icybm1
    integer(kind=wip) icxbm0, icxbm1
    integer(kind=wip) iczb
    integer(kind=wip) nitr, bidirection
    integer(kind=wip) icybp0prior, icxbp0prior
    integer(kind=wip) icybm0prior, icxbm0prior
    integer(kind=wip) ncc, ncc2
    integer(kind=wip) ibs, ibr

    integer(kind=wip) itr
    integer(kind=wip) istatus(mpi_status_size, 4), ierr
#ifndef SYNC_COM
    integer(kind=wip),dimension(4) :: irq
    integer(kind=wip) nrq
#endif

    mylm=mylm0*1

    !---- common parameters for upper level moment communication. ----
    ! ---- 3D rank order rule. ----
    !     ipx = mod(myrank,npx)
    !     ipy = mod((myrank - ipx)/npx,npy)
    !     ipz = mod((myrank - ipx - ipy*npx)/(npx*npy),npz)
    !      mylm = (nmax+1)*(nmax+1)
    !      nsczdiv = (ncellz / mcell_size_z-1) / npz + 1
    !      nscydiv = (ncelly / mcell_size_y-1) / npy + 1
    !      nscxdiv = (ncellx / mcell_size_x-1) / npx + 1
    !      nscellz = (ncellz - 1) / mcell_size_z + 1
    !      nscelly = (ncelly - 1) / mcell_size_y + 1
    !      nscellx = (ncellx - 1) / mcell_size_x + 1
    !...  global address of process cells.
    icz0 = (ncellz / mcell_size_z * ipz) / npz + 1
    icz1 = icz0 + nsczdiv - 1
    icy0 = (ncelly / mcell_size_y * ipy) / npy + 1
    icy1 = icy0 + nscydiv - 1
    icx0 = (ncellx / mcell_size_x * ipx) / npx + 1
    icx1 = icx0 + nscxdiv - 1

    !---- upper level moment communication starts here. ----

    ! ULmoment +Z
    np_supercell = (npz * mcell_size_z-1) / ncellz + 1
    ipz_pdest = mod(ipz + np_supercell - 1/npz + npz, npz)*npy*npx &
         &     + ipy*npx + ipx
    ipz_psrc  = mod(ipz - np_supercell + 1/npz + npz, npz)*npy*npx &
         &     + ipy*npx + ipx
    ! ULmoment -Z
    ipz_mdest = mod(ipz - np_supercell + 1/npz + npz, npz)*npy*npx &
         &     + ipy*npx + ipx
    ipz_msrc  = mod(ipz + np_supercell - 1/npz + npz, npz)*npy*npx &
         &     + ipy*npx + ipx

    nitr = nscellz / 2
    bidirection = mod(nscellz, 2)   ! when bidirection=1 always bidirectional.
    if(nscellz > npz) then
       nitr = npz / 2
       bidirection = mod(npz, 2)
    endif

    DO itr = 1, nitr
       if (itr == 1) then        ! when npz=2 this code sends to same address.
          iczbp0 = mod(icz0 - nsczdiv - 1 + nscellz, nscellz) + 1
          iczbp1 = mod(icz1 - nsczdiv - 1 + nscellz, nscellz) + 1
          iczbm0 = mod(icz0 + nsczdiv - 1, nscellz) + 1
          iczbm1 = mod(icz1 + nsczdiv - 1, nscellz) + 1
          ncc = 0
          DO icx = icx0, icx1
             DO icy = icy0, icy1
                DO icz = icz0, icz1
                   DO m = 1, mylm
                      ncc = ncc + 1
                      ccbuf(ncc       ) = wm(m, 1, icz, icy, icx )
                      ccbuf(ncc + mylm) = wm(m, 2, icz, icy, icx )
                   END DO
                   ncc = ncc + mylm
                END DO
             END DO
          END DO

#ifdef PROFILE_COMM
    call mpi_barrier(mpi_comm_world, ierr)
    call start_profile(comm_fmm_local_top_id)
#endif
#ifdef MCOM_TIMER
          te = MPI_WTIME()
          tcomul_lcl = tcomul_lcl + te - ts
          call mpi_barrier(mpi_comm_world,ierr)
          ts = MPI_WTIME()
          tcomul_pbar = tcomul_pbar + ts - te
#endif
#ifdef SYNC_COM
          call mpi_sendrecv(ccbuf, ncc, MPI_PREC_MOMVEC_COM, &
               &              ipz_pdest, myrank, &
               &              rccbufp(1,1), ncc, MPI_PREC_MOMVEC_COM, &
               &              ipz_psrc, ipz_psrc, mpi_comm_world, istatus, ierr)
          call mpi_sendrecv(ccbuf, ncc, MPI_PREC_MOMVEC_COM, &
               &              ipz_mdest, myrank, &
               &              rccbufm(1,1), ncc, MPI_PREC_MOMVEC_COM, &
               &              ipz_msrc, ipz_msrc, mpi_comm_world, istatus, ierr)
#else
          call mpi_irecv(rccbufp(1,1), ncc, &
               &              MPI_PREC_MOMVEC_COM, ipz_psrc,   1, &
               &              mpi_comm_world, irq(1), ierr)
          call mpi_isend(ccbuf, ncc, &
               &              MPI_PREC_MOMVEC_COM, ipz_pdest,  1, &
               &              mpi_comm_world, irq(2), ierr)
          call mpi_irecv(rccbufm(1,1), ncc, &
               &              MPI_PREC_MOMVEC_COM, ipz_msrc,   2, &
               &              mpi_comm_world, irq(3), ierr)
          call mpi_isend(ccbuf, ncc, &
               &              MPI_PREC_MOMVEC_COM, ipz_mdest,  2, &
               &              mpi_comm_world, irq(4), ierr)
          nrq = 4
          call mpi_waitall(nrq, irq, istatus, ierr)
#endif

#ifdef PROFILE_COMM
          call end_profile(comm_fmm_local_top_id)
#endif
#ifdef MCOM_TIMER
          te = MPI_WTIME()
          tcomul_com = tcomul_com + te - ts
          nccomul_com = nccomul_com + 1
#ifdef MCOM_TIMER_DETAIL
          write(myrank+NBFU,*) "Rank: ", myrank,
          &                    " TcomUL-Z-1st: ", te - ts
#endif
         ts = te
#endif
#ifdef DBG_COMUL
          write(myrank+NBFU,*) "+Z ipz_pdest,ipz_psrc: ",ipz_pdest,ipz_psrc
          write(myrank+NBFU,*) "+Z itr,iczbp0,iczbp1 : ",itr,iczbp0,iczbp1, &
               &                    " ncc: ",ncc
          write(myrank+NBFU,*) "-Z ipz_mdest,ipz_msrc: ",ipz_mdest,ipz_msrc
          write(myrank+NBFU,*) "-Z itr,iczbm0,iczbm1 : ",itr,iczbm0,iczbm1, &
               &                    " ncc: ",ncc
#endif

          ncc2 = 0
          DO icx = icx0, icx1
             DO icy = icy0, icy1
                DO iczb = iczbp0, iczbp1
                   DO m = 1, mylm
                      ncc2 = ncc2 + 1
                      wm(m,1,iczb,icy,icx) = rccbufp(ncc2,       1)
                      wm(m,2,iczb,icy,icx) = rccbufp(ncc2 + mylm,1)
                   END DO
                   ncc2 = ncc2 + mylm
                END DO
             END DO
          END DO

          ncc2 = 0
          DO icx = icx0, icx1
             DO icy = icy0, icy1
                DO iczb = iczbm0, iczbm1
                   DO m = 1, mylm
                      ncc2 = ncc2 + 1
                      wm(m,1,iczb,icy,icx) = rccbufm(ncc2,       1)
                      wm(m,2,iczb,icy,icx) = rccbufm(ncc2 + mylm,1)
                   END DO
                   ncc2 = ncc2 + mylm
                END DO
             END DO
          END DO

       elseif(itr < nitr .or. bidirection == 1) then

          ibs = mod(itr, 2) + 1
          ibr = mod(itr+1, 2) + 1
          iczbp0 = mod(icz0 - nsczdiv*itr - 1+nscellz, nscellz) + 1
          iczbp1 = mod(iczbp0 + nsczdiv - 1-1+nscellz,  nscellz) +1
          iczbm0 = mod(icz0 + nsczdiv*itr - 1, nscellz) + 1
          iczbm1 = mod(iczbm0 + nsczdiv - 1 - 1, nscellz) + 1

#ifdef PROFILE_COMM
    call mpi_barrier(mpi_comm_world, ierr)
    call start_profile(comm_fmm_local_top_id)
#endif
#ifdef MCOM_TIMER
          te = MPI_WTIME()
          tcomul_lcl = tcomul_lcl + te - ts
          call mpi_barrier(mpi_comm_world,ierr)
          ts = MPI_WTIME()
          tcomul_pbar = tcomul_pbar + ts - te
#endif
#ifdef SYNC_COM
          call mpi_sendrecv(rccbufp(1,ibs), ncc, MPI_PREC_MOMVEC_COM, &
               &              ipz_pdest, myrank, &
               &              rccbufp(1,ibr), ncc, MPI_PREC_MOMVEC_COM, &
               &              ipz_psrc, ipz_psrc, mpi_comm_world, istatus, ierr)
          call mpi_sendrecv(rccbufm(1,ibs), ncc, MPI_PREC_MOMVEC_COM, &
               &              ipz_mdest, myrank, &
               &              rccbufm(1,ibr), ncc, MPI_PREC_MOMVEC_COM, &
               &              ipz_msrc, ipz_msrc, mpi_comm_world, istatus, ierr)
#else
          call mpi_irecv(rccbufp(1,ibr), ncc, &
               &              MPI_PREC_MOMVEC_COM, ipz_psrc,   1, &
               &              mpi_comm_world, irq(1), ierr)
          call mpi_isend(rccbufp(1,ibs), ncc, &
               &              MPI_PREC_MOMVEC_COM, ipz_pdest,  1, &
               &              mpi_comm_world, irq(2), ierr)
          call mpi_irecv(rccbufm(1,ibr), ncc, &
               &              MPI_PREC_MOMVEC_COM, ipz_msrc,   2, &
               &              mpi_comm_world, irq(3), ierr)
          call mpi_isend(rccbufm(1,ibs), ncc, &
               &              MPI_PREC_MOMVEC_COM, ipz_mdest,  2, &
               &              mpi_comm_world, irq(4), ierr)
          nrq = 4
          call mpi_waitall(nrq, irq, istatus, ierr)
#endif

#ifdef PROFILE_COMM
          call end_profile(comm_fmm_local_top_id)
#endif
#ifdef MCOM_TIMER
          te = MPI_WTIME()
          tcomul_com = tcomul_com + te - ts
          nccomul_com = nccomul_com + 1
#ifdef MCOM_TIMER_DETAIL

          write(myrank+NBFU,*) "Rank: ", myrank, &
               &                    " TcomUL-Z-2nd&later: ", te - ts
#endif
          ts = te
#endif
#ifdef DBG_COMUL
          write(myrank+NBFU,*) "+Z itr,iczbp0,iczbp1 : ",itr,iczbp0,iczbp1, &
               &                    " ncc: ",ncc
          write(myrank+NBFU,*) "-Z itr,iczbm0,iczbm1 : ",itr,iczbm0,iczbm1, &
               &                    " ncc: ",ncc
#endif

          ncc2 = 0
          DO icx = icx0, icx1
             DO icy = icy0, icy1
                DO iczb = iczbp0, iczbp1
                   DO m = 1, mylm
                      ncc2 = ncc2 + 1
                      wm(m,1,iczb,icy,icx) = rccbufp(ncc2,       ibr)
                      wm(m,2,iczb,icy,icx) = rccbufp(ncc2 + mylm,ibr)
                   END DO
                   ncc2 = ncc2 + mylm
                END DO
             END DO
          END DO

          ncc2 = 0
          DO icx = icx0, icx1
             DO icy = icy0, icy1
                DO iczb = iczbm0, iczbm1
                   DO m = 1, mylm
                      ncc2 = ncc2 + 1
                      wm(m,1,iczb,icy,icx) = rccbufm(ncc2,       ibr)
                      wm(m,2,iczb,icy,icx) = rccbufm(ncc2 + mylm,ibr)
                   END DO
                   ncc2 = ncc2 + mylm
                END DO
             END DO
          END DO

       else   ! final iteration for (+Z).

          ibs = mod(itr, 2) + 1
          ibr = mod(itr+1, 2) + 1
          iczbp0 = mod(icz0 - nsczdiv*itr - 1+ nscellz, nscellz) + 1
          iczbp1 = mod(iczbp0 + nsczdiv -1-1+ nscellz,  nscellz) +1

#ifdef PROFILE_COMM
    call mpi_barrier(mpi_comm_world, ierr)
    call start_profile(comm_fmm_local_top_id)
#endif
#ifdef MCOM_TIMER
          te = MPI_WTIME()
          tcomul_lcl = tcomul_lcl + te - ts
          call mpi_barrier(mpi_comm_world,ierr)
          ts = MPI_WTIME()
          tcomul_pbar = tcomul_pbar + ts - te
#endif

          call mpi_sendrecv(rccbufp(1,ibs), ncc, MPI_PREC_MOMVEC_COM, &
               &             ipz_pdest, &
               &             myrank, rccbufp(1,ibr), ncc, MPI_PREC_MOMVEC_COM, &
               &             ipz_psrc, ipz_psrc, mpi_comm_world,istatus,ierr)

#ifdef PROFILE_COMM
          call end_profile(comm_fmm_local_top_id)
#endif
#ifdef MCOM_TIMER
          te = MPI_WTIME()
          tcomul_com = tcomul_com + te - ts
          nccomul_com = nccomul_com + 1
#ifdef MCOM_TIMER_DETAIL
          write(myrank+NBFU,*) "Rank: ", myrank, &
               &                    " TcomUL-Z-final-for(+Z): ", te - ts
#endif
          ts = te
#endif
#ifdef DBG_COMUL
          write(myrank+NBFU,*) "+Z itr,iczbp0,iczbp1 : ",itr,iczbp0,iczbp1, &
               &                    " ncc: ",ncc
#endif

          ncc2 = 0
          DO icx = icx0, icx1
             DO icy = icy0, icy1
                DO iczb = iczbp0, iczbp1
                   DO m = 1, mylm
                      ncc2 = ncc2 + 1
                      wm(m,1,iczb,icy,icx) = rccbufp(ncc2,       ibr)
                      wm(m,2,iczb,icy,icx) = rccbufp(ncc2 + mylm,ibr)
                   END DO
                   ncc2 = ncc2 + mylm
                END DO
             END DO
          END DO

       end if
    END DO

    ! ULmoment +Y
    np_supercell = (npy * mcell_size_y-1) / ncelly + 1
    ipy_pdest = ipz*npy*npx &
         &     + mod(ipy + np_supercell - 1/npy + npy, npy)*npx + ipx
    ipy_psrc = ipz*npy*npx &
         &     + mod(ipy - np_supercell + 1/npy + npy, npy)*npx + ipx
    ! ULmoment -Y
    ipy_mdest = ipz*npy*npx &
         &     + mod(ipy - np_supercell + 1/npy + npy, npy)*npx + ipx
    ipy_msrc = ipz*npy*npx &
         &     + mod(ipy + np_supercell - 1/npy + npy, npy)*npx + ipx

    nitr = nscelly / 2
    bidirection = mod(nscelly, 2)   ! when bidirection=1 always bidirectional.
    if(nscelly > npy) then
       nitr = npy / 2
       bidirection = mod(npy, 2)
    endif

    DO icx = icx0, icx1
       DO itr = 1, nitr
          if (itr == 1) then
             icybp0 = mod(icy0 - nscydiv - 1+ nscelly, nscelly) + 1
             icybp1 = mod(icy1 - nscydiv - 1+ nscelly, nscelly) + 1
             icybm0 = mod(icy0 + nscydiv - 1, nscelly) + 1
             icybm1 = mod(icy1 + nscydiv - 1, nscelly) + 1
             ncc = 2 * mylm * nscellz * (icy1 - icy0 + 1)

#ifdef PROFILE_COMM
    call mpi_barrier(mpi_comm_world, ierr)
    call start_profile(comm_fmm_local_top_id)
#endif
#ifdef MCOM_TIMER
             te = MPI_WTIME()
             tcomul_lcl = tcomul_lcl + te - ts
             call mpi_barrier(mpi_comm_world,ierr)
             ts = MPI_WTIME()
             tcomul_pbar = tcomul_pbar + ts - te
#endif
#ifdef SYNC_COM
             call mpi_sendrecv(wm(1,1,1,icy0,icx), ncc, &
                  &                MPI_PREC_MOMVEC_COM, ipy_pdest, myrank, &
                  &                wm(1,1,1,icybp0,icx), ncc, MPI_PREC_MOMVEC_COM, &
                  &                ipy_psrc, ipy_psrc, mpi_comm_world, &
                  &                istatus, ierr)
             call mpi_sendrecv(wm(1,1,1,icy0,icx), ncc, &
                  &                MPI_PREC_MOMVEC_COM, ipy_mdest, myrank, &
                  &                wm(1,1,1,icybm0,icx), ncc, MPI_PREC_MOMVEC_COM, &
                  &                ipy_msrc, ipy_msrc, mpi_comm_world, &
                  &                istatus, ierr )
#else
             call mpi_irecv(wm(1,1,1,icybp0,icx), ncc, &
                  &                MPI_PREC_MOMVEC_COM, ipy_psrc,     1, &
                  &                mpi_comm_world, irq(1), ierr)
             call mpi_isend(wm(1,1,1,icy0,icx), ncc, &
                  &                MPI_PREC_MOMVEC_COM, ipy_pdest,    1, &
                  &                mpi_comm_world, irq(2), ierr)
             call mpi_irecv(wm(1,1,1,icybm0,icx), ncc, &
                  &                MPI_PREC_MOMVEC_COM, ipy_msrc,     2, &
                  &                mpi_comm_world, irq(3), ierr)
             call mpi_isend(wm(1,1,1,icy0,icx), ncc, &
                  &                MPI_PREC_MOMVEC_COM, ipy_mdest,    2, &
                  &                mpi_comm_world, irq(4), ierr)
             nrq = 4
             call mpi_waitall(nrq, irq, istatus, ierr)
#endif

#ifdef PROFILE_COMM
          call end_profile(comm_fmm_local_top_id)
#endif
#ifdef MCOM_TIMER
             te = MPI_WTIME()
             tcomul_com = tcomul_com + te - ts
             nccomul_com = nccomul_com + 1
#ifdef MCOM_TIMER_DETAIL
             write(myrank+NBFU,*) "Rank: ", myrank, &
                  &                    " TcomUL-Y-1st: ", te - ts
#endif
             ts = te
#endif
#ifdef DBG_COMUL
             write(myrank+NBFU,*) "+Y ipy_pdest,ipy_psrc: ",ipy_pdest,ipy_psrc
             write(myrank+NBFU,*) "+Y itr,icybp0,icybp1 : ",itr,icybp0,icybp1, &
                  &                    " ncc: ",ncc," icx: ",icx
             write(myrank+NBFU,*) "-Y ipy_mdest,ipy_msrc: ",ipy_mdest,ipy_msrc
             write(myrank+NBFU,*) "-Y itr,icybm0,icybm1 : ",itr,icybm0,icybm1, &
                  &                    " ncc: ",ncc," icx: ",icx
#endif

          elseif(itr < nitr .or. bidirection == 1) then

             icybp0prior = icybp0
             icybp0 = mod(icy0 - nscydiv * itr - 1 + nscelly, &
                  &                     nscelly) + 1
             icybm0prior = icybm0
             icybm0 = mod(icy0 + nscydiv * itr - 1, nscelly) + 1

#ifdef PROFILE_COMM
    call mpi_barrier(mpi_comm_world, ierr)
    call start_profile(comm_fmm_local_top_id)
#endif
#ifdef MCOM_TIMER
             te = MPI_WTIME()
             tcomul_lcl = tcomul_lcl + te - ts
             call mpi_barrier(mpi_comm_world,ierr)
             ts = MPI_WTIME()
             tcomul_pbar = tcomul_pbar + ts - te
#endif
#ifdef SYNC_COM
             call mpi_sendrecv(wm(1,1,1,icybp0prior,icx), ncc, &
                  &                   MPI_PREC_MOMVEC_COM, ipy_pdest, myrank, &
                  &                   wm(1,1,1,icybp0,icx), ncc, MPI_PREC_MOMVEC_COM, &
                  &                   ipy_psrc, ipy_psrc, mpi_comm_world, &
                  &                   istatus, ierr )
             call mpi_sendrecv(wm(1,1,1,icybm0prior,icx), ncc, &
                  &                   MPI_PREC_MOMVEC_COM, ipy_mdest, myrank, &
                  &                   wm(1,1,1,icybm0,icx), ncc, MPI_PREC_MOMVEC_COM, &
                  &                   ipy_msrc, ipy_msrc, mpi_comm_world, &
                  &                   istatus, ierr )
#else
             call mpi_irecv(wm(1,1,1,icybp0,icx), ncc, &
                  &                   MPI_PREC_MOMVEC_COM, ipy_psrc,  1, &
                  &                   mpi_comm_world, irq(1), ierr)
             call mpi_isend(wm(1,1,1,icybp0prior,icx), ncc, &
                  &                   MPI_PREC_MOMVEC_COM, ipy_pdest, 1, &
                  &                   mpi_comm_world, irq(2), ierr)
             call mpi_irecv(wm(1,1,1,icybm0,icx), ncc, &
                  &                   MPI_PREC_MOMVEC_COM, ipy_msrc,  2, &
                  &                   mpi_comm_world, irq(3), ierr)
             call mpi_isend(wm(1,1,1,icybm0prior,icx), ncc, &
                  &                   MPI_PREC_MOMVEC_COM, ipy_mdest, 2, &
                  &                   mpi_comm_world, irq(4), ierr)
             nrq = 4
             call mpi_waitall(nrq, irq, istatus, ierr)
#endif

#ifdef PROFILE_COMM
          call end_profile(comm_fmm_local_top_id)
#endif
#ifdef MCOM_TIMER
             te = MPI_WTIME()
             tcomul_com = tcomul_com + te - ts
             nccomul_com = nccomul_com + 1
#ifdef MCOM_TIMER_DETAIL
             write(myrank+NBFU,*) "Rank: ", myrank, &
                  &                    " TcomUL-Y-2nd&later: ", te - ts
#endif
             ts = te
#endif
#ifdef DBG_COMUL
             write(myrank+NBFU,*) "+Y itr,icybp0prior,icybp0,ncc : ", &
                  &                    itr,icybp0prior,icybp0,ncc," icx: ",icx
             write(myrank+NBFU,*) "-Y itr,icybm0prior,icybm0,ncc : ", &
                  &                    itr,icybm0prior,icybm0,ncc," icx: ",icx
#endif

          else     ! final iteration for (+Y).

             icybp0prior = icybp0
             icybp0 = mod(icy0 - nscydiv * itr - 1 + nscelly, &
                  &                     nscelly) + 1

#ifdef PROFILE_COMM
    call mpi_barrier(mpi_comm_world, ierr)
    call start_profile(comm_fmm_local_top_id)
#endif
#ifdef MCOM_TIMER
             te = MPI_WTIME()
             tcomul_lcl = tcomul_lcl + te - ts
             call mpi_barrier(mpi_comm_world,ierr)
             ts = MPI_WTIME()
             tcomul_pbar = tcomul_pbar + ts - te
#endif

             call mpi_sendrecv(wm(1,1,1,icybp0prior,icx), ncc, &
                  &                MPI_PREC_MOMVEC_COM, &
                  &                ipy_pdest, myrank, wm(1,1,1,icybp0,icx), ncc, &
                  &                MPI_PREC_MOMVEC_COM, &
                  &                ipy_psrc, ipy_psrc, mpi_comm_world,istatus,ierr)

#ifdef PROFILE_COMM
          call end_profile(comm_fmm_local_top_id)
#endif
#ifdef MCOM_TIMER
             te = MPI_WTIME()
             tcomul_com = tcomul_com + te - ts
             nccomul_com = nccomul_com + 1
#ifdef MCOM_TIMER_DETAIL
             write(myrank+NBFU,*) "Rank: ", myrank, &
                  &                    " TcomUL-Y-final-for(+Y): ", te - ts
#endif
             ts = te
#endif
#ifdef DBG_COMUL
             write(myrank+NBFU,*) "+Y itr,icybp0prior,icybp0,ncc : ", &
                  &                    itr,icybp0prior,icybp0,ncc," icx: ",icx
#endif

          endif
       END DO
    END DO

    ! ULmoment +X
    np_supercell = (npx * mcell_size_x-1) / ncellx + 1
    ipx_pdest = ipz*npy*npx + ipy*npx &
         &     + mod(ipx + np_supercell - 1/npx + npx, npx)
    ipx_psrc = ipz*npy*npx + ipy*npx &
         &     + mod(ipx - np_supercell + 1/npx + npx, npx)
    ! ULmoment -X
    ipx_mdest = ipz*npy*npx + ipy*npx &
         &     + mod(ipx - np_supercell + 1/npx + npx, npx)
    ipx_msrc = ipz*npy*npx + ipy*npx &
         &     + mod(ipx + np_supercell - 1/npx + npx, npx)

    nitr = nscellx / 2
    bidirection = mod(nscellx, 2)   ! when bidirection=1 always bidirectional.
    if(nscellx > npx) then
       nitr = npx / 2
       bidirection = mod(npx, 2)
    endif

    DO itr = 1, nitr
       if (itr ==1) then
          icxbp0 = mod(icx0 - nscxdiv - 1 + nscellx, nscellx) + 1
          icxbp1 = mod(icx1 - nscxdiv - 1 + nscellx, nscellx) + 1
          icxbm0 = mod(icx0 + nscxdiv - 1, nscellx) + 1
          icxbm1 = mod(icx1 + nscxdiv - 1, nscellx) + 1
          ncc = 2 * mylm * nscellz * nscelly * (icx1 - icx0 + 1)

#ifdef PROFILE_COMM
    call mpi_barrier(mpi_comm_world, ierr)
    call start_profile(comm_fmm_local_top_id)
#endif
#ifdef MCOM_TIMER
          te = MPI_WTIME()
          tcomul_lcl = tcomul_lcl + te - ts
          call mpi_barrier(mpi_comm_world,ierr)
          ts = MPI_WTIME()
          tcomul_pbar = tcomul_pbar + ts - te
#endif
#ifdef SYNC_COM
          call mpi_sendrecv(wm(1,1,1,1,icx0), ncc, &
               &             MPI_PREC_MOMVEC_COM, ipx_pdest, myrank, &
               &             wm(1,1,1,1,icxbp0), ncc, MPI_PREC_MOMVEC_COM, &
               &             ipx_psrc, ipx_psrc, mpi_comm_world, istatus, ierr)
          call mpi_sendrecv(wm(1,1,1,1,icx0), ncc, &
               &             MPI_PREC_MOMVEC_COM, ipx_mdest, myrank, &
               &             wm(1,1,1,1,icxbm0), ncc, MPI_PREC_MOMVEC_COM, &
               &             ipx_msrc, ipx_msrc, mpi_comm_world, istatus, ierr)
#else
          call mpi_irecv(wm(1,1,1,1,icxbp0), ncc, &
               &             MPI_PREC_MOMVEC_COM, ipx_psrc,    1, &
               &             mpi_comm_world, irq(1), ierr)
          call mpi_isend(wm(1,1,1,1,icx0), ncc, &
               &             MPI_PREC_MOMVEC_COM, ipx_pdest,   1, &
               &             mpi_comm_world, irq(2), ierr)
          call mpi_irecv(wm(1,1,1,1,icxbm0), ncc, &
               &             MPI_PREC_MOMVEC_COM, ipx_msrc,    2, &
               &             mpi_comm_world, irq(3), ierr)
          call mpi_isend(wm(1,1,1,1,icx0), ncc, &
               &             MPI_PREC_MOMVEC_COM, ipx_mdest,   2, &
               &             mpi_comm_world, irq(4), ierr)
          nrq = 4
          call mpi_waitall(nrq, irq, istatus, ierr)
#endif

#ifdef PROFILE_COMM
          call end_profile(comm_fmm_local_top_id)
#endif
#ifdef MCOM_TIMER
          te = MPI_WTIME()
          tcomul_com = tcomul_com + te - ts
          nccomul_com = nccomul_com + 1
#ifdef MCOM_TIMER_DETAIL
          write(myrank+NBFU,*) "Rank: ", myrank, &
               &                    " TcomUL-X-1st: ", te - ts
#endif
          ts = te
#endif
#ifdef DBG_COMUL
          write(myrank+NBFU,*) "+X ipx_pdest,ipx_psrc: ",ipx_pdest,ipx_psrc
          write(myrank+NBFU,*) "+X itr,icxbp0,icxbp1 : ",itr,icxbp0,icxbp1, &
               &                    " ncc: ",ncc
          write(myrank+NBFU,*) "-X ipx_mdest,ipx_msrc: ",ipx_mdest,ipx_msrc
          write(myrank+NBFU,*) "-X itr,icxbm0,icxbm1 : ",itr,icxbm0,icxbm1, &
               &                    " ncc: ",ncc
#endif

       elseif(itr < nitr .or. bidirection == 1) then

          icxbp0prior = icxbp0
          icxbp0 = mod(icx0 - nscxdiv*itr -1+ nscellx, nscellx) + 1
          icxbm0prior = icxbm0
          icxbm0 = mod(icx0 + nscxdiv * itr -1, nscellx) + 1

#ifdef PROFILE_COMM
    call mpi_barrier(mpi_comm_world, ierr)
    call start_profile(comm_fmm_local_top_id)
#endif
#ifdef MCOM_TIMER
          te = MPI_WTIME()
          tcomul_lcl = tcomul_lcl + te - ts
          call mpi_barrier(mpi_comm_world,ierr)
          ts = MPI_WTIME()
          tcomul_pbar = tcomul_pbar + ts - te
#endif
#ifdef SYNC_COM
          call mpi_sendrecv(wm(1,1,1,1,icxbp0prior), ncc, &
               &                MPI_PREC_MOMVEC_COM, ipx_pdest, myrank,  &
               &                wm(1,1,1,1,icxbp0), ncc, MPI_PREC_MOMVEC_COM, &
               &                ipx_psrc, ipx_psrc, mpi_comm_world, &
               &                istatus, ierr )
          call mpi_sendrecv(wm(1,1,1,1,icxbm0prior), ncc, &
               &                MPI_PREC_MOMVEC_COM, ipx_mdest, myrank, &
               &                wm(1,1,1,1,icxbm0), ncc, MPI_PREC_MOMVEC_COM, &
               &                ipx_msrc, ipx_msrc, mpi_comm_world, &
               &                istatus, ierr )
#else
          call mpi_irecv(wm(1,1,1,1,icxbp0), ncc, &
               &                MPI_PREC_MOMVEC_COM, ipx_psrc,   1, &
               &                mpi_comm_world, irq(1), ierr)
          call mpi_isend(wm(1,1,1,1,icxbp0prior), ncc, &
               &                MPI_PREC_MOMVEC_COM, ipx_pdest,  1, &
               &                mpi_comm_world, irq(2), ierr)
          call mpi_irecv(wm(1,1,1,1,icxbm0), ncc, &
               &                MPI_PREC_MOMVEC_COM, ipx_msrc,   2, &
               &                mpi_comm_world, irq(3), ierr)
          call mpi_isend(wm(1,1,1,1,icxbm0prior), ncc, &
               &                MPI_PREC_MOMVEC_COM, ipx_mdest,  2, &
               &                mpi_comm_world, irq(4), ierr)
          nrq = 4
          call mpi_waitall(nrq, irq, istatus, ierr)
#endif

#ifdef PROFILE_COMM
          call end_profile(comm_fmm_local_top_id)
#endif
#ifdef MCOM_TIMER
          te = MPI_WTIME()
          tcomul_com = tcomul_com + te - ts
          nccomul_com = nccomul_com + 1
#ifdef MCOM_TIMER_DETAIL
          write(myrank+NBFU,*) "Rank: ", myrank, &
               &                    " TcomUL-X-2nd&later: ", te - ts
#endif
          ts = te
#endif
#ifdef DBG_COMUL
          write(myrank+NBFU,*) "+X itr,icxbp0prior,icxbp0,ncc : ", &
               &                    itr,icxbp0prior,icxbp0,ncc
          write(myrank+NBFU,*) "-X itr,icxbm0prior,icxbm0,ncc : ", &
               &                    itr,icxbm0prior,icxbm0,ncc
#endif

       else    ! final iteration for (+X).

          icxbp0prior = icxbp0
          icxbp0 = mod(icx0 - nscxdiv*itr -1+ nscellx, nscellx) + 1

#ifdef PROFILE_COMM
    call mpi_barrier(mpi_comm_world, ierr)
    call start_profile(comm_fmm_local_top_id)
#endif
#ifdef MCOM_TIMER
          te = MPI_WTIME()
          tcomul_lcl = tcomul_lcl + te - ts
          call mpi_barrier(mpi_comm_world,ierr)
          ts = MPI_WTIME()
          tcomul_pbar = tcomul_pbar + ts - te
#endif

          call mpi_sendrecv(wm(1,1,1,1,icxbp0prior), ncc, &
               &             MPI_PREC_MOMVEC_COM, &
               &             ipx_pdest, myrank, wm(1,1,1,1,icxbp0), ncc, &
               &             MPI_PREC_MOMVEC_COM, &
               &             ipx_psrc, ipx_psrc, mpi_comm_world,istatus,ierr)

#ifdef PROFILE_COMM
          call end_profile(comm_fmm_local_top_id)
#endif
#ifdef MCOM_TIMER
          te = MPI_WTIME()
          tcomul_com = tcomul_com + te - ts
          nccomul_com = nccomul_com + 1
#ifdef MCOM_TIMER_DETAIL
          write(myrank+NBFU,*) "Rank: ", myrank, &
               &                    " TcomUL-X-final-for(+X): ", te - ts
#endif
          ts = te
#endif
#ifdef DBG_COMUL
          write(myrank+NBFU,*) "+X itr,icxbp0prior,icxbp0,ncc : ", &
               &                    itr,icxbp0prior,icxbp0,ncc
#endif

       end if
    END DO

  end subroutine comm_fmm_local_top

!---------------------------------------------------------------------
!>
!! \brief  Subroutine to perform MPI communication to gather 
!!         multipole expansion coefficient before the M2M. \n
!!         (local communications in the lower levels)
!! \author Yoshimichi Andoh, Shin-ichi Ichikawa
!<
  subroutine comm_fmm_local_multi_dr(mylm0, nlz, nly, nlx, &
       &           wm,nbsize,nsczdiv,nscydiv,nscxdiv, &
       &           nscellz, nscelly,nscellx, &
       &           nbound_zm, nbound_ym, nbound_xm, &
       &           nbound_zp, nbound_yp, nbound_xp, il)
!---------------------------------------------------------------------
!.... Adjacent communication module for FMM lower level moment.    ....
!.... Parallel bucket relay method between processes, whose        ....
!.... relative position within super-cell are same, is employed.   ....
!.... Case of multiple cell for each process in each direction is  ....
!.... one major case of parallel bucket relay. At the same time,   ....
!.... case of sigle cell for each process in each direction is     ....
!.... also assumed.                                                ....
!.... The super-cell tends to be shared by multiple processes.     ....
!.... The moment calculation is assumed to be carried out          ....
!.... redundantly among the processes within the same super-cell.  ....
!.... The number of super-cell per process is one for shared case. ....
!.... For these cases, adjacent communication is established with  ....
!.... far process. This means communication for large super-cell   ....
!.... is inefficient for directly linked mesh network.             ....
!.... Since single cell case tends to be dominant for upper level  ....
!.... moment, this subroutine LLmoment is inadequate. Alternative  ....
!.... is ULmoment.                                                 ....
!.... Values of the number of extra-process adjacent cells, which  ....
!.... are referred for moment calculation, are expressed by the    ....
!.... variables nbound_zm, nbound_zp, and so on.                   ....
!.... When number of cells includes factor 3 three cells are merged....
!.... to construct upper level cell. The level to apply this merge ....
!.... number three is assumed to be upper level, and it is assumed ....
!.... that the unification subroutine ULmoment is applied to these ....
!.... upper levels. This subroutine LLmoment is assumed to be used ....
!.... for the levels. Even that, localized data structure can be   ....
!.... advatageous for the level where cell merge number for upper  ....
!.... level is three for cell size factor of larger power of 3.    ....
!.... This is reflected to boundary cell number 6, 7, and 8.       ....
!.... In some cases of process number irregular pairing of         ....
!.... bucket-relay communication is required. This is controled by ....
!.... relative position of beginning cell or ending cell within    ....
!.... merging cell group. Also, global address of beginning cell   ....
!.... or ending cell is considered to determine destination rank   ....
!.... boundary cells.                                              ....
!.... When the cell size is 2^n*3^m the number of processes assumed....
!.... are 1, 3, , 3^m, 2*3^m, 2^2*3^m, 2^3*3^m, , 2^n*3^m in each  ....
!.... dimension.                                                   ....
    use dr_cntl, only : mtbd, mtbdp2
    use mpi_3d_grid
    implicit none
    integer(kind=wip) mylm,mylm0
    integer(kind=wip), intent(in) :: nlz, nly, nlx
    integer(kind=wip) nbsize
    integer(kind=wip) nscxdiv, nscydiv, nsczdiv
    integer(kind=wip), intent(in) :: nbound_zm, nbound_ym, nbound_xm ! the number of minus boundary cells.
    integer(kind=wip), intent(in) :: nbound_zp, nbound_yp, nbound_xp ! the number of plus boundary cells.
    integer(kind=wip), intent(in) :: il
    real(kind=wvp) wm(mylm0,2, nbound_zm + nsczdiv + nbound_zp, &
   &                           nbound_ym + nscydiv + nbound_yp, &
   &                           nbound_xm + nscxdiv + nbound_xp)
#ifndef FJ_RDMA
    real(kind=wvp) ccbufp (2*mylm0*nbsize*nscydiv*nscxdiv)
    real(kind=wvp) rccbufp(2*mylm0*nbsize*nscydiv*nscxdiv, 2)
    real(kind=wvp) ccbufm (2*mylm0*nbsize*nscydiv*nscxdiv)
    real(kind=wvp) rccbufm(2*mylm0*nbsize*nscydiv*nscxdiv, 2)
#endif
    integer(kind=wip) nscellx, nscelly, nscellz

    integer(kind=wip) np_supercell
    integer(kind=wip) ipz_pdest, ipy_pdest, ipx_pdest
    integer(kind=wip) ipz_psrc, ipy_psrc, ipx_psrc
    integer(kind=wip) ipz_mdest, ipy_mdest, ipx_mdest
    integer(kind=wip) ipz_msrc, ipy_msrc, ipx_msrc

    integer(kind=wip) nb_destp, nb_destm
    integer(kind=wip) npow
    integer(kind=wip) nitr
    integer(kind=wip) itr

    integer(kind=wip) icz, icy, icx
    integer(kind=wip) iczp0 , iczp1
    integer(kind=wip) iczbp0, iczbp1
    integer(kind=wip) iczm0 , iczm1
    integer(kind=wip) iczbm0, iczbm1
    integer(kind=wip) iczb
    integer(kind=wip) icyp0 , icyp1
    integer(kind=wip) icybp0, icybp1 ! rdma
    integer(kind=wip) icym0 , icym1
    integer(kind=wip) icybm0, icybm1 ! rdma
    integer(kind=wip) icxp0 , icxp1
    integer(kind=wip) icxbp0, icxbp1
    integer(kind=wip) icxm0 , icxm1
    integer(kind=wip) icxbm0, icxbm1

    integer(kind=wip) iczg0, iczg1    ! global cell address of starting and final cell.
    integer(kind=wip) icyg0, icyg1    ! global cell address of starting and final cell.
    integer(kind=wip) icxg0, icxg1    ! global cell address of starting and final cell.

    integer(kind=wip) ncc
    integer(kind=wip) nccps, nccpr
    integer(kind=wip) nccms, nccmr

    integer(kind=wip) m
    integer(kind=wip) ibs, ibr

    integer(kind=wip) istatus(mpi_status_size, 4), ierr
    integer(kind=wip),dimension(4) :: irq
    integer(kind=wip) nrq

#ifdef DEBUG_RDMAYA
    integer hoge
#endif
#ifdef FJ_RDMA
    real(kind=wvp),pointer :: ccbufp (:)
    real(kind=wvp),pointer :: rccbufp(:,:)
    real(kind=wvp),pointer :: ccbufm (:)
    real(kind=wvp),pointer :: rccbufm(:,:)
    integer(8), pointer :: heap_raddr(:)
    integer(8) :: heap_laddr
    integer(8) :: wm_laddr, wm_4th_offset, wm_3th_offset

    integer(8), pointer :: wm_raddr(:)
    type(c_ptr) :: cptr
    integer(8)  :: ccbufp_offset
    integer(8)  :: rccbufp_offset
    integer(8)  :: ccbufm_offset
    integer(8)  :: rccbufm_offset
    integer(4) ids   ! data size in integer
    integer(4) ibds  ! data size for one data-block
#if defined(PRECISION_M2L_SP) || defined(PRECISION_M2L_MIX)
    ibds=8    ! 4*2 SP for wvp
    ids=4     ! byte real(4),   used as nccp*ids, nccm*ids
#else
    ibds=16   ! 8*2 DP for wvp
    ids=8     ! byte real(8)
#endif

    cptr = rdma_allocate_coarray((2*mylm0*nbsize*nscydiv*nscxdiv)*ibds, ccbufp_offset)
    call c_f_pointer(cptr, fptr=ccbufp, shape=[2*mylm0*nbsize*nscydiv*nscxdiv])
    cptr = rdma_allocate_coarray((2*mylm0*nbsize*nscydiv*nscxdiv*2)*ibds, rccbufp_offset)
    call c_f_pointer(cptr, fptr=rccbufp, shape=[2*mylm0*nbsize*nscydiv*nscxdiv, 2])
    cptr = rdma_allocate_coarray((2*mylm0*nbsize*nscydiv*nscxdiv)*ibds, ccbufm_offset)
    call c_f_pointer(cptr, fptr=ccbufm, shape=[2*mylm0*nbsize*nscydiv*nscxdiv])
    cptr = rdma_allocate_coarray((2*mylm0*nbsize*nscydiv*nscxdiv*2)*ibds, rccbufm_offset)
    call c_f_pointer(cptr, fptr=rccbufm, shape=[2*mylm0*nbsize*nscydiv*nscxdiv, 2])
    !
    heap_laddr = rdma_get_heap_laddr()
    cptr       = rdma_get_heap_raddr()
    call c_f_pointer(cptr, fptr=heap_raddr, shape=[nprocs])
    !
    wm_laddr = rdma_get_laddr(wm)
    cptr     = rdma_get_raddr(wm)
    call c_f_pointer(cptr, fptr=wm_raddr, shape=[nprocs])
    wm_4th_offset = mylm0*(nbound_zm+nsczdiv+nbound_zp)*(nbound_ym+nscydiv+nbound_yp)
    wm_3th_offset = mylm0*(nbound_zm+nsczdiv+nbound_zp)
#endif  /** FJ_RDMA **/

    mylm=mylm0*1
    ! ---- 3D rank order rule. ----
    !     ipx=mod(myrank, npx)
    !     ipy=mod((myrank - ipx) / npx, npy)
    !     ipz=mod((myrank - ipx - ipy*npx) / (npx*npy), npz)

#ifdef DBG_COMLL
    write(myrank+NBFU,*) "npz,npy,npx:",npz,npy,npx
    write(myrank+NBFU,*) "myrank,ipz,ipy,ipx:",myrank,ipz,ipy,ipx
#endif

    !----- common parameters for lower level moment communication. -----
    !      nsczdiv = (ncellz / mcell_size_z - 1) / npz + 1  ! super cells/proc.
    !      nscydiv = (ncelly / mcell_size_y - 1) / npy + 1  ! super cells/proc.
    !      nscxdiv = (ncellx / mcell_size_x - 1) / npx + 1  ! super cells/proc.
    ! ---------- calculation of extra process boundary cell size. ----------
    ! ... mtbd = 3, mtbd+1 = 4, and mtbdp2 = 5 for DR method of 3 cells merging.
    !      nscellz = ncellz / mcell_size_z
    ! ... npow=2 means number of merging cells for upper level is 2, whereas
    ! ... npow=3 means number of merging cells for upper level is 3.
    ! ... the number of super-cells are always even number when 2 cells are merged
    ! ... to make upper level cell, while it is power of 3 number when 3 cells are
    ! ... merged.
    !      npow = nlz
    ! ... global cell address of starting cell; iczg0 and ending cell; iczg1.
    !      mbd = 2
    !      mtbd = 3
    !      iczg0 = (nscellz * ipz) / npz + 1
    !      iczg1 = (nscellz * (ipz + 1) - 1) / npz + 1
    !      if(npow == 2) then
    !         nbound_zm = mbd + 1 - mod(iczg0, npow)      ! number of minus boundary cells.
    !         nbound_zp = mbd     + mod(iczg1, npow)      ! number of plus boundary cells.
    !      else
    !         nbound_zm = mtbd     + mod(iczg0 + 2, npow) ! number of minus boundary cells.
    !         nbound_zp = mtbd + 2 - mod(iczg1 + 2, npow) ! number of plus boundary cells.
    !      endif
    !     likewise for Y-direction and X-direction.

    !----- lower level moment communication starts here. -----

    ! LLmoment +Z
    ! the following is the number of procceses per cell constructing process group in a
    ! super-cell. which is equal to 1 when the number of global cell is larger than the
    ! number of processes. 
    !   npz          : the number of process in Z-direction.
    !   nscellz      : the global number of super cells in Z-direction.

      np_supercell = (npz - 1) / nscellz + 1

    ! destination/source 1-dimensional rank for neighbor communications in plus-direction.
    ! neighbor rank means the rank relatively same position in the neighbor process group
    ! in neighbor super-cell.

      ipz_pdest = mod(ipz + np_supercell - 1/npz + npz, npz)*npy*npx + ipy*npx + ipx
      ipz_psrc  = mod(ipz - np_supercell + 1/npz + npz, npz)*npy*npx + ipy*npx + ipx

    ! LLmoment -Z
    ! destination/source 1-dimensional rank for neighbor communications in minus-direction.

      ipz_mdest = mod(ipz - np_supercell + 1/npz  + npz, npz)*npy*npx + ipy*npx + ipx
      ipz_msrc  = mod(ipz + np_supercell - 1/npz  + npz, npz)*npy*npx + ipy*npx + ipx

    !     The number of boundary cells in destination process is irregular
    !     when the number of cells is the product of power-of-two number and
    !     power-of-three number. The position within merging cell of first cell
    !     or last cell in a process is used to determine the boundary size.
    !     Irregular communication pairing is realized by calculating residual
    !     boundary cells which are not received yet on self rank and destination rank.


    call preprocess_comm(ipz,nbound_zm,nbound_zp, &
     &                   npz,nlz,nscellz,nsczdiv, &
     &                   nb_destm,nb_destp,nitr)

   !npow = nlz      ! number of merging cells on this FMM level.

   !if(npow == 2) then  ! cell merge number for upper level is two.

       !  ...   number of super-cells always has factor 2 when number of merging cells for
       !  ...   upper level is 2. boundary size in the destination process is always own
       !  ...   plus-side boundary size for plus direction communication, or is minus-side
       !  ...   boundary size for minus direction communication.

   !   nb_destm = nbound_zm      ! plus boundary cells on neighbor rank in minus-direction.
   !   nb_destp = nbound_zp      ! minus boundary cells on neighbor rank in plus-direction.

   !else                ! cell merge number for upper level is three.

       !  ...   number of super-cells is a multiple of power of 3 number.
       !  ...   this means super-cells consist of multiple set of merging 3 cells.
       !  ...   when the position within merging cells of the last cell on neighbor rank
       !  ...   is 1st, 2nd, or 3rd the number of boundary cells is 5, 4, or 3, respectively.
       !  ...   when the position within merging cells of the 1st cell on neighbor rank
       !  ...   is 1st, 2nd, or 3rd the number of boundary cells is 3, 4, or 5, respectively.

   !   iczg0 = (nscellz * ipz) / npz + 1            ! beginning global cell address of this rank.
   !   iczg1 = (nscellz * (ipz + 1) - 1) / npz + 1  ! final global cell ddress of this rank.
       ! plus boundary cells on neighbor rank in minus direction.
   !   nb_destm = mtbdp2 - mod(iczg0 + 1, npow)
       ! minus boundary cells on neighbor rank in plus direction.
   !   nb_destp = mtbd   + mod(iczg1    , npow)

   !endif
    ! the number of bucket-relay repetition. boundary cells devided by cells per process.
   !nitr = max( (nbound_zm - 1) / nsczdiv + 1,  (nbound_zp - 1) / nsczdiv + 1    )
   !nitr = max( (nb_destm - 1) / nsczdiv + 1, nitr)
   !nitr = max( (nb_destp - 1) / nsczdiv + 1, nitr)

#ifdef DBG_COMLL
    write(myrank+NBFU,*) "+Z: ipz_pdest,ipz_psrc,np_supercell: ", &
         &            ipz_pdest,ipz_psrc,np_supercell
    write(myrank+NBFU,*) "-Z: ipz_mdest,ipz_msrc,np_supercell: ", &
         &            ipz_mdest,ipz_msrc,np_supercell
    write(myrank+NBFU,*) "+-Z: nitr;",nitr
    if(npow /= 2) write(myrank+NBFU,*) "+-Z: iczg0,iczg1,nscellz;",iczg0,iczg1, nscellz
    write(myrank+NBFU,*) "+-Z: nb_destp ; ",nb_destp
    write(myrank+NBFU,*) "+-Z: nb_destm ; ",nb_destm
#endif

    do itr = 1, nitr
       ! the note "p" or "m" denotes plus or minus direction communication, respectively.
       ! the note "b" denotes cell local cell address of receiving cells.
       if (itr == 1) then        ! first iteration
!Z+
          iczp0 = nbound_zm + nsczdiv - nb_destp + 1      ! beginning local cell address of sending cells.
          if(iczp0 < nbound_zm + 1) iczp0 = nbound_zm + 1 ! boundary cells > local cells.
          iczp1 = nbound_zm + nsczdiv                     ! final local cell address of sending cells.
          iczbp0 = nbound_zm - nsczdiv + 1                ! beginning local cell address of receiving cells.
          if(iczbp0 < 1) iczbp0 = 1                       ! cells per process is larger than boundary cells.
          iczbp1 = nbound_zm                              ! final local cell address of receiving cells.
          nccps = 2*nscxdiv*nscydiv*(iczp1 - iczp0 + 1)*mylm   ! the number of buffered sending cells.
          nccpr = 2*nscxdiv*nscydiv*(iczbp1 - iczbp0 + 1)*mylm ! the number of buffered receiving cells.

          ncc = 0
          do icz = iczp1, iczp0, -1
             do icx = nbound_xm + 1, nbound_xm + nscxdiv
                do icy = nbound_ym + 1, nbound_ym + nscydiv
                   do m = 1, mylm
                      ncc = ncc + 1
                      ccbufp(ncc       ) = wm(m, 1, icz, icy, icx )
                      ccbufp(ncc + mylm) = wm(m, 2, icz, icy, icx )
                   end do
                   ncc = ncc + mylm
                end do
             end do
          end do
!Z-
          iczm0 = nbound_zm + 1                           ! beginning local cell address of sending cells.
          iczm1 = nbound_zm + nb_destm                    ! final local cell address of sending cells.
          if(iczm1 > nbound_zm + nsczdiv) iczm1 = nbound_zm + nsczdiv  ! boundary cells > local cells.
          iczbm0 = nbound_zm + nsczdiv + 1                ! beginning local cell address of receiving cells.
          iczbm1 = nbound_zm + 2*nsczdiv                  ! final local cell address of receiving cells.
          if(iczbm1 > nbound_zm + nsczdiv + nbound_zp) &  ! cells per process is larger than boundary cells.
               &         iczbm1 = nbound_zm + nsczdiv + nbound_zp
          nccms = 2*nscxdiv*nscydiv*(iczm1 - iczm0 + 1)*mylm   ! the number of buffered sending cells.
          nccmr = 2*nscxdiv*nscydiv*(iczbm1 - iczbm0 + 1)*mylm ! the number of buffered receiving cells.

          ncc = 0
          do icz = iczm0, iczm1
             do icx = nbound_xm + 1, nbound_xm + nscxdiv
                do icy = nbound_ym + 1, nbound_ym + nscydiv
                   do m = 1, mylm
                      ncc = ncc + 1
                      ccbufm(ncc       ) = wm(m, 1, icz, icy, icx )
                      ccbufm(ncc + mylm) = wm(m, 2, icz, icy, icx )
                   end do
                   ncc = ncc + mylm
                end do
             end do
          end do

#ifdef MCOM_TIMER
          te = MPI_WTIME()
          tcomll_lcl = tcomll_lcl + te - ts
          call mpi_barrier(mpi_comm_world,ierr)
          ts = MPI_WTIME()
          tcomll_pbar = tcomll_pbar + ts - te
#endif

#ifdef PROFILE_COMM
          call mpi_barrier(mpi_comm_world, ierr)
          call start_profile(comm_fmm_local_multi_dr_id)
#endif
#ifdef FJ_RDMA
          call rdma_post(ipz_psrc)
          if (ipz_psrc .ne. ipz_msrc) then
             call rdma_post(ipz_msrc)
             call rdma_wait(ipz_mdest)
          endif
          call rdma_wait(ipz_pdest)
          !
          if (ipz_pdest .ne. ipz_mdest) then
             call rdma_put_post(ipz_pdest, heap_raddr(ipz_pdest+1)+rccbufp_offset, heap_laddr+ccbufp_offset, nccps*ids)
             call rdma_put_post(ipz_mdest, heap_raddr(ipz_mdest+1)+rccbufm_offset, heap_laddr+ccbufm_offset, nccms*ids)
             call rdma_wait(ipz_psrc)
             call rdma_wait(ipz_msrc)
          else
             call rdma_put     (ipz_pdest, heap_raddr(ipz_pdest+1)+rccbufp_offset, heap_laddr+ccbufp_offset, nccps*ids)
             call rdma_put_post(ipz_mdest, heap_raddr(ipz_mdest+1)+rccbufm_offset, heap_laddr+ccbufm_offset, nccms*ids)
             call rdma_wait(ipz_msrc)
          endif
#else
          call mpi_irecv(rccbufp(1,1), nccpr, &
               &              MPI_PREC_MOMVEC_COM, ipz_psrc,   1, &
               &              mpi_comm_world, irq(1), ierr)
          call mpi_isend(ccbufp, nccps, &
               &              MPI_PREC_MOMVEC_COM, ipz_pdest,  1, &
               &              mpi_comm_world, irq(2), ierr)
          call mpi_irecv(rccbufm(1,1), nccmr, &
               &              MPI_PREC_MOMVEC_COM, ipz_msrc,   2, &
               &              mpi_comm_world, irq(3), ierr)
          call mpi_isend(ccbufm, nccms, &
               &              MPI_PREC_MOMVEC_COM, ipz_mdest,  2, &
               &              mpi_comm_world, irq(4), ierr)
          nrq = 4
          call mpi_waitall(nrq, irq, istatus, ierr)
#endif
#ifdef PROFILE_COMM
          call end_profile(comm_fmm_local_multi_dr_id)
#endif
          
#ifdef MCOM_TIMER
          te = MPI_WTIME()
          tcomll_com = tcomll_com + te - ts
          nccomll_com = nccomll_com + 1
#ifdef MCOM_TIMER_DETAIL
          write(myrank+NBFU,*) "Rank: ", myrank, &
               &                    " TcomLL-Z-1st: ", te - ts
#endif
          ts = te
#endif

          ncc = 0
          do iczb = iczbp1, iczbp0, -1
             do icx = nbound_xm + 1, nbound_xm + nscxdiv
                do icy = nbound_ym + 1, nbound_ym + nscydiv
                   do m = 1, mylm
                      ncc = ncc + 1
                      wm(m,1,iczb,icy,icx) = rccbufp(ncc,       1)
                      wm(m,2,iczb,icy,icx) = rccbufp(ncc + mylm,1)
                   end do
                   ncc = ncc + mylm
                end do
             end do
          end do

          ncc = 0
          do iczb = iczbm0, iczbm1
             do icx = nbound_xm + 1, nbound_xm + nscxdiv
                do icy = nbound_ym + 1, nbound_ym + nscydiv
                   do m = 1, mylm
                      ncc = ncc + 1
                      wm(m,1,iczb,icy,icx) = rccbufm(ncc,       1)
                      wm(m,2,iczb,icy,icx) = rccbufm(ncc + mylm,1)
                   end do
                   ncc = ncc + mylm
                end do
             end do
          end do

#ifdef DBG_COMLL
          write(myrank+NBFU,*) "+Z 1st        : ", &
               &            "iczp0,iczp1,iczbp0,iczbp1: ", &
               &             iczp0,iczp1,iczbp0,iczbp1
          write(myrank+NBFU,*) "+Z 1st        : ", &
               &            "nb_destp, nccps, nccpr: ", &
               &             nb_destp, nccps, nccpr
          write(myrank+NBFU,*) "-Z 1st        : ", &
               &            "iczm0,iczm1,iczbm0,iczbm1: ", &
               &             iczm0,iczm1,iczbm0,iczbm1
          write(myrank+NBFU,*) "-Z 1st        : ", &
               &            "nb_destm, nccms, nccmr: ", &
               &             nb_destm, nccms, nccmr
#endif

       else      ! iteration follows == irregular&regular pairing. =

          ibs = mod(itr, 2) + 1            ! send buffer index.
          ibr = mod(itr + 1, 2) + 1        ! receive buffer index.
          iczp0 = iczbp0                   ! beginning local cell address of sending cells. which is
                                           ! equal to that of receiving cells in previous iteration.
          if(nsczdiv*itr > nb_destp) &     ! destination boundary cells < total bucket-relay cells.
               &            iczp0 = nbound_zm + nsczdiv - nb_destp + 1
          iczp1 = iczbp1                   ! final local cell address of sending cells. which is
                                           ! equal to that of receiving cells in previous iteration.
          iczbp0 = nbound_zm - nsczdiv*itr + 1   ! beginning local cell address of receiving cells.
          if(iczbp0 < 1) iczbp0 = 1              ! necessary cells is smaller than cell per process.
          iczbp1 = nbound_zm - nsczdiv*(itr - 1) ! final local cell address of receiving cells.
          nccps = 2*nscxdiv*nscydiv*(iczp1 - iczp0 + 1)*mylm    ! the number of buffered sending cells.
          nccpr = 2*nscxdiv*nscydiv*(iczbp1 - iczbp0 + 1)*mylm  ! the number of buffered sending cells.

          iczm0 = iczbm0                   ! beginning local cell address of sending cells. which is
                                           ! equal to that of receiving cells in previous iteration.
          iczm1 = iczbm1                   ! final local cell address of sending cells. which is
                                           ! equal to that of receiving cells in previous iteration.
          if(nsczdiv*itr > nb_destm) &     ! destination boundary cells < total bucket-relay cells.
               &            iczm1 = nbound_zm + nb_destm
          iczbm0 = nbound_zm + nsczdiv*itr + 1   ! beginning local cell address of receiving cells.
          iczbm1 = nbound_zm + nsczdiv*(itr + 1) ! final local cell address of receiving cells.
          if(iczbm1 > nbound_zm + nsczdiv + nbound_zp) &      ! necessary cells < cell per process.
               &            iczbm1 = nbound_zm + nsczdiv + nbound_zp
          nccms = 2*nscxdiv*nscydiv*(iczm1 - iczm0 + 1)*mylm    ! the number of buffered sending cells.
          nccmr = 2*nscxdiv*nscydiv*(iczbm1 - iczbm0 + 1)*mylm  ! the number of buffered receiving cells.

#ifdef DBG_COMLL
          write(myrank+NBFU,*) "+Z 2nd & later: itr; ",itr, &
               &              " iczp0,iczp1: ",iczp0,iczp1," nccps: ",nccps
          write(myrank+NBFU,*) "+Z 2nd & later: itr; ",itr, &
               &              " iczbp0,iczbp1: ",iczbp0,iczbp1," nccpr: ",nccpr
          write(myrank+NBFU,*) "-Z 2nd & later: itr; ",itr, &
               &              " iczm0,iczm1: ",iczm0,iczm1," nccms: ",nccms
          write(myrank+NBFU,*) "-Z 2nd & later: itr; ",itr, &
               &              " iczbm0,iczbm1: ",iczbm0,iczbm1," nccmr: ",nccmr
#endif
#ifdef MCOM_TIMER
          te = MPI_WTIME()
          tcomll_lcl = tcomll_lcl + te - ts
          call mpi_barrier(mpi_comm_world,ierr)
          ts = MPI_WTIME()
          tcomll_pbar = tcomll_pbar + ts - te
#endif

#ifdef PROFILE_COMM
          call mpi_barrier(mpi_comm_world, ierr)
          call start_profile(comm_fmm_local_multi_dr_id)
#endif
#ifdef FJ_RDMA
          if(nccpr > 0) call rdma_post(ipz_psrc)
          if(nccmr > 0) call rdma_post(ipz_msrc)
          if(nccps > 0) call rdma_wait(ipz_pdest)
          if(nccms > 0) call rdma_wait(ipz_mdest)
          !
          if(nccps > 0) then
             call rdma_put_post(ipz_pdest, heap_raddr(ipz_pdest+1)+rccbufp_offset+(ibr-1)*(mylm0*nbsize*nscydiv*nscxdiv)*ibds, &
                  & heap_laddr+rccbufp_offset+(ibs-1)*(mylm0*nbsize*nscydiv*nscxdiv)*ibds, nccps*ids)
          end if
          if(nccms > 0) then
             call rdma_put_post(ipz_mdest, heap_raddr(ipz_mdest+1)+rccbufm_offset+(ibr-1)*(mylm0*nbsize*nscydiv*nscxdiv)*ibds, &
                  & heap_laddr+rccbufm_offset+(ibs-1)*(mylm0*nbsize*nscydiv*nscxdiv)*ibds, nccms*ids)
          end if
          !
          if(nccpr > 0) call rdma_wait(ipz_psrc)
          if(nccmr > 0) call rdma_wait(ipz_msrc)
#else
          nrq = 0
          if(nccpr > 0) then
             nrq = nrq + 1
             call mpi_irecv(rccbufp(1,ibr), nccpr, &
                  &                   MPI_PREC_MOMVEC_COM, ipz_psrc,     1, &
                  &                   mpi_comm_world, irq(nrq), ierr)

#ifdef DBG_COMLL
             write(myrank+NBFU,*) "+Z 2nd & later: recipient itr; ",itr, &
                  &              " iczbp0,iczbp1: ",iczbp0,iczbp1," nccpr: ",nccpr
#endif

          endif     ! recipient (p)
          if(nccps > 0) then
             nrq = nrq + 1
             call mpi_isend(rccbufp(1,ibs), nccps, &
                  &                   MPI_PREC_MOMVEC_COM, ipz_pdest,    1, &
                  &                   mpi_comm_world, irq(nrq), ierr)

#ifdef DBG_COMLL
             write(myrank+NBFU,*) "+Z 2nd & later: sender itr; ",itr, &
                  &              " iczp0,iczp1: ",iczp0,iczp1," nccps: ",nccps
#endif

          endif     ! sender (p)
          if(nccmr > 0) then
             nrq = nrq + 1
             call mpi_irecv(rccbufm(1,ibr), nccmr, &
                  &                   MPI_PREC_MOMVEC_COM, ipz_msrc,     2, &
                  &                   mpi_comm_world, irq(nrq), ierr)

#ifdef DBG_COMLL
             write(myrank+NBFU,*) "-Z 2nd & later: recipient itr; ",itr, &
                  &              " iczbm0,iczbm1: ",iczbm0,iczbm1," nccmr: ",nccmr
#endif

          endif     ! recipient (m)
          if(nccms > 0) then
             nrq = nrq + 1
             call mpi_isend(rccbufm(1,ibs), nccms, &
                  &                   MPI_PREC_MOMVEC_COM, ipz_mdest,    2, &
                  &                   mpi_comm_world, irq(nrq), ierr)

#ifdef DBG_COMLL
             write(myrank+NBFU,*) "-Z 2nd & later: sender itr; ",itr, &
                  &              " iczm0,iczm1: ",iczm0,iczm1," nccms: ",nccms
#endif

          endif     ! sender (m)

          call mpi_waitall(nrq, irq, istatus, ierr)
#endif
#ifdef PROFILE_COMM
          call end_profile(comm_fmm_local_multi_dr_id)
#endif
#ifdef MCOM_TIMER
          te = MPI_WTIME()
          tcomll_com = tcomll_com + te - ts
          nccomll_com = nccomll_com + 1
#ifdef MCOM_TIMER_DETAIL
          write(myrank+NBFU,*) "Rank: ", myrank, &
               &                    " TcomLL-Z-2nd&later: ", te - ts
#endif
         ts = te
#endif

          if(nccpr > 0) then
             ncc = 0
             do iczb = iczbp1, iczbp0, -1
                do icx = nbound_xm + 1, nbound_xm + nscxdiv
                   do icy = nbound_ym + 1, nbound_ym + nscydiv
                      do m = 1, mylm
                         ncc = ncc + 1
                         wm(m,1,iczb,icy,icx) = rccbufp(ncc,       ibr)
                         wm(m,2,iczb,icy,icx) = rccbufp(ncc + mylm,ibr)
                      end do
                      ncc = ncc + mylm
                   end do
                end do
             end do
          endif     ! recipient (p)

          if(nccmr > 0) then
             ncc = 0
             do iczb = iczbm0, iczbm1
                do icx = nbound_xm + 1, nbound_xm + nscxdiv
                   do icy = nbound_ym + 1, nbound_ym + nscydiv
                      do m = 1, mylm
                         ncc = ncc + 1
                         wm(m,1,iczb,icy,icx) = rccbufm(ncc,       ibr)
                         wm(m,2,iczb,icy,icx) = rccbufm(ncc + mylm,ibr)
                      end do
                      ncc = ncc + mylm
                   end do
                end do
             end do
          endif     ! recipient (m)

          !            endif     ! irregular&regular pairing.

       endif   ! iteration
    end do     ! iteration

    ! LLmoment +Y
    np_supercell = (npy - 1) / nscelly + 1
    ipy_pdest = ipz*npy*npx + mod(ipy + np_supercell - 1/npy + npy, npy)*npx + ipx
    ipy_psrc = ipz*npy*npx + mod(ipy - np_supercell + 1/npy + npy, npy)*npx + ipx
    ! LLmoment -Y
    ipy_mdest = ipz*npy*npx + mod(ipy - np_supercell + 1/npy + npy, npy)*npx + ipx
    ipy_msrc = ipz*npy*npx + mod(ipy + np_supercell - 1/npy + npy, npy)*npx + ipx

    !     The number of boundary cells in destination process is irregular
    !     when number of cells is product of power-of-two number and
    !     power-of-three number. The position within merging cell of first cell
    !     or last cell in a process is used to know the boundary size.
    !     Irregular communication pairing is realized by calculating residual
    !     boundary cells which are not received yet on self rank and
    !     destination rank.

! input    : ipy,nbound_ym,nbound_yp
! constant : npy,nly,nscelly,nscydiv
! output   : nb_destm,nb_destp,nitr
    call preprocess_comm(ipy,nbound_ym,nbound_yp, &
     &                   npy,nly,nscelly,nscydiv, &
     &                   nb_destm,nb_destp,nitr)    

   !npow = nly

   !if(npow == 2) then  ! cell merge number for upper level is two.

       !  ...   number of super-cells always has factor 2 when number of merging cells for
       !  ...   upper level is 2. boundary size in the destination process is always own
       !  ...   plus-side boundary size for plus direction communication or minus-side
       !  ...   boundary size for minus direction communication.
   !   nb_destm = nbound_ym
   !   nb_destp = nbound_yp

   !else                ! cell merge number for upper level is three.

       !  ...   number of super-cells is a multiple of power of 3 number.
       !  ...   this means super-cells consist of multiple set of merging 3 cells.
       !  ...   when the position within merging cells of the last cell on neighbor rank
       !  ...   is 1st, 2nd, or 3rd the number of boundary cells is 5, 4, or 3, respectively.
       !  ...   when the position within merging cells of the 1st cell on neighbor rank
       !  ...   is 1st, 2nd, or 3rd the number of boundary cells is 3, 4, or 5, respectively.

   !   icyg0 = (nscelly * ipy) / npy + 1
   !   icyg1 = (nscelly * (ipy + 1) - 1) / npy + 1
   !   nb_destm = mtbdp2 - mod(icyg0 + 1, npow)
   !   nb_destp = mtbd   + mod(icyg1    , npow)

   !endif

   !nitr = max( (nbound_ym - 1) / nscydiv + 1, (nbound_yp - 1) / nscydiv + 1 )
   !nitr = max( (nb_destm - 1) / nscydiv + 1, nitr)
   !nitr = max( (nb_destp - 1) / nscydiv + 1, nitr)

#ifdef DBG_COMLL
    write(myrank+NBFU,*) "+Y: ipy_pdest,ipy_psrc,np_supercell: ", &
         &            ipy_pdest,ipy_psrc,np_supercell
    write(myrank+NBFU,*) "+Y: nitr;",nitr
    write(myrank+NBFU,*) "-Y: ipy_pdest,ipy_psrc,np_supercell: ", &
         &            ipy_mdest,ipy_msrc,np_supercell
    write(myrank+NBFU,*) "-Y: nitr;",nitr
    if(npow /= 2) write(myrank+NBFU,*) "+-Y: icyg0,icyg1,nscelly;",icyg0,icyg1, &
         &                    nscelly
    write(myrank+NBFU,*) "+-Y: nb_destp ; ",nb_destp
    write(myrank+NBFU,*) "+-Y: nb_destm ; ",nb_destm
#endif

    do icx = nbound_xm + 1, nbound_xm + nscxdiv
       do itr = 1, nitr
          if (itr == 1) then                ! first iteration
!
             icyp0 = nbound_ym + nscydiv - nb_destp + 1
             if(icyp0 < nbound_ym + 1) icyp0 = nbound_ym + 1
             icyp1 = nbound_ym + nscydiv
             icybp0 = nbound_ym - nscydiv + 1
             if(icybp0 < 1) icybp0 = 1
             icybp1 = nbound_ym
#ifdef DEBUG_RDMAYA
!debug y+,1
!            icybp0=icbp0_test(il)%iaddr(itr)  ! OK
!            icybp1=icbp1_test(il)%iaddr(itr)  ! OK
!debug
#endif
             nccps = 2*(nbound_zm + nsczdiv + nbound_zp)*(icyp1 - icyp0 + 1)*mylm
             nccpr = 2*(nbound_zm + nsczdiv + nbound_zp)*(icybp1 - icybp0 + 1)*mylm
!Y-
             icym0 = nbound_ym + 1
             icym1 = nbound_ym + nb_destm
             if(icym1 > nbound_ym + nscydiv) &
                  &                 icym1 = nbound_ym + nscydiv
             icybm0 = nbound_ym + nscydiv + 1
             icybm1 = nbound_ym + 2*nscydiv
             if(icybm1 > nbound_ym + nscydiv + nbound_yp) &
                  &                 icybm1 = nbound_ym + nscydiv + nbound_yp
             nccms = 2*(nbound_zm + nsczdiv + nbound_zp) &
                  &                *(icym1 - icym0 + 1)*mylm
             nccmr = 2*(nbound_zm + nsczdiv + nbound_zp) &
                  &                *(icybm1 - icybm0 + 1)*mylm
#ifdef DEBUG_RDMAYA
!debug y-,1
!            icybm0=icbm0_test(il)%iaddr(itr)  ! OK
!            icybm1=icbm1_test(il)%iaddr(itr)  ! OK
!debug
#endif

#ifdef MCOM_TIMER
             te = MPI_WTIME()
             tcomll_lcl = tcomll_lcl + te - ts
             call mpi_barrier(mpi_comm_world,ierr)
             ts = MPI_WTIME()
             tcomll_pbar = tcomll_pbar + ts - te
#endif
#ifdef DEBUG_RDMAYA
! nakao
 call mpi_irecv(hoge,                   1, MPI_INTEGER, ipy_psrc,  1, &
 mpi_comm_world, irq(1), ierr)
 call mpi_isend(icbp0_y(il)%iaddr(itr), 1, MPI_INTEGER, ipy_pdest, 1, &
 mpi_comm_world, irq(2), ierr)
 call mpi_waitall(2, irq, istatus, ierr)
 if (hoge .ne. icybp0) then
   write(*,*) "Error",hoge,icybp0,itr,' ,il',il
 endif
! nakao
#endif
#ifdef PROFILE_COMM
             call mpi_barrier(mpi_comm_world, ierr)
             call start_profile(comm_fmm_local_multi_dr_id)
#endif
#ifdef FJ_RDMA
             call rdma_post(ipy_psrc)
             if (ipy_psrc .ne. ipy_msrc) then
                call rdma_post(ipy_msrc)
                call rdma_wait(ipy_mdest)
             endif
             call rdma_wait(ipy_pdest)
             !
             if (ipy_pdest .ne. ipy_mdest) then
                call rdma_put_post(ipy_pdest, &
                wm_raddr(ipy_pdest+1) + ((icx-1)*wm_4th_offset + (icbp0_y(il)%iaddr(itr)-1)*wm_3th_offset)*ibds, &
                wm_laddr              + ((icx-1)*wm_4th_offset + (icyp0-1)                 *wm_3th_offset)*ibds, &
                nccps*ids)
                call rdma_put_post(ipy_mdest, &
                wm_raddr(ipy_mdest+1) + ((icx-1)*wm_4th_offset + (icbm0_y(il)%iaddr(itr)-1)*wm_3th_offset)*ibds, &
                wm_laddr              + ((icx-1)*wm_4th_offset + (icym0-1)                 *wm_3th_offset)*ibds, &
                nccms*ids)
                call rdma_wait(ipy_psrc)
                call rdma_wait(ipy_msrc)
             else
                call rdma_put(ipy_pdest, &
                wm_raddr(ipy_pdest+1) + ((icx-1)*wm_4th_offset + (icbp0_y(il)%iaddr(itr)-1)*wm_3th_offset)*ibds, &
                wm_laddr              + ((icx-1)*wm_4th_offset + (icyp0-1)                 *wm_3th_offset)*ibds, &
                nccps*ids)
                call rdma_put_post(ipy_mdest, &
                wm_raddr(ipy_mdest+1) + ((icx-1)*wm_4th_offset + (icbm0_y(il)%iaddr(itr)-1)*wm_3th_offset)*ibds, &
                wm_laddr              + ((icx-1)*wm_4th_offset + (icym0-1)                 *wm_3th_offset)*ibds, &
                nccms*ids)
                call rdma_wait(ipy_msrc)
             end if
#else
             call mpi_irecv(wm(1,1,1,icybp0,icx), nccpr, &
                  &                MPI_PREC_MOMVEC_COM, ipy_psrc,     1, &
                  &                mpi_comm_world, irq(1), ierr)
             call mpi_isend(wm(1,1,1,icyp0,icx), nccps, &
                  &                MPI_PREC_MOMVEC_COM, ipy_pdest,    1, &
                  &                mpi_comm_world, irq(2), ierr)
             call mpi_irecv(wm(1,1,1,icybm0,icx), nccmr, &
                  &                MPI_PREC_MOMVEC_COM, ipy_msrc,     2, &
                  &                mpi_comm_world, irq(3), ierr)
             call mpi_isend(wm(1,1,1,icym0,icx), nccms, &
                  &                MPI_PREC_MOMVEC_COM, ipy_mdest,    2, &
                  &                mpi_comm_world, irq(4), ierr)
             nrq = 4
             call mpi_waitall(nrq, irq, istatus, ierr)
#endif
#ifdef PROFILE_COMM
             call end_profile(comm_fmm_local_multi_dr_id)
#endif
#ifdef MCOM_TIMER
             te = MPI_WTIME()
             tcomll_com = tcomll_com + te - ts
             nccomll_com = nccomll_com + 1
#ifdef MCOM_TIMER_DETAIL
             write(myrank+NBFU,*) "Rank: ", myrank, &
                  &                    " TcomLL-Y-1st: ", te - ts
#endif
             ts = te
#endif
#ifdef DBG_COMLL
             write(myrank+NBFU,*) "+Y 1st        : ", &
                  &                    "icyp0,icyp1,icybp0,icybp1: ", &
                  &                     icyp0,icyp1,icybp0,icybp1
             write(myrank+NBFU,*) "+Y 1st        : ", &
                  &            "nb_destp, nccps, nccpr: ", &
                  &             nb_destp, nccps, nccpr
             write(myrank+NBFU,*) "-Y 1st        : ", &
                  &                    "icym0,icym1,icybm0,icybm1: ", &
                  &                     icym0,icym1,icybm0,icybm1
             write(myrank+NBFU,*) "-Y 1st        : ", &
                  &            "nb_destm, nccms, nccmr: ", &
                  &             nb_destm, nccms, nccmr
#endif

          else      ! iteration follows == irregular&regular pairing. =
!
             icyp0 = icybp0
             if(nscydiv*itr > nb_destp) &
                  &               icyp0 = nbound_ym + nscydiv - nb_destp + 1
             icyp1 = icybp1
             icybp0 = nbound_ym - nscydiv*itr + 1
             if(icybp0 < 1) icybp0 = 1
             icybp1 = nbound_ym - nscydiv*(itr - 1)
             nccps = 2*(nbound_zm + nsczdiv + nbound_zp)*(icyp1 - icyp0 + 1)*mylm
             nccpr = 2*(nbound_zm + nsczdiv + nbound_zp)*(icybp1 - icybp0 + 1)*mylm
!Y-
             icym0 = icybm0
             icym1 = icybm1
             if(nscydiv*itr > nb_destm) &
                  &               icym1 = nbound_ym + nb_destm
             icybm0 = nbound_ym + nscydiv*itr + 1
             icybm1 = nbound_ym + nscydiv*(itr + 1)
             if(icybm1 > nbound_ym + nscydiv + nbound_yp) &
                  &               icybm1 = nbound_ym + nscydiv + nbound_yp
             nccms = 2*(nbound_zm + nsczdiv + nbound_zp)*(icym1 - icym0 + 1)*mylm
             nccmr = 2*(nbound_zm + nsczdiv + nbound_zp)*(icybm1 - icybm0 + 1)*mylm

#ifdef DBG_COMLL
             write(myrank+NBFU,*) "+Y 2nd & later: itr: ",itr, &
                  &              " icyp0,icyp1: ",icyp0,icyp1," nccps: ",nccps
             write(myrank+NBFU,*) "+Y 2nd & later: itr: ",itr, &
                  &              " icybp0,icybp1: ",icybp0,icybp1," nccpr: ",nccpr
             write(myrank+NBFU,*) "-Y 2nd & later: itr: ",itr, &
                  &              " icym0,icym1: ",icym0,icym1," nccms: ",nccms
             write(myrank+NBFU,*) "-Y 2nd & later: itr: ",itr, &
                  &              " icybm0,icybm1: ",icybm0,icybm1," nccmr: ",nccmr
#endif
#ifdef MCOM_TIMER
             te = MPI_WTIME()
             tcomll_lcl = tcomll_lcl + te - ts
             call mpi_barrier(mpi_comm_world,ierr)
             ts = MPI_WTIME()
             tcomll_pbar = tcomll_pbar + ts - te
#endif
#ifdef PROFILE_COMM
             call mpi_barrier(mpi_comm_world, ierr)
             call start_profile(comm_fmm_local_multi_dr_id)
#endif
#ifdef FJ_RDMA
             if(nccpr > 0) call rdma_post(ipy_psrc)
             if(nccmr > 0) call rdma_post(ipy_msrc)
             if(nccps > 0) call rdma_wait(ipy_pdest)
             if(nccms > 0) call rdma_wait(ipy_mdest)
             !
             if(nccps > 0) then
                call rdma_put_post(ipy_pdest, &
                     wm_raddr(ipy_pdest+1) + ((icx-1)*wm_4th_offset + (icbp0_y(il)%iaddr(itr)-1)*wm_3th_offset)*ibds, &
                     wm_laddr              + ((icx-1)*wm_4th_offset + (icyp0-1)                 *wm_3th_offset)*ibds, &
                     nccps*ids)
             end if
             if(nccms > 0) then
                call rdma_put_post(ipy_mdest, &
                     wm_raddr(ipy_mdest+1) + ((icx-1)*wm_4th_offset + (icbm0_y(il)%iaddr(itr)-1)*wm_3th_offset)*ibds, &
                     wm_laddr              + ((icx-1)*wm_4th_offset + (icym0-1)                 *wm_3th_offset)*ibds, &
                     nccms*ids)
             end if
             !
             if(nccpr > 0) call rdma_wait(ipy_psrc)
             if(nccmr > 0) call rdma_wait(ipy_msrc)
#else
             nrq = 0
             if(nccpr > 0) then
                nrq = nrq + 1
                call mpi_irecv(wm(1,1,1,icybp0,icx), nccpr, &
                     &                      MPI_PREC_MOMVEC_COM, ipy_psrc,   1, &
                     &                      mpi_comm_world, irq(nrq), ierr)

#ifdef DBG_COMLL
                write(myrank+NBFU,*) "+Y 2nd & later: recipient itr: ",itr, &
                     &              " icybp0,icybp1: ",icybp0,icybp1," nccpr: ",nccpr
#endif

             endif         ! recipient (p)
             if(nccps > 0) then
                nrq = nrq + 1
                call mpi_isend(wm(1,1,1,icyp0,icx), nccps, &
                     &                      MPI_PREC_MOMVEC_COM, ipy_pdest,  1, &
                     &                      mpi_comm_world, irq(nrq), ierr)

#ifdef DBG_COMLL
                write(myrank+NBFU,*) "+Y 2nd & later: sender itr: ",itr, &
                     &              " icyp0,icyp1: ",icyp0,icyp1," nccps: ",nccps
#endif

             endif     ! sender (p)
             if(nccmr > 0) then
                nrq = nrq + 1
                call mpi_irecv(wm(1,1,1,icybm0,icx), nccmr, &
                     &                      MPI_PREC_MOMVEC_COM, ipy_msrc,   2, &
                     &                      mpi_comm_world, irq(nrq), ierr)

#ifdef DBG_COMLL
                write(myrank+NBFU,*) "-Y 2nd & later: recipient itr: ",itr, &
                     &              " icybm0,icybm1: ",icybm0,icybm1," nccmr: ",nccmr
#endif

             endif         ! recipient (m)
             if(nccms > 0) then
                nrq = nrq + 1
                call mpi_isend(wm(1,1,1,icym0,icx), nccms, &
                     &                      MPI_PREC_MOMVEC_COM, ipy_mdest,  2, &
                     &                      mpi_comm_world, irq(nrq), ierr)

#ifdef DBG_COMLL
                write(myrank+NBFU,*) "-Y 2nd & later: sender itr: ",itr, &
                     &              " icym0,icym1: ",icym0,icym1," nccms: ",nccms
#endif
             endif         ! sender (m)

             call mpi_waitall(nrq, irq, istatus, ierr)
#endif
#ifdef PROFILE_COMM
             call end_profile(comm_fmm_local_multi_dr_id)
#endif
#ifdef MCOM_TIMER
             te = MPI_WTIME()
             tcomll_com = tcomll_com + te - ts
             nccomll_com = nccomll_com + 1
#ifdef MCOM_TIMER_DETAIL
             write(myrank+NBFU,*) "Rank: ", myrank, &
                  &                    " TcomLL-Y-2nd&later: ", te - ts
#endif
             ts = te
#endif

             !               endif         ! irregular&regular pairing.

          endif         ! iteration
       end do           ! iteration
    end do            ! ipx

    ! LLmoment +X
    np_supercell = (npx - 1) / nscellx + 1
    ipx_pdest = ipz*npy*npx + ipy*npx + mod(ipx + np_supercell - 1/npx + npx, npx)
    ipx_psrc = ipz*npy*npx + ipy*npx + mod(ipx - np_supercell + 1/npx + npx, npx)
    ! LLmoment -X
    ipx_mdest = ipz*npy*npx + ipy*npx + mod(ipx - np_supercell + 1/npx + npx, npx)
    ipx_msrc = ipz*npy*npx + ipy*npx + mod(ipx + np_supercell - 1/npx + npx, npx)

    !     The number of boundary cells in destination process is irregular
    !     when number of cells is product of power-of-two number and
    !     power-of-three number. The position within merging cell of first cell
    !     or last cell in a process is used to know the boundary size.
    !     Irregular communication pairing is realized by calculating residual
    !     boundary cells which are not received yet on self rank and
    !     destination rank.

    call preprocess_comm(ipx,nbound_xm,nbound_xp, &
     &                   npx,nlx,nscellx,nscxdiv, &
     &                   nb_destm,nb_destp,nitr)

   !npow = nlx

   !if(npow == 2) then  ! cell merge number for upper level is two.

       !  ...   number of super-cells always has factor 2 when number of merging cells for
       !  ...   upper level is 2. boundary size in the destination process is always own
       !  ...   plus-side boundary size for plus direction communication or minus-side
       !  ...   boundary size for minus direction communication.
   !   nb_destm = nbound_xm
   !   nb_destp = nbound_xp

   !else                ! cell merge number for upper level is three.

       !  ...   number of super-cells is a multiple of power of 3 number.
       !  ...   this means super-cells consist of multiple set of merging 3 cells.
       !  ...   when the position within merging cells of the last cell on neighbor rank
       !  ...   is 1st, 2nd, or 3rd the number of boundary cells is 5, 4, or 3, respectively.
       !  ...   when the position within merging cells of the 1st cell on neighbor rank
       !  ...   is 1st, 2nd, or 3rd the number of boundary cells is 3, 4, or 5, respectively.

   !   icxg0 = (nscellx * ipx) / npx + 1
   !   icxg1 = (nscellx * (ipx + 1) - 1) / npx + 1
   !   nb_destm = mtbdp2 - mod(icxg0 + 1, npow)
   !   nb_destp = mtbd   + mod(icxg1    , npow)

   !endif

   !nitr = max( (nbound_xm - 1) / nscxdiv + 1, (nbound_xp - 1) / nscxdiv + 1  )
   !nitr = max( (nb_destm - 1) / nscxdiv + 1, nitr)
   !nitr = max( (nb_destp - 1) / nscxdiv + 1, nitr)

#ifdef DBG_COMLL
    write(myrank+NBFU,*) "+X: ipx_pdest,ipx_psrc,np_supercell: ", &
         &            ipx_pdest,ipx_psrc,np_supercell
    write(myrank+NBFU,*) "+X: nitr;",nitr
    write(myrank+NBFU,*) "-X: ipx_mdest,ipx_msrc,np_supercell: ", &
         &            ipx_mdest,ipx_msrc,np_supercell
    write(myrank+NBFU,*) "-X: nitr;",nitr
    if(npow /= 2) write(myrank+NBFU,*) "+-X: icxg0,icxg1,nscellx;",icxg0,icxg1, nscellx
    write(myrank+NBFU,*) "+-X: nb_destp ; ",nb_destp
    write(myrank+NBFU,*) "+-X: nb_destm ; ",nb_destm
#endif

    do itr = 1, nitr
       if (itr == 1) then    ! first iteration
!x+
          icxp0 = nbound_xm + nscxdiv - nb_destp + 1
          if(icxp0 < nbound_xm + 1) icxp0 = nbound_xm + 1
          icxp1 = nbound_xm + nscxdiv
          icxbp0 = nbound_xm - nscxdiv + 1
          if(icxbp0 < 1) icxbp0 = 1
          icxbp1 = nbound_xm
          nccps = 2*(nbound_zm + nsczdiv + nbound_zp)*(nbound_ym + nscydiv + nbound_yp)*(icxp1 - icxp0 + 1)*mylm
          nccpr = 2*(nbound_zm + nsczdiv + nbound_zp)*(nbound_ym + nscydiv + nbound_yp)*(icxbp1 - icxbp0 + 1)*mylm
!x-
          icxm0 = nbound_xm + 1
          icxm1 = nbound_xm + nb_destm
          if(icxm1 > nbound_xm + nscxdiv) &
               &              icxm1 = nbound_xm + nscxdiv
          icxbm0 = nbound_xm + nscxdiv + 1
          icxbm1 = nbound_xm + 2*nscxdiv
          if(icxbm1 > nbound_xm + nscxdiv + nbound_xp) &
               &              icxbm1 = nbound_xm + nscxdiv + nbound_xp
          nccms = 2*(nbound_zm + nsczdiv + nbound_zp)*(nbound_ym + nscydiv + nbound_yp)*(icxm1 - icxm0 + 1)*mylm
          nccmr = 2*(nbound_zm + nsczdiv + nbound_zp)*(nbound_ym + nscydiv + nbound_yp)*(icxbm1 - icxbm0 + 1)*mylm

#ifdef MCOM_TIMER
          te = MPI_WTIME()
          tcomll_lcl = tcomll_lcl + te - ts
          call mpi_barrier(mpi_comm_world,ierr)
          ts = MPI_WTIME()
          tcomll_pbar = tcomll_pbar + ts - te
#endif

#ifdef PROFILE_COMM
          call mpi_barrier(mpi_comm_world, ierr)    
          call start_profile(comm_fmm_local_multi_dr_id)
#endif

#ifdef FJ_RDMA
          call rdma_post(ipx_psrc)
          if (ipx_psrc .ne. ipx_msrc) then
             call rdma_post(ipx_msrc)
             call rdma_wait(ipx_mdest)
          endif
          call rdma_wait(ipx_pdest)
          !
          if (ipx_pdest .ne. ipx_mdest) then
             call rdma_put_post(ipx_pdest, &
                  wm_raddr(ipx_pdest+1) + (icbp0_x(il)%iaddr(itr)-1)*wm_4th_offset*ibds, &
                  wm_laddr              + (icxp0-1) *wm_4th_offset*ibds, &
                  nccps*ids)
             call rdma_put_post(ipx_mdest, &
                  wm_raddr(ipx_mdest+1) + (icbm0_x(il)%iaddr(itr)-1)*wm_4th_offset*ibds, &
                  wm_laddr              + (icxm0-1) *wm_4th_offset*ibds, &
                  nccms*ids)
             call rdma_wait(ipx_psrc)
             call rdma_wait(ipx_msrc)
          else
             call rdma_put(ipx_pdest, &
                  wm_raddr(ipx_pdest+1) + (icbp0_x(il)%iaddr(itr)-1)*wm_4th_offset*ibds, &
                  wm_laddr              + (icxp0-1) *wm_4th_offset*ibds, &
                  nccps*ids)
             call rdma_put_post(ipx_mdest, &
                  wm_raddr(ipx_mdest+1) + (icbm0_x(il)%iaddr(itr)-1)*wm_4th_offset*ibds, &
                  wm_laddr              + (icxm0-1) *wm_4th_offset*ibds, &
                  nccms*ids)
             call rdma_wait(ipx_msrc)
          endif

#else /** not RDMA **/

          call mpi_irecv(wm(1,1,1,1,icxbp0), nccpr, &
               &             MPI_PREC_MOMVEC_COM, ipx_psrc,      1, &
               &             mpi_comm_world, irq(1), ierr)
          call mpi_isend(wm(1,1,1,1,icxp0), nccps, &
               &             MPI_PREC_MOMVEC_COM, ipx_pdest,     1, &
               &             mpi_comm_world, irq(2), ierr)
          call mpi_irecv(wm(1,1,1,1,icxbm0), nccmr, &
               &             MPI_PREC_MOMVEC_COM, ipx_msrc,      2, &
               &             mpi_comm_world, irq(3), ierr)
          call mpi_isend(wm(1,1,1,1,icxm0), nccms, &
               &             MPI_PREC_MOMVEC_COM, ipx_mdest,     2, &
               &             mpi_comm_world, irq(4), ierr)
          nrq = 4
          call mpi_waitall(nrq, irq, istatus, ierr)
#endif /** RDMA or MPI **/

#ifdef PROFILE_COMM
          call end_profile(comm_fmm_local_multi_dr_id)
#endif
#ifdef MCOM_TIMER
          te = MPI_WTIME()
          tcomll_com = tcomll_com + te - ts
          nccomll_com = nccomll_com + 1
#ifdef MCOM_TIMER_DETAIL
          write(myrank+NBFU,*) "Rank: ", myrank, &
               &                    " TcomLL-X-1st: ", tle - tls
#endif
          ts = te
#endif

#ifdef DBG_COMLL
          write(myrank+NBFU,*) "+X 1st        : ", &
               &                    "icxp0,icxp1,icxbp0,icxbp1: ", &
               &                     icxp0,icxp1,icxbp0,icxbp1
          write(myrank+NBFU,*) "+X 1st        : ", &
               &            "nb_destp, nccps, nccpr: ", &
               &             nb_destp, nccps, nccpr
          write(myrank+NBFU,*) "-X 1st        : ", &
               &                    "icxm0,icxm1,icxbm0,icxbm1: ", &
               &                     icxm0,icxm1,icxbm0,icxbm1
          write(myrank+NBFU,*) "-X 1st        : ", &
               &            "nb_destm, nccms, nccmr: ", &
               &             nb_destm, nccms, nccmr
#endif

       else      ! iteration follows == irregular&regular pairing. =
!x+
          icxp0 = icxbp0
          if(nscxdiv*itr > nb_destp) &
               &            icxp0 = nbound_xm + nscxdiv - nb_destp + 1
          icxp1 = icxbp1
          icxbp0 = nbound_xm - nscxdiv*itr + 1
          if(icxbp0 < 1) icxbp0 = 1
          icxbp1 = nbound_xm - nscxdiv*(itr - 1)
          nccps = 2*(nbound_zm + nsczdiv + nbound_zp)*(nbound_ym + nscydiv + nbound_yp)*(icxp1 - icxp0 + 1)*mylm
          nccpr = 2*(nbound_zm + nsczdiv + nbound_zp)*(nbound_ym + nscydiv + nbound_yp)*(icxbp1 - icxbp0 + 1)*mylm
!x-
          icxm0 = icxbm0
          icxm1 = icxbm1
          if(nscxdiv*itr > nb_destm) &
               &            icxm1 = nbound_xm + nb_destm
          icxbm0 = nbound_xm + nscxdiv*itr + 1
          icxbm1 = nbound_xm + nscxdiv*(itr + 1)
          if(icxbm1 > nbound_xm + nscxdiv + nbound_xp) &
               &            icxbm1 = nbound_xm + nscxdiv + nbound_xp
          nccms = 2*(nbound_zm + nsczdiv + nbound_zp)*(nbound_ym + nscydiv + nbound_yp)*(icxm1 - icxm0 + 1)*mylm
          nccmr = 2*(nbound_zm + nsczdiv + nbound_zp)*(nbound_ym + nscydiv + nbound_yp)*(icxbm1 - icxbm0 + 1)*mylm

#ifdef DBG_COMLL
          write(myrank+NBFU,*) "+X 2nd & later: itr: ",itr, &
               &              " icxp0,icxp1: ",icxp0,icxp1," nccps: ",nccps
          write(myrank+NBFU,*) "+X 2nd & later: itr: ",itr, &
               &              " icxbp0,icxbp1: ",icxbp0,icxbp1," nccpr: ",nccpr
          write(myrank+NBFU,*) "-X 2nd & later: itr: ",itr, &
               &              " icxm0,icxm1: ",icxm0,icxm1," nccms: ",nccms
          write(myrank+NBFU,*) "-X 2nd & later: itr: ",itr, &
               &              " icxbm0,icxbm1: ",icxbm0,icxbm1," nccmr: ",nccmr
#endif
#ifdef MCOM_TIMER
          te = MPI_WTIME()
          tcomll_lcl = tcomll_lcl + te - ts
          call mpi_barrier(mpi_comm_world,ierr)
          ts = MPI_WTIME()
          tcomll_pbar = tcomll_pbar + ts - te
#endif
#ifdef PROFILE_COMM
          call mpi_barrier(mpi_comm_world, ierr)
          call start_profile(comm_fmm_local_multi_dr_id)
#endif
#ifdef FJ_RDMA
          if(nccpr > 0) call rdma_post(ipx_psrc)
          if(nccmr > 0) call rdma_post(ipx_msrc)
          if(nccps > 0) call rdma_wait(ipx_pdest)
          if(nccms > 0) call rdma_wait(ipx_mdest)
          !
          if(nccps > 0) then
             call rdma_put_post(ipx_pdest, &
                  wm_raddr(ipx_pdest+1) + (icbp0_x(il)%iaddr(itr)-1)*wm_4th_offset*ibds, &
                  wm_laddr              + (icxp0-1) *wm_4th_offset*ibds, &
                  nccps*ids)
          end if
          if(nccms > 0) then
             call rdma_put_post(ipx_mdest, &
                  wm_raddr(ipx_mdest+1) + (icbm0_x(il)%iaddr(itr)-1)*wm_4th_offset*ibds, &
                  wm_laddr              + (icxm0-1) *wm_4th_offset*ibds, &
                  nccms*ids)
          end if
          !
          if(nccpr > 0) call rdma_wait(ipx_psrc)
          if(nccmr > 0) call rdma_wait(ipx_msrc)
#else
          nrq = 0
          if(nccpr > 0) then
             nrq = nrq + 1
             call mpi_irecv(wm(1,1,1,1,icxbp0), nccpr, &
                  &                   MPI_PREC_MOMVEC_COM, ipx_psrc,     1, &
                  &                   mpi_comm_world, irq(nrq), ierr)

#ifdef DBG_COMLL
             write(myrank+NBFU,*) "+X 2nd & later: recipient itr: ",itr, &
                       &              " icxbp0,icxbp1: ",icxbp0,icxbp1," nccpr: ",nccpr
#endif

          endif                ! recipient (p)
          if(nccps > 0) then
             nrq = nrq + 1
             call mpi_isend(wm(1,1,1,1,icxp0), nccps, &
                  &                   MPI_PREC_MOMVEC_COM, ipx_pdest,    1, &
                  &                   mpi_comm_world, irq(nrq), ierr)

#ifdef DBG_COMLL
             write(myrank+NBFU,*) "+X 2nd & later: sender itr: ",itr, &
                  &              " icxp0,icxp1: ",icxp0,icxp1," nccps: ",nccps
#endif

          endif                ! sender (p)
          if(nccmr > 0) then
             nrq = nrq + 1
             call mpi_irecv(wm(1,1,1,1,icxbm0), nccmr, &
                  &                   MPI_PREC_MOMVEC_COM, ipx_msrc,     2, &
                  &                   mpi_comm_world, irq(nrq), ierr)

#ifdef DBG_COMLL
             write(myrank+NBFU,*) "-X 2nd & later: recipient itr: ",itr, &
                  &              " icxbm0,icxbm1: ",icxbm0,icxbm1," nccmr: ",nccmr
#endif

          endif                ! recipient (m)
          if(nccms > 0) then
             nrq = nrq + 1
             call mpi_isend(wm(1,1,1,1,icxm0), nccms, &
                  &                   MPI_PREC_MOMVEC_COM, ipx_mdest,    2, &
                  &                   mpi_comm_world, irq(nrq), ierr)

#ifdef DBG_COMLL
             write(myrank+NBFU,*) "-X 2nd & later: sender itr: ",itr, &
                  &              " icxm0,icxm1: ",icxm0,icxm1," nccms: ",nccms
#endif

          endif                ! sender (m)

          call mpi_waitall(nrq, irq, istatus, ierr)
#endif
#ifdef PROFILE_COMM
          call end_profile(comm_fmm_local_multi_dr_id)
#endif
#ifdef MCOM_TIMER
          te = MPI_WTIME()
          tcomll_com = tcomll_com + te - ts
          nccomll_com = nccomll_com + 1
#ifdef MCOM_TIMER_DETAIL
          write(myrank+NBFU,*) "Rank: ", myrank, &
               &                    " TcomLL-X-2nd&later: ", tle - tls
#endif
          ts = te
#endif

          !            endif                ! irregular&regular pairing.

       endif                ! iteration
    end do                  ! iteration

#ifdef FJ_RDMA
    call rdma_deallocate(4)
#endif
  end subroutine comm_fmm_local_multi_dr

!---------------------------------------------------------------------
  subroutine preprocess_comm(ip0,nbound_m,nbound_p, &
     &                       np,nl,nscell,nscdiv,  &
     &                       nb_destm,nb_destp,nitr)
!---------------------------------------------------------------------
  use dr_cntl, only : mtbd, mtbdp2
  implicit none
  integer(4),intent(in) :: ip0,nbound_m,nbound_p
  integer(4),intent(in) :: np,nl,nscell,nscdiv
  integer(4),intent(out) :: nb_destm,nb_destp,nitr
  integer(4) npow,icg0,icg1
  integer(4) ip

    if(    ip0 .gt. np)then
      ip=ip0-np  ! periodic boundary +
    elseif(ip0 .lt.  0)then
      ip=ip0+np  ! periodic boundary -
    else
      ip=ip0
    endif

    npow = nl

    if(npow == 2) then  ! cell merge number for upper level is two.

       !  ...   number of super-cells always has factor 2 when number of merging cells for
       !  ...   upper level is 2. boundary size in the destination process is always own
       !  ...   plus-side boundary size for plus direction communication or minus-side
       !  ...   boundary size for minus direction communication.
       nb_destm = nbound_m
       nb_destp = nbound_p

    else                ! cell merge number for upper level is three.

       !  ...   number of super-cells is a multiple of power of 3 number.
       !  ...   this means super-cells consist of multiple set of merging 3 cells.
       !  ...   when the position within merging cells of the last cell on neighbor rank
       !  ...   is 1st, 2nd, or 3rd the number of boundary cells is 5, 4, or 3, respectively.
       !  ...   when the position within merging cells of the 1st cell on neighbor rank
       !  ...   is 1st, 2nd, or 3rd the number of boundary cells is 3, 4, or 5, respectively.

       icg0 = (nscell * ip) / np + 1
       icg1 = (nscell * (ip + 1) - 1) / np + 1
       nb_destm = mtbdp2 - mod(icg0 + 1, npow)
       nb_destp = mtbd   + mod(icg1    , npow)

    endif

    nitr = max( (nbound_m - 1) / nscdiv + 1, (nbound_p - 1) / nscdiv + 1 )
    nitr = max( (nb_destm - 1) / nscdiv + 1, nitr)
    nitr = max( (nb_destp - 1) / nscdiv + 1, nitr)

  end subroutine preprocess_comm
!----------------------------------------------------------------------
  subroutine get_icbp01(nbd_m,nbd_p,nitr, &
     &                  nscdiv,nb_dest,icbp0,icbp1)
    implicit none
    integer(4),intent(in)  :: nbd_m,nbd_p,nitr,nscdiv,nb_dest
    integer(4),intent(out) :: icbp0(nitr),icbp1(nitr)
    integer(4) itr,icp0,icp1
    integer(4) icbp0_tmp,icbp1_tmp
    integer(4) ip1

       do itr = 1, nitr
          if (itr == 1) then                ! first iteration
             icp0 = nbd_m + nscdiv - nb_dest + 1
             if(icp0 < nbd_m + 1) icp0 = nbd_m + 1
             icp1 = nbd_m + nscdiv
             icbp0_tmp = nbd_m - nscdiv + 1
             if(icbp0_tmp < 1) icbp0_tmp = 1
             icbp1_tmp = nbd_m
          else      ! iteration follows == irregular&regular pairing. =
             icp0 = icbp0_tmp
             if(nscdiv*itr > nb_dest) &
                  &               icp0 = nbd_m + nscdiv - nb_dest + 1
             icp1 = icbp1_tmp
             icbp0_tmp = nbd_m - nscdiv*itr + 1
             if(icbp0_tmp < 1) icbp0_tmp = 1
             icbp1_tmp = nbd_m - nscdiv*(itr - 1)
          endif         ! iteration
          icbp0(itr)=icbp0_tmp
          icbp1(itr)=icbp1_tmp
       enddo

  end subroutine get_icbp01
!----------------------------------------------------------------------
  subroutine get_icbm01(nbd_m,nbd_p,nitr, &
     &                  nscdiv,nb_dest,icbm0,icbm1)
    implicit none
    integer(4),intent(in)  :: nbd_m,nbd_p,nitr,nb_dest,nscdiv
    integer(4),intent(out) :: icbm0(nitr),icbm1(nitr)
    integer(4) itr,icm0,icm1
    integer(4) icbm0_tmp,icbm1_tmp

       do itr = 1, nitr
          if (itr == 1) then                ! first iteration
             icm0 = nbd_m + 1
             icm1 = nbd_m + nb_dest
             if(icm1 > nbd_m + nscdiv) &
                  &                 icm1 = nbd_m + nscdiv
             icbm0_tmp = nbd_m + nscdiv + 1
             icbm1_tmp = nbd_m + 2*nscdiv
             if(icbm1_tmp > nbd_m + nscdiv + nbd_p) &
                  &                 icbm1_tmp = nbd_m + nscdiv + nbd_p
          else      ! iteration follows == irregular&regular pairing. =
             icm0 = icbm0_tmp
             icm1 = icbm1_tmp
             if(nscdiv*itr > nb_dest) &
                  &               icm1 = nbd_m + nb_dest
             icbm0_tmp = nbd_m + nscdiv*itr + 1
             icbm1_tmp = nbd_m + nscdiv*(itr + 1)
             if(icbm1_tmp > nbd_m + nscdiv + nbd_p) &
                  &               icbm1_tmp = nbd_m + nscdiv + nbd_p
          endif         ! iteration
          icbm0(itr)=icbm0_tmp
          icbm1(itr)=icbm1_tmp
       enddo

  end subroutine get_icbm01
!----------------------------------------------------------------------
end module comm_fmm_mod
