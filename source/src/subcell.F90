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
!! \brief  Module and subroutines to set up and to maintain meta-deta
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to set up and to maintain meta-deta
!! \author Yoshimichi Andoh, Shin-ichi Ichikawa
!<
!----------------------------------------------------------------------
module subcell

  use mpi_tool
  use mpi_3d_grid, only : npx, npy, npz
  use domain, only : lxdiv, lydiv, lzdiv, ixmin, iymin, izmin, &
       & ncellx, ncelly, ncellz, max_nsegments_per_cell
  implicit none

  integer(4) :: nselfatm=0 !< number of atoms in myrank region
  integer(4) :: nselfseg=0 !< number of segments in myrank region
  integer(4) :: ndatm=0 !< number of atoms in myrank region + its halo region
  integer(4) :: max_seg !< maximum number of sements, used only in allocation

  integer(4),allocatable :: lsegtop(:)     !< store local segment number
  integer(4),allocatable :: lseg_natoms(:) !< store natom in segment

  integer(4),allocatable :: i2m(:) !< converter of atom id from global to local
  integer(4),allocatable :: m2i(:) !< converter of atom id from local to global
  integer(4),allocatable :: tag(:,:,:) !< store address of the first atom in subcell
  integer(4),allocatable :: na_per_cell(:,:,:) !< store atom mumber in subcell

  integer(4),allocatable :: ndseg_fmmn(:,:,:)  !< number of segment in subcell

  real(8),allocatable :: wseg_cx(:) !< mass center of segmemt x (local)
  real(8),allocatable :: wseg_cy(:) !< mass center of segment y (local)
  real(8),allocatable :: wseg_cz(:) !< mass center of segment z (local)

contains

!-----------------------------------------------------------------------
  subroutine alloc_i2m(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    integer(4) :: i
    allocate(i2m(ivalue))

!$omp parallel do default(shared)
    do i = 1,ivalue
       i2m(i) = -1
    end do
  end subroutine alloc_i2m

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to update segment mass centers
!! \author Yoshimichi Andoh
!<
  subroutine update_wsegc()
!----------------------------------------------------------------------
    use atom_mass
    use trajectory_mpi
    use param
    implicit none
    integer(4) :: i,ipar,i0,k0
    real(8) :: totm
    real(8) :: tmpx,tmpy,tmpz

!$omp parallel do default(none) &
!$omp& private(i0,k0) &
!$omp& private(i,ipar,tmpx,tmpy,tmpz,totm) &
!$omp& shared(nselfseg,mass) &
!$omp& shared(paranum) &
!$omp& shared(wseg_cx,wseg_cy,wseg_cz,wkxyz) &
!$omp& shared(m2i,lsegtop,lseg_natoms)
    do k0=1,nselfseg
       tmpx=0.0d0
       tmpy=0.0d0
       tmpz=0.0d0
       totm=0.0d0
       do i0 = lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          i=m2i(i0)
          ipar=paranum(i)
#ifdef TIP4
          if(mass(ipar) < 0) cycle
#endif
          totm=totm+mass(ipar)
          tmpx=tmpx+wkxyz(1,i0)*mass(ipar)
          tmpy=tmpy+wkxyz(2,i0)*mass(ipar)
          tmpz=tmpz+wkxyz(3,i0)*mass(ipar)
       enddo ! i0
       tmpx=tmpx/totm
       tmpy=tmpy/totm
       tmpz=tmpz/totm
       wseg_cx(k0)=tmpx
       wseg_cy(k0)=tmpy
       wseg_cz(k0)=tmpz
    enddo ! k0
  end subroutine update_wsegc

!---------------------------------------------------------------------
!>
!! \brief  Subroutine to assign segments to subcell (global -> local)
!! \author Yoshimichi Andoh
!<
  subroutine calc_ia2c()  !! called once in main.F90
!---------------------------------------------------------------------
    use mpi_3d_grid, only : get_ip
    use domain, only : is_my_cell, is_cell_for_me
    use trajectory_mpi, only : nadirect
    use segments
    use unit_cell
    use cell_shape
    implicit none
    real(8) :: x0, y0, z0
    integer(4) :: icx, icy, icz
    integer(4) :: icx0, icy0, icz0, icxyz0
    integer(4) :: i, k, ll
    integer(4) :: k0,i0
    real(8) :: rdcellx,rdcelly,rdcellz
    integer(4) :: k00,k000,ii,ic
    real(8) :: xfrac,yfrac,zfrac
    real(8) :: box(3,3),av(3),bv(3),cv(3)
! closed in this subrotine
    integer(4),allocatable :: jtok(:)
    integer(4),allocatable :: ktoj(:),ktoj2(:)

    allocate(jtok(nsegments))
    allocate(ktoj(max_seg))
    allocate(ktoj2(max_seg))

    rdcellx=dble(ncellx)/cellx
    rdcelly=dble(ncelly)/celly
    rdcellz=dble(ncellz)/cellz
    call cell_convert1(cellx,celly,cellz,alpha,beta,gamma,box)
    av(:)=box(:,1)
    bv(:)=box(:,2)
    cv(:)=box(:,3)
    call cell_convert3(cellx,celly,cellz,alpha,beta,gamma,cellxh,cellyh, &
         &  cellzh,cellvol,sinbeta,cosbeta, &
         &  singamma,cosgamma,b_factor,g_factor)

    ndseg_fmmn = 0  !! allocated as ndseg_fmmn(lzdiv,lydiv,lxdiv)
    ktoj =-1
    jtok =-1
    m2i = 0
!debug
!       write(myrank+10000,*) ixmin,iymin,izmin; call flush(myrank+10000)
!debug
    i0 = 0 ! local atom id
    ll = 0 ! local segment id
    do k = 1, nsegments   ! nsegment is total number of segment in system
       zfrac= seg_cz(k)/g_factor
       yfrac=(seg_cy(k)-zfrac*b_factor)/singamma
       xfrac= seg_cx(k)-yfrac*cosgamma-zfrac*cosbeta
       x0 = xfrac+cellxh
       y0 = yfrac+cellyh
       z0 = zfrac+cellzh
       icx=min0( int(x0*rdcellx),ncellx-1 )+1
       icy=min0( int(y0*rdcelly),ncelly-1 )+1
       icz=min0( int(z0*rdcellz),ncellz-1 )+1

       if(is_cell_for_me(icx,icy,icz)) then
          if(is_my_cell(icx, icy, icz)) then
             icx0=icx-ixmin+1 ! global -> local
             icy0=icy-iymin+1 ! global -> local
             icz0=icz-izmin+1 ! global -> local
             icxyz0 = (icx0-1)*(lydiv*lzdiv) + (icy0-1)*lzdiv + (icz0-1)
             ndseg_fmmn(icz0,icy0,icx0) = ndseg_fmmn(icz0,icy0,icx0) + 1
             if (ndseg_fmmn(icz0,icy0,icx0) .gt. max_nsegments_per_cell) then
                write(0,'(a)') 'ERROR: ndseg_fmmn overflow'
                write(0,*) ndseg_fmmn(icz0,icy0,icx0),max_nsegments_per_cell
                call modylas_abort()
             endif
             k0=max_nsegments_per_cell*icxyz0+ndseg_fmmn(icz0,icy0,icx0)
             wseg_cx(k0) = seg_cx(k)
             wseg_cy(k0) = seg_cy(k)
             wseg_cz(k0) = seg_cz(k)
             ktoj2(k0)=k
             ll=ll+1

             do i = segtop(k)+1,segtop(k)+seg_natoms(k)
                i0 = i0 + 1
                if(i0.gt.nadirect) then
                   write(0,*) 'ERROR: i0 overflow in calc_ia2c'
                   call modylas_abort()
                endif
                i2m(i)=i0   ! this i0 start with 1, not same as i0 in main
                m2i(i0)=i
             end do

          end if
       end if
    enddo
    nselfatm = i0
    nselfseg = ll

    !^^^ Left-align wseg ^^^!
    k000=0
    do ii=1,lxdiv*lydiv*lzdiv
       icz0=mod(ii-1,lzdiv)     +1
       icy0=mod(ii-1,lzdiv*lydiv)
       icy0=icy0/lzdiv          +1
       icx0=(ii-1)/(lzdiv*lydiv)+1
       icxyz0 = (icx0-1)*(lzdiv*lydiv) + (icy0-1)*lzdiv + (icz0-1)
       k00 = max_nsegments_per_cell*icxyz0
       do ic=1,ndseg_fmmn(icz0,icy0,icx0)
          k000=k000+1
          k0=k00+ic
          seg_cx(k000)=wseg_cx(k0)
          seg_cy(k000)=wseg_cy(k0)
          seg_cz(k000)=wseg_cz(k0)
          k=ktoj2(k0)
          ktoj(k000)=k
          jtok(k)=k000
       enddo ! ic
    enddo ! ii

    !^^^ reset wseg_cx ^^^!
    wseg_cx=0d0
    wseg_cy=0d0
    wseg_cz=0d0
!$omp parallel default(none) &
!$omp& private(k0,k,i,i0) &
!$omp& shared(nselfseg,wseg_cx,wseg_cy,wseg_cz,i2m,ktoj) &
!$omp& shared(seg_cx,seg_cy,seg_cz) &
!$omp& shared(segtop,seg_natoms,lsegtop,lseg_natoms)
!$omp do
    do k0=1,nselfseg
       wseg_cx(k0)=seg_cx(k0)
       wseg_cy(k0)=seg_cy(k0)
       wseg_cz(k0)=seg_cz(k0)
       k = ktoj(k0)
       i = segtop(k)+1
       i0 = i2m(i)
       lsegtop(k0)    =i0
       lseg_natoms(k0)=seg_natoms(k)
    enddo
!$omp end do
!$omp end parallel

    deallocate(ktoj,ktoj2)
    deallocate(jtok)
  end subroutine calc_ia2c

end module subcell
