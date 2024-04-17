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
!! \brief Module and subroutine to apply periodic boundary condition.
!<
!----------------------------------------------------------------------
!>
!! \brief Module containing subroutines to apply PBC.
!! \author Yoshimichi Andoh
!<
module boundary
  use subcell
  implicit none

contains
!----------------------------------------------------------------------
!<
!! \brief Subroutines to apply PBC for segment mass center.
!! \author Yoshimichi Andoh
!<
  subroutine cell_edge
!----------------------------------------------------------------------
    use atom_mass
    use trajectory_org
    use unit_cell
    use cell_shape
    use segments
    use param
    use matrix_inverse_mod
    implicit none
    integer(4) :: i,isegstart,isegend,k,ipar
    real(8) :: totm,tmpx,tmpy,tmpz,segrx,segry,segrz
    real(8) :: segsx,segsy,segsz
    real(8) :: box(3,3),r_box(3,3)

    call cell_convert1(cellx,celly,cellz,alpha,beta,gamma,box)
    call matrix_inverse(box,r_box)

!$omp parallel do default(none) &
!$omp& private(k,i,isegstart,isegend,ipar) &
!$omp& private(tmpx,tmpy,tmpz,totm) &
!$omp& private(segsx,segsy,segsz,segrx,segry,segrz) &
!$omp& shared(nsegments,segtop,seg_natoms,box,r_box) &
!$omp& shared(paranum,mass,xyz,seg_cx,seg_cy,seg_cz)
    do k = 1, nsegments
       isegstart=segtop(k)+1
       isegend  =isegstart+seg_natoms(k)-1
       tmpx=0.0d0;tmpy=0.0d0;tmpz=0.0d0;totm=0d0
       do i = isegstart, isegend
          ipar=paranum(i)
#ifdef TIP4
          if(mass(ipar) < 0) cycle
#endif
          totm=totm+mass(ipar)
          tmpx=tmpx+xyz(1,i)*mass(ipar)
          tmpy=tmpy+xyz(2,i)*mass(ipar)
          tmpz=tmpz+xyz(3,i)*mass(ipar)
       enddo ! i
       tmpx=tmpx/totm
       tmpy=tmpy/totm
       tmpz=tmpz/totm
       !       ### normalize
       segsx=r_box(1,1)*tmpx+r_box(1,2)*tmpy+r_box(1,3)*tmpz
       segsy=r_box(2,1)*tmpx+r_box(2,2)*tmpy+r_box(2,3)*tmpz
       segsz=r_box(3,1)*tmpx+r_box(3,2)*tmpy+r_box(3,3)*tmpz
       !       ### check periodic boundary condition
       segsx=segsx-floor(segsx+0.5d0)  ! -1,0,+1
       segsy=segsy-floor(segsy+0.5d0)  ! -1,0,+1
       segsz=segsz-floor(segsz+0.5d0)  ! -1,0,+1
       !       ### recover unit
       segrx=box(1,1)*segsx+box(1,2)*segsy+box(1,3)*segsz
       segry=box(2,1)*segsx+box(2,2)*segsy+box(2,3)*segsz
       segrz=box(3,1)*segsx+box(3,2)*segsy+box(3,3)*segsz
       !       ### segment
       seg_cx(k)=segrx
       seg_cy(k)=segry
       seg_cz(k)=segrz
       !       ### atom
       do i = isegstart, isegend
          xyz(1,i)=xyz(1,i)+segrx-tmpx
          xyz(2,i)=xyz(2,i)+segry-tmpy
          xyz(3,i)=xyz(3,i)+segrz-tmpz
       enddo ! i
    enddo ! k
  end subroutine cell_edge
!----------------------------------------------------------------------
!>
!! \brief Subroutines to apply PBC for i-j atom distance. 
!!        (-DONEPROC_AXIS pressumed.)
!! \author Yoshimichi Andoh
!<
  subroutine pbc_pair(dlx,dly,dlz)
!----------------------------------------------------------------------
    use unit_cell
    implicit none
    real(8), intent(inout) :: dlx,dly,dlz

    if(dlx.ge.+cellxh) dlx=dlx-cellx;if(dlx.lt.-cellxh) dlx=dlx+cellx
    if(dly.ge.+cellyh) dly=dly-celly;if(dly.lt.-cellyh) dly=dly+celly
    if(dlz.ge.+cellzh) dlz=dlz-cellz;if(dlz.lt.-cellzh) dlz=dlz+cellz
  end subroutine pbc_pair
!----------------------------------------------------------------------
  subroutine pbc_pair_dr(dlx,dly,dlz)
!>
!! \brief Subroutines to apply PBC for i-j atom distance. 
!!        (-DONEPROC_AXIS pressumed.)
!! \author Yoshimichi Andoh
!<
!----------------------------------------------------------------------
    use unit_cell
    use precision  !precision module
    implicit none
    real(kind=wrp), intent(inout) :: dlx,dly,dlz

    if(dlx.ge.+cellxh) dlx=dlx-cellx;if(dlx.lt.-cellxh) dlx=dlx+cellx
    if(dly.ge.+cellyh) dly=dly-celly;if(dly.lt.-cellyh) dly=dly+celly
    if(dlz.ge.+cellzh) dlz=dlz-cellz;if(dlz.lt.-cellzh) dlz=dlz+cellz
  end subroutine pbc_pair_dr
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to apply the periodic boundary condition
!! \author Yoshimichi Andoh
!<
  subroutine apply_pbc()
!-----------------------------------------------------------------------
    use trajectory_mpi
    use param
    use mpi_tool
    use unit_cell
    use cell_shape
    use mpi_3d_grid, only : npx, npy, npz
    use domain, only : lxdiv, lydiv, lzdiv, ixmin, iymin, izmin, ixmax, iymax, &
         & izmax, ncellx, ncelly, ncellz
    use dr_cntl, only : nbd,nbd2
    implicit none
    integer(4)::j0
    real(8) :: box(3,3),av(3),bv(3),cv(3)
    integer(4)::jxb,jyb,iz,jzb

    call cell_convert1(cellx,celly,cellz,alpha,beta,gamma,box)
    av(:)=box(:,1)
    bv(:)=box(:,2)
    cv(:)=box(:,3)

!default: MTD
    !------------ -Z direction
    if(    izmin==1)then

!$omp parallel private(j0)
      do jxb=0,lxdiv+nbd
        do jyb=0,lydiv+nbd
!$omp do
          do j0=tag(0,jyb,jxb), &
     &          tag(0,jyb,jxb) + na_per_cell(0,jyb,jxb)-1
            wkxyz(:,j0)=wkxyz(:,j0)-cv(:)
          enddo
!$omp end do nowait
        enddo
      enddo
!$omp end parallel

    endif
    !------------ +Z direction
    if(izmax==ncellz)then

!$omp parallel private(j0)
      do jxb=0,lxdiv+nbd
        do jyb=0,lydiv+nbd
!$omp do
          do j0=tag(lzdiv+1,jyb,jxb), &
     &          tag(lzdiv+1,jyb,jxb) + na_per_cell(lzdiv+1,jyb,jxb)-1
            wkxyz(:,j0)=wkxyz(:,j0)+cv(:)
          enddo
!$omp end do nowait
        enddo
      enddo
!$omp end parallel

    endif

    !------------ -Y direction
    if(    iymin==1)then

!$omp parallel private(j0)
      do jxb=0,lxdiv+nbd
        do jyb=0,0
!$omp do
          do j0=tag(0,jyb,jxb), &
     &          tag(lzdiv+nbd,jyb,jxb) + na_per_cell(lzdiv+nbd,jyb,jxb)-1
            wkxyz(:,j0)=wkxyz(:,j0)-bv(:)
          enddo
!$omp end do nowait
        enddo
      enddo
!$omp end parallel

    endif
    !------------ +Y direction
    if(iymax==ncelly)then

!$omp parallel private(j0)
      do jxb=0,lxdiv+nbd
        do jyb=lydiv+nbd,lydiv+nbd
!$omp do
          do j0=tag(0,jyb,jxb), &
     &          tag(lzdiv+nbd,jyb,jxb) + na_per_cell(lzdiv+nbd,jyb,jxb)-1
            wkxyz(:,j0)=wkxyz(:,j0)+bv(:)
          enddo
!$omp end do nowait
        enddo
      enddo
!$omp end parallel

    endif

    !------------ -X direction
    if(    ixmin==1)then

!$omp parallel private(j0)
      do jxb=0,0
        do jyb=0,lydiv+nbd
!$omp do
          do j0=tag(0,jyb,jxb), &
     &          tag(lzdiv+nbd,jyb,jxb) + na_per_cell(lzdiv+nbd,jyb,jxb)-1
            wkxyz(:,j0)=wkxyz(:,j0)-av(:)
          enddo
!$omp end do nowait
        enddo
      enddo
!$omp end parallel

    endif
    !------------ +X direction
    if(ixmax==ncellx)then

!$omp parallel private(j0)
      do jxb=lxdiv+nbd,lxdiv+nbd
        do jyb=0,lydiv+nbd
!$omp do
          do j0=tag(0,jyb,jxb), &
     &          tag(lzdiv+nbd,jyb,jxb) + na_per_cell(lzdiv+nbd,jyb,jxb)-1
            wkxyz(:,j0)=wkxyz(:,j0)+av(:)
          enddo
!$omp end do nowait
        enddo
      enddo
!$omp end parallel

    endif

  end subroutine apply_pbc

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to apply PBC to dx,dy,dz
!! \author Yoshimichi Andoh
!<
  subroutine PBCij(dx,dy,dz)   !! not used with npt_pr
    use unit_cell
    implicit none
    real(8), intent(inout) :: dx,dy,dz

    if(cuboid)then

       if(dx.gt.+cellx/2d0) dx=dx-cellx
       if(dx.lt.-cellx/2d0) dx=dx+cellx
       if(dy.gt.+celly/2d0) dy=dy-celly
       if(dy.lt.-celly/2d0) dy=dy+celly
       if(dz.gt.+cellz/2d0) dz=dz-cellz
       if(dz.lt.-cellz/2d0) dz=dz+cellz

    else
       write(0,*) 'ERROR: cell shape must be cuboid'
       call modylas_abort()
    endif
  end subroutine PBCij
!----------------------------------------------------------------------
!>
!! \brief Subroutines to scale cell geometry and apply PBC.
!! \author Yoshimichi Andoh
!<
  subroutine apply_cellgeometry_scaling
    use unit_cell
    use extended_system
    use trajectory_mpi
    use matrix_inverse_mod
    use cell_shape
    implicit none
    integer(4) :: i0,k0
    real(8) :: tmpx,tmpy,tmpz,segrx,segry,segrz
    real(8) :: segsx,segsy,segsz
    real(8) :: box0(3,3), r_box0(3,3)
    real(8) :: box(3,3),r_box(3,3)

!target
    call cell_convert1(cellx0,celly0,cellz0,alpha0,beta0,gamma0,box0)
    call matrix_inverse(box0,r_box0)
!present
    call cell_convert1(cellx,celly,cellz,alpha,beta,gamma,box)
    call matrix_inverse(box,r_box)

    do k0=1,nselfseg
       tmpx=wseg_cx(k0)
       tmpy=wseg_cy(k0)
       tmpz=wseg_cz(k0)
       !       ### normalize
       segsx=r_box(1,1)*tmpx+r_box(1,2)*tmpy+r_box(1,3)*tmpz
       segsy=r_box(2,1)*tmpx+r_box(2,2)*tmpy+r_box(2,3)*tmpz
       segsz=r_box(3,1)*tmpx+r_box(3,2)*tmpy+r_box(3,3)*tmpz
       !       ### recover unit
       segrx=box0(1,1)*segsx+box0(1,2)*segsy+box0(1,3)*segsz
       segry=box0(2,1)*segsx+box0(2,2)*segsy+box0(2,3)*segsz
       segrz=box0(3,1)*segsx+box0(3,2)*segsy+box0(3,3)*segsz
       !       ### segment
       wseg_cx(k0)=segrx
       wseg_cy(k0)=segry
       wseg_cz(k0)=segrz
       !       ### atom
       do i0 = lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          wkxyz(1,i0)=wkxyz(1,i0)+segrx-tmpx
          wkxyz(2,i0)=wkxyz(2,i0)+segry-tmpy
          wkxyz(3,i0)=wkxyz(3,i0)+segrz-tmpz
       enddo ! i
    enddo ! k

    cellx=cellx0; celly=celly0; cellz=cellz0
    alpha=alpha0; beta =beta0 ; gamma=gamma0
    call cell_convert3(cellx0,celly0,cellz0,alpha0,beta0,gamma0, &
         &  cellxh,cellyh, cellzh, &
         &  cellvol,sinbeta,cosbeta, &
         &  singamma,cosgamma,b_factor,g_factor)
    call cell_convert1(cellx0,celly0,cellz0,alpha0,beta0,gamma0,box)
    call apply_pbc()

  end subroutine apply_cellgeometry_scaling
!----------------------------------------------------------------------

end module boundary

