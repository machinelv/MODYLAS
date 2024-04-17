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
!! \brief  Module and subroutines to distribute atom to subcells
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to distribute atom to subcells
!! \author Yoshimichi Andoh
!<
module dist_atom2cell
implicit none

contains 

!----------------------------------------------------------------------
!>
!! \brief  Wrapper subroutine to atom2cell
!! \author Kazushi Fujimoto  (Yoshimichi Andoh mod.)
!<
  subroutine atom2cell_wrap
!----------------------------------------------------------------------
    use center_of_mass
    implicit none

    call atom2cell     ! default
#ifdef SEGSHAKE
    call list_own_i0
#endif

  end subroutine atom2cell_wrap

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to initialize meta-deta
!! \author Yoshimichi Andoh
!<
  subroutine atom2cell
!----------------------------------------------------------------------
    use subcell, only : tag,na_per_cell,m2i,i2m,lsegtop,lseg_natoms, &
   &                    ndseg_fmmn
    use domain, only : lxdiv, lydiv, lzdiv
    use trajectory_org
    use trajectory_mpi
    use dr_cntl, only : nbd,nbd2
    implicit none
    integer(4) :: i, ii,ic, i0, i00, iz,iy,ix
    integer(4) :: j0x,j0y,j0z
    integer(4) :: icx0,icy0,icz0
    integer(4) :: jx,jy,jz
    include 'mpif.h'
    integer(4) :: k0
! closed in this subrotine
    integer(4),allocatable :: m2i2(:)

    allocate(m2i2(nadirect))
    !### initialize ###
!$omp parallel default(shared) &
!$omp& private(ix,iy,iz,i,i0)
!$omp do
!default: MTD
    do ix=0,lxdiv+nbd 
       do iy=0,lydiv+nbd 
          do iz=0,lzdiv+nbd 
             tag(iz,iy,ix)=0
             na_per_cell(iz,iy,ix)=0
          enddo
       enddo
    enddo
!$omp end do nowait
!$omp do
    do i=1,n
       i2m(i) =-1       ! initialize
    enddo
!$omp end do
!$omp do
    do i0=1,nadirect
       m2i2(i0)=m2i(i0) ! store m2i
       m2i(i0)=-1       ! initialize
    enddo
!$omp end do
!$omp end parallel

    !### self range ###
    k0=0   !! a count for segment number
!default: MTD
    do j0x=0,lxdiv+nbd 
       do j0y=0,lydiv+nbd 
          do j0z=0,lzdiv+nbd 
             IF(j0x.ge.1.and.j0x.le.lxdiv .and. &
          &     j0y.ge.1.and.j0y.le.lydiv .and. &
          &     j0z.ge.1.and.j0z.le.lzdiv )THEN
                if(j0z==nbd)THEN
                   i00=(j0x)*na1cell*(lydiv+nbd2)*(lzdiv+nbd2) &
              &       +(j0y)*na1cell*(lzdiv+nbd2) &
              &       +(j0z)*na1cell
                endif
             ELSE
                cycle ! outside of myrank
             ENDIF

!default: MTD
             icx0 = j0x-nbd+1 ! address in myrank 
             icy0 = j0y-nbd+1 !
             icz0 = j0z-nbd+1 !
             tag(j0z,j0y,j0x)=i00+1
             do ic = 1, ndseg_fmmn(icz0,icy0,icx0)
                k0=k0+1 ! local
                do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
                   i00=i00+1 ! local
                   m2i(i00)=m2i2(i0)
                enddo ! i0
                lsegtop(k0)=i00-lseg_natoms(k0)+1 ! renew lsegtop
             enddo ! ic

             IF(ndseg_fmmn(icz0,icy0,icx0)==0)then
                na_per_cell(j0z,j0y,j0x)=0
             ELSE
                na_per_cell(j0z,j0y,j0x)=i00-tag(j0z,j0y,j0x)+1
                if(na_per_cell(j0z,j0y,j0x).gt.na1cell)then
!! note this error will be relaxed, only when i00 > na5cell
                   write(0,*) 'ERROR: na_per_cell over flowed!', &
                        &      na_per_cell(j0z,j0y,j0x),na1cell
                   call modylas_abort()
                endif
             ENDIF

          enddo ! jz
       enddo ! jy
    enddo ! jx

!$omp parallel default(shared) &
!$omp& private(ii,i0,i,icx0,icy0,icz0)
!$omp do
    do ii=1,lxdiv*lydiv*lzdiv
!default: MTD
       icz0=mod(ii-1,lzdiv)     +nbd
       icy0=mod(ii-1,lzdiv*lydiv)
       icy0=icy0/lzdiv          +nbd
       icx0=(ii-1)/(lzdiv*lydiv)+nbd
       do i0=tag(icz0,icy0,icx0), &
      &      tag(icz0,icy0,icx0)+na_per_cell(icz0,icy0,icx0)-1
          i=m2i(i0)
if(i.le.0)then
   write(0,*) 'ERROR: i is less than 0 in atom2cell.'
   call modylas_abort
endif
          i2m(i)=i0
          wkxyz(1,i0)=xyz(1,i)
          wkxyz(2,i0)=xyz(2,i)
          wkxyz(3,i0)=xyz(3,i)
          wkv(1:3,i0)=v(1:3,i)
       enddo ! i0
    enddo ! ii
!$omp end do
!$omp end parallel

    deallocate(m2i2)
end subroutine atom2cell

end module dist_atom2cell
