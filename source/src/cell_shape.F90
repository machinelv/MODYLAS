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
!! \brief Module and subroutine which converts lengths and angles of 
!!        unit cell into the geometry matrix.
!<
!----------------------------------------------------------------------
!>
!! \brief Module which converts lengths and angles of unit cell into
!!        the geometry matrix.
!! \author Noriyuki Yoshii
!<
!----------------------------------------------------------------------
module cell_shape

contains

!----------------------------------------------------------------------
!>
!! \brief Module which converts lengths and angles of unit cell into
!!        cell geometry matrix box(:,:).
!! \author Noriyuki Yoshii
!<
  subroutine cell_convert1(wlx,wly,wlz,alpha,beta,gamma,box)
    use md_const  !! module+, wkx,wly,wlz
    implicit none
    real(8) :: wlx,wly,wlz,alpha,beta,gamma
    real(8) :: box(3,3)
    real(8) :: frac

    frac=dcos(alpha*deg2rad)-dcos(beta*deg2rad) &  !!! dcos,dsin check
         &    *dcos(gamma*deg2rad)

    box(1,1)=wlx  !! lengtha
    box(1,2)=wly*dcos(gamma*deg2rad)
    box(1,3)=wlz*dcos(beta*deg2rad)
    box(2,1)=0.d0
    box(2,2)=wly*dsin(gamma*deg2rad)
    box(2,3)=wlz*frac/dsin(gamma*deg2rad)
    box(3,1)=0.d0
    box(3,2)=0.d0
    box(3,3)=dsqrt(wlz**2-(box(1,3)**2+box(2,3)**2))
  end subroutine cell_convert1
!----------------------------------------------------------------------
!>
!! \brief Subroutine which converts the matrix representing shape of 
!!        simulation box into the lengths and angles of unit cell
!! \author Noriyuki Yoshii
!<
  subroutine cell_convert2(wlx,wly,wlz,alpha,beta,gamma,box)
!----------------------------------------------------------------------
    use md_const
    implicit none
    real(8) :: wlx,wly,wlz,alpha,beta,gamma
    real(8) :: box(3,3)
    real(8) :: ab,bc,ca

    wlx=box(1,1)
    wly=dsqrt(box(1,2)**2+box(2,2)**2)
    wlz=dsqrt(box(1,3)**2+box(2,3)**2+box(3,3)**2)

    ab=box(1,1)*box(1,2)
    bc=box(1,2)*box(1,3)+box(2,2)*box(2,3)
    ca=box(1,3)*box(1,1)

    alpha = dacos(bc/wly/wlz)*rad2deg
    beta  = dacos(ca/wlz/wlx)*rad2deg
    gamma = dacos(ab/wlx/wly)*rad2deg
  end subroutine cell_convert2
!----------------------------------------------------------------------
!>
!! \brief Subroutine which prepares the parameters of simulation box
!! \author Noriyuki Yoshii
!<
  subroutine cell_convert3(wlx,wly,wlz,alpha,beta,gamma,wlxh,wlyh,wlzh,cellvol, &
       &  sinbeta,cosbeta,singamma,cosgamma,b_factor,g_factor)
!----------------------------------------------------------------------
    use md_const
    implicit none
    real(8) :: wlx,wly,wlz,alpha,beta,gamma
    real(8) :: wlxh,wlyh,wlzh
    real(8) :: cosalpha
    real(8) :: sinbeta,cosbeta
    real(8) :: singamma,cosgamma
    real(8) :: b_factor,g_factor
    real(8) :: cellvol

    wlxh = wlx/2d0
    wlyh = wly/2d0
    wlzh = wlz/2d0

    cosalpha = cos(alpha*deg2rad)
    sinbeta  = sin(beta *deg2rad)
    cosbeta  = cos(beta *deg2rad)
    singamma = sin(gamma*deg2rad)
    cosgamma = cos(gamma*deg2rad)
    b_factor = (cosalpha-cosbeta*cosgamma)/singamma
    g_factor = sqrt(1.0d0-cosbeta**2-b_factor**2)

    cellvol  = (singamma*g_factor) * wlx * wly * wlz
  end subroutine cell_convert3
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to initilize periodic boundary condition
!!         for constant pressure.
!! \author Kensuke Iwahashi
!<
  subroutine init_md_periodic
    use unit_cell
    use md_periodic, only : md_periodic__type, CLUSTER
    implicit none
    if (md_periodic__type == CLUSTER)  return

    call cell_convert3(cellx,celly,cellz,alpha,beta,gamma,cellxh,cellyh, &
         &  cellzh,cellvol,sinbeta,cosbeta, &
         &  singamma,cosgamma,b_factor,g_factor)
  end subroutine init_md_periodic
!----------------------------------------------------------------------
!>
!! \brief Subroutine to apply scaling of unit cell volume.
!! \author Yoshimichi Andoh
!<
  subroutine apply_volume_scaling
    use unit_cell
    use extended_system
    implicit none
    real(8) :: volscale_factor=1d0
    volscale_factor=sysvol0/cellvol
    volscale_factor=volscale_factor**(1d0/3d0)
    cellx=cellx*volscale_factor
    celly=celly*volscale_factor
    cellz=cellz*volscale_factor
    call cell_convert3(cellx,celly,cellz,alpha,beta,gamma,cellxh,cellyh, &
         &  cellzh,cellvol,sinbeta,cosbeta, &
         &  singamma,cosgamma,b_factor,g_factor)
  end subroutine apply_volume_scaling

end module cell_shape

