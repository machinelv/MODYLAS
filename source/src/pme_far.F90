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
!! \brief Modules and subroutines which relates to the reciprocal part
!!        of the particle mesh Ewald (PME) method
!<
!----------------------------------------------------------------------
!>
!! \brief Module to store mathmatical functions for spline.
!! \author Kazushi Fujimoto and Shinichi Ichikawa
!<
!------------------------------------------------------------------------------
module PmeMdMath
  ! The mathmatical functions for spline.
  !    Bspline1: Spline function.
  !    Spline: Spline function by recursive cord.
  !    Dspline: Derivative of spline functions
  use precision  !precision module is loaded by default to declare PME subroutine
  real(kind=wrrp),parameter :: PI=3.1415926535897932384626433832795029d0
  real(kind=wrrp),parameter :: PI_sqrt=1.77245385090551602729816748334d0
  real(kind=wrrp),parameter :: r_PI_sqrt=1.0d0/PI_sqrt
contains
  subroutine Bspline1(x, n, work)
    implicit none
    integer(4)     :: n, i, j
    real(kind=dp)  :: x, coeff, work(n)

    work(1) = 1.0d0 - x
    work(2) = x

    do i = 3, n
        coeff = 1.0d0 / dble(i - 1)
        work(i) = x * work(i-1) * coeff
        do j = 1, i - 2
            work(i - j) = coeff * ((dble(j) + x) * work(i - j - 1) + (dble(i - j) - x) * work(i - j))
        enddo
        work(1) = (1.0d0 - x) * work(1) * coeff
    enddo
  end subroutine
  ! M. Griebel et al. Numerical Simulation in Molecular Dynamics Springer 
  recursive function Spline(x, p) result(spl)
    implicit none
    real(kind=dp)         :: x, spl
    integer(4), intent(in) :: p 
    !write(*,*) x, spl, p
    if((x <= 0.0d0) .or. (x >= dble(p))) then
        spl =  0.0d0
    else
      if(p == 1) then
        spl = 0.0d0
      elseif(p == 2) then
        spl = 1.0d0 - dabs(x - 1.0d0);
      else
        spl =  (x * spline(x, p - 1) + (dble(p) - x) * spline(x - 1.0d0, p - 1)) / (dble(p) - 1.0d0);
      endif
    endif
  end function
  real(8) function Dspline(x, p)
    real(kind=dp)  :: x
    integer(4)     :: p
    Dspline = spline(x, p - 1) - spline(x - 1.0d0, p - 1)
  end function
end module PmeMdMath
!------------------------------------------------------------------------------
!>
!! \brief Module to store variables for PME method.
!! \author Kazushi Fujimoto and Shinichi Ichikawa
!<
module PmeVariable
  use precision  !precision module is loaded by default to declare PME subroutine
  implicit none
  integer(4) :: maxorder, bsorder
  integer(4) :: nfft(3), nffthalf(3), ngrid
  real(kind=dp) :: ewald2Inverse, ewald_alpha
  real(kind=dp) :: qfact
  real(kind=wrrp), parameter:: kCHARGECONSTANT = 33.2063709667747d0
  real(kind=dp ),  allocatable:: bk2(:, :)
  real(kind=wrrp), allocatable:: coefficient(:,:,:)
  real(kind=wrrp), allocatable:: msv1(:), msv2(:), msv3(:)
end module PmeVariable
!------------------------------------------------------------------------------
!>
!! \brief Module to store cell parameter.
!! \author Kazushi Fujimoto and Shinichi Ichikawa
!<
module MdCellParameter
  use precision  !precision module is loaded by default to declare PME subroutine
  implicit none
  real(kind=dp):: volume, reciprocal_lattice_vector(6)
contains
!------------------------------------------------------------------------------
  subroutine GetCellParameter()
    ! Calculate reciprocal lattice vectors
    use cell_shape,  only : cell_convert3
    use unit_cell,   only : cellx, celly, cellz
    use unit_cell,   only : cellxh, cellyh, cellzh, cellvol
    use unit_cell,   only : sinbeta, cosbeta, singamma, cosgamma
    use unit_cell,   only : b_factor, g_factor
    use unit_cell,   only : alpha, beta, gamma
    use PmeVariable, only : nfft
    real(kind=wrrp) :: ar1,br1,br2,cr1,cr2,cr3
    real(kind=wrrp) :: rcellvol
    integer(4) :: i

    call cell_convert3(cellx,celly,cellz,alpha,beta,gamma,cellxh,cellyh, &
         &  cellzh,cellvol,sinbeta,cosbeta, &
         &  singamma,cosgamma,b_factor,g_factor)
    ar1 = cellx
    br1 = celly * cosgamma
    br2 = celly * singamma
    cr1 = cellz * cosbeta
    cr2 = cellz * b_factor
    cr3 = cellz * g_factor
    volume = cellvol
    rcellvol = 1.0d0/cellvol
    ! Calculate reciprocal lattice vectors
    reciprocal_lattice_vector(1) = (br2*cr3          ) * rcellvol ! 11
    reciprocal_lattice_vector(2) = (        - cr3*br1) * rcellvol ! 21
    reciprocal_lattice_vector(3) = (br1*cr2 - cr1*br2) * rcellvol ! 31
    reciprocal_lattice_vector(4) = (cr3*ar1          ) * rcellvol ! 22
    reciprocal_lattice_vector(5) = (        - ar1*cr2) * rcellvol ! 32
    reciprocal_lattice_vector(6) = (ar1*br2          ) * rcellvol ! 33
  end subroutine GetCellParameter
!------------------------------------------------------------------------------
end module MdCellParameter
!------------------------------------------------------------------------------
!>
!! \brief Module to reciprocal functions in PME method.
!! \author Kazushi Fujimoto and Shinichi Ichikawa
!<
module PmeReciprocalFunctions
  use precision  !precision module is loaded by default to declare PME subroutine
  implicit none
contains
!------------------------------------------------------------------------------
  subroutine CalculateBk2(bsord, maxo)
    ! bk(m_i) = 1.0d0 / 
    !           sum_{k=0}^{n-2} M_n(k+1) * exp(2 * pi * i * m_i * k / K_i) 
    ! n: bsorder
    ! K_i: The number of grids of i-axis
    use precision  !precision module is loaded by default to declare PME subroutine
    use PmeMdMath,   only : PI, Bspline1
    use PmeVariable, only : nfft, nffthalf
    use PmeVariable, only : bsorder, maxorder, bk2
    implicit none
    integer(4) :: bsord, maxo
    integer(4) :: nfftmax
    integer(4) ::  k, ix, m_i
    real(kind=dp) :: work(maxo), cos_v, sin_v, angle
  
    bsorder = bsord
    maxorder = maxo
    do ix = 1, 3
      nffthalf(ix) = nfft(ix) / 2
    enddo
    nfftmax = max(nfft(1), nfft(2), nfft(3))
    allocate(bk2(3, nfftmax))
    bk2(:, :) = 0.0d0
  
    ! IF n (bsorder) is ODD.
    if( mod(bsorder, 2) == 1 ) then
        do ix = 1, 3
            do m_i = 1, nfft(ix)
                ! IF n is odd and 2|m_i| = K_i this interpolation result fails but, 
                ! since it occurs in the tail of the reciprocal sum, we can set b(m_i) 
                ! arbitrarily to zero in this case.(P.8581 in JCP 103,8577(1995))
                if(m_i == int( 0.5d0 * dble(nfft(ix)) )) then
                    bk2(ix, m_i) = 0.0d0
                endif
                cos_v = 0.0d0
                sin_v = 0.0d0
                call Bspline1(0.0d0, bsorder, work)
                do k = 1, bsorder -1
                    angle = 2.0d0 * PI * dble((m_i - 1) * (k - 1)) / dble(nfft(ix))
                    cos_v = cos_v + work(k) * cos(angle)
                    sin_v = sin_v + work(k) * sin(angle)
                enddo
                bk2(ix, m_i) = 1.0d0 / (cos_v**2 + sin_v**2)
            enddo
        enddo
    ! IF n is EVEN.
    else
        do ix = 1, 3
            do m_i = 1, nfft(ix)
                cos_v = 0.0d0
                sin_v = 0.0d0
                call Bspline1(0.0d0, bsorder, work)
                do k = 1, bsorder - 1
                    angle = 2.0d0 * PI * dble((m_i - 1) * (k - 1)) / dble(nfft(ix))
                    cos_v = cos_v + work(k) * cos(angle)
                    sin_v = sin_V + work(k) * sin(angle)
                enddo
                bk2(ix, m_i) = 1.0d0 / (cos_v**2 + sin_v**2)
            enddo
        enddo
    endif

  end subroutine CalculateBk2
!------------------------------------------------------------------------------
  subroutine CalculateMsv()
    ! msv = m_i for 0 <= m_i <= K_i / 2
    !     = m_i - K_i for otherwise
    ! K_i: the number of grids of i-axis
    use PmeVariable, only : nfft, nffthalf, ewald2Inverse
    use PmeVariable, only : msv1, msv2, msv3
#ifndef PME_XYDIV
    use comm_pme_mod, only : ngxdivy, ngydivz, iminpy, jminpz
#else
    use comm_pme_mod, only : ngzdivy, ngydivx, kminpy, jminpx
#endif
    implicit none
    real(kind=dp) :: m1, m2, m3
    integer(4) :: k1, k2, k3
    integer(4) :: kd1, kd2, kd3

#ifndef PME_XYDIV
    do kd3 = 1, nfft(3)
       k3 = kd3 - 1
       if(k3 < nffthalf(3)) then
          m3 = dble(k3)
       else
          m3 = dble(k3 - nfft(3)) 
       endif
       msv3(kd3) = m3
    end do
    do kd1 = 1, ngxdivy
       k1 = iminpy + kd1 - 1
       if(k1 < nffthalf(1)) then
          m1 = dble(k1)
       else
          m1 = dble(k1 - nfft(1))
       endif
       msv1(kd1) = m1
    end do
    do kd2 = 1, ngydivz
       k2 = jminpz + kd2 - 1
       if(k2 < nffthalf(2)) then
          m2 = dble(k2)
       else
          m2 = dble(k2 - nfft(2))
       endif
       msv2(kd2) = m2
    end do
#else
    do kd1 = 1, nfft(1)
       k1 = kd1 - 1
       if(k1 < nffthalf(1)) then
          m1 = dble(k1)
       else
          m1 = dble(k1 - nfft(1))
       endif
       msv1(kd1) = m1
    end do
    do kd3 = 1, ngzdivy
       k3 = kminpy + kd3 - 1
       if(k3 < nffthalf(3)) then
          m3 = dble(k3)
       else
          m3 = dble(k3 - nfft(3))
       endif
       msv3(kd3) = m3
    end do
    do kd2 = 1, ngydivx
       k2 = jminpx + kd2 - 1
       if(k2 < nffthalf(2)) then
          m2 = dble(k2)
       else
          m2 = dble(k2 - nfft(2))
       endif
       msv2(kd2) = m2
    end do
#endif

  end subroutine CalculateMsv
!------------------------------------------------------------------------------
  subroutine CalculateCoefficient()
    ! C(m_1, m_2, m_3) * B(m_1, m_2, m_3) = exp(-pi * pi * m**2 / beta**2) / (m**2 * pi * V)
    !                                     * |b_1(m_1)|**2 * |b_2(m_2)|**2 * |b_3(m_3)|**2
    !     m = msv1 * a*_1 + msv2 * a*_2 + msv3 * a*_3
    !       a: reciprocal_lattice_vector
    !       msv = m_i for 0 <= m_i < K_i (K_i is # of girds of i-axis)
    !           = 0 for otherwise
    !       msv and b(m) were already calcualted in CalculateMsv and CalculateBk2, respectively.
    !     msq = |m|**2
    !     V: volume or the MD cell
    ! C(0, 0, 0) = 0
    ! Eq. 3.9 and 4.8 in JCP 103, 8577 (1995)
    use physics_const, only : md_QQ_4PiE
    use PmeMdMath,     only : PI
    use PmeVariable,   only : bk2, coefficient
    use PmeVariable,   only : nfft, nffthalf, ewald2Inverse
    use PmeVariable,   only : msv1, msv2, msv3
    use PmeVariable,   only : qfact
    use MdCellParameter, only : volume, reciprocal_lattice_vector
#ifndef PME_XYDIV
    use comm_pme_mod, only : ngxdivy, ngydivz, iminpy, jminpz
#else
    use comm_pme_mod, only : ngzdivy, ngydivx, kminpy, jminpx
#endif
    implicit none
    real(kind=dp) :: m1, m2, m3, m2b, m3b, m3t, msqn
    integer(4):: k1, k2, k3
    integer(4):: kd1, kd2, kd3

#ifndef PME_XYDIV

!$omp  parallel default(none) &
!$omp& private(k1,k2,k3,kd1,kd2,m1,m2,m3,m2b,m3b,m3t,msqn) &
!$omp& shared(ngydivz,ngxdivy,nfft,jminpz,iminpy) &
!$omp& shared(bk2,reciprocal_lattice_vector,volume,ewald2Inverse) &
!$omp& shared(msv1,msv2,msv3) &
!$omp& shared(coefficient)
!$omp  do collapse(2)
    do kd2 = 1, ngydivz
       do kd1 = 1, ngxdivy
          k2 = jminpz + kd2
          k1 = iminpy + kd1
          do k3 = 1, nfft(3)

             m1 = msv1(kd1) * reciprocal_lattice_vector(1)
             m2 = msv1(kd1) * reciprocal_lattice_vector(2) + msv2(kd2) * reciprocal_lattice_vector(4)
             m3 = msv1(kd1) * reciprocal_lattice_vector(3) + msv2(kd2) * reciprocal_lattice_vector(5) &
                 + msv3(k3) * reciprocal_lattice_vector(6)

             msqn = m1**2 + m2**2 + m3**2

             if(msqn > Cdzero) then
                 coefficient(k3, kd1, kd2) = &
                &                  dexp(-msqn * PI * PI * ewald2Inverse)  &
                &   / (msqn * PI * volume) * bk2(1, k1) * bk2(2, k2) * bk2(3, k3)
             !.. for wave number 0 with no operation in reciprocal space.
             else
                coefficient(k3, kd1, kd2) = Czero
             endif
          enddo
       enddo
    enddo
!$omp end do
!$omp end parallel

#else

!$omp  parallel default(none) &
!$omp& private(k1,k2,k3,kd1,kd2,m1,m2,m3,m2b,m3b,m3t,msqn) &
!$omp& shared(ngydivx,ngzdivy,nfft,jminpx,kminpy) &
!$omp& shared(bk2,reciprocal_lattice_vector,volume,ewald2Inverse) &
!$omp& shared(msv1,msv2,msv3) &
!$omp& shared(coefficient)
!$omp  do collapse(2)
    do kd2 = 1, ngydivx
       do kd3 = 1, ngzdivy
          k2 = jminpx + kd2
          k3 = kminpy + kd3
          do k1 = 1, nfft(1)

             m1 = msv1(k1) * reciprocal_lattice_vector(1)
             m2 = msv1(k1) * reciprocal_lattice_vector(2) + msv2(kd2) * reciprocal_lattice_vector(4)
             m3 = msv1(k1) * reciprocal_lattice_vector(3) + msv2(kd2) * reciprocal_lattice_vector(5) &
                + msv3(kd3) * reciprocal_lattice_vector(6)

             msqn = m1**2 + m2**2 + m3**2

             if(msqn > Cdzero) then
                 coefficient(k1, kd3, kd2) = &
                &                  dexp(-msqn * PI * PI * ewald2Inverse)  &
                &   / (msqn * PI * volume) * bk2(1, k1) * bk2(2, k2) * bk2(3, k3)
             !.. for wave number 0 with no operation in reciprocal space.
             else
                coefficient(k1, kd3, kd2) = Czero
             endif
          enddo
       enddo
    enddo
!$omp end do
!$omp end parallel

#endif

  end subroutine CalculateCoefficient
!------------------------------------------------------------------------------
  subroutine BsplineFill(thetai)
    ! Calcute the Cardinal B-splines M_2p(u_i) of i-th atom (thetai).
    ! u_i is the scaled fractional coordinate of i-the atom.
    !    u_i = K_1(or 2 or 3) * a*_1(or 2 or 3) * r_i
    !      a*: reciprocal lattice vector
    !      r_i: coordinate of i-the atom.
    !      0 <= u_i < K_1(or 2 or 3)

    use trajectory_mpi, only : wkxyz, nadirect
    use subcell,        only : tag, na_per_cell
    use domain,         only : lxdiv, lydiv, lzdiv
    use dr_cntl,        only : nbd
    use PmeVariable, only : bsorder
    use PmeVariable, only : nfft
    use MdCellParameter , only : reciprocal_lattice_vector
    implicit none
    real(kind=wrrp), intent(inout)   :: thetai(3, 2, bsorder, nadirect)
    integer(4) :: icx, icy, icz
    integer(4) :: i0, ix
    real(kind=dp) :: wr, rx(3), xi, yi, zi

!$omp  parallel default(none) &
!$omp& private(icx, icy, icz, i0, ix) &
!$omp& private(xi, yi, zi, rx, wr) &
!$omp& shared(reciprocal_lattice_vector, nfft) &
!$omp& shared(wkxyz, thetai, nadirect) &
!$omp& shared(lxdiv, lydiv, lzdiv, tag, na_per_cell)
!default: MTD
    do icx = 1-nbd, lxdiv+nbd
    do icy = 1-nbd, lydiv+nbd
!$omp  do
      do i0=tag(1-nbd,icy,icx), &
         &      tag(lzdiv+nbd,icy,icx)+na_per_cell(lzdiv+nbd,icy,icx)-1

        xi = wkxyz(1, i0)
        yi = wkxyz(2, i0)
        zi = wkxyz(3, i0)
        ! rx = a*_1(or 2 or 3) * r_i
        rx(1) = xi*reciprocal_lattice_vector(1) + yi*reciprocal_lattice_vector(2) + zi*reciprocal_lattice_vector(3) ! ri*kn
        rx(2) = yi*reciprocal_lattice_vector(4) + zi*reciprocal_lattice_vector(5)
        rx(3) = zi*reciprocal_lattice_vector(6)
        do ix = 1, 3
          ! wr = u_i - [u_i]
          wr = (rx(ix) - floor(rx(ix))) * dble(nfft(ix))
          wr = wr - floor(wr)      ! by SI.
          ! The only arguments of the fourth spline interpolation that do NOT have a value of zero are
          !    u_i - [u_i]       (0 <= u_i - [u_i]       <= 1)
          !    u_i - ([u_i] - 1) (1 <= u_i - ([u_i] - 1) <= 2)
          !    u_i - ([u_i] - 2) (2 <= u_i - ([u_i] - 2) <= 3)
          !    u_i - ([u_i] - 3) (3 <= u_i - ([u_i] - 3) <= 4)
          ! since M_2(u) = 0 for u < 0 or u > 2.
          call Bsplgen(wr, ix, i0, nadirect, thetai)
        enddo

      enddo  ! i0.
!$omp end do
    enddo
    enddo
!$omp end parallel

  end subroutine BsplineFill
!------------------------------------------------------------------------------
  subroutine Bsplgen(w,  dimension, atomID, natom, thetai)
    !theta(*, 1, *, *) : Spline function
    !theta(*, 2, *, *) : Differentiation of spline function
    use PmeMdMath,   only: Spline, Dspline
    use PmeVariable, only: bsorder
    implicit none
    integer(4),      intent(in) :: dimension, atomID, natom
    real(kind=dp),   intent(in) :: w
    real(kind=wrrp), intent(inout) :: thetai(3, 2, bsorder, natom)
    integer(4) :: i

    do i = 1, bsorder
        thetai(dimension, 1, bsorder + 1 - i, atomID) = Spline(dble(i-1) + w, bsorder)
        thetai(dimension, 2, bsorder + 1 - i, atomID) = Dspline(dble(i-1) + w, bsorder)
    enddo

  end subroutine Bsplgen
!------------------------------------------------------------------------------
  subroutine CalculateQgrid(natom, bxyz, bcharge, bthetai, wqgrid, nomp)
    ! Q(m_1, m_2, m_3) = sum_{i=1}^N sum_{n_1, n_2, n_3}  q_i * M_n(u_{1,i} -k_1 - n_1*K_1)
    !                                                         * M_n(u_{2,i} -k_2 - n_2*K_2)
    !                                                         * M_n(u_{3,i} -k_3 - n_3*K_3)
    ! M_n: n-th spline function. thetai in this program
    ! u_i: scaled fractional coordinate of i-th atom
    ! K_a: # of grids along a-th axis
    ! n_a: all integer

    ! note: Atom data bxyz, bcharge, and bthetai MUST INCLUDE that of cells, thickness 1 cell,
    !       outside of process boundary. This is to avoid communication of qgrid for superposing
    !       contribution from atoms in adjacent process.

    use PmeVariable, only : nfft, nffthalf, bsorder
    use PmeVariable, only : qfact
    use comm_pme_mod, only : ngxdiv, ngydiv, ngzdiv, imin, jmin, kmin
#ifndef PME_XYDIV
    use comm_pme_mod, only : ngxdivmax
#else
    use comm_pme_mod, only : ngzdivmax
#endif
    use omp_lib
    implicit none
    integer(4), intent(in)    :: natom
    integer(4)                :: imax, jmax, kmax
    integer(4), intent(in)    :: nomp
#ifndef PME_XYDIV
    real(kind=wrrp), intent(inout) :: wqgrid(ngxdivmax,ngydiv,ngzdiv,0:nomp-1)
#else
    ! difference with no PME_XYDIV case in declaration order of dimension is appropriate
    ! for FFTE code.
    real(kind=wrrp), intent(inout) :: wqgrid(ngzdivmax,ngydiv,ngxdiv,0:nomp-1)
#endif
    real(kind=dp),  intent(in)  :: bxyz(3, natom)
    real(kind=dp),  intent(in)  :: bcharge(natom)
    real(kind=wrrp), intent(in)  :: bthetai(3, 2, bsorder, natom)
    real(kind=wrrp) :: real_part
    real(kind=wrrp) :: charge
    integer(4) :: kd, ig1, ig2, ig3
    integer(4) :: igx, igy, igz
    integer(4) :: m1, m2, m3
    integer(4) :: i, iam

    imax = imin + ngxdiv
    jmax = jmin + ngydiv
    kmax = kmin + ngzdiv
    iam = 0

!$omp  parallel default(none) &
!$omp& private(i, kd, m1, m2, m3, charge, real_part, iam) &
!$omp& private(ig1, ig2, ig3, igx, igy, igz) &
!$omp& shared(natom, bxyz, nfft, nffthalf) &
!$omp& shared(imin, jmin, kmin, imax, jmax, kmax) &
!$omp& shared(qfact, wqgrid, bsorder, bthetai, bcharge)
!$  iam = omp_get_thread_num()

!! to superpose different qgrid contribution initialization is moved to caller side.
!!    wqgrid(:,:,:,iam)=Czero
!$omp do
    do i = 1, natom

       charge = qfact * bcharge(i)   ! for wave number 0 with no operation in reciprocal space.

       m1 = floor((bxyz(1, i) + 0.5) * dble(nfft(1))) + 1
       m2 = floor((bxyz(2, i) + 0.5) * dble(nfft(2))) + 1
       m3 = floor((bxyz(3, i) + 0.5) * dble(nfft(3))) + 1

#ifndef PME_XYDIV
!ocl novrec
       do ig3 = 0, bsorder - 1
          igz = m3 - ig3
          if(igz <= kmin .or. kmax < igz) cycle
          igz = igz - kmin
          do ig2 = 0, bsorder - 1
              igy = m2 - ig2
              if(igy <= jmin .or. jmax < igy) cycle
              igy = igy - jmin
              do ig1 = 0, bsorder - 1
                  igx = m1 - ig1
                  if(igx <= imin .or. imax < igx) cycle
                  igx = igx - imin
                  real_part = charge * bthetai(1, 1, bsorder - ig1, i) &
                                     * bthetai(2, 1, bsorder - ig2, i) &
                                     * bthetai(3, 1, bsorder - ig3, i) 
                  wqgrid(igx, igy, igz, iam) = wqgrid(igx, igy, igz, iam) + real_part
              enddo
          enddo
       enddo
#else
!ocl novrec
       do ig1 = 0, bsorder - 1
          igx = m1 - ig1
          if(igx <= imin .or. imax < igx) cycle
          igx = igx - imin
          do ig2 = 0, bsorder - 1
              igy = m2 - ig2
              if(igy <= jmin .or. jmax < igy) cycle
              igy = igy - jmin
              do ig3 = 0, bsorder - 1
                  igz = m3 - ig3
                  if(igz <= kmin .or. kmax < igz) cycle
                  igz = igz - kmin
                  real_part = charge * bthetai(1, 1, bsorder - ig1, i) &
                                     * bthetai(2, 1, bsorder - ig2, i) &
                                     * bthetai(3, 1, bsorder - ig3, i) 
                  wqgrid(igz, igy, igx, iam) = wqgrid(igz, igy, igx, iam) + real_part
              enddo
          enddo
       enddo
#endif
    enddo  ! i.
!$omp end do
!$omp end parallel

!$omp parallel default(shared) &
!$omp& private(igx,igy,igz,iam)
!$omp do collapse(2)
#ifndef PME_XYDIV
    do igz = 1, ngzdiv
       do igy = 1, ngydiv
          do iam = 1, nomp-1
             do igx = 1, ngxdiv
                wqgrid(igx,igy,igz,0)=wqgrid(igx,igy,igz,0) + wqgrid(igx,igy,igz,iam)
             enddo
          enddo
       enddo
    enddo
#else
    do igx = 1, ngxdiv
       do igy = 1, ngydiv
          do iam = 1, nomp-1
             do igz = 1, ngzdiv
                wqgrid(igz,igy,igx,0)=wqgrid(igz,igy,igx,0) + wqgrid(igz,igy,igx,iam)
             enddo
          enddo
       enddo
    enddo
#endif
!$omp end do
!$omp end parallel

  end subroutine CalculateQgrid
!------------------------------------------------------------------------------
  subroutine CalculatePmeReciprocalEnergy(qgrid, energy, virial)
    ! Calculate reciprocal energy and virials of PME.
    ! E = 0.5 * sum_{m_1 = 0}^{K_1 - 1} sum_{m_2 = 0}^{K_2 - 1} sum_{m_3 = 0}^{K_3 - 1}
    !     exp(-pi * pi * msv**2 / beta**2) / (msv**2 * pi * V) 
    !                                     * |b_1(m_1)|**2 * |b_2(m_2)|**2 * |b_3(m_3)|**2
    !     * F[Q](m_1, m_2, m_3) * F[Q](-m_1, -m_2, -m_3)
    ! F[Q]: Fourie transform of qgrid.

    use PmeMdmath,        only : PI
    use PmeVariable,      only : nfft, coefficient, ewald2Inverse
    use PmeVariable,      only : msv1, msv2, msv3
    use MdCellParameter , only : reciprocal_lattice_vector
#ifndef PME_XYDIV
    use comm_pme_mod,     only : ngxdivy, ngydivz, iminpy, jminpz
#else
    use comm_pme_mod,     only : ngzdivy, ngydivx, kminpy, jminpx
#endif
#ifdef DBG_PME_WNCOMP
    use comm_pme_mod,     only : work_dbg
#endif
    implicit none
#ifndef PME_XYDIV
    real(kind=wrrp), intent(inout) :: qgrid(2, nfft(3), ngxdivy, ngydivz)
#else
    real(kind=wrrp), intent(inout) :: qgrid(2, nfft(1), ngzdivy, ngydivx)
#endif
    integer(4) :: k1, k2, k3
    real(kind=wrrp) :: norm, grid_energy, m2, virial_coef
    real(kind=wrrp) :: mk1, mk2, mk3, mk2b, mk3b, mk3t, mm1, mm2, mm3
    real(kind=wrrp) :: qre, qim
    real(kind=wrrp), intent(inout) :: energy, virial(3,3)

    energy = 0.0d0
    virial(:, :) = 0.0d0

!$omp  parallel default(none) &
!$omp& private(k1, k2, k3) &
!$omp& private(m2, mk1, mk2, mk3, mk2b, mk3b, mk3t, mm1, mm2, mm3) &
!$omp& private(qre, qim, norm, grid_energy, virial_coef) &
!$omp& shared(qgrid, coefficient, ewald2Inverse) &
!$omp& shared(nfft) &
!$omp& shared(msv1, msv2, msv3, reciprocal_lattice_vector) &
#ifndef PME_XYDIV
!$omp& shared(ngxdivy, ngydivz, iminpy, jminpz) &
#else
!$omp& shared(ngzdivy, ngydivx, kminpy, jminpx) &
#endif 
#ifdef DBG_PME_WNCOMP
!$omp& shared(work_dbg) &
#endif 
!$omp& reduction(+:virial, energy)

#ifndef PME_XYDIV   /*--------------------------*/

!$omp do collapse(2)
    do k2=1,ngydivz
      do k1=1,ngxdivy
        do k3=1,nfft(3)
          if(k3*(k1+iminpy)*(k2+jminpz) == 1) cycle
          qre = qgrid(1, k3, k1, k2)
          qim = qgrid(2, k3, k1, k2)
          norm =  qre*qre + qim*qim
          grid_energy = Chalf * coefficient(k3,k1,k2) * norm
          energy = energy + grid_energy

          !! for backwards fft
          qgrid(1, k3, k1, k2) = qgrid(1, k3, k1, k2) * coefficient(k3,k1,k2)
          qgrid(2, k3, k1, k2) = qgrid(2, k3, k1, k2) * coefficient(k3,k1,k2)
          ! virial

#else     /*------------------------------------*/

!$omp do collapse(2)
    do k2=1,ngydivx
      do k3=1,ngzdivy
        do k1=1,nfft(1)
          if(k1*(k3+kminpy)*(k2+jminpx) == 1) cycle

          !! norm =  real(qgrid(k1, k2, k3) * conjg(qgrid(k1, k2, k3)))
          qre = qgrid(1, k1, k3, k2)
          qim = qgrid(2, k1, k3, k2)
          norm =  qre*qre + qim*qim
          grid_energy = Chalf * coefficient(k1,k3,k2) * norm
          energy = energy + grid_energy

          !! for backwards fft
          qgrid(1, k1, k3, k2) = qgrid(1, k1, k3, k2) * coefficient(k1,k3,k2)
          qgrid(2, k1, k3, k2) = qgrid(2, k1, k3, k2) * coefficient(k1,k3,k2)
          ! virial

#endif    /*------------------------------------*/

          mk1 = msv1(k1) * reciprocal_lattice_vector(1)
          mk2 = msv1(k1) * reciprocal_lattice_vector(2) + msv2(k2) * reciprocal_lattice_vector(4)
          mk3 = msv1(k1) * reciprocal_lattice_vector(3) + msv2(k2) * reciprocal_lattice_vector(5) &
              + msv3(k3) * reciprocal_lattice_vector(6)
          mk1 = Ctwo * PI * mk1
          mk2 = Ctwo * PI * mk2
          mk3 = Ctwo * PI * mk3

          mm1 = mk1*mk1; mm2 = mk2*mk2; mm3 = mk3*mk3
          m2 = mm1 + mm2 + mm3
          virial_coef = (Ctwo / m2) * (Cone + Cquarter * m2 * ewald2Inverse)
          virial(1, 1) = virial(1, 1) - grid_energy * (virial_coef * mm1 - Cone)
          virial(2, 2) = virial(2, 2) - grid_energy * (virial_coef * mm2 - Cone)
          virial(3, 3) = virial(3, 3) - grid_energy * (virial_coef * mm3 - Cone)
          virial(2, 1) = virial(2, 1) - grid_energy * virial_coef * mk2 * mk1
          virial(3, 1) = virial(3, 1) - grid_energy * virial_coef * mk3 * mk1
          virial(3, 2) = virial(3, 2) - grid_energy * virial_coef * mk3 * mk2

#ifdef DBG_PME_WNCOMP
#ifndef PME_XYDIV
          work_dbg(1,k3,k1+iminpy,k2+jminpz)= -grid_energy * (virial_coef * mm1 - Cone)
          work_dbg(2,k3,k1+iminpy,k2+jminpz)= -grid_energy * (virial_coef * mm2 - Cone)
          work_dbg(3,k3,k1+iminpy,k2+jminpz)= -grid_energy * (virial_coef * mm3 - Cone)
          work_dbg(4,k3,k1+iminpy,k2+jminpz)= -grid_energy * virial_coef * mk2 * mk1
          work_dbg(5,k3,k1+iminpy,k2+jminpz)= -grid_energy * virial_coef * mk3 * mk1
          work_dbg(6,k3,k1+iminpy,k2+jminpz)= -grid_energy * virial_coef * mk3 * mk2
          work_dbg(7,k3,k1+iminpy,k2+jminpz)=  grid_energy
#else
          work_dbg(1,k1,k3+kminpy,k2+jminpx)= -grid_energy * (virial_coef * mm1 - Cone)
          work_dbg(2,k1,k3+kminpy,k2+jminpx)= -grid_energy * (virial_coef * mm2 - Cone)
          work_dbg(3,k1,k3+kminpy,k2+jminpx)= -grid_energy * (virial_coef * mm3 - Cone)
          work_dbg(4,k1,k3+kminpy,k2+jminpx)= -grid_energy * virial_coef * mk2 * mk1
          work_dbg(5,k1,k3+kminpy,k2+jminpx)= -grid_energy * virial_coef * mk3 * mk1
          work_dbg(6,k1,k3+kminpy,k2+jminpx)= -grid_energy * virial_coef * mk3 * mk2
          work_dbg(7,k1,k3+kminpy,k2+jminpx)=  grid_energy
#endif
#endif

        enddo  ! k3/k1
      enddo  ! k1/k3
    enddo  ! k2
!$omp end do
!$omp end parallel
  end subroutine CalculatePmeReciprocalEnergy
!------------------------------------------------------------------------------
  subroutine CalculateForceReciprocal(bso, qgrid, &
                                      & natom, thetai, m2i, xyz, nakind, charge, w3_f, nomp)
    use PmeVariable,     only : nfft, nffthalf, bsorder
    use PmeVariable,     only : qfact
    use MdCellParameter, only : reciprocal_lattice_vector
    use subcell,         only : tag, na_per_cell
    use domain,          only : lxdiv, lydiv, lzdiv
    use param,           only : paranum
    use comm_pme_mod,    only : ngxdiv, ngydiv, ngzdiv, imin, jmin, kmin
    implicit none
    integer(4), intent(in) :: bso
    integer(4), intent(in) :: natom, nakind
    integer(4), intent(in) :: m2i(natom)
    integer(4), intent(in) :: nomp
#ifndef PME_XYDIV
    real(kind=wrrp), intent(in)  :: qgrid(2*bso+ngxdiv,2*bso+ngydiv,2*bso+ngzdiv)
#elif PME_XYDIV
    real(kind=wrrp), intent(in)  :: qgrid(2*bso+ngzdiv,2*bso+ngydiv,2*bso+ngxdiv)
#endif
    real(kind=wrrp), intent(in)  :: thetai(3, 2, bsorder, natom)
    real(kind=dp),  intent(in)  :: xyz(3, natom)
    real(kind=dp),  intent(in)  :: charge(nakind)
    real(kind=dp),  intent(inout) :: w3_f(3, natom, 0:nomp-1)
    integer(4) :: icx, icy, icz
    integer(4) :: i0, ipar
    integer(4) :: m1, m2, m3
    integer(4) :: ig1, ig2, ig3
    integer(4) :: igx, igy, igz
    real(kind=wrrp)  :: xi, yi, zi
    real(kind=wrrp)  :: rx, ry, rz
    real(kind=wrrp)  :: t11, t12, t21, t22, t31, t32
    real(kind=wrrp)  :: q, qq
    real(kind=wrap) :: fterm(3), fx, fy, fz

    ! bso=bsorder+bso_mgn
!$omp  parallel default(none) &
!$omp& private(q, qq) &
!$omp& private(icx, icy, icz) &
!$omp& private(i0, ipar, m1, m2, m3) &
!$omp& private(xi, yi, zi, rx, ry, rz) &
!$omp& private(ig1, ig2, ig3) &
!$omp& private(igx, igy, igz) &
!$omp& private(t11, t12, t21, t22, t31, t32) &
!$omp& private(fterm, fx, fy, fz) &
!$omp& shared(lxdiv, lydiv, lzdiv, tag, na_per_cell) &
!$omp& shared(qgrid, nfft, nffthalf) &
!$omp& shared(bsorder, bso) &
!$omp& shared(reciprocal_lattice_vector, thetai) &
!$omp& shared(xyz, charge, qfact, m2i, paranum, w3_f) &
!$omp& shared(imin, jmin, kmin)

!default: MTD
    do icx=1,lxdiv
       do icy=1,lydiv
!$omp do
          do i0=tag(    1,icy,icx),&
         &      tag(lzdiv,icy,icx)+na_per_cell(lzdiv,icy,icx)-1
            ipar=paranum(m2i(i0))
            q = charge(ipar) * qfact   ! for wave number 0 with no operation in reciprocal space.

            xi = xyz(1, i0)
            yi = xyz(2, i0)
            zi = xyz(3, i0)
            rx = xi*reciprocal_lattice_vector(1) + yi*reciprocal_lattice_vector(2) + zi*reciprocal_lattice_vector(3)
            ry = yi*reciprocal_lattice_vector(4) + zi*reciprocal_lattice_vector(5)
            rz = zi*reciprocal_lattice_vector(6)
            m1 = floor((rx + 0.5) * dble(nfft(1))) + 1
            m2 = floor((ry + 0.5) * dble(nfft(2))) + 1
            m3 = floor((rz + 0.5) * dble(nfft(3))) + 1
            fterm(:) = Czero
#ifndef PME_XYDIV
            do ig3 = 0, bsorder - 1
              igz = m3 - kmin - ig3 + bso
              t31 = thetai(3, 1, bsorder - ig3, i0)
              t32 = thetai(3, 2, bsorder - ig3, i0)
              do ig2 = 0, bsorder - 1
                igy = m2 - jmin - ig2 + bso
                t21 = thetai(2, 1, bsorder - ig2, i0)
                t22 = thetai(2, 2, bsorder - ig2, i0)
                do ig1 = 0, bsorder - 1
                  igx = m1 - imin - ig1 + bso
                  t12 = thetai(1, 2, bsorder - ig1, i0)
                  t11 = thetai(1, 1, bsorder - ig1, i0)
                  
                  qq = q * qgrid(igx, igy, igz)
                  fterm(1) = fterm(1) - qq * t12 * t21 * t31 * nfft(1)
                  fterm(2) = fterm(2) - qq * t11 * t22 * t31 * nfft(2)
                  fterm(3) = fterm(3) - qq * t11 * t21 * t32 * nfft(3)
                enddo
              enddo
            enddo
#elif PME_XYDIV
            do ig1 = 0, bsorder - 1
              igx = m1 - imin - ig1 + bso
              t12 = thetai(1, 2, bsorder - ig1, i0)
              t11 = thetai(1, 1, bsorder - ig1, i0)
              do ig2 = 0, bsorder - 1
                igy = m2 - jmin - ig2 + bso
                t21 = thetai(2, 1, bsorder - ig2, i0)
                t22 = thetai(2, 2, bsorder - ig2, i0)
                do ig3 = 0, bsorder - 1
                  igz = m3 - kmin - ig3 + bso
                  t31 = thetai(3, 1, bsorder - ig3, i0)
                  t32 = thetai(3, 2, bsorder - ig3, i0)
                  
                  qq = q * qgrid(igz, igy, igx)
                  fterm(1) = fterm(1) - qq * t12 * t21 * t31 * nfft(1)
                  fterm(2) = fterm(2) - qq * t11 * t22 * t31 * nfft(2)
                  fterm(3) = fterm(3) - qq * t11 * t21 * t32 * nfft(3)
                enddo
              enddo
            enddo
#endif
            fx = fterm(1) * reciprocal_lattice_vector(1)
            fy = fterm(1) * reciprocal_lattice_vector(2) + fterm(2) * reciprocal_lattice_vector(4)
            fz = fterm(1) * reciprocal_lattice_vector(3) + fterm(2) * reciprocal_lattice_vector(5) &
               + fterm(3) * reciprocal_lattice_vector(6)
            w3_f(1, i0, 0) = w3_f(1, i0, 0) + fx
            w3_f(2, i0, 0) = w3_f(2, i0, 0) + fy
            w3_f(3, i0, 0) = w3_f(3, i0, 0) + fz

          enddo  ! i0
!$omp end do
       end do
    end do
!$omp end parallel

  end subroutine CalculateForceReciprocal

end module PmeReciprocalFunctions
!------------------------------------------------------------------------------
!>
!! \brief Module to calculate far part of PME.
!! \author Kazushi Fujimoto and Shinichi Ichikawa
!<
module pme_far
  !precision module is loaded by default to declare PME subroutine
  use precision, only : dp, wrrp, wrap, MPI_PREC_PME_COM,  &
                      & Czero, Cdzero, Cone, Cdone, Ctwo, Cquarter
  use omp_lib
!!! nomp should have value thread num at initialize_openmp in opnemp_tool.
!!! but, nomp had not value threadnum when init_md_pmewald was called. it was 1.
!*  use openmp_tool, only : nomp
  use subcell
  use PmeReciprocalFunctions, only : CalculateBk2, CalculateCoefficient, BsplineFill,&
                                   & CalculateQgrid, CalculatePmeReciprocalEnergy, &
                                   & CalculateForceReciprocal, &
                                   & CalculateMsv

#include "timing.h90"

  implicit none

  integer(4) :: npoint
  real(kind=wrrp), allocatable :: thetai(:, :, :, :)
  real(kind=dp),  allocatable :: bfcharge(:)
  real(kind=dp),  allocatable :: bfxyz(:,:)
  real(kind=wrrp), allocatable :: bfthetai(:,:,:,:)

  private nomp
  integer(4) :: nomp=1  ! default

contains

!---------------------------------------------------------------------
!>
!! \brief  Subroutine to allocate arrays in meta-data style (PME)
!! \author Yoshimichi Andoh
!<
  subroutine fmod_alloc_pme_arraies
!----------------------------------------------------------------------
    use comm_pme_mod, only : maxorder
    use trajectory_mpi
    implicit none

  end subroutine fmod_alloc_pme_arraies
!----------------------------------------------------------------------
  subroutine init_md_pmewald()
!----------------------------------------------------------------------
    use trajectory_org
    use trajectory_mpi, only : nadirect
    use unit_cell, only      : cellx, celly, cellz
    use md_const
    use md_periodic
    use ewald_mod
    use physics_const, only : md_QQ_4PiE
    use PmeVariable, only : qfact
    use PmeVariable, only : ewald2Inverse
    use PmeVariable, only : nfft
    use PmeVariable, only : coefficient
    use PmeVariable, only : msv1, msv2, msv3

    use domain, only : ncellx, ncelly, ncellz
    use mpi_3d_grid
    use ewald_variables, only : ewald_alpha
    use comm_pme_mod
    use mpi_tool
    implicit none
    integer(4) :: i,infft
    integer(4) :: leng, icnt
    integer(4),allocatable :: kkx(:),kky(:),kkz(:)
#ifndef PME_XYDIV
!.. for yz-div of PME resulting in small number of cells along z-direction, hence shorter 
!.. inner-most loop length, in real space calc by MTD method.
    integer(4),allocatable :: kkxpy(:),kkypz(:)
#else
!.. for xy-div of PME resulting in maximum number of cells along z-direction, hence maximum
!.. inner-most loop length, in real space calc by MTD method.
    integer(4),allocatable :: kkzpy(:),kkypx(:)
#endif
    integer(4) :: bso
    include 'mpif.h'

!$  nomp = omp_get_max_threads()
    !
    if (md_periodic__type /= PMEWALD)  return
    !
    !     Set parameters for particle mesh mewald.
    !
#ifdef ONEPROC_AXIS
    !!! This restriction was caused by only old code:
    !!! periodic boundary condition to PME grid point derived from atom coordinate together with
    !!! conttribution of atoms in Halo result in doublefold amplitude of qgrid.
    !!!! This restriction can be removed. !!!!!!!
    !if(myrank==0) write(0,*) "ERROR:: ONEPROC_AXIS with PME is not", &
    !     &              " supported yet"
    !call modylas_abort()
#endif
    !!! This restriction has no reason. Number of sub-cells itsself has no relation with PME:
    !!! Both of FFTE and PZFFT3DV can handle FFT length multiple of power of 2, 3, and 5 numbers. 
    !!! Only restriction is that FFT grid point must be the multiple of number of processes.
    !!!! This restriction can be modefied. !!!!!!!
    !!!if(mod(ncellx,2)/=0 .or. mod(ncelly,2)/=0 .or. mod(ncellz,2)/=0)then
    !!!   if(myrank==0) write(0,*) "ERROR:: ncell=2powers*3powers is not supported yet"
    !!!   call modylas_abort()
    !!!endif
    if(mod(nfft1,npx)/=0 .or. mod(nfft2,npy)/=0 .or. mod(nfft3,npz)/=0 &
   &                     .or. mod(nfft1,npy)/=0 .or. mod(nfft2,npz)/=0 &   ! YZ-division in reciprocal.
   &                     .or. mod(nfft3,npy)/=0 .or. mod(nfft2,npx)/=0) then  ! YX-division.
       if(myrank==0) write(0,*) "ERROR:: FFT grid is not the multiple of number of processes"
       call modylas_abort()
    endif
#ifdef PRECISION_PME_SP
    !!! The reason of this restriction is only in "pzfft3dv" code.
    !!! The precision of "pzfft3dv" is double precision, its old standard fortran makes 
    !!! use of "use statement" for precision impossible.
    if(myrank==0) write(0,*) "ERROR:: PRECISION_PME_SP with PME is not", &
         &              " supported yet. Please consult development group of Modylas."
    call modylas_abort()
#endif
    
    ! ewald
    ewald2Inverse = Cdone / (ewald_alpha * ewald_alpha)
    ! unit of ewald_alpha in modylas is 1/m, whereas its 1/nm in PME.
    ! ewald2Inverse = Cone / (ewald_alpha * ewald_alpha * 1.0d+18) when ewald_alpha of PME_Math.
    !
    npoint = nfft1 * nfft2 * nfft3
    
    bso = bsorder+bso_mgn                ! bso_mgn=6, in module margin_sizes, through comm_pme_mod.
    !
    !     Calculate MPI parallel division variables.
    !
    infft=nfft1; call value_check(infft)
    infft=nfft2; call value_check(infft)
    infft=nfft3; call value_check(infft)

    ! process decomposition.
    allocate( kkx(0:npx-1) ) ; kkx=0
    allocate( kky(0:npy-1) ) ; kky=0
    allocate( kkz(0:npz-1) ) ; kkz=0
#ifndef PME_XYDIV
    allocate( kkxpy(0:npy-1) ) ; kkxpy=0
    allocate( kkypz(0:npz-1) ) ; kkypz=0
#else
    allocate( kkzpy(0:npy-1) ) ; kkzpy=0
    allocate( kkypx(0:npx-1) ) ; kkypx=0
#endif

    icnt=0
    do i=1,nfft1
       kkx(icnt)=kkx(icnt)+1
       icnt=icnt+1
       if(icnt==npx) icnt=0
    enddo

    icnt=0
    do i=1,nfft2
       kky(icnt)=kky(icnt)+1
       icnt=icnt+1
       if(icnt==npy) icnt=0
    enddo

    icnt=0
    do i=1,nfft3
       kkz(icnt)=kkz(icnt)+1
       icnt=icnt+1
       if(icnt==npz) icnt=0
    enddo
#ifndef PME_XYDIV
    icnt=0
    do i=1,nfft1
       kkxpy(icnt)=kkxpy(icnt)+1
       icnt=icnt+1
       if(icnt==npy) icnt=0
    enddo

    icnt=0
    do i=1,nfft2
       kkypz(icnt)=kkypz(icnt)+1
       icnt=icnt+1
       if(icnt==npz) icnt=0
    enddo
#else
    icnt=0
    do i=1,nfft3
       kkzpy(icnt)=kkzpy(icnt)+1
       icnt=icnt+1
       if(icnt==npy) icnt=0
    enddo

    icnt=0
    do i=1,nfft1
       kkypx(icnt)=kkypx(icnt)+1
       icnt=icnt+1
       if(icnt==npx) icnt=0
    enddo
#endif

    ngxdiv = kkx(ipx)
    ngydiv = kky(ipy)
    ngzdiv = kkz(ipz)
    ngxmin = minval(kkx(:))
    ngymin = minval(kky(:))
    ngzmin = minval(kkz(:))
#ifndef PME_XYDIV
    ngxdivmax = (nfft1-1)/npx + 1
    ngxdivy = kkxpy(ipy)
    ngydivz = kkypz(ipz)
#else
    ngzdivmax = (nfft3-1)/npz + 1
    ngzdivy = kkzpy(ipy)
    ngydivx = kkypx(ipx)
#endif

    imin=0; jmin=0; kmin=0
    do i=0,ipx-1
       imin=imin+kkx(i)
    enddo
    do i=0,ipy-1
       jmin=jmin+kky(i)
    enddo
    do i=0,ipz-1
       kmin=kmin+kkz(i)
    enddo
#ifndef PME_XYDIV
    iminpy=0; jminpz=0
    do i=0,ipy-1
       iminpy=iminpy+kkxpy(i)
    enddo
    do i=0,ipz-1
       jminpz=jminpz+kkypz(i)
    enddo
#else
    kminpy=0; jminpx=0
    do i=0,ipy-1
       kminpy=kminpy+kkzpy(i)
    enddo
    do i=0,ipx-1
       jminpx=jminpx+kkypx(i)
    enddo
#endif

    !
    !     Allocate work arrays.
    !

#ifndef PME_XYDIV
    allocate(work0((nfft1-1)/npx+1,ngydiv,ngzdiv))
    allocate(work0_omp((nfft1-1)/npx+1,ngydiv,ngzdiv,0:nomp-1))
    allocate(mpitmp((nfft1-1)/npx+1,ngydiv,ngzdiv,0:npx-1))
    allocate(work1(2,((nfft1-1)/npx+1)*npx,ngydiv,ngzdiv))
    allocate(work2(2,((nfft3-1)/npx+1)*npx,ngxdivy,ngydivz))
    leng = max(2*bso*ngydiv*ngzdiv,(2*bso+ngxdiv)*2*bso*ngzdiv)
#else
    allocate(work0((nfft3-1)/npz+1,ngydiv,ngxdiv))
    allocate(work0_omp((nfft3-1)/npz+1,ngydiv,ngxdiv,0:nomp-1))
    allocate(mpitmp ((nfft3-1)/npz+1,ngydiv,ngxdiv,0:npz-1))
    allocate(work1(2,((nfft3-1)/npz+1)*npz,ngydiv,ngxdiv))
    allocate(work2(2,((nfft1-1)/npz+1)*npz,ngzdivy,ngydivx))
    leng = max(2*bso*ngydiv*ngxdiv,(2*bso+ngzdiv)*2*bso*ngxdiv)
#endif
    allocate(buff(leng,2))
    allocate(buff2(leng,2))
    work0=Czero
    work0_omp=Czero
    mpitmp=Czero
    work1=Czero
    work2=Czero
#ifdef DBG_PME_WNCOMP
#ifndef PME_XYDIV
    allocate(work_dbg(7,nfft3,nfft1,nfft2))
#else
    allocate(work_dbg(7,nfft1,nfft3,nfft2))
#endif
#endif

    ! ... PME grid size variables for reciprocal calc.
    nfft(1) = nfft1
    nfft(2) = nfft2
    nfft(3) = nfft3

    !
    !     initialize communicator.
    !
    call comm_create()
    !
    !     allocate PME arrays
    !
    !--- charge constatnt factor for use in qgrid mapping.
    !--- for wave number 0 with no operation in reciprocal space.
    qfact = dsqrt(md_QQ_4PiE)

#ifndef PME_XYDIV
    allocate(coefficient(nfft(3), ngxdivy, ngydivz))
    allocate(msv3(nfft(3) ))
    allocate(msv1(ngxdivy))
    allocate(msv2(ngydivz))
#else
    allocate(coefficient(nfft(1), ngzdivy, ngydivx))
    allocate(msv1(nfft(1) ))
    allocate(msv3(ngzdivy))
    allocate(msv2(ngydivx))
#endif
    coefficient(:,:,:) = Czero
    msv1 = Czero; msv2 = Czero; msv3 = Czero

    ! qgrid for qgrid communication and force calculation in which imaginary part is not necessary.
    ! qgrid is also in calculate_energy_reciprocal. In this case qgrid has imaginary part.
#ifndef PME_XYDIV
    allocate(qgrid(2*bso+ngxdiv,2*bso+ngydiv,2*bso+ngzdiv))   ! qgrid is in module comm_pme and is allocated in init_md_pmewald.
#else
    allocate(qgrid(2*bso+ngzdiv,2*bso+ngydiv,2*bso+ngxdiv))   ! qgrid is in module comm_pme and is allocated in init_md_pmewald.
#endif
    qgrid = Czero

    allocate(bfcharge(nadirect))
    allocate(bfxyz(3,nadirect))
    bfcharge = Czero
    bfxyz = Cdzero

    allocate(thetai(3, 2, bsorder, nadirect))
    thetai = Czero
    allocate(bfthetai(3, 2, bsorder, nadirect))
    bfthetai = Czero

    !
    !     Calculate coefficients for reciprocal calculation.
    !
    call CalculateBk2(bsorder, maxorder)  ! maxorder=8 and  bsorder=4 in module comm_pme_mod.
    !
    call CalculateMsv()
    !
    deallocate(kkx,kky,kkz)
#ifndef PME_XYDIV
    deallocate(kkxpy)
    deallocate(kkypz)
#else
    deallocate( kkzpy )
    deallocate( kkypx )
#endif
    !
    !     calculate constant part of potential energy in Ewald's method
    !
    call md_calculate_ewald_e_constant()
    !
  end subroutine init_md_pmewald

!-----------------------------------------------------------------------
  subroutine md_add_pmewald()

    use trajectory_mpi,  only : wkxyz, nadirect
    use subcell,         only : tag, na_per_cell, m2i
    use domain,          only : lxdiv, lydiv,lzdiv
    use MdCellParameter, only : GetCellParameter, reciprocal_lattice_vector
    use atom_virial
    use md_const        !!!! several doubled declarations with PME modules
    use coulomb_mod,     only : chgv
    use forces,          only : w3_f
    use md_monitors
    use domain,          only : lxdiv, lydiv, lzdiv
    use param
    use mpi_3d_grid
    use unit_cell
    use cell_shape
    use comm_pme_mod,    only : bsorder, bso_mgn
    use comm_pme_mod,    only : nfft1, nfft2, nfft3
    use comm_pme_mod,    only : imin, jmin, kmin
    use comm_pme_mod,    only : ngxdiv,ngydiv,ngzdiv
    use comm_pme_mod,    only : ngxmin,ngymin,ngzmin
    use comm_pme_mod,    only : work0_omp
    use comm_pme_mod,    only : work1, work2
    use comm_pme_mod,    only : mpitmp
    use comm_pme_mod,    only : qgrid
#ifndef PME_XYDIV
    use comm_pme_mod,    only : ngxdivmax
    use comm_pme_mod,    only : iminpy, jminpz, ngxdivy, ngydivz
    use comm_pme_mod,    only : icommx, icommy_ffte, icommz_ffte
    use comm_pme_mod,    only : myrankx, nrankx
#else
    use comm_pme_mod,    only : ngzdivmax
    use comm_pme_mod,    only : kminpy, jminpx, ngzdivy, ngydivx
    use comm_pme_mod,    only : icommz, icommy_ffte, icommx_ffte
    use comm_pme_mod,    only : myrankz, nrankz
#endif
    use comm_pme_mod,    only : comm_qgrid
    use mpi_tool
    use dr_cntl,         only : nbd
#ifdef DBG_PME_WNCOMP
    use md_condition,    only : mdstep
    use comm_pme_mod,    only : work_dbg
#endif

    implicit none
    integer(4) :: ipar
    integer(4) :: lr,icz,icy,icx
    integer(4) :: i0, noff, na, i0s, i0fin
    integer(4) :: igz,igy,igx
    integer(4) :: bso
    integer(4) :: nakind
    integer(4) :: natom, nblk, leng, nbmax, iam
    real(kind=wrrp) :: energy, virial(3, 3)
!!    real(kind=wrap) :: energy, virial(3, 3)
    ! FOR DEBUG
    integer(4) :: i, j, k
    real(kind=dp) :: xi, yi, zi
    include 'mpif.h'
    integer(4) :: ierr
#ifdef DEBUGFCE
    real(kind=wrap) :: epme
    epme = Czero
#endif
#ifdef DBG_PME_WNCOMP
    real(kind=wrp) :: ecpme, vl_cpme, rd_vij(6)
    integer(4)     :: k1,k2,k3
    ecpme = Czero
    vl_cpme = Czero
#endif

    TIME_START(TM_ADD_PMEWALD)

    energy = Cdzero
    virial(:,:) = Cdzero

    TIMED_START(TM_PME_COEFF)
    ! ... reciprocal vector and cell parameters.
    call GetCellParameter()

    ! ... coefficient for reciprocal energy.
    call CalculateCoefficient()

    TIMED_STOP(TM_PME_COEFF)
    TIMED_START(TM_PME_BSPFILL)
    ! ... generate thetai.
    CALL BsplineFill(thetai)

    TIMED_STOP(TM_PME_BSPFILL)
    TIMED_START(TM_PME_CALC_QGRID)
    ! ... pack atom data for qgrid calculation.
    ! note: Atom data buffer bfxyz, bfcharge, and bfthetai etc. MUST INCLUDE that of cells, 
    !       thickness 1 cell, outside of process boundary. This is to avoid communication
    !       of qgrid for superposing contribution from atoms in adjacent process.
    !
    ! find maximum number per thread of atoms in Z-directional cell blocks.
    nbmax = 0; natom = 0
!default: MTD
    do icx = 1-nbd, lxdiv+nbd
    do icy = 1-nbd, lydiv+nbd
      leng = tag(lzdiv+nbd,icy,icx)+na_per_cell(lzdiv+nbd,icy,icx) - tag(1-nbd,icy,icx)
      nblk = (leng - 1)/nomp + 1
      nbmax = max(nbmax,nblk)
      natom = natom + leng
    enddo
    enddo
    nblk = nbmax

    ! buffering.
!$omp  parallel default(none) &
!$omp& private(icx, icy, i0s, i0, i0fin, noff, na, ipar, iam) &
!$omp& private(xi, yi, zi) &
!$omp& shared(nblk, lxdiv, lydiv, lzdiv, tag, na_per_cell) &
!$omp& shared(paranum, m2i, chgv, wkxyz, thetai, reciprocal_lattice_vector) &
!$omp& shared(bfcharge, bfxyz, bfthetai)
!$  iam = omp_get_thread_num()
    noff = 0
!default: MTD
    do icx = 1-nbd, lxdiv+nbd
    do icy = 1-nbd, lydiv+nbd
      i0fin = tag(lzdiv+nbd,icy,icx)+na_per_cell(lzdiv+nbd,icy,icx)-1
!$omp do
      do i0s=tag(1-nbd,icy,icx), &
         &      tag(lzdiv+nbd,icy,icx)+na_per_cell(lzdiv+nbd,icy,icx)-1, nblk
      na = noff + nblk*iam
      do i0=i0s, min(i0s+nblk-1, i0fin)
        na = na + 1
        ipar=paranum(m2i(i0))
        bfcharge(na) = chgv(ipar)
        xi = wkxyz(1,i0)
        yi = wkxyz(2,i0)
        zi = wkxyz(3,i0)
        bfxyz(1, na) = xi*reciprocal_lattice_vector(1) + yi*reciprocal_lattice_vector(2) + zi*reciprocal_lattice_vector(3) ! ri*kn
        bfxyz(2, na) = yi*reciprocal_lattice_vector(4) + zi*reciprocal_lattice_vector(5)
        bfxyz(3, na) = zi*reciprocal_lattice_vector(6)
        bfthetai(:,:,:,na) = thetai(:,:,:,i0)
      enddo  ! i0.
      enddo  ! i0s.
!$omp end do
!default: MTD
      noff = noff + (tag(lzdiv+nbd,icy,icx)+na_per_cell(lzdiv+nbd,icy,icx) &
      &       - tag(1-nbd,icy,icx))
    enddo
    enddo
!$omp end parallel

    ! initiaization of qgrid working array.
    work0_omp = Czero

    !CALL calculate_qgrid(natom, bfxyz, bfcharge, bfthetai, work0_omp, nomp)
    CALL CalculateQgrid(natom, bfxyz, bfcharge, bfthetai, work0_omp, nomp)

    TIMED_STOP(TM_PME_CALC_QGRID)
    TIME_STOP(TM_ADD_PMEWALD)
#ifdef WTIMER
    time2=MPI_WTIME()
    write(on,*) '_qreductn', time2-time1
    time1=MPI_WTIME()
#endif
    TIME_BARRIER(TMB_ADD_PMEWALD_GATHER)
    TIME_START(TM_ADD_PMEWALD_GATHER)

#ifndef PME_XYDIV
!... next does not work correctly for the case ngxdiv != (nfft1-1)/npx+1 except for last rank.
!... must be ngxdiv=((nfft1-1)/npx+1).
    CALL MPI_GATHER(work0_omp, ((nfft1-1)/npx+1)*ngydiv*ngzdiv, &
         &                MPI_PREC_PME_COM, &
         &                mpitmp, ((nfft1-1)/npx+1)*ngydiv*ngzdiv, &
         &                MPI_PREC_PME_COM, 0, icommx, ierr)
#else
!... next does not work correctly for the case ngzdiv != (nfft3-1)/npz+1 except for last rank.
!... must be ngzdiv=((nfft3-1)/npz+1).
    CALL MPI_GATHER(work0_omp, ((nfft3-1)/npz+1)*ngydiv*ngxdiv, &
         &                MPI_PREC_PME_COM, &
         &                mpitmp,  ((nfft3-1)/npz+1)*ngydiv*ngxdiv, &
         &                MPI_PREC_PME_COM, 0, icommz, ierr)
#endif

    TIME_STOP(TM_ADD_PMEWALD_GATHER)
#ifdef WTIMER
    time2=MPI_WTIME()
    write(on,*) '_mpigathr', time2-time1
    time1=MPI_WTIME()
#endif
    TIME_BARRIER(TMB_PMEWALD_PZFFT3DV)

#ifndef PME_XYDIV
    IF (myrankx == 0) THEN
#else
    IF (myrankz == 0) THEN
#endif

       TIME_START(TM_ADD_PMEWALD)
!$omp  parallel default(none) &
!$omp& private(lr,igx,igy,igz) &
!$omp& shared(mpitmp,work1) &
#ifndef PME_XYDIV
!$omp& shared(ngydiv,ngzdiv,nfft1,npx)
!... next does not work correctly for the case ngxdiv != (nfft1-1)/npx+1 except for last rank.
       DO lr = 0, npx - 1
!$omp do
          DO igz = 1, ngzdiv
             DO igy = 1, ngydiv
                DO igx = 1, (nfft1-1)/npx+1
                   work1(1,((nfft1-1)/npx+1)*lr+igx,igy,igz) = mpitmp(igx,igy,igz,lr)
                   work1(2,((nfft1-1)/npx+1)*lr+igx,igy,igz) = Czero
                END DO
             END DO
          END DO
!$omp end do
       END DO
#else
!$omp& shared(ngxdiv,ngydiv,nfft3,npz)
!... next does not work correctly for the case ngzdiv != (nfft3-1)/npz+1 except for last rank.
       DO lr = 0, npz - 1
!$omp do
          DO igx = 1, ngxdiv
             DO igy = 1, ngydiv
                DO igz = 1, (nfft3-1)/npz+1
                   ! IF (((nfft3-1)/npz+1)*lr+igz > nfft3) cycle
                   work1(1,((nfft3-1)/npz+1)*lr+igz,igy,igx) = mpitmp(igz,igy,igx,lr)
                   work1(2,((nfft3-1)/npz+1)*lr+igz,igy,igx) = Czero
                END DO
             END DO
          END DO
!$omp end do
       END DO
#endif
!$omp end parallel

       TIME_STOP(TM_ADD_PMEWALD)
       TIME_START(TM_PMEWALD_PZFFT3DV)
       !
       !.... FFT forward transform
       !
#ifndef PME_XYDIV
       CALL PZFFT3DV(work1, work2, nfft1, nfft2, nfft3, &
            &                 icommy_ffte, icommz_ffte, npy ,npz, -2)
#else
       CALL PZFFT3DV(work1, work2, nfft3, nfft2, nfft1, &
            &                 icommy_ffte, icommx_ffte, npy ,npx, -2)
#endif

       TIME_STOP(TM_PMEWALD_PZFFT3DV)

#ifdef DBG_PME_WNCOMP
       work_dbg=0.0d0
#endif

      !
      !.... reciprocal space calculation
      !
      TIME_START(TM_ADD_PMEWALD)
      TIMED_START(TM_PME_ENERGY)

      CALL CalculatePmeReciprocalEnergy(work2, energy, virial)

      TIMED_STOP(TM_PME_ENERGY)
      TIME_STOP(TM_ADD_PMEWALD)
      TIME_START(TM_PMEWALD_PZFFT3DV)
      !
      !.... FFT inverse transform
      !
#ifndef PME_XYDIV
       CALL PZFFT3DV(work2, work1, nfft1, nfft2, nfft3, &
            &                 icommy_ffte, icommz_ffte, npy, npz, +2)
#else
       CALL PZFFT3DV(work2, work1, nfft3, nfft2, nfft1, &
            &                 icommy_ffte, icommx_ffte, npy, npx, +2)
#endif
       TIME_STOP(TM_PMEWALD_PZFFT3DV)
       TIME_START(TM_ADD_PMEWALD)
       !
       !
!$omp  parallel default(none) &
!$omp& private(lr,igx,igy,igz) &
!$omp& shared(mpitmp,work1) &
#ifndef PME_XYDIV
!$omp& shared(ngydiv,ngzdiv,nfft1,npoint,npx)
       !.... SCATTER X-direction
       !... next does not work correctly for the case ngxdiv != (nfft1-1)/npx+1 except for last rank.
       DO lr = 0, npx - 1
!$omp do
          DO igz = 1, ngzdiv
             DO igy = 1, ngydiv
                DO igx = 1, (nfft1-1)/npx+1
                   mpitmp(igx,igy,igz,lr) = work1(1,((nfft1-1)/npx+1)*lr+igx,igy,igz)*npoint
#else
!$omp& shared(ngxdiv,ngydiv,nfft3,npoint,npz)
       !.... SCATTER Z-direction
       !... next does not work correctly for the case ngzdiv != (nfft3-1)/npz+1 except for last rank.
       DO lr = 0, npz - 1
!$omp do
          DO igx = 1, ngxdiv
             DO igy = 1, ngydiv
                DO igz = 1, (nfft3-1)/npz+1
                   ! IF (((nfft3-1)/npz+1)*lr+igz > nfft3) cycle
                   mpitmp (igz,igy,igx,lr) = work1(1,((nfft3-1)/npz+1)*lr+igz,igy,igx)*npoint
#endif
                END DO
             END DO
          END DO
!$omp end do
       END DO
!$omp end parallel

       TIME_STOP(TM_ADD_PMEWALD)

    END IF !! myrankx==0/myrankz==0

#ifdef DBG_PME_WNCOMP
#ifndef PME_XYDIV
    if(myrankx==0) then
!correct only when nfft1/nfft3 is the multiple of npy, nfft2 is the multiple of npz/npx
    do k2=jminpz+1, jminpz+ngydivz
      call mpi_gather(work_dbg(1,1,iminpy+1,k2), 7*nfft3*ngxdivy, &
         &            MPI_PREC_PME_COM,                       &
         &            work_dbg(1,1,1       ,k2), 7*nfft3*ngxdivy, &
         &            MPI_PREC_PME_COM, 0,   icommy_ffte, ierr )
    end do
      call mpi_gather(work_dbg(1,1,1,jminpz+1), 7*nfft3*nfft1*ngydivz, &
         &            MPI_PREC_PME_COM,                              &
         &            work_dbg(1,1,1       ,1), 7*nfft3*nfft1*ngydivz, &
         &            MPI_PREC_PME_COM, 0,  icommz_ffte, ierr         )
    endif  ! myrankx=0.
#else
    if(myrankz==0) then
!correct only when nfft1/nfft3 is the multiple of npy, nfft2 is the multiple of npz/npx
    do k2=jminpx+1, jminpx+ngydivx
      call mpi_gather(work_dbg(1,1,kminpy+1,k2), 7*nfft1*ngzdivy, &
         &            MPI_PREC_PME_COM,                       &
         &            work_dbg(1,1,1       ,k2), 7*nfft1*ngzdivy, &
         &            MPI_PREC_PME_COM, 0,   icommy_ffte, ierr )
    end do
      call mpi_gather(work_dbg(1,1,1,jminpx+1), 7*nfft1*nfft3*ngydivx, &
         &            MPI_PREC_PME_COM,                              &
         &            work_dbg(1,1,1       ,1), 7*nfft1*nfft3*ngydivx, &
         &            MPI_PREC_PME_COM, 0,  icommx_ffte, ierr         )
    endif  ! myrankz=0.
#endif

    if(myrank == 0) then
      ecpme=0.0d0; rd_vij=0.0d0
!ocl loop_nointerchange
!ocl nosimd
      do k3=1, nfft3
      do k2=1, nfft2
      do k1=1, nfft1
#ifndef PME_XYDIV
        rd_vij(1) = rd_vij(1) + work_dbg(1,k3,k1,k2)
        rd_vij(2) = rd_vij(2) + work_dbg(2,k3,k1,k2)
        rd_vij(3) = rd_vij(3) + work_dbg(3,k3,k1,k2)
        rd_vij(4) = rd_vij(4) + work_dbg(4,k3,k1,k2)
        rd_vij(5) = rd_vij(5) + work_dbg(5,k3,k1,k2)
        rd_vij(6) = rd_vij(6) + work_dbg(6,k3,k1,k2)
        ecpme = ecpme + work_dbg(7,k3,k1,k2)
#else
        rd_vij(1) = rd_vij(1) + work_dbg(1,k1,k3,k2)
        rd_vij(2) = rd_vij(2) + work_dbg(2,k1,k3,k2)
        rd_vij(3) = rd_vij(3) + work_dbg(3,k1,k3,k2)
        rd_vij(4) = rd_vij(4) + work_dbg(4,k1,k3,k2)
        rd_vij(5) = rd_vij(5) + work_dbg(5,k1,k3,k2)
        rd_vij(6) = rd_vij(6) + work_dbg(6,k1,k3,k2)
        ecpme = ecpme + work_dbg(7,k1,k3,k2)
#endif
      end do
      end do
      end do

      write(*,*) "DBG_PME_WNSERI mdstep: ",mdstep," global EC: ",ecpme
      write(*,*) "DBG_PME_WNSERI mdstep: ",mdstep," global Vl_c_PME: ",vl_cpme
      write(*,*) "DBG_PME_WNSERI mdstep: ",mdstep," global wk_v11: ",rd_vij(1)
      write(*,*) "DBG_PME_WNSERI mdstep: ",mdstep," global wk_v22: ",rd_vij(2)
      write(*,*) "DBG_PME_WNSERI mdstep: ",mdstep," global wk_v33: ",rd_vij(3)
      write(*,*) "DBG_PME_WNSERI mdstep: ",mdstep," global wk_v21: ",rd_vij(4)
      write(*,*) "DBG_PME_WNSERI mdstep: ",mdstep," global wk_v31: ",rd_vij(5)
      write(*,*) "DBG_PME_WNSERI mdstep: ",mdstep," global wk_v32: ",rd_vij(6)
    endif  ! myrank=/ca0.

#endif   /* DBG_PME_WNCOMP */

#ifdef WTIMER
    time2=MPI_WTIME()
    write(on,*) '_pzfft3+-', time2-time1
    call mpi_barrier(mpi_comm_world,ierr)     
    time1=MPI_WTIME()
#endif

    TIME_BARRIER(TMB_COMM_QGRID)
    TIME_START(TM_COMM_QGRID)

#ifndef PME_XYDIV
    CALL comm_qgrid(ngxdivmax,ngxdiv,ngydiv,ngzdiv,ngxmin,ngymin,ngzmin,icommx)
#else
    CALL comm_qgrid(ngzdivmax,ngzdiv,ngydiv,ngxdiv,ngzmin,ngymin,ngxmin,icommz)
#endif

    TIME_STOP(TM_COMM_QGRID)
    TIME_START(TM_ADD_PMEWALD)
    TIMED_START(TM_PME_FORCE)

#ifdef WTIMER
    time2=MPI_WTIME()
    write(on,*) '_commqgrd', time2-time1
    time1=MPI_WTIME()
#endif

    bso=bsorder+bso_mgn        ! bso_mgn=6, in module margin_sizes, through comm_pme_mod.
    nakind = size(chgv)
    call CalculateForceReciprocal(bso, qgrid, &
                                  & nadirect, thetai, m2i, wkxyz, nakind, chgv, w3_f, nomp)

    TIMED_STOP(TM_PME_FORCE)
#ifdef WTIMER
    time2=MPI_WTIME()
    write(on,*) '_calforce', time2-time1
#endif
    wk_p_energy = wk_p_energy + energy
    wk_vir2(1,0)=wk_vir2(1,0)+virial(1,1) ! wk_v11
    wk_vir2(2,0)=wk_vir2(2,0)+virial(2,2) ! wk_v22
    wk_vir2(3,0)=wk_vir2(3,0)+virial(3,3) ! wk_v33
    wk_vir2(4,0)=wk_vir2(4,0)+virial(2,1) ! wk_v21
    wk_vir2(5,0)=wk_vir2(5,0)+virial(3,1) ! wk_v31
    wk_vir2(6,0)=wk_vir2(6,0)+virial(3,2) ! wk_v32
    TIME_STOP(TM_ADD_PMEWALD)
#ifdef DEBUGFCE
    epme = 0.0d0
    call mpi_allreduce(energy,epme,1, &
         &     MPI_PREC_PME_COM,mpi_sum,mpi_comm_world,ierr)
#ifdef KCAL
    if(myrank==0) write(*,*) 'Pot(PME  )=',epme*kJ_mol/4.184d0,'[kcal/mol]'
#else
    if(myrank==0) write(*,*) 'Pot(PME  )=',epme*kJ_mol,'[kJ/mol]'
#endif
#endif

  end subroutine md_add_pmewald
!----------------------------------------------------------------------
#ifdef WTIMER
  real function etime(time)
    real, intent(out) :: time(2)
    call cpu_time(etime)
    time(1) = etime
    time(2) = 0
  end function etime
!!
!     SUBROUTINE CPU_TIME(TIME)
!     REAL, INTENT(OUT) :: TIME

!     INTEGER :: COUNT, COUNT_RATE

!     CALL SYSTEM_CLOCK(COUNT,COUNT_RATE)
!     IF (COUNT_RATE /= 0) THEN
!       TIME = REAL(COUNT)/COUNT_RATE
!     ELSE
!       TIME = -1
!     END IF
!     END SUBROUTINE CPU_TIME
#endif
!	real function dtime(time)
!	real function etime(time)
!	real time(2)
!	double precision,save :: last_time = 0
!	double precision this_time
!	call cpu_time(this_time)
!	time(1) = this_time - last_time
!	time(2) = 0
!	dtime = time(1)
!	last_time = this_time
!	end
!----------------------------------------------------------------------
!#endif
!-------------------------------------------------------------------------
  subroutine value_check(ivalue)
    implicit none
    integer, intent(inout) :: ivalue
    integer :: iret, ivalue0

    ivalue0=ivalue

    call factor_check(2,ivalue,iret)
    if(iret.lt.0) call factor_check(3,ivalue,iret)
    if(iret.lt.0) call factor_check(5,ivalue,iret)
    if(iret.lt.0) then
       write(0,*) 'ERROR: nfft=',ivalue0,' must be 2**x or 3**x or 5**x'
       call modylas_abort()
    endif
  end subroutine value_check
!-------------------------------------------------------------------------
  subroutine factor_check(idiv,inumber,iflag)
    implicit none
    integer, intent(in) :: idiv
    integer, intent(inout) :: inumber
    integer, intent(inout) :: iflag

    do
       iflag=mod(inumber,idiv)
       if(iflag.eq.0) then
          inumber=inumber/idiv
       else
          if(inumber.eq.1) then
             iflag=0
          else
             iflag=-1
          endif
          exit
       endif
    enddo
  end subroutine factor_check
!----------------------------------------------------------------------
  subroutine check_pmewald_alpha(val)
    use ewald_variables
    use mpi_tool
    implicit none
    real(kind=dp), intent(in) :: val
    if(val.eq.0)then
       write(0,*) &
  &    'ERROR: alpha for Ewald/PME should not be 0. Check your mddef.'
       call modylas_abort()
    endif
  end subroutine check_pmewald_alpha

end module pme_far

