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
!! \brief  Module to calculate coefficients of potential and force in L2P transformation of FMM.
!<
!----------------------------------------------------------------------
!>
!! \brief  Calculating coefficients of potential and force (from regular solid harmonics in cartesian coordinates computed by utilizing recurrence relation).
!! \author Tatsuya Sakashita, Noriyuki Yoshii.
!<
module regular_solid_harmonics_cartesian
  use math_functions
  implicit none

contains

  subroutine calculate_regular_harmonics_array(n_max, factorial_array, minus_double_factorial_array, &
       & x, y, z, array)
    implicit none
    integer, intent(in) :: n_max
    real(8), intent(in) :: factorial_array(0:2*n_max), minus_double_factorial_array(0:n_max)
    real(8), intent(in) :: x, y, z
    complex(8), intent(out) :: array(0:n_max, 0:n_max)
    real(8) :: r_square
    real(8) :: fact,pll,pmm,pmmp1
    integer :: l, m
    real(8) :: re, im
    complex(8) :: pow_xy

    r_square = x**2 + y**2 + z**2
    
    pow_xy = dcmplx(1d0, 0)

    do m=0,n_max-1
       pmm = minus_double_factorial_array(m)
       array(m,m) = pmm * pow_xy / factorial_array(2*m) ! When l=m
       pmmp1=z*(2*m+1)*pmm
       array(m+1,m) = pmmp1 * pow_xy / factorial_array(2*m+1) ! When l=m+1
       do l=m+2, n_max ! for first index of array
          pll=(z*(2*l-1)*pmmp1-(l+m-1)*r_square*pmm)/(l-m)
          pmm=pmmp1
          pmmp1=pll
          array(l,m) = pll * pow_xy / factorial_array(l+m)
       enddo
       pow_xy = pow_xy * dcmplx(x,y)
    enddo

    pmm = minus_double_factorial_array(n_max)
    array(n_max,n_max) = pmm * pow_xy / factorial_array(2*n_max) ! When l=m
  end subroutine calculate_regular_harmonics_array

  real(8) function double_factorial(m) result(pmm)
    implicit none
    integer, intent(in) :: m
    integer :: i, fact
    
    fact=1
    pmm = 1d0
    do i=1,m
       pmm = pmm * fact
       fact = fact + 2
    enddo
  end function double_factorial

  real(8) pure elemental function minus_double_factorial(m) result(pmm)
    implicit none
    integer, intent(in) :: m
    integer :: i, fact
    
    fact=1
    pmm = 1d0
    do i=1,m
       pmm = pmm * fact
       fact = fact + 2
    enddo
    pmm = (-1)**m * pmm
  end function minus_double_factorial

  ! calculating array of the following
  !  (-1)^m (2m-1)!!  for m=0,1
  !  (2m-1)!!  for m=max_m
  subroutine calculate_minus_double_factorial_array(m, array)
    implicit none
    integer, intent(in) :: m
    real(8), intent(out) :: array(0:m)
    integer :: i, fact

    fact=1
    array(0) = 1
    do i=1,m-1
       array(i) = - array(i-1) * fact
       fact = fact + 2
    enddo
    !array(m) = double_factorial(m)
    array(m) = (-1)**(m-1) * array(m-1) * fact
  end subroutine calculate_minus_double_factorial_array

  subroutine calculate_inverse_factorial_array(n, array)
    implicit none
    integer, intent(in) :: n
    real(8), intent(out) :: array(0:n)
    real(8) :: array_tmp(0:n)
    integer :: i

    array_tmp(0) = 1
    do i=1,n
       array_tmp(i) = array_tmp(i-1) * i
    enddo
    do i=0,n
       array(i) = 1 / array_tmp(i)
    enddo
  end subroutine calculate_inverse_factorial_array
  
  subroutine prepare_factorial_arrays(n_max, factorial_inverse_array, minus_double_factorial_array)
    implicit none
    integer, intent(in) :: n_max
    real(8), intent(out) :: factorial_inverse_array(0:2*n_max), minus_double_factorial_array(0:n_max)

    call calculate_minus_double_factorial_array(n_max, minus_double_factorial_array)
    call calculate_inverse_factorial_array(2*n_max, factorial_inverse_array)
  end subroutine prepare_factorial_arrays

  subroutine calculate_regular_harmonics_1dim_array(n_max, factorial_inverse_array, minus_double_factorial_array, &
       & x, y, z, array)
    use fmm_l_m_index
    implicit none
    integer, intent(in) :: n_max
    real(8), intent(in) :: x, y, z
    real(8), intent(in) :: factorial_inverse_array(0:2*n_max), minus_double_factorial_array(0:n_max)
    complex(8), intent(out) :: array((n_max+1)*(n_max+2)/2)
    real(8) :: r_square
    real(8) :: fact,pll,pmm,pmmp1
    integer :: l, m, m1
    real(8) :: re, im
    complex(8) :: pow_xy

    r_square = x**2 + y**2 + z**2
    
    pow_xy = dcmplx(1d0, 0)

    do m=0,n_max-1
       pmm = minus_double_factorial_array(m)
       m1 = translate_l_m_to_1dim(m,m)
       array(m1) = pmm * pow_xy * factorial_inverse_array(2*m) ! When l=m
       pmmp1=z*(2*m+1)*pmm
       m1 = translate_l_m_to_1dim(m+1,m)
       array(m1) = pmmp1 * pow_xy * factorial_inverse_array(2*m+1) ! When l=m+1
       do l=m+2, n_max ! for first index of array
          pll=(z*(2*l-1)*pmmp1-(l+m-1)*r_square*pmm)/(l-m)
          pmm=pmmp1
          pmmp1=pll
          m1 = translate_l_m_to_1dim(l,m)
          array(m1) = pll * pow_xy * factorial_inverse_array(l+m)
       enddo
       pow_xy = pow_xy * dcmplx(x,y)
    enddo

    pmm = minus_double_factorial_array(n_max)
    m1 = translate_l_m_to_1dim(n_max,n_max)
    array(m1) = pmm * pow_xy * factorial_inverse_array(2*n_max) ! When l=m
  end subroutine calculate_regular_harmonics_1dim_array
  
  subroutine calculate_regular_harmonics_force_1dim_array(n_max, factorial_inverse_array, minus_double_factorial_array, &
       & x, y, z, array_re, array_im, f_coeff)
    use fmm_l_m_index
    implicit none
    integer, intent(in) :: n_max
    real(8), intent(in) :: x, y, z
    real(8), intent(in) :: factorial_inverse_array(0:2*n_max), minus_double_factorial_array(0:n_max)
    real(8), intent(out) :: array_re((n_max+1)*(n_max+2)/2)
    real(8), intent(out) :: array_im((n_max+1)*(n_max+2)/2)
    complex(8), intent(out) :: f_coeff(3,(n_max+1)*(n_max+2)/2)
    real(8) :: r_square
    real(8) :: fact,pll,pmm,pmmp1
    integer :: l, m, m1, m_store
    real(8) :: re, im
    complex(8) :: pow_xy
    real(8) :: d_plus_re, d_plus_im

    f_coeff(:,:) = dcmplx(0d0, 0d0)
    r_square = x**2 + y**2 + z**2
    
    pow_xy = dcmplx(1d0, 0)

    ! When m = 0
    m = 0
    pmm = minus_double_factorial_array(m)
    m1 = translate_l_m_to_1dim(m,m)
    re = pmm * dreal(pow_xy) * factorial_inverse_array(2*m) ! When l=m
    array_re(m1) = 0.5d0 * re
    array_im(m1) = 0
    m_store = translate_l_m_to_1dim(m+1,m+1) ! m_store is m+1 (i.e. m1=m-1)
    f_coeff(1,m_store) = dcmplx( re, -0d0)
    f_coeff(2,m_store) = dcmplx(-0d0, -re)
    f_coeff(3,m_store-1) = dcmplx(-re, 0d0)
    if (m >= 1) then
       m_store = translate_l_m_to_1dim(m+1,m-1) ! m_store is m-1 (i.e. m1=m+1)
       f_coeff(1,m_store) = f_coeff(1,m_store) + dcmplx(-re,  0d0)
       f_coeff(2,m_store) = f_coeff(2,m_store) + dcmplx(-0d0, -re)
    endif
    
    pmmp1=z*(2*m+1)*pmm
    m1 = translate_l_m_to_1dim(m+1,m)
    re = pmmp1 * dreal(pow_xy) * factorial_inverse_array(2*m+1) ! When l=m+1
    array_re(m1) = 0.5d0 * re
    array_im(m1) = 0d0
    if ((m+2) <= n_max) then
       m_store = translate_l_m_to_1dim(m+2,m+1) ! m_store is m+1 (i.e. m1=m-1)
       f_coeff(1,m_store) = dcmplx( re, -0d0)
       f_coeff(2,m_store) = dcmplx(-0d0, -re)
       f_coeff(3,m_store-1) = dcmplx(-re,  0d0)
       if (m >= 1) then
          m_store = translate_l_m_to_1dim(m+2,m-1) ! m_store is m-1 (i.e. m1=m+1)
          f_coeff(1,m_store) = f_coeff(1,m_store) + dcmplx(-re,  0d0)
          f_coeff(2,m_store) = f_coeff(2,m_store) + dcmplx(-0d0, -re)
       endif
    endif
    do l=m+2, n_max ! for first index of array
       pll=(z*(2*l-1)*pmmp1-(l+m-1)*r_square*pmm)/(l-m)
       pmm=pmmp1
       pmmp1=pll
       m1 = translate_l_m_to_1dim(l,m)
       re = pll * dreal(pow_xy) * factorial_inverse_array(l+m)
       array_re(m1) = 0.5d0 * re
       array_im(m1) = 0d0
       if (l < n_max) then
          !write(*,*) "$l,m=", l+1,m+1
          m_store = translate_l_m_to_1dim(l+1,m+1) ! m_store is m+1 (i.e. m1=m-1)
          f_coeff(1,m_store) = dcmplx( re, -0d0)
          f_coeff(2,m_store) = dcmplx(-0d0, -re)
          f_coeff(3,m_store-1) = dcmplx(-re, 0d0)
          
          if (((l+1) >= (m-1)) .and. (m>=1)) then
             !write(*,*) "@l,m=", l+1,m-1
             m_store = translate_l_m_to_1dim(l+1,m-1) ! m_store is m-1 (i.e. m1=m+1)
             f_coeff(1,m_store) = f_coeff(1,m_store) + dcmplx(-re,  0d0)
             f_coeff(2,m_store) = f_coeff(2,m_store) + dcmplx(-0d0, -re)
          endif
       endif
    enddo
    pow_xy = pow_xy * dcmplx(x,y)
    
    ! When m >= 1
    do m=1,n_max-1
       pmm = minus_double_factorial_array(m)
       m1 = translate_l_m_to_1dim(m,m)
       re = pmm * dreal(pow_xy) * factorial_inverse_array(2*m) ! When l=m
       im = pmm * dimag(pow_xy) * factorial_inverse_array(2*m) ! When l=m
       array_re(m1) = re
       array_im(m1) = im
       m_store = translate_l_m_to_1dim(m+1,m+1) ! m_store is m+1 (i.e. m1=m-1)
       f_coeff(1,m_store) = dcmplx( re, -im)
       f_coeff(2,m_store) = dcmplx(-im, -re)
       f_coeff(3,m_store-1) = dcmplx(-re, im)
       if (m >= 1) then
          m_store = translate_l_m_to_1dim(m+1,m-1) ! m_store is m-1 (i.e. m1=m+1)
          f_coeff(1,m_store) = f_coeff(1,m_store) + dcmplx(-re,  im)
          f_coeff(2,m_store) = f_coeff(2,m_store) + dcmplx(-im, -re)
       endif
       
       pmmp1=z*(2*m+1)*pmm
       m1 = translate_l_m_to_1dim(m+1,m)
       re = pmmp1 * dreal(pow_xy) * factorial_inverse_array(2*m+1) ! When l=m+1
       im = pmmp1 * dimag(pow_xy) * factorial_inverse_array(2*m+1) ! When l=m+1
       array_re(m1) = re
       array_im(m1) = im
       if ((m+2) <= n_max) then
          m_store = translate_l_m_to_1dim(m+2,m+1) ! m_store is m+1 (i.e. m1=m-1)
          f_coeff(1,m_store) = dcmplx( re, -im)
          f_coeff(2,m_store) = dcmplx(-im, -re)
          f_coeff(3,m_store-1) = dcmplx(-re,  im)
          if (m >= 1) then
             m_store = translate_l_m_to_1dim(m+2,m-1) ! m_store is m-1 (i.e. m1=m+1)
             f_coeff(1,m_store) = f_coeff(1,m_store) + dcmplx(-re,  im)
             f_coeff(2,m_store) = f_coeff(2,m_store) + dcmplx(-im, -re)
          endif
       endif
       do l=m+2, n_max ! for first index of array
          pll=(z*(2*l-1)*pmmp1-(l+m-1)*r_square*pmm)/(l-m)
          pmm=pmmp1
          pmmp1=pll
          m1 = translate_l_m_to_1dim(l,m)
          re = pll * dreal(pow_xy) * factorial_inverse_array(l+m)
          im = pll * dimag(pow_xy) * factorial_inverse_array(l+m)
          array_re(m1) = re
          array_im(m1) = im
          if (l < n_max) then
             !write(*,*) "$l,m=", l+1,m+1
             m_store = translate_l_m_to_1dim(l+1,m+1) ! m_store is m+1 (i.e. m1=m-1)
             f_coeff(1,m_store) = dcmplx( re, -im)
             f_coeff(2,m_store) = dcmplx(-im, -re)
             f_coeff(3,m_store-1) = dcmplx(-re, im)

             if (((l+1) >= (m-1)) .and. (m>=1)) then
                !write(*,*) "@l,m=", l+1,m-1
                m_store = translate_l_m_to_1dim(l+1,m-1) ! m_store is m-1 (i.e. m1=m+1)
                f_coeff(1,m_store) = f_coeff(1,m_store) + dcmplx(-re,  im)
                f_coeff(2,m_store) = f_coeff(2,m_store) + dcmplx(-im, -re)
             endif
          endif
       enddo
       pow_xy = pow_xy * dcmplx(x,y)
    enddo

    f_coeff(3,2) = dcmplx(-array_re(1), 0d0)  ! coeff for f_z when l=0, m=1
    do l=1, n_max-1
       m1 = translate_l_m_to_1dim(l,1)
       m_store = translate_l_m_to_1dim(l+1,0)
       f_coeff(1,m_store) = dcmplx(-array_re(m1), 0d0)
       f_coeff(2,m_store) = dcmplx(-array_im(m1), 0d0)
       f_coeff(3,m_store) = dcmplx(-array_re(m1-1), 0d0)
    enddo

    pmm = minus_double_factorial_array(n_max)
    m1 = translate_l_m_to_1dim(n_max,n_max)
    array_re(m1) = pmm * dreal(pow_xy) * factorial_inverse_array(2*n_max) ! When l=m
    array_im(m1) = pmm * dimag(pow_xy) * factorial_inverse_array(2*n_max) ! When l=m
  end subroutine calculate_regular_harmonics_force_1dim_array

end module regular_solid_harmonics_cartesian

