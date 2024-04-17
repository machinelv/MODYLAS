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
!! \brief  Module for rotation matrices to be used in M2L transformation of FMM.
!<
!----------------------------------------------------------------------
!>
!! \brief  Calculating direct (inverse) rotation matrix with respect to the basis of singular (regular) solid harmonics in 1-dimensional index, respectively.
!! \author Tatsuya Sakashita, Noriyuki Yoshii.
!<
module regular_singular_rotation_1dim
  use precision
  implicit none
  
contains

  subroutine generate_d_matrix(beta, lmax, d_mat)
    USE ieee_arithmetic
    use math_functions, only : factorial
    implicit none
    integer, intent(in) :: lmax
    real(8), intent(in) :: beta
    real(kind=wmp), intent(out) :: d_mat(-lmax:lmax,-lmax:lmax,0:lmax)
    real(8) :: cos_b, cos_b_half, cos_b_half_square, tan_b_half
    integer :: l, m, n
    real(8), parameter :: eps=1d-5

    cos_b = cos(beta)
    cos_b_half = cos(beta/2)
    cos_b_half_square = cos(beta/2)**2
    tan_b_half = tan(beta/2)
    
    d_mat(0,0,0) = 1
    d_mat(0,0,1) = cos_b
    d_mat(1,0,1) = - sin(beta)/sqrt(2d0)
    d_mat(1,-1,1) = (1-cos_b) / 2;    d_mat(1,1,1) = (1+cos_b) / 2
    
    do l=2,lmax
       do m=0,l-2
          do n=-m,m
             if (l==3 .and. m==2 .and. n==-2) then
                write(*,*) "icchi"
                stop
             endif
             d_mat(m,n,l) = dble(l*(2*l-1)) / sqrt(dble(l**2-m**2)*(l**2-n**2)) * &
                  & ( (cos_b-dble(m*n)/(l*(l-1))) * d_mat(m,n,l-1) &
                  &   - sqrt(dble((l-1)**2-m**2)*((l-1)**2-n**2)) / (l-1) / (2*l-1) * d_mat(m,n,l-2) )
             if (.not.ieee_is_finite(d_mat(m,n,l))) then
                write(0,*) "notfinite=", "m=", m, " n=", n, "l=", l, &
                     & "d=", d_mat(m,n,l), d_mat(m,n,l-1), d_mat(m,n,l-2)
                stop
             endif
          enddo
       enddo
       
       if (abs(cos_b_half)>=eps) then
          d_mat(l,l,l) = cos_b_half_square**l
          if (.not.ieee_is_finite(d_mat(l,l,l))) then
             write(0,*) "nan00=", d_mat(l,l,l)
             stop
          endif
          d_mat(l-1,l-1,l) = (l*cos_b-l+1) * d_mat(l-1,l-1,l-1)
          if (.not.ieee_is_finite(d_mat(l-1,l-1,l))) then
             write(0,*) "nan0=", d_mat(l-1,l-1,l)
             stop
          endif
          do n=l,-l+1,-1
             d_mat(l,n-1,l) = - sqrt(dble(l+n)/(l-n+1)) * tan_b_half * d_mat(l,n,l)
             if (.not.ieee_is_finite(d_mat(l,n-1,l))) then
                write(0,*) "nan1=", d_mat(l,n-1,l)
                stop
             endif
          enddo
          do n=l-1,-l+2,-1
             if(abs(l*cos_b - n) <= eps) then
                d_mat(l-1,n-1,l) = sqrt(dble(factorial(2*l-1))/(factorial(l-n+1)*factorial(l+n-1))) &
                     & * cos_b_half**(l+n-2) * (-sin(beta/2))**(l-n)
             else
                d_mat(l-1,n-1,l) = - sqrt(dble(l+n)/(l-n+1)) * (l*cos_b-n+1) / (l*cos_b-n) * tan_b_half * d_mat(l-1,n,l)
                if (.not.ieee_is_finite(d_mat(l-1,n-1,l))) then
                   d_mat(l-1,n-1,l) = - sqrt(dble(l+n)/(l-n+1)) * (l*cos_b-n+1) / (l*cos_b-n) * tan_b_half * d_mat(l-1,n,l)
                endif
             endif
          enddo
       else
          d_mat(l,:,l) = 0d0
          d_mat(l-1,:,l) = 0d0
          d_mat(l,-l,l) = 1d0;  d_mat(l-1,1-l,l) = -1d0
       endif
    enddo

    ! copying by using mirror symmetry
    do l=1,lmax
       do m=0,l
          do n=-m,m
             d_mat(n,m,l) = (-1)**(m+n) * d_mat(m,n,l)
             d_mat(-m,-n,l) = (-1)**(m+n) * d_mat(m,n,l)
             d_mat(-n,-m,l) = d_mat(m,n,l)
          enddo
       enddo
    enddo

  end subroutine generate_d_matrix

  integer function factorial_cancel(l,n)
    use math_functions, only : factorial
    implicit none
    integer, intent(in) :: l, n
    integer :: i
    factorial_cancel = 1
    do i=l+n, 2*l-1
       factorial_cancel = factorial_cancel * i
    enddo
    factorial_cancel = factorial_cancel 
  end function factorial_cancel

  real(8) function factorial_denom(l)
    use math_functions, only : factorial
    implicit none
    integer, intent(in) :: l
    integer :: i
    factorial_denom = 1.d0
    do i=l+2, 2*l-1
       factorial_denom = factorial_denom * i
    enddo
    factorial_denom = factorial_denom / factorial(l-1)
  end function factorial_denom
  
  integer(4) pure elemental function translate_rot_index_to_1dim(l, p, q)
    implicit none
    integer, intent(in) :: l, p, q
    
    translate_rot_index_to_1dim = l*(l+1)*(2*l+1)/6 + (l+1)*q + p + 1

  end function translate_rot_index_to_1dim
  
  ! For M2L  direct:singular,  inverse:regular
  subroutine generate_rotation_matrix_1dim(beta, alpha, lmax, &
       & sh_rot1, sh_rot2, &
       & sh_inv_rot1, sh_inv_rot2)
    use math_functions, only : factorial
    implicit none
    real(8), intent(in) :: beta, alpha
    integer, intent(in) :: lmax
    real(kind=wmp), intent(out) :: sh_rot1((lmax+1)*(lmax+2)*(2*lmax+3)/6,2), &
   &                               sh_rot2((lmax+1)*(lmax+2)*(2*lmax+3)/6,2)
    real(kind=wmp), intent(out) :: sh_inv_rot1((lmax+1)*(lmax+2)*(2*lmax+3)/6,2), &
   &                               sh_inv_rot2((lmax+1)*(lmax+2)*(2*lmax+3)/6,2)
    real(8) :: angle, coeff
    complex(8) :: sh_d(0:lmax,0:lmax,0:lmax)
    integer :: l, m, n, ind_1dim

    call generate_sh_d_matrix(-beta, lmax, sh_d)  ! -beta is direct rotation

    ! direct: singular type coeff
    do l=0,lmax
       do m=0,l
          do n=0,l
             ind_1dim = translate_rot_index_to_1dim(l, m, n)
             angle = n*alpha   ! Assume exp(-angle). So angle = n*alpha
             coeff = sqrt( dble(factorial(l-n) * factorial(l+n)) / (factorial(l-m) * factorial(l+m)) )
             if (m /= 0) then
                sh_rot1(ind_1dim,1) = coeff * &
             &      (-1)**m *   cos(angle)*dreal(sh_d(m,n,l))
                sh_rot1(ind_1dim,2) = coeff * &
             &      (-1)**m * (+sin(angle)*dimag(sh_d(m,n,l)))
                sh_rot2(ind_1dim,1) = coeff * &
             &      (-1)**m *   sin(angle)*dreal(sh_d(m,n,l))
                sh_rot2(ind_1dim,2) = coeff * &
             &      (-1)**m * (-cos(angle)*dimag(sh_d(m,n,l)))
             else
                sh_rot1(ind_1dim,1) = coeff * (-1)**m * cos(angle)*dreal(sh_d(m,n,l))
                sh_rot1(ind_1dim,2) = 0d0
                sh_rot2(ind_1dim,1) = coeff * (-1)**m * sin(angle)*dreal(sh_d(m,n,l))
                sh_rot2(ind_1dim,2) = 0d0 
             endif
          enddo
       enddo
    enddo

    call generate_sh_d_matrix(beta, lmax, sh_d)  ! +beta is inverse rotation

    ! inverse: regular type coeff
    do l=0,lmax
       do m=0,l
          angle = -m*alpha
          do n=0,l
             ind_1dim = translate_rot_index_to_1dim(l, m, n)
             coeff = sqrt( dble(factorial(l-m) * factorial(l+m)) / (factorial(l-n) * factorial(l+n)) )
             if (m /= 0) then
                sh_inv_rot1(ind_1dim,1) = coeff * &
                     & dreal(sh_d(m,n,l)) * cos(angle)
                sh_inv_rot1(ind_1dim,2) = coeff * &
                     & dreal(sh_d(m,n,l)) * sin(angle)
                sh_inv_rot2(ind_1dim,1) = coeff * dreal( &
                     & dimag(sh_d(m,n,l)) * dcmplx(-sin(angle), cos(angle)))
                sh_inv_rot2(ind_1dim,2) = coeff * dimag( &
                     & dimag(sh_d(m,n,l)) * dcmplx(-sin(angle), cos(angle)))
             else
                sh_inv_rot1(ind_1dim,1) = coeff * dreal(sh_d(m,n,l)) * cos(angle)
                sh_inv_rot1(ind_1dim,2) = 0d0
                sh_inv_rot2(ind_1dim,1) = coeff * dreal(dimag(sh_d(m,n,l)) * dcmplx(-sin(angle), 0d0))
                sh_inv_rot2(ind_1dim,2) = coeff * dimag(dimag(sh_d(m,n,l)) * dcmplx(-sin(angle), 0d0))
             endif
          enddo
       enddo
    enddo
  end subroutine generate_rotation_matrix_1dim

  subroutine generate_sh_d_matrix(beta, lmax, sh_d)
    USE ieee_arithmetic
    implicit none
    integer, intent(in) :: lmax
    real(8), intent(in) :: beta
    complex(8), intent(out) :: sh_d(0:lmax,0:lmax,0:lmax)
    real(kind=wmp) :: d_mat(-lmax:lmax,-lmax:lmax,0:lmax)
    integer :: l, m, n

    call generate_d_matrix(beta, lmax, d_mat)

    do l=0,lmax
       do m=0,l
          sh_d(m,0,l) = dcmplx( d_mat(m,0,l), d_mat(m,0,l) )
          do n=1,l
             sh_d(m,n,l) = dcmplx( d_mat(m,n,l) + (-1)**n * d_mat(m,-n,l), d_mat(m,n,l) - (-1)**n * d_mat(m,-n,l) )
          enddo
       enddo
    enddo
  end subroutine generate_sh_d_matrix

end module regular_singular_rotation_1dim
