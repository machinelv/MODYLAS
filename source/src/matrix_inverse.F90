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
!! \brief Module and subroutine which calculate the inverse matrix
!<
!----------------------------------------------------------------------
!>
!! \brief Module to calculate the inverse matrix
!! \author Noriyuki Yoshii
!<
module matrix_inverse_mod
  implicit none

contains

  subroutine matrix_inverse (m,m_inv)
!----------------------------------------------------------------------
    implicit none
    real(8), intent(in) :: m(3,3)
    real(8), intent(out) :: m_inv(3,3)
    real(8) :: x,rx

    x = m(1,1)*m(2,2)*m(3,3) + m(1,2)*m(2,3)*m(3,1) &
    &  + m(1,3)*m(2,1)*m(3,2) - m(1,3)*m(2,2)*m(3,1) &
    &  - m(1,1)*m(2,3)*m(3,2) - m(1,2)*m(2,1)*m(3,3)

    rx=1.d0/x

    m_inv(1,1)=(m(2,2)*m(3,3)-m(2,3)*m(3,2))*rx
    m_inv(1,2)=(m(3,2)*m(1,3)-m(3,3)*m(1,2))*rx
    m_inv(1,3)=(m(1,2)*m(2,3)-m(1,3)*m(2,2))*rx
    m_inv(2,1)=(m(2,3)*m(3,1)-m(2,1)*m(3,3))*rx
    m_inv(2,2)=(m(3,3)*m(1,1)-m(3,1)*m(1,3))*rx
    m_inv(2,3)=(m(1,3)*m(2,1)-m(1,1)*m(2,3))*rx
    m_inv(3,1)=(m(2,1)*m(3,2)-m(2,2)*m(3,1))*rx
    m_inv(3,2)=(m(3,1)*m(1,2)-m(3,2)*m(1,1))*rx
    m_inv(3,3)=(m(1,1)*m(2,2)-m(1,2)*m(2,1))*rx
  end subroutine matrix_inverse

end module matrix_inverse_mod
