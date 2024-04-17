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
!! \brief Module to create index table to be used in M2M transformation of FMM.
!<
!----------------------------------------------------------------------
!>
!! \brief Creating index table to be used in M2M transformation of FMM.
!! \author Noriyuki Yoshii, Tatsuya Sakashita
!<
module fmm_l_m_index
  implicit none
  integer(4) :: len_m2m

contains

!>
!! \brief Function to translate pair of indices (l,m) into 1-dimensional index.
!! \author Noriyuki Yoshii, Tatsuya Sakashita
!<
  integer(4) pure elemental function translate_l_m_to_1dim(l, m)
    implicit none
    integer, intent(in) :: l, m
    
    translate_l_m_to_1dim = l*(l+1)/2 + 1 + m
  end function translate_l_m_to_1dim

!>
!! \brief Subroutine to create table to translate into 1-dimensional index for M2M.
!! \author Noriyuki Yoshii, Tatsuya Sakashita
!<
  subroutine init_m2m_indices(nmax, ind_m2m)
!---------------------------------------------------------------------
    implicit none
    integer(4), intent(in) :: nmax
    integer(4), intent(out), allocatable :: ind_m2m(:,:)
    integer(4) :: j,k,m,l, m1,m2, i1dim

    len_m2m = 0
    do j=0,nmax
       do k=0,j
          do l=0,j
             do m=max(0,-(j-l)+k), min(l,(j-l)+k)
                len_m2m = len_m2m + 1
             enddo
          enddo
       enddo
    enddo
    allocate(ind_m2m(2,len_m2m))

    i1dim = 0
    do j=0,nmax
       do k=0,j
          do l=0,j
             do m=max(0,-(j-l)+k), min(l,(j-l)+k)
                i1dim = i1dim + 1
                ind_m2m(1,i1dim) = translate_l_m_to_1dim(j, k)
                ind_m2m(2,i1dim) = translate_l_m_to_1dim(l, m)
             enddo
          enddo
       enddo
    enddo
  end subroutine init_m2m_indices

end module fmm_l_m_index
