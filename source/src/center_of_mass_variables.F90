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
!! \brief Module to store information of distance constrain
!!        between two mass centers.
!<
!----------------------------------------------------------------------
!>
!! \brief Module to store information of distance constraints.
!! \author Yoshimichi Andoh
!<
module center_of_mass_variables
  implicit none
  integer(4) :: groupAtop, groupAend, groupBtop, groupBend
  real(8) :: KratCOM = 0.0d0,KratCOM2
  real(8) :: d_COMR=0d0  !! initial value
  logical :: constrain_COM=.false.
  logical :: change_distCOM=.false.

contains
!---------------------------------------------------------------------
!>
!! \brief Subroutine to assign top and end atom numbers in groupA.
!! \author Yoshimichi Andoh
!<
  subroutine fmod_set_groupA(ivalue1, ivalue2)
    implicit none
    integer(4), intent(in) :: ivalue1, ivalue2
    groupAtop = ivalue1
    groupAend = ivalue2
  end subroutine fmod_set_groupA
!---------------------------------------------------------------------
!>
!! \brief Subroutine to assign top and end atom numbers in groupB.
!! \author Yoshimichi Andoh
!<
  subroutine fmod_set_groupB(ivalue1, ivalue2)
    implicit none
    integer(4), intent(in) :: ivalue1, ivalue2
    groupBtop = ivalue1
    groupBend = ivalue2
  end subroutine fmod_set_groupB
!---------------------------------------------------------------------
!>
!! \brief Subroutine to set constrain distance between groupA and B.
!! \author Yoshimichi Andoh
!<
  subroutine fmod_set_dist_COM(value1)
    implicit none
    real(8), intent(inout) :: value1
    ! convert [A] ==> [m]
    value1 = value1 * 1.0d-10
    KratCOM = value1
    KratCOM2 = value1*value1
  end subroutine fmod_set_dist_COM
!---------------------------------------------------------------------
!>
!! \brief Subroutine to set increment of constrain distance.
!! \author Yoshimichi Andoh
!<
  subroutine fmod_set_d_COMR(value1)
    implicit none
    real(8) :: value1
    ! convert [A] ==> [m]
    value1 = value1 * 1.0d-10
    d_COMR = value1
  end subroutine fmod_set_d_COMR
!---------------------------------------------------------------------
!>
!! \brief Subroutine to read constrain_COM in .mddef.
!! \author Yoshimichi Andoh
!<
  subroutine fmod_constrain_COM(ivalue)
    use mpi_tool
    implicit none
    integer(4), intent(in) :: ivalue
    if (ivalue /= 0) then
       constrain_COM = .true.
       if(myrank==0)write(*,*) 'constrain_COM = yes'
    else
       constrain_COM = .false.  !! default
       if(myrank==0)write(*,*) 'constrain_COM = no, or not set yes/no'
    endif
  end subroutine fmod_constrain_COM
!---------------------------------------------------------------------
!>
!! \brief Subroutine to read change_distCOM in .mddef.
!! \author Yoshimichi Andoh
!<
  subroutine fmod_change_distCOM(ivalue)
    use mpi_tool
    implicit none
    integer(4), intent(in) :: ivalue
    if (ivalue /= 0) then
       change_distCOM = .true.
       if(myrank==0)write(*,*) 'change_distCOM = yes'
    else
       change_distCOM = .false.  !! default
       if(myrank==0)write(*,*) 'change_distCOM = no, or not set yes/no'
    endif
  end subroutine fmod_change_distCOM

end module center_of_mass_variables
