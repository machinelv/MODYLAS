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
!! \brief  Module and subroutines to store original atom coordinates 
!!         and velocities.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to store original atom coordinates and velocities.
!! \author Kensuke Iwahashi
!<
module trajectory_org
  use mpi_tool
  implicit none
  real(8),allocatable :: xyz(:,:)
  real(8),allocatable :: v(:,:)
  integer(4) :: n=0,nmols=0

contains

  subroutine fmod_n(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    n = ivalue
  end subroutine fmod_n

  subroutine fmod_alloc_xyz(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(xyz(3,ivalue))
  end subroutine fmod_alloc_xyz

  subroutine fmod_set_xyz(i0,value1,value2,value3,line_number)
    implicit none
    integer(4), intent(in) :: i0, line_number
    real(8), intent(in) :: value1,value2,value3
    integer(4) :: i
    i = i0 + 1
    if (i > n) then
       write(0,'(a,i0)')  'ERROR (line ', line_number, ') : '
       write(0,'(a,i0,a)')  &
            &    'The number of md_gegeric.positions is out of bounds.  '// &
            &    'It must be less than ', n, '.'
       call modylas_abort()
    endif
    xyz(1,i) = value1
    xyz(2,i) = value2
    xyz(3,i) = value3
  end subroutine fmod_set_xyz
!----------------------------------------------------------------------
  subroutine fmod_alloc_v(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(v(3,ivalue))
  end subroutine fmod_alloc_v
!----------------------------------------------------------------------
  subroutine fmod_set_v(i0, value_x, value_y, value_z, line_number)
    implicit none
    integer(4), intent(in) :: i0, line_number
    real(8), intent(in) :: value_x, value_y, value_z
    integer(4) :: i
    i = i0 + 1
    if (i > n) then
       write(0,'(a,i0)')  'ERROR:'! line ', line_number
       write(0,'(a,i0,a)') &
            &    'The number of md_gegeric.velocities is out of bounds.  '// &
            &    'It must be less than ', n, '.'
       call modylas_abort()
    endif
    v(1,i) = value_x
    v(2,i) = value_y
    v(3,i) = value_z
  end subroutine fmod_set_v

end module trajectory_org
