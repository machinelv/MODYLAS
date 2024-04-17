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
!! \brief  Module and subroutines to store segment information.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to set and store segment information.
!! \author Kazushi FUJIMOTO, Kensuke Iwahashi
!<
module segments
  implicit none
  integer(4) :: nsegments=-1
  integer(4),allocatable :: segtop(:), seg_natoms(:)
  real(8),allocatable :: seg_cx(:)
  real(8),allocatable :: seg_cy(:)
  real(8),allocatable :: seg_cz(:)

  integer(4),allocatable :: irearrange(:)
  logical :: seg_ready = .FALSE.

contains

  subroutine fmod_alloc_segtop(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(segtop(ivalue))
  end subroutine fmod_alloc_segtop
!----------------------------------------------------------------------
  subroutine fmod_alloc_seg_natoms(ivalue)
    implicit none
    integer(4) :: ivalue
    allocate(seg_natoms(ivalue))
    seg_natoms = -1
  end subroutine fmod_alloc_seg_natoms
!----------------------------------------------------------------------
  subroutine fmod_alloc_seg_center_x(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(seg_cx(ivalue))
  end subroutine fmod_alloc_seg_center_x
!----------------------------------------------------------------------
  subroutine fmod_alloc_seg_center_y(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(seg_cy(ivalue))
  end subroutine fmod_alloc_seg_center_y
!----------------------------------------------------------------------
  subroutine fmod_alloc_seg_center_z(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    allocate(seg_cz(ivalue))
  end subroutine fmod_alloc_seg_center_z
!----------------------------------------------------------------------
  subroutine fmod_howmany_segments(ivalue)
    implicit none
    integer(4), intent(in) :: ivalue
    nsegments = ivalue
  end subroutine fmod_howmany_segments
!----------------------------------------------------------------------
  subroutine fmod_set_seg_natoms(iseg0, ivalue, line_number)
    use mpi_tool
    implicit none
    integer(4), intent(in) :: iseg0, ivalue, line_number
    integer(4) :: iseg
    iseg = iseg0 + 1
    if (iseg > nsegments) then
       write(0,'(a,i0)')  'ERROR (line ', line_number, ') : '
       write(0,'(a,i0,a)')  &
            &    'The number of md_periodic.segments[].howmany_atoms ' // &
            &    'is out of bounds.  '// &
            &    'It must be less than ', nsegments, '.'
       call modylas_abort()
    endif
    seg_natoms(iseg) = ivalue
  end subroutine fmod_set_seg_natoms
!----------------------------------------------------------------------
  subroutine init_segment_check
    !     check atom order:
    !           atom numbers must be in order along segment id number
    use trajectory_org, only : n
    use mpi_tool
    implicit none
    integer(4) :: im,j,k
    integer(4) :: isegstart, isegend
    include 'mpif.h'

    if( nsegments.le.0 ) then
       if(myrank==0) write(0,'("ERROR: nsegment <= 0",i8)') nsegments
       call modylas_abort()
    endif

    k = segtop(1)
    if( k.ne.0 ) then
       if(myrank==0) write(0,9000) k
       call modylas_abort()
    endif

    do im= 1,nsegments

       isegstart = segtop(im)+1
       isegend   = isegstart + seg_natoms(im)-1

       do j= isegstart, isegend
          if(k+1 .ne. j) then
             write(0,9000) k
             !              stop
          endif
          k=k+1
       enddo

    enddo

    if( k.ne.n ) then
       if(myrank==0) write(0,9000) k-1
       call modylas_abort()
    endif

9000 format('ERROR: the order of atoms is wrong:',i10)
  end subroutine init_segment_check
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to read segments.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_segment
    use device_numbers, only : f_mdff
    use mpi_tool
    implicit none
    include 'mpif.h'
    integer(4) :: ierr

    if(myrank==0) read(f_mdff) nsegments
    call MPI_Bcast(nsegments, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
    call fmod_alloc_seg_natoms(nsegments)
    call fmod_alloc_seg_center_x(nsegments)
    call fmod_alloc_seg_center_y(nsegments)
    call fmod_alloc_seg_center_z(nsegments)
    call fmod_alloc_segtop(nsegments)
    if(myrank==0) then
       read(f_mdff) seg_natoms
       read(f_mdff) segtop
    endif
    call MPI_Bcast(seg_natoms, nsegments, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(segtop, nsegments, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
  end subroutine read_segment
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write segments to mdff.bin.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_segment
    use device_numbers, only : f_mdff
    implicit none
    write(f_mdff) nsegments
    write(f_mdff) seg_natoms
    write(f_mdff) segtop
  end subroutine write_segment
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write segments.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_segment
    implicit none
    write(*,*) '[write_segment]'
    write(*,*) nsegments
    write(*,*) seg_natoms
    write(*,*) segtop
  end subroutine write_memory_segment

  subroutine segment_is_not_ready
    use mpi_tool
    implicit none
    write(0,*) 'ERROR: Segment must be read before.'
    call modylas_abort()
  end subroutine segment_is_not_ready

  subroutine fmod_set_seg2atom(iseg0,j0,ivalue0,line_number)
    use trajectory_org, only : n
    use mpi_tool
    implicit none
    integer(4), intent(in) :: iseg0, j0, ivalue0, line_number
    integer(4) :: iseg, j, ivalue
    integer(4),save :: ic=0

    iseg = iseg0 + 1
    ivalue = ivalue0 + 1
    j = j0 + 1

    if (j > seg_natoms(iseg)) then
       write(0,'(a,i0)')  'ERROR (line ', line_number, ') : '
       write(0,'(a,i0,a)')  &
            &    'The number of md_periodic.segments.[1].atoms[] is ' // &
            &    'out of bounds.  '// &
            &    'It must be less than ', seg_natoms(iseg), '.'
       call modylas_abort()
    endif

    if(j == 1 ) then
       segtop(iseg) = ivalue0
    endif
    ic = ivalue

    if (ic == n) then
       seg_ready = .TRUE.
    endif
  end subroutine fmod_set_seg2atom

  subroutine fmod_alloc_irearrange(ivalue)
    use trajectory_org, only : n
    implicit none
    integer(4), intent(in) :: ivalue
    integer(4) :: i
    allocate(irearrange(ivalue))
    do i = 1, n
       irearrange(i) = i
    enddo
  end subroutine fmod_alloc_irearrange

end module segments
