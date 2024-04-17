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
!! \brief  Module and subroutines to control .mdxyz.bin file I/O.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to control .mdxyz.bin file I/O.
!! \author Kensuke Iwahashi, Kazushi Fujimoto
!<
module file_mdxyz_bin
  implicit none

contains

!>
!! \brief  Subroutine to read .mdxyz.bin.
!! \author Kensuke Iwahashi, Kazushi Fujimoto
!<
  subroutine read_mdxyzbin
    use trajectory_org
    use thermostat
    use file_restart
    use device_numbers
    use session_name_mod
    use mpi_tool
    use unit_cell
    use subcell, only : alloc_i2m
    use file_utility, only : open_binary_file
    implicit none
    integer(4) :: io, n_nhc
    !KF
    integer(4) :: ierr
    include 'mpif.h'
!KF end
    if(myrank==0) then
       call open_binary_file(f_restart_bin, trim(session_name)// '.mdxyz.bin')
       read(f_restart_bin,end=100) n
    endif
    call MPI_Bcast(n, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
    if (.not. allocated(xyz)) then
       call alloc_i2m(n)
       call fmod_alloc_xyz(n)
       call fmod_alloc_v(n)
    endif
    if(myrank.eq.0) then
       !       Read coordinates and velocities of atoms.
       !       Trajectory must be arranged in the same order as input.
       read(f_restart_bin) xyz(1:3,:), v(1:3,:)
       !       Read positions and velocities of thermostats.
       read(f_restart_bin) n_nhc
       if (n_nhc .ne. nnhc) then
          write(0,*) 'ERROR: nnhc'
          call modylas_abort()
       endif
       read(f_restart_bin) rss, vss
       !       Read positions and velocities of barostats.
       read(f_restart_bin) n_nhc
       if (n_nhc .ne. nnhc) then
          write(0,*) 'ERROR: nnhc'
          call modylas_abort()
       endif
       read(f_restart_bin) rssb, vssb
       !       Read cell parameters (length and angles).
       read(f_restart_bin) cellx,celly,cellz, alpha,beta,gamma,vboxg
       cellxh = 0.5d0 * cellx
       cellyh = 0.5d0 * celly
       cellzh = 0.5d0 * cellz
       if (abs(alpha-90.0) .gt. 1.0d-5)  cuboid = .false.
       if (abs(beta -90.0) .gt. 1.0d-5)  cuboid = .false.
       if (abs(gamma-90.0) .gt. 1.0d-5)  cuboid = .false.
    endif ! myrank==0
100 continue
    close(f_restart_bin)
    !KF
    call MPI_Bcast(xyz  ,  3*n, MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(v    ,  3*n, MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(n_nhc,  1  , MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(rss  ,  5  , MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(vss  ,  5  , MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(rssb ,  5  , MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(vssb ,  5  , MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(cellx,  1  , MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(celly,  1  , MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(cellz,  1  , MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(alpha,  1  , MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(beta ,  1  , MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(gamma,  1  , MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(vboxg,  9  , MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(cellxh, 1  , MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(cellyh, 1  , MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(cellzh, 1  , MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(cuboid, 1  ,MPI_LOGICAL  , 0, MPI_COMM_WORLD, ierr)
    !KF end
  end subroutine read_mdxyzbin
!--------------------------------------------------------------------------
!>
!! \brief  Subroutine to write mdxyz in binary format.
!! \author Kensuke Iwahashi
!<
  subroutine write_mdxyzbin
    use trajectory_org
    use file_restart
    use device_numbers
    use session_name_mod
    use mpi_tool
    use thermostat
    use unit_cell
    use file_utility, only : create_new_binary_file
    implicit none
    if (myrank==0) then
       call create_new_binary_file(f_restart_bin, trim(session_name)//'.mdxyz.bin')
       write(f_restart_bin) n
       write(f_restart_bin) xyz(1:3,:), v(1:3,:)
       write(f_restart_bin) mnhc
       write(f_restart_bin) rss, vss
       write(f_restart_bin) mnhc
       write(f_restart_bin) rssb, vssb
       write(f_restart_bin) cellx,celly,cellz, alpha,beta,gamma,vboxg
       close(f_restart_bin)
       write(*,*) "mdxyz.bin was succesfully created!"
    endif
  end subroutine write_mdxyzbin

end module file_mdxyz_bin
