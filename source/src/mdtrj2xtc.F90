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
!! \brief  Program to convert file from mdtrj to xtc.
!!         THis only works with Single CPU. When you use more than 1cpu,
!!         this waits for at MPI barrier eternally.
!! \author Ryo Urano
!<
program mdtrj2xtc
  use trajectory_org
  use atom_mass
  use thermostat
  use md_const
  use device_numbers
  use session_name_mod
  use param
  use mpi_tool
  use unit_cell
  use commandline_args
  use input_mod
  use file_xtc
  implicit none
  integer(4),parameter :: f_xyz=97, f_cel=98
  integer(4) :: mdstep, natom, io, i, n_nhc, ndisplay, nintvl, k
  integer(4),allocatable :: atomtype(:)
  real(8) :: dt,massi
  character(2) :: str_atom

  call mpistart()

#ifdef XTC

  if (myrank==0) then
     write(6,*) 'Please specify the interval of the output from mdtrj to xtc?'
     read(5,*) nintvl

     ndisplay = n

     call parse_args
     call parse_input(.false.)

     xtc_output=.true.
     call open_xtc
     open(f_trj, file=trim(session_name)// '.mdtrj', &
          &         iostat=io, status='old', &
          &         access='sequential',form='unformatted')

     k = 0
     do while (.true.)
        k = k + 1
        read(f_trj, end=100) mdstep, dt
        read(f_trj) natom
        read(f_trj) xyz(1:3,:), v(1:3,:)
        read(f_trj) n_nhc
        read(f_trj) rss, vss
        read(f_trj) n_nhc
        read(f_trj) rssb,vssb
        read(f_trj) cellx,celly,cellz,alpha,beta,gamma,vboxg

        if (mod(k,nintvl) .ne. 0)  cycle
        write(*,*) '#mdstep ', mdstep


    xtcf % time= mdstep*dt *1e+12
    xtcf % STEP= mdstep-xtc_start
    xtcf % box(1,1)=cellx*1e+9
    xtcf % box(2,2)=celly*1e+9
    xtcf % box(3,3)=cellz*1e+9
    xtcf%pos = real(xyz)*1e+9 !! nm 
    xtcf % prec=1000.000000
    call xtc_out % write(n, xtcf % STEP, xtcf % time, xtcf % box, xtcf%pos, xtcf % prec)

     enddo
100  continue
     close(f_trj)
     call xtc_out % close

 endif
#else

 write(*,*) " ########################################################"
 write(*,*) " This program depends on the externaly library: xdrfile"
 write(*,*) " Without the library, nothing happens in this program"
 write(*,*) " ########################################################"
 
#endif

  call mpiend()

end program mdtrj2xtc
