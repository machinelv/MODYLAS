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
!! \brief  Program to convert file from mdtrj to xyz format.
!! \author Kensuke Iwahashi
!<
!----------------------------------------------------------------------
!>
!! \brief  Program to convert file from mdtrj to xyz format.
!! \author Kensuke Iwahashi
!<
program mdtrj2xyz
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
  implicit none
  integer(4),parameter :: f_xyz=97, f_cel=98
  integer(4) :: mdstep, natom, io, i, n_nhc, ndisplay, nintvl, k
  integer(4),allocatable :: atomtype(:)
  real(8) :: dt,massi
  real(8),allocatable :: vdum(:,:)
  character(2) :: str_atom

  call mpistart()
  if (myrank==0) then
     write(6,'(a,i0)') 'myrank in C=', myrank
     write(6,*) 'How many atoms to display?'
     write(6,*) '(if negative value is given, the all atoms is used)'
     read(5,*) ndisplay
     write(6,*) 'How many interval steps to display?'
     read(5,*) nintvl
     call parse_args
     call parse_input(.false.)
     if(ndisplay .le. 0) ndisplay = n
     !       - - - - - - - - - - - - - - - - - - - - -
     allocate(atomtype(n))
     allocate(vdum(3,n))
     !       - - - - - - - - - - - - - - - - - - - - -
     mass = mass / md_ATOMIC_MASS_UNIT
     do i=1,n
        massi = mass(paranum(i))
        if (massi .gt. 1.0e30) then
           massi = massi * 1.0e-30
        endif
        if (dabs(massi-1.00800d0) .lt. 1.0d-2) then
           atomtype(i) = 1
        else if (dabs(massi-12.01100d0) .lt. 1.0d-2) then
           atomtype(i) = 2
        else if (dabs(massi-14.00700d0) .lt. 1.0d-2) then
           atomtype(i) = 3
        else if (dabs(massi-15.99900d0) .lt. 1.0d-2) then
           atomtype(i) = 4
        else if (dabs(massi-30.974000d0) .lt. 1.0d-2) then
           atomtype(i) = 5
        else if (dabs(massi-32.06000d0) .lt. 1.0d-2) then
           atomtype(i) = 6
        else if (dabs(massi-4.00260d0) .lt. 1.0d-2) then
           atomtype(i) = 7
        else if (dabs(massi-20.17970d0) .lt. 1.0d-2) then
           atomtype(i) = 8
        else if (dabs(massi-40.08000d0) .lt. 1.0d-2) then
           atomtype(i) = 9
        else if (dabs(massi-65.37000d0) .lt. 1.0d-2) then
           atomtype(i) = 10
        else if (dabs(massi-55.84700d0) .lt. 1.0d-2) then
           atomtype(i) = 11
        else if (dabs(massi-22.98977d0) .lt. 1.0d-2) then
           atomtype(i) = 12
        else if (dabs(massi-35.45000d0) .lt. 1.0d-2) then
           atomtype(i) = 13
        else if (dabs(massi-39.102000d0) .lt. 1.0d-2) then
           atomtype(i) = 14
        else
           atomtype(i) = 0
        endif
     enddo
     !     - - - - - - - - - - - - - - - - - - - - -
     open(f_trj, file=trim(session_name)// '.mdtrj', &
          &         iostat=io, status='old', &
          &         access='sequential',form='unformatted')
     open(f_xyz, file=trim(session_name)// '.xyz', &
          &         iostat=io, status='replace', &
          &         access='sequential',form='formatted')
     open(f_cel, file=trim(session_name)// '.cellsize', &
          &         iostat=io, status='replace', &
          &         access='sequential',form='formatted')
     k = 0
     do while (.true.)
        k = k + 1
        read(f_trj, end=100) mdstep, dt
        !          write(*,*) mdstep, dt  ! for debug
        read(f_trj) natom
        read(f_trj) xyz(1:3,:), v(1:3,:)
        read(f_trj) n_nhc
        read(f_trj) rss, vss
        read(f_trj) n_nhc
        read(f_trj) rssb,vssb
        read(f_trj) cellx,celly,cellz,alpha,beta,gamma,vboxg

        write(f_cel,*) cellx*1.0d10, celly*1.0d10, cellz*1.0d10

        if (mod(k,nintvl) .ne. 0)  cycle
        ndisplay = min(ndisplay,natom)
        write(f_xyz,*) ndisplay
        write(f_xyz,*) '# ', mdstep
        xyz = xyz * 1.0d10
        do i=1,ndisplay
           if (atomtype(i) .eq. 1) then
              str_atom = 'H '
           else if (atomtype(i) .eq. 2) then
              str_atom = 'C '
           else if (atomtype(i) .eq. 3) then
              str_atom = 'N '
           else if (atomtype(i) .eq. 4) then
              str_atom = 'O '
           else if (atomtype(i) .eq. 5) then
              str_atom = 'P '
           else if (atomtype(i) .eq. 6) then
              str_atom = 'S '
           else if (atomtype(i) .eq. 7) then
              str_atom = 'He'
           else if (atomtype(i) .eq. 8) then
              str_atom = 'Ne'
           else if (atomtype(i) .eq. 9) then
              str_atom = 'Ca'
           else if (atomtype(i) .eq. 10) then
              str_atom = 'Zn'
           else if (atomtype(i) .eq. 11) then
              str_atom = 'Fe'
           else if (atomtype(i) .eq. 12) then
              str_atom = 'Na'
           else if (atomtype(i) .eq. 13) then
              str_atom = 'Cl'
           else if (atomtype(i) .eq. 14) then
              str_atom = 'K '
           else if (atomtype(i) .eq. 0) then
              str_atom = 'XX'
           endif
           write(f_xyz,'(a2,3e20.10)') str_atom, xyz(1:3,i)
        enddo
        !         write(f_xyz,'(2x,3e20.10)') alpha,beta,gamma
        !         write(f_xyz,'(2x,3e20.10)') cellx,celly,cellz
     enddo
100  continue
     close(f_trj)
     close(f_xyz)
     !     - - - - - - - - - - - - - - - - - - - - -
     deallocate(atomtype)
  endif
  call mpiend()
end program mdtrj2xyz
