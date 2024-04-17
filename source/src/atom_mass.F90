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
!! \brief  Module and subroutines to set atom mass.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to read/write mass information.
!! \author Kazushi FUJIMOTO
!<
module atom_mass
  implicit none
  real(8),allocatable :: mass(:),r_mass(:)
  real(8) :: totalmass

contains

!>
!! \brief  Subroutine to set mass value.
!! \author Kazushi FUJIMOTO
!<
  subroutine fmod_set_mol_mass(value)
    use md_const
    implicit none
    real(8) :: value
    integer(4),save :: i=0
    i = i + 1

    mass(i) = value * md_ATOMIC_MASS_UNIT
    r_mass(i) = 1.0d0/mass(i)
  end subroutine fmod_set_mol_mass
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to read molecule.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_mass
    use device_numbers, only : f_mdff
    use param
    use mpi_tool
    implicit none
    include 'mpif.h'
    integer(4) :: ierr

    if(myrank==0) read(f_mdff) mass
    call MPI_Bcast(mass, npara, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    r_mass(:) = 1.0d0 / mass(:)
  end subroutine read_mass
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write mass in mdff.bin.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_mass
    use device_numbers, only : f_mdff
    implicit none
    write(f_mdff) mass
  end subroutine write_mass
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write mass.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_mass
    implicit none
    write(*,*) '[write_mass]'
    write(*,*) mass
  end subroutine write_memory_mass

end module atom_mass
