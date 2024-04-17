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
!! \brief  Module and subroutines for converting from mdff to mdff.bin
!<
!>
!! \brief  Module for converting from mdff to mdff.bin
!! \author  Kazushi FUJIMOTO
!<
module file_mdff
  implicit none

contains

subroutine write_mdffbin_restart
  use angle
  use device_numbers
  use session_name_mod
  use mpi_tool
  use bond, only : write_md_charmm_a_bond
  use bond_morse, only: write_md_bond_morse
  use bond_morse2, only: write_md_bond_morse2
  !use bond_pcff, only : write_md_pcff_a_bond
  use angle_morse, only : write_md_angleA_morse, write_md_angleB_morse
  !use angle_pcff, only : write_md_pcff_angleA,write_md_pcff_angleB
  !use bond_bond, only : write_md_pcff_bbA, write_md_pcff_bbB
  !use bond_angle, only : write_md_pcff_baA, write_md_pcff_baB
  !use angle_angle, only : write_md_pcff_aaA, write_md_pcff_aaB
  use UB, only : write_md_charmm_a_ub
  use dihedral, only : write_md_charmm_dihedralA, write_md_charmm_dihedralB
  use CMAP
  use improper_torsion, only : write_md_charmm_itorsionA, write_md_charmm_itorsionB
  use void123, only : write_md_void
  use special14, only : write_md_LJ_special, write_md_coulomb_special
  use segments, only : write_segment
  use atom_mass, only : write_mass
  use param, only : write_parameter_number
  use coulomb_mod, only : write_chgv
  use lj_mod, only : write_epsilon_sqrt, write_R_half
  use file_utility, only : create_binary_file
  implicit none
  integer(4) :: dummyInt = 0

  if(myrank==0) then
     call create_binary_file(f_mdff, trim(session_name)//'_restart.mdff.bin')

     call write_forcefield
     call write_parameter_number
     call write_segment
     call write_molecule
     call write_mass
     call write_chgv
     call write_epsilon_sqrt
     call write_R_half
     call write_md_shake
     call write_md_void
     call write_md_LJ_special
     call write_md_coulomb_special
     call write_md_charmm_a_bond
     call write_md_bond_morse
     call write_md_bond_morse2
     !call write_md_pcff_a_bond
     write(f_mdff) dummyInt
     call write_md_charmm_angleA
     call write_md_charmm_angleB
     call write_md_angleA_morse
     call write_md_angleB_morse
     !call write_md_pcff_angleA
     write(f_mdff) dummyInt
     !call write_md_pcff_angleB
     write(f_mdff) dummyInt
     call write_md_charmm_a_ub
     !call write_md_pcff_bbA
     write(f_mdff) dummyInt
     !call write_md_pcff_bbB
     write(f_mdff) dummyInt
     !call write_md_pcff_baA
     write(f_mdff) dummyInt
     !call write_md_pcff_baB
     write(f_mdff) dummyInt
     !call write_md_pcff_aaA
     write(f_mdff) dummyInt
     !call write_md_pcff_aaB
     write(f_mdff) dummyInt
     call write_md_charmm_dihedralA
     call write_md_charmm_dihedralB
     call write_md_charmm_CMAPA
     call write_md_charmm_CMAPB
     call write_md_charmm_CMAPC
     call write_md_charmm_CMAPD
     call write_md_charmm_CMAPE
     call write_md_charmm_itorsionA
     call write_md_charmm_itorsionB
     !       call write_position_constrain   !NOT support from binary input
     call write_CMAPVersion
     close(f_mdff)
     write(*,*) 'restart mdff.bin was succesfully created!'
  endif
end subroutine write_mdffbin_restart

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to write mdff firle
!! \author Kazushi FUJIMOTO
!<
  subroutine write_mdffbin
    use angle
    use device_numbers
    use session_name_mod
    use mpi_tool
    use bond, only : write_md_charmm_a_bond
    use bond_morse, only: write_md_bond_morse
    use bond_morse2, only: write_md_bond_morse2
    !use bond_pcff, only : write_md_pcff_a_bond
    use angle_morse, only : write_md_angleA_morse, write_md_angleB_morse
    !use angle_pcff, only : write_md_pcff_angleA,write_md_pcff_angleB
    !use bond_bond, only : write_md_pcff_bbA, write_md_pcff_bbB
    !use bond_angle, only : write_md_pcff_baA, write_md_pcff_baB
    !use angle_angle, only : write_md_pcff_aaA, write_md_pcff_aaB
    use UB, only : write_md_charmm_a_ub
    use dihedral, only : write_md_charmm_dihedralA, write_md_charmm_dihedralB
    use CMAP
    use improper_torsion, only : write_md_charmm_itorsionA, write_md_charmm_itorsionB
    use void123, only : write_md_void
    use special14, only : write_md_LJ_special, write_md_coulomb_special
    use segments, only : write_segment
    use atom_mass, only : write_mass
    use param, only : write_parameter_number
    use coulomb_mod, only : write_chgv
    use lj_mod, only : write_epsilon_sqrt, write_R_half
    use file_utility, only : create_new_binary_file
    implicit none
    integer(4) :: dummyInt = 0

    if(myrank==0) then
       call create_new_binary_file(f_mdff, trim(session_name)//'.mdff.bin')

       call write_forcefield
       call write_parameter_number
       call write_segment
       call write_molecule
       call write_mass
       call write_chgv
       call write_epsilon_sqrt
       call write_R_half
       call write_md_shake
       call write_md_void
       call write_md_LJ_special
       call write_md_coulomb_special
       call write_md_charmm_a_bond
       call write_md_bond_morse
       call write_md_bond_morse2
       !call write_md_pcff_a_bond
       write(f_mdff) dummyInt
       call write_md_charmm_angleA
       call write_md_charmm_angleB
       call write_md_angleA_morse
       call write_md_angleB_morse
       !call write_md_pcff_angleA
       write(f_mdff) dummyInt
       !call write_md_pcff_angleB
       write(f_mdff) dummyInt
       call write_md_charmm_a_ub
       !call write_md_pcff_bbA
       write(f_mdff) dummyInt
       !call write_md_pcff_bbB
       write(f_mdff) dummyInt
       !call write_md_pcff_baA
       write(f_mdff) dummyInt
       !call write_md_pcff_baB
       write(f_mdff) dummyInt
       !call write_md_pcff_aaA
       write(f_mdff) dummyInt
       !call write_md_pcff_aaB
       write(f_mdff) dummyInt
       call write_md_charmm_dihedralA
       call write_md_charmm_dihedralB
       call write_md_charmm_CMAPA
       call write_md_charmm_CMAPB
       call write_md_charmm_CMAPC
       call write_md_charmm_CMAPD
       call write_md_charmm_CMAPE
       call write_md_charmm_itorsionA
       call write_md_charmm_itorsionB
       !       call write_position_constrain   !NOT support from binary input
       call write_CMAPVersion
       close(f_mdff)
       write(*,*) 'mdff.bin was succesfully created!'
    endif
  end subroutine write_mdffbin
!-------------------------------------------------------------------------
  subroutine write_memory_mdff   !!! not-used now, cut
    use angle
    use device_numbers
    use mpi_tool
    use bond, only : write_memory_md_charmm_a_bond
    use bond_morse, only: write_memory_md_bond_morse
    use bond_morse2, only: write_memory_md_bond_morse2
    !use bond_pcff, only : write_memory_md_pcff_a_bond
    !use angle_pcff, only : write_memory_md_pcff_angleA,write_memory_md_pcff_angleB
    use angle_morse, only : write_memory_md_angleA_morse,write_memory_md_angleB_morse
    use UB, only : write_memory_md_charmm_a_ub
    !use bond_bond, only : write_memory_md_pcff_bbA,write_memory_md_pcff_bbB
    !use bond_angle, only : write_memory_md_pcff_baA,write_memory_md_pcff_baB
    !use angle_angle, only : write_memory_md_pcff_aaA,write_memory_md_pcff_aaB
    use dihedral, only : write_memory_md_charmm_dihedralA, write_memory_md_charmm_dihedralB
    use CMAP
    use improper_torsion, only : write_memory_md_charmm_itorsionA, write_memory_md_charmm_itorsionB
    use void123, only : write_memory_md_void
    use special14, only : write_memory_md_LJ_special, write_memory_md_coulomb_special
    use segments, only : write_memory_segment
    use atom_mass, only : write_memory_mass
    use param, only : write_memory_parameter_number
    use coulomb_mod, only : write_memory_chgv
    use lj_mod, only : write_memory_epsilon_sqrt, write_memory_R_half
    implicit none
    integer(4) :: dummyInt = 0

    if(myrank==0) then
       call write_memory_forcefield
       call write_memory_parameter_number
       call write_memory_segment
       call write_memory_molecule
       call write_memory_mass
       call write_memory_chgv
       call write_memory_epsilon_sqrt
       call write_memory_R_half
       call write_memory_md_shake
       call write_memory_md_void
       call write_memory_md_LJ_special
       call write_memory_md_coulomb_special
       call write_memory_md_charmm_a_bond
       call write_memory_md_bond_morse
       call write_memory_md_bond_morse2
       !call write_memory_md_pcff_a_bond
       write(*,*) dummyInt
       call write_memory_md_charmm_angleA
       call write_memory_md_charmm_angleB
       call write_memory_md_angleA_morse
       call write_memory_md_angleB_morse
       !call write_memory_md_pcff_angleA
       write(*,*) dummyInt
       !call write_memory_md_pcff_angleB
       write(*,*) dummyInt
       call write_memory_md_charmm_a_ub
       !call write_memory_md_pcff_bbA
       write(*,*) dummyInt
       !call write_memory_md_pcff_bbB
       write(*,*) dummyInt
       !call write_memory_md_pcff_baA
       write(*,*) dummyInt
       !call write_memory_md_pcff_baB
       write(*,*) dummyInt
       !call write_memory_md_pcff_aaA
       write(*,*) dummyInt
       !call write_memory_md_pcff_aaB
       write(*,*) dummyInt
       call write_memory_md_charmm_dihedralA
       call write_memory_md_charmm_dihedralB
       call write_memory_md_charmm_CMAPA
       call write_memory_md_charmm_CMAPB
       call write_memory_md_charmm_CMAPC
       call write_memory_md_charmm_CMAPD
       call write_memory_md_charmm_CMAPE
       call write_memory_md_charmm_itorsionA
       call write_memory_md_charmm_itorsionB
       !       call write_memory_position_constrain   !NOT supported yet
       call write_memory_CMAPVersion
    endif
  end subroutine write_memory_mdff
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to read force field
!! \author Kazushi FUJIMOTO
!<
  subroutine read_forcefield
    use md_condition
    use device_numbers
    use CMAP
    use mpi_tool
    use force_field_numbers
    implicit none
    include 'mpif.h'
    integer(4) :: ierr

    if(myrank==0) then
       read(f_mdff) md_condition__force_field
       if(md_condition__force_field.ne.CHARMM)then
          write(*,*) 'Warrning: non-CHARMM potential is used!'
       endif
    endif
    call MPI_Bcast(md_condition__force_field,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)

    return
  end subroutine read_forcefield
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write force field mdff.bin
!! \author Kazushi FUJIMOTO
!<
  subroutine write_forcefield
    use md_condition
    use device_numbers
    implicit none
    write(f_mdff) md_condition__force_field
    return
  end subroutine write_forcefield
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write memory forcefiled
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_forcefield
    use md_condition
    implicit none
    write(*,*) '[write_forcefield]'
    write(*,*) md_condition__force_field
    return
  end subroutine write_memory_forcefield
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to read parameter number
!! \author Kazushi FUJIMOTO
!<
  subroutine read_parameter_number
    use trajectory_org
    use param
    use device_numbers
    use mpi_tool
    implicit none
    include 'mpif.h'
    integer(4) :: ierr

    if(myrank==0) then
       read(f_mdff) npara
    endif
    call MPI_Bcast(npara,  1, MPI_INTEGER4, 0,  MPI_COMM_WORLD, ierr)
    allocate(paranum(n))

    if(myrank==0) read(f_mdff) paranum
    call MPI_Bcast(paranum, n, MPI_INTEGER4, 0,  MPI_COMM_WORLD, ierr)
  end subroutine read_parameter_number
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to read molecule.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_molecule
    use trajectory_org
    use molecules
    use device_numbers
    use mpi_tool
    implicit none
    include 'mpif.h'
    integer(4) :: ierr

    if(myrank==0) read(f_mdff) md_periodic__howmany_molecules
    call MPI_Bcast(md_periodic__howmany_molecules, 1, &
         &               MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)

    call fmod_alloc_mol_natoms(md_periodic__howmany_molecules)
    call fmod_alloc_moltop(md_periodic__howmany_molecules)
    if(myrank==0) then
       read(f_mdff) molecules__howmany_atoms
       read(f_mdff) moltop
       !       read(f_mdff) mol2atom
       !       read(f_mdff) atom2mol
    endif
    call MPI_Bcast(molecules__howmany_atoms, &
         &               md_periodic__howmany_molecules, MPI_INTEGER4, 0, &
         &               MPI_COMM_WORLD, ierr)
    call MPI_Bcast(moltop, &
         &               md_periodic__howmany_molecules, MPI_INTEGER4, 0, &
         &               MPI_COMM_WORLD, ierr)
  end subroutine read_molecule
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write molecule to mdff.bin File.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_molecule
    use molecules
    use device_numbers
    implicit none
    write(f_mdff) md_periodic__howmany_molecules
    write(f_mdff) molecules__howmany_atoms
    write(f_mdff) moltop
  end subroutine write_molecule
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write molecule.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_molecule
    use molecules
    implicit none
    write(*,*) '[write_molecule]'
    write(*,*) md_periodic__howmany_molecules
    write(*,*) molecules__howmany_atoms
    write(*,*) moltop
    return
  end subroutine write_memory_molecule
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to read shake.
!! \author Kazushi FUJIMOTO
!<
  subroutine read_md_shake
    use trajectory_org
    use shake_rattle_roll
    use md_condition
    use device_numbers
    use param
    !pshake
    use pshake_init
    use mpi_tool
    implicit none
    include 'mpif.h'
    integer(4) :: ierr
    integer(4) :: rtotnconst, rngrp
    integer(4) :: i,j

    allocate(ShakeGroupLeader(npara))
    if(myrank==0) read(f_mdff) rtotnconst
    call MPI_Bcast(rtotnconst, 1, MPI_INTEGER4, &
         &               0, MPI_COMM_WORLD, ierr)
    call fmod_totnconst(rtotnconst)
    if(rtotnconst==0)  then
       rngrp_ps = 0
       allocate(nconstraints(rngrp_ps))
       allocate(atom1S(rngrp_ps,10), atom2S(rngrp_ps,10))
       allocate(slength(rngrp_ps,10))
       return
    endif
    if(myrank==0) read(f_mdff) rngrp
    call MPI_Bcast(rngrp, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)

    if(myrank==0)  then
       do i = 1, npara
          read(f_mdff) ShakeGroupLeader(i)
       enddo
    endif
    call MPI_Bcast(ShakeGroupLeader, npara, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)

    allocate(nconstraints(rngrp))
    if(myrank==0)  then
       do i = 1, rngrp
          read(f_mdff) nconstraints(i)
       enddo
    endif
    call MPI_Bcast(nconstraints, rngrp, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)

    allocate(atom1S(rngrp,10), atom2S(rngrp,10))
    allocate(slength(rngrp,10))
    atom1S(:,:) = -1
    atom2S(:,:) = -1
    slength(:,:) = -1.0d0
    if(myrank==0)  then
       do i = 1, rngrp
          do j = 1, nconstraints(i)
             read(f_mdff) atom1S(i,j), atom2S(i,j), slength(i,j)
          enddo
       enddo
    endif
    call MPI_Bcast(atom1S, rngrp*10, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(atom2S, rngrp*10, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(slength, rngrp*10, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    !pshake
    rngrp_ps = rngrp
  end subroutine read_md_shake
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write shake in mdff.bin File.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_md_shake
    use shake_rattle_roll
    use device_numbers
    use param
    use pshake_init
    implicit none
    integer(4) :: rngrp, i, j
    write(f_mdff) totnconst
    if(totnconst==0)  then
       return
    endif
    rngrp = size(nconstraints(:))
    write(f_mdff) rngrp
    do i = 1, npara
       write(f_mdff) ShakeGroupLeader(i)
    enddo
    do i = 1, rngrp
       write(f_mdff) nconstraints(i)
    enddo
    do i = 1, rngrp
       do j = 1, nconstraints(i)
          write(f_mdff) atom1S(i,j), atom2S(i,j), slength(i,j)
       enddo
    enddo
  end subroutine write_md_shake
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to write shake.
!! \author Kazushi FUJIMOTO
!<
  subroutine write_memory_md_shake
    use shake_rattle_roll
    use param
    use pshake_init
    implicit none
    integer(4) :: rngrp, i, j
    write(*,*) '[write_md_shake]'
    write(*,*) totnconst
    if(totnconst==0)  then
       return
    endif
    rngrp = size(nconstraints(:))
    write(*,*) rngrp
    do i = 1, npara
       write(*,*) ShakeGroupLeader(i)
    enddo
    do i = 1, rngrp
       write(*,*) nconstraints(i)
    enddo
    do i = 1, rngrp
       do j = 1, nconstraints(i)
          write(*,*) atom1S(i,j), atom2S(i,j), slength(i,j)
       enddo
    enddo
  end subroutine write_memory_md_shake

end module file_mdff
