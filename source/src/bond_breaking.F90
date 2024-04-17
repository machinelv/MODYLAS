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
!! \brief  Module and subroutines to calculate bond-breaking potential.
!<
!>
!! \brief  Module to control bond-breaking events.
!! \author Zhiye Tang
!<
!----------------------------------------------------------------------
module bond_breaking_event
    use bond_morse2, only: nbondmorse2_as
    use angle_morse, only : nanglemorseAs
    use dihedral, only : ndihedralAs
    use improper_torsion, only : nitorsionAs
    implicit none
    ! for omp local
    integer(4),allocatable :: brokenBondListLocal(:)
    ! for node
    integer(4),allocatable :: brokenBondList(:)
    ! global
    integer(4),allocatable :: brokenBondListG(:)
    integer(4),allocatable :: brokenBondNodeNbondMPI(:)
    integer(4),allocatable :: brokenBondMPIBuffer(:)
    
    integer(4) :: MAX_BROKEN_BOND_LOCAL=20
    integer(4) :: MAX_BROKEN_BOND=100
    integer(4) :: MAX_BROKEN_BOND_G=200

    real(8) :: BROKENMASK = -1.0d0
    real(8) :: BROKENMASKDIHERAL = 0.0d0
    contains
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to allocate memory for bond breaking detection
!! \author Zhiye Tang
!< 
    subroutine bond_breaking_event_allocate()
        use mpi_tool, only : nprocs
        implicit none
        allocate(brokenBondListLocal(MAX_BROKEN_BOND_LOCAL*2))
        allocate(brokenBondList(MAX_BROKEN_BOND*2))
        allocate(brokenBondListG(MAX_BROKEN_BOND_G*2))
        allocate(brokenBondNodeNbondMPI(nprocs))
        allocate(brokenBondMPIBuffer(MAX_BROKEN_BOND_LOCAL*2))
    end subroutine
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to detect bond breaking
!! \author Zhiye Tang
!< 
    subroutine detect_bond_breaking()
        use subcell
        use bond_morse2
        use param
        use mpi_tool
        use trajectory_mpi
        use md_condition, only : allow_bond_breaking
        use boundary
        implicit none
        integer(4) :: iam, atomIDJ, k0, i0, j0, iA, iB, ipar, ii, jj
        real(8) :: drAB__x,drAB__y,drAB__z, b2,cutoff2
    
        integer(4) :: localNBrokenBond, memPos
        ! num of bond broken in this node
        integer(4) :: nBrokenBond
        ! global
        integer(4) :: nBrokenBondG
        ! ntotalBondBrokenNode = 1 if bond broken in this node, 0 not
        ! ntotalBondBrokenNodeMPI, for reduce
        integer(4) :: ntotalBondBrokenNode,ntotalBondBrokenNodeMPI
        ! in case of bond breaking on only node, rank
        integer(4) :: oneNodeID, oneNodeIDMPI
        ! atomID
        integer(4) :: atomID1, atomID2
        include 'mpif.h'

        if(allow_bond_breaking<=0.0d0) return
! skip this subroutine if no bond morse2 is defined in the system
        if(nbondmorse2_as.eq.0) return
        iam = 0

        cutoff2 = allow_bond_breaking * allow_bond_breaking

        localNBrokenBond = 0
        nBrokenBond = 0
        memPos = 0

!$omp parallel default(none) &
!$omp& private(k0,i0,j0,iA,iB,ipar,atomIDJ) &
!$omp& shared(nselfseg,lsegtop,lseg_natoms) &
!$omp& shared(m2i,paranum,nA,atom2A,i2m,wkxyz,cutoff2) &
!$omp& private(drAB__x,drAB__y,drAB__z, b2, ii) &
!$omp& private(brokenBondListLocal) &
!$omp& firstprivate(localNBrokenBond) &
!$omp& shared(memPos, brokenBondList, DA) &
!$omp& shared(MAX_BROKEN_BOND, MAX_BROKEN_BOND_LOCAL) &
!$omp& reduction(+:nBrokenBond) &
!$omp& shared(myrank)

!$omp do
        do k0=1,nselfseg
            do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
                iA   = m2i(i0) ! iA is the global atomID
                ipar = paranum(iA)
                do iB = 1, nA(ipar)
                    atomIDJ=atom2A(ipar,iB)+(iA-ipar)
                    j0 = i2m(atomIDJ)
!default: MTD
                    if(iA < atomIDJ)cycle ! avoid double count
                    if(DA(ipar, iB)<0.0d0) then
                        ! already removed
                        !write(*,*) iA, atomIDJ, "is already removed"
                        cycle
                    endif

                    drAB__x = wkxyz(1,i0) - wkxyz(1,j0)
                    drAB__y = wkxyz(2,i0) - wkxyz(2,j0)
                    drAB__z = wkxyz(3,i0) - wkxyz(3,j0)
#ifdef ONEPROC_AXIS
                    call pbc_pair(drAB__x,drAB__y,drAB__z)
#endif
                    b2 = drAB__x*drAB__x+drAB__y*drAB__y+drAB__z*drAB__z
                    if(b2 > cutoff2) then
                        localNBrokenBond = localNBrokenBond + 1
                        if(localNBrokenBond>MAX_BROKEN_BOND_LOCAL)then
                            write(0,*) "Error: localNBrokenBond>MAX_BROKEN_BOND_LOCAL"
                        endif
                        brokenBondListLocal(localNBrokenBond*2-1) = iA
                        brokenBondListLocal(localNBrokenBond*2  ) = atomIDJ
                        !write(*,*) "Debug ", myrank, iA, atomIDJ
                    endif
                enddo
            enddo
        enddo
!$omp end do
        nBrokenBond = nBrokenBond + localNBrokenBond

!$omp critical
        do ii=1, localNBrokenBond
            memPos = memPos + 1
            if(memPos>MAX_BROKEN_BOND)then
                write(0,*) "Error: memPos>MAX_BROKEN_BOND"
            endif
            brokenBondList(memPos*2-1) = brokenBondListLocal(ii*2-1)
            brokenBondList(memPos*2  ) = brokenBondListLocal(ii*2  )
        enddo
!$omp end critical
!$omp end parallel

        ! nBrokenBond
        ! brokenBondList
        if( nBrokenBond > 0 ) then
            ntotalBondBrokenNode = 1
        else
            ntotalBondBrokenNode = 0
        endif

        ! communicate the total number of nodes that have bond breaking
        ! into ntotalBondBrokenNodeMPI
        call mpi_allreduce(ntotalBondBrokenNode, ntotalBondBrokenNodeMPI, &
      &      1, MPI_INT, MPI_SUM, mpi_comm_world, ipar)
        
        ! communicate bond breaking information into
        ! nBrokenBondG
        ! brokenBondListG
        if( ntotalBondBrokenNodeMPI == 0 )then
            nBrokenBondG = 0
        else if( ntotalBondBrokenNodeMPI == 1 )then
            ! bcast the node rank
            if(nBrokenBond==0)then
                ! receiving nodes
                oneNodeID = 0
            else
                ! sending nodes
                oneNodeID = myrank + 1
            endif
            
            call mpi_allreduce(oneNodeID, oneNodeIDMPI, &
      &         1, MPI_INT, MPI_SUM, mpi_comm_world, ipar)

            oneNodeIDMPI = oneNodeIDMPI - 1


            nBrokenBondG = nBrokenBond
            ! bcast nBrokenBondG
            call MPI_Bcast(nBrokenBondG, 1, MPI_INT, oneNodeIDMPI, mpi_comm_world, ipar)

            ! on rank = myrank
            ! copy brokenBondList to brokenBondMPIBuffer
            if( oneNodeIDMPI == myrank ) then
                do jj=1, nBrokenBondG
                    brokenBondMPIBuffer(jj*2-1) = brokenBondList(jj*2-1)
                    brokenBondMPIBuffer(jj*2  ) = brokenBondList(jj*2)
                enddo
            endif
            ! bcast brokenBondList
            call MPI_Bcast(brokenBondMPIBuffer, nBrokenBondG * 2, MPI_INT, oneNodeIDMPI, &
      &         mpi_comm_world, ipar)
            ! copy brokenBondMPIBuffer to brokenBondListG
            if(nBrokenBondG > MAX_BROKEN_BOND_G) then
                write(0,*) "Error: nBrokenBondG>MAX_BROKEN_BOND_G"
            endif
            do jj=1, nBrokenBondG
                brokenBondListG(jj*2-1) = brokenBondMPIBuffer(jj*2-1)
                brokenBondListG(jj*2  ) = brokenBondMPIBuffer(jj*2)
            enddo
        else if( ntotalBondBrokenNodeMPI > 1 ) then
            ! all reduce nBrokenBondG
            call mpi_allreduce(nBrokenBond, nBrokenBondG, 1, &
      &         MPI_INT, MPI_SUM, mpi_comm_world, ipar)
            if(nBrokenBondG > MAX_BROKEN_BOND_G) then
                write(0,*) "Error: nBrokenBondG>MAX_BROKEN_BOND_G"
            endif
            ! gather nBrokenBond into brokenBondNodeNbondMPI
            call mpi_allgather(nBrokenBond, 1, MPI_INT, brokenBondNodeNbondMPI, &
      &         1, MPI_INT, mpi_comm_world, ipar)
            memPos = 0
            do ii=1, nprocs
                if(brokenBondNodeNbondMPI(ii)>0) then
                    ! on rank = myrank
                    ! copy brokenBondList to brokenBondMPIBuffer
                    ! note that MPI rank is ii-1
                    if( ii-1 == myrank ) then
                        do jj=1, nBrokenBond
                            brokenBondMPIBuffer(jj*2-1) = brokenBondList(jj*2-1)
                            brokenBondMPIBuffer(jj*2  ) = brokenBondList(jj*2)
                        enddo
                    endif
                    ! bcast brokenBondMPIBuffer
                    call MPI_Bcast(brokenBondMPIBuffer, brokenBondNodeNbondMPI(ii) * 2, &
      &                 MPI_INT, ii-1, mpi_comm_world, ipar)
                    ! copy brokenBondMPIBuffer to brokenBondListG
                    ! copy brokenBondMPIBuffer to 
                    do jj=1, brokenBondNodeNbondMPI(ii)
                        brokenBondListG(memPos+jj*2-1) = brokenBondMPIBuffer(jj*2-1)
                        brokenBondListG(memPos+jj*2  ) = brokenBondMPIBuffer(jj*2)
                    enddo
                    memPos = memPos + brokenBondNodeNbondMPI(ii) * 2
                endif
            enddo
        endif

        ! now all nodes should be on the same page
        ! let them update the topology individually
        ! nBrokenBondG
        ! brokenBondListG
        do ii=1,nBrokenBondG
            atomID1 = brokenBondListG(ii*2-1)
            atomID2 = brokenBondListG(ii*2)
            if(myrank==0) then
                call bond_breaking_event_print(atomID1, atomID2)
            endif
            if(nbondmorse2_as.ne.0) then
                call bond_breaking_event_bond_morse2(atomID1, atomID2)
            endif
            if(nanglemorseAs.ne.0) then
                call bond_breaking_event_angle_morse(atomID1, atomID2)
            endif
            if(ndihedralAs.ne.0) then
                call bond_breaking_event_dihedral(atomID1, atomID2)
            endif
            if(nitorsionAs.ne.0) then
                call bond_breaking_event_improper_torsion(atomID1, atomID2)
            endif
        enddo

    end subroutine

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to detect bond breaking of bond morse2
!! \author Zhiye Tang
!<  
    subroutine bond_breaking_event_bond_morse2(atomID1, atomID2)
        implicit none
        integer(4) :: atomID1, atomID2
        ! atomID1 and atomID2 are global atom ID
  
        call bond_breaking_event_bond_morse2_direct(atomID1, atomID2)
        call bond_breaking_event_bond_morse2_direct(atomID2, atomID1)
  
    end subroutine
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to detect bond breaking of bond morse2 direct
!! \author Zhiye Tang
!<   
    subroutine bond_breaking_event_bond_morse2_direct(atomID1, atomID2)
        use bond_morse2
        use param
        implicit none
        integer(4) :: atomID1, atomID2
        ! atomID1 and atomID2 are global atom ID
        integer(4) :: i1, atomloop2, ipar1
  
        ipar1 = paranum(atomID1)
  
        do i1=1, nA(ipar1)
            ! atomloop2 is the global atom ID of the interacting atom
            atomloop2 = atom2A(ipar1,i1)+(atomID1-ipar1)
            if(atomloop2 == atomID2) then
                ! this is the broken bond
                call bond_breaking_event_bond_morse2_setA(ipar1, i1)
                call bond_breaking_event_bond_morse2_print(atomID1, atomID1, atomID2)
            endif
        enddo
  
    end subroutine
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to detect bond breaking of bond morse2 setA
!! \author Zhiye Tang
!< 
    subroutine bond_breaking_event_bond_morse2_setA(ipar1, i1)
        use bond_morse2
        implicit none
        integer(4) :: ipar1, i1
  
        DA(ipar1, i1) = BROKENMASK
  
    end subroutine
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to print bond breaking information of bond morse2
!! \author Zhiye Tang
!< 
    subroutine bond_breaking_event_bond_morse2_print(atomID, atomID1, atomID2)
        use mpi_tool
        implicit none
        integer(4) :: atomID1, atomID2, atomID

#ifdef PRINTBROKEN
        if(myrank==0) write(*,'(A17I8A1I8I8)') "Bond morse2 Atom[", atomID, "]", atomID1, atomID2
#endif
    end subroutine
  
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to process bond breaking of angle morse
!! \author Zhiye Tang
!<  
    subroutine bond_breaking_event_angle_morse(atomID1, atomID2)
        implicit none
        integer(4) :: atomID1, atomID2
        ! atomID1 and atomID2 are global atom ID
  
        call bond_breaking_event_angle_morse_direct(atomID1, atomID2)
        call bond_breaking_event_angle_morse_direct(atomID2, atomID1)
  
    end subroutine
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to process bond breaking of angle morse direct
!! \author Zhiye Tang
!<  
    subroutine bond_breaking_event_angle_morse_direct(atomID1, atomID2)
        use angle_morse
        use param
        implicit none
        integer(4) :: atomID1, atomID2
        integer(4) :: ipar1
        integer(4) :: i1, atomloopA2, atomloopA3, atomloopB1, atomloopB3
  
        ipar1 = paranum(atomID1)
  
        ! do atomID1, A series
        do i1=1, nA(ipar1)
            atomloopA2 = atom2A(ipar1,i1)+(atomID1-ipar1)
            atomloopA3 = atom3A(ipar1,i1)+(atomID1-ipar1)
            if(atomloopA2 == atomID2) then
                KmA1(ipar1, i1) = BROKENMASK
                KmA2(ipar1, i1) = BROKENMASK
                KmA3(ipar1, i1) = BROKENMASK
                KmA4(ipar1, i1) = BROKENMASK
                call bond_breaking_event_angle_morse_print(atomID1, atomID1, atomloopA2, atomloopA3)
                ! need to do atomloopA3
                call bond_breaking_event_angle_morse_indirect(atomloopA3, atomID1,atomID2)
            else if( atomloopA3 == atomID2) then
                KmA1(ipar1, i1) = BROKENMASK
                KmA2(ipar1, i1) = BROKENMASK
                KmA3(ipar1, i1) = BROKENMASK
                KmA4(ipar1, i1) = BROKENMASK
                call bond_breaking_event_angle_morse_print(atomID1, atomID1, atomloopA2, atomloopA3)
                ! need to do atomloopA2
                call bond_breaking_event_angle_morse_indirect(atomloopA2, atomID1,atomID2)
            endif
        enddo
  
        ! do atomID1, B series
        do i1=1, nB(ipar1)
            atomloopB1 = atom1B(ipar1,i1)+(atomID1-ipar1)
            atomloopB3 = atom3B(ipar1,i1)+(atomID1-ipar1)
            if(atomloopB1 == atomID2) then
                KmB1(ipar1, i1) = BROKENMASK
                KmB2(ipar1, i1) = BROKENMASK
                KmB3(ipar1, i1) = BROKENMASK
                KmB4(ipar1, i1) = BROKENMASK
                call bond_breaking_event_angle_morse_print(atomID1, atomloopB1, atomID1, atomloopB3)
                ! need to do atomloopB3
                call bond_breaking_event_angle_morse_indirect(atomloopB3, atomID1,atomID2)
            else if( atomloopB3 == atomID2) then
                KmB1(ipar1, i1) = BROKENMASK
                KmB2(ipar1, i1) = BROKENMASK
                KmB3(ipar1, i1) = BROKENMASK
                KmB4(ipar1, i1) = BROKENMASK
                call bond_breaking_event_angle_morse_print(atomID1, atomloopB1, atomID1, atomloopB3)
                ! need to do atomloopB1
                call bond_breaking_event_angle_morse_indirect(atomloopB1, atomID1,atomID2)
            endif
        enddo
        
    end subroutine
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to process bond breaking of angle morse indirect
!! \author Zhiye Tang
!<  
    subroutine bond_breaking_event_angle_morse_indirect(atomID, atomID1, atomID2)
        use angle_morse
        use param
        implicit none
        integer(4) :: atomID1, atomID2, atomID
        integer(4) :: ipar1
        integer(4) :: i1, atomloopA2, atomloopA3, atomloopB1, atomloopB3
  
        ipar1 = paranum(atomID)
  
        ! do atomID, A series
        do i1=1, nA(ipar1)
            atomloopA2 = atom2A(ipar1,i1)+(atomID-ipar1)
            atomloopA3 = atom3A(ipar1,i1)+(atomID-ipar1)
            if( (atomloopA2 == atomID1 .and. atomloopA3 == atomID2) .or. &
      &         (atomloopA2 == atomID2 .and. atomloopA3 == atomID1)) then
                KmA1(ipar1, i1) = BROKENMASK
                KmA2(ipar1, i1) = BROKENMASK
                KmA3(ipar1, i1) = BROKENMASK
                KmA4(ipar1, i1) = BROKENMASK
                call bond_breaking_event_angle_morse_print(atomID, atomID, atomloopA2, atomloopA3)
            endif
        enddo
  
        ! do atomID, B series
        do i1=1, nB(ipar1)
            atomloopB1 = atom1B(ipar1,i1)+(atomID-ipar1)
            atomloopB3 = atom3B(ipar1,i1)+(atomID-ipar1)
            if( (atomloopB1 == atomID1 .and. atomloopB3 == atomID2) .or. &
      &         (atomloopB1 == atomID2 .and. atomloopB3 == atomID1)) then
                KmB1(ipar1, i1) = BROKENMASK
                KmB2(ipar1, i1) = BROKENMASK
                KmB3(ipar1, i1) = BROKENMASK
                KmB4(ipar1, i1) = BROKENMASK
                call bond_breaking_event_angle_morse_print(atomID, atomloopB1, atomID, atomloopB3)
            endif
        enddo
  
    end subroutine
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to print bond breaking information of angle morse
!! \author Zhiye Tang
!<  
    subroutine bond_breaking_event_angle_morse_print(atomID, atomID1, atomID2, atomID3)
        use mpi_tool
        implicit none
        integer(4) :: atomID1, atomID2, atomID3, atomID

#ifdef PRINTBROKEN
        if(myrank==0) write(*,'(A17I8A1I8I8I8)') "Angle morse Atom[", atomID, "]", atomID1, atomID2, atomID3
#endif
  
    end subroutine
  
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to process bond breaking of dihedral
!! \author Zhiye Tang
!<   
    subroutine bond_breaking_event_dihedral(atomID1, atomID2)
        implicit none
        integer(4) :: atomID1, atomID2
        ! atomID1 and atomID2 are global atom ID
  
        call bond_breaking_event_dihedral_direct(atomID1, atomID2)
        call bond_breaking_event_dihedral_direct(atomID2, atomID1)
  
    end subroutine
  
#ifdef DIHEDRAL_TABLE
    subroutine bond_breaking_event_dihedral_direct(atomID1, atomID2)
        use dihedral
        use param
        implicit none
        integer(4) :: atomID1, atomID2
        integer(4) :: ipar1
        integer(4) :: i1, atomloopA2, atomloopA3, atomloopA4
        integer(4) :: atomloopB1, atomloopB3, atomloopB4
  
        ipar1 = paranum(atomID1)
  
        ! do atomID1, A series
        do i1=1, nA_t(ipar1)
            atomloopA2 = atom2A_t(ipar1,i1)+(atomID1-ipar1)
            atomloopA3 = atom3A_t(ipar1,i1)+(atomID1-ipar1)
            atomloopA4 = atom4A_t(ipar1,i1)+(atomID1-ipar1)
            if(atomloopA2 == atomID2) then
                indexA_t(ipar1, i1) = -1
                call bond_breaking_event_dihedral_print(atomID1, atomID1, atomloopA2, atomloopA3, atomloopA4)
                ! need to do atomloopA3
                call bond_breaking_event_dihedral_indirect(atomloopA3, atomID1,atomID2)
                ! need to do atomloopA4
                call bond_breaking_event_dihedral_indirect(atomloopA4, atomID1,atomID2)
            else if( atomloopA3 == atomID2) then
                indexA_t(ipar1, i1) = -1
                call bond_breaking_event_dihedral_print(atomID1, atomID1, atomloopA2, atomloopA3, atomloopA4)
                ! need to do atomloopA2
                call bond_breaking_event_dihedral_indirect(atomloopA2, atomID1,atomID2)
                ! need to do atomloopA4
                call bond_breaking_event_dihedral_indirect(atomloopA4, atomID1,atomID2)
            else if( atomloopA4 == atomID2) then
                indexA_t(ipar1, i1) = -1
                call bond_breaking_event_dihedral_print(atomID1, atomID1, atomloopA2, atomloopA3, atomloopA4)
                ! need to do atomloopA2
                call bond_breaking_event_dihedral_indirect(atomloopA2, atomID1,atomID2)
                ! need to do atomloopA3
                call bond_breaking_event_dihedral_indirect(atomloopA3, atomID1,atomID2)
            endif
        enddo
  
        ! do atomID1, B series
        do i1=1, nB_t(ipar1)
            atomloopB1 = atom1B_t(ipar1,i1)+(atomID1-ipar1)
            atomloopB3 = atom3B_t(ipar1,i1)+(atomID1-ipar1)
            atomloopB4 = atom4B_t(ipar1,i1)+(atomID1-ipar1)
            if(atomloopB1 == atomID2) then
                indexB_t(ipar1, i1) = -1
                call bond_breaking_event_dihedral_print(atomID1, atomloopB1, atomID1, atomloopB3, atomloopB4)
                ! need to do atomloopB3
                call bond_breaking_event_dihedral_indirect(atomloopB3, atomID1,atomID2)
                ! need to do atomloopB4
                call bond_breaking_event_dihedral_indirect(atomloopB4, atomID1,atomID2)
            else if( atomloopB3 == atomID2) then
                indexB_t(ipar1, i1) = -1
                call bond_breaking_event_dihedral_print(atomID1, atomloopB1, atomID1, atomloopB3, atomloopB4)
                ! need to do atomloopB1
                call bond_breaking_event_dihedral_indirect(atomloopB1, atomID1,atomID2)
                ! need to do atomloopB4
                call bond_breaking_event_dihedral_indirect(atomloopB4, atomID1,atomID2)
            else if( atomloopB4 == atomID2) then
                indexB_t(ipar1, i1) = -1
                call bond_breaking_event_dihedral_print(atomID1, atomloopB1, atomID1, atomloopB3, atomloopB4)
                ! need to do atomloopB1
                call bond_breaking_event_dihedral_indirect(atomloopB1, atomID1,atomID2)
                ! need to do atomloopB3
                call bond_breaking_event_dihedral_indirect(atomloopB3, atomID1,atomID2)
            endif
        enddo
        
    end subroutine
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to process bond breaking of dihedral indirect
!! \author Zhiye Tang
!<  
    subroutine bond_breaking_event_dihedral_indirect(atomID, atomID1, atomID2)
        use dihedral
        use param
        implicit none
        integer(4) :: atomID1, atomID2, atomID
        integer(4) :: ipar1
        integer(4) :: i1, atomloopA2, atomloopA3, atomloopA4
        integer(4) :: atomloopB1, atomloopB3, atomloopB4
  
        ipar1 = paranum(atomID)
  
        ! do atomID1, A series
        do i1=1, nA_t(ipar1)
            atomloopA2 = atom2A_t(ipar1,i1)+(atomID-ipar1)
            atomloopA3 = atom3A_t(ipar1,i1)+(atomID-ipar1)
            atomloopA4 = atom4A_t(ipar1,i1)+(atomID-ipar1)
            if( (atomloopA2 == atomID1 .and. atomloopA3 == atomID2) .or. &
      &         (atomloopA2 == atomID2 .and. atomloopA3 == atomID1) .or. &
      &         (atomloopA3 == atomID1 .and. atomloopA4 == atomID2) .or. &
      &         (atomloopA3 == atomID2 .and. atomloopA4 == atomID1) ) then
                indexA_t(ipar1, i1) = -1
                call bond_breaking_event_dihedral_print(atomID, atomID, atomloopA2, atomloopA3, atomloopA4)
            endif
        enddo
  
        ! do atomID1, B series
        do i1=1, nB_t(ipar1)
            atomloopB1 = atom1B_t(ipar1,i1)+(atomID-ipar1)
            atomloopB3 = atom3B_t(ipar1,i1)+(atomID-ipar1)
            atomloopB4 = atom4B_t(ipar1,i1)+(atomID-ipar1)
            if( (atomloopB3 == atomID1 .and. atomloopB4 == atomID2) .or. &
      &         (atomloopB3 == atomID2 .and. atomloopB4 == atomID1) ) then
                indexB_t(ipar1, i1) = -1
                call bond_breaking_event_dihedral_print(atomID, atomloopB1, atomID, atomloopB3, atomloopB4)
            endif
        enddo
        
    end subroutine

#else /* DIHEDRAL_TABLE */
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to process bond breaking of dihedral direct
!! \author Zhiye Tang
!<  
    subroutine bond_breaking_event_dihedral_direct(atomID1, atomID2)
        use dihedral
        use param
        implicit none
        integer(4) :: atomID1, atomID2
        integer(4) :: ipar1
        integer(4) :: i1, atomloopA2, atomloopA3, atomloopA4
        integer(4) :: atomloopB1, atomloopB3, atomloopB4
  
        ipar1 = paranum(atomID1)
  
        ! do atomID1, A series
        do i1=1, nA(ipar1)
            atomloopA2 = atom2A(ipar1,i1)+(atomID1-ipar1)
            atomloopA3 = atom3A(ipar1,i1)+(atomID1-ipar1)
            atomloopA4 = atom4A(ipar1,i1)+(atomID1-ipar1)
            if(atomloopA2 == atomID2) then
                KdA(ipar1, i1) = BROKENMASKDIHERAL
                call bond_breaking_event_dihedral_print(atomID1, atomID1, atomloopA2, atomloopA3, atomloopA4)
                ! need to do atomloopA3
                call bond_breaking_event_dihedral_indirect(atomloopA3, atomID1,atomID2)
                ! need to do atomloopA4
                call bond_breaking_event_dihedral_indirect(atomloopA4, atomID1,atomID2)
            else if( atomloopA3 == atomID2) then
                KdA(ipar1, i1) = BROKENMASKDIHERAL
                call bond_breaking_event_dihedral_print(atomID1, atomID1, atomloopA2, atomloopA3, atomloopA4)
                ! need to do atomloopA2
                call bond_breaking_event_dihedral_indirect(atomloopA2, atomID1,atomID2)
                ! need to do atomloopA4
                call bond_breaking_event_dihedral_indirect(atomloopA4, atomID1,atomID2)
            else if( atomloopA4 == atomID2) then
                KdA(ipar1, i1) = BROKENMASKDIHERAL
                call bond_breaking_event_dihedral_print(atomID1, atomID1, atomloopA2, atomloopA3, atomloopA4)
                ! need to do atomloopA2
                call bond_breaking_event_dihedral_indirect(atomloopA2, atomID1,atomID2)
                ! need to do atomloopA3
                call bond_breaking_event_dihedral_indirect(atomloopA3, atomID1,atomID2)
            endif
        enddo
  
        ! do atomID1, B series
        do i1=1, nB(ipar1)
            atomloopB1 = atom1B(ipar1,i1)+(atomID1-ipar1)
            atomloopB3 = atom3B(ipar1,i1)+(atomID1-ipar1)
            atomloopB4 = atom4B(ipar1,i1)+(atomID1-ipar1)
            if(atomloopB1 == atomID2) then
                KdB(ipar1, i1) = BROKENMASKDIHERAL
                call bond_breaking_event_dihedral_print(atomID1, atomloopB1, atomID1, atomloopB3, atomloopB4)
                ! need to do atomloopB3
                call bond_breaking_event_dihedral_indirect(atomloopB3, atomID1,atomID2)
                ! need to do atomloopB4
                call bond_breaking_event_dihedral_indirect(atomloopB4, atomID1,atomID2)
            else if( atomloopB3 == atomID2) then
                KdB(ipar1, i1) = BROKENMASKDIHERAL
                call bond_breaking_event_dihedral_print(atomID1, atomloopB1, atomID1, atomloopB3, atomloopB4)
                ! need to do atomloopB1
                call bond_breaking_event_dihedral_indirect(atomloopB1, atomID1,atomID2)
                ! need to do atomloopB4
                call bond_breaking_event_dihedral_indirect(atomloopB4, atomID1,atomID2)
            else if( atomloopB4 == atomID2) then
                KdB(ipar1, i1) = BROKENMASKDIHERAL
                call bond_breaking_event_dihedral_print(atomID1, atomloopB1, atomID1, atomloopB3, atomloopB4)
                ! need to do atomloopB1
                call bond_breaking_event_dihedral_indirect(atomloopB1, atomID1,atomID2)
                ! need to do atomloopB3
                call bond_breaking_event_dihedral_indirect(atomloopB3, atomID1,atomID2)
            endif
        enddo
        
    end subroutine
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to process bond breaking of dihedral indirect
!! \author Zhiye Tang
!<  
    subroutine bond_breaking_event_dihedral_indirect(atomID, atomID1, atomID2)
        use dihedral
        use param
        implicit none
        integer(4) :: atomID1, atomID2, atomID
        integer(4) :: ipar1
        integer(4) :: i1, atomloopA2, atomloopA3, atomloopA4
        integer(4) :: atomloopB1, atomloopB3, atomloopB4
  
        ipar1 = paranum(atomID)
  
        ! do atomID1, A series
        do i1=1, nA(ipar1)
            atomloopA2 = atom2A(ipar1,i1)+(atomID-ipar1)
            atomloopA3 = atom3A(ipar1,i1)+(atomID-ipar1)
            atomloopA4 = atom4A(ipar1,i1)+(atomID-ipar1)
            if( (atomloopA2 == atomID1 .and. atomloopA3 == atomID2) .or. &
      &         (atomloopA2 == atomID2 .and. atomloopA3 == atomID1) .or. &
      &         (atomloopA3 == atomID1 .and. atomloopA4 == atomID2) .or. &
      &         (atomloopA3 == atomID2 .and. atomloopA4 == atomID1) ) then
                KdA(ipar1, i1) = BROKENMASKDIHERAL
                call bond_breaking_event_dihedral_print(atomID, atomID, atomloopA2, atomloopA3, atomloopA4)
            endif
        enddo
  
        ! do atomID1, B series
        do i1=1, nB(ipar1)
            atomloopB1 = atom1B(ipar1,i1)+(atomID-ipar1)
            atomloopB3 = atom3B(ipar1,i1)+(atomID-ipar1)
            atomloopB4 = atom4B(ipar1,i1)+(atomID-ipar1)
            if( (atomloopB3 == atomID1 .and. atomloopB4 == atomID2) .or. &
      &         (atomloopB3 == atomID2 .and. atomloopB4 == atomID1) ) then
                KdB(ipar1, i1) = BROKENMASKDIHERAL
                call bond_breaking_event_dihedral_print(atomID, atomloopB1, atomID, atomloopB3, atomloopB4)
            endif
        enddo
        
    end subroutine
#endif /* DIHEDRAL_TABLE */
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to print bond breaking information of dihedral
!! \author Zhiye Tang
!<  
    subroutine bond_breaking_event_dihedral_print(atomID, atomID1, atomID2, atomID3, atomID4)
        use mpi_tool
        implicit none
        integer(4) :: atomID1, atomID2, atomID3, atomID4, atomID
  
#ifdef PRINTBROKEN
        if(myrank==0) write(*,'(A17I8A1I8I8I8I8)') "Dihedral    Atom[", atomID, "]", atomID1, atomID2, atomID3, atomID4
#endif
    end subroutine
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to process bond breaking of improper
!! \author Zhiye Tang
!<   
    subroutine bond_breaking_event_improper_torsion(atomID1, atomID2)
        implicit none
        integer(4) :: atomID1, atomID2
        ! atomID1 and atomID2 are global atom ID
  
        call bond_breaking_event_improper_torsion_direct(atomID1, atomID2)
        call bond_breaking_event_improper_torsion_direct(atomID2, atomID1)
  
    end subroutine
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to process bond breaking of improper direct
!! \author Zhiye Tang
!<  
    subroutine bond_breaking_event_improper_torsion_direct(atomID1, atomID2)
        use improper_torsion
        use param
        implicit none
        integer(4) :: atomID1, atomID2
        integer(4) :: ipar1
        integer(4) :: i1, atomloopA2, atomloopA3, atomloopA4
        integer(4) :: atomloopB1, atomloopB3, atomloopB4
  
        ipar1 = paranum(atomID1)
  
        ! do atomID1, A series
        do i1=1, nA(ipar1)
            atomloopA2 = atom2A(ipar1,i1)+(atomID1-ipar1)
            atomloopA3 = atom3A(ipar1,i1)+(atomID1-ipar1)
            atomloopA4 = atom4A(ipar1,i1)+(atomID1-ipar1)
            if(atomloopA2 == atomID2) then
                call bond_breaking_event_improper_torsion_maskA(ipar1, i1)
                call bond_breaking_event_improper_torsion_print(atomID1, atomID1, atomloopA2, atomloopA3, atomloopA4)
                ! need to do atomloopA3
                call bond_breaking_event_improper_torsion_indirect(atomloopA3, atomID1,atomID2)
                ! need to do atomloopA4
                call bond_breaking_event_improper_torsion_indirect(atomloopA4, atomID1,atomID2)
            else if( atomloopA3 == atomID2) then
                call bond_breaking_event_improper_torsion_maskA(ipar1, i1)
                call bond_breaking_event_improper_torsion_print(atomID1, atomID1, atomloopA2, atomloopA3, atomloopA4)
                ! need to do atomloopA2
                call bond_breaking_event_improper_torsion_indirect(atomloopA2, atomID1,atomID2)
                ! need to do atomloopA4
                call bond_breaking_event_improper_torsion_indirect(atomloopA4, atomID1,atomID2)
            else if( atomloopA4 == atomID2) then
                call bond_breaking_event_improper_torsion_maskA(ipar1, i1)
                call bond_breaking_event_improper_torsion_print(atomID1, atomID1, atomloopA2, atomloopA3, atomloopA4)
                ! need to do atomloopA2
                call bond_breaking_event_improper_torsion_indirect(atomloopA2, atomID1,atomID2)
                ! need to do atomloopA3
                call bond_breaking_event_improper_torsion_indirect(atomloopA3, atomID1,atomID2)
            endif
        enddo
  
        ! do atomID1, B series
        do i1=1, nB(ipar1)
            atomloopB1 = atom1B(ipar1,i1)+(atomID1-ipar1)
            atomloopB3 = atom3B(ipar1,i1)+(atomID1-ipar1)
            atomloopB4 = atom4B(ipar1,i1)+(atomID1-ipar1)
            if(atomloopB1 == atomID2) then
                call bond_breaking_event_improper_torsion_maskB(ipar1, i1)
                call bond_breaking_event_improper_torsion_print(atomID1, atomloopB1, atomID1, atomloopB3, atomloopB4)
                ! need to do atomloopB3
                call bond_breaking_event_improper_torsion_indirect(atomloopB3, atomID1,atomID2)
                ! need to do atomloopB4
                call bond_breaking_event_improper_torsion_indirect(atomloopB4, atomID1,atomID2)
            else if( atomloopB3 == atomID2) then
                call bond_breaking_event_improper_torsion_maskB(ipar1, i1)
                call bond_breaking_event_improper_torsion_print(atomID1, atomloopB1, atomID1, atomloopB3, atomloopB4)
                ! need to do atomloopB1
                call bond_breaking_event_improper_torsion_indirect(atomloopB1, atomID1,atomID2)
                ! need to do atomloopB4
                call bond_breaking_event_improper_torsion_indirect(atomloopB4, atomID1,atomID2)
            else if( atomloopB4 == atomID2) then
                call bond_breaking_event_improper_torsion_maskB(ipar1, i1)
                call bond_breaking_event_improper_torsion_print(atomID1, atomloopB1, atomID1, atomloopB3, atomloopB4)
                ! need to do atomloopB1
                call bond_breaking_event_improper_torsion_indirect(atomloopB1, atomID1,atomID2)
                ! need to do atomloopB3
                call bond_breaking_event_improper_torsion_indirect(atomloopB3, atomID1,atomID2)
            endif
        enddo
        
    end subroutine
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to process bond breaking of improper indirect
!! \author Zhiye Tang
!<  
    subroutine bond_breaking_event_improper_torsion_indirect(atomID, atomID1, atomID2)
        use improper_torsion
        use param
        implicit none
        integer(4) :: atomID1, atomID2, atomID
        integer(4) :: ipar1
        integer(4) :: i1, atomloopA2, atomloopA3, atomloopA4
        integer(4) :: atomloopB1, atomloopB3, atomloopB4
  
        ipar1 = paranum(atomID)
  
        ! do atomID1, A series
        do i1=1, nA(ipar1)
            atomloopA2 = atom2A(ipar1,i1)+(atomID-ipar1)
            atomloopA3 = atom3A(ipar1,i1)+(atomID-ipar1)
            atomloopA4 = atom4A(ipar1,i1)+(atomID-ipar1)
            if( (atomloopA2 == atomID1 .and. atomloopA3 == atomID2) .or. &
      &         (atomloopA2 == atomID2 .and. atomloopA3 == atomID1) .or. &
  
      &         (atomloopA3 == atomID1 .and. atomloopA4 == atomID2) .or. &
      &         (atomloopA3 == atomID2 .and. atomloopA4 == atomID1) .or. &
  
      &         (atomloopA2 == atomID1 .and. atomloopA4 == atomID2) .or. &
      &         (atomloopA2 == atomID2 .and. atomloopA4 == atomID1) ) then
                call bond_breaking_event_improper_torsion_maskA(ipar1, i1)
                call bond_breaking_event_improper_torsion_print(atomID, atomID, atomloopA2, atomloopA3, atomloopA4)
            endif
        enddo
  
        ! do atomID1, B series
        do i1=1, nB(ipar1)
            atomloopB1 = atom1B(ipar1,i1)+(atomID-ipar1)
            atomloopB3 = atom3B(ipar1,i1)+(atomID-ipar1)
            atomloopB4 = atom4B(ipar1,i1)+(atomID-ipar1)
            if( (atomloopB1 == atomID1 .and. atomloopB3 == atomID2) .or. &
      &         (atomloopB1 == atomID2 .and. atomloopB3 == atomID1) .or. &
  
      &         (atomloopB3 == atomID1 .and. atomloopB4 == atomID2) .or. &
      &         (atomloopB3 == atomID2 .and. atomloopB4 == atomID1) .or. &
  
      &         (atomloopB1 == atomID1 .and. atomloopB4 == atomID2) .or. &
      &         (atomloopB1 == atomID2 .and. atomloopB4 == atomID1) ) then
                call bond_breaking_event_improper_torsion_maskB(ipar1, i1)
                call bond_breaking_event_improper_torsion_print(atomID, atomloopB1, atomID, atomloopB3, atomloopB4)
            endif
        enddo
        
    end subroutine
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to process bond breaking of dihedral, mask force constant
!! \author Zhiye Tang
!<  
    subroutine bond_breaking_event_improper_torsion_maskA(ipar1, i1)
        use improper_torsion
        implicit none
        integer(4) :: ipar1, i1
#if defined(GROEXT_OPLS) || defined(GAFF)
        KdA(ipar1, i1) = BROKENMASK
#else
        KitA(ipar1, i1) = BROKENMASK
#endif
    end subroutine
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to process bond breaking of dihedral, mask force constant
!! \author Zhiye Tang
!<
    subroutine bond_breaking_event_improper_torsion_maskB(ipar1, i1)
        use improper_torsion
        implicit none
        integer(4) :: ipar1, i1
#if defined(GROEXT_OPLS) || defined(GAFF)
        KdB(ipar1, i1) = BROKENMASK
#else
        KitB(ipar1, i1) = BROKENMASK
#endif
    end subroutine
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to print bond breaking information of improper
!! \author Zhiye Tang
!<
    subroutine bond_breaking_event_improper_torsion_print(atomID, atomID1, atomID2, atomID3, atomID4)
        use mpi_tool
        implicit none
        integer(4) :: atomID1, atomID2, atomID3, atomID4, atomID
  
#ifdef PRINTBROKEN
        if(myrank==0) write(*,'(A17,I8,A1,I8,I8,I8,I8)') "Improper    Atom[", atomID, "]", atomID1, atomID2, atomID3, atomID4
#endif
  
    end subroutine
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to print bond breaking information
!! \author Zhiye Tang
!<
    subroutine bond_breaking_event_print(atomID1, atomID2)
        use md_condition, only : mdstep
        use file_mdbb, only : write_mdbb_entry
        implicit none
        integer(4) :: atomID1, atomID2
  
        write(*,'(A3,I10,A24,I8,A5,I8,A9)') "At ", mdstep, " bond breaking between ", atomID1, " and ", atomID2, " (1-base)"
  
        call write_mdbb_entry(mdstep, atomID1, atomID2)
    end subroutine
  
  end module
  
  
