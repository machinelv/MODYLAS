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
!! \brief  Module and subroutines to constrain atom position.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to constran atom position.
!! \author  Atsushi Yamada
!<
module position_constrain
  use omp_lib
  use subcell
  implicit none
  !     type_p_constrain = 0 (no position constrain) / 1 (harmonic) /
  !                        2 (huge mass) / 3 (fix (for optimize))
  !               (*) "huge_mass" cannot work well fixing atomic coordinates
  !                    in NPT
  !               (*) "huge_mass" cannot work if you specified individual
  !                    molecules among the same molecular species
  !     atom_type_p_constrain = 0 (all types of atoms are constrained)
  !                             1 (heavy atoms are constrained)
  !     molecules_p_constrain = 0 (all molecules are listed)
  !                             1 (listed molecules are specified)
  !     posi_type_p_constrain = 0 (initial coordinates in mdxyz are used)
  !                             1 (the coordinates are specified)
  !     npconstrain : the number of atoms to be constrained
  integer(4) :: type_p_constrain=0
  integer(4) :: atom_type_p_constrain, molecules_p_constrain
  integer(4) :: posi_type_p_constrain
  integer(4) :: npconstrain=0, natom_allow_uncons=0
  integer(4) :: howmany_mol_p_constrain
  integer(4),allocatable :: atom_p_constrain(:), mol_p_constrain(:)
  integer(4),allocatable :: atoms_allow_uncons(:)
  integer(4),allocatable :: one_or_zero_fix_atom(:)
  real(8) :: K_p_constrain
  real(8),allocatable :: ref_crd0_p_constrain(:,:)

contains

!----------------------------------------------------------------------
!  Initial procedures for position constrain
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to initialize for position constrain.
!! \author  Atsushi Yamada
!<
  subroutine init_position_constrain
    use atom_mass
    use device_numbers
    use session_name_mod
    use trajectory_org
    use md_const
    use ensemble_numbers
    use md_condition
    use molecules
    use param
    use mpi_tool
    implicit none
    logical :: flg_allow_unconst
    integer(4) :: i,j,k,l, im, jm, num, icnt, ioerr
    real(8) ::  mass_heavy, mass_huge
    character(1024) :: line
    character(30) :: key
    integer(4) :: ierr
    include 'mpif.h'

    ! === In the case that posicion constrain is specified in mdff ===
    if( type_p_constrain .ne. 0 ) then

       ! -- log --
       if(myrank==0) then
          if(     type_p_constrain == 0 ) then
!            write(*,*) "No position constrain"
          else if(type_p_constrain == 1 ) then
             write(*,*) "Harmonic positiotn constrain"
             write(*,*) "Harmonic force constant= ", K_p_constrain, " [J/m^2]"
          else if(type_p_constrain == 2 ) then
             if(md_condition__ensemble == NPT_A .or. md_condition__ensemble == NPT_PR ) then
                write(0,9300)
                call modylas_abort()
             else if(md_condition__ensemble == OPT) then
                write(0,*) "ERROR: You can not use huge_mass option in opt"
                call modylas_abort()
             endif
             write(*,*) "Huge Mass positiotn constrain"
             write(*,9200)
          else if(type_p_constrain == 3 ) then
             write(*,*) "Fix positiotn constrain in opt"
          else
             write(0,*) "ERROR: Input Error for Position Constrain"
             call modylas_abort()
          endif
       endif

       ! -- allocate --
       if( md_condition__ensemble==OPT .or. type_p_constrain==1) then
          allocate( one_or_zero_fix_atom(n) )
          one_or_zero_fix_atom(:) = 1
       endif

       mass_heavy = 1.5d0 * md_ATOMIC_MASS_UNIT
       mass_huge  = 9.99999d99

       ! -- Check the number of atoms to be constrained --
       num = 0
       do i= 1,n
          if(      molecules_p_constrain==0) then        ! all molecules

             if(      atom_type_p_constrain==0) then     !   all types of atoms
                call check_atom_allow_unconstrain(i,flg_allow_unconst)
                if(.not.flg_allow_unconst) num = num + 1
             else if( atom_type_p_constrain==1) then     !   heavy atoms
                if( mass(paranum(i)) .gt. mass_heavy ) then
                   call check_atom_allow_unconstrain(i,flg_allow_unconst)
                   if(.not.flg_allow_unconst) num=num+1
                endif
             endif

          else if( molecules_p_constrain==1) then        ! specified molecules

             im = atom2mol(i)
             do  jm= 1,howmany_mol_p_constrain

                if(mol_p_constrain(jm)==im) then

                   if(      atom_type_p_constrain==0) then     !   all types of atoms
                      call check_atom_allow_unconstrain(i,flg_allow_unconst)
                      if(.not.flg_allow_unconst) num = num + 1
                   else if( atom_type_p_constrain==1) then     !   heavy atoms
                      if( mass(paranum(i)) .gt. mass_heavy ) then
                         call check_atom_allow_unconstrain(i,flg_allow_unconst)
                         if(.not.flg_allow_unconst) num = num + 1
                      endif
                   endif

                endif
             enddo

          endif
       enddo

       npconstrain = num
       if(myrank==0) write(*,9100) npconstrain  !(log)

       ! allocate
       allocate( atom_p_constrain(npconstrain) )

       ! Make atoms list to be constrained
       k = 0
       do i= 1,n
          if(      molecules_p_constrain==0) then        ! all molecules
             if(      atom_type_p_constrain==0) then     !   all types of atoms
                call check_atom_allow_unconstrain(i,flg_allow_unconst)
                if(.not.flg_allow_unconst) then
                   k = k + 1
                   atom_p_constrain(k) = i
                endif
             else if( atom_type_p_constrain==1) then     !   heavy atoms
                if( mass(paranum(i)) .gt. mass_heavy ) then
                   call check_atom_allow_unconstrain(i,flg_allow_unconst)
                   if(.not.flg_allow_unconst) then
                      k = k + 1
                      atom_p_constrain(k) = i
                   endif
                endif
             endif

          else if( molecules_p_constrain==1) then        ! specified molecules

             im = atom2mol(i)
             do  jm= 1,howmany_mol_p_constrain
                if(mol_p_constrain(jm)==im) then

                   if(      atom_type_p_constrain==0) then  !   all types of atoms
                      call check_atom_allow_unconstrain(i,flg_allow_unconst)
                      if(.not.flg_allow_unconst) then
                         k = k + 1
                         atom_p_constrain(k) = i
                      endif
                   else if( atom_type_p_constrain==1) then  !   heavy atoms
                      if( mass(paranum(i)) .gt. mass_heavy ) then
                         call check_atom_allow_unconstrain(i,flg_allow_unconst)
                         if(.not.flg_allow_unconst) then
                            k = k + 1
                            atom_p_constrain(k) = i
                         endif
                      endif
                   endif

                endif
             enddo

          endif
       enddo

       ! Set parameters
       do k= 1,npconstrain

          i = atom_p_constrain(k)

          ! Set mass to huge and set velocity to zero
          ! (the velocity of the center of mass can not be zero)
          ! NOTE1: if conbining SHAKE, this program does not run
          ! NOTE2: huge_mass option can not be applied to a specified molecule.
          !        this option changes the mass of all molecules included in the species
          if( type_p_constrain==2 ) then        ! huge mass
             mass(paranum(i))   = mass_huge
             r_mass(paranum(i)) = 1.0d0/mass_huge
             v(:,i) = 0.0d0
          else if( type_p_constrain==3 .or. type_p_constrain==1) then
             ! fix for optimize and harmonic
             one_or_zero_fix_atom(i) = 0
          endif

       enddo


       ! === In the case that posicion constrain is specified
       !                          by sessionname.posiconst input file ===
    else if(npconstrain==0 .and. type_p_constrain==0) then


       !! Now, this version read extra input file, "posi_const.dat"
       !! which gives the atoms list being constrained
       !! The format of the file posi_const.dat
       !! line 1  : keyword
       !!         : "no-psition-constrain",
       !!         : "huge_mass", "harmonic", "fix"(for opt)
       !!         : In the case of "harmonic",
       !!           write force constant in the same line
       !!           e.x.   harmonic   1.0
       !! line 2- : 0 or 1  (0: to be constrained, 1: free)
       !!           for all atoms
       !!           in the case of harmonic, write coordinates parameters
       !!           e.x.
       !!           0   1.00  -12.123  2.234
       !!           0   0.10    2.345  7.999
       !!           1   0.30    1.345  9.999
       !!           0   0.20    3.345  2.999
       !!           1   23.00   43.90 -23.00
       !!           ....  (total n)

9000   format(a)

    if(myrank==0) then

       open(f_pconst,file=trim(session_name)//'.posiconst', status="old",iostat=ioerr)

       ! --support old version : file name = posi_const.dat--
       if(ioerr .ne. 0) then
          open(f_pconst, file="posi_const.dat", status="old",err=100)
       endif

       ! -- read keyword (fist line) --
       read(f_pconst,9000) line
       read(line,*) key
       write(*,*) "Input Keyword= ", key
       if(      key == "no-position-constrain" ) then
          type_p_constrain = 0
!         write(*,*)"No position constrain"
       else if( key == "harmonic" ) then
          type_p_constrain = 1
          write(*,*)"Harmonic positiotn constrain"
          read(line,*) key, K_p_constrain
          !convert unit : [kcal /mol A^2] => [J / m^2]
          K_p_constrain = K_p_constrain*1000.0d0*md_CALORIE/md_AVOGADRO*1.0d20
          write(*,*) "Harmonic force constant= ", K_p_constrain, " [J/m^2]"
       else if( key == "huge_mass" ) then
          type_p_constrain = 2
          if(md_condition__ensemble == NPT_A .or. md_condition__ensemble == NPT_PR ) then
             write(0,9300)
             call modylas_abort()
          else if(md_condition__ensemble == OPT)then
             write(0,*) "ERROR: You can not use huge_mass option in opt"
             call modylas_abort()
          endif
          write(*,*)"Huge Mass positiotn constrain"
          write(*,9200)
       else if( key == "fix" ) then
          type_p_constrain = 3
          write(*,*)"Fix positiotn constrain in opt"
       else
          write(0,*) "ERROR: Input Error for Position Constrain"
          call modylas_abort()
       endif

    endif  ! myrank=0

    call MPI_Bcast(type_p_constrain,1,MPI_INTEGER4,0, MPI_COMM_WORLD, ierr)
    if(type_p_constrain == 1) then
       call MPI_Bcast(K_p_constrain,1,MPI_REAL8,0, MPI_COMM_WORLD, ierr)
    endif

    mass_heavy = 1.5d0 * md_ATOMIC_MASS_UNIT
    mass_huge  = 9.99999d99

    !!(opt needs to allocate one_or_zero_fix_atom
    !!     even for no position constrain )
    if( md_condition__ensemble==OPT .or. type_p_constrain==1) then
       allocate( one_or_zero_fix_atom(n) )
       one_or_zero_fix_atom(:) = 1
    endif

    if( type_p_constrain == 0 ) return

    if(myrank==0) then

       icnt=0
       do i=1,n
          read(f_pconst,*) j
          if(j==0) icnt = icnt + 1
          if(md_condition__ensemble==OPT.or.type_p_constrain==1)then
             one_or_zero_fix_atom(i) = j
          endif
       enddo
       npconstrain = icnt

    endif

    call MPI_Bcast(npconstrain,1,MPI_INTEGER4,0, MPI_COMM_WORLD, ierr)
    if( md_condition__ensemble==OPT .or. type_p_constrain==1) then
       call MPI_Bcast(one_or_zero_fix_atom,n,MPI_INTEGER4,0, MPI_COMM_WORLD, ierr)
    endif


    ! allocate
    allocate( atom_p_constrain(npconstrain) )
    allocate( ref_crd0_p_constrain(3,n) )


    if(myrank==0) then
       rewind f_pconst
       read(f_pconst,9000) line
    endif

    k = 0
    do i= 1,n
       if(myrank==0) then
          if( key == "harmonic" ) then
             read(f_pconst,*) j,(ref_crd0_p_constrain(l,i),l=1,3)
             ref_crd0_p_constrain(:,i) = ref_crd0_p_constrain(:,i) * 1.0d-10
          else
             read(f_pconst,*) j
          endif
       endif
       call MPI_Bcast(j,1,MPI_INTEGER4,0,MPI_COMM_WORLD, ierr)

       if(j==0) then
          k = k + 1
          atom_p_constrain(k) = i

          ! huge mass
          if( type_p_constrain==2 ) then
             mass(paranum(i))   = mass_huge
             r_mass(paranum(i)) = 1.0d0/mass_huge
             v(:,i) = 0.0d0
          endif
       endif
    enddo

    call MPI_Bcast(ref_crd0_p_constrain,3*n,MPI_REAL8,0, MPI_COMM_WORLD, ierr)


    if(myrank==0) then
       close(f_pconst)
       write(*,9100) npconstrain  !(log)
    endif

 endif

 return

   ! in the case of no sessionname.posiconst or posi_const.dat file
100 continue
 type_p_constrain = 0
!write(*,*)"No position constrain"

   call MPI_Bcast(type_p_constrain,1,MPI_INTEGER4,0, MPI_COMM_WORLD, ierr)

9100 format("The total number of atoms to be constrained :",i8)
9200 format("WARNING of Huge Mass option:", / &
          &  "  This option can not be applied to a specified molecule.", / &
          &  "  This changes the atomic masses in all molecules ", &
          &  "  among the same species", / &
          &  "  Check it by yourself" )
9300 format("ERROR of Huge Mass option:", / &
          &       "   This does not work well in NPT ensemble")
 end subroutine init_position_constrain

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to check if the atom is unconstrained.
!! \author  Atsushi Yamada
!<
  subroutine  check_atom_allow_unconstrain(k,flg_allow_unconst)
    implicit none
    logical ::  flg_allow_unconst
    integer(4) :: i,k

    flg_allow_unconst = .false.
    do i=1,natom_allow_uncons
       if(k==atoms_allow_uncons(i))then
          flg_allow_unconst = .true.
          exit
       endif
    enddo
  end subroutine check_atom_allow_unconstrain

!======================================================================
!>
!! \brief  Subroutine to calculate force and energy for position constrain
!! \author  Atsushi Yamada
!<
  subroutine add_position_constrain_forces()
    !
    !      V = K*r^2   (r=|r-r0|)
    !
    use atom_mass
    use md_const
    use forces
    use md_monitors
    use trajectory_org
    use trajectory_mpi
    use atom_virial
    implicit none
    integer(4) :: i, k0,i0, iam
    real(8) :: r, r2, Fr__x, Fr__y, Fr__z, dr__x, dr__y, dr__z
    real(8) :: U_p_constrain

    if(type_p_constrain.ne.1) return

    iam = 0
    U_p_constrain = 0.0d0
!$omp parallel default(shared) &
!$omp& private(k0,i0,i,iam) &
!$omp& private(dr__x,dr__y,dr__z,r2,r,Fr__x,Fr__y,Fr__z) &
!$omp& reduction(+:U_p_constrain)
!$    iam = omp_get_thread_num()
!$omp do
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          i = m2i(i0)
          dr__x = wkxyz(1,i0) - ref_crd0_p_constrain(1,i)
          dr__y = wkxyz(2,i0) - ref_crd0_p_constrain(2,i)
          dr__z = wkxyz(3,i0) - ref_crd0_p_constrain(3,i)
          dr__x = dr__x * dble(1-one_or_zero_fix_atom(i))
          dr__y = dr__y * dble(1-one_or_zero_fix_atom(i))
          dr__z = dr__z * dble(1-one_or_zero_fix_atom(i))
          r2    = dr__x*dr__x + dr__y*dr__y + dr__z*dr__z
          r     = dsqrt(r2)
          Fr__x = -2.0d0 * K_p_constrain * dr__x
          Fr__y = -2.0d0 * K_p_constrain * dr__y
          Fr__z = -2.0d0 * K_p_constrain * dr__z
          w3_f(1,i0,iam) = w3_f(1,i0,iam) + Fr__x
          w3_f(2,i0,iam) = w3_f(2,i0,iam) + Fr__y
          w3_f(3,i0,iam) = w3_f(3,i0,iam) + Fr__z
          U_p_constrain = U_p_constrain + K_p_constrain * r2
          !! virial .... now ignored
       enddo ! i0
    enddo ! k0
!$omp end do
!$omp end parallel

    wk_p_energy = wk_p_energy + U_p_constrain
  end subroutine add_position_constrain_forces
!======================================================================
!>
!! \brief  Subroutine to set atomic coordinate to reference position for position constrain.
!! \author  Atsushi Yamada
!<
  subroutine position_constrain_fix_atom()
!
!   This is for NPT because atomic coordinates change by changing cell size
!   For "huge_mass" keyword in position constrain
    use atom_mass
    use md_const
    use md_forces
    use md_monitors
    use trajectory_org
    use trajectory_mpi
    use param
    use atom_virial
    implicit none
    integer(4) :: i,i0,k0

    if(type_p_constrain.ne.2) return  !! huge_mass keyword

    !! NEED to check later
    write(*,*) "check program: position_constrain_fix_atom"

    !! Now, this subroutine does not work (->core dumped) ay 2012.02.09
!$omp parallel do default(shared) &
!$omp& private(k0,i0,i)
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          i = m2i(i0)
          if(mass(paranum(i)) .ge. 9.9d99) then
             wkxyz(1,i0) = ref_crd0_p_constrain(1,i)
             wkxyz(2,i0) = ref_crd0_p_constrain(2,i)
             wkxyz(3,i0) = ref_crd0_p_constrain(3,i)
          endif
       enddo
    enddo
  end subroutine position_constrain_fix_atom

!>
!! \brief Subroutine to set variables of /position constrain/type. 
!! \author  Atsushi Yamada
!<
  subroutine fmod_p_cons__type(ival)
    implicit none
    integer(4), intent(in) :: ival
    type_p_constrain = ival
  end subroutine fmod_p_cons__type
!---------------------------------------------------------------------
!>
!! \brief Subroutine to set variables of /position constrain/force_constant.
!! \author  Atsushi Yamada
!<
  subroutine fmod_p_cons__fc(val)
    implicit none
    real(8), intent(in) :: val
    K_p_constrain = val
  end subroutine fmod_p_cons__fc
!---------------------------------------------------------------------
!>
!! \brief Subroutine to set variables of /position constrain/atom/atom_type.
!! \author  Atsushi Yamada
!<
  subroutine fmod_p_cons_atom__atom_type(ival)
    implicit none
    integer(4), intent(in) :: ival
    atom_type_p_constrain = ival
  end subroutine fmod_p_cons_atom__atom_type
!---------------------------------------------------------------------
!>
!! \brief Subroutine to set variables of /position constrain/atom/molecules.
!! \author  Atsushi Yamada
!<
  subroutine fmod_p_cons_atom__molecules(ival)
    implicit none
    integer(4), intent(in) :: ival
    molecules_p_constrain = ival
  end subroutine fmod_p_cons_atom__molecules
!---------------------------------------------------------------------
!>
!! \brief Subroutine to set variables of /position constrain/atom/nmolecule.
!! \author  Atsushi Yamada
!<
  subroutine fmod_p_cons_atom__nmolecule(ival)
    implicit none
    integer(4), intent(in) :: ival
    howmany_mol_p_constrain = ival
  end subroutine fmod_p_cons_atom__nmolecule
!---------------------------------------------------------------------
!>
!! \brief Subroutine to allocate arrays for position constrain.
!! \author  Atsushi Yamada
!<
  subroutine fmod_alloc_mol_p_constrain(ival)
    implicit none
    integer(4), intent(in) :: ival
    allocate( mol_p_constrain(ival) )
  end subroutine fmod_alloc_mol_p_constrain
!---------------------------------------------------------------------
!>
!! \brief Subroutine to set variables for /position constrain/atom/nmolecule.
!! \author  Atsushi Yamada
!<
  subroutine add_to_mol_p_constrain(i,im)
    implicit none
    integer(4), intent(in) :: i,im
    mol_p_constrain(i) = im
  end subroutine add_to_mol_p_constrain
!---------------------------------------------------------------------
!>
!! \brief  Subroutine to set variables of /position constrain/atom/natom_allow_unconstrain.
!! \author  Atsushi Yamada
!<
  subroutine fmod_p_cons_atom__natm_alw_ucn(ival)
    implicit none
    integer(4), intent(in) :: ival
    natom_allow_uncons = ival
  end subroutine fmod_p_cons_atom__natm_alw_ucn
!---------------------------------------------------------------------
!>
!! \brief  Subroutine to allocate arrays of /position constrain/atom/natom_allow_unconstrain.
!! \author  Atsushi Yamada
!<
  subroutine fmod_alloc_natom_allow_uncons(ival)
    implicit none
    integer(4), intent(in) :: ival
    allocate( atoms_allow_uncons(ival) )
  end subroutine fmod_alloc_natom_allow_uncons
!---------------------------------------------------------------------
!>
!! \brief  Subroutine to set variables of /position constrain/atom/allow atom unconstrain.
!! \author  Atsushi Yamada
!<
  subroutine add_to_atom_allow_uncons(i,ival)
    implicit none
    integer(4), intent(in) :: i,ival
    atoms_allow_uncons(i) = ival
  end subroutine add_to_atom_allow_uncons
!---------------------------------------------------------------------
!>
!! \brief  Subroutine to set variables of /position constrain/position/type.
!! \author  Atsushi Yamada
!<
  subroutine fmod_p_cons_position__type(ival)
    implicit none
    integer(4), intent(in) :: ival
    posi_type_p_constrain = ival
  end subroutine fmod_p_cons_position__type
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to allocate arrays of position constrain/position/coordinate.
!! \author  Atsushi Yamada
!<
  subroutine fmod_alloc_ref_crd0_p_cons
    use trajectory_org
    implicit none
    allocate(ref_crd0_p_constrain(3,n))
    ref_crd0_p_constrain(:,:) = 9.999d99
  end subroutine fmod_alloc_ref_crd0_p_cons
!---------------------------------------------------------------------
!>
!! \brief  Subroutine to set variables of position constrain/position/coordinate.
!! \author  Atsushi Yamada
!<
  subroutine fmod_p_cons_position__ref_crd0(ival,v1,v2,v3)
    implicit none
    integer(4), intent(in) :: ival
    real(8), intent(in) :: v1, v2, v3
    ref_crd0_p_constrain(1,ival) = v1
    ref_crd0_p_constrain(2,ival) = v2
    ref_crd0_p_constrain(3,ival) = v3
  end subroutine fmod_p_cons_position__ref_crd0
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine for reading position constrain part in mdff.bin file
!! \author  Atsushi Yamada
!<
!-------------------------------------------------------------------------
  subroutine read_position_constrain
    use trajectory_org
    use device_numbers, only : f_mdff
    use mpi_tool
    implicit none
    include 'mpif.h'
    integer(4) :: ierr

    if(myrank==0) read(f_mdff) type_p_constrain
    call MPI_Bcast(type_p_constrain,1,MPI_INTEGER4,0,MPI_COMM_WORLD, ierr)
    if(type_p_constrain.eq.0) return
    if(type_p_constrain.eq.1)  then
       if(myrank==0) then
          read(f_mdff) K_p_constrain
       endif
       call MPI_Bcast(K_p_constrain,1,MPI_REAL8,0,MPI_COMM_WORLD, ierr)
    endif
    if(myrank==0) then
       read(f_mdff) atom_type_p_constrain
       read(f_mdff) molecules_p_constrain
       read(f_mdff) howmany_mol_p_constrain
    endif
    call MPI_Bcast(atom_type_p_constrain,1,MPI_INTEGER4,0,MPI_COMM_WORLD, ierr)
    call MPI_Bcast(molecules_p_constrain,1,MPI_INTEGER4,0,MPI_COMM_WORLD, ierr)
    call MPI_Bcast(howmany_mol_p_constrain,1,MPI_INTEGER4,0,MPI_COMM_WORLD, ierr)

    allocate( mol_p_constrain(howmany_mol_p_constrain) )

    if(molecules_p_constrain.eq.1) then
       if(myrank==0) then
          read(f_mdff) mol_p_constrain
       endif
       call MPI_Bcast(mol_p_constrain,howmany_mol_p_constrain, &
            &  MPI_INTEGER4,0,MPI_COMM_WORLD, ierr)
    endif
    if(myrank==0) then
       read(f_mdff) natom_allow_uncons
    endif
    call MPI_Bcast(natom_allow_uncons,1,MPI_INTEGER4,0,MPI_COMM_WORLD, ierr)
    allocate( atoms_allow_uncons(natom_allow_uncons) )

    if(natom_allow_uncons.ne.0) then
       if(myrank==0) read(f_mdff) atoms_allow_uncons
       call MPI_Bcast(atoms_allow_uncons,natom_allow_uncons, &
            &  MPI_INTEGER4,0,MPI_COMM_WORLD, ierr)
    endif
    if(myrank==0) read(f_mdff) posi_type_p_constrain
    call MPI_Bcast(posi_type_p_constrain,1, &
         &               MPI_INTEGER4,0,MPI_COMM_WORLD, ierr)
    if(type_p_constrain==1) then
       allocate( ref_crd0_p_constrain(3,n) )
       if(myrank==0) read(f_mdff) ref_crd0_p_constrain
       call MPI_Bcast(ref_crd0_p_constrain,3*n, &
            &                   MPI_REAL8,0,MPI_COMM_WORLD, ierr)
    endif

    return
  end subroutine read_position_constrain
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine for writing position constrain part in mdff.bin file
!! \author  Atsushi Yamada
!<
!-------------------------------------------------------------------------
  subroutine write_position_constrain
    use device_numbers, only : f_mdff
    implicit none
    write(f_mdff) type_p_constrain
    if(type_p_constrain.eq.0) return
    if(type_p_constrain.eq.1)  then
       write(f_mdff) K_p_constrain
    endif
    write(f_mdff) atom_type_p_constrain
    write(f_mdff) molecules_p_constrain
    write(f_mdff) howmany_mol_p_constrain
    if(molecules_p_constrain.eq.1) then
       write(f_mdff) mol_p_constrain
    endif
    write(f_mdff) natom_allow_uncons
    if(natom_allow_uncons.ne.0) then
       write(f_mdff) atoms_allow_uncons
    endif
    write(f_mdff) posi_type_p_constrain
    if(type_p_constrain==1) then
       write(f_mdff) ref_crd0_p_constrain
    endif
  end subroutine write_position_constrain
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine for writing position constrain part
!! \author  Atsushi Yamada
!<
  subroutine write_memory_position_constrain
    implicit none
    write(*,*) '[write_position_constrain]'
    write(*,*) type_p_constrain
    if(type_p_constrain.eq.0) return
    if(type_p_constrain.eq.1)  then
       write(*,*) K_p_constrain
    endif
    write(*,*) atom_type_p_constrain
    write(*,*) molecules_p_constrain
    write(*,*) howmany_mol_p_constrain
    if(molecules_p_constrain.eq.1) then
       write(*,*) mol_p_constrain
    endif
    write(*,*) natom_allow_uncons
    if(natom_allow_uncons.ne.0) then
       write(*,*) atoms_allow_uncons
    endif
    write(*,*) posi_type_p_constrain
    if(type_p_constrain==1) then
       write(*,*) ref_crd0_p_constrain
    endif
   end subroutine write_memory_position_constrain

end module position_constrain
