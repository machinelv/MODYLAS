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
!! \brief  Module and subroutines to read inputs.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to read inputs (.mddef, .mdff, .mdxyz).
!! \author Kensuke Iwahashi
!<
!----------------------------------------------------------------------
module input_mod
  implicit none

contains

!>
!! \brief  Subroutine to read input files.
!! \author Kensuke Iwahashi
!<
  subroutine parse_input(read_mddef)
    use version
    use parse
    use device_numbers
    use session_name_mod
    use mpi_tool 
#ifdef GROEXT
    use file_groext
#endif
   implicit none
    logical :: read_mddef, lerr
    character(LEN=10) :: chara
    ! get input version number from mddef
    call parse_open(trim(session_name)//'.mddef', lerr)
    call parse_get('/input/version', chara, lerr)

    input_version = chara

    if(lerr) then ! NO keyword of input version: default=older version
       input_version = "0.9.0"
       if(myrank==0) then
          write(*,*) "No input version number"
          write(*,*) " --> Oldest input version 0.9.0 was used"
       endif
       lerr=.false.
    endif
#ifdef GROEXT
    call parse_groext()
#endif
    call parse_close(lerr)

!   if(trim(MODYLAS_version).ne.trim(input_version)) then
!      if(myrank==0)  write(*,*) "Warning: input version differs from MODYLAS version"
!   endif

    ! read input files depending on input version
    if(     trim(input_version) == '0.9.0') then
       call  parse_input_v0_9_0(read_mddef)
    else if(trim(input_version) == '1.0.0') then
       call  parse_input_v1_0_0(read_mddef)
    else
       if(myrank==0) write(0,*)"ERROR: wrong version number in mddef"
       call modylas_abort
    endif


  end subroutine parse_input
!------------------------------------------------------------------------
!     Read input for ver. 0.9.0 (or older version)
!------------------------------------------------------------------------
!>
!! \brief  Subroutine to read input files in version 0.9.0 or older.
!! \author Kensuke Iwahashi
!<
  subroutine parse_input_v0_9_0(read_mddef)
    use parse
    use file_mdmntr
    use file_mdtrj
    use file_restart
    use file_dcd
    use file_xtc
    use device_numbers
    use session_name_mod
    use file_utility
    use md_const
    use ensemble_numbers
    use md_condition
    use mpi_3d_grid, only : set_npx, set_npy, set_npz, mpi_manual_division_flg
    use md_multiplestep
    use thermostat
    use param
    use mol_info
    use md_periodic
    use center_of_mass_variables
    use shake_rattle_roll
    use mpi_tool
    use bond
    use bond_morse
    !use bond_pcff
    use angle
    !use angle_pcff
    !use bond_bond
    !use parse_angle_pcff
    use UB
    use dihedral
    use improper_torsion
    use CMAP
    use void123
    use special14
    use opt_mod
    use position_constrain
    use comm_pme_mod
    use pme_far
    use ewald_variables
    use ewald_mod
    use unit_cell
    use segments
    use trajectory_org
    use extended_system
    use coulomb_mod, only : fmod_set_mol_chgv, read_chgv
    use lj_mod, only : fmod_set_mol_ljs__epsilon_sqrt, fmod_set_mol_ljs__r_half, &
         & read_epsilon_sqrt, read_r_half
    use molecules, only : fmod_alloc_atom2mol, fmod_alloc_mol2atom
    use atom_mass, only : fmod_set_mol_mass, read_mass
    use fmm_parameters, only : fmod_set_nmax
    use fmm_far, only : fmod_set_lgflg, fmod_set_fmm_nlevel
    use parse_shake, only : fmod_alloc_yyparse_shake
    use cutoff_radius, only : fmod_cutrad, fmod_cutrad2
    use force_field_numbers
    use random_seed
    use maxwell_distribution, only : fmod_maxwell_temperature, fmod_reset_maxwell
    use file_force, only : fmod_force_start, fmod_force_interval
    use file_mdff, only : read_forcefield, read_md_shake, read_molecule, read_parameter_number
    use file_mdxyz_bin, only : read_mdxyzbin
    use md_oplsaa_special_divide, only : fmod_md_oplsaa_lj_sp_divide, fmod_md_oplsaa_cl_sp_divide
    use system_dipole, only: set_up_system_dipole
#ifdef TIP4
    use tip4p, only : fmod_set_tip4p_geometry, fmod_set_msite_position
#endif
    use trajectory_mpi
    use center_of_mass_variables
    use domain
    implicit none
    integer(4) :: i, j, k, io, im
    integer(4) :: natom, nthermostat, nbarostat
    integer(4) :: ivalue1, ivalue2, ivalue3
    integer(4) :: ivalue4, ivalue5, ivalue6
    real(8)    :: value1, value2, value3, value4, value5
    real(8)    :: value6, value7, value8, value9
    real(8)   :: qvalue1, qvalue2, qvalue3, qvalue4
    character(LEN=100) :: getname,segname
    character(LEN=10) :: chara
    logical :: read_mddef
    logical :: lerr=.false., lerr1, lerr2, lerr3, lerr4, lerr5
    logical :: lerr6, lerr7, lerr8, lerr9, lerr_bk
    logical :: velocity_scaling_flag = .false.
    logical :: volume_scaling_flag = .false.
    logical :: cellgeometry_scaling_flag = .false.
    logical :: read_velocity
    integer(4), allocatable :: natommol(:), nmolmol(:), molnseg(:)
    integer(4) :: hws, sid, nas, nall, temp, nvc, nvl
    integer(4) :: nsl, nsct, nbond, nag, nub, nit, ncm, ndi
    integer(4) :: dummyInt


    call  set_default_keyword_0_9_0

    !
    !     mdxyz
    !
    if (can_read_binary_file(trim(session_name)//'.mdxyz.bin')) then
       if(myrank==0) write(*,*) 'Reading .mdxyz.bin file started.'
       call read_mdxyzbin()
       if(myrank==0) write(*,*) 'Reading .mdxyz.bin file ended successfully!'
    else
       if(myrank==0) write(*,*) 'Reading .mdxyz file started.'
       call parse_ignore_key('/atom/positions')
       call parse_ignore_key('/atom/velocities')
       call parse_set_debug(0)
       call parse_open(trim(session_name)//'.mdxyz', lerr)
       if (lerr) then
          write(0,*) parse_show_error_message()
          write(0,*) 'ERROR in mdxyz file.'
          call modylas_abort()
       else
          call parse_get('/atom/natom',natom,lerr)
          if (.not. lerr) then
             call fmod_n(natom)
             if(natom .ne. 0) then
                call alloc_i2m(natom)
                call fmod_alloc_xyz(natom)
                call fmod_alloc_v(natom)
             endif
          endif
          getname='/atom/maxwell velocities/temperature'
          call parse_get(trim(getname),value1,lerr)
          if(.not.lerr) then
             call fmod_maxwell_temperature(value1)
             read_velocity = .false.
          else
             read_velocity = .true.
          endif
          call parse_get('/thermostat/nthermostat', nthermostat,lerr)
          do i = 1, nthermostat
             call parse_get('/thermostat/positions',value1,lerr)
             if (.not. lerr)  call fmod_set_rss(i-1, value1, 1)
          enddo
          do i =1, nthermostat
             call parse_get('/thermostat/velocities',value1,lerr)
             if (.not. lerr)  call fmod_set_vss(i-1, value1, 1)
          enddo
          call parse_get('/barostat/nbarostat', nbarostat,lerr)
          do i = 1, nbarostat
             call parse_get('/barostat/positions',value1,lerr)
             if (.not. lerr)  call fmod_set_rssb(i-1, value1, 1)
          enddo
          do i = 1, nbarostat
             call parse_get('/barostat/velocities',value1,lerr)
             if (.not. lerr)  call fmod_set_vssb(i-1, value1, 1)
          enddo
          call parse_get('/periodic/cell/length/x',qvalue1,lerr)
          if (.not. lerr) then
             qvalue1 = qvalue1 * 1.0d-10
             value1 = qvalue1
             value2 = qvalue1 * 0.5d0
             call fmod_cellx(value1)
             call fmod_cellxh(value2)
          endif
          call parse_get('/periodic/cell/length/y',qvalue1,lerr)
          if (.not. lerr) then
             qvalue1 = qvalue1 * 1.0d-10
             value1 = qvalue1
             value2 = qvalue1 * 0.5d0
             call fmod_celly(value1)
             call fmod_cellyh(value2)
          endif
          call parse_get('/periodic/cell/length/z',qvalue1,lerr)
          if (.not. lerr) then
             qvalue1 = qvalue1 * 1.0d-10
             value1 = qvalue1
             value2 = qvalue1 * 0.5d0
             call fmod_cellz(value1)
             call fmod_cellzh(value2)
          endif
          call parse_get('/periodic/cell/angle/alpha',value1,lerr)
          if (.not. lerr)  call fmod_set_alpha(value1)
          call parse_get('/periodic/cell/angle/beta',value1,lerr)
          if (.not. lerr)  call fmod_set_beta(value1)
          call parse_get('/periodic/cell/angle/gamma',value1,lerr)
          if (.not. lerr)  call fmod_set_gamma(value1)
          call parse_get('/periodic/cell/vboxg',value1,lerr1)
          call parse_get('/periodic/cell/vboxg',value2,lerr2)
          call parse_get('/periodic/cell/vboxg',value3,lerr3)
          call parse_get('/periodic/cell/vboxg',value4,lerr4)
          call parse_get('/periodic/cell/vboxg',value5,lerr5)
          call parse_get('/periodic/cell/vboxg',value6,lerr6)
          call parse_get('/periodic/cell/vboxg',value7,lerr7)
          call parse_get('/periodic/cell/vboxg',value8,lerr8)
          call parse_get('/periodic/cell/vboxg',value9,lerr9)
          if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. &
               &        lerr5 .or. lerr6 .or. lerr7 .or. lerr8 .or. lerr9)) then
             call fmod_set_vboxg(value1,value2,value3,value4,value5,value6,value7,value8,value9)
          endif
          call parse_close(lerr)
          if(myrank==0) write(*,*) 'Reading .mdxyz file ended successfully!'
       endif
       if (lerr) then
          call parse_abort(parse_show_error_message())
       endif
       !
       !       Direct reading of coordinates and velocities
       !
       call direct_read(trim(session_name)//'.mdxyz', read_velocity)
    endif


    !
    !     mdff.bin
    !
    if (.not.can_read_binary_file(trim(session_name)//'.mdff.bin')) then
!      if(myrank==0)write(6,*) "No ", trim(session_name), ".mdff.bin"
       goto 100 !GOTO READING MDFF FILE
    endif

    if(myrank==0) then
       open(f_mdff, file=trim(session_name)//'.mdff.bin', iostat=io, &
            &       status='old', access='sequential', form='unformatted')
    endif
    write(*,*) 'Reading .mdff.bin file started.'

    call read_forcefield
    call read_parameter_number
    call alloc
    call read_segment
    call read_molecule
    call read_mass
    call read_chgv
    call read_epsilon_sqrt
    call read_R_half
    call read_md_shake
    call read_md_void    ! void (Common)
    call read_md_LJ_special ! special (LJ)
    call read_md_coulomb_special  ! special (CL)
    call read_md_charmm_a_bond
    call read_md_bond_morse
    !call read_md_pcff_a_bond
    if(myrank==0) read(f_mdff) dummyInt
    call read_md_charmm_angleA
    call read_md_charmm_angleB
    !call read_md_pcff_angleA
    if(myrank==0) read(f_mdff) dummyInt
    !call read_md_pcff_angleB
    if(myrank==0) read(f_mdff) dummyInt
    call read_md_charmm_a_ub
    !call read_md_pcff_bbA
    if(myrank==0) read(f_mdff) dummyInt
    !call read_md_pcff_bbB
    if(myrank==0) read(f_mdff) dummyInt
    call read_md_charmm_dihedralA
    call read_md_charmm_dihedralB
    call read_md_charmm_CMAPA
    call read_md_charmm_CMAPB
    call read_md_charmm_CMAPC
    call read_md_charmm_CMAPD
    call read_md_charmm_CMAPE
    call read_md_charmm_itorsionA
    call read_md_charmm_itorsionB
    !ay+
    !      call read_position_constrain   !NOT support from binary input
    !ay+
    call read_CMAPVersion
    close(f_mdff)

    if(myrank==0) then
       write(*,*) 'Reading mdff.bin file ended successfully!'
       !       For debug force field.
       !       call write_memory_mdff
    endif
    goto 200 !GOTO READING MDDEF FILE

100 continue
    write(*,*) 'Reading .mdff file started.'

    !K. FUJIMOTO 2013.08.23
    call parse_open(trim(session_name)//'.mdff', lerr)

    call parse_get('/forcefield/type',chara,lerr)
    if(trim(chara)=='charmm')then
#ifdef OPLSAMBER
       if(myrank==0) then
          write(0,*) 'ERROR: -DOPLSAMBER is set ', &
               &     'in -DCMAKE_Fortran_FLAGS=" " for cmake. ',  &
               &     'It must be eliminated to handle input with type=charmm in .mdff file'
          call modylas_abort()
       endif
#endif
       md_condition__force_field=100 ! CHARMM
       call parse_get('/forcefield/CMAPVersion',ivalue1,lerr)
       if (lerr) then
          if (myrank == 0) then
             write(*,*) "###########################################"
             write(*,*) "CMAP Version is not set in your mdff file."
             write(*,*) "So, CMAP Version was set to be 36."
             write(*,*) "If you do not use CMAP Version 36,"
             write(*,*) "please write key word 'CMAPVersion=22'"
             write(*,*) "between <forcefiled> and </forcefiled> "
             write(*,*) "in your mdff file as follws."
             write(*,*)
             write(*,*) "<forcefiled>"
             write(*,*) "  type=charmm"
             write(*,*) "  CMAPVersion=22"
             write(*,*) "</forcefiled>"
             write(*,*) "###########################################"
          endif
       else
          call fmod_set_CMAPVersion(ivalue1)
       endif
    elseif(trim(chara)=='oplsaa')then
#ifndef OPLSAMBER
       if(myrank==0) then
          write(0,*)'ERROR: -DOPLSAMBER is necessary ', &
               &     'in -DCMAKE_Fortran_FLAGS=" " for cmake, ',  &
               &     'to handle input with type=oplsaa in .mdff file'
          call modylas_abort()
       endif
#endif
       md_condition__force_field=OPLSAA ! OPLS
       call parse_get('/forcefield/special_divide_lj',value1,lerr)
       if (.not. lerr)  call fmod_md_oplsaa_lj_sp_divide(value1)
       call parse_get('/forcefield/special_divide_coulomb', &
            &       value1,lerr)
       if (.not. lerr)  call fmod_md_oplsaa_cl_sp_divide(value1)
    elseif(trim(chara)=='amber')then
#ifndef OPLSAMBER
       if(myrank==0) then
          write(0,*)'ERROR: -DOPLSAMBER is necessary ', &
               &     'in -DCMAKE_Fortran_FLAGS=" " for cmake, ',  &
               &     'to handle input with type=amber in .mdff file'
          call modylas_abort()
       endif
#endif
       md_condition__force_field=AMBER ! AMBER
       call parse_get('/forcefield/special_divide_lj',value1,lerr)
       if (.not. lerr)  call fmod_md_oplsaa_lj_sp_divide(value1)
       call parse_get('/forcefield/special_divide_coulomb', &
            &       value1,lerr)
       if (.not. lerr)  call fmod_md_oplsaa_cl_sp_divide(value1)
!gaff
    elseif(trim(chara)=='gaff')then
#ifndef OPLSAMBER
       if(myrank==0) then
          write(0,*)'ERROR: -DOPLSAMBER is necessary ', &
               &     'in -DCMAKE_Fortran_FLAGS=" " for cmake, ',  &
               &     'to handle input with type=gaff in .mdff file'
#ifndef GAFF
          write(0,*)'ERROR: -DGAFF is necessary ', &
               &     'in -DCMAKE_Fortran_FLAGS=" " for cmake, ',  &
               &     'to handle input with type=gaff in .mdff file'
          call modylas_abort()
#endif
          call modylas_abort()
       endif
#endif
       md_condition__force_field=GAmbFF ! General Amber FF (GAFF)
       call parse_get('/forcefield/special_divide_lj',value1,lerr)
       if (.not. lerr)  call fmod_md_oplsaa_lj_sp_divide(value1)
       call parse_get('/forcefield/special_divide_coulomb', &
            &       value1,lerr)
       if (.not. lerr)  call fmod_md_oplsaa_cl_sp_divide(value1)
    elseif(trim(chara)=='kremer')then
#ifdef OPLSAMBER
       if(myrank==0) then
          write(0,*) 'ERROR: -DOPLSAMBER is set ', &
     &     'in -DCMAKE_Fortran_FLAGS=" " for cmake. ',  &
     &     'It must be eliminated to handle input with type=kremer in .mdff file'
          call modylas_abort()
       endif
#endif
       md_condition__force_field=KREMER ! Kremer's potential
    !elseif(trim(chara)=='pcff')  then
    !   md_condition__force_field=PCFF 
    endif  !! trim(chara)
    
    call parse_get('/molecules/nmolecule',nmol,lerr)
    call fmod_alloc_molinfo(nmol)
    if(.not. lerr) allocate(natommol(nmol))
    if(.not. lerr) natommol=0
    if(.not. lerr) allocate(nmolmol(nmol))
    if(.not. lerr) nmolmol=0
    do i = 1, nmol
       call parse_get('/atom/molecule',ivalue1, lerr)
       call parse_get('/atom/molecule',ivalue2, lerr)
       if(.not. lerr) nmolmol(i) = ivalue2
    enddo
    
    npara=0
    nall =0
    do i = 1, nmol
       call parse_specify_key('/molecules/molecule',i)
       call parse_get('/molecules/molecule/id',ivalue1,lerr)
       call parse_get('/molecules/molecule/natom',ivalue1,lerr)
       if(.not. lerr) then
          call fmod_set_molinfo(i-1,nmolmol(i),ivalue1)
          natommol(i) = ivalue1
       endif
       npara = npara + natommol(i)
       nall = nall + natommol(i)*nmolmol(i)
    enddo
    call fmod_set_param_num(npara,nmol,nmolmol,natommol)
    call alloc
    call fmod_alloc_atom2mol(nall)
    call fmod_alloc_mol2atom(nall)
    call fmod_alloc_irearrange(nall)
    
    hws = 0
    allocate(molnseg(nmol)) ; molnseg=0
    do i = 1, nmol
       call parse_specify_key('/molecules/molecule',i)
       getname = '/molecules/molecule/segments/nsegment'
       call parse_get(trim(getname),molnseg(i),lerr)
       if (.not. lerr) then
          hws = hws + molnseg(i) * nmolmol(i)
       endif
    enddo
    if(hws .ne. 0) then
       nsegments = hws
       call fmod_alloc_seg_natoms(nsegments)
       call fmod_alloc_seg_center_x(nsegments)
       call fmod_alloc_seg_center_y(nsegments)
       call fmod_alloc_seg_center_z(nsegments)
       call fmod_alloc_segtop(nsegments)
    endif
    
    do i = 1, nmol
       call fmod_set_mol_ID(i-1)
       call parse_specify_key('/molecules/molecule',i)
       call fmod_set_mol_seg_num(molnseg(i))
       segname='/molecules/molecule/segments/segment'
       do j = 1, molnseg(i)
          call parse_specify_key(segname,j)
          getname='/molecules/molecule/segments/segment/ID'
          call parse_get(trim(getname),sid,lerr1)
          getname='/molecules/molecule/segments/segment/natom'
          call parse_get(trim(getname),nas,lerr2)
          if (.not. (lerr1 .or. lerr2)) then
             call fmod_set_mol_seg_natom(sid, nas)
          endif
          do k = 1, nas
             getname='/molecules/molecule/segments/segment/atom'
             call parse_get(trim(getname),ivalue1,lerr)
             if (.not. (lerr1 .or. lerr)) then
                call fmod_set_mol_seg2atom(sid,k-1,ivalue1)
             endif
          enddo
       enddo
       call fmod_set_mol_mol_natoms()
       call fmod_set_mol_mol2atom()
    enddo
    
    
    do i = 1, nmol
       call parse_specify_key('/molecules/molecule',i)
       !READING MASS AND R_MASS
       do j = 1, natommol(i)
          call parse_get('/molecules/molecule/mass',value1, lerr)
          if (.not. lerr)  call fmod_set_mol_mass(value1)
       enddo
       !READING CHARGE
       do j = 1, natommol(i)
          call parse_get('/molecules/molecule/charge',value1,lerr)
          if (.not. lerr)  call fmod_set_mol_chgv(value1)
       enddo
       !READING LJ EPSILON
       do j = 1, natommol(i)
          call parse_get('/molecules/molecule/epsilon',qvalue1,lerr)
          if (.not. lerr) then
             !convert unit:[kcal/mol] => [J]
             value1 = sqrt(qvalue1*1000.0d0 * md_CALORIE / md_AVOGADRO)
             call fmod_set_mol_ljs__epsilon_sqrt(value1)
          endif
       enddo
       !READING LJ R-HALF
       do j = 1, natommol(i)
          call parse_get('/molecules/molecule/r',qvalue1,lerr)
          if (.not. lerr) then
             if( md_condition__force_field == CHARMM .or. &
                  &              md_condition__force_field == KREMER  )then
                !convert unit[A] => [m]
                value1 = qvalue1 * 1.0d-10 * 0.5d0
                call fmod_set_mol_ljs__r_half(value1)
             elseif(md_condition__force_field == OPLSAA .or. &
                  & md_condition__force_field == AMBER  .or. &
                  & md_condition__force_field == GAmbFF  )then
                !convert unit[A] => [m]
                qvalue1 = qvalue1 * 1.0d-10 * 0.5d0
                value1 = qvalue1 * 2.0d0 ** (1.0d0/6.0d0)
                call fmod_set_mol_ljs__r_half(value1)
             endif
          endif
       enddo
    enddo
    
    call fmod_alloc_yyparse_shake(npara)
    getname = '/molecules/molecule/shake pair'
    do i = 1, nmol
       call parse_specify_key('/molecules/molecule',i)
       temp = 0
       do while(.true.)
          call parse_get(trim(getname), ivalue1,lerr)
          if(lerr) exit
          call parse_get(trim(getname), ivalue2,lerr)
          call parse_get(trim(getname), qvalue1,lerr)
          value1 = qvalue1 * 1.0d-10
          call add_to_mol_shake_list(i,ivalue1,ivalue2,value1)
          temp = temp + 1
       enddo
       totnconst = totnconst + temp * nmolmol(i)
    enddo
    call set_shake()
    
    ! READING VOID PAIR
    call fmod_alloc_yyparse_lj_top(npara)
    call fmod_alloc_yyparse_lj_n(npara)
    call fmod_alloc_yyparse_coulomb_top(npara)
    call fmod_alloc_yyparse_coulomb_n(npara)
    do i =1, nmol
       !! LJ void !!
       call parse_specify_key('/molecules/molecule',i)
       getname = '/molecules/molecule/nvoidpair_lj'
       call parse_get(trim(getname), nvl, lerr)
       getname = '/molecules/molecule/lj void pair'
       do j = 1, nvl
          call parse_get(trim(getname), ivalue1,lerr1)
          call parse_get(trim(getname), ivalue2,lerr2)
          if (.not. (lerr1 .or. lerr2)) then
             call add_to_mol_lj_void(i, ivalue1, ivalue2)
          endif
       enddo
       !! COULOMB VOID PAIR
       getname = '/molecules/molecule/nvoidpair_coulomb'
       call parse_get(trim(getname),nvc, lerr)
       getname = '/molecules/molecule/coulomb void pair'
       do j = 1, nvc
          call parse_get(trim(getname),ivalue1, lerr1)
          call parse_get(trim(getname),ivalue2, lerr2)
          if (.not. (lerr1 .or. lerr2)) then
             call add_to_mol_coulomb_void(i,ivalue1,ivalue2)
          endif
       enddo
    enddo
    call set_md_void()
    
    !READING SPECAIL PAIR
    do i = 1, nmol
       !! LJ special PAIR !!
       call parse_specify_key('/molecules/molecule',i)
       getname = '/molecules/molecule/nspecialpair_lj'
       call parse_get(trim(getname), nsl, lerr)
       getname = '/molecules/molecule/lj special pair'
       do j = 1, nsl
          call parse_get(trim(getname),ivalue1,lerr1)
          call parse_get(trim(getname),ivalue2,lerr2)
          call parse_get(trim(getname),qvalue1,lerr3)
          call parse_get(trim(getname),qvalue2,lerr4)
          if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4)) then
             if( md_condition__force_field == CHARMM .or. md_condition__force_field == KREMER )then
                !converting unit: [kcal/mol] => [J]
                !converting unit: [A] => [m]
                !epsilon needs not to be square-rooted and R not halfed
                value1 = qvalue1 * 1000d0 * md_CALORIE / md_AVOGADRO
                value2 = qvalue2 * 1.0d-10
                call add_to_mol_lj_list(i,ivalue1,ivalue2,value1,value2)
             elseif( md_condition__force_field == OPLSAA .or. md_condition__force_field == AMBER )then
                continue !! nothing to do in OPLS/AMBER
             endif
          endif
       enddo
       !! COUBLOMB special PAIR !!
       getname = '/molecules/molecule/nspecialpair_coulomb'
       call parse_get(trim(getname),nsct,lerr)
       getname = '/molecules/molecule/coulomb special pair'
       do j = 1, nsct
          call parse_get(trim(getname),ivalue1,lerr1)
          call parse_get(trim(getname),ivalue2,lerr2)
          call parse_get(trim(getname),value1,lerr3)
          if (.not. (lerr1 .or. lerr2 .or. lerr3)) then
             call add_to_mol_coulomb_list(i,ivalue1,ivalue2,value1)
          endif
       enddo
    enddo
    call set_md_lj_special()
    call set_md_coulomb_special()
    
    
    !READING BOND PAIR
    call fmod_alloc_yyparse_bond_top(npara)
    call fmod_alloc_yyparse_bond_n(npara)
    do i = 1, nmol
       call parse_specify_key('/molecules/molecule',i)
       call parse_get('/molecules/molecule/nbond',nbond,lerr)
       do j = 1, nbond
          call parse_get('/molecules/molecule/bond',ivalue1,lerr1)
          call parse_get('/molecules/molecule/bond',ivalue2,lerr2)
          call parse_get('/molecules/molecule/bond',qvalue1,lerr3)
          call parse_get('/molecules/molecule/bond',qvalue2,lerr4)
          if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4)) then
             !convert unit : [kcal /mol A^2] => [J / m^2]
             value1 = qvalue1*1000.0d0*md_CALORIE/md_AVOGADRO*1.0d20
             !convert unit : [A] => [m]
             value2 = qvalue2 * 1.0d-10
             call  add_to_mol_bond_list(i,ivalue1,ivalue2,value1,value2)
          endif
       enddo
    enddo
    call set_md_charmm_a_bond()

    !READING BOND MORSE PAIR
    call fmod_alloc_yyparse_bond_morse_top(npara)
    call fmod_alloc_yyparse_bond_morse_n(npara)
    do i = 1, nmol
       call parse_specify_key('/molecules/molecule',i)
       call parse_get('/molecules/molecule/nbondmorse',nbond,lerr)
       do j = 1, nbond
          call parse_get('/molecules/molecule/bond morse',ivalue1,lerr1)
          call parse_get('/molecules/molecule/bond morse',ivalue2,lerr2)
          call parse_get('/molecules/molecule/bond morse',qvalue1,lerr3)
          call parse_get('/molecules/molecule/bond morse',qvalue2,lerr4)
          call parse_get('/molecules/molecule/bond morse',qvalue3,lerr5)
          if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. lerr5)) then
             !convert unit : [kcal /mol] => [J / m^2]
             value1 = qvalue1*1000.0d0*md_CALORIE/md_AVOGADRO
             !convert unit : [A-1] => [m-1]
             value2 = qvalue2 * 1.0d10
             !convert unit : [A] => [m]
             value3 = qvalue3 * 1.0d-10
             call  add_to_mol_bond_morse_list(i,ivalue1,ivalue2,value1,value2, value3)
          endif
       enddo
    enddo
    call set_md_bond_morse()
    !READING BOND PCFF PAIR
    !call fmod_alloc_yyparse_bond_pcff_top(npara)
    !call fmod_alloc_yyparse_bond_pcff_n(npara)
    !do i = 1, nmol
    !   call parse_specify_key('/molecules/molecule',i)
    !   call parse_get('/molecules/molecule/nbond_pcff',nbond,lerr)
    !   do j = 1, nbond
    !      call parse_get('/molecules/molecule/bond_pcff',ivalue1,lerr1)
    !      call parse_get('/molecules/molecule/bond_pcff',ivalue2,lerr2)
    !      call parse_get('/molecules/molecule/bond_pcff',qvalue1,lerr3)
    !      call parse_get('/molecules/molecule/bond_pcff',qvalue2,lerr4)
    !      call parse_get('/molecules/molecule/bond_pcff',qvalue3,lerr5)
    !      call parse_get('/molecules/molecule/bond_pcff',qvalue4,lerr6)
    !      if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 &
    !            .or. lerr5 .or. lerr6)) then
    !         !convert unit : [kcal /mol A^2] => [J / m^2]
    !         value1 = qvalue1*1000.0d0*md_CALORIE/md_AVOGADRO*1.0d20
    !         !convert unit : [kcal /mol A^3] => [J / m^3]
    !         value2 = qvalue2*1000.0d0*md_CALORIE/md_AVOGADRO*1.0d30
    !         !convert unit : [kcal /mol A^4] => [J / m^4]
    !         value3 = qvalue3*1000.0d0*md_CALORIE/md_AVOGADRO*1.0d40
    !         !convert unit : [A] => [m]
    !         value4 = qvalue4 * 1.0d-10
    !         call add_to_mol_bond_pcff_list &
    !              (i,ivalue1,ivalue2,value1,value2,value3,value4)
    !      endif
    !   enddo
    !enddo
    !call set_md_pcff_a_bond()
    
    
    !READING ANGLE PAIR
    call fmod_alloc_yyparse_angleA_top(npara)
    call fmod_alloc_yyparse_angleA_n(npara)
    call fmod_alloc_yyparse_angleB_top(npara)
    call fmod_alloc_yyparse_angleB_n(npara)
    do i = 1, nmol
       call parse_specify_key('/molecules/molecule',i)
       call parse_get('/molecules/molecule/nangle',nag,lerr)
       do j = 1, nag
          call parse_get('/molecules/molecule/angle',ivalue1,lerr1)
          call parse_get('/molecules/molecule/angle',ivalue2,lerr2)
          call parse_get('/molecules/molecule/angle',ivalue3,lerr3)
          call parse_get('/molecules/molecule/angle',qvalue1,lerr4)
          call parse_get('/molecules/molecule/angle',qvalue2,lerr5)
          if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. lerr5)) then
             !convert unit : [kcal/mol.rad**2] => [J/rad**2]
             value1 = qvalue1 * 1000.0d0 * md_CALORIE/md_AVOGADRO
             !convert unit : [degree]=>[rad]
             value2 = qvalue2 * md_DEGREE
             call  add_to_mol_angle_list(i,ivalue1,ivalue2,ivalue3,value1,value2)
          endif
       enddo
    enddo
    call set_md_charmm_angleA
    call set_md_charmm_angleB
    
    !READING ANGLE PCFF PAIR
    !call fmod_alloc_yyparse_angle_pcffA_top(npara)
    !call fmod_alloc_yyparse_angle_pcffA_n(npara)
    !call fmod_alloc_yyparse_angle_pcffB_top(npara)
    !call fmod_alloc_yyparse_angle_pcffB_n(npara)
    !do i = 1, nmol
    !   call parse_specify_key('/molecules/molecule',i)
    !   call parse_get('/molecules/molecule/nangle_pcff',nag,lerr)
    !   do j = 1, nag
    !      call parse_get('/molecules/molecule/angle_pcff',ivalue1,lerr1)
    !      call parse_get('/molecules/molecule/angle_pcff',ivalue2,lerr2)
    !      call parse_get('/molecules/molecule/angle_pcff',ivalue3,lerr3)
    !      call parse_get('/molecules/molecule/angle_pcff',qvalue1,lerr4)
    !      call parse_get('/molecules/molecule/angle_pcff',qvalue2,lerr5)
    !      call parse_get('/molecules/molecule/angle_pcff',qvalue3,lerr6)
    !      call parse_get('/molecules/molecule/angle_pcff',qvalue4,lerr7)
    !      if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. lerr5 &
    !            .or. lerr6 .or. lerr7 )) then
    !         !convert unit : [kcal/mol.rad**2] => [J/rad**2]
    !         value1 = qvalue1 * 1000.0d0 * md_CALORIE/md_AVOGADRO
    !         !convert unit : [kcal/mol.rad**3] => [J/rad**3]
    !         value2 = qvalue2 * 1000.0d0 * md_CALORIE/md_AVOGADRO
    !         !convert unit : [kcal/mol.rad**4] => [J/rad**4]
    !         value3 = qvalue3 * 1000.0d0 * md_CALORIE/md_AVOGADRO
    !         !convert unit : [degree]=>[rad]
    !         value4 = qvalue4 * md_DEGREE
    !         call  add_to_mol_angle_pcff_list&
    !               (i,ivalue1,ivalue2,ivalue3,value1,value2,value3,value4)
    !      endif
    !   enddo
    !enddo
    !call set_md_pcff_angleA
    !call set_md_pcff_angleB
    
    !READING UB PAIR
    call fmod_alloc_yyparse_ub_top(npara)
    call fmod_alloc_yyparse_ub_n(npara)
    do i = 1, nmol
       call parse_specify_key('/molecules/molecule',i)
       call parse_get('/molecules/molecule/nub',nub,lerr)
       do  j = 1, nub
          call parse_get('/molecules/molecule/ub',ivalue1,lerr1)
          call parse_get('/molecules/molecule/ub',ivalue2,lerr2)
          call parse_get('/molecules/molecule/ub',ivalue3,lerr3)
          call parse_get('/molecules/molecule/ub',qvalue1,lerr4)
          call parse_get('/molecules/molecule/ub',qvalue2,lerr5)
          if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. lerr5)) then
             !convert unit : pkcal/mol.A**2] => [J/m**2]
             value1 = qvalue1*1000.0d0*md_CALORIE/md_AVOGADRO*1.0d20
             !convert unit : [A] => [m]
             value2 = qvalue2 * 1.0d-10
             call  add_to_mol_ub_list(i,ivalue1,ivalue2,ivalue3,value1,value2)
          endif
       enddo
    enddo
    call set_md_charmm_a_ub()
    
    !READING DIHEDRAL PAIRS
    call fmod_alloc_yyparse_dihedA_top(npara)
    call fmod_alloc_yyparse_dihedA_n(npara)
    call fmod_alloc_yyparse_dihedB_top(npara)
    call fmod_alloc_yyparse_dihedB_n(npara)
#ifdef DIHEDRAL_TABLE
       call init_dihedral_table()
#endif
    do i = 1, nmol
       call parse_specify_key('/molecules/molecule',i)
       call parse_get('/molecules/molecule/ndihedral',ndi,lerr)
       do j = 1, ndi
          call parse_get('/molecules/molecule/dihedral', ivalue1,lerr1)
          call parse_get('/molecules/molecule/dihedral', ivalue2,lerr2)
          call parse_get('/molecules/molecule/dihedral', ivalue3,lerr3)
          call parse_get('/molecules/molecule/dihedral', ivalue4,lerr4)
          call parse_get('/molecules/molecule/dihedral', qvalue1,lerr5)
          call parse_get('/molecules/molecule/dihedral', qvalue2,lerr6)
          call parse_get('/molecules/molecule/dihedral', qvalue3,lerr7)
          if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. &
               &                   lerr4 .or. lerr5 .or. lerr6 .or. lerr7)) then
             
             if(    md_condition__force_field == CHARMM .or. &
                  & md_condition__force_field == KREMER .or. &
                  & md_condition__force_field == GAmbFF  ) then

                !convert unit : [kcal/mol] => [J]
                value1 = qvalue1 * 1000.0d0 * md_CALORIE/md_AVOGADRO
                value2 = qvalue2
                !convert unit : [degree] => [rad]
                value3 = qvalue3 * md_DEGREE
                call add_to_mol_dihedral_list(i, ivalue1, ivalue2, ivalue3, ivalue4, value1, value2, value3)
                
             elseif(md_condition__force_field == OPLSAA .or. md_condition__force_field == AMBER  )then
                
                !convert unit : [kcal/mol] => [J]
                qvalue1 = qvalue1 * 1000.0d0 * md_CALORIE/md_AVOGADRO
                value1 = qvalue1 * 0.5d0
                value2 = qvalue2
                !convert unit : [degree] => [rad]
                value3 = qvalue3 * md_DEGREE
                value3 =         - value3
                call add_to_mol_dihedral_list(i, ivalue1, ivalue2, ivalue3, ivalue4, value1, value2, value3)
             else
                write(0,*) 'ERROR: dihedral reading'
                call modylas_abort()
             endif
          endif
       enddo
    enddo
    call set_md_charmm_dihedralA()
    call set_md_charmm_dihedralB()
#ifdef DIHEDRAL_TABLE
       call finalize_dihedral_table()
#endif
    !	stop 'DIHED'
    
    if(    md_condition__force_field==CHARMM) then
       !READING CMAP PAIRS
       call fmod_alloc_yyparse_CMAP_topA(npara)
       call fmod_alloc_yyparse_CMAP_nA(npara)
       call fmod_alloc_yyparse_CMAP_topB(npara)
       call fmod_alloc_yyparse_CMAP_nB(npara)
       call fmod_alloc_yyparse_CMAP_topC(npara)
       call fmod_alloc_yyparse_CMAP_nC(npara)
       call fmod_alloc_yyparse_CMAP_topD(npara)
       call fmod_alloc_yyparse_CMAP_nD(npara)
       call fmod_alloc_yyparse_CMAP_topE(npara)
       call fmod_alloc_yyparse_CMAP_nE(npara)
       do i =1, nmol
          call parse_specify_key('/molecules/molecule',i)
          call parse_get('/molecules/molecule/ncmap',ncm,lerr)
          do j = 1, ncm
             call parse_get('/molecules/molecule/CMAP',ivalue1,lerr1)
             call parse_get('/molecules/molecule/CMAP',ivalue2,lerr2)
             call parse_get('/molecules/molecule/CMAP',ivalue3,lerr3)
             call parse_get('/molecules/molecule/CMAP',ivalue4,lerr4)
             call parse_get('/molecules/molecule/CMAP',ivalue5,lerr5)
             call parse_get('/molecules/molecule/CMAP',ivalue6,lerr6)
             if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. lerr5 .or. lerr6 )) then
                call add_to_mol_cmap_list(i, ivalue1, ivalue2, ivalue3, ivalue4, ivalue5, ivalue6)
             endif
          enddo
       enddo
       call set_md_charmm_CMAPA()
       call set_md_charmm_CMAPB()
       call set_md_charmm_CMAPC()
       call set_md_charmm_CMAPD()
       call set_md_charmm_CMAPE()
    endif
    
    !READING ITORSION PAIRS
    call fmod_alloc_yyparse_itorsA_top(npara)
    call fmod_alloc_yyparse_itorsA_n(npara)
    call fmod_alloc_yyparse_itorsB_top(npara)
    call fmod_alloc_yyparse_itorsB_n(npara)
    do i =1, nmol
       call parse_specify_key('/molecules/molecule',i)
       call parse_get('/molecules/molecule/nitorsion',nit,lerr)
       do j = 1, nit
#ifdef GROEXT_OPLS
          call parse_get('/molecules/molecule/itorsion', ivalue1,lerr1)
          call parse_get('/molecules/molecule/itorsion', ivalue2,lerr2)
          call parse_get('/molecules/molecule/itorsion', ivalue3,lerr3)
          call parse_get('/molecules/molecule/itorsion', ivalue4,lerr4)
          call parse_get('/molecules/molecule/itorsion', qvalue1,lerr5)
          call parse_get('/molecules/molecule/itorsion', qvalue2,lerr6)
          call parse_get('/molecules/molecule/itorsion', qvalue3,lerr7)
          !convert unit : [kcal/mol.rad**2] => [J/rad**2]
          qvalue1 = qvalue1 * 1000.0d0 * md_CALORIE / md_AVOGADRO
          value1 = qvalue1 * 0.5d0
          value2 = qvalue2
          !convert unit:[degree]=>[rad]
          value3 = qvalue3 * md_DEGREE
          if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. lerr5 .or. lerr6 .or. lerr7)) then
             call add_to_mol_itorsion_list(i, ivalue1, ivalue2, ivalue3, ivalue4, value1, value2, value3)
          endif
#elif defined(GAFF)
          call parse_get('/molecules/molecule/itorsion', ivalue1,lerr1)
          call parse_get('/molecules/molecule/itorsion', ivalue2,lerr2)
          call parse_get('/molecules/molecule/itorsion', ivalue3,lerr3)
          call parse_get('/molecules/molecule/itorsion', ivalue4,lerr4)
          call parse_get('/molecules/molecule/itorsion', qvalue1,lerr5)
          call parse_get('/molecules/molecule/itorsion', ivalue5,lerr6)
          call parse_get('/molecules/molecule/itorsion', ivalue6,lerr7)
          !E=K[1+d*cos(n*phi)]
          !convert unit : [kcal/mol] => [J]
          value1 = qvalue1 * 1000.0d0 * md_CALORIE / md_AVOGADRO ! K
          value2 = ivalue5   ! d
          value3 = ivalue6   ! n
          if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. lerr5 .or. lerr6 .or. lerr7)) then
             call add_to_mol_itorsion_list(i, ivalue1, ivalue2, ivalue3, ivalue4, value1, value2, value3)
          endif
#else
          call parse_get('/molecules/molecule/itorsion', ivalue1,lerr1)
          call parse_get('/molecules/molecule/itorsion', ivalue2,lerr2)
          call parse_get('/molecules/molecule/itorsion', ivalue3,lerr3)
          call parse_get('/molecules/molecule/itorsion', ivalue4,lerr4)
          call parse_get('/molecules/molecule/itorsion', qvalue1,lerr5)
          call parse_get('/molecules/molecule/itorsion', qvalue2,lerr6)
          !convert unit : [kcal/mol.rad**2] => [J/rad**2]
          value1 = qvalue1 * 1000.0d0 * md_CALORIE / md_AVOGADRO
          !convert unit:[degree]=>[rad]
          value2 = qvalue2 * md_DEGREE
          if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. lerr5 .or. lerr6)) then
             call add_to_mol_itorsion_list(i, ivalue1, ivalue2, ivalue3, ivalue4, value1, value2)
          endif
#endif
       enddo
    enddo
    call set_md_charmm_itorsionA()
    call set_md_charmm_itorsionB()
    
    !--- read position constrain --
    call parse_get('/position constrain/type', chara, lerr)
    if (.not. lerr) then
       if(trim(chara) =='harmonic') then
          call fmod_p_cons__type(1)
          getname = '/position constrain/force_constant'
          call parse_get(trim(getname), value1, lerr)
          if (.not. lerr)  then
             !convert unit : [kcal /mol A^2] => [J / m^2]
             value1 = value1*1000.0d0*md_CALORIE/md_AVOGADRO*1.0d20
             call fmod_p_cons__fc(value1)
          endif
       else if(trim(chara) =='huge_mass') then
          call fmod_p_cons__type(2)
       else if(trim(chara) =='fix') then
          call fmod_p_cons__type(3)
       endif
    endif
    getname = '/position constrain/atom/atom_type'
    call parse_get(trim(getname), chara, lerr)
    if (.not. lerr) then
       if(trim(chara) =='all') then
          call fmod_p_cons_atom__atom_type(0)
       else if(trim(chara) =='heavy') then
          call fmod_p_cons_atom__atom_type(1)
       endif
    endif
    getname = '/position constrain/atom/molecules'
    call parse_get(trim(getname), chara, lerr)
    if (.not. lerr) then
       if(trim(chara) =='all') then
          call fmod_p_cons_atom__molecules(0)
       else if(trim(chara) =='specify') then
          call fmod_p_cons_atom__molecules(1)
          getname = '/position constrain/atom/nmolecule'
          call parse_get(trim(getname), ivalue1, lerr)
          if (.not. lerr)  then
             call fmod_p_cons_atom__nmolecule(ivalue1)
             call fmod_alloc_mol_p_constrain(ivalue1)
          endif
          k = ivalue1
          do  i= 1,k
             getname = '/position constrain/atom/molecule'
             call parse_get(trim(getname),ivalue1,lerr1)
             call parse_get(trim(getname),ivalue2,lerr2)
             if (.not. (lerr1 .or. lerr2 )) then
                im = 0
                do j= 1,ivalue1
                   im = im + nmolmol(j)
                enddo
                im = im + ivalue2 + 1
                call add_to_mol_p_constrain(i,im)
             endif
          enddo
       endif
    endif
    getname = '/position constrain/atom/natom_allow_unconstrain'
    call parse_get(trim(getname), ivalue1, lerr)
    if (.not. lerr)  then
       call fmod_p_cons_atom__natm_alw_ucn(ivalue1)
       call fmod_alloc_natom_allow_uncons(ivalue1)
       k = ivalue1
       do  i= 1,k
          getname = '/position constrain/atom/allow atom unconstrain'
          call parse_get(trim(getname),ivalue1, lerr)
          if(.not.lerr)call add_to_atom_allow_uncons(i,ivalue1+1)
       enddo
    endif
    getname = '/position constrain/position/type'
    call parse_get(trim(getname), chara, lerr)
    if (.not. lerr) then
       if(trim(chara) =='initial') then
          call fmod_p_cons_position__type(0)
       else if(trim(chara) =='specify') then
          call fmod_p_cons_position__type(1)
          getname = '/position constrain/position/coordinate'
          call fmod_alloc_ref_crd0_p_cons
          do while(.true.)
             call parse_get(trim(getname),ivalue1,lerr)
             if (lerr) exit
             call parse_get(trim(getname),qvalue1,lerr1)
             call parse_get(trim(getname),qvalue2,lerr2)
             call parse_get(trim(getname),qvalue3,lerr3)
             if ( lerr1 .or. lerr2 .or. lerr3 ) then
                write(0,*)  'ERROR:/position constrain/position/coordinate'
                call modylas_abort()
             endif
             value1 = qvalue1 * 1.0d-10
             value2 = qvalue2 * 1.0d-10
             value3 = qvalue3 * 1.0d-10
             call fmod_p_cons_position__ref_crd0(ivalue1+1, value1,value2,value3 )
          enddo
       endif
    endif

    call parse_close(lerr)

    if(myrank==0) write(*,*) 'Reading .mdff file ended successfully!'

    !
    !     mddef  0.9.0b
    !
200 continue
    if (read_mddef) then
       call parse_open(trim(session_name)//'.mddef', lerr)
       if(myrank==0) write(*,*) 'Reading .mddef file started.'

       if (.not. lerr) then
          call parse_get('/output/ascii', chara, lerr)
          if(trim(chara) == 'yes') then
             call fmod_ascii_output(1)
          else if(trim(chara) == 'no') then
             call fmod_ascii_output(0)
          endif
          !backup files starting with #
          call parse_get('/output/backup', chara, lerr)
          IF(trim(chara) == 'yes') THEN
             call fmod_backup_output(1)
          else if(trim(chara) == 'no') then
             call fmod_backup_output(0)
          ENDIF
          call parse_get('/output/restart/start', ivalue1, lerr)
          if (.not. lerr)  call fmod_restart_start(ivalue1)
          call parse_get('/output/restart/interval', ivalue1, lerr)
          if (.not. lerr)  call fmod_restart_interval(ivalue1)
          call parse_get('/output/monitor/start', ivalue1, lerr)
          if (.not. lerr)  call fmod_mntr_start(ivalue1)
          call parse_get('/output/monitor/interval', ivalue1, lerr)
          if (.not. lerr)  call fmod_mntr_interval(ivalue1)
          call parse_get('/output/force/start', ivalue1, lerr)
          if (.not. lerr)  call fmod_force_start(ivalue1)
          call parse_get('/output/force/interval', ivalue1, lerr)
          if (.not. lerr)  call fmod_force_interval(ivalue1)
          !trj
          call parse_get('/output/mdtrj', chara, lerr)
          IF(trim(chara) == 'yes') THEN
             call fmod_mdtrj_output(1)
          ENDIF
          call parse_get('/output/trajectory/start', ivalue1, lerr)
          if (.not. lerr)  call fmod_trj_start(ivalue1)
          call parse_get('/output/trajectory/interval', ivalue1, lerr)
          if (.not. lerr)  call fmod_trj_interval(ivalue1)
          !dcd
          call parse_get('/output/dcd', chara, lerr)
          IF(trim(chara) == 'yes') THEN
             call fmod_dcd_output(1)
          ENDIF
          call parse_get('/output/trjdcd/start', ivalue1, lerr)
          if (.not. lerr)  call fmod_dcd_start(ivalue1)
          call parse_get('/output/trjdcd/interval', ivalue1, lerr)
          if (.not. lerr)  call fmod_dcd_interval(ivalue1)
#ifdef XTC
          !xtc
          call parse_get('/output/xtc', chara, lerr)
          IF(trim(chara) == 'yes') THEN
             call fmod_xtc_output(1)
          ENDIF
          call parse_get('/output/trjxtc/start', ivalue1, lerr)
          if (.not. lerr)  call fmod_xtc_start(ivalue1)
          call parse_get('/output/trjxtc/interval', ivalue1, lerr)
          if (.not. lerr)  call fmod_xtc_interval(ivalue1)
#endif
          !dcd
!#ifdef XTC
!          !xtc
!          call parse_get('/output/xtc', chara, lerr)
!          IF(trim(chara) == 'yes') THEN
!             call fmod_xtc_output(1)
!          ENDIF
!          call parse_get('/output/trjxtc/start', ivalue1, lerr)
!          if (.not. lerr)  call fmod_xtc_start(ivalue1)
!          call parse_get('/output/trjxtc/interval', ivalue1, lerr)
!          if (.not. lerr)  call fmod_xtc_interval(ivalue1)
!#endif
          call parse_get('/output/system_dipole',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
            call set_up_system_dipole(1)
          endif
          call parse_get('/randomseed/seed', ivalue1, lerr)
          if (.not. lerr)  call fmod_randomseed(ivalue1)
! TZYBEGIN
          call parse_get('/condition/dt', value1, lerr)
         !  if (.not. lerr)  call fmod_md_condition__dt(value1)
          if (.not. lerr) then
              call fmod_md_condition__dt(value1)
              if( myrank == 0 ) write(*,*) "dt is deprecated. Use dt_long instead"
          endif
          call parse_get('/condition/dt_long', value1, lerr)
          if (.not. lerr)  call fmod_md_condition__dt(value1)
! TZYEND
          call parse_get('/condition/initial_step', ivalue1, lerr)
          if(.not. lerr) then
             call fmod_mdstep(ivalue1) !generic step
          else
             call fmod_mdstep(0) !generic step
          endif
          !
          call parse_get('/condition/velocity_scaling',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             velocity_scaling_flag = .true.
          endif
          !
          call parse_get('/condition/volume_scaling',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             volume_scaling_flag = .true.
          endif
          call parse_get('/condition/volume', value1, lerr)
          if (.not. lerr)  call fmod_sysvolume(value1)
          ! 
          call parse_get('/condition/cellgeometry_scaling',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             cellgeometry_scaling_flag = .true.
          endif
          call parse_get('/condition/cellx', value1, lerr)
          if (.not. lerr)  call fmod_syscellx(value1)
          call parse_get('/condition/celly', value1, lerr)
          if (.not. lerr)  call fmod_syscelly(value1)
          call parse_get('/condition/cellz', value1, lerr)
          if (.not. lerr)  call fmod_syscellz(value1)
          call parse_get('/condition/alpha', value1, lerr)
          if (.not. lerr)  call fmod_sysalpha(value1)
          call parse_get('/condition/beta', value1, lerr)
          if (.not. lerr)  call fmod_sysbeta(value1)
          call parse_get('/condition/gamma', value1, lerr)
          if (.not. lerr)  call fmod_sysgamma(value1)
          !
          if (volume_scaling_flag)    call fmod_volume_scaling(1)
          if (cellgeometry_scaling_flag)    call fmod_cellgeometry_scaling(1)
          !
          call parse_get('/condition/ensemble', chara, lerr)
          if (trim(chara) == 'nve') then
             call fmod_md_condition__ensemble(NVE)
             if (velocity_scaling_flag)  call fmod_velocity_scaling(.true., .false.)
          elseif (trim(chara) == 'nvt') then
             call fmod_md_condition__ensemble(NVT)
             getname='/thermostat/tau_Q'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauqtherm(value1)
             !ya
             call parse_get('/thermostat/initialize',chara,lerr)
             if (.not. lerr .and. trim(chara) == 'yes') then
                initialize_thermostat_flag = .true.
             endif
             if (velocity_scaling_flag) then
                call parse_abort('Velocity scaling should be used with NVE ensemble!')
             endif
          elseif (trim(chara) == 'nvlxlylzt') then
             call fmod_md_condition__ensemble(NVLXLYLZT)
             getname='/thermostat/tau_Q'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauqtherm(value1)
             !ya
             call parse_get('/thermostat/initialize',chara,lerr)
             if (.not. lerr .and. trim(chara) == 'yes') then
                initialize_thermostat_flag = .true.
             endif

          elseif (trim(chara) == 'npt_a') then

             call fmod_md_condition__ensemble(NPT_A)
             if (velocity_scaling_flag) then
                call parse_abort('Velocity scaling cannot be used in NPT.')
             endif
             getname='/thermostat/tau_Q'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauqtherm(value1)
             call parse_get('/thermostat/initialize',chara,lerr)
             if (.not. lerr .and. trim(chara) == 'yes') then
                initialize_thermostat_flag = .true.
             endif
             getname='/barostat/tau_Q'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauqbaro(value1)
             getname='/barostat/tau_W'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauwpres(value1)
             call parse_get('/barostat/initialize',chara,lerr)
             if (.not. lerr .and. trim(chara) == 'yes') then
                initialize_barostat_flag = .true.
             endif

          elseif (trim(chara) == 'npt_z') then

             call fmod_md_condition__ensemble(NPT_Z)
             if (velocity_scaling_flag) then
                call parse_abort('Velocity scaling cannot be used in NPT.')
             endif
             getname='/thermostat/tau_Q'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauqtherm(value1)
             call parse_get('/thermostat/initialize',chara,lerr)
             if (.not. lerr .and. trim(chara) == 'yes') then
                initialize_thermostat_flag = .true.
             endif
             getname='/barostat/tau_Q'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauqbaro(value1)
             getname='/barostat/tau_W'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauwpres(value1)
             call parse_get('/barostat/initialize',chara,lerr)
             if (.not. lerr .and. trim(chara) == 'yes') then
                initialize_barostat_flag = .true.
             endif

          elseif (trim(chara) == 'nptlz') then

             call fmod_md_condition__ensemble(NPTLZ)
             if (velocity_scaling_flag) then
                call parse_abort('Velocity scaling cannot be used in NPT.')
             endif
             getname='/thermostat/tau_Q'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauqtherm(value1)
             call parse_get('/thermostat/initialize',chara,lerr)
             if (.not. lerr .and. trim(chara) == 'yes') then
                initialize_thermostat_flag = .true.
             endif
             getname='/barostat/tau_Q'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauqbaro(value1)
             getname='/barostat/tau_W'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauwpres(value1)
             call parse_get('/barostat/initialize',chara,lerr)
             if (.not. lerr .and. trim(chara) == 'yes') then
                initialize_barostat_flag = .true.
             endif
!
          elseif (trim(chara) == 'nlxlypzt') then

             call fmod_md_condition__ensemble(NLXLYPZT)
             if (velocity_scaling_flag) then
                call parse_abort('Velocity scaling cannot be used in NPT.')
             endif
             getname='/thermostat/tau_Q'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauqtherm(value1)
             call parse_get('/thermostat/initialize',chara,lerr)
             if (.not. lerr .and. trim(chara) == 'yes') then
                initialize_thermostat_flag = .true.
             endif
             getname='/barostat/tau_Q'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauqbaro(value1)
             getname='/barostat/tau_W'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauwpres(value1)
             call parse_get('/barostat/initialize',chara,lerr)
             if (.not. lerr .and. trim(chara) == 'yes') then
                initialize_barostat_flag = .true.
             endif

          elseif (trim(chara) == 'npxpypzt') then

             call fmod_md_condition__ensemble(NPXPYPZT)
             if (velocity_scaling_flag) then
                call parse_abort('Velocity scaling cannot be used in NPT.')
             endif
             getname='/thermostat/tau_Q'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauqtherm(value1)
             call parse_get('/thermostat/initialize',chara,lerr)
             if (.not. lerr .and. trim(chara) == 'yes') then
                initialize_thermostat_flag = .true.
             endif
             getname='/barostat/tau_Q'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauqbaro(value1)
             getname='/barostat/tau_W'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauwpres(value1)
             call parse_get('/barostat/initialize',chara,lerr)
             if (.not. lerr .and. trim(chara) == 'yes') then
                initialize_barostat_flag = .true.
             endif

          elseif (trim(chara) == 'nlxlylzt') then

             call fmod_md_condition__ensemble(NLXLYLZT)
             if (velocity_scaling_flag) then
                call parse_abort('Velocity scaling cannot be used in NPT.')
             endif
             getname='/thermostat/tau_Q'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauqtherm(value1)
             call parse_get('/thermostat/initialize',chara,lerr)
             if (.not. lerr .and. trim(chara) == 'yes') then
                initialize_thermostat_flag = .true.
             endif
             getname='/barostat/tau_Q'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauqbaro(value1)
             getname='/barostat/tau_W'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauwpres(value1)
             call parse_get('/barostat/initialize',chara,lerr)
             if (.not. lerr .and. trim(chara) == 'yes') then
                initialize_barostat_flag = .true.
             endif

          elseif (trim(chara) == 'nptlzxy') then

             call fmod_md_condition__ensemble(NPTLZxy)
             if (velocity_scaling_flag) then
                call parse_abort('Velocity scaling cannot be used in NPT.')
             endif
             getname='/thermostat/tau_Q'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauqtherm(value1)
             call parse_get('/thermostat/initialize',chara,lerr)
             if (.not. lerr .and. trim(chara) == 'yes') then
                initialize_thermostat_flag = .true.
             endif
             getname='/barostat/tau_Q'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauqbaro(value1)
             getname='/barostat/tau_W'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauwpres(value1)
             call parse_get('/barostat/initialize',chara,lerr)
             if (.not. lerr .and. trim(chara) == 'yes') then
                initialize_barostat_flag = .true.
             endif

          elseif (trim(chara) == 'npt_pr') then
             call fmod_md_condition__ensemble(NPT_PR)
             if (velocity_scaling_flag) then
                call parse_abort('Velocity scaling cannot be used in NPT.')
             endif
             getname='/thermostat/tau_Q'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauqtherm(value1)
             call parse_get('/thermostat/initialize',chara,lerr)
             if (.not. lerr .and. trim(chara) == 'yes') then
                initialize_thermostat_flag = .true.
             endif
             getname='/barostat/tau_Q'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauqbaro(value1)
             getname='/barostat/tau_W'
             call parse_get(trim(getname), value1,lerr)
             if (.not. lerr)  call fmod_tauwpres(value1)
             call parse_get('/barostat/initialize',chara,lerr)
             if (.not. lerr .and. trim(chara) == 'yes') then
                initialize_barostat_flag = .true.
             endif
             !NI (+ YA suggestion)-->
             ! For "npt_pr", read pressure before ensemble tag.
             call parse_get('/condition/pressure', value1, lerr)
             if (.not. lerr)  call fmod_syspres(value1)
             !<--NI
             !NI-->
             ! NtT ensemble (constant external stress)
             ! Get reference cell (if it exists in mdff)
             call parse_get('/barostat/lx', value1, lerr1)
             call parse_get('/barostat/ly', value2, lerr2)
             call parse_get('/barostat/lz', value3, lerr3)
             call parse_get('/barostat/alpha', value4, lerr4)
             call parse_get('/barostat/beta', value5, lerr5)
             call parse_get('/barostat/gamma', value6, lerr6)
             if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. lerr5 .or. lerr6) ) then
                call fmod_referece_box(value1,value2,value3,value4,value5,value6)
             endif
             ! Get external stress (if it exists in mdff)
             call parse_get('/barostat/stress1', value1, lerr1) !!xx
             call parse_get('/barostat/stress2', value2, lerr2) !!xy
             call parse_get('/barostat/stress3', value3, lerr3) !!xz
             call parse_get('/barostat/stress4', value4, lerr4) !!yy
             call parse_get('/barostat/stress5', value5, lerr5) !!yz
             call parse_get('/barostat/stress6', value6, lerr6) !!zz
             if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. lerr5 .or. lerr6) ) then
                call fmod_sigmaS(value1,value2,value3,value4,value5,value6)
             else
                sigmaS = 0.d0
             endif
             !<--NI
          elseif (trim(chara) == 'opt') then
             call fmod_md_condition__ensemble(OPT)
             if (velocity_scaling_flag) then
                call parse_abort('Velocity scaling cannot be used in OPT.')
             endif
             getname='/condition/optimize/step_length'
             call parse_get(trim(getname),value1,lerr)
             value1 = value1 * 1.0d-10
             if (.not. lerr) then
                call fmod_opt_condition__step_length(value1)
             endif
             getname='/condition/optimize/convergence'
             call parse_get(trim(getname),value1,lerr)
             if (.not. lerr) then
                call fmod_opt_condition__convergence(value1)
             endif
             getname='/condition/optimize/up_rate'
             call parse_get(trim(getname),value1,lerr)
             if (.not. lerr) then
                call fmod_opt_condition__step_lurate(value1)
             endif
             getname='/condition/optimize/down_rate'
             call parse_get(trim(getname),value1,lerr)
             if (.not. lerr) then
                call fmod_opt_condition__step_ldrate(value1)
             endif
          endif
          call parse_get('/condition/steps', ivalue1, lerr)
          if (.not. lerr)  call fmod_howmany_steps(ivalue1)
          call parse_get('/condition/temperature', value1, lerr)
          if (.not. lerr)  call fmod_systemp(value1)
          call parse_get('/condition/pressure', value1, lerr)
          if (.not. lerr)  call fmod_syspres(value1)
          call parse_get('/condition/respa/style', chara, lerr)
          if (.not. lerr .and. trim(chara) == 'XO') then
             XORESPA=1; XIRESPA=0
          elseif (.not. lerr .and. trim(chara) == 'XI') then
             XORESPA=0; XIRESPA=1
          endif
! TZYBEGIN
          call parse_get('/condition/respa/nstep_skip_middle', ivalue1, lerr)
         !  if (.not. lerr)  call fmod_set_maxMTm(ivalue1)
          if (.not. lerr) then
              call fmod_set_maxMTm(ivalue1)
              if( myrank == 0 ) write(*,*) "nstep_skip_middle is deprecated. Use nshort_per_middle instead"
          endif
          call parse_get('/condition/respa/nstep_skip_long', ivalue1, lerr)
         !  if (.not. lerr)  call fmod_set_maxMTl(ivalue1)
          if (.not. lerr) then
              call fmod_set_maxMTl(ivalue1)
              if( myrank == 0 ) write(*,*) "nstep_skip_long is deprecated. Use nmiddle_per_long instead"
          endif

          call parse_get('/condition/respa/nshort_per_middle', ivalue1, lerr)
          if (.not. lerr)  call fmod_set_maxMTm(ivalue1)
          call parse_get('/condition/respa/nmiddle_per_long', ivalue1, lerr)
          if (.not. lerr)  call fmod_set_maxMTl(ivalue1)
! TZYEND

          call parse_get('/condition/respa/Pref_inner', value1, lerr)
          if (.not. lerr)  call fmod_set_Pinner(value1)
          call parse_get('/condition/respa/Pref_outer', value1, lerr)
          if (.not. lerr)  call fmod_set_Pouter(value1)
          call parse_get('/condition/respa/Pref_outermost', value1, lerr)
          if (.not. lerr)  call fmod_set_Poutermost(value1)
          call parse_get('/condition/maxwell_velocities', chara, lerr)
          if(trim(chara)=='yes') then
             call fmod_reset_maxwell(1)
          endif
          call parse_get('/mpi/division', chara, lerr)
          if (.not. lerr .and. trim(chara) == 'manual') then
             mpi_manual_division_flg = .true.
          endif
          call parse_get('/mpi/nxdiv', ivalue1, lerr)
          if (.not. lerr)  call set_npx(ivalue1)
          call parse_get('/mpi/nydiv', ivalue1, lerr)
          if (.not. lerr)  call set_npy(ivalue1)
          call parse_get('/mpi/nzdiv', ivalue1, lerr)
          if (.not. lerr)  call set_npz(ivalue1)
          call parse_get('/shake/maxiteration', ivalue1, lerr)
          if (.not. lerr)  call fmod_shake_max_iteration(ivalue1)
          call parse_get('/shake/shake_tolerance', value1, lerr)
          if (.not. lerr)  call fmod_md_shake__shake_tolerance(value1)
          call parse_get('/shake/rattle_tolerance', value1, lerr)
          if (.not. lerr)  call fmod_md_shake__rattle_tolerance(value1)
          call parse_get('/periodic/force/cutoff', value1, lerr)
          value1 = value1 * 1.0d-10
          value2 = value1 * value1
          if (.not. lerr) then
             call fmod_cutrad(value1)
             call fmod_cutrad2(value2)
          endif
          !
          call parse_get('/periodic/force/LJcorrection', chara, lerr)
          if(trim(chara)=='yes') then
             call fmod_LJ_correction(1)
          elseif(trim(chara)=='no') then
             call fmod_LJ_correction(0)
          elseif(lerr) then
             call fmod_LJ_correction(-1)
          endif
          !
          call parse_get('/periodic/type', chara, lerr)
          if(trim(chara) == 'cluster') then
             write(0,*) 'ERROR: type=cluster is not supported'
             call modylas_abort()
             call fmod_md_periodic__type(CLUSTER)
          elseif(trim(chara) == 'simple') then
             write(0,*) 'ERROR: type=simple is not supported'
             call modylas_abort()
             call fmod_md_periodic__type(SIMPLE)
          elseif(trim(chara) == 'ewald') then
             call fmod_md_periodic__type(EWALD)
             call parse_get('/periodic/ewald/alpha',value1,lerr)
             if (.not. lerr) then
               call fmod_md_ewald__alpha(value1)
               call check_pmewald_alpha(value1)
             endif
             call parse_get('/periodic/ewald/h2max',ivalue1,lerr)
             if (.not. lerr)  call fmod_md_ewald__max_h2(ivalue1)

          elseif(trim(chara) == 'pmewald') then

             call fmod_md_periodic__type(PMEWALD)
             call parse_get('/periodic/pme/alpha',value1,lerr)
             if (.not. lerr) then
                call fmod_md_ewald__alpha(value1)
                call check_pmewald_alpha(value1)
             endif
             call parse_get('/periodic/pme/nlevel',ivalue1,lerr)
             if (.not. lerr)    call fmod_set_fmm_nlevel(ivalue1)
             call parse_get('/periodic/pme/ncellx',ivalue1,lerr)
             if (.not. lerr)    call fmod_set_fmm_ncellx(ivalue1)
             call parse_get('/periodic/pme/ncelly',ivalue1,lerr)
             if (.not. lerr)    call fmod_set_fmm_ncelly(ivalue1)
             call parse_get('/periodic/pme/ncellz',ivalue1,lerr)
             if (.not. lerr)    call fmod_set_fmm_ncellz(ivalue1)
             call parse_get('/periodic/pme/ncell',ivalue1,lerr)
             if (.not. lerr)    call fmod_set_ncell(ivalue1)
             call parse_get('/periodic/pme/bsorder',value1,lerr)
             if (.not. lerr)  call fmod_pmewald_bsorder(value1)
             call parse_get('/periodic/pme/nfft1',ivalue1,lerr)
             if (.not. lerr)  call fmod_pmewald_nfft1(ivalue1)
             call parse_get('/periodic/pme/nfft2',ivalue2,lerr)
             if (.not. lerr)  call fmod_pmewald_nfft2(ivalue2)
             call parse_get('/periodic/pme/nfft3',ivalue3,lerr)
             if (.not. lerr)  call fmod_pmewald_nfft3(ivalue3)

          elseif(trim(chara) == 'cmm') then

             call fmod_md_periodic__type(FMM)
             call parse_get('/periodic/cmm/ndirect',ivalue1,lerr)
             if (.not. lerr)    call fmod_set_fmm_ndirect(ivalue1)
             call parse_get('/periodic/cmm/nmax',ivalue1,lerr)
             if (.not. lerr)    call fmod_set_nmax(ivalue1)
             call parse_get('/periodic/cmm/ULswitch',ivalue1,lerr)
             if (.not. lerr)    call fmod_set_lgflg(ivalue1)
             call parse_get('/periodic/cmm/sterm',chara,lerr)
             if(trim(chara)=='yes') then
                call fmod_ewald_sterm(1)
             elseif(trim(chara)=='no'.or.lerr) then
                call fmod_ewald_sterm(0)
             endif

             call parse_get('/periodic/cmm/nlevel',ivalue1,lerr)
             if (.not. lerr)    call fmod_set_fmm_nlevel(ivalue1)
             call parse_get('/periodic/cmm/ncellx',ivalue1,lerr)
             if (.not. lerr)    call fmod_set_fmm_ncellx(ivalue1)
             call parse_get('/periodic/cmm/ncelly',ivalue1,lerr)
             if (.not. lerr)    call fmod_set_fmm_ncelly(ivalue1)
             call parse_get('/periodic/cmm/ncellz',ivalue1,lerr)
             if (.not. lerr)    call fmod_set_fmm_ncellz(ivalue1)
             call parse_get('/periodic/cmm/ncell',ivalue1,lerr)
             if (.not. lerr) then
                call warn_old_ncell
                call fmod_set_ncell(ivalue1)
                call fmod_set_fmm_nlevel( calculate_fmm_nlevel(ivalue1) )
             endif
             do while(.true.)
                call parse_get('/periodic/cmm/mincell',ivalue1,lerr)
                if(lerr) exit
             enddo
          endif

          !ay     hidden command
          lerr_bk = lerr
          call parse_get('/forcefield/CMAPVersion',ivalue1,lerr)
          if(.not.lerr) then
             call fmod_set_CMAPVersion(ivalue1)
          endif
          lerr = lerr_bk

          !ya     hidden command
          call parse_get('/debug/0step_stop',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             call fmod_set_0step_stop
          endif

          !ya     hidden command
          call parse_get('/allocation parameters/totnconstL',ivalue1,lerr)
          if (.not. lerr) then
             call fmod_set_totnconstL(ivalue1)
          endif
          call parse_get('/allocation parameters/na1cell',ivalue1,lerr)
          if (.not. lerr) then
             call fmod_set_na1cell(ivalue1)
          endif
          call parse_get('/allocation parameters/max_nsegments_per_cell', &
     &                    ivalue1,lerr)
          if (.not. lerr) then
             call fmod_set_maxnsegmentspercell(ivalue1)
          endif

#ifdef SEGSHAKE
          call parse_get('/COM/constrain_COM',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             call fmod_constrain_COM(1)
          else
             call fmod_constrain_COM(0)
          endif
          call parse_get('/COM/groupAtop', ivalue1, lerr1)
          call parse_get('/COM/groupAend', ivalue2, lerr2)
          if(.not. lerr1 .and. .not. lerr2) then
             call fmod_set_groupA(ivalue1, ivalue2)
          endif
          call parse_get('/COM/groupBtop', ivalue1, lerr1)
          call parse_get('/COM/groupBend', ivalue2, lerr2)
          if(.not. lerr1 .and. .not. lerr2) then
             call fmod_set_groupB(ivalue1, ivalue2)
          endif
          call parse_get('/COM/dist_COM', qvalue1, lerr1)
          if(.not. lerr1) then
             value1=qvalue1
             call fmod_set_dist_COM(value1)
          endif
          call parse_get('/COM/change_distCOM',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             call fmod_change_distCOM(1)
          endif
          call parse_get('/COM/deltaR', qvalue1, lerr1)
          if(.not. lerr1) then
             value1=qvalue1
             call fmod_set_d_COMR(value1)
          endif
#endif

#ifdef TIP4
        call parse_get('/tip4p/type',chara,lerr)
        if (.not. lerr)then
        if(trim(chara) == 'original') then
          call fmod_set_tip4p_geometry(0)
        elseif(trim(chara) == '2005') then
          call fmod_set_tip4p_geometry(1)
        else
          write(*,*) 'keyword in <tip4p> type= </tip4p> is wrong', &
     &    ', it should be original or 2005.'
          call mpiend
        endif
        endif
!
        call parse_get('/tip4p/msite',chara,lerr)
        if (.not. lerr)then
          if(trim(chara) == 'modylas') then
            call fmod_set_msite_position(0)
          elseif(trim(chara) == 'gromacs') then
            call fmod_set_msite_position(1)
          else
            write(*,*) 'keyword in <tip4p> msite= </tip4p> is wrong', &
     &      ', it should be modylas or gromacs.'
            call mpiend
          endif
        endif
#endif 

!
! closing .mddef
!
          call parse_close(lerr)
          if(myrank==0) write(*,*) 'Reading .mddef file ended successfully!'


          if (lerr) then
             call parse_abort(parse_show_error_message())
          endif
       endif
    endif
    !     Show error message except parse_close().
    if (lerr) then
       call parse_abort(parse_show_error_message())
    endif
  end subroutine parse_input_v0_9_0
!-------------------------------------------------------------------------
  subroutine set_default_keyword_0_9_0
    use file_utility, only : ascii_output
    implicit none

    ! /output/ascii=no
    ascii_output=.false.

  end subroutine set_default_keyword_0_9_0
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to abort reading inputs.
!! \author Kensuke Iwahashi
!<
  subroutine parse_abort(string)
    use mpi_tool
    implicit none
    character(len=*), intent(in) :: string
#ifdef MPI
    include 'mpif.h'
    integer(4) :: ierr, myrank
    call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)
    if(myrank == 0 ) then
       write(0,*) trim(string)
    endif
#else /* MPI */
    write(0,*) trim(string)
#endif /* MPI */
    call modylas_abort()
  end subroutine parse_abort
!-------------------------------------------------------------------------
  subroutine alloc
    !fujimoto
    use atom_mass
    use lj_mod
    use coulomb_mod
    use param
    use mpi_tool
    implicit none
    allocate(mass(npara))         ; mass=0d0
    allocate(r_mass(npara))       ; r_mass=0d0
    allocate(chgv(npara))         ; chgv=0d0
    allocate(epsilon_sqrt(npara)) ; epsilon_sqrt=0d0
    allocate(R_half(npara))       ; R_half=0d0
  end subroutine alloc
!-------------------------------------------------------------------------

!------------------------------------------------------------------------
!     Read input for ver. 1.0.0 (2013.10.24)
!------------------------------------------------------------------------
!>
!! \brief  Subroutine for reading input file of version 1.0.0
!! \author  Atsushi Yamada
!<
!-------------------------------------------------------------------------
  subroutine parse_input_v1_0_0(read_mddef)
    use parse
    use file_mdmntr
    use file_mdtrj
    use file_dcd
    use file_xtc
    use file_restart
    use device_numbers
    use session_name_mod
    use file_utility
    use md_const
    use md_condition
    use ensemble_numbers
    use mpi_3d_grid, only : set_npx, set_npy, set_npz, mpi_manual_division_flg
    use md_multiplestep
    use thermostat
    use param
    use mol_info
    use md_periodic
    use shake_rattle_roll
    use mpi_tool
    use bond
    use bond_morse
    use bond_morse2
    !use bond_pcff
    use angle
    use angle_morse
    !use angle_pcff
    !use bond_bond
    !use parse_angle_pcff
    !use bond_bond
    !use parse_bb 
    !use bond_angle
    !use parse_ba
    !use angle_angle
    !use parse_angle_angle
    use UB
    use dihedral
    use improper_torsion
    use CMAP
    use void123
    use special14
    use opt_mod
    use position_constrain
    use comm_pme_mod
    use pme_far
    use ewald_variables
    use ewald_mod
    use unit_cell
    use segments
    use trajectory_org
    use extended_system
    use coulomb_mod, only : fmod_set_mol_chgv, read_chgv
    use lj_mod, only : fmod_set_mol_ljs__epsilon_sqrt, fmod_set_mol_ljs__r_half, &
         & read_epsilon_sqrt, read_r_half
    use molecules, only : fmod_alloc_atom2mol, fmod_alloc_mol2atom
    use atom_mass, only : fmod_set_mol_mass, read_mass
    use fmm_parameters, only : fmod_set_nmax
    use fmm_far, only : fmod_set_lgflg, fmod_set_fmm_nlevel
    use parse_shake, only : fmod_alloc_yyparse_shake
    use cutoff_radius, only : fmod_cutrad, fmod_cutrad2
    use force_field_numbers
    use maxwell_distribution, only : fmod_reset_maxwell
    use file_force, only : fmod_force_start, fmod_force_interval
    use file_mdff, only : read_forcefield, read_md_shake, read_molecule, read_parameter_number
    use file_mdxyz_bin, only : read_mdxyzbin
    use md_oplsaa_special_divide, only : fmod_md_oplsaa_lj_sp_divide, fmod_md_oplsaa_cl_sp_divide, &
        & fmod_md_oplsaa_lj_sp_divide_mddef, fmod_md_oplsaa_coulomb_sp_divide_mddef, &
        & fmod_md_check_scaling_factor_is_set
    use system_dipole, only: set_up_system_dipole
#ifdef GROEXT
    use file_groext
#endif
#ifdef TIP4
    use tip4p, only : fmod_set_tip4p_geometry, fmod_set_msite_position
#endif
    use trajectory_mpi
    use domain
    use center_of_mass_variables
    implicit none
    integer(4) :: i, j, k, io, im
    integer(4) :: natom, nthermostat, nbarostat
    integer(4) :: ivalue1, ivalue2, ivalue3
    integer(4) :: ivalue4, ivalue5, ivalue6
    real(8)    :: value1, value2, value3, value4, value5
    real(8)    :: value6, value7, value8, value9
    real(8)    :: value10, value11
    real(8)   :: qvalue1, qvalue2, qvalue3, qvalue4, qvalue5
    real(8)   :: qvalue6, qvalue7, qvalue8, qvalue9, qvalue10
    real(8)   :: qvalue11
    character(LEN=100) :: getname,segname,name_TopPar_Spc
    character(LEN=10) :: chara
    logical :: read_mddef
    logical :: lerr, lerr1, lerr2, lerr3, lerr4, lerr5
    logical :: lerr6, lerr7, lerr8, lerr9, lerr_bk
    logical :: lerr10, lerr11, lerr12, lerr13, lerr14
    logical :: velocity_scaling_flag = .false.
    logical :: velocity_scaling_region_wise_flag = .false.
    logical :: volume_scaling_flag = .false.
    logical :: cellgeometry_scaling_flag = .false.
    integer(4), allocatable :: natommol(:), nmolmol(:), molnseg(:)
    integer(4) :: hws, sid, nas, nall, temp, nvc, nvl
    integer(4) :: nsl, nsct, nbond, nag, nub, nit, ncm, ndi
    integer(4) :: dummyInt

#ifdef GROEXT
    if(gmx_input_groext) then
       call convert_gmx_to_modylas()
    end if
    if(gmx_output_groext) then
       call convert_modylas_to_gmx()
    end if
    if(gmx_convert_only) then
  call mpiend()
    end if
#endif

    !
    !     mdxyz
    !
#ifdef GROEXT
    if (.not. gmx_input_groext .and. can_read_binary_file(trim(session_name)//'.mdxyz.bin')) then
#else
    if (can_read_binary_file(trim(session_name)//'.mdxyz.bin')) then       
#endif
       if(myrank==0) write(*,*) 'Reading mdxyz.bin file started.'
       call read_mdxyzbin()
       if(myrank==0) write(*,*) 'Reading mdxyz.bin file ended successfully!'

    else
#ifdef GROEXT
       if (gmx_input_groext) then
          call parse_open_string(str_mdxyz, lerr)
       else
          call parse_open(trim(session_name)//'.mdxyz', lerr)
       end if
#else
       call parse_open(trim(session_name)//'.mdxyz', lerr)
#endif
       if(myrank==0) write(*,*) 'Reading .mdxyz file started.'
       if (lerr) then
          write(0,*) parse_show_error_message()
          write(0,*) 'ERROR in mdxyz file.'
          call modylas_abort()
       else
          call parse_get('/atom/natom',natom,lerr)
          if (.not. lerr) then
             call fmod_n(natom)
             if(natom .ne. 0) then
                call alloc_i2m(natom)
                call fmod_alloc_xyz(natom)
                call fmod_alloc_v(natom)
             endif
          endif
          ! positions
          do i = 1, natom
             call parse_get('/atom/positions',qvalue1,lerr1)
             call parse_get('/atom/positions',qvalue2,lerr2)
             call parse_get('/atom/positions',qvalue3,lerr3)
             if (.not. (lerr1 .or. lerr2 .or. lerr3)) then
                value1 = qvalue1 * 1.0d-10
                value2 = qvalue2 * 1.0d-10
                value3 = qvalue3 * 1.0d-10
                call fmod_set_xyz(i-1,value1,value2,value3,1)
             endif
          enddo
          ! velocities
          getname='/atom/velocities'
          do i = 1, natom
             call parse_get(trim(getname), qvalue1,lerr1)
             call parse_get(trim(getname), qvalue2,lerr2)
             call parse_get(trim(getname), qvalue3,lerr3)
             if (.not. (lerr1 .or. lerr2 .or. lerr3)) then
                value1 = qvalue1 * 1.0d-10
                value2 = qvalue2 * 1.0d-10
                value3 = qvalue3 * 1.0d-10
                call fmod_set_v(i-1,value1,value2,value3,1)
             endif
          enddo
          call parse_get('/thermostat/nthermostat', nthermostat,lerr)
          do i = 1, nthermostat
             call parse_get('/thermostat/positions',value1,lerr)
             if (.not. lerr)  call fmod_set_rss(i-1, value1, 1)
          enddo
          do i =1, nthermostat
             call parse_get('/thermostat/velocities',value1,lerr)
             if (.not. lerr)  call fmod_set_vss(i-1, value1, 1)
          enddo
          call parse_get('/barostat/nbarostat', nbarostat,lerr)
          do i = 1, nbarostat
             call parse_get('/barostat/positions',value1,lerr)
             if (.not. lerr)  call fmod_set_rssb(i-1, value1, 1)
          enddo
          do i = 1, nbarostat
             call parse_get('/barostat/velocities',value1,lerr)
             if (.not. lerr)  call fmod_set_vssb(i-1, value1, 1)
          enddo
          call parse_get('/periodic cell/length',qvalue1,lerr)  ! x
          if (.not. lerr) then
             qvalue1 = qvalue1 * 1.0d-10
             value1 = qvalue1
             value2 = qvalue1 * 0.5d0
             call fmod_cellx(value1)
             call fmod_cellxh(value2)
          endif
          call parse_get('/periodic cell/length',qvalue1,lerr)  ! y
          if (.not. lerr) then
             qvalue1 = qvalue1 * 1.0d-10
             value1 = qvalue1
             value2 = qvalue1 * 0.5d0
             call fmod_celly(value1)
             call fmod_cellyh(value2)
          endif
          call parse_get('/periodic cell/length',qvalue1,lerr)  ! z
          if (.not. lerr) then
             qvalue1 = qvalue1 * 1.0d-10
             value1 = qvalue1
             value2 = qvalue1 * 0.5d0
             call fmod_cellz(value1)
             call fmod_cellzh(value2)
          endif
          call parse_get('/periodic cell/angle',value1,lerr) ! alpha
          if (.not. lerr)  call fmod_set_alpha(value1)
          call parse_get('/periodic cell/angle',value1,lerr) ! beta
          if (.not. lerr)  call fmod_set_beta(value1)
          call parse_get('/periodic cell/angle',value1,lerr) ! gamma
          if (.not. lerr)  call fmod_set_gamma(value1)
          call parse_get('/periodic cell/vboxg',value1,lerr1)
          call parse_get('/periodic cell/vboxg',value2,lerr2)
          call parse_get('/periodic cell/vboxg',value3,lerr3)
          call parse_get('/periodic cell/vboxg',value4,lerr4)
          call parse_get('/periodic cell/vboxg',value5,lerr5)
          call parse_get('/periodic cell/vboxg',value6,lerr6)
          call parse_get('/periodic cell/vboxg',value7,lerr7)
          call parse_get('/periodic cell/vboxg',value8,lerr8)
          call parse_get('/periodic cell/vboxg',value9,lerr9)
          if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. &
               &        lerr5 .or. lerr6 .or. lerr7 .or. lerr8 .or. lerr9)) then
             call fmod_set_vboxg(value1,value2,value3,value4,value5, &
                  &                         value6,value7,value8,value9)
          endif
          call parse_close(lerr)
          if(myrank==0) write(*,*) 'Reading .mdxyz file ended successfully!'

       endif
       if (lerr) then
          call parse_abort(parse_show_error_message())
       endif
    endif


    !
    !     mdff
    !
    if(.not.can_read_binary_file(trim(session_name)//'.mdff.bin'))then
!      if(myrank==0)write(6,*) "No ", trim(session_name), ".mdff.bin"
       goto 100 !GOTO READING MDFF FILE (ascii format)
    endif

    !===  mdff.bin -- read binary format ===
    if(myrank==0) then
       open(f_mdff, file=trim(session_name)//'.mdff.bin', iostat=io, &
            &       status='old', access='sequential', form='unformatted')
    endif
    if(myrank==0) write(*,*) 'Reading .mdff file started.'

    call read_forcefield
    call read_parameter_number
    call alloc
    call read_segment
    call read_molecule
    call read_mass
    call read_chgv
    call read_epsilon_sqrt
    call read_R_half
    call read_md_shake
    call read_md_void    ! void (Common)
    call read_md_LJ_special ! special (LJ)
    call read_md_coulomb_special  ! special (CL)
    call read_md_charmm_a_bond
    call read_md_bond_morse
    call read_md_bond_morse2
    !call read_md_pcff_a_bond
    if(myrank==0) read(f_mdff) dummyInt
    call read_md_charmm_angleA
    call read_md_charmm_angleB
    call read_md_angleA_morse
    call read_md_angleB_morse
    !call read_md_pcff_angleA
    if(myrank==0) read(f_mdff) dummyInt
    !call read_md_pcff_angleB
    if(myrank==0) read(f_mdff) dummyInt
    call read_md_charmm_a_ub
    !call read_md_pcff_bbA
    if(myrank==0) read(f_mdff) dummyInt
    !call read_md_pcff_bbB
    if(myrank==0) read(f_mdff) dummyInt
    !call read_md_pcff_baA
    if(myrank==0) read(f_mdff) dummyInt
    !call read_md_pcff_baB
    if(myrank==0) read(f_mdff) dummyInt
    !call read_md_pcff_aaA
    if(myrank==0) read(f_mdff) dummyInt
    !call read_md_pcff_aaB
    if(myrank==0) read(f_mdff) dummyInt
    call read_md_charmm_dihedralA
    call read_md_charmm_dihedralB
    call read_md_charmm_CMAPA
    call read_md_charmm_CMAPB
    call read_md_charmm_CMAPC
    call read_md_charmm_CMAPD
    call read_md_charmm_CMAPE
    call read_md_charmm_itorsionA
    call read_md_charmm_itorsionB
    !ay+
    !      call read_position_constrain   !NOT support from binary input
    !ay+
    call read_CMAPVersion
    close(f_mdff)

    if(myrank==0) then
       write(*,*) 'Reading .mdff.bin file ended successfully!'
       !       For debug force field.
       !       call write_memory_mdff
    endif

    goto 200 !GOTO READING MDDEF FILE

    !===  mdff -- read ascii format ver.1.0.0 ===
100 continue
    if(myrank==0) write(*,*) 'Reading .mdff file started.'

#ifdef GROEXT
    if (gmx_input_groext) then
       call parse_open_string(str_mdff, lerr)
    else
       call parse_open(trim(session_name)//'.mdff', lerr)
    end if
#else
    call parse_open(trim(session_name)//'.mdff', lerr)
#endif

    if (lerr) then
       write(0,*) parse_show_error_message()
       write(0,*) 'ERROR in mdff file.'
       call modylas_abort()
    else
       call parse_get('/forcefield/type',chara,lerr)
       if(trim(chara)=='charmm')then
#ifdef OPLSAMBER
          if(myrank==0) then
             write(0,*) 'ERROR: -DOPLSAMBER is set ', &
               &     'in -DCMAKE_Fortran_FLAGS=" " for cmake. ',  &
               &     'It must be eliminated to handle input with type=charmm in .mdff file'
!                 &   'in src/Makefile incompatible with type=charmm'
             call modylas_abort()
          endif
#endif
          md_condition__force_field=100
          call parse_get('/forcefield/cmap_version',ivalue1,lerr)
          if (lerr) then
             if (myrank == 0) then
                write(*,*) "###########################################"
                write(*,*) "CMAP Version is not set in your mdff file."
                write(*,*) "So, CMAP Version was set to be 36."
                write(*,*) "If you do not use CMAP Version 36,"
                write(*,*) "please write key word 'CMAPVersion=22'"
                write(*,*) "between <forcefiled> and </forcefiled> "
                write(*,*) "in your mdff file as follws."
                write(*,*)
                write(*,*) "<forcefiled>"
                write(*,*) "  type=charmm"
                write(*,*) "  cmap_version=22"
                write(*,*) "</forcefiled>"
                write(*,*) "###########################################"
             endif
          else
             call fmod_set_CMAPVersion(ivalue1)
          endif
       elseif(trim(chara)=='oplsaa')then
#ifndef OPLSAMBER
          if(myrank==0) then
             write(0,*) 'ERROR: -DOPLSAMBER is necessary ', &
               &     'in -DCMAKE_Fortran_FLAGS=" " for cmake, ',  &
               &     'to handle input with type=oplsaa in .mdff file'
             call modylas_abort()
          endif
#endif
          md_condition__force_field=OPLSAA ! OPLS
          call parse_get('/forcefield/special_divide_lj',value1,lerr)
          if (.not. lerr)  call fmod_md_oplsaa_lj_sp_divide(value1)
          call parse_get('/forcefield/special_divide_coulomb',value1,lerr)
          if (.not. lerr)  call fmod_md_oplsaa_cl_sp_divide(value1)
       elseif(trim(chara)=='amber')then
#ifndef OPLSAMBER
          if(myrank==0) then
             write(0,*) 'ERROR: -DOPLSAMBER is necessary ', &
               &     'in -DCMAKE_Fortran_FLAGS=" " for cmake, ',  &
               &     'to handle input with type=amber in .mdff file'
             call modylas_abort()
          endif
#endif
          md_condition__force_field=AMBER ! AMBER
          call parse_get('/forcefield/special_divide_lj',value1,lerr)
          if (.not. lerr)  call fmod_md_oplsaa_lj_sp_divide(value1)
          call parse_get('/forcefield/special_divide_coulomb',value1,lerr)
          if (.not. lerr)  call fmod_md_oplsaa_cl_sp_divide(value1)
!gaff
       elseif(trim(chara)=='gaff')then
#ifndef OPLSAMBER
          if(myrank==0) then
             write(0,*) 'ERROR: -DOPLSAMBER is necessary ', &
               &     'in -DCMAKE_Fortran_FLAGS=" " for cmake, ',  &
               &     'to handle input with type=gaff in .mdff file'
#ifndef GAFF
             write(0,*) 'ERROR: -DGAFF is necessary ', &
               &     'in -DCMAKE_Fortran_FLAGS=" " for cmake, ',  &
               &     'to handle input with type=gaff in .mdff file'
             call modylas_abort()
#endif
             call modylas_abort()
          endif
#endif
          md_condition__force_field=GAmbFF ! General Amber FF (GAFF)
          call parse_get('/forcefield/special_divide_lj',value1,lerr)
          if (.not. lerr)  call fmod_md_oplsaa_lj_sp_divide(value1)
          call parse_get('/forcefield/special_divide_coulomb',value1,lerr)
          if (.not. lerr)  call fmod_md_oplsaa_cl_sp_divide(value1)
       elseif(trim(chara)=='kremer')then
#ifdef OPLSAMBER
          if(myrank==0) then
             write(0,*) 'ERROR: -DOPLSAMBER is set ', &
     &     'in -DCMAKE_Fortran_FLAGS=" " for cmake. ',  &
     &     'It must be eliminated to handle input with type=kremer in .mdff file'
             call modylas_abort()
          endif
#endif
          md_condition__force_field=KREMER ! Kremer's potential
       !elseif(trim(chara) == 'pcff') then
       !   md_condition__force_field=PCFF
       endif  !! trim(chara)

       call parse_get('/topology and parameters/nspecies',nmol,lerr)
       call fmod_alloc_molinfo(nmol)
       if(.not. lerr) allocate(natommol(nmol))
       if(.not. lerr) natommol=0
       if(.not. lerr) allocate(nmolmol(nmol))
       if(.not. lerr) nmolmol=0
       do i = 1, nmol
          call parse_get('/system',ivalue1, lerr)
          call parse_get('/system',ivalue2, lerr)
          if(.not. lerr) nmolmol(i) = ivalue2
       enddo

       npara=0
       nall =0
       name_TopPar_Spc='/topology and parameters/species'
       do i = 1, nmol
          call parse_specify_key(trim(name_TopPar_Spc),i)
          !         '/topology and parameters/species/id' is not used now.
          call parse_get(trim(name_TopPar_Spc)//'/id',ivalue1,lerr)
          call parse_get(trim(name_TopPar_Spc)//'/natom',ivalue1,lerr)
          if(.not. lerr) then
             call fmod_set_molinfo(i-1,nmolmol(i),ivalue1)
             natommol(i) = ivalue1
          endif
          npara = npara + natommol(i)
          nall = nall + natommol(i)*nmolmol(i)
       enddo
       call fmod_set_param_num(npara,nmol,nmolmol,natommol)
       call alloc
       call fmod_alloc_atom2mol(nall)
       call fmod_alloc_mol2atom(nall)
       call fmod_alloc_irearrange(nall)

       hws = 0
       allocate(molnseg(nmol)) ; molnseg=0
       do i = 1, nmol
          call parse_specify_key('/topology and parameters/species',i)
          getname = trim(name_TopPar_Spc)//'/segments/nsegment'
          call parse_get(trim(getname),molnseg(i),lerr)
          if (.not. lerr) then
             hws = hws + molnseg(i) * nmolmol(i)
          endif
       enddo
       if(hws .ne. 0) then
          nsegments = hws
          call fmod_alloc_seg_natoms(nsegments)
          call fmod_alloc_seg_center_x(nsegments)
          call fmod_alloc_seg_center_y(nsegments)
          call fmod_alloc_seg_center_z(nsegments)
          call fmod_alloc_segtop(nsegments)
       endif

       do i = 1, nmol
          call fmod_set_mol_ID(i-1)
          call parse_specify_key('/topology and parameters/species',i)
          call fmod_set_mol_seg_num(molnseg(i))
          segname=trim(name_TopPar_Spc)//'/segments/segment'
          do j = 1, molnseg(i)
             call parse_specify_key(segname,j)
             getname=trim(name_TopPar_Spc)//'/segments/segment/ID'
             call parse_get(trim(getname),sid,lerr1)
             if(lerr1) then
                getname=trim(name_TopPar_Spc)//'/segments/segment/id' !(small letter is better)
                call parse_get(trim(getname),sid,lerr1)
             endif
             getname=trim(name_TopPar_Spc)//'/segments/segment/natom'
             call parse_get(trim(getname),nas,lerr2)
             if (.not. (lerr1 .or. lerr2)) then
                call fmod_set_mol_seg_natom(sid, nas)
             endif
             do k = 1, nas
                getname=trim(name_TopPar_Spc)//'/segments/segment/atom'
                call parse_get(trim(getname),ivalue1,lerr)
                if (.not. (lerr1 .or. lerr)) then
                   call fmod_set_mol_seg2atom(sid,k-1,ivalue1)
                endif
             enddo
          enddo
          call fmod_set_mol_mol_natoms()
          call fmod_set_mol_mol2atom()
       enddo


       do i = 1, nmol
          call parse_specify_key('/topology and parameters/species',i)
          !READING MASS AND R_MASS
          do j = 1, natommol(i)
             getname=trim(name_TopPar_Spc)//'/mass'
             call parse_get(trim(getname),value1,lerr)
             if (.not. lerr)  call fmod_set_mol_mass(value1)
          enddo
          !READING CHARGE
          do j = 1, natommol(i)
             getname=trim(name_TopPar_Spc)//'/charge'
             call parse_get(trim(getname),value1,lerr)
             if (.not. lerr)  call fmod_set_mol_chgv(value1)
          enddo
          !READING LJ EPSILON
          do j = 1, natommol(i)
             getname=trim(name_TopPar_Spc)//'/epsilon'
             call parse_get(trim(getname),qvalue1,lerr)
             if (.not. lerr) then
                !convert unit:[kcal/mol] => [J]
                value1 = sqrt(qvalue1*1000.0d0 * md_CALORIE / md_AVOGADRO)
                call fmod_set_mol_ljs__epsilon_sqrt(value1)
             endif
          enddo
          !READING LJ R-HALF
          do j = 1, natommol(i)
             getname=trim(name_TopPar_Spc)//'/r'
             call parse_get(trim(getname),qvalue1,lerr)
             if (.not. lerr) then
                if( md_condition__force_field == CHARMM .or. &
                 &  md_condition__force_field == KREMER  )then
                   !convert unit[A] => [m]
                   value1 = qvalue1 * 1.0d-10 * 0.5d0
                   call fmod_set_mol_ljs__r_half(value1)
                elseif(md_condition__force_field == OPLSAA .or. &
                 &     md_condition__force_field == AMBER  .or. &
                 &     md_condition__force_field == GAmbFF )then
                   !convert unit[A] => [m]
                   qvalue1 = qvalue1 * 1.0d-10 * 0.5d0
                   value1 = qvalue1 * 2.0d0 ** (1.0d0/6.0d0)
                   call fmod_set_mol_ljs__r_half(value1)
                !elseif(md_condition__force_field == PCFF  )then
                !   !convert unit[A] => [m]
                !   value1 = qvalue1*1e-10
                !   call fmod_set_mol_ljs__r_half(value1)
                endif
             endif
          enddo
       enddo

       call fmod_alloc_yyparse_shake(npara)
       getname = trim(name_TopPar_Spc)//'/shake pair'
       do i = 1, nmol
          call parse_specify_key('/topology and parameters/species',i)
          temp = 0
          do while(.true.)
             call parse_get(trim(getname), ivalue1,lerr)
             if(lerr) exit
             call parse_get(trim(getname), ivalue2,lerr)
             call parse_get(trim(getname), qvalue1,lerr)
             value1 = qvalue1 * 1.0d-10
             call add_to_mol_shake_list(i,ivalue1,ivalue2,value1)
             temp = temp + 1
          enddo
          totnconst = totnconst + temp * nmolmol(i)
       enddo
       call set_shake()

       ! READING VOID PAIR
       call fmod_alloc_yyparse_lj_top(npara)
       call fmod_alloc_yyparse_lj_n(npara)
       call fmod_alloc_yyparse_coulomb_top(npara)
       call fmod_alloc_yyparse_coulomb_n(npara)
       do i =1, nmol
          !! LJ void !!
          call parse_specify_key('/topology and parameters/species',i)
          getname = trim(name_TopPar_Spc)//'/nvoidpair_lj'
          call parse_get(trim(getname), nvl, lerr)
          getname = trim(name_TopPar_Spc)//'/lj void pair'
          do j = 1, nvl
             call parse_get(trim(getname), ivalue1,lerr1)
             call parse_get(trim(getname), ivalue2,lerr2)
             if (.not. (lerr1 .or. lerr2)) then
                call add_to_mol_lj_void(i, ivalue1, ivalue2)
             endif
          enddo
          !! COULOMB VOID PAIR
          getname = trim(name_TopPar_Spc)//'/nvoidpair_coulomb'
          call parse_get(trim(getname),nvc, lerr)
          getname = trim(name_TopPar_Spc)//'/coulomb void pair'
          do j = 1, nvc
             call parse_get(trim(getname),ivalue1, lerr1)
             call parse_get(trim(getname),ivalue2, lerr2)
             if (.not. (lerr1 .or. lerr2)) then
                call add_to_mol_coulomb_void(i,ivalue1,ivalue2)
             endif
          enddo
       enddo
       call set_md_void()

       !READING SPECAIL PAIR
       do i = 1, nmol
          !! LJ special PAIR !!
          call parse_specify_key('/topology and parameters/species',i)
          getname = trim(name_TopPar_Spc)//'/nspecialpair_lj'
          call parse_get(trim(getname), nsl, lerr)
          getname = trim(name_TopPar_Spc)//'/lj special pair'
          do j = 1, nsl
             call parse_get(trim(getname),ivalue1,lerr1)
             call parse_get(trim(getname),ivalue2,lerr2)
             call parse_get(trim(getname),qvalue1,lerr3)
             call parse_get(trim(getname),qvalue2,lerr4)
             if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4)) then
                if( md_condition__force_field == CHARMM .or. md_condition__force_field == KREMER )then
                   !converting unit: [kcal/mol] => [J]
                   !converting unit: [A] => [m]
                   !epsilon needs not to be square-rooted and R not halfed
                   value1 = qvalue1 * 1000d0 * md_CALORIE / md_AVOGADRO
                   value2 = qvalue2 * 1.0d-10
                   call add_to_mol_lj_list(i,ivalue1,ivalue2,value1,value2)
                elseif( md_condition__force_field == OPLSAA .or. &
                   &    md_condition__force_field == AMBER  .or. &
                   &    md_condition__force_field == GAmbFF )then
                   continue !! nothing to do in OPLS/AMBER
                endif
             endif
          enddo
          !! COUBLOMB special PAIR !!
          getname = trim(name_TopPar_Spc)//'/nspecialpair_coulomb'
          call parse_get(trim(getname),nsct,lerr)
          getname = trim(name_TopPar_Spc)//'/coulomb special pair'
          do j = 1, nsct
             call parse_get(trim(getname),ivalue1,lerr1)
             call parse_get(trim(getname),ivalue2,lerr2)
             call parse_get(trim(getname),value1,lerr3)
             if (.not. (lerr1 .or. lerr2 .or. lerr3)) then
                call add_to_mol_coulomb_list(i,ivalue1,ivalue2,value1)
             endif
          enddo
       enddo
       call set_md_lj_special()
       call set_md_coulomb_special()


       !READING BOND PAIR
       call fmod_alloc_yyparse_bond_top(npara)
       call fmod_alloc_yyparse_bond_n(npara)
       do i = 1, nmol
          call parse_specify_key('/topology and parameters/species',i)
          call parse_get(trim(name_TopPar_Spc)//'/nbond',nbond,lerr)
          do j = 1, nbond
             call parse_get(trim(name_TopPar_Spc)//'/bond',ivalue1,lerr1)
             call parse_get(trim(name_TopPar_Spc)//'/bond',ivalue2,lerr2)
             call parse_get(trim(name_TopPar_Spc)//'/bond',qvalue1,lerr3)
             call parse_get(trim(name_TopPar_Spc)//'/bond',qvalue2,lerr4)
             if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4)) then
                !convert unit : [kcal /mol A^2] => [J / m^2]
                value1 = qvalue1*1000.0d0*md_CALORIE/md_AVOGADRO*1.0d20
                !convert unit : [A] => [m]
                value2 = qvalue2 * 1.0d-10
                call  add_to_mol_bond_list(i,ivalue1,ivalue2,value1,value2)
             endif
          enddo
       enddo
       call set_md_charmm_a_bond()

       !READING BOND MORSE PAIR
       call fmod_alloc_yyparse_bond_morse_top(npara)
       call fmod_alloc_yyparse_bond_morse_n(npara)
       do i = 1, nmol
          call parse_specify_key('/topology and parameters/species',i)
          call parse_get(trim(name_TopPar_Spc)//'/nbondmorse',nbond,lerr)
          do j = 1, nbond
             call parse_get(trim(name_TopPar_Spc)//'/bond morse',ivalue1,lerr1)
             call parse_get(trim(name_TopPar_Spc)//'/bond morse',ivalue2,lerr2)
             call parse_get(trim(name_TopPar_Spc)//'/bond morse',qvalue1,lerr3)
             call parse_get(trim(name_TopPar_Spc)//'/bond morse',qvalue2,lerr4)
             call parse_get(trim(name_TopPar_Spc)//'/bond morse',qvalue3,lerr5)
             if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. lerr5)) then
                !convert unit : [kcal /mol] => [J / m^2]
                value1 = qvalue1*1000.0d0*md_CALORIE/md_AVOGADRO
                !convert unit : [A-1] => [m-1]
                value2 = qvalue2 * 1.0d10
                !convert unit : [A] => [m]
                value3 = qvalue3 * 1.0d-10
                call  add_to_mol_bond_morse_list(i,ivalue1,ivalue2,value1,value2, value3)
             endif
          enddo
       enddo
       call set_md_bond_morse()

       !READING BOND MORSE2 PAIR
       call fmod_alloc_yyparse_bond_morse2_top(npara)
       call fmod_alloc_yyparse_bond_morse2_n(npara)
       do i = 1, nmol
         call parse_specify_key('/topology and parameters/species',i)
         call parse_get(trim(name_TopPar_Spc)//'/nbondmorse2',nbond,lerr)
         do j = 1, nbond
           call parse_get(trim(name_TopPar_Spc)//'/bond morse2', ivalue1,lerr1)
           call parse_get(trim(name_TopPar_Spc)//'/bond morse2', ivalue2,lerr2)
           call parse_get(trim(name_TopPar_Spc)//'/bond morse2', qvalue1,lerr3)
           call parse_get(trim(name_TopPar_Spc)//'/bond morse2', qvalue2,lerr4)
           call parse_get(trim(name_TopPar_Spc)//'/bond morse2', qvalue3,lerr5)
           call parse_get(trim(name_TopPar_Spc)//'/bond morse2', qvalue4,lerr6)
           if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. lerr5 .or. lerr6)) then
             !convert unit : [kcal /mol] => [J]
             value1 = qvalue1*1000.0d0*md_CALORIE/md_AVOGADRO
             !convert unit : [A^-1] => [m^-1]
             value2 = qvalue2 * 1d10
             !convert unit : [A^-1] => [m^-1]
             value3 = qvalue3 * 1d20
             !convert unit : [A] => [m]
             value4 = qvalue4 * 1.0d-10
             call  add_to_mol_bond_morse2_list &
                   (i,ivalue1,ivalue2,value1,value2,value3,value4)
           endif
         enddo
       enddo
       call set_md_bond_morse2()

               !READING ANGLE MORSE PAIR
       call fmod_alloc_yyparse_angleA_morse_top(npara)
       call fmod_alloc_yyparse_angleA_morse_n(npara)
       call fmod_alloc_yyparse_angleB_morse_top(npara)
       call fmod_alloc_yyparse_angleB_morse_n(npara)
       do i = 1, nmol
         call parse_specify_key(trim(name_TopPar_Spc),i)
         call parse_get(trim(name_TopPar_Spc)//'/nanglemorse',nag,lerr)
         getname = trim(name_TopPar_Spc)//'/angle morse'
         do j = 1, nag
           call parse_get(trim(getname),ivalue1,lerr1)
           call parse_get(trim(getname),ivalue2,lerr2)
           call parse_get(trim(getname),ivalue3,lerr3)
           call parse_get(trim(getname),qvalue1,lerr4)
           call parse_get(trim(getname),qvalue2,lerr5)
           call parse_get(trim(getname),qvalue3,lerr6)
           call parse_get(trim(getname),qvalue4,lerr7)
           call parse_get(trim(getname),qvalue5,lerr8)
           call parse_get(trim(getname),qvalue6,lerr9)
           call parse_get(trim(getname),qvalue7,lerr10)
           call parse_get(trim(getname),qvalue8,lerr11)
           call parse_get(trim(getname),qvalue9,lerr12)
           call parse_get(trim(getname),qvalue10,lerr13)
           call parse_get(trim(getname),qvalue11,lerr14)
           if (.not. (lerr1  .or. lerr2  .or. lerr3  .or. &
                      lerr4  .or. lerr5  .or. lerr6  .or. &
                      lerr7  .or. lerr8  .or. lerr9  .or. &
                      lerr10 .or. lerr11 .or. lerr12 .or. &
                      lerr13 .or. lerr14)) then
              !convert unit : [kcal/mol.rad] => [J/rad]
              value1 = qvalue1 * 1000.0d0 * md_CALORIE/md_AVOGADRO
              !convert unit : [kcal/mol.rad**2] => [J/rad**2]
              value2 = qvalue2 * 1000.0d0 * md_CALORIE/md_AVOGADRO
              !convert unit : [kcal/mol.rad**3] => [J/rad**3]
              value3 = qvalue3 * 1000.0d0 * md_CALORIE/md_AVOGADRO
              !convert unit : [kcal/mol.rad**4] => [J/rad**4]
              value4 = qvalue4 * 1000.0d0 * md_CALORIE/md_AVOGADRO
              !convert unit : [degree]=>[rad]
              value5 = qvalue5 * md_DEGREE
              !convert unit : [A^-1] => [m^-1]
              value6 = qvalue6 * 1d10
              !convert unit : [A^-1] => [m^-1]
              value7 = qvalue7 * 1d20
              !convert unit : [A] => [m]
              value8 = qvalue8 * 1.0d-10
              !convert unit : [A^-1] => [m^-1]
              value9 = qvalue9 * 1d10
              !convert unit : [A^-1] => [m^-1]
              value10 = qvalue10 * 1d20
              !convert unit : [A] => [m] 
              value11 = qvalue11 * 1.0d-10
              call  add_to_mol_angle_morse_list &
                       (i,ivalue1,ivalue2,ivalue3,value1,value2, &
                       value3, value4, value5, value6, value7, &
                       value8, value9, value10, value11)
            endif
          enddo
        enddo
        call set_md_angleA_morse()
        call set_md_angleB_morse()

    !READING BOND PCFF PAIR
    !call fmod_alloc_yyparse_bond_pcff_top(npara)
    !call fmod_alloc_yyparse_bond_pcff_n(npara)
    !do i = 1, nmol
    !    call parse_specify_key('/topology and parameters/species',i)
    !    call parse_get(trim(name_TopPar_Spc)//'/nbond_pcff',nbond,lerr)
    !   do j = 1, nbond
    !      call parse_get(trim(name_TopPar_Spc)//'/bond_pcff',ivalue1,lerr1)
    !      call parse_get(trim(name_TopPar_Spc)//'/bond_pcff',ivalue2,lerr2)
    !      call parse_get(trim(name_TopPar_Spc)//'/bond_pcff',qvalue1,lerr3)
    !      call parse_get(trim(name_TopPar_Spc)//'/bond_pcff',qvalue2,lerr4)
    !      call parse_get(trim(name_TopPar_Spc)//'/bond_pcff',qvalue3,lerr5)
    !      call parse_get(trim(name_TopPar_Spc)//'/bond_pcff',qvalue4,lerr6)
    !      if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 &
    !            .or. lerr5 .or. lerr6)) then
    !         !convert unit : [kcal /mol A^2] => [J / m^2]
    !         value1 = qvalue1*1000.0d0*md_CALORIE/md_AVOGADRO*1.0d20
    !         !convert unit : [kcal /mol A^3] => [J / m^3]
    !         value2 = qvalue2*1000.0d0*md_CALORIE/md_AVOGADRO*1.0d30
    !         !convert unit : [kcal /mol A^4] => [J / m^4]
    !         value3 = qvalue3*1000.0d0*md_CALORIE/md_AVOGADRO*1.0d40
    !         !convert unit : [A] => [m]
    !         value4 = qvalue4 * 1.0d-10
    !         call add_to_mol_bond_pcff_list &
    !              (i,ivalue1,ivalue2,value1,value2,value3,value4)
    !      endif
    !   enddo
    !enddo
    !call set_md_pcff_a_bond()


       !READING ANGLE PAIR
       call fmod_alloc_yyparse_angleA_top(npara)
       call fmod_alloc_yyparse_angleA_n(npara)
       call fmod_alloc_yyparse_angleB_top(npara)
       call fmod_alloc_yyparse_angleB_n(npara)
       do i = 1, nmol
          call parse_specify_key(trim(name_TopPar_Spc),i)
          call parse_get(trim(name_TopPar_Spc)//'/nangle',nag,lerr)
          getname = trim(name_TopPar_Spc)//'/angle'
          do j = 1, nag
             call parse_get(trim(getname),ivalue1,lerr1)
             call parse_get(trim(getname),ivalue2,lerr2)
             call parse_get(trim(getname),ivalue3,lerr3)
             call parse_get(trim(getname),qvalue1,lerr4)
             call parse_get(trim(getname),qvalue2,lerr5)
             if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. lerr5)) then
                !convert unit : [kcal/mol.rad**2] => [J/rad**2]
                value1 = qvalue1 * 1000.0d0 * md_CALORIE/md_AVOGADRO
                !convert unit : [degree]=>[rad]
                value2 = qvalue2 * md_DEGREE
                call  add_to_mol_angle_list(i,ivalue1,ivalue2,ivalue3,value1,value2)
             endif
          enddo
       enddo
       call set_md_charmm_angleA
       call set_md_charmm_angleB

       !READING ANGLE PCFF PAIR
       !call fmod_alloc_yyparse_angle_pcffA_top(npara)
       !call fmod_alloc_yyparse_angle_pcffA_n(npara)
       !call fmod_alloc_yyparse_angle_pcffB_top(npara)
       !call fmod_alloc_yyparse_angle_pcffB_n(npara)
       !do i = 1, nmol
       !   call parse_specify_key(trim(name_TopPar_Spc),i)
       !   call parse_get(trim(name_TopPar_Spc)//'/nangle_pcff',nag,lerr)
       !   getname = trim(name_TopPar_Spc)//'/angle_pcff'
       !   do j = 1, nag
       !      call parse_get(trim(getname),ivalue1,lerr1)
       !      call parse_get(trim(getname),ivalue2,lerr2)
       !      call parse_get(trim(getname),ivalue3,lerr3)
       !      call parse_get(trim(getname),qvalue1,lerr4)
       !      call parse_get(trim(getname),qvalue2,lerr5)
       !      call parse_get(trim(getname),qvalue3,lerr6)
       !      call parse_get(trim(getname),qvalue4,lerr7)
       !      if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. lerr5 &
       !            .or. lerr6 .or. lerr7 )) then
       !         !convert unit : [kcal/mol/rad**2] => [J/rad**2]
       !         value1 = qvalue1 * 1000.0d0 * md_CALORIE/md_AVOGADRO
       !         !convert unit : [kcal/mol/rad**3] => [J/rad**3]
       !         value2 = qvalue2 * 1000.0d0 * md_CALORIE/md_AVOGADRO
       !         !convert unit : [kcal/mol/rad**4] => [J/rad**4]
       !         value3 = qvalue3 * 1000.0d0 * md_CALORIE/md_AVOGADRO
       !         !convert unit : [degree]=>[rad]
       !         value4 = qvalue4 * md_DEGREE
       !         call  add_to_mol_angle_pcff_list&
       !               (i,ivalue1,ivalue2,ivalue3,value1,value2,value3,value4)
       !      endif
       !   enddo
       !enddo
       !call set_md_pcff_angleA
       !call set_md_pcff_angleB

       !READING PCFF BOND BOND PAIR
       !call fmod_alloc_yyparse_bond_bondA_top(npara)
       !call fmod_alloc_yyparse_bond_bondA_n(npara)
       !call fmod_alloc_yyparse_bond_bondB_top(npara)
       !call fmod_alloc_yyparse_bond_bondB_n(npara)
       !do i = 1, nmol
       !   call parse_specify_key(trim(name_TopPar_Spc),i)
       !   call parse_get(trim(name_TopPar_Spc)//'/nbond_bond',nag,lerr)
       !   getname = trim(name_TopPar_Spc)//'/bond_bond'
       !   do j = 1, nag
       !      call parse_get(trim(getname),ivalue1,lerr1)
       !      call parse_get(trim(getname),ivalue2,lerr2)
       !      call parse_get(trim(getname),ivalue3,lerr3)
       !      call parse_get(trim(getname),qvalue1,lerr4)
       !      call parse_get(trim(getname),qvalue2,lerr5)
       !      call parse_get(trim(getname),qvalue3,lerr6)
       !      if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. lerr5 &
       !            .or. lerr6)) then
       !         !convert unit : [kcal/mol/A**2] => [J/m**2]
       !         value1 = qvalue1*1000.0d0*md_CALORIE/md_AVOGADRO*1.0d20
       !         !convert unit : [A] => [m]
       !         value2 = qvalue2 * 1.0d-10
       !         !convert unit : [A] => [m]
       !         value3 = qvalue3 * 1.0d-10
       !         call  add_to_mol_bond_bond_list&
       !               (i,ivalue1,ivalue2,ivalue3,value1,value2,value3)
       !      endif
       !   enddo
       !enddo
       !call set_md_pcff_bbA
       !call set_md_pcff_bbB

       !READING PCFF BOND ANGLE PAIR
       !call fmod_alloc_yyparse_bond_angleA_top(npara)
       !call fmod_alloc_yyparse_bond_angleA_n(npara)
       !call fmod_alloc_yyparse_bond_angleB_top(npara)
       !call fmod_alloc_yyparse_bond_angleB_n(npara)
       !do i = 1, nmol
       !   call parse_specify_key(trim(name_TopPar_Spc),i)
       !   call parse_get(trim(name_TopPar_Spc)//'/nbond_angle',nag,lerr)
       !   getname = trim(name_TopPar_Spc)//'/bond_angle'
       !   do j = 1, nag
       !      call parse_get(trim(getname),ivalue1,lerr1)
       !      call parse_get(trim(getname),ivalue2,lerr2)
       !      call parse_get(trim(getname),ivalue3,lerr3)
       !      call parse_get(trim(getname),qvalue1,lerr4)
       !      call parse_get(trim(getname),qvalue2,lerr5)
       !      call parse_get(trim(getname),qvalue3,lerr6)
       !      call parse_get(trim(getname),qvalue4,lerr7)
       !      call parse_get(trim(getname),qvalue5,lerr8)
       !      if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. lerr5 &
       !          .or. lerr6   .or. lerr7 .or. lerr8)) then
       !         !convert unit : [kcal/mol/A/rad] => [J/m/rad]
       !         value1 = qvalue1*1000.0d0*md_CALORIE/md_AVOGADRO*1.0d10
       !         !convert unit : [kcal/mol/A/rad] => [J/m/rad]
       !         value2 = qvalue2*1000.0d0*md_CALORIE/md_AVOGADRO*1.0d10
       !         !convert unit : [A] => [m]
       !         value3 = qvalue3 * 1.0d-10
       !         !convert unit : [A] => [m]
       !         value4 = qvalue4 * 1.0d-10
       !         !convert unit : [degree]=>[rad]
       !         value5 = qvalue5 * md_DEGREE
       !         call  add_to_mol_bond_angle_list&
       !              (i,ivalue1,ivalue2,ivalue3,value1,value2,value3,value4,value5)
       !      endif
       !   enddo
       !enddo
       !call set_md_pcff_baA
       !call set_md_pcff_baB

       !READING PCFF ANGLE ANGLE PAIR
       !call fmod_alloc_yyparse_angle_angleA_top(npara)
       !call fmod_alloc_yyparse_angle_angleA_n(npara)
       !call fmod_alloc_yyparse_angle_angleB_top(npara)
       !call fmod_alloc_yyparse_angle_angleB_n(npara)
       !do i = 1, nmol
       !   call parse_specify_key(trim(name_TopPar_Spc),i)
       !   call parse_get(trim(name_TopPar_Spc)//'/nangle_angle',nag,lerr)
       !   getname = trim(name_TopPar_Spc)//'/angle_angle'
       !   do j = 1, nag
       !      call parse_get(trim(getname),ivalue1,lerr1)
       !      call parse_get(trim(getname),ivalue2,lerr2)
       !      call parse_get(trim(getname),ivalue3,lerr3)
       !      call parse_get(trim(getname),ivalue4,lerr4)
       !      call parse_get(trim(getname),qvalue1,lerr5)
       !      call parse_get(trim(getname),qvalue2,lerr6)
       !      call parse_get(trim(getname),qvalue3,lerr7)
       !      call parse_get(trim(getname),qvalue4,lerr8)
       !      if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. lerr5 &
       !          .or. lerr6   .or. lerr7 .or. lerr8)) then
       !         !convert unit : [kcal/mol/rad**2] => [J/rad**2]
       !         value1 = qvalue1 * 1000.0d0 * md_CALORIE/md_AVOGADRO
       !         !convert unit : [degree]=>[rad]
       !         value2 = qvalue2 * md_DEGREE
       !         !convert unit : [degree]=>[rad]
       !         value3 = qvalue3 * md_DEGREE
       !         !convert unit : [degree]=>[rad]
       !         value4 = qvalue4 * md_DEGREE
       !         call  add_to_mol_angle_angle_list&
       !          (i,ivalue1,ivalue2,ivalue3,ivalue4,value1,value2,value3,value4)
       !      endif
       !   enddo
       !enddo
       !call set_md_pcff_aaA
       !call set_md_pcff_aaB

       !READING UB PAIR
       call fmod_alloc_yyparse_ub_top(npara)
       call fmod_alloc_yyparse_ub_n(npara)
       do i = 1, nmol
          call parse_specify_key('/topology and parameters/species',i)
          call parse_get(trim(name_TopPar_Spc)//'/nub',nub,lerr)
          do  j = 1, nub
             call parse_get(trim(name_TopPar_Spc)//'/ub',ivalue1,lerr1)
             call parse_get(trim(name_TopPar_Spc)//'/ub',ivalue2,lerr2)
             call parse_get(trim(name_TopPar_Spc)//'/ub',ivalue3,lerr3)
             call parse_get(trim(name_TopPar_Spc)//'/ub',qvalue1,lerr4)
             call parse_get(trim(name_TopPar_Spc)//'/ub',qvalue2,lerr5)
             if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. lerr5)) then
                !convert unit : pkcal/mol.A**2] => [J/m**2]
                value1 = qvalue1*1000.0d0*md_CALORIE/md_AVOGADRO*1.0d20
                !convert unit : [A] => [m]
                value2 = qvalue2 * 1.0d-10
                call  add_to_mol_ub_list(i,ivalue1,ivalue2,ivalue3,value1,value2)
             endif
          enddo
       enddo
       call set_md_charmm_a_ub()

       !READING DIHEDRAL PAIRS
       call fmod_alloc_yyparse_dihedA_top(npara)
       call fmod_alloc_yyparse_dihedA_n(npara)
       call fmod_alloc_yyparse_dihedB_top(npara)
       call fmod_alloc_yyparse_dihedB_n(npara)
#ifdef DIHEDRAL_TABLE
       call init_dihedral_table()
#endif
       do i = 1, nmol
          call parse_specify_key('/topology and parameters/species',i)
          call parse_get(trim(name_TopPar_Spc)//'/ndihedral',ndi,lerr)
          getname = trim(name_TopPar_Spc)//'/dihedral'
          do j = 1, ndi
             call parse_get(trim(getname),ivalue1,lerr1)
             call parse_get(trim(getname),ivalue2,lerr2)
             call parse_get(trim(getname),ivalue3,lerr3)
             call parse_get(trim(getname),ivalue4,lerr4)
             call parse_get(trim(getname),qvalue1,lerr5)
             call parse_get(trim(getname),qvalue2,lerr6)
             call parse_get(trim(getname),qvalue3,lerr7)
             if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. lerr5 .or. lerr6 .or. lerr7)) then

                if( md_condition__force_field == CHARMM .or. md_condition__force_field == KREMER .or.  &
                    !md_condition__force_field == PCFF   .or. md_condition__force_field == GAmbFF) then
                                                             md_condition__force_field == GAmbFF) then

                   !convert unit : [kcal/mol] => [J]
                   value1 = qvalue1 * 1000.0d0 * md_CALORIE/md_AVOGADRO
                   value2 = qvalue2
                   !convert unit : [degree] => [rad]
                   value3 = qvalue3 * md_DEGREE
                   call add_to_mol_dihedral_list(i, ivalue1, ivalue2, ivalue3, ivalue4, value1, value2, value3)

                elseif(md_condition__force_field == OPLSAA .or. md_condition__force_field == AMBER  )then

                   !convert unit : [kcal/mol] => [J]
                   qvalue1 = qvalue1 * 1000.0d0 * md_CALORIE/md_AVOGADRO
                   value1 = qvalue1 * 0.5d0
                   value2 = qvalue2
                   !convert unit : [degree] => [rad]
                   value3 = qvalue3 * md_DEGREE
                   value3 =         - value3
                   call add_to_mol_dihedral_list(i, ivalue1, ivalue2, ivalue3, ivalue4, value1, value2, value3)
                else
                   write(0,*) 'ERROR: dihedral reading'
                   call modylas_abort()
                endif
             endif
          enddo
       enddo
       call set_md_charmm_dihedralA()
       call set_md_charmm_dihedralB()
#ifdef DIHEDRAL_TABLE
       call finalize_dihedral_table()
#endif
       if(    md_condition__force_field==CHARMM) then
          !READING CMAP PAIRS
          call fmod_alloc_yyparse_CMAP_topA(npara)
          call fmod_alloc_yyparse_CMAP_nA(npara)
          call fmod_alloc_yyparse_CMAP_topB(npara)
          call fmod_alloc_yyparse_CMAP_nB(npara)
          call fmod_alloc_yyparse_CMAP_topC(npara)
          call fmod_alloc_yyparse_CMAP_nC(npara)
          call fmod_alloc_yyparse_CMAP_topD(npara)
          call fmod_alloc_yyparse_CMAP_nD(npara)
          call fmod_alloc_yyparse_CMAP_topE(npara)
          call fmod_alloc_yyparse_CMAP_nE(npara)
          do i =1, nmol
             call parse_specify_key('/topology and parameters/species',i)
             call parse_get(trim(name_TopPar_Spc)//'/ncmap',ncm,lerr)
             do j = 1, ncm
                call parse_get(trim(name_TopPar_Spc)//'/cmap',ivalue1,lerr1)
                call parse_get(trim(name_TopPar_Spc)//'/cmap',ivalue2,lerr2)
                call parse_get(trim(name_TopPar_Spc)//'/cmap',ivalue3,lerr3)
                call parse_get(trim(name_TopPar_Spc)//'/cmap',ivalue4,lerr4)
                call parse_get(trim(name_TopPar_Spc)//'/cmap',ivalue5,lerr5)
                call parse_get(trim(name_TopPar_Spc)//'/cmap',ivalue6,lerr6)
                if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. lerr5 .or. lerr6 )) then
                   call add_to_mol_cmap_list(i, ivalue1, ivalue2, ivalue3, ivalue4, ivalue5, ivalue6)
                endif
             enddo
          enddo
          call set_md_charmm_CMAPA()
          call set_md_charmm_CMAPB()
          call set_md_charmm_CMAPC()
          call set_md_charmm_CMAPD()
          call set_md_charmm_CMAPE()
       endif

       !READING ITORSION PAIRS
       call fmod_alloc_yyparse_itorsA_top(npara)
       call fmod_alloc_yyparse_itorsA_n(npara)
       call fmod_alloc_yyparse_itorsB_top(npara)
       call fmod_alloc_yyparse_itorsB_n(npara)
       do i =1, nmol
          call parse_specify_key('/topology and parameters/species',i)
          call parse_get(trim(name_TopPar_Spc)//'/nitorsion',nit,lerr)
          getname = trim(name_TopPar_Spc)//'/itorsion'
          do j = 1, nit
#if defined(GROEXT_OPLS) && defined(GAFF)
#error Please do not define GROEXT_OPLS and GAFF at the same time in cmake
#error This is because the improper parameter will not not read properly by doing so
#endif
#ifdef GROEXT_OPLS
             call parse_get(trim(getname), ivalue1,lerr1)
             call parse_get(trim(getname), ivalue2,lerr2)
             call parse_get(trim(getname), ivalue3,lerr3)
             call parse_get(trim(getname), ivalue4,lerr4)
             call parse_get(trim(getname), qvalue1,lerr5)
             call parse_get(trim(getname), qvalue2,lerr6)
             call parse_get(trim(getname), qvalue3,lerr7)
             !convert unit : [kcal/mol.rad**2] => [J/rad**2]
             qvalue1 = qvalue1 * 1000.0d0 * md_CALORIE / md_AVOGADRO
             value1 = qvalue1 * 0.5d0
             value2 = qvalue2
             !convert unit:[degree]=>[rad]
             value3 = qvalue3 * md_DEGREE
             if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. lerr5 .or. lerr6 .or. lerr7)) then
                call add_to_mol_itorsion_list(i, ivalue1, ivalue2, ivalue3, ivalue4, value1, value2, value3)
             endif
#elif defined(GAFF)
             call parse_get(trim(getname), ivalue1,lerr1)
             call parse_get(trim(getname), ivalue2,lerr2)
             call parse_get(trim(getname), ivalue3,lerr3)
             call parse_get(trim(getname), ivalue4,lerr4)
             call parse_get(trim(getname), qvalue1,lerr5)
             call parse_get(trim(getname), ivalue5,lerr6)
             call parse_get(trim(getname), ivalue6,lerr7)
             !E=K[1+d*cos(n*phi)]
             !convert unit : [kcal/mol] => [J]
             value1 = qvalue1 * 1000.0d0 * md_CALORIE / md_AVOGADRO ! K
             value2 = ivalue5   ! d
             value3 = ivalue6   ! n
             if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. lerr5 .or. lerr6 .or. lerr7)) then
               call add_to_mol_itorsion_list(i, ivalue1, ivalue2, ivalue3, ivalue4, value1, value2, value3)
             endif
#else
             call parse_get(trim(getname),ivalue1,lerr1)
             call parse_get(trim(getname),ivalue2,lerr2)
             call parse_get(trim(getname),ivalue3,lerr3)
             call parse_get(trim(getname),ivalue4,lerr4)
             call parse_get(trim(getname),qvalue1,lerr5)
             call parse_get(trim(getname),qvalue2,lerr6)
             !convert unit : [kcal/mol.rad**2] => [J/rad**2]
             value1 = qvalue1 * 1000.0d0 * md_CALORIE / md_AVOGADRO
             !convert unit:[degree]=>[rad]
             value2 = qvalue2 * md_DEGREE
             if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. lerr5 .or. lerr6)) then
                call add_to_mol_itorsion_list(i, ivalue1, ivalue2, ivalue3, ivalue4, value1, value2)
             endif
#endif
          enddo
       enddo
       call set_md_charmm_itorsionA()
       call set_md_charmm_itorsionB()

       !--- read position constrain --
       call parse_get('/position constrain/type', chara, lerr)
       if (.not. lerr) then
          if(trim(chara) =='harmonic') then
             call fmod_p_cons__type(1)
             getname = '/position constrain/force_constant'
             call parse_get(trim(getname), value1, lerr)
             if (.not. lerr)  then
                !convert unit : [kcal /mol A^2] => [J / m^2]
                value1 = value1*1000.0d0*md_CALORIE/md_AVOGADRO*1.0d20
                call fmod_p_cons__fc(value1)
             endif
          else if(trim(chara) =='huge_mass') then
             call fmod_p_cons__type(2)
          else if(trim(chara) =='fix') then
             call fmod_p_cons__type(3)
          endif
       endif
       getname = '/position constrain/atom/atom_type'
       call parse_get(trim(getname), chara, lerr)
       if (.not. lerr) then
          if(trim(chara) =='all') then
             call fmod_p_cons_atom__atom_type(0)
          else if(trim(chara) =='heavy_atoms') then
             call fmod_p_cons_atom__atom_type(1)
          endif
       endif
       getname = '/position constrain/atom/molecules'
       call parse_get(trim(getname), chara, lerr)
       if (.not. lerr) then
          if(trim(chara) =='all') then
             call fmod_p_cons_atom__molecules(0)
          else if(trim(chara) =='specified') then
             call fmod_p_cons_atom__molecules(1)
             getname = '/position constrain/atom/nmolecule'
             call parse_get(trim(getname), ivalue1, lerr)
             if (.not. lerr)  then
                call fmod_p_cons_atom__nmolecule(ivalue1)
                call fmod_alloc_mol_p_constrain(ivalue1)
             endif
             k = ivalue1
             do  i= 1,k
                getname = '/position constrain/atom/molecule'
                call parse_get(trim(getname),ivalue1,lerr1)
                call parse_get(trim(getname),ivalue2,lerr2)
                if (.not. (lerr1 .or. lerr2 )) then
                   !if(ivalue1==0) then
                   !   im = ivalue2+1
                   !else
                   im = 0
                   do j= 1,ivalue1
                      im = im + nmolmol(j)
                   enddo
                   im = im + ivalue2 + 1
                   !endif
                   call add_to_mol_p_constrain(i,im)
                endif
             enddo
          endif
       endif
       getname = '/position constrain/atom/natom_allow_unconstrain'
       call parse_get(trim(getname), ivalue1, lerr)
       if (.not. lerr)  then
          call fmod_p_cons_atom__natm_alw_ucn(ivalue1)
          call fmod_alloc_natom_allow_uncons(ivalue1)
          k = ivalue1
          do  i= 1,k
             getname = '/position constrain/atom/allow atom unconstrain'
             call parse_get(trim(getname),ivalue1, lerr)
             if(.not.lerr)call add_to_atom_allow_uncons(i,ivalue1+1)
          enddo
       endif
       getname = '/position constrain/position/type'
       call parse_get(trim(getname), chara, lerr)
       if (.not. lerr) then
          if(trim(chara) =='initial') then
             call fmod_p_cons_position__type(0)
          else if(trim(chara) =='specified') then
             call fmod_p_cons_position__type(1)
             getname = '/position constrain/position/coordinate'
             call fmod_alloc_ref_crd0_p_cons
             do while(.true.)
                call parse_get(trim(getname),ivalue1,lerr)
                if (lerr) exit
                call parse_get(trim(getname),qvalue1,lerr1)
                call parse_get(trim(getname),qvalue2,lerr2)
                call parse_get(trim(getname),qvalue3,lerr3)
                if ( lerr1 .or. lerr2 .or. lerr3 ) then
                   write(0,*) 'ERROR:/position constrain/position/coordinate'
                   call modylas_abort()
                endif
                value1 = qvalue1 * 1.0d-10
                value2 = qvalue2 * 1.0d-10
                value3 = qvalue3 * 1.0d-10
                call fmod_p_cons_position__ref_crd0(ivalue1+1,value1,value2,value3 )
             enddo
          endif
       endif

       call parse_close(lerr)
    endif

    if(myrank==0) write(*,*) 'Reading .mdff file ended successfully!'


    !
    !     mddef 1.0.0
    !
200 continue
    if (read_mddef) then
       call parse_open(trim(session_name)//'.mddef', lerr)
       if(myrank==0) write(*,*) 'Reading .mddef file started.'

       call parse_get('/input/version', chara, lerr)
#ifdef GROEXT
       call parse_groext()
#endif
       ! output
       call parse_get('/output/ascii', chara, lerr)
       if(trim(chara) == 'yes') then
          call fmod_ascii_output(1)
       else if(trim(chara) == 'no') then
          call fmod_ascii_output(0)
       endif
       !backup files starting with #
       call parse_get('/output/backup', chara, lerr)
       IF(trim(chara) == 'yes') THEN
          call fmod_backup_output(1)
       else if(trim(chara) == 'no') then
          call fmod_backup_output(0)
       ENDIF
       call parse_get('/output/restart/start', ivalue1, lerr)
       if (.not. lerr)  call fmod_restart_start(ivalue1)
       call parse_get('/output/restart/interval', ivalue1, lerr)
       if (.not. lerr)  call fmod_restart_interval(ivalue1)
       call parse_get('/output/monitor/start', ivalue1, lerr)
       if (.not. lerr)  call fmod_mntr_start(ivalue1)
       call parse_get('/output/monitor/interval', ivalue1, lerr)
       if (.not. lerr)  call fmod_mntr_interval(ivalue1)
       call parse_get('/output/force/start', ivalue1, lerr)
       if (.not. lerr)  call fmod_force_start(ivalue1)
       call parse_get('/output/force/interval', ivalue1, lerr)
       if (.not. lerr)  call fmod_force_interval(ivalue1)
       call parse_get('/output/system_dipole',chara,lerr)
       if (.not. lerr .and. trim(chara) == 'yes') then
         call set_up_system_dipole(1)
       endif

       !mdtrj
       call parse_get('/output/mdtrj', chara, lerr)
       IF(trim(chara) == 'yes') THEN
          call fmod_mdtrj_output(1)
       ENDIF
       call parse_get('/output/trajectory/start', ivalue1, lerr)
       if (.not. lerr)  call fmod_trj_start(ivalue1)
       call parse_get('/output/trajectory/interval', ivalue1, lerr)
       if (.not. lerr)  call fmod_trj_interval(ivalue1)
       !dcd
       call parse_get('/output/dcd', chara, lerr)
       IF(trim(chara) == 'yes') THEN
          call fmod_dcd_output(1)
       ENDIF
       call parse_get('/output/trjdcd/start', ivalue1, lerr)
       if (.not. lerr)  call fmod_dcd_start(ivalue1)
       call parse_get('/output/trjdcd/interval', ivalue1, lerr)
       if (.not. lerr)  call fmod_dcd_interval(ivalue1)

#ifdef XTC
          !xtc
          call parse_get('/output/xtc', chara, lerr)
          IF(trim(chara) == 'yes') THEN
             call fmod_xtc_output(1)
          ENDIF
          call parse_get('/output/trjxtc/start', ivalue1, lerr)
          if (.not. lerr)  call fmod_xtc_start(ivalue1)
          call parse_get('/output/trjxtc/interval', ivalue1, lerr)
          if (.not. lerr)  call fmod_xtc_interval(ivalue1)
#endif
       !dcd
!#ifdef XTC
!          !xtc
!          call parse_get('/output/xtc', chara, lerr)
!          IF(trim(chara) == 'yes') THEN
!             call fmod_xtc_output(1)
!          ENDIF
!          call parse_get('/output/trjxtc/start', ivalue1, lerr)
!          if (.not. lerr)  call fmod_xtc_start(ivalue1)
!          call parse_get('/output/trjxtc/interval', ivalue1, lerr)
!          if (.not. lerr)  call fmod_xtc_interval(ivalue1)
!#endif

! TZYBEGIN
       call parse_get('/integrator/dt', value1, lerr)
      !  if (.not. lerr)  call fmod_md_condition__dt(value1)
       if (.not. lerr) then
           call fmod_md_condition__dt(value1)
           if( myrank == 0 ) write(*,*) "dt is deprecated. Use dt_long instead"
       endif
       call parse_get('/integrator/dt_long', value1, lerr)
       if (.not. lerr)  call fmod_md_condition__dt(value1)
! TZTEND
       call parse_get('/integrator/initial_step', ivalue1, lerr)
       if(.not. lerr) then
          call fmod_mdstep(ivalue1) !generic step
       else
          call fmod_mdstep(0) !generic step
       endif

       call parse_get('/integrator/initial_step', value1, lerr)
       
       call parse_get('/scaling14/special_divide_lj', value1, lerr)
       if(.not. lerr) call fmod_md_oplsaa_lj_sp_divide_mddef(value1)

       call parse_get('/scaling14/special_divide_coulomb', value1, lerr)
       if(.not. lerr) call fmod_md_oplsaa_coulomb_sp_divide_mddef(value1)

       ! ensemble
       call parse_get('/ensemble/velocity_scaling',chara,lerr)
       if (.not. lerr .and. trim(chara) == 'yes') then
          velocity_scaling_flag = .true.

          call parse_get('/ensemble/velocity_scaling_region_wise',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             velocity_scaling_region_wise_flag = .true.
          endif
       endif
       !
       call parse_get('/ensemble/volume_scaling',chara,lerr)
       if (.not. lerr .and. trim(chara) == 'yes') then
          volume_scaling_flag = .true.
       endif
       call parse_get('/ensemble/volume', value1, lerr)
       if (.not. lerr)  call fmod_sysvolume(value1)
       !
       call parse_get('/ensemble/cellgeometry_scaling',chara,lerr)
       if (.not. lerr .and. trim(chara) == 'yes') then
          cellgeometry_scaling_flag = .true.
       endif
       call parse_get('/ensemble/cellx', value1, lerr)
       if (.not. lerr)  call fmod_syscellx(value1)
       call parse_get('/ensemble/celly', value1, lerr)
       if (.not. lerr)  call fmod_syscelly(value1)
       call parse_get('/ensemble/cellz', value1, lerr)
       if (.not. lerr)  call fmod_syscellz(value1)
       call parse_get('/ensemble/alpha', value1, lerr)
       if (.not. lerr)  call fmod_sysalpha(value1)
       call parse_get('/ensemble/beta', value1, lerr)
       if (.not. lerr)  call fmod_sysbeta(value1)
       call parse_get('/ensemble/gamma', value1, lerr)
       if (.not. lerr)  call fmod_sysgamma(value1)
       ! 
       if (volume_scaling_flag)    call fmod_volume_scaling(1)
       if (cellgeometry_scaling_flag)    call fmod_cellgeometry_scaling(1)
       !

       call parse_get('/ensemble/deformx', value1, lerr)
       if(.not. lerr) call fmod_md_condition__deformx(value1)

       call parse_get('/ensemble/deformy', value1, lerr)
       if(.not. lerr) call fmod_md_condition__deformy(value1)

       call parse_get('/ensemble/deformz', value1, lerr)
       if(.not. lerr) call fmod_md_condition__deformz(value1)

       call parse_get('/ensemble/bond_breaking_distance', value1, lerr)
       if(.not. lerr) call fmod_set_allow_bond_breaking(value1)
       
       call parse_get('/ensemble/ensemble', chara, lerr)
       if (trim(chara) == 'nve') then
          call fmod_md_condition__ensemble(NVE)
          if (velocity_scaling_flag)  call fmod_velocity_scaling(velocity_scaling_flag, velocity_scaling_region_wise_flag)
       elseif (trim(chara) == 'nvt') then
          call fmod_md_condition__ensemble(NVT)
          if (velocity_scaling_flag) then
             call parse_abort('Velocity scaling should be used with NVE.')
          endif
          getname='/ensemble/thermostat/tau_Q'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauqtherm(value1)
          call parse_get('/ensemble/thermostat/initialize',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             initialize_thermostat_flag = .true.
          endif
!
       elseif (trim(chara) == 'nvlxlylzt') then
          call fmod_md_condition__ensemble(NVLXLYLZT)
          getname='/ensemble/thermostat/tau_Q'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauqtherm(value1)
          call parse_get('/ensemble/thermostat/initialize',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             initialize_thermostat_flag = .true.
          endif
          
       elseif (trim(chara) == 'npt_a') then
          call fmod_md_condition__ensemble(NPT_A)
          if (velocity_scaling_flag) then
             call parse_abort('Velocity scaling cannot be used in NPT.')
          endif
          getname='/ensemble/thermostat/tau_Q'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauqtherm(value1)
          call parse_get('/ensemble/thermostat/initialize',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             initialize_thermostat_flag = .true.
          endif
          getname='/ensemble/barostat/tau_Q'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauqbaro(value1)
          getname='/ensemble/barostat/tau_W'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauwpres(value1)
          call parse_get('/ensemble/barostat/initialize',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             initialize_barostat_flag = .true.
          endif

       elseif (trim(chara) == 'npt_z') then
          call fmod_md_condition__ensemble(NPT_Z)
          if (velocity_scaling_flag) then
            call parse_abort('Velocity scaling cannot be used in NPT.')
          endif
          getname='/ensemble/thermostat/tau_Q'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauqtherm(value1)
          call parse_get('/ensemble/thermostat/initialize',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             initialize_thermostat_flag = .true.
          endif
          getname='/ensemble/barostat/tau_Q'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauqbaro(value1)
          getname='/ensemble/barostat/tau_W'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauwpres(value1)
          call parse_get('/ensemble/barostat/initialize',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             initialize_barostat_flag = .true.
          endif
          
       elseif (trim(chara) == 'nptlz') then
          call fmod_md_condition__ensemble(NPTLZ)
          if (velocity_scaling_flag) then
            call parse_abort('Velocity scaling cannot be used in NPT.')
          endif
          getname='/ensemble/thermostat/tau_Q'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauqtherm(value1)
          call parse_get('/ensemble/thermostat/initialize',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             initialize_thermostat_flag = .true.
          endif
          getname='/ensemble/barostat/tau_Q'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauqbaro(value1)
          getname='/ensemble/barostat/tau_W'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauwpres(value1)
          call parse_get('/ensemble/barostat/initialize',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             initialize_barostat_flag = .true.
          endif

       elseif (trim(chara) == 'nlxlypzt') then
          call fmod_md_condition__ensemble(NLXLYPZT)
          if (velocity_scaling_flag) then
            call parse_abort('Velocity scaling cannot be used in NPT.')
          endif
          getname='/ensemble/thermostat/tau_Q'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauqtherm(value1)
          call parse_get('/ensemble/thermostat/initialize',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             initialize_thermostat_flag = .true.
          endif
          getname='/ensemble/barostat/tau_Q'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauqbaro(value1)
          getname='/ensemble/barostat/tau_W'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauwpres(value1)
          call parse_get('/ensemble/barostat/initialize',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             initialize_barostat_flag = .true.
          endif

       elseif (trim(chara) == 'npxpypzt') then
          call fmod_md_condition__ensemble(NPXPYPZT)
          if (velocity_scaling_flag) then
            call parse_abort('Velocity scaling cannot be used in NPT.')
          endif
          getname='/ensemble/thermostat/tau_Q'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauqtherm(value1)
          call parse_get('/ensemble/thermostat/initialize',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             initialize_thermostat_flag = .true.
          endif
          getname='/ensemble/barostat/tau_Q'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauqbaro(value1)
          getname='/ensemble/barostat/tau_W'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauwpres(value1)
          call parse_get('/ensemble/barostat/initialize',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             initialize_barostat_flag = .true.
          endif

       elseif (trim(chara) == 'nlxlylzt') then
          call fmod_md_condition__ensemble(NLXLYLZT)
          if (velocity_scaling_flag) then
            call parse_abort('Velocity scaling cannot be used in NPT.')
          endif
          getname='/ensemble/thermostat/tau_Q'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauqtherm(value1)
          call parse_get('/ensemble/thermostat/initialize',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             initialize_thermostat_flag = .true.
          endif
          getname='/ensemble/barostat/tau_Q'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauqbaro(value1)
          getname='/ensemble/barostat/tau_W'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauwpres(value1)
          call parse_get('/ensemble/barostat/initialize',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             initialize_barostat_flag = .true.
          endif

       elseif (trim(chara) == 'nptlzxy') then
          call fmod_md_condition__ensemble(NPTLZxy)
          if (velocity_scaling_flag) then
            call parse_abort('Velocity scaling cannot be used in NPT.')
          endif
          getname='/ensemble/thermostat/tau_Q'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauqtherm(value1)
          call parse_get('/ensemble/thermostat/initialize',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             initialize_thermostat_flag = .true.
          endif
          getname='/ensemble/barostat/tau_Q'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauqbaro(value1)
          getname='/ensemble/barostat/tau_W'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauwpres(value1)
          call parse_get('/ensemble/barostat/initialize',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             initialize_barostat_flag = .true.
          endif
          
       elseif (trim(chara) == 'npt_pr') then
          
          call fmod_md_condition__ensemble(NPT_PR)
          if (velocity_scaling_flag) then
             call parse_abort('Velocity scaling cannot be used in NPT.')
          endif
          getname='/ensemble/thermostat/tau_Q'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauqtherm(value1)
          call parse_get('/ensemble/thermostat/initialize',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             initialize_thermostat_flag = .true.
          endif
          getname='/ensemble/barostat/tau_Q'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauqbaro(value1)
          getname='/ensemble/barostat/tau_W'
          call parse_get(trim(getname), value1,lerr)
          if (.not. lerr)  call fmod_tauwpres(value1)
          call parse_get('/ensemble/barostat/initialize',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             initialize_barostat_flag = .true.
          endif
          ! For "npt_pr", read pressure before ensemble tag.
          call parse_get('/ensemble/pressure', value1, lerr)
          if (.not. lerr)  call fmod_syspres(value1)
          ! NtT ensemble (constant external stress)
          ! Get reference cell (if it exists in mdff)
          call parse_get('/ensemble/barostat/lx', value1, lerr1)
          call parse_get('/ensemble/barostat/ly', value2, lerr2)
          call parse_get('/ensemble/barostat/lz', value3, lerr3)
          call parse_get('/ensemble/barostat/alpha', value4, lerr4)
          call parse_get('/ensemble/barostat/beta', value5, lerr5)
          call parse_get('/ensemble/barostat/gamma', value6, lerr6)
          if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. lerr5 .or. lerr6) ) then
             call fmod_referece_box(value1,value2,value3,value4,value5,value6)
          endif
          ! Get external stress (if it exists in mdff)
          call parse_get('/ensemble/barostat/stress1', value1, lerr1)	!!xx
          call parse_get('/ensemble/barostat/stress2', value2, lerr2)	!!xy
          call parse_get('/ensemble/barostat/stress3', value3, lerr3)	!!xz
          call parse_get('/ensemble/barostat/stress4', value4, lerr4)	!!yy
          call parse_get('/ensemble/barostat/stress5', value5, lerr5)	!!yz
          call parse_get('/ensemble/barostat/stress6', value6, lerr6)	!!zz
          if (.not. (lerr1 .or. lerr2 .or. lerr3 .or. lerr4 .or. lerr5 .or. lerr6) ) then
             call fmod_sigmaS(value1,value2,value3,value4,value5,value6)
          else
             sigmaS = 0.d0
          endif
          
       elseif (trim(chara) == 'opt') then
          call fmod_md_condition__ensemble(OPT)
          if (velocity_scaling_flag) then
             call parse_abort('Velocity scaling cannot be used in OPT.')
          endif
          getname='/integrator/optimize/step_length'
          call parse_get(trim(getname),value1,lerr)
          value1 = value1 * 1.0d-10
          if (.not. lerr) then
             call fmod_opt_condition__step_length(value1)
          endif
          getname='/integrator/optimize/convergence'
          call parse_get(trim(getname),value1,lerr)
          if (.not. lerr) then
             call fmod_opt_condition__convergence(value1)
          endif
          getname='/integrator/optimize/up_rate'
          call parse_get(trim(getname),value1,lerr)
          if (.not. lerr) then
             call fmod_opt_condition__step_lurate(value1)
          endif
          getname='/integrator/optimize/down_rate'
          call parse_get(trim(getname),value1,lerr)
          if (.not. lerr) then
             call fmod_opt_condition__step_ldrate(value1)
          endif
       endif
       
       call parse_get('/integrator/steps', ivalue1, lerr)
       if (.not. lerr)  call fmod_howmany_steps(ivalue1)
       call parse_get('/ensemble/temperature', value1, lerr)
       if (.not. lerr)  call fmod_systemp(value1)
       call parse_get('/ensemble/pressure', value1, lerr)
       if (.not. lerr)  call fmod_syspres(value1)
! TZYBEGIN
       getname='/integrator/multiple time step/nstep_skip_middle'
       call parse_get(trim(getname), ivalue1, lerr)
      !  if (.not. lerr)  call fmod_set_maxMTm(ivalue1)
       if (.not. lerr) then
           call fmod_set_maxMTm(ivalue1)
           if( myrank == 0 ) write(*,*) "nstep_skip_middle is deprecated. Use nshort_per_middle instead"
       endif
       getname='/integrator/multiple time step/nstep_skip_long'
       call parse_get(trim(getname), ivalue1, lerr)
      !  if (.not. lerr)  call fmod_set_maxMTl(ivalue1)
       if (.not. lerr) then
           call fmod_set_maxMTl(ivalue1)
           if( myrank == 0 ) write(*,*) "nstep_skip_long is deprecated. Use nmiddle_per_long instead"
       endif


       getname='/integrator/multiple time step/nshort_per_middle'
       call parse_get(trim(getname), ivalue1, lerr)
       if (.not. lerr)  call fmod_set_maxMTm(ivalue1)
       getname='/integrator/multiple time step/nmiddle_per_long'
       call parse_get(trim(getname), ivalue1, lerr)
       if (.not. lerr)  call fmod_set_maxMTl(ivalue1)
! TZYEND
       call parse_get('/ensemble/Pref_inner', value1, lerr)
       if (.not. lerr)  call fmod_set_Pinner(value1)
       call parse_get('/ensemble/Pref_outer', value1, lerr)
       if (.not. lerr)  call fmod_set_Pouter(value1)
       call parse_get('/ensemble/Pref_outermost', value1, lerr)
       if (.not. lerr)  call fmod_set_Poutermost(value1)
       call parse_get('/ensemble/maxwell_velocities', chara, lerr)
       if(trim(chara)=='yes') then
          call fmod_reset_maxwell(1)
       endif
       ! mpi
       call parse_get('/mpi/division', chara, lerr)
       if (.not. lerr .and. trim(chara) == 'manual') then
          mpi_manual_division_flg = .true.
       endif
       call parse_get('/mpi/nxdiv', ivalue1, lerr)
       if (.not. lerr)  call set_npx(ivalue1)
       call parse_get('/mpi/nydiv', ivalue1, lerr)
       if (.not. lerr)  call set_npy(ivalue1)
       call parse_get('/mpi/nzdiv', ivalue1, lerr)
       if (.not. lerr)  call set_npz(ivalue1)
       ! shake
       call parse_get('/integrator/shake/maxiteration', ivalue1,lerr)
       if (.not. lerr)  call fmod_shake_max_iteration(ivalue1)
       call parse_get('/integrator/shake/shake_tolerance', value1,lerr)
       if (.not. lerr)  call fmod_md_shake__shake_tolerance(value1)
       
       ! intermolecular interaction
       getname='/intermolecular interaction/twobody/cutoff'
       call parse_get(trim(getname), value1, lerr)
       value1 = value1 * 1.0d-10
       value2 = value1 * value1
       if (.not. lerr) then
          call fmod_cutrad(value1)
          call fmod_cutrad2(value2)
       endif
! TZYBEGIN
       getname='/intermolecular interaction/non_even_dihedral_potential'
       call parse_get(trim(getname), chara, lerr)
       if(trim(chara)=='yes') then
          call fmod_md_condition__bNonEvenDihedralPotential(1)
       elseif(trim(chara)=='no') then
          call fmod_md_condition__bNonEvenDihedralPotential(0)
       elseif(lerr) then
          call fmod_md_condition__bNonEvenDihedralPotential(-1)
       endif
! TZYEND
       getname='/intermolecular interaction/twobody/LJcorrection'
       call parse_get(trim(getname), chara, lerr)
       if(trim(chara)=='yes') then
          call fmod_LJ_correction(1)
       elseif(trim(chara)=='no') then
          call fmod_LJ_correction(0)
       elseif(lerr) then
          call fmod_LJ_correction(-1)
       endif
       call parse_get('/intermolecular interaction/type', chara, lerr)
       if(trim(chara) == 'cluster') then
          write(0,*) 'ERROR: type=cluster is not supported'
          call modylas_abort()
          call fmod_md_periodic__type(CLUSTER)
       elseif(trim(chara) == 'simple') then
          write(0,*) 'ERROR: type=simple is not supported'
          call modylas_abort()
          call fmod_md_periodic__type(SIMPLE)
       elseif(trim(chara) == 'ewald') then
          call fmod_md_periodic__type(EWALD)
          getname='/intermolecular interaction/ewald/alpha'
          call parse_get(trim(getname),value1,lerr)
          if (.not. lerr) then
             call fmod_md_ewald__alpha(value1)
             call check_pmewald_alpha(value1)
          endif
          getname='/intermolecular interaction/ewald/h2max'
          call parse_get(trim(getname),ivalue1,lerr)
          if (.not. lerr)  call fmod_md_ewald__max_h2(ivalue1)
       elseif(trim(chara) == 'pme') then
          call fmod_md_periodic__type(PMEWALD)
          getname='/intermolecular interaction/pme/alpha'
          call parse_get(trim(getname),value1,lerr)
          if (.not. lerr) then
             call fmod_md_ewald__alpha(value1)
             call check_pmewald_alpha(value1)
          endif
          getname='/intermolecular interaction/pme/bsorder'
          call parse_get(trim(getname),value1,lerr)
          if (.not. lerr)  call fmod_pmewald_bsorder(value1)
          getname='/intermolecular interaction/pme/nfft1'
          call parse_get(trim(getname),ivalue1,lerr)
          if (.not. lerr)  call fmod_pmewald_nfft1(ivalue1)
          getname='/intermolecular interaction/pme/nfft2'
          call parse_get(trim(getname),ivalue2,lerr)
          if (.not. lerr)  call fmod_pmewald_nfft2(ivalue2)
          getname='/intermolecular interaction/pme/nfft3'
          call parse_get(trim(getname),ivalue3,lerr)
          if (.not. lerr)  call fmod_pmewald_nfft3(ivalue3)
          
       elseif(trim(chara) == 'fmm') then
          call fmod_md_periodic__type(FMM)
          getname='/intermolecular interaction/fmm/nmax'
          call parse_get(trim(getname),ivalue1,lerr)
          if (.not. lerr)    call fmod_set_nmax(ivalue1)
          getname='/intermolecular interaction/fmm/ULswitch'
          call parse_get(trim(getname),ivalue1,lerr)
          if (.not. lerr)    call fmod_set_lgflg(ivalue1)
          getname='/intermolecular interaction/fmm/sterm'
          call parse_get(trim(getname),chara,lerr)
          if(trim(chara)=='yes') then
             call fmod_ewald_sterm(1)
          elseif(trim(chara)=='no'.or.lerr) then
             call fmod_ewald_sterm(0)
          endif
       endif
       
       getname='/intermolecular interaction/fmm/nlevel'
       call parse_get(trim(getname),ivalue1,lerr)
       if (.not. lerr)    call fmod_set_fmm_nlevel(ivalue1)
       getname='/intermolecular interaction/ncellx'
       call parse_get(trim(getname),ivalue1,lerr)
       if (.not. lerr)    call fmod_set_fmm_ncellx(ivalue1)
       getname='/intermolecular interaction/ncelly'
       call parse_get(trim(getname),ivalue1,lerr)
       if (.not. lerr)    call fmod_set_fmm_ncelly(ivalue1)
       getname='/intermolecular interaction/ncellz'
       call parse_get(trim(getname),ivalue1,lerr)
       if (.not. lerr)    call fmod_set_fmm_ncellz(ivalue1)
       !old-style
       getname='/intermolecular interaction/ncell'
       call parse_get(trim(getname),ivalue1,lerr)
       if (.not. lerr) then
          call warn_old_ncell
          call fmod_set_ncell(ivalue1)
          call fmod_set_fmm_nlevel( calculate_fmm_nlevel(ivalue1) )
       endif
       !old-style
       
       
       !ay     hidden command
       lerr_bk = lerr
       call parse_get('/forcefield/cmap_version',ivalue1,lerr)
       if(.not.lerr) then
          call fmod_set_CMAPVersion(ivalue1)
       endif
       lerr = lerr_bk

       !ya     hidden command
       call parse_get('/debug/0step_stop',chara,lerr)
       if (.not. lerr .and. trim(chara) == 'yes') then
           call fmod_set_0step_stop
       endif
       
          !ya     hidden command
          call parse_get('/allocation parameters/totnconstL',ivalue1,lerr)
          if (.not. lerr) then
             call fmod_set_totnconstL(ivalue1)
          endif
          call parse_get('/allocation parameters/na1cell',ivalue1,lerr)
          if (.not. lerr) then
             call fmod_set_na1cell(ivalue1)
          endif
          call parse_get('/allocation parameters/max_nsegments_per_cell', &
     &                    ivalue1,lerr)
          if (.not. lerr) then
             call fmod_set_maxnsegmentspercell(ivalue1)
          endif

#ifdef SEGSHAKE
          call parse_get('/COM/constrain_COM',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             call fmod_constrain_COM(1)
          else
             call fmod_constrain_COM(0)
          endif
          call parse_get('/COM/groupAtop', ivalue1, lerr1)
          call parse_get('/COM/groupAend', ivalue2, lerr2)
          if(.not. lerr1 .and. .not. lerr2) then
             call fmod_set_groupA(ivalue1, ivalue2)
          endif
          call parse_get('/COM/groupBtop', ivalue1, lerr1)
          call parse_get('/COM/groupBend', ivalue2, lerr2)
          if(.not. lerr1 .and. .not. lerr2) then
             call fmod_set_groupB(ivalue1, ivalue2)
          endif
          call parse_get('/COM/dist_COM', qvalue1, lerr1)
          if(.not. lerr1) then
             value1=qvalue1
             call fmod_set_dist_COM(value1)
          endif
          call parse_get('/COM/change_distCOM',chara,lerr)
          if (.not. lerr .and. trim(chara) == 'yes') then
             call fmod_change_distCOM(1)
          endif
          call parse_get('/COM/deltaR', qvalue1, lerr1)
          if(.not. lerr1) then
             value1=qvalue1
             call fmod_set_d_COMR(value1)
          endif
#endif
#ifdef TIP4
        call parse_get('/tip4p/type',chara,lerr)
        if (.not. lerr)then
        if(trim(chara) == 'original') then
          call fmod_set_tip4p_geometry(0)
        elseif(trim(chara) == '2005') then
          call fmod_set_tip4p_geometry(1)
        else
          write(*,*) 'keyword in <tip4p> type= </tip4p> is wrong', &
     &    ', it should be original or 2005.'
          call mpiend
        endif
        endif
!
        call parse_get('/tip4p/msite',chara,lerr)
        if (.not. lerr)then
          if(trim(chara) == 'modylas') then
            call fmod_set_msite_position(0)
          elseif(trim(chara) == 'gromacs') then
            call fmod_set_msite_position(1)
          else
            write(*,*) 'keyword in <tip4p> msite= </tip4p> is wrong', &
     &      ', it should be modylas or gromacs.'
            call mpiend
          endif
        endif
#endif 
 
!
! closing .mddef
!

       call parse_close(lerr)
       call fmod_md_check_scaling_factor_is_set()
       if(myrank==0) write(*,*) 'Reading .mddef file ended successfully!'
       
       if (lerr) then
          call parse_abort(parse_show_error_message())
       endif
    endif
    !     Show error message except parse_close().
    if (lerr) then
       call parse_abort(parse_show_error_message())
    endif
  end subroutine parse_input_v1_0_0
!-------------------------------------------------------------------------
!>
!! \brief  Subroutine to read mdxyz without parse module.
!! \author Kensuke Iwahashi
!<
  subroutine direct_read(filename, read_velocity)
    use trajectory_org
    implicit none
    character(len=*):: filename
    logical:: read_velocity
    character(len=80):: line
    real(8):: val1, val2, val3
    integer(4):: istatus, ia, iv
    logical:: atom, positions, velocities
    include 'mpif.h'
    integer(4) :: ierr
    atom = .false.
    positions = .false.
    velocities = .false.
    ia = 0
    iv = 0
    open(99, file=filename, status='old')
    istatus = 0
    do while(istatus == 0)
       read(99, '(a80)', iostat=istatus) line
       if (index(line, '<atom>') > 0) then
          atom = .true.
       else if (index(line, '</atom>') > 0) then
          atom = .false.
       else if (atom .and. index(line, '<positions>') > 0) then
          positions = .true.
       else if (atom .and. index(line, '</positions>') > 0 ) then
          positions = .false.
       else if (atom .and. index(line, '<velocities>') > 0 ) then
          if (read_velocity)  velocities = .true.
       else if (atom .and. index(line, '</velocities>') > 0) then
          velocities = .false.
       else if (positions) then
          ia = ia + 1
          read(line, *) val1, val2, val3
          xyz(1,ia) = val1 * 1.0d-10
          xyz(2,ia) = val2 * 1.0d-10
          xyz(3,ia) = val3 * 1.0d-10
       else if (velocities) then
          iv = iv + 1
          read(line, *) val1, val2, val3
          v(1,iv) = val1 * 1.0d-10
          v(2,iv) = val2 * 1.0d-10
          v(3,iv) = val3 * 1.0d-10
       endif
    enddo
    close(99)
    if (ia > 0) then
       call MPI_Bcast(xyz, 3*ia, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    endif
    if (iv > 0) then
       call MPI_Bcast(v,   3*iv, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    endif
  end subroutine direct_read
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to validate MD inputs.
!! \author Kensuke Iwahashi
!<
  subroutine init_md_check
    use trajectory_org
    use cutoff_radius
    use md_condition
    use md_periodic
    use molecules
    use mpi_tool
    use unit_cell
    use ensemble_numbers
    use file_restart, only : restart_interval
    use domain, only : ncellx,ncelly,ncellz
    implicit none
    integer(4) :: type
    include 'mpif.h'
    if(myrank.eq.mpiout) then
       ! check MPI process number
       if(nprocs .gt. ncellx*ncelly*ncellz)then
          write(0,*) 'ERROR: NPROCS=',nprocs,  &
     &               '.gt. total number of subcells:', &
     &                ncellx*ncelly*ncellz
          write(0,*) 'Select smalller NPROCS, or greater ncell.   '
          call modylas_abort
       endif

       if (md_periodic__howmany_molecules < 0) then
          write(0,*) 'ERROR: md_periodic.howmany_molecules must be described.'
          call modylas_abort
       endif
       !
       type = 0
       if(md_periodic__type == CLUSTER) then
          if(md_condition__ensemble == NVE) then
             type = 1
          endif
       else if (md_periodic__type == FMM) then
          if(md_condition__ensemble /= NPT_PR) then
             type = 1
          endif
       else
          if(md_condition__ensemble == NVE) then
             type = 1
          else if(md_condition__ensemble == NVT) then
             type = 1
          else if(md_condition__ensemble == NVLXLYLZT) then
             type = 1
          else if(md_condition__ensemble == NPT_A) then
             type = 1
          else if(md_condition__ensemble == NPT_Z) then
             type = 1
          else if(md_condition__ensemble == NPTLZ) then
             type = 1
          else if(md_condition__ensemble == NLXLYPZT) then
             type = 1
          else if(md_condition__ensemble == NPXPYPZT) then
             type = 1
          else if(md_condition__ensemble == NLXLYLZT) then
             type = 1
          else if(md_condition__ensemble == NPTLZxy) then
             type = 1
          else if(md_condition__ensemble == NPT_PR) then
             type = 1
          else if(md_condition__ensemble == OPT) then
             type = 1
          endif
       endif
       if (type == 0) then
          write(0,*) 'ERROR: Bad condition.'
          call modylas_abort
       endif
       !
       !     check if the last trajectory will be saved
       if (mod(md_condition__howmany_steps,restart_interval) /= 0) then
          write(0,*) 'ERROR: the last trajectory will not be saved'
          call modylas_abort
       endif
       !
       !     check periodic boundary
       if(md_periodic__type /= CLUSTER) then

          !       check validity of cut-off length
          if(cutrad > cellxh) then
             write(0,*)  'ERROR: cut-off length for force is greater than a half of cell'
             call modylas_abort
          endif
          if(cutrad > cellyh) then
             write(0,*)  'ERROR: cut-off length for force is greater than a half of cell'
             call modylas_abort
          endif
          if(cutrad > cellzh) then
             write(0,*)  'ERROR: cut-off length for force is greater than a half of cell'
             call modylas_abort
          endif
       endif
    endif
  end subroutine init_md_check

  subroutine warn_old_ncell
    use domain, only : ncellx_input, ncelly_input, ncellz_input, nlevel_input
    use mpi_tool, only : myrank
    implicit none

!   if(.not.(ncellx_input .and. ncelly_input .and. ncellz_input .and. nlevel_input)) then
!      if(myrank==0) then
!         write(*,*) 'WARNING: Only ncell was read from *.mddef, ', &
!              & 'as old modylas for NPROCS=2powers'
!      endif
!   endif
  end subroutine warn_old_ncell

  ! This function finds the positive exponent m such that power^m * const = n.
  pure elemental integer(4) function find_exponent(n_in, power)
    integer(4), intent(in) :: n_in, power
    integer(4) :: n

    n = n_in
    find_exponent = 0

    do while (mod(n,power) == 0)
       n = n / power
       find_exponent = find_exponent + 1
    enddo
  end function find_exponent

  pure elemental integer(4) function calculate_fmm_nlevel(ncell)
    implicit none
    integer(4), intent(in) :: ncell
    
    calculate_fmm_nlevel = find_exponent(ncell, 2) + find_exponent(ncell, 3)
  end function calculate_fmm_nlevel
  
end module input_mod
