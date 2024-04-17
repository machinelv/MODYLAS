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
!! \brief  Subroutines for prolog, epilog and file I/O of MODYLAS.
!! \author Kensuke Iwahashi
!<
!----------------------------------------------------------------------
!>
!! \brief  Subroutine for initialization of MODYLAS.
!! \author Kensuke Iwahashi
!<
  subroutine initialize_application
    use md_condition
    use device_numbers
    use md_periodic
    use lj_mod
    use mpi_tool
    use openmp_tool, only : nomp
    use comm_bound_mod
    use mpi_3d_grid
    use domain
    use fmm_far
    use ewald_mod
    use pme_far
    use kinetic_energy
    use CMAP, only : PreCMAP, ncmap
    use shake_rattle_roll
    use thermostat
    use barostat
    use position_constrain
    use ensemble_numbers
    use segments, only : init_segment_check
    use cell_shape, only : init_md_periodic
    use atom_virial
    use file_application
    use input_mod, only : init_md_check
    use session_name_mod
    use file_utility
    use dr_cntl
    use comm_direct2_dr
    use comm_direct3_dr
    use fmm_near_dr
#ifdef TIP4
    use tip4p, only : init_tip4p
#endif
    use center_of_mass
    use table_functions, only : init_table
    implicit none

    ! check preprecessors
#ifdef CHARMMFSW
    if(myrank==0) write(*,*) "-DCHARMMFSW is set"
#endif
#ifdef OPLSAMBER
    if(myrank==0) write(*,*) "-DOPLSAMBER is set"
#endif
#ifdef GAFF
    if(myrank==0) write(*,*) "-DGAFF is set"
#endif
#ifdef ONEPROC_AXIS
    if(myrank==0) write(*,*) "-DONEPROC_AXIS is out of support."
    call mpiend
#endif
#ifdef TIP4
    if(myrank==0) write(*,*) "-DTIP4 is set"
#endif
#ifdef SYNC_COM
    if(myrank==0) write(*,*) "-DSYNC_COM is set"
#endif
#ifdef GMXFORT
    if(myrank==0) write(*,*) "-DGMXFORT is set"
#endif
#ifdef XTC
    if(myrank==0) write(*,*) "-DXTC is set"
#endif
#ifdef FJ_RDMA
    if(myrank==0) write(*,*) "-DFJ_RDMA is set"
#endif
#ifdef PROFILE_COMM
    if(myrank==0) write(*,*) "-DPROFILE_COMM is set"
#endif

!>> for debug
#ifdef NOTABLE
    if(myrank==0) write(*,*) "-DNOTABLE is set"
#endif
#ifdef DEBUGFCE
    if(myrank==0) write(*,*) "-DDEBUGFCE is set"
#ifdef KCAL
    if(myrank==0) write(*,*) "-DKCAL is set"
#endif
#ifdef ATM
    if(myrank==0) write(*,*) "-DATM is set"
#endif
#endif
#ifdef DEBUG_MTDFMM
    if(myrank==0) write(*,*) "-DDEBUG_MTDFMM is set"
#endif
#ifdef DEBUG_RDMAYA
    if(myrank==0) write(*,*) "-DDEBUG_RDMAYA is set"
#endif
!<< for debug

#ifdef PRECISION_P2P_MIX
    if(myrank==0) write(*,*) "-DPRECISION_P2P_MIX is set"
#elif defined (PRECISION_P2P_SP)
    if(myrank==0) write(*,*) "-DPRECISION_P2P_SP is set"
#endif
#ifdef MOMENT
    if(myrank==0) write(*,*) "-DMOMENT is set"
#endif

#ifdef PRECISION_M2L_MIX
    if(myrank==0) write(*,*) "-DPRECISION_M2L_MIX is set"
#elif defined(PRECISION_M2L_SP)
    if(myrank==0) write(*,*) "-DPRECISION_M2L_SP is set"
#endif

#ifdef GROEXT
    if(myrank==0) write(*,*) "-DGROEXT is set"
#endif

#ifdef PROF_TIMER
    if(myrank==0) write(*,*) "-DPROF_TIMER is set"
#endif

    !     initialize MD objects
    call init_g_main
    call init_md_check
    call init_segment_check
#ifdef TIP4
    call init_tip4p
#endif
#ifdef SEGSHAKE
    call init_center_of_mass
#endif
    if (.not. exist_mdffbin) then
       call init_md_condition
    endif
    call init_md_velocity
    call init_md_periodic
    call init_mpi_3d_grid
    call init_fmm_domain_div
    call check_parallel_condition !check nprocs and nomp
    call check_cutofflength
    call fmod_set_maxsegments
    call fmod_alloc_metadata
    call fmod_alloc_kinetic_energy
    call fmod_alloc_atom_virial
    call fmod_alloc_smp_wk

    if (   md_periodic__type == FMM) then
       call fmod_alloc_multipole
       call init_fmm
    elseif(md_periodic__type == PMEWALD)then
       call fmod_alloc_pme_arraies
       call init_md_pmewald
    elseif(md_periodic__type == EWALD)then
       call fmod_alloc_ewald_arraies
       call init_md_ewald
    endif
    call init_position_constrain
    call init_comm_bound()
    call init_drtbl
    call init_comm_direct_2_dr() !! DR method
    call init_comm_direct_3_dr() !! DR method
    call init_energy_direct_dr() !! DR method
    call init_shake_local()
    call pshake_initialize1
    call pshake_initialize2
    call pshake_finish1
    if(ncmap .ne. 0) call PreCMAP()
#ifdef CHARMMFSW
    call init_charmmfsw
#endif
    call init_table

    if(md_condition__ensemble==NVT  .or. &
         &   md_condition__ensemble==NPT_A.or. &
         &   md_condition__ensemble==NPT_Z.or. &
         &   md_condition__ensemble==NPTLZ.or. &
         &   md_condition__ensemble==NLXLYPZT.or. &
         &   md_condition__ensemble==NPXPYPZT.or. &
         &   md_condition__ensemble==NLXLYLZT.or. &
         &   md_condition__ensemble==NVLXLYLZT.or. &
         &   md_condition__ensemble==NPTLZxy.or. &
         &   md_condition__ensemble==NPT_PR) then
       call init_nhc()
       if(initialize_thermostat_flag) call zero_thermostat
    endif
    if(md_condition__ensemble==NPT_A.or. &
         &   md_condition__ensemble==NPT_Z.or. &
         &   md_condition__ensemble==NPTLZ.or. &
         &   md_condition__ensemble==NLXLYPZT.or. &
         &   md_condition__ensemble==NPXPYPZT.or. &
         &   md_condition__ensemble==NLXLYLZT.or. &
         &   md_condition__ensemble==NPTLZxy.or. &
         &   md_condition__ensemble==NPT_PR) then
       call init_nhca()
       if(initialize_thermostat_flag) call zero_thermostat
       if(initialize_barostat_flag)   call zero_barostat
    endif
    call md_calc_corrector_constant()

! limitation in modylas ver 1.1.0
    if (    md_periodic__type == FMM) then
       if(  md_condition__ensemble==NPT_Z.or. &
         &   md_condition__ensemble==NPTLZ.or. &
         &   md_condition__ensemble==NLXLYPZT.or. &
         &   md_condition__ensemble==NPXPYPZT.or. &
         &   md_condition__ensemble==NLXLYLZT.or. &
         &   md_condition__ensemble==NPTLZxy.or. &
         &   md_condition__ensemble==NPT_PR) then
    if(myrank==0)then
       write(0,'(a,i0,a)')  &
       &    'ERROR: When type=fmm is selected,   ' // &
       &    'the selectable ensemble=(opt, nve, nvt, npt_a).'
    endif
       call mpiend()
       endif
    endif
!

  end subroutine initialize_application
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to allocate arrays in meta-data style (FMM)
!! \author Yoshimichi Andoh
!<
  subroutine fmod_alloc_metadata
!----------------------------------------------------------------------
    use omp_lib
    use trajectory_org
    use trajectory_mpi
    use forces
    use md_forces
    use subcell
    use domain
    use shake_rattle_roll
    use md_multiplestep
    use opt_mod, only : opt_direction
    use ensemble_numbers
    use md_condition
    use bond_breaking_event, only: bond_breaking_event_allocate
    use margin_sizes
    use openmp_tool, only : nomp
    use void123, only : void_numL,void_topL,void_list_j0, &
                        maxj0,maxj0_per_i0
    use dr_cntl, only : nbd,nbd2
#ifdef DIHEDRAL_TABLE
    use dihedral, only : kblk, k0_w, k0_n, atom1B, jsize,jsize2, &
        w3_f1,w3_f2,w3_f3,fi0_2,fi0_3,jlist,i0a_2,i0c_2,i0d_2,i0_2 &
    & , atom1B_t
#else
    use dihedral, only : kblk, k0_w, k0_n, atom1B, jsize,jsize2, &
        w3_f1,w3_f2,w3_f3,fi0_2,fi0_3,jlist,i0a_2,i0c_2,i0d_2,i0_2
#endif

    implicit none
    integer(4) :: itmp
    include 'mpif.h'
#ifdef FJ_RDMA
    type(c_ptr) :: wkxyz_cptr
#endif

!   if(myrank==0) write(*,*) "Starting allocation of many arrays."
!   flush(6)
    !############
    !  metadata
    !############
    allocate(tag(0:lzdiv+nbd,0:lydiv+nbd,0:lxdiv+nbd))
    allocate(na_per_cell(0:lzdiv+nbd,0:lydiv+nbd,0:lxdiv+nbd))
    !############
    !  segment
    !############
    max_seg = max_nsegments_per_cell*lzdiv*lydiv*lxdiv
    allocate(wseg_cz(max_seg))
    allocate(wseg_cy(max_seg))
    allocate(wseg_cx(max_seg))
    allocate(ndseg_fmmn(lzdiv,lydiv,lxdiv))
    itmp=max_nsegments_per_cell*(lxdiv+nbd2)*(lydiv+nbd2)*(lzdiv+nbd2)
    allocate(lsegtop(itmp))
    allocate(lseg_natoms(itmp))
    !############
    !  atom
    !############
    itmp=ncellx*ncelly*ncellz
    na1cell=max( int(int(n/itmp)*na1cellmargin),na1cell_init)
    na5cell=na1cell*(lzdiv+nbd2)
    nalongcell=na1cell*(lzdiv)
    nadirect=na1cell*(lxdiv+nbd2)*(lydiv+nbd2)*(lzdiv+nbd2)
    narea   =na1cell*(lzdiv+nbd2)*(lydiv+nbd2)
    naline  =na1cell*(lzdiv+nbd2)
    !Coordinate & Velocity
    allocate(wkxyz(3,nadirect))
#ifdef FJ_RDMA
    call rdma_register_addr(wkxyz, (3*nadirect)*8)
    wkxyz_laddr = rdma_get_laddr(wkxyz)
    wkxyz_cptr  = rdma_get_raddr(wkxyz)
    call c_f_pointer(wkxyz_cptr, fptr=wkxyz_raddr, shape=[nprocs])
#endif
    allocate(wkv(3,nadirect))
    allocate(m2i(nadirect))
    !Force
    allocate(wk_f(3,nadirect))
#ifdef FJ_RDMA
    call rdma_register_addr(wk_f, (3*nadirect)*8)
#endif
    allocate(w3_f(3,nadirect,0:nomp-1))
    !void
    allocate(void_numL(nadirect))
    allocate(void_topL(nadirect))
    maxj0=nadirect*maxj0_per_i0
    allocate(void_list_j0(maxj0))
    void_numL=0 !! initial value is very important. Never delete these lines
    void_topL=0
    void_list_j0=0
    !TABLE for Force_calculation
    allocate(chgv_table(na5cell,0:nomp-1))
    allocate(epsilon_sqrt_table(na5cell,0:nomp-1))
    allocate(R_half_table(na5cell,0:nomp-1)) !! CHARMM
    !MT
    allocate(fshort(3,nadirect))
    allocate(fmiddle(3,nadirect))
#ifdef FJ_RDMA
    call rdma_register_addr(fshort, (3*nadirect)*8)
    call rdma_register_addr(fmiddle, (3*nadirect)*8)
#endif
    allocate(flong(3,nadirect))
    !SHAKE
    allocate(wk_dfc(3,nadirect))
    allocate(xyzstr(3,nadirect))
    allocate(kshake(nadirect))
    allocate(vstr(3,nadirect))
    !Optimize
    if(md_condition__ensemble==OPT) then  ! optimize
       allocate(opt_direction(3,nadirect))
    endif

!dihedral
    itmp=max_nsegments_per_cell*(lxdiv+nbd2)*(lydiv+nbd2)*(lzdiv+nbd2)
    kblk=itmp
    k0_w=(kblk-1)/nomp + 1     
    k0_n=(kblk-1)/k0_w + 1 
!
#ifdef DIHEDRAL_TABLE
    jsize=size(atom1B_t(:,:),2) 
#else
    jsize=size(atom1B(:,:),2) 
#endif
    jsize2=nadirect
    kblk=k0_n  
    jsize=jsize*jsize2 / k0_n + 1
!
!   if(myrank==0) write(*,*) "Allocating potentially huge arrays called w3_f1, w3_f2, and w3_f3"
!   flush(6)
! TZYBEGIN
    ! allocate(w3_f1(3,jsize,kblk,0:nomp-1)) 
    ! allocate(w3_f2(3,jsize,kblk,0:nomp-1)) 
    ! allocate(w3_f3(3,jsize,kblk,0:nomp-1)) 
    allocate(w3_f1(3,jsize,0:nomp-1)) 
    allocate(w3_f2(3,jsize,0:nomp-1)) 
    allocate(w3_f3(3,jsize,0:nomp-1)) 
! TZYEND
!   if(myrank==0) write(*,*) "Allocation of  w3_f1, w3_f2, and w3_f3 has been performed"
!   flush(6)
    allocate(fi0_2(3,jsize2)) 
    allocate(fi0_3(3,jsize,kblk)) 
    allocate(jlist(jsize,kblk)) 
    allocate(i0a_2(jsize,kblk)) 
    allocate(i0c_2(jsize,kblk)) 
    allocate(i0d_2(jsize,kblk)) 
    allocate(i0_2(jsize,kblk)) 
!dihedral

if(allow_bond_breaking>0.0d0) call bond_breaking_event_allocate()
  end subroutine fmod_alloc_metadata
