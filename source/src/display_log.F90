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
!! \brief  Module and subroutines to output calculation log.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to output log.
!! \author Yoshimichi Andoh
!<
module display_log_mod
  use omp_lib

contains

!>
!! \brief  Subroutine to output log.
!! \author Yoshimichi Andoh
!<
  subroutine output_log   ! output log
    use maxwell_distribution, only : reset_maxwell
    use version
    use mpi_tool
    use trajectory_org
    use param
    use pshake_init
    use CMAP
    use cutoff_radius
    use shake_rattle_roll
    use extended_system
    use file_dcd
    use device_numbers
    use session_name_mod
    use void123, only : nvoid
    use special14, only : nljsp, ncoulombsp
    use md_multiplestep
    use md_condition
    use md_periodic
    use ewald_mod
    use ewald_variables
    use fmm_parameters, only : nmax
    use fmm_far, only : lgflg
    use mpi_3d_grid, only : npx, npy, npz, mpi_manual_division_flg
    use domain, only : lxdiv, lydiv, lzdiv, ncellx, ncelly, ncellz, nlevel
    use comm_pme_mod
    use md_oplsaa_special_divide
    use file_mdmntr
    use file_mdtrj
    use ensemble_numbers
    use force_field_numbers
    use file_restart, only : restart_interval
#ifdef CHARMMFSW
    use lj_mod, only : Ron, Roff
#endif
#ifdef TIP4
    use tip4p, only : msite, type_of_msiteposition
#endif
    implicit none
    include 'mpif.h'
    integer(4) :: ierr
    real(8)::derfc
    integer(4)::nomp=1

!$  nomp = omp_get_max_threads()

    IF(myrank==0)THEN
       write(*,'(a)') ' '
       write(*,'(a,a,a)') &
            & '## ', trim(session_name), '.log   ' // &
            & ' -- MD conditions output from MD calculation by modylas'
       !version
       write(*,'(a)') '#'
       write(*,'(a)') '# VERSION'
       write(*,'(a)') '#'
       write(*,'(a,a10)') 'MODYLAS version  :', trim(MODYLAS_version)
       write(*,'(a,a10)') 'input   version  :', trim(input_version)
       !parameters
       write(*,'(a)') '#'
       write(*,'(a)') '# PARAMETERS'
       write(*,'(a)') '#'
       write(*,'(a,i10)') 'natom            :', n
#ifdef TIP4
       write(*,'(a,i10)') 'nmsite(TIP4P)    :', msite
       if(    type_of_msiteposition==0)then
       write(*,'(a,a13)') 'place of msite   :','MODYLAS-style'
       elseif(type_of_msiteposition==1)then
       write(*,'(a,a13)') 'place of msite   :','GROMACS-style'
       endif
#endif
       write(*,'(a,i10)') 'nmolecule        :', nmols
       write(*,'(a,i10)') 'npara            :', npara
       write(*,'(a,i10)') 'nvoid            :', nvoid
       write(*,'(a,i10)') 'nljsp            :', nljsp
       write(*,'(a,i10)') 'nclsp            :', ncoulombsp
       write(*,'(a,i10)') 'nshake           :', totnconst
       write(*,'(a,i10)') 'nshake_group     :', n_type_ps
       write(*,'(a,i10)') 'pshake_group     :', rngrp_ps
       write(*,'(a,i10)') 'Deg. of freedom  :', degree_of_freedom
       !integration
       write(*,'(a)') '#'
       write(*,'(a)') '# INTEGRATOR'
       write(*,'(a)') '#'
       write(*,'(a,i10)') 'MD_steps         :', md_condition__howmany_steps
!      write(*,'(a,f10.1,a)') 'delta t          :', dt*1d+15, ' [fs]'
       write(*,'(a,f10.1,a)') 'delta t (long)   :', dt*1d+15, ' [fs]'
       write(*,'(a,f10.1,a)') 'delta t (middle) :', dt*1d+15/maxMTl, ' [fs]'
       write(*,'(a,f10.1,a)') 'delta t (short)  :', dt*1d+15/maxMTm/maxMTl, ' [fs]'
!      write(*,'(a,i10)') 'nskip_middle     :', maxMTm
!      write(*,'(a,i10)') 'nskip_long       :', maxMTl
       write(*,'(a,es10.1)') 'shake_torelance  :', shake_tolerance
       write(*,'(a,es10.1)') 'roll_torelance   :', roll_tolerance
       write(*,'(a,i10)') 'mntr_interval    :',mntr_interval
       write(*,'(a,i10)') 'trj_interval     :',trj_interval
       if(dcd_output)then
          write(*,'(a,i10)') 'dcd_interval     :',dcd_interval
       endif
       write(*,'(a,i10)') 'restart_interval :',restart_interval
       !ensemble
       write(*,'(a)') '#'
       write(*,'(a)') '# ENSEMBLE'
       write(*,'(a)') '#'
       if(    md_condition__ensemble==NVE) then
          write(*,'(a,a10)') 'ensemble         :', 'NVE'
       elseif(md_condition__ensemble==NVT) then
          write(*,'(a,a10)') 'ensemble         :', 'NVT'
       elseif(md_condition__ensemble==NVLXLYLZT) then
          write(*,'(a,a10)') 'ensemble         :', 'NVLXLYLZT'
       elseif(md_condition__ensemble==NPT_A) then
          write(*,'(a,a10)') 'ensemble         :', 'NPT_A'
       elseif(md_condition__ensemble==NPT_Z) then
          write(*,'(a,a10)') 'ensemble         :', 'NPT_Z'
       elseif(md_condition__ensemble==NPTLZ) then
          write(*,'(a,a10)') 'ensemble         :', 'NPTLZ'
       elseif(md_condition__ensemble==NLXLYPZT) then
          write(*,'(a,a10)') 'ensemble         :', 'NLXLYPZT'
       elseif(md_condition__ensemble==NPXPYPZT) then
          write(*,'(a,a10)') 'ensemble         :', 'NPXPYPZT'
       elseif(md_condition__ensemble==NLXLYLZT) then
          write(*,'(a,a10)') 'ensemble         :', 'NLXLYLZT'
       elseif(md_condition__ensemble==NPTLZxy) then
          write(*,'(a,a10)') 'ensemble         :', 'NPTLZxy'
       elseif(md_condition__ensemble==NPT_PR) then
          write(*,'(a,a10)') 'ensemble         :', 'NPT_PR'
       elseif(md_condition__ensemble==OPT) then
          write(*,'(a,a10)') 'ensemble         :', 'OPT'
       endif
       write(*,'(a,f10.2,a)') 'temperature      :', systemp,' [K]'
       write(*,'(a,es10.3,a)')'pressure         :', syspres,' [Pa]'
       write(*,'(a,es10.3,a)')' _inner          :', Pinner ,' [Pa]'
       write(*,'(a,es10.3,a)')' _outer          :', Pouter ,' [Pa]'
       write(*,'(a,es10.3,a)')' _outermost      :', Poutermost,' [Pa]'
       if(reset_maxwell)then
          write(*,'(a,a10)')'maxwell_velocity :', 'yes'
       else
          write(*,'(a,a10)')'maxwell_velocity :', 'no'
       endif
       if(velocity_scaling)then
          write(*,'(a,a10)')'velocity_scaling :', 'yes'
       else
          write(*,'(a,a10)')'velocity_scaling :', 'no'
       endif
       if(velocity_scaling_region_wise)then
          write(*,'(a,a10)')'velocity_scaling_region_region_wise :', 'yes'
!      else
!         write(*,'(a,a10)')'velocity_scaling_region_region_wise :', 'no'
       endif
       if(volume_scaling)then
          write(*,'(a,a10)')'volume_scaling   :', 'yes'
       else
          write(*,'(a,a10)')'volume_scaling   :', 'no'
       endif
       if(cellgeometry_scaling)then
          write(*,'(a,a10)')'cellgeo_scaling  :', 'yes'
       else
          write(*,'(a,a10)')'cellgeo_scaling  :', 'no'
       endif
       write(*,'(a,f10.1,a)')'tau_Q (thermo)   :', tautherm*1d+12,' [ps]'
       write(*,'(a,f10.1,a)')'tau_Q (baro)     :', taubaro*1d+12,' [ps]'
       write(*,'(a,f10.1,a)')'tau_W (baro)     :', taubaroW*1d+12,' [ps]'
       !FF
       write(*,'(a)') '#'
       write(*,'(a)') '# INTERACTION'
       write(*,'(a)') '#'
       if(     md_condition__force_field == CHARMM)then
          write(*,'(a,a10)') 'force field      :', 'CHARMM'
          write(*,'(a,i10)') 'CMAP_version     :', CMAPVersion
       else if(md_condition__force_field == OPLSAA)then
          write(*,'(a,a10)') 'force field      :', 'OPLSAA'
          write(*,'(a,f10.2)') 'LJ 1-4 scaling   :', 1d0/md_oplsaa_lj_sp_divide
          write(*,'(a,f10.2)') 'CL 1-4 scaling   :', 1d0/md_oplsaa_coulomb_sp_divide
       else if(md_condition__force_field == AMBER )then
          write(*,'(a,a10)') 'force field      :', 'AMBER '
       else if(md_condition__force_field == GAmbFF)then
          write(*,'(a,a10)') 'force field      :', '  GAFF'
          write(*,'(a,f10.2)') 'LJ 1-4 scaling   :', 1d0/md_oplsaa_lj_sp_divide
          write(*,'(a,f10.2)') 'CL 1-4 scaling   :', 1d0/md_oplsaa_coulomb_sp_divide
       else if(md_condition__force_field == KREMER )then
          write(*,'(a,a10)') 'force field      :', 'KREMER'
       endif
       if(allow_bond_breaking > 0.0d0) then
          write(*,'(a,es10.3,a)') 'bond breaking at :', allow_bond_breaking * 1.0d10, " [A]"
       else
          write(*,'(a)') "bond breaking    :  disabled"
       endif
       !LJ
       write(*,'(a)')     'pairwise method  :       MTD'
       write(*,'(a,f10.1,a)') 'LJ cutoff radius :', cutrad*1d+10, ' [A]'
       write(*,'(a)')     'Table functions  :        on'
#ifdef CHARMMFSW
       write(*,'(a)')     'LJ force switch  :        on'
       write(*,'(a,f10.1,a)')'LJ switch on     :', Ron*1d+10, ' [A]'
       write(*,'(a,f10.1,a)')'LJ switch off    :', Roff*1d+10, ' [A]'
#else
       write(*,'(a)')     'LJ force switch  :       off'
#endif
       if(LJ_LRC)then
          write(*,'(a)')     'LJ correction    :        on'
       else
          write(*,'(a)')     'LJ correction    :       off'
       endif
       !Coulomb
       if (md_periodic__type == FMM) then
          write(*,'(a)')     'Coulomb method   :       FMM'
          write(*,'(a,i10)') 'nmax             :', nmax
          write(*,'(a,i10)') 'nlevel           :', nlevel
          write(*,'(a,i10)') 'ncellx           :', ncellx
          write(*,'(a,i10)') 'ncelly           :', ncelly
          write(*,'(a,i10)') 'ncellz           :', ncellz
          write(*,'(a,i10)') 'ULswitch         :', lgflg
          if(ewald_sterm)then
             write(*,'(a)')     'Ewald surf. term :        on'
          else
             write(*,'(a)')     'Ewald surf. term :       off'
          endif
       elseif (md_periodic__type == EWALD) then
          write(*,'(a)')     'Coulomb method   :     Ewald'
          write(*,'(a,es10.3,a)') 'Ewald alpha      :', ewald_alpha,' [1/m]'
          write(*,'(a,es10.1,a)') 'Ewald tolerance  :', derfc(ewald_alpha*cutrad)
          write(*,'(a,i10)') 'ncellx           :', ncellx
          write(*,'(a,i10)') 'ncelly           :', ncelly
          write(*,'(a,i10)') 'ncellz           :', ncellz
          write(*,'(a,i10)') '|hmax|**2        :', max_h2
       elseif (md_periodic__type == PMEWALD) then
          write(*,'(a)')     'Coulomb method   :       PME'
          write(*,'(a,es10.3,a)') 'Ewald alpha      :', ewald_alpha,' [1/m]'
          write(*,'(a,es10.1,a)') 'Ewald tolerance  :', derfc(ewald_alpha*cutrad)
          write(*,'(a,i10)') 'ncellx           :', ncellx
          write(*,'(a,i10)') 'ncelly           :', ncelly
          write(*,'(a,i10)') 'ncellz           :', ncellz
          write(*,'(a,i10)') 'bspline order    :', bsorder
          write(*,'(a,i10)') 'fftgridx         :', nfft1
          write(*,'(a,i10)') 'fftgridy         :', nfft2
          write(*,'(a,i10)') 'fftgridz         :', nfft3
       else
          write(0,*) 'ERROR: This modylas is FMM/PME/Ewald version.'
          write(0,*) '        Only FMM/PME/Ewald is supported.'
          call modylas_abort()
       endif
       !MPI & OpenMP
       write(*,'(a)') '#'
       write(*,'(a)') '# PARALLEL CONDITION'
       write(*,'(a)') '#'
       write(*,'(a,i10)') 'nprocs           :', nprocs
       write(*,'(a,i10)') 'nomp             :', nomp
       if(mpi_manual_division_flg)THEN
          write(*,'(a)')     'MPI division     :    manual'
       else
          write(*,'(a)')     'MPI division     :      auto'
       endif
       write(*,'(a,i10)') 'nprocx           :', npx
       write(*,'(a,i10)') 'nprocy           :', npy
       write(*,'(a,i10)') 'nprocz           :', npz
       write(*,'(a,i10)') 'lxdiv            :', lxdiv
       write(*,'(a,i10)') 'lydiv            :', lydiv
       write(*,'(a,i10)') 'lzdiv            :', lzdiv
       if (md_periodic__type == PMEWALD) then
          write(*,'(a,i10)') 'ngxdiv           :', ngxdiv
          write(*,'(a,i10)') 'ngydiv           :', ngydiv
          write(*,'(a,i10)') 'ngzdiv           :', ngzdiv
       endif
       write(*,'(a)')     'MPI comm. method :       MTD'
       !Precision
       write(*,'(a)') '#'
       write(*,'(a)') '# PRECISION SETTING'
       write(*,'(a)') '#'
#if defined(PRECISION_P2P_MIX)
       write(*,'(a)') 'P2P              : DP/SP mix'
#elif defined(PRECISION_P2P_SP)
       write(*,'(a)') 'P2P              :   full-SP'
#else
       write(*,'(a)') 'P2P              :   full-DP'
#endif
#if defined(PRECISION_M2L_MIX)
       write(*,'(a)') 'M2L              : DP/SP mix'
#elif defined(PRECISION_M2L_SP)
       write(*,'(a)') 'M2L              :   full-SP'
#else /** default **/
       write(*,'(a)') 'M2L              :   full-DP'
#endif

    ENDIF
    call mpi_barrier(mpi_comm_world, ierr)
  end subroutine output_log

end module display_log_mod

