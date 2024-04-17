!----------------------------------------------------------------------
!MODYLAS ver. 1.1.0 
!
!Copyright (c) 2014-2019 Nagoya University
!              2020-2023 The University of Tokyo
!
!Released under the MIT license.
!see https://opensource.org/licenses/MIT
!---------------------------------------------------------------------
!MODYLAS Developers:
!Yoshimichi Andoh, Kazushi Fujimoto, Tatsuya Sakashita, Noriyuki Yoshii, 
!Zhiye Tang, Jiachao Zhang, Yuta Asano, Ryo Urano, Tetsuro Nagai, 
!Atsushi Yamada, Hidekazu Kojima, Kensuke Iwahashi, Fumiyasu Mizutani, 
!Shin-ichi Ichikawa, and Susumu Okazaki.
!----------------------------------------------------------------------
!>
!! \file
!! \brief  Module and subroutines to profile whole program.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to profile run time.
!! \author Shin-ichi Ichikawa
!<
#ifdef PROF_TIMER

#ifdef __INTEL_COMPILER
#define gettod clockx
#endif
      module m_timer_simple
#include "timing.h90"
      implicit none
      integer, parameter :: n_tim = 140   ! num of timers. PLEASE MODIFY.
      integer, parameter :: iunit = 7     ! output unit number (6 means stdout)
      integer, parameter :: ign_zero = 1  ! 1: timer result is not printed if ncount[i] equals 0
      integer, parameter :: l_name = 30   ! length of name

      double precision :: timers(n_tim), tsta(n_tim)
      double precision :: ttot(3)
      integer :: ncount(n_tim), ncnt_e(n_tim)
      integer :: iregion_id, in_region
      integer :: ext_target(n_tim), n_ext_target=0
      integer :: len_name(n_tim)
      character(l_name) name(n_tim)
      integer :: prof_start=-1
      end module m_timer_simple

!-----------------------------------------------------------------------
      subroutine timer_init
      use m_timer_simple
      implicit none
      include 'mpif.h'
      real(4), external :: etime
      integer i, l, id_env, iostat
      character(1024) env
      integer n, iwk(n_tim)
      real(4) rdum, rarray(2)

      prof_start = -1
      iregion_id = 0

      ! clear counters
      do i = 1, n_tim
        timers(i) = 0.0d0
        ncount(i) = 0
        ncnt_e(i) = 0
        name(i)   = ''
      enddo

      ! set timer name. PLEASE MODIFY.
      name(TM_ALL)  = 'TMI_ALL'      ! includes all.
      name(TM_NPT_A)= 'TMI_NPT_A'    ! includes cost under npt_a. rec_curr_state and unmeasured local code around npt_a are not included.
      name(TM_NVE)  = 'TMI_NVE'      ! includes cost under npt_a. rec_curr_state and unmeasured local code around nve are not included.
      name(TM_P2M)  = 'TMI_P2M'
      name(TM_M2M)  = 'TMI_M2M'
      name(TM_M2L)  = 'TMI_M2L'
      name(TM_L2L)  = 'TMI_L2L'
      name(TM_L2P)  = 'TMI_L2P'
      name(TM_FMM_EWALD)               = 'TMI_FMM_EWALD'
      name(TM_FMM_COMM_GLOBAL)         = 'TMI_FMM_COMM_GLOBAL'
      name(TM_FMM_COMM_LOCAL)          = 'TMI_FMM_COMM_LOCAL'
      name(TM_COMM_FMM_SUPER)          = 'TMI_COMM_FMM_SUPER'
      name(TM_BG_ENE)                  = 'TMI_BG_ENE'
      name(TM_BCG_SURF)                = 'TMI_BCG_SURF'
      name(TM_RM_EWALD_SURF)           = 'TMI_RM_EWALD_SURF'
      name(TM_FLONG_POST)              = 'TMI_FLONG_POST'
      name(TM_M2L_TRANSP)              = 'TMI_M2L_TRANSP'
      name(TM_M2L_REDUCT)              = 'TMI_M2L_REDUCT'

      name(TM_ALLREDUCE_VIR)           = 'TMI_ALLREDUCE_VIR'            ! included in TM_PRE_MTLOOP, TM_MTLOOP_SHAKEITR_LOCAL, and TM_MTLOOP_RATTLEITR_LOCAL.
      name(TM_MTLOOP_SEGPREP)          = 'TMI_MTLOOP_SEGPREP'
      name(TM_MTLOOP_SHAKEITR_LOCAL)   = 'TMI_MTLOOP_SHAKEITR_LOCAL'
      name(TM_SHAKE_ROLL)              = 'TMI_SHAKE_ROLL'
      name(TM_MTLOOP_SHAKEITR_POST)    = 'TMI_MTLOOP_SHAKEITR_POST'
      name(TM_MTLOOP_LOCAL)            = 'TMI_MTLOOP_LOCAL'             ! for NVE.
      name(TM_ROLL_CONV)               = 'TMI_ROLL_CONV'
      name(TM_ALLREDUCE_ROLL)          = 'TMI_ALLREDUCE_ROLL'           ! not included in TM_SHAKE_ROLL nor in TM_RATTLE_ROLL nor in TM_ROLL_CONV.
      name(TM_COMM_BOUND)              = 'TMI_COMM_BOUND'
      name(TM_COMM_DIRECT3)            = 'TMI_COMM_DIRECT3'
      name(TM_APPL_PBC)                = 'TMI_APPL_PBC'
      name(TM_FSHORT)                  = 'TMI_FSHORT'
      name(TM_COMM_DIRECT2_FSHORT)     = 'TMI_COMM_DIRECT2_FSHORT'
      name(TM_RMVOID_ETC)              = 'TMI_RMVOID_ETC'
      name(TM_FORCES_MIDDLE_BUFINIT)   = 'TMI_FORCES_MIDDLE_BUFINIT'
      name(TM_WRAP_FORCES_LOCAL)       = 'TMI_WRAP_FORCES_LOCAL'         ! for NVE.
      name(TM_FMIDDL)                  = 'TMI_FMIDDL'                   ! for NVE.
      name(TM_COMM_DIRECT2_FSHT_MDL_LNG) = 'TMI_COMM_DIRECT2_FSHT_MDL_LNG' ! for NVE.

      name(TM_ENERGY_DIRECT)           = 'TMI_ENERGY_DIRECT'
      name(TM_FMM_NEAR)                = 'TMI_FMM_NEAR'                 ! included in TM_ENERGY_DIRECT, used when Not MTD.
      name(TM_COMM_DIRECT2_FMIDDL)     = 'TMI_COMM_DIRECT2_FMIDDL'
      name(TM_MTLOOP_SEGPREP2)         = 'TMI_MTLOOP_SEGPREP2'
      name(TM_MTLOOP_RATTLEITR_LOCAL)  = 'TMI_MTLOOP_RATTLEITR_LOCAL'
      name(TM_RATTLE_ROLL)             = 'TMI_RATTLE_ROLL'
      name(TM_RM_MOM)                  = 'TMI_RM_MOM'
      name(TM_ALLREDUCE_RM_MOM)        = 'TMI_ALLREDUCE_RM_MOM'         ! independent.
      name(TM_THERMO)                  = 'TMI_THERMO'                   ! called also from out of TM_ALL. Cared by TIME_PROFILE_START.
      name(TM_ALLREDUCE_THERMO)        = 'TMI_ALLREDUCE_THERMO'         ! included in TM_THERMO. called also from out of TM_ALL. 
      name(TM_PRE_MTLOOP)              = 'TMI_PRE_MTLOOP'
      name(TM_CHARMM_DIHEDRAL)         = 'TMI_CHARMM_DIHEDRAL'          ! included in TM_FSHORT.
      name(TM_CHARMM_ANGLE)            = 'TMI_CHARMM_ANGLE'             ! included in TM_FSHORT.
      name(TM_CHARMM_BOND)             = 'TMI_CHARMM_BOND'              ! included in TM_FSHORT.
      name(TM_ALLREDUCE_SCALERVIRIAL)  = 'TMI_ALLREDUCE_SCALERVIRIAL'   ! included in TM_MTLOOP_SHAKEITR_LOCAL and TM_MTLOOP_RATTLEITR_LOCAL.

      name(TM_PME)             = 'TMI_PME'   ! includes TM_PME_ss_ww. TM_ADD_PMEWALD, TM_ADD_PMEWALD_GATHER, TM_PMEWALD_PZFFT3DV, TM_COMM_QGRID.
      name(TM_PME_ss_ww)       = 'TMI_PME_ss_ww'
      name(TM_PME_NEAR)        = 'TMI_PME_NEAR'
      name(TM_FMM_IO)          = 'TMI_FMM_IO'
      name(TM_FLONG_PRE)       = 'TMI_FLONG_PRE'                      ! for PME flong.
      name(TM_ADD_PMEWALD)        = 'TMI_ADD_PMEWALD'                 ! for PME flong.
      name(TM_ADD_PMEWALD_GATHER) = 'TMI_ADD_PMEWALD_GATHER'          ! for PME flong
      name(TM_PMEWALD_PZFFT3DV)   = 'TMI_PMEWALD_PZFFT3DV'            ! for PME flong. includes alltoall comm.
      name(TM_COMM_QGRID)         = 'TMI_COMM_QGRID'                  ! for PME flong.
      name(TM_ADD_PMEWALD_TI0)    = 'TMI_ADD_PMEWALD_TI0'             ! for PME flong.
      name(TM_ADD_PMEWALD_GATHER_TI0) = 'TMI_ADD_PMEWALD_GATHER_TI0'  ! for PME flong
      name(TM_PMEWALD_PZFFT3DV_TI0)   = 'TMI_PMEWALD_PZFFT3DV_TI0'    ! for PME flong. includes alltoall comm.
      name(TM_COMM_QGRID_TI0)         = 'TMI_COMM_QGRID_TI0'          ! for PME flong.

      name(TM_IO_RECCURSTATE)  = 'TMI_IO_RECCURSTATE'

      name(TMB_PRE_MAINLOOP)           = 'TMBI_PRE_MAINLOOP'           ! unnecessary for performance evaluation.
      name(TMB_MDSTEP_TOP)             = 'TMBI_MDSTEP_TOP'             ! unnecessary for performance evaluation.
      name(TMB_FMM_COMM_GLOBAL)        = 'TMBI_FMM_COMM_GLOBAL'
      name(TMB_FMM_COMM_LOCAL)         = 'TMBI_FMM_COMM_LOCAL'
      name(TMB_COMM_BOUND)             = 'TMBI_COMM_BOUND'
      name(TMB_COMM_DIRECT3)           = 'TMBI_COMM_DIRECT3'
      name(TMB_COMM_DIRECT2_FSHORT)    = 'TMBI_COMM_DIRECT2_FSHORT'
      name(TMB_COMM_DIRECT2_FMIDDL)    = 'TMBI_COMM_DIRECT2_FMIDDL'
      name(TMB_COMM_DIRECT2_FSHT_MDL_LNG) = 'TMBI_COMM_DIRECT2_FSHT_MDL_LNG'  ! for PME.
      name(TMB_COMM_FMM_SUPER)         = 'TMBI_COMM_FMM_SUPER'
      name(TMB_ALLREDUCE_VIR)          = 'TMBI_ALLREDUCE_VIR'         ! included in TM_PRE_MTLOOP, TM_MTLOOP_RATTLEITR_LOCAL, TM_MTLOOP_SHAKEITR_LOCAL, and TM_MTLOOP_RATTLEITR_LOCAL.
      name(TMB_ALLREDUCE_ROLL)         = 'TMBI_ALLREDUCE_ROLL'        ! not included in TM_SHAKE_ROLL nor in TM_RATTLE_ROLL nor in TM_ROLL_CONV.
      name(TMB_ALLREDUCE_RM_MOM)       = 'TMBI_ALLREDUCE_RM_MOM'      ! independent.
      name(TMB_ALLREDUCE_THERMO)       = 'TMBI_ALLREDUCE_THERMO'      ! included in TM_THERMO. called also from out of TM_ALL. 
      name(TMB_FMM_IO)                 = 'TMBI_FMM_IO'
      name(TMB_IO_RECCURSTATE)         = 'TMBI_IO_RECCURSTATE'        ! independent.
      name(TMB_ALLREDUCE_SCALERVIRIAL) = 'TMBI_ALLREDUCE_SCALERVIRIAL'   ! included in TM_MTLOOP_SHAKEITR_LOCAL and TM_MTLOOP_RATTLEITR_LOCAL.
      name(TMB_ADD_PMEWALD_GATHER)     = 'TMBI_ADD_PMEWALD_GATHER'     ! for PME flong.
      name(TMB_PMEWALD_PZFFT3DV)       = 'TMBI_PMEWALD_PZFFT3DV'       ! for PME flong.
      name(TMB_COMM_QGRID)             = 'TMBI_COMM_QGRID'             ! for PME flong.

      name(TM_EVB_MD2EVB)              = 'TMI_EVB_MD2EVB'
      name(TM_COMEVB_PBC_ATOM)         = 'TMI_COMEVB_PBC_ATOM'
      name(TM_COMEVB_PBC_SEGMENT)      = 'TMI_COMEVB_PBC_SEGMENT'
      name(TM_EVB_COMPOSIT_LIST)       = 'TMI_EVB_COMPOSIT_LIST'
      name(TM_EVB_IO_CMPSITLIST)       = 'TMI_EVB_IO_CMPSITLIST'      /* included in TM_EVB_COMPOSIT_LIST. */
      name(TM_COMEVB_MERGE_EVBGRP)     = 'TMI_COMEVB_MERGE_EVBGRP'
      name(TM_COMEVB_PBC_EVBGRP)       = 'TMI_COMEVB_PBC_EVBGRP'
      name(TM_EVB_ENERGY)              = 'TMI_EVB_ENERGY'
      name(TM_EVB_EIGENV)              = 'TMI_EVB_EIGENV'
      name(TM_COMEVB_PBC_EIGENV)       = 'TMI_COMEVB_PBC_EIGENV'
      name(TM_EVB_SEGPURGE)            = 'TMI_EVB_SEGPURGE'
      name(TM_COMEVB_MERGE_SEGPURGE)   = 'TMI_COMEVB_MERGE_SEGPURGE'
      name(TM_EVB_UPD_WSEGC)           = 'TMI_EVB_UPD_WSEGC'
      name(TM_COMEVB_SEGSORT)          = 'TMI_COMEVB_SEGSORT'
      name(TM_EVB2MD)                  = 'TMI_EVB2MD'

      name(TMB_COMEVB_PBC_ATOM)        = 'TMBI_COMEVB_PBC_ATOM'
      name(TMB_COMEVB_PBC_SEGMENT)     = 'TMBI_COMEVB_PBC_SEGMENT'
      name(TMB_EVB_IO_CMPSITLIST)      = 'TMBI_EVB_IO_CMPSITLIST'     /* included in TM_EVB_COMPOSIT_LIST. barrier after IO. */
      name(TMB_COMEVB_MERGE_EVBGRP)    = 'TMBI_COMEVB_MERGE_EVBGRP'
      name(TMB_COMEVB_PBC_EVBGRP)      = 'TMBI_COMEVB_PBC_EVBGRP'
      name(TMB_COMEVB_PBC_EIGENV)      = 'TMBI_COMEVB_PBC_EIGENV'
      name(TMB_COMEVB_MERGE_SEGPURGE)  = 'TMBI_COMEVB_MERGE_SEGPURGE'
      name(TMB_COMEVB_SEGSORT)         = 'TMBI_COMEVB_SEGSORT'

      name(TM_M2L_LEV0)                = 'TMI_M2L_LEV0'
      name(TM_M2L_LEV1)                = 'TMI_M2L_LEV1'
      name(TM_M2L_LEV2)                = 'TMI_M2L_LEV2'
      name(TM_M2L_LEV3)                = 'TMI_M2L_LEV3'
      name(TM_M2L_LEV4)                = 'TMI_M2L_LEV4'
      name(TM_M2L_LEV5)                = 'TMI_M2L_LEV5'
      name(TM_M2L_LEV6)                = 'TMI_M2L_LEV6'
      name(TM_M2L_LEV7)                = 'TMI_M2L_LEV7'
      name(TM_M2L_LEV8)                = 'TMI_M2L_LEV8'

      name(TM_FFT_F1)                  = 'TMI_FFT_F1'
      name(TM_ALLTOALL_F1)             = 'TMI_ALLTOALL_F1'
      name(TM_FFTFCOPY_F1)             = 'TMI_FFTFCOPY_F1'
      name(TM_FFT_F2)                  = 'TMI_FFT_F2'
      name(TM_ALLTOALL_F2)             = 'TMI_ALLTOALL_F2'
      name(TM_FFT_F3)                  = 'TMI_FFT_F3'
      name(TM_FFT_B1)                  = 'TMI_FFT_B1'
      name(TM_ALLTOALL_B1)             = 'TMI_ALLTOALL_B1'
      name(TM_FFT_B2)                  = 'TMI_FFT_B2'
      name(TM_FFTFCOPY_B1)             = 'TMI_FFTFCOPY_B1'
      name(TM_ALLTOALL_B2)             = 'TMI_ALLTOALL_B2'
      name(TM_FFT_B3)                  = 'TMI_FFT_B3'
      name(TMB_ALLTOALL_F1)            = 'TMBI_ALLTOALL_F1'
      name(TMB_ALLTOALL_F2)            = 'TMBI_ALLTOALL_F2'
      name(TMB_ALLTOALL_B1)            = 'TMBI_ALLTOALL_B1'
      name(TMB_ALLTOALL_B2)            = 'TMBI_ALLTOALL_B2'

      name(TM_PME_COEFF)               = 'TMI_PME_COEFF'              /* included in TM_ADD_PMEWALD. */
      name(TM_PME_BSPFILL)             = 'TMI_PME_BSPFILL'            /* included in TM_ADD_PMEWALD. */
      name(TM_PME_CALC_QGRID)          = 'TMI_PME_CALC_QGRID'         /* included in TM_ADD_PMEWALD. */
      name(TM_PME_ENERGY)              = 'TMI_PME_ENERGY'             /* included in TM_ADD_PMEWALD. */
      name(TM_PME_FORCE)               = 'TMI_PME_FORCE'              /* included in TM_ADD_PMEWALD. */


!unused      name(3)  = 'TM_MAIN_LOOP'
!force retire      name(3)  = 'TM_FMM'
!unused      name(3)  = 'TM_SHAKE'
!unused      name(4)  = 'TM_RATTLE'
!unused      name(5)  = 'TM_COMM_DIRECT'
!unused      name(6)  = 'TM_MIGRATION'
!unused      name(7)  = 'TM_ENE_REDUCTION'
!unused      name(8)  = 'TM_OUTPUT'
!unused      name(3)  = 'TM_COMM_FMM'

      ! keep len_trim(name(i)) for "!PA"
      do i = 1, n_tim
        len_name(i) = len_trim(name(i))
      enddo

      ! real, user, sys of total
      rdum = etime(rarray)
      ttot(1) = MPI_WTIME()  ! real
      ttot(2) = rarray(1)    ! user
      ttot(3) = rarray(2)    ! sys

      return
      end

!-----------------------------------------------------------------------
      subroutine timer_fin
      use m_timer_simple
      implicit none
      include 'mpif.h'
      real(4), external :: etime
      character(MPI_MAX_PROCESSOR_NAME) procname
      integer i, n_err
      integer m, n_err_g
      integer myrank, nprocs, ierr, resultlen
      real(4) rdum, rarray(2)
      double precision, allocatable :: ttot_g(:,:), timers_g(:,:)
      integer, allocatable :: ncount_g(:,:)
      character(MPI_MAX_PROCESSOR_NAME), allocatable :: procname_g(:)

      ! convert usec to sec
      do i = 1, n_tim
        timers(i) = timers(i) * 1.0d-6
      enddo
      ! user and system cpu time
      rdum = etime(rarray)
      ttot(1) = MPI_WTIME() - ttot(1)  ! real
      ttot(2) = rarray(1)   - ttot(2)  ! user
      ttot(3) = rarray(2)   - ttot(3)  ! sys

      ! check count of timer call
      n_err = 0
      do i = 1, n_tim
        if (ncount(i) .ne. ncnt_e(i)) then
          write(*,*) 'Error: timer_fin: not pair timer_(sta|end), ',    &
     &         'ID = ', i
          n_err = n_err + 1
        endif
      enddo
      call MPI_ALLREDUCE(n_err, n_err_g, 1, MPI_INTEGER, MPI_MAX,       &
     &                   MPI_COMM_WORLD, ierr)
      n_err = n_err_g
      if (n_err .gt. 0) then
        return
      endif

      call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
      call MPI_GET_PROCESSOR_NAME(procname, resultlen, ierr)

      if (myrank .eq. 0) then
        allocate(ttot_g(3,0:nprocs-1))
        allocate(procname_g(0:nprocs-1))
        allocate(ncount_g(n_tim,0:nprocs-1))
        allocate(timers_g(n_tim,0:nprocs-1))
      else
        allocate(ttot_g(3,0:0))
        allocate(procname_g(0:0))
        allocate(ncount_g(n_tim,0:0))
        allocate(timers_g(n_tim,0:0))
      endif

      call MPI_GATHER(ttot,   3, MPI_DOUBLE_PRECISION,                  &
     &                ttot_g, 3, MPI_DOUBLE_PRECISION,                  &
     &                0, MPI_COMM_WORLD, ierr)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_GATHER(procname,   MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER,&
     &                procname_g, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER,&
     &                0, MPI_COMM_WORLD, ierr)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_GATHER(ncount,   n_tim, MPI_INTEGER,                     &
     &                ncount_g, n_tim, MPI_INTEGER,                     &
     &                0, MPI_COMM_WORLD, ierr)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_GATHER(timers,   n_tim, MPI_DOUBLE_PRECISION,            &
     &                timers_g, n_tim, MPI_DOUBLE_PRECISION,            &
     &                0, MPI_COMM_WORLD, ierr)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      if (myrank .eq. 0) then
        if (iunit .ne. 6) open(iunit, file='timer.out')

        do m = 0, nprocs - 1
          ! print header
          write(iunit,'("cfj # rank = ",i6,                             &
     &                 2x,"real,user,sys = ",3f13.3,                    &
     &                 2x,"procname = ",a)')                            &
     &         m, ttot_g(1:3,m), procname_g(m)(:len_trim(procname_g(m)))
          if (iregion_id .ne. 0) then
            write(iunit,'("cfj # FJ_TIMER_REGION = ",i3)') iregion_id
          endif
          if (n_ext_target .ne. 0) then
            write(iunit,'("cfj # FJ_TIMER_EXT_TARGETS =")',advance='no')
            do i = 1, n_tim
              if (ext_target(i) .eq. 1) then
                write(iunit,'(1x,i0)',advance='no') i
              endif
            enddo
            write(iunit,'()')
          endif

          ! print timers
          do i = 1, n_tim
            if (ign_zero .ne. 1 .or. ncount_g(i,m) .ne. 0) then
              write(iunit,'("cfj ",i3,1x,a,1x,f16.6,1x,i8)')            &
     &             i, name(i), timers_g(i,m), ncount_g(i,m)
            endif
          enddo
          write(iunit,'()')
        enddo

        if (iunit .ne. 6) close(iunit)
      endif

      deallocate(ttot_g)
      deallocate(procname_g)
      deallocate(ncount_g)
      deallocate(timers_g)

      return
      end

!-----------------------------------------------------------------------
      subroutine timer_bar(id)
      use m_timer_simple
      implicit none
      include 'mpif.h'
      integer id
      integer ierr
      double precision tend

      if(prof_start > 0) then
!$omp master
        call gettod(tsta(id))
        ncount(id) = ncount(id) + 1
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call gettod(tend)
        timers(id) = timers(id) + (tend - tsta(id))
        ncnt_e(id) = ncnt_e(id) + 1
!$omp end master
      endif

      return
      end

! -------------------------------------------
      subroutine timer_bar_com(id,icomm)
      use m_timer_simple
      implicit none
      include 'mpif.h'
      integer id
      integer icomm
      integer ierr
      double precision tend

      if(prof_start > 0) then
!$omp master
        call gettod(tsta(id))
        ncount(id) = ncount(id) + 1
        call MPI_BARRIER(icomm, ierr)
        call gettod(tend)
        timers(id) = timers(id) + (tend - tsta(id))
        ncnt_e(id) = ncnt_e(id) + 1
!$omp end master
      endif

      return
      end

!-----------------------------------------------------------------------
      subroutine timer_sta(id)
      use m_timer_simple
      implicit none
      integer id

      if(prof_start > 0) then
#ifdef FJ_HWPROF
!!      call start_collection(name(id))
        call fapp_start(trim(name(id)), 0, 0)
#endif

!$omp master
        call gettod(tsta(id))
        ncount(id) = ncount(id) + 1
!$omp end master
      endif

      return
      end

!-----------------------------------------------------------------------
      subroutine timer_end(id)
      use m_timer_simple
      implicit none
      integer id
      double precision tend

      if(prof_start > 0) then
!$omp master
        call gettod(tend)
        timers(id) = timers(id) + (tend - tsta(id))
        ncnt_e(id) = ncnt_e(id) + 1
!$omp end master

#ifdef FJ_HWPROF
!!      call stop_collection(name(id))
        call fapp_stop(trim(name(id)), 0, 0)
#endif
      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine timer_profsta()
      use m_timer_simple
      implicit none

!$omp master
      prof_start = 1
!$omp end master

      return
      end
!-----------------------------------------------------------------------
      subroutine timer_profstop()
      use m_timer_simple
      implicit none

!$omp master
      prof_start = -1
!$omp end master

      return
      end
!-----------------------------------------------------------------------

#endif /* PROF_TIMER */
