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
!! \brief  Module and subroutines to control file output.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to control file output.
!! \author Kensuke Iwahashi
!<
module file_application
  implicit none
contains
!-----------------------------------------------------------------------
!> 
!! \brief  Subroutine to write all of various outputs.
!! \author Kensuke Iwahashi
!< 
  subroutine record_current_state
    use md_condition
    use mpi_tool
    use comm_bound_mod, only : pre_record_data
    use file_mdrun
    use file_mdmntr
    use file_mdtrj
    use file_restart
    use file_force
    use file_dcd
#ifdef DIVIDEDCD
    use file_dcd_divided
#endif
    use file_xtc
    use file_utility, only : ascii_output
#ifdef GROEXT
    use file_groext
#endif
    implicit none
    logical :: flg_pre_record_data=.false. ! default

!### record mdrun & mdmntr ###
    IF (myrank.eq.mpiout) THEN
       call record_current_monitors
       call record_current_runtime_info
    ENDIF

#ifdef DIVIDEDCD
    IF(trj_output)THEN
       if(mod((mdstep-trj_start),trj_interval)==0) then
          flg_pre_record_data=.true.
       endif
    ENDIF
    IF(rst_output)THEN
       if (mod((mdstep-restart_start),restart_interval)==0) then
          flg_pre_record_data=.true.
       endif
    ENDIF
#else
    IF(    dcd_output)THEN
       if( mod((mdstep-dcd_start),dcd_interval)==0) then
          flg_pre_record_data=.true.
       endif
    ENDIF
    IF(trj_output)THEN
       if(mod((mdstep-trj_start),trj_interval)==0) then
          flg_pre_record_data=.true.
       endif
    ENDIF
    IF(rst_output)THEN
       if (mod((mdstep-restart_start),restart_interval)==0) then
          flg_pre_record_data=.true.
       endif
    ENDIF
#ifdef XTC
    IF(xtc_output)THEN
       if (mod((mdstep-restart_start),restart_interval)==0) then
          flg_pre_record_data=.true.
       endif
    ENDIF
#endif        
#endif

    if(flg_pre_record_data) call pre_record_data

    !### record mdtrj ###
    if (mod((mdstep-trj_start),trj_interval)==0) then
       IF (myrank.eq.mpiout) THEN
          if(trj_output) call record_current_trajectory
       ENDIF
    endif

#ifdef DIVIDEDCD
    if (mod((mdstep-dcd_start),dcd_interval)==0) then
       IF (myrank.eq.mpiout) THEN
          if(dcd_output) call record_current_trajectory_dcd_info
       ENDIF
       if(dcd_output) call record_current_trajectory_dcd_xyz
    endif
#else
    !### record .dcd ###
    if (mod((mdstep-dcd_start),dcd_interval)==0) then
       IF (myrank.eq.mpiout) THEN
          if(dcd_output) call record_current_trajectory_dcd
       ENDIF
    endif
#endif

#ifdef XTC
    !### record .xtc ###
    if (mod((mdstep-xtc_start),xtc_interval)==0) then
       IF (myrank.eq.mpiout) THEN
          if(xtc_output) call record_current_trajectory_xtc
       ENDIF
    endif
#endif

    !### record restart.asc ###
    if (mod((mdstep-restart_start),restart_interval)==0)then
       IF (myrank.eq.mpiout) THEN
          if(ascii_output) call record_restart_ascii
       ENDIF
    endif

    !### record restart.bin ###
    if (mod((mdstep-restart_start),restart_interval)==0)then
       IF (myrank.eq.mpiout) THEN
          if(rst_output) call record_restart_binary
       ENDIF
    endif

#ifdef GROEXT
    if (mod((mdstep-restart_start),restart_interval)==0)then
       IF (myrank.eq.mpiout) THEN
          if(gmx_output_groext) then
             if( gmx_write_gro ) then
                call update_gro_positions
                call record_gro
             end if
          end if
       ENDIF
    endif
#endif

#ifdef SEGSHAKE
    !### record calculated mean force ###
    if (mod((mdstep-force_start),force_interval)==0) then
       IF (myrank.eq.mpiout) THEN
          call record_current_force
       ENDIF
    endif
#endif

  flg_pre_record_data=.false. ! Reset value for next cycle. 

  end subroutine record_current_state
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to open output files and write their header.
!! \author Kensuke Iwahashi
!<
  subroutine init_g_main
    use mpi_tool
    use file_mdrun
    use file_mdmntr
    use file_mdtrj
    use file_restart
    use file_force
    use file_dcd
    use file_mdbb
    use md_condition, only : allow_bond_breaking
#ifdef DIVIDEDCD
    use file_dcd_divided
#endif
    use file_xtc
    implicit none

#ifdef DIVIDEDCD
    call open_dcd_xyz
#endif

    if (myrank.eq.mpiout) then
       !       open *.mdtrj file to which trajectory of MD calculation
       !       will be written
       call open_trj
       call write_trajectory_header
       !
#ifdef DIVIDEDCD
       !       open *.dcd.header file to which cell info of MD calculation
       !       will be written in dcd format
       call open_dcd_info
       call write_dcd_info_header
#else
       !       open *.dcd file to which trajectory of MD calculation
       !       will be written in dcd format
       call open_dcd
       call write_dcd_header
#ifdef XTC
       !       open *.xtc file to which trajectory of MD calculation
       !       will be written in dcd format
       call open_xtc
#endif /**XTC**/
#endif /**DIVIDEDCD**/
       !
       !       open *.mdmntr file to which monitor variables (Hamiltonian etc.)
       !       will be written
       call open_mntr
       call write_monitors_header
       !
       !       open *.mdrun file to which runtime information will be written
       call open_run
       !
#ifdef SEGSHAKE
       !       open *.force file to which mean-force will be written
       call open_force
       call write_force_header
#endif
       !       open *.restart.bin file to which force log will be written
       call open_restart
       !

if(allow_bond_breaking>0.0d0) then
   call open_mdbb
endif

    endif
  end subroutine init_g_main
!-----------------------------------------------------------------------
!>
!! \brief  Subroutine to close all output files.
!! \author Kensuke Iwahashi
!<
  subroutine cleanup
    use device_numbers
    use md_condition
    use file_restart
    use file_dcd
    use file_xtc
    use file_mdtrj

    close(f_rvf)
    if(trj_output) close(f_trj)
    close(f_run)
    close(f_mntr)
    if(dcd_output) close(f_dcd)

#ifdef XTC
    if(xtc_output) then
       ! call xtcf % close
       call xtc_out % close
    endif
#endif    

#ifdef SEGSHAKE
    close(f_force)
#endif
    if(rst_output) close(f_restart_bin)

    if(allow_bond_breaking>0.0d0) then
      close(f_mdbb)
    endif

  end subroutine cleanup

end module file_application
