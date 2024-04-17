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
!! \brief  Module and subroutines to set scaling parameter of 1-4 interaction.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to set 1-4 scaling parameter.
!! \author Kensuke Iwahashi, Yoshimichi Andoh
!<
module md_oplsaa_special_divide
  implicit none
  integer(4) :: nclsp=0  !! must be zero as default (note by YA)
  real(8) :: md_oplsaa_lj_sp_divide=2.0
  real(8) :: md_oplsaa_coulomb_sp_divide=2.0
  logical :: md_oplsaa_lj_sp_divide_is_set=.false.
  logical :: md_oplsaa_coulomb_sp_divide_is_set=.false.
contains

  subroutine fmod_md_oplsaa_lj_sp_divide(value)
    implicit none
    real(8), intent(in) :: value
    md_oplsaa_lj_sp_divide = value
    md_oplsaa_lj_sp_divide_is_set = .true.
  end subroutine fmod_md_oplsaa_lj_sp_divide

  subroutine fmod_md_oplsaa_cl_sp_divide(value)
    implicit none
    real(8), intent(in) :: value
    md_oplsaa_coulomb_sp_divide = value
    md_oplsaa_coulomb_sp_divide_is_set = .true.
  end subroutine fmod_md_oplsaa_cl_sp_divide

  subroutine fmod_md_oplsaa_lj_sp_divide_mddef(value)
    use force_field_numbers
    use md_condition, only : md_condition__force_field
    use mpi_tool, only : myrank
    implicit none
  
    real(8), intent(in) :: value

    if( md_condition__force_field == CHARMM .or. &
    &   md_condition__force_field == KREMER ) then
        if(value >= 0.0d0 .and. value /= 1.0d0) then
          if(myrank==0) write(0,*) "Warning: CHARMM/KREMER aaserts that LJ&
          & 1-4 scaling factor is 1.0, override special_divide_lj in mddef by 1.0"
          !flush(0)
        endif
    else if( md_condition__force_field == OPLSAA ) then
        if(value >= 0.0d0 .and. value /= 0.5d0) then
          if(myrank==0) write(0,*) "Warning: OPLSAA aaserts that LJ 1-4 scaling&
          & factor is 0.5, override special_divide_lj in mddef by 0.5"
          !flush(0)
        endif
!    else if( md_condition__force_field == PCFF ) then
!        if(value >= 0.0d0 .and. value /= 1.0d0) then
!          if(myrank==0) write(0,*) "Warning: PCFF aaserts that LJ 1-4 scaling&
!          & factor is 1.0, override special_divide_lj in mddef by 1.0"
!          !flush(0)
!        endif
    else if( md_condition__force_field == AMBER .or. &
    &        md_condition__force_field == GAmbFF ) then
        if( md_oplsaa_lj_sp_divide_is_set ) then
            if( md_oplsaa_lj_sp_divide /= value ) then
                if(myrank==0) write(0,*) "Warning: special_divide_lj in mddef&
                & is inconsistent with the value in mdff, use value in mdff ", md_oplsaa_lj_sp_divide
                !flush(0)
            endif
        else
          md_oplsaa_lj_sp_divide = value
          md_oplsaa_lj_sp_divide_is_set = .true.
          if(myrank==0) write(0,*) "md_oplsaa_lj_sp_divide = ", md_oplsaa_lj_sp_divide, " is read from mddef for GAFF"
        endif
    endif
    !flush(0)

  end subroutine fmod_md_oplsaa_lj_sp_divide_mddef

  subroutine fmod_md_oplsaa_coulomb_sp_divide_mddef(value)
    use force_field_numbers
    use md_condition, only : md_condition__force_field
    use mpi_tool, only : myrank
    implicit none
    real(8), intent(in) :: value

    if( md_condition__force_field == CHARMM .or. &
    &   md_condition__force_field == KREMER ) then
        if(value >= 0.0d0 .and. value /= 1.0d0) then
          if(myrank==0) write(0,*) "Warning: CHARMM/KREMER aaserts that Coulomb&
          & 1-4 scaling factor is 1.0, override special_divide_coulomb in mddef by 1.0"
          !flush(0)
        endif
    else if( md_condition__force_field == OPLSAA ) then
        if(value >= 0.0d0 .and. value /= 0.5d0) then
          if(myrank==0) write(0,*) "Warning: OPLSAA aaserts that Coulomb 1-4&
          & scaling factor is 0.5, override special_divide_coulomb in mddef by 0.5"
          !flush(0)
        endif
!    else if( md_condition__force_field == PCFF ) then
!        if(value >= 0.0d0 .and. value /= 1.0d0) then
!          if(myrank==0) write(0,*) "Warning: PCFF aaserts that Coulomb 1-4&
!          & scaling factor is 1.0, override special_divide_coulomb in mddef by 1.0"
!          !flush(0)
!        endif
    else if( md_condition__force_field == AMBER .or. &
    &        md_condition__force_field == GAmbFF ) then
        if( md_oplsaa_coulomb_sp_divide_is_set ) then
            if( md_oplsaa_coulomb_sp_divide /= value ) then
                if(myrank==0) write(0,*) "Warning: special_divide_coulomb in mddef&
                & is inconsistent with the value in mdff, use value in mdff ", md_oplsaa_coulomb_sp_divide
                !flush(0)
            endif
        else
          md_oplsaa_coulomb_sp_divide = value
          md_oplsaa_coulomb_sp_divide_is_set = .true.
          if(myrank==0) write(0,*) "md_oplsaa_coulomb_sp_divide = ", md_oplsaa_coulomb_sp_divide, " is read from mddef for GAFF"
        endif
    endif
    !flush(0)

  end subroutine fmod_md_oplsaa_coulomb_sp_divide_mddef

  subroutine fmod_md_check_scaling_factor_is_set()
    use force_field_numbers
    use md_condition, only : md_condition__force_field
    use mpi_tool, only : modylas_abort
    use mpi_tool, only : myrank
    implicit none

    if( md_condition__force_field == CHARMM .or. &
    &   md_condition__force_field == KREMER ) then
        ! do nothing
    else if( md_condition__force_field == OPLSAA ) then
        ! do nothing
!    else if( md_condition__force_field == PCFF ) then
!        ! do nothing
    else if( md_condition__force_field == AMBER .or. &
    &        md_condition__force_field == GAmbFF ) then
        if( .not. md_oplsaa_lj_sp_divide_is_set ) then
            if(myrank==0) then
                write(0,*) "Error: special_divide_lj must be set either in mdff"
                write(0,*) "       or (/scaling14/special_divide_lj) in mddef for GAFF"
                !flush(0)
                call modylas_abort()
            endif
        endif
        if( .not. md_oplsaa_coulomb_sp_divide_is_set ) then
            if(myrank==0) then
                write(0,*) "Error: special_divide_coulomb must be set either in mdff "
                write(0,*) "       or in mddef (/scaling14/special_divide_coulomb) for GAFF"
                !flush(0)
                call modylas_abort()
            endif
        endif
    endif
  
  end subroutine fmod_md_check_scaling_factor_is_set

end module md_oplsaa_special_divide
