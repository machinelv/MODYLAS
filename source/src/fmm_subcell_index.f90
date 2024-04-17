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
!! \brief Module and subroutines to check relationship between subcell/nprocs.
!<
!----------------------------------------------------------------------
!>
!! \brief Module and subroutines to check relationship between subcell/nprocs.
!! \author Tatsuya Sakashita, Shin-ichi Ichikawa
!<
module fmm_subcell_index
  implicit none

contains

  ! NOT same as lzdiv
  integer pure elemental function get_cell_size(ncell, mcell_size, nproc)
    implicit none
    integer, intent(in) :: ncell, mcell_size
    integer, intent(in) :: nproc

    get_cell_size = (ncell / mcell_size - 1) / nproc + 1
  end function get_cell_size

  ! Assume npow = 2 or 3
  integer pure elemental function get_nbound_minus(nsupercell, np, ip)
    implicit none
    integer, intent(in) :: nsupercell
    integer, intent(in) :: np, ip
    integer :: npow       ! number of power element (2 or 3)
    integer :: icg_start
    
    npow=2
    if(mod(nsupercell,2) /= 0) npow=3
    icg_start = (nsupercell * ip) / np  ! 0-started
    get_nbound_minus = 2 * npow + mod(icg_start, npow)
  end function get_nbound_minus

  integer pure elemental function get_nbound_plus(nsupercell, np, ip)
    implicit none
    integer, intent(in) :: nsupercell
    integer, intent(in) :: np, ip
    integer :: npow       ! number of power element (2 or 3)
    integer :: icg_end

    npow=2
    if(mod(nsupercell,2) /= 0) npow=3
    icg_end = (nsupercell * (ip+1) - 1) / np  ! 0-started
    get_nbound_plus = 3 * npow - 1 - mod(icg_end, npow)
  end function get_nbound_plus

  integer pure elemental function get_ncell_with_bounds(ncell, mcell_size, np, ip)
    implicit none
    integer, intent(in) :: ncell, mcell_size
    integer, intent(in) :: np, ip
    integer :: nsupercell

    nsupercell = (ncell - 1) / mcell_size + 1
  
    get_ncell_with_bounds = get_cell_size(ncell, mcell_size, np) &
         & + get_nbound_minus(nsupercell, np, ip) &
         & + get_nbound_plus(nsupercell, np, ip)

  end function get_ncell_with_bounds

  !-----------------------------------------------------------------------
  !>
  !! \brief Factorize the number of cells to power of 2 or 3
  !! \author Tatsuya Sakashita
  !! \param[in] ncell the number of cells in some axis
  !! \param[in] nlevel total number of FMM levels
  !! \param[out] array_nl
  !! \param[in] axis character which expresses axis to display error message
  !! \note
  !! If the following conditions are not satisfied, display error message, and abort
  !! \li \p ncell should not inculde other factors except 2 or 3.
  !! \li For factorization \p ncell \f$=2^m 3^n\f$, m + n should be equal to \p nlevel.
  !! The argument \p axis is needed, because this subroutine itself cannot distingush axis.
  !! \todo
  !! Issuing abort in this subroutine may not be a good idea.
  !! However, Fortran grammer does not have "exception".
  subroutine factor_ncell(ncell, nlevel, array_nl, axis)
    use mpi_tool
    implicit none
    integer, intent(in) :: ncell, nlevel
    integer, intent(out) :: array_nl(0:nlevel)
    character, intent(in) :: axis  ! "x", "y", "z" for displaying error
    integer :: dividend
    integer :: index, count_two, count_three

    count_two = 0;  count_three = 0
    dividend = ncell
    index = 0
    do while(mod(dividend,2).eq.0)
       array_nl(index)=2
       dividend=dividend/2
       index = index + 1
       count_two = count_two + 1
    enddo
    do while(mod(dividend,3).eq.0)
       array_nl(index)=3
       dividend=dividend/3
       index = index + 1
       count_three = count_three + 1
    enddo

    if (dividend /= 1) then
       if(myrank==0) then
          write(0,'(3a,i0,4a,i0,a)') "ERROR: Number of sub-cells ncell", axis, "=", ncell, &
               & " should have only prime factors 2 or 3.", &
               & " The specified ncell", axis, " has extra factor ", dividend, "."
       endif
       call modylas_abort()
    endif

    if (index /= nlevel) then
       if (myrank==0) then
          write(0,'(6a,i0,a,i0,a,i0)') "ERROR: The sum of exponents of power of two and three in factors of ncell", axis, &
               & " should be greater than or equal to nlevel.", &
               & " ncell", axis, "=2^", count_two, "*3^", count_three, &
               & "; nlevel=", nlevel
       endif
       call modylas_abort()
    endif
    
    array_nl(nlevel)=1  ! for FMM Ewald
  end subroutine factor_ncell

  !----------------------------------------------------------------------
  !###
  !### check subcell and supercell structure ###!
  !### by Ichikawa @ 2012/12/26
  !###
  !----------------------------------------------------------------------
  subroutine check_supercell_versus_nprocs(nscellx, nscelly, nscellz)
    use mpi_3d_grid, only : npx, npy, npz
    use mpi_tool
    implicit none
    integer, intent(in) :: nscellx, nscelly, nscellz

    if (myrank==0) then
       if(nscellx > npx) then
          if(mod(nscellx,npx) /= 0) then
             write(0,'(a,i0,a,i0,a)') "ERROR: number of super-cells nscellx=", nscellx, &
                  & " cannot be divided evenly by number of processes npx=", npx, "."
             call modylas_abort()
          endif
       else if(nscellx < npx) then
          if(mod(npx,nscellx) /= 0) then
             write(0,'(a,i0,a,i0,a)') "ERROR: number of processes npx=", npx, &
                  & " cannot be divided evenly by cells number of super-cells nscellx=", nscellx, "."
             call modylas_abort()
          endif
       endif
       
       if(nscelly > npy) then
          if(mod(nscelly,npy) /= 0) then
             write(0,'(a,i0,a,i0,a)') "ERROR: number of super-cells nscelly=", nscelly, &
                  & " cannot be divided evenly by number of processes npy=", npy, "."
             call modylas_abort()
          endif
       else if(nscelly < npy) then
          if(mod(npy,nscelly) /= 0) then
             write(0,'(a,i0,a,i0,a)') "ERROR: number of processes npy=", npy, &
                  & " cannot be divided evenly by cells number of super-cells nscelly=", nscelly, "."
             call modylas_abort()
          endif
       endif
       
       if(nscellz > npz) then
          if(mod(nscellz,npz) /= 0) then
             write(0,'(a,i0,a,i0,a)') "ERROR: number of super-cells nscellz=", nscellz, &
                  & " cannot be divided evenly by number of processes npz=", npz, "."
             call modylas_abort()
          endif
       else if(nscellz < npz) then
          if(mod(npz,nscellz) /= 0) then
             write(0,'(a,i0,a,i0,a)') "ERROR: number of processes npz=", npz, &
                  & " cannot be divided evenly by cells number of super-cells nscellz=", nscellz, "."
             call modylas_abort()
          endif
       endif
    endif
  end subroutine check_supercell_versus_nprocs
  
end module fmm_subcell_index
