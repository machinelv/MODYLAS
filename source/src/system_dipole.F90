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
!! \brief  Module and subroutines to calculate dipole moment of system.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to calculate dipole moment of system.
!! \author Ryo Urano
!<
module system_dipole
use mpi_tool, only: myrank, mpiend
implicit none
real(8) :: dM(3), wk_dM(3)
real(8) :: M(3),  wk_M(3)
logical :: bSystemDipole = .false.
integer(4) :: file_number = 1757

contains
  subroutine set_up_system_dipole(i)
    use session_name_mod, only: session_name
    implicit none
    integer(4) :: i, io
    if(i == 1) bSystemDipole = .true.

    if(bSystemDipole) then
      if(myrank==0) then
        open(file_number, file=trim(session_name)//'_dipole.dat', iostat=io, &
             status='replace', access='sequential', form='formatted')

        write(file_number,'(A)') "#################################################################################"
        write(file_number,'(A)') "# Calculated system dipole of "//trim(session_name)
        write(file_number,'(A)') "# M(t) = sum_i^N qi*ri(t)"
        write(file_number,'(A)') "# dM(t)/dt = sum_i^N qi*vi(t)"
        write(file_number,'(A)') "# M(t): System dipole"
        write(file_number,'(A)') "# qi: charge of i-th atom"
        write(file_number,'(A)') "# ri: coordinate of i-th atom"
        write(file_number,'(A)') "# vi: velocitiy of i-th atom"
        write(file_number,'(A)') "# N: total number of atoms in the system"
        write(file_number,'(A)') "#################################################################################"
        write(file_number, '(A10, 6A20)') "#    istep",  "M(x)", "M(y)", "M(z)", "dM/dt(x)",  "dM/dt(y)", "dM/dt(z)"
        call flush(file_number)
      endif
    endif

  end subroutine set_up_system_dipole

  subroutine close_system_dipole()
    implicit none
    if(myrank==0) then
      close(file_number)
    endif
  end subroutine close_system_dipole

  subroutine calculating_system_dipole_in_each_process(istep)
    use omp_lib
    use trajectory_mpi, only: wkv, wkxyz
    use subcell, only: nselfseg, lsegtop, lseg_natoms, m2i
    use coulomb_mod, only: chgv
    use param, only: paranum
    implicit none
    include 'mpif.h'
    integer(4) :: ierr
    integer(4) :: i0, k0, ipar, istep
    real(8) ::  Mx,  My,  Mz
    real(8) :: dMx, dMy, dMz

    dMx = 0.0d0
    dMy = 0.0d0
    dMz = 0.0d0
     Mx = 0.0d0
     My = 0.0d0
     Mz = 0.0d0
!$omp parallel default(shared) &
!$omp& private(k0,i0,ipar) &
!$omp& reduction(+:dMx, dMy, dMz, Mx, My, Mz)    
!$omp do
    do k0=1,nselfseg
      do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
         ipar=paranum(m2i(i0))
         dMx = dMx + chgv(ipar) * wkv(1, i0)
         dMy = dMy + chgv(ipar) * wkv(2, i0)
         dMz = dMz + chgv(ipar) * wkv(3, i0)
         Mx  = Mx + chgv(ipar) * wkxyz(1, i0)
         My  = My + chgv(ipar) * wkxyz(2, i0)
         Mz  = Mz + chgv(ipar) * wkxyz(3, i0)
      enddo ! i0
    enddo ! k0
!$omp end do
!$omp end parallel

    wk_dM(1) = dMx
    wk_dM(2) = dMy
    wk_dM(3) = dMz
    wk_M(1) = Mx
    wk_M(2) = My
    wk_M(3) = Mz

    dM(:) = 0.0d0
    call MPI_Reduce(wk_dM, dM, 3, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(wk_M, M, 3, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if(myrank == 0) then
      if(ierr /= 0) then
        write(*,*) "ERROR(calculating_system_dipole_in_each_process): MPI_REDUCE failed"
        call flush(6)
        call mpiend()
      else
        write(file_number, '(i10, 6ES20.12)') istep, M(1), M(2), M(3), dM(1), dM(2), dM(3)
        call flush(file_number)
      endif
    endif
  end subroutine

end module system_dipole
