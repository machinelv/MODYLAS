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
!! \brief  Module and subroutines to calculate Ewald surface term.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to calculate Ewald surface term.
!! \author Noriyuki Yoshii, Yoshimichi Andoh, Ryo Urano
!<
module surface_term
  use omp_lib
  use mpi_tool
  implicit none
  real(8)::sysdpl(3),wk_sysdpl(3)

contains

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to calculate the dipole of unit cell
!! \author Yoshimichi Andoh
!<
  subroutine calc_system_dipole
!----------------------------------------------------------------------
    use domain, only : nlevel
    use fmm_far
    use fmm_l_m_index
    use unit_cell
    implicit none
    integer(4) :: j0,k0,m1
    complex(8) :: M1m1    ! M_1^(-1)
    complex(8) :: M10     ! M_1^0
    complex(8) :: M1p1    ! M_1^(+1)

    if(nmax.lt.1)then
       write(0,*) 'ERROR: nmax.lt.1, then dipole cannot be defined!'
       call modylas_abort()
    endif

    !### calc M_1^(-1) ###!
    j0=+1
    k0=1 ! for -1
    m1 = translate_l_m_to_1dim(j0, k0)
    M1m1=cmplx(wm(nlevel)%array(m1,1,1,1,1),-wm(nlevel)%array(m1,2,1,1,1)) * (-1)**k0

    !### calc M_1^( 0) ###!
    j0=+1
    k0= 0
    m1 = translate_l_m_to_1dim(j0, k0)
    M10 =dcmplx(wm(nlevel)%array(m1,1,1,1,1),wm(nlevel)%array(m1,2,1,1,1))

    !### calc M_1^(+1) ###!
    j0=+1
    k0=+1
    m1 = translate_l_m_to_1dim(j0, k0)
    M1p1=dcmplx(wm(nlevel)%array(m1,1,1,1,1),wm(nlevel)%array(m1,2,1,1,1))

    !### calc system dipole ###!
    sysdpl(1)=2*dreal(M1p1)
    sysdpl(2)=2*dimag(M1p1)
    sysdpl(3)=-M10

#if defined(PRECISION_M2L_MIX) || defined(PRECISION_M2L_SP)
      sysdpl = sysdpl * 1d-10
#endif

    ! debug
    !     if(myrank==0)then
    !      write(*,*) M1m1
    !      write(*,*) M10
    !      write(*,*) M1p1
    !      write(*,*) sysdpl
    !      write(*,*)
    !      do j0=0,nmax
    !      do k0=-j0,j0
    !      m1 = j0*j0+j0+1+k0
    !      write(*,*) wm_global(nlevel)%array(m1,1,1,1)!/cellx**j0
    !      enddo
    !      enddo
    !     endif
    !     call modylas_abort()
  end subroutine calc_system_dipole
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to remove the interaction by the Ewald surface term
!! \author Yoshimichi Andoh
!<
  subroutine remove_ewald_surface_term
!----------------------------------------------------------------------
    use trajectory_mpi
    use coulomb_mod
    use md_const
    use forces
    use md_monitors
    use param
    use atom_virial
    use unit_cell
    use subcell
    use nonneutral, only: qsum
    implicit none
    real(8) :: DD, coef0, coefi
    real(8) :: ewald_st, Cv
    integer(4) :: k0,i0,ipar
    integer(4) :: iam=0
    integer(4) :: k
    real(8) ::   xi_st,yi_st,zi_st,r2_st
    real(8) :: Pot_chg_st,temp,temp_Pot,coef1
    real(8) :: tempind_x,tempind_y,tempind_z
    real(8) :: sysdpl_ind_x,sysdpl_ind_y,sysdpl_ind_z

    integer(4)::ierr
    include 'mpif.h'

    Pot_chg_st=0.0d0
    temp=0.0d0
    temp_Pot=0.0d0
    
    tempind_x=0.0d0
    tempind_y=0.0d0
    tempind_z=0.0d0


    !### calc. system-dipole ###
    sysdpl=sysdpl*md_ELEMENTARY_CHARGE
    DD=sysdpl(1)*sysdpl(1)+sysdpl(2)*sysdpl(2) +sysdpl(3)*sysdpl(3)

    !### calc. force to be removed ###
    coef0=md_ELEMENTARY_CHARGE / (3d0*cellvol*md_VACUUM_DIELECTRIC_CONSTANT)
    coef1=md_ELEMENTARY_CHARGE / (6d0*cellvol*md_VACUUM_DIELECTRIC_CONSTANT)          

!$omp parallel default(none) &
!$omp& private(iam,k0,i0,ipar,coefi) &
!$omp& private(xi_st,yi_st,zi_st,r2_st) &
!$omp& shared(qsum) &
!$omp& shared(nselfseg,lsegtop,lseg_natoms) &
!$omp& shared(chgv,paranum,m2i,coef0) &
!$omp& shared(wkxyz,coef1) &
!$omp& reduction(+:Pot_chg_st) &
!$omp& reduction(+:tempind_x,tempind_y,tempind_z) &
!$omp& shared(w3_f,sysdpl)
!$  iam=omp_get_thread_num()
!$omp do
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          ipar=paranum(m2i(i0))
          coefi=chgv(ipar)*coef0
          w3_f(1,i0,iam)=w3_f(1,i0,iam)+coefi*sysdpl(1)
          w3_f(2,i0,iam)=w3_f(2,i0,iam)+coefi*sysdpl(2)
          w3_f(3,i0,iam)=w3_f(3,i0,iam)+coefi*sysdpl(3)
       enddo
    enddo
!$omp end do
!$omp end parallel

    if(myrank==0)then

       !### calc. poteintail to be removed ###
       ewald_st=DD/(6d0*cellvol*md_VACUUM_DIELECTRIC_CONSTANT)
       wk_p_energy = wk_p_energy - ewald_st
       !### calc. virial to be removed ###
       Cv=-1d0/(6d0*cellvol*md_VACUUM_DIELECTRIC_CONSTANT)
       wk_vir2(1,0)=wk_vir2(1,0)+Cv*(DD-2d0*sysdpl(1)*sysdpl(1)) ! Pxx
       wk_vir2(2,0)=wk_vir2(2,0)+Cv*(DD-2d0*sysdpl(2)*sysdpl(2)) ! Pyy
       wk_vir2(3,0)=wk_vir2(3,0)+Cv*(DD-2d0*sysdpl(3)*sysdpl(3)) ! Pzz
       wk_vir2(4,0)=wk_vir2(4,0)+Cv*(  -2d0*sysdpl(1)*sysdpl(2)) ! Pxy=Pyx
       wk_vir2(5,0)=wk_vir2(5,0)+Cv*(  -2d0*sysdpl(1)*sysdpl(3)) ! Pxz=Pzx
       wk_vir2(6,0)=wk_vir2(6,0)+Cv*(  -2d0*sysdpl(2)*sysdpl(3)) ! Pyz=Pzy
#ifdef DEBUGFCE
#ifdef KCAL
       write(*,*) 'Pot(surf) =', -ewald_st*kJ_mol/4.184d0,'[kcal/mol]'
#else
       write(*,*) 'Pot(surf) =', -ewald_st*kJ_mol,'[kJ/mol]'
#endif
#endif
    endif

  end subroutine remove_ewald_surface_term
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to remove the interaction by the Ewald surface term
!!         for non-charge neutral system.
!! \author Ryo Urano
!<
  subroutine remove_charged_ewaldsurfaceterm
!----------------------------------------------------------------------
    use trajectory_mpi
    use coulomb_mod
    use md_const
    use forces
    use md_monitors
    use param
    use atom_virial
    use unit_cell
    use subcell
    use nonneutral
    implicit none
    real(8) :: DD, coef0, coefi
    real(8) :: ewald_st, Cv
    integer(4) :: k0,i0,ipar
    integer(4) :: iam=0
    integer(4) :: k
    real(8) ::   xi_st,yi_st,zi_st,r2_st
    real(8) :: Pot_chg_st,temp,temp_Pot,coef1
    real(8) :: tempind_x,tempind_y,tempind_z
    real(8) :: sysdpl_ind_x,sysdpl_ind_y,sysdpl_ind_z

    integer(4)::ierr
    include 'mpif.h'

    Pot_chg_st=0.0d0
    temp=0.0d0
    temp_Pot=0.0d0
    
    tempind_x=0.0d0
    tempind_y=0.0d0
    tempind_z=0.0d0


    !### calc. system-dipole ###
    sysdpl=sysdpl*md_ELEMENTARY_CHARGE
    DD=sysdpl(1)*sysdpl(1)+sysdpl(2)*sysdpl(2) +sysdpl(3)*sysdpl(3)

    !### calc. force to be removed ###
    coef0=md_ELEMENTARY_CHARGE / (3d0*cellvol*md_VACUUM_DIELECTRIC_CONSTANT)
    coef1=md_ELEMENTARY_CHARGE / (6d0*cellvol*md_VACUUM_DIELECTRIC_CONSTANT)          

!$omp parallel default(none) &
!$omp& private(iam,k0,i0,ipar,coefi) &
!$omp& private(xi_st,yi_st,zi_st,r2_st) &
!$omp& shared(qsum) &
!$omp& shared(nselfseg,lsegtop,lseg_natoms) &
!$omp& shared(chgv,paranum,m2i,coef0) &
!$omp& shared(wkxyz,coef1) &
!$omp& reduction(+:Pot_chg_st) &
!$omp& reduction(+:tempind_x,tempind_y,tempind_z) &
!$omp& shared(w3_f,sysdpl)
!$  iam=omp_get_thread_num()
!$omp do
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          ipar=paranum(m2i(i0))
          coefi=chgv(ipar)*coef0

!!!!!!!!!!!!!!!! charged surface term
!!! note that surface term is subtraction term not addition term
         coefi=chgv(ipar)*coef0  !! why they donot use 1/3 instead of  1/6
         !  ! temp=temp+chgv(ipar)
          xi_st=wkxyz(1,i0)                                                    
          yi_st=wkxyz(2,i0)
          zi_st=wkxyz(3,i0)
         
          w3_f(1,i0,iam)=w3_f(1,i0,iam)-qsum*coefi*xi_st
          w3_f(2,i0,iam)=w3_f(2,i0,iam)-qsum*coefi*yi_st
          w3_f(3,i0,iam)=w3_f(3,i0,iam)-qsum*coefi*zi_st

          coefi=chgv(ipar)*coef1 !! why they donot use 1/3 instead of  1/6
          r2_st=xi_st*xi_st+yi_st*yi_st+zi_st*zi_st

!! - (- 2pi \sumi qi /3V qj rj^2 ) /(4 pi epislon) 
          Pot_chg_st=Pot_chg_st+coefi*qsum*r2_st

!! for virial calculation          
          tempind_x=tempind_x+ qsum*chgv(ipar)*xi_st*xi_st
          tempind_y=tempind_y+ qsum*chgv(ipar)*yi_st*yi_st
          tempind_z=tempind_z+ qsum*chgv(ipar)*zi_st*zi_st
       enddo
    enddo
!$omp end do
!$omp end parallel


!!! charged surface term
       wk_p_energy = wk_p_energy +Pot_chg_st

       !### calc. virial to be removed ###
       Cv=-1d0/(6d0*cellvol*md_VACUUM_DIELECTRIC_CONSTANT)

!!! assuming isotropic box
!! for qi in Qtot*\sum_i qi ! qsum is already multiplied by md_ELEMENTARY_CHARGE
       Cv=Cv*md_ELEMENTARY_CHARGE  ! for q_i in sysdpl_ind
!   endif

!!! surface charged virial 
       wk_vir2(1,0)=wk_vir2(1,0)+Cv*(tempind_x) ! Pxx
       wk_vir2(2,0)=wk_vir2(2,0)+Cv*(tempind_y) ! Pyy
       wk_vir2(3,0)=wk_vir2(3,0)+Cv*(tempind_z) ! Pzz

#ifdef DEBUGFCE
!! surface charged potential 
    call mpi_allreduce(Pot_chg_st,temp_Pot,1, &
         &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

    if(myrank==0)then
       ! write(*,*) 'DD =', DD
#ifdef KCAL
       write(*,*) 'Pot(surf_charged) =', temp_Pot*kJ_mol/4.184d0,'[kcal/mol]'
#else
       write(*,*) 'Pot(surf_charged) =', temp_Pot*kJ_mol,'[kJ/mol]'
#endif
    endif
#endif
  end subroutine remove_charged_ewaldsurfaceterm

end module surface_term

