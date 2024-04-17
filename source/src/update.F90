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
!! \brief  Module and subroutines to update coordinates and velocities
!!         by inner/outer propagator.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to update coordinates and velocities.
!! \author Yoshimichi Andoh and Noriyuki Yoshii
!<
!----------------------------------------------------------------------
module update
  use subcell
  implicit none
  real(8), parameter :: e2=1.d0/6.d0, e4=e2/20.d0, e6=e4/42.d0, e8=e6/72.d0

contains

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to update coordinates
!! \author Yoshimichi Andoh
!<
  subroutine update_coordinates(vtemp,veigv,aa2)
!----------------------------------------------------------------------
    use md_multiplestep
    use thermostat
    use trajectory_mpi
    use matrix_inverse_mod
    use atom_mass
    use param
    implicit none
    integer(4) :: i
    integer(4) :: i0,k0
    real(8) :: aa,aa2(3),arg2,poly,bb(3)
    real(8) :: veigv(3,3),veig(3),vtemp(3,3)
    real(8) :: u1,u2,u3,uv1,uv2,uv3

    vtemp = vboxg
    call diagonal(vtemp,veig,veigv)
    call matrix_inverse(veigv,vtemp)

    do i=1,3
       aa    =dexp(dthL*veig(i))
       aa2(i)=aa*aa
       arg2  =(veig(i)*dthL)*(veig(i)*dthL)
       poly  =(((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.d0
       bb(i) =aa*poly*dtL
    enddo
!$omp parallel default(none) &
!$omp& private(k0,i0,u1,u2,u3,uv1,uv2,uv3) &
!$omp& shared(nselfseg,lsegtop,lseg_natoms) &
!$omp& shared(wkxyz,wkv,vtemp,aa2,bb,veigv,mass,paranum,m2i)
!$omp do
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
#ifdef TIP4
          if(mass(paranum(m2i(i0))) .lt. 0) cycle
#endif
          u1 = wkxyz(1,i0)*vtemp(1,1) + wkxyz(2,i0)*vtemp(1,2) + wkxyz(3,i0)*vtemp(1,3)
          u2 = wkxyz(1,i0)*vtemp(2,1) + wkxyz(2,i0)*vtemp(2,2) + wkxyz(3,i0)*vtemp(2,3)
          u3 = wkxyz(1,i0)*vtemp(3,1) + wkxyz(2,i0)*vtemp(3,2) + wkxyz(3,i0)*vtemp(3,3)
          uv1 = wkv(1,i0)*vtemp(1,1)+wkv(2,i0)*vtemp(1,2)+wkv(3,i0)*vtemp(1,3)
          uv2 = wkv(1,i0)*vtemp(2,1)+wkv(2,i0)*vtemp(2,2)+wkv(3,i0)*vtemp(2,3)
          uv3 = wkv(1,i0)*vtemp(3,1)+wkv(2,i0)*vtemp(3,2)+wkv(3,i0)*vtemp(3,3)
          u1=u1*aa2(1)+uv1*bb(1)
          u2=u2*aa2(2)+uv2*bb(2)
          u3=u3*aa2(3)+uv3*bb(3)
          wkxyz(1,i0)=u1*veigv(1,1)+u2*veigv(1,2)+u3*veigv(1,3)
          wkxyz(2,i0)=u1*veigv(2,1)+u2*veigv(2,2)+u3*veigv(2,3)
          wkxyz(3,i0)=u1*veigv(3,1)+u2*veigv(3,2)+u3*veigv(3,3)
       enddo ! i0
    enddo ! k0
!$omp end do
!$omp end parallel

  end subroutine update_coordinates
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to update velocities at inner loop
!! \author Yoshimichi Andoh
!<
  subroutine update_velocities_inner()
!----------------------------------------------------------------------
    use atom_mass
    use extended_system
    use thermostat
    use md_const
    use md_condition
    use trajectory_org
    use trajectory_mpi
    use atom_virial
    use md_multiplestep
    use md_forces
    use param
    use shake_rattle_roll
    use unit_cell
    use matrix_inverse_mod
    implicit none
    integer(4) :: i,i0,k0,ipar
    real(8) :: aa,trvg
    real(8) :: uv1,uv2,uv3,u1,u2,u3
    real(8) :: veigv(3,3),veig(3),vtemp(3,3)
    real(8) :: aa2(3),arg2,poly,bb(3)

    trvg=(vboxg(1,1)+vboxg(2,2)+vboxg(3,3))*degree_of_freedom_inverse
    vtemp(1,1)=vboxg(1,1)+trvg
    vtemp(1,2)=vboxg(1,2)
    vtemp(1,3)=vboxg(1,3)
    vtemp(2,1)=0d0
    vtemp(2,2)=vboxg(2,2)+trvg
    vtemp(2,3)=vboxg(2,3)
    vtemp(3,1)=0d0
    vtemp(3,2)=0d0
    vtemp(3,3)=vboxg(3,3)+trvg
    call diagonal (vtemp,veig,veigv)
    call matrix_inverse (veigv,vtemp)

    do i=1,3
       aa    =dexp(-veig(i)*dthL/2d0)
       aa2(i)=aa*aa
       arg2  =(-veig(i)*dthL/2d0)*(-veig(i)*dthL/2d0)
       poly  =(((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.d0
       bb(i) =aa*poly*dtL/2d0
    enddo

!$omp parallel do default(none) &
!$omp& private(k0,i0,u1,u2,u3,uv1,uv2,uv3,ipar) &
!$omp& shared(nselfseg,lsegtop,lseg_natoms,m2i,paranum) &
!$omp& shared(wkv,vtemp,veigv) &
!$omp& shared(r_mass,fshort,aa2,bb,totnconst,wk_dfc,mass)
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          ipar=paranum(m2i(i0))
#ifdef TIP4
          if(mass(ipar) .lt. 0) cycle
#endif
          u1 = wkv(1,i0)*vtemp(1,1) + wkv(2,i0)*vtemp(1,2) + wkv(3,i0)*vtemp(1,3)
          u2 = wkv(1,i0)*vtemp(2,1) + wkv(2,i0)*vtemp(2,2) + wkv(3,i0)*vtemp(2,3)
          u3 = wkv(1,i0)*vtemp(3,1) + wkv(2,i0)*vtemp(3,2) + wkv(3,i0)*vtemp(3,3)
          uv1= &
               &     fshort(1,i0)*r_mass(ipar)*vtemp(1,1) &
               &    +fshort(2,i0)*r_mass(ipar)*vtemp(1,2) &
               &    +fshort(3,i0)*r_mass(ipar)*vtemp(1,3)
          uv2= &
               &     fshort(1,i0)*r_mass(ipar)*vtemp(2,1) &
               &    +fshort(2,i0)*r_mass(ipar)*vtemp(2,2) &
               &    +fshort(3,i0)*r_mass(ipar)*vtemp(2,3)
          uv3= &
               &     fshort(1,i0)*r_mass(ipar)*vtemp(3,1) &
               &    +fshort(2,i0)*r_mass(ipar)*vtemp(3,2) &
               &    +fshort(3,i0)*r_mass(ipar)*vtemp(3,3)
          if(totnconst > 0)then
             uv1=uv1+ &
                  &     wk_dfc(1,i0)*r_mass(ipar)*vtemp(1,1) &
                  &    +wk_dfc(2,i0)*r_mass(ipar)*vtemp(1,2) &
                  &    +wk_dfc(3,i0)*r_mass(ipar)*vtemp(1,3)
             uv2=uv2+ &
                  &     wk_dfc(1,i0)*r_mass(ipar)*vtemp(2,1) &
                  &    +wk_dfc(2,i0)*r_mass(ipar)*vtemp(2,2) &
                  &    +wk_dfc(3,i0)*r_mass(ipar)*vtemp(2,3)
             uv3=uv3+ &
                  &     wk_dfc(1,i0)*r_mass(ipar)*vtemp(3,1) &
                  &    +wk_dfc(2,i0)*r_mass(ipar)*vtemp(3,2) &
                  &    +wk_dfc(3,i0)*r_mass(ipar)*vtemp(3,3)
          endif
          u1=u1*aa2(1)+uv1*bb(1)
          u2=u2*aa2(2)+uv2*bb(2)
          u3=u3*aa2(3)+uv3*bb(3)
          wkv(1,i0)=u1*veigv(1,1)+u2*veigv(1,2)+u3*veigv(1,3)
          wkv(2,i0)=u1*veigv(2,1)+u2*veigv(2,2)+u3*veigv(2,3)
          wkv(3,i0)=u1*veigv(3,1)+u2*veigv(3,2)+u3*veigv(3,3)
       enddo ! i0
    enddo ! k0

  end subroutine update_velocities_inner
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to update velocities at outer loop
!! \author Yoshimichi Andoh
!<
  subroutine update_velocities_outer(dtin)
!----------------------------------------------------------------------
    use atom_mass
    use trajectory_mpi
    use md_multiplestep
    use md_forces
    use param
    implicit none
    integer(4) :: i0,k0,ipar
    real(8) :: dtin

!$omp parallel do default(none) &
!$omp& private(k0,i0,ipar) &
!$omp& shared(nselfseg,lsegtop,lseg_natoms,m2i,paranum) &
!$omp& shared(wkv,fmiddle,dtin,r_mass,mass)
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          ipar=paranum(m2i(i0))
#ifdef TIP4
          if(mass(ipar) .lt. 0) cycle
#endif
          wkv(:,i0) = wkv(:,i0) + 0.5d0*dtin*fmiddle(:,i0)*r_mass(ipar)
       enddo ! i0
    enddo ! k0

  end subroutine update_velocities_outer
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to update velocities at outermost loop
!! \author Yoshimichi Andoh
!<
  subroutine update_velocities_outermost(dtin)
!----------------------------------------------------------------------
    use atom_mass
    use trajectory_mpi
    use md_multiplestep
    use md_forces
    use param
    implicit none
    integer(4) :: i0,k0,ipar
    real(8) :: dtin

!$omp parallel do default(none) &
!$omp& private(k0,i0,ipar) &
!$omp& shared(nselfseg,lsegtop,lseg_natoms,m2i,paranum) &
!$omp& shared(wkv,flong,dtin,r_mass,mass)
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          ipar=paranum(m2i(i0))
#ifdef TIP4
          if(mass(ipar) .lt. 0) cycle
#endif
          wkv(:,i0) = wkv(:,i0) + 0.5d0*dtin*flong(:,i0)*r_mass(ipar)
       enddo ! i0
    enddo ! k0
  end subroutine update_velocities_outermost

end module update
