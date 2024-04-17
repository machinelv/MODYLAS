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
!! \brief  Subroutines to integrate the equation of motion
!!         of thermostat
!! \author Yoshimichi Andoh
!<
!----------------------------------------------------------------------
!>
!! \brief  Subroutine for thermostat
!! \author Yoshimichi Andoh
!! \details This code bases on G.J.Martyna et al.,
!!          Mol.Phys., 87, 1117 (1996).
!<
module thermostat
  use mpi_tool
  implicit none
  integer(4), parameter :: mnhc=5,mnys=5,nsc=1,nnhc=5,nys=5
  real(8) :: rss(mnhc),vss(mnhc),gms(mnhc),wmts(mnhc)
  real(8) :: rssb(mnhc),vssb(mnhc),bms(mnhc),wmbs(mnhc)
  real(8) :: ttss(mnhc),bttss(mnhc),ppss(mnhc),bppss(mnhc),st=0.0d0
  real(8) :: wdti2(mnys),wdti4(mnys),wdti8(mnys)
  real(8) :: rssstr(mnhc),rssbstr(mnhc),vssstr(mnhc),vssbstr(mnhc)
  real(8) :: box(3,3),vboxg(3,3)
  ! For NtT ensemble
  real(8) :: box0(3,3)=0.d0,sigmaS(3,3)=0.d0,StressTerm(3,3)=0.d0
  logical :: use_thermo_scale=.false.,use_baro_scale=.false.
  real(8) :: okt,gnd2kt,www,dnsc
  logical :: initialize_thermostat_flag=.false.
  logical :: initialize_barostat_flag=.false.

  ! variables related to nhc_a
  integer(4) :: nbaro
  real(8) :: bkt,dnd2kt,baromass
  real(8) :: rbaromass,onethird

contains

!----------------------------------------------------------------------
  subroutine nhchain()
!----------------------------------------------------------------------
    use extended_system
    use trajectory_mpi
    use md_const
    use md_condition
    use md_multiplestep
    use kinetic_energy
    use atom_mass
    use param
    implicit none
    integer(4) :: i,j,k,i0,k0
    real(8) :: aa,totke,kinT
    real(8) :: totke_t, totke_r
    real(8) :: dkin
    real(8) :: scalev
    call k_energy_scaler(totke)

    kinT = 2.0d0 * totke &
         &       * degree_of_freedom_inverse * rvkbolz
    dkin=kinT*degree_of_freedom

    if (kinT .eq. 0.0d0)  kinT = 0.1d0

    scalev=1d0

    !     ### main part of NHC ###
    do i=1,nsc
       do j=1,nys

          gms(nnhc) =(wmts(nnhc-1)*vss(nnhc-1)**2-okt)/wmts(nnhc)
          vss(nnhc) =vss(nnhc) +gms(nnhc)*wdti4(j)

          do k=nnhc-1,2,-1
             aa     =dexp(-wdti8(j)*vss(k+1))
             gms(k) =(wmts(k-1)*vss(k-1)**2-okt)/wmts(k)
             vss(k) =vss(k)*aa*aa+gms(k)*wdti4(j)*aa
          enddo

          aa     =exp(-wdti8(j)*vss(2))
          gms(1) =(dkin-gnd2kt)/wmts(1)
          vss(1) =vss(1) *aa*aa+gms(1)*wdti4(j)*aa

          aa=exp(-wdti2(j)*vss(1))
          scalev=scalev*aa
          dkin=dkin*aa*aa

          do k=1,nnhc
             rss(k)=rss(k)+vss(k)*wdti2(j)
          enddo

          aa=exp(-wdti8(j)*vss(2))
          gms(1)=(dkin-gnd2kt)/wmts(1)
          vss(1)=vss(1)*aa*aa+gms(1)*wdti4(j)*aa
          do k=2,nnhc-1,+1
             aa     =dexp(-wdti8(j)*vss(k+1))
             gms(k) =(wmts(k-1)*vss(k-1)**2-okt)/wmts(k)
             vss(k) =vss(k)*aa*aa+gms(k)*wdti4(j)*aa
          enddo

          gms(nnhc) =(wmts(nnhc-1)*vss(nnhc-1)**2-okt)/wmts(nnhc)
          vss(nnhc) =vss(nnhc) +gms(nnhc)*wdti4(j)

       enddo
    enddo

!$omp parallel do default(shared) &
!$omp& private(k0,i0)
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
#ifdef TIP4
       if(mass(paranum(m2i(i0))) .lt. 0) cycle
#endif
          wkv(1,i0) = scalev*wkv(1,i0)
          wkv(2,i0) = scalev*wkv(2,i0)
          wkv(3,i0) = scalev*wkv(3,i0)
       end do ! i0
    end do ! k0

    !     ### Hamiltonian calculation ###
    ttss(1) =0.5d0*wmts(1)*vss(1) *vss(1)
    ppss(1) =gnd2kt*rss(1)
    do i = 2,nnhc
       ttss(i) =0.5d0*wmts(i)*vss(i) *vss(i)
       ppss(i) =okt*rss(i)
    enddo
    st=0.d0
    do i=1,nnhc
       st=st+ttss(i)+ppss(i)
    enddo
    st = st * md_BOLTZMANN
  end subroutine nhchain

!! [Ref.] G.J.Martyna et al., Mol.Phys., 87, 1117 (1996).
!! :: NOTE :: There are many misprints in the above article.
!!            So, your own check is essential.
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to initialize constants
!! \author Yoshimichi Andoh
!<
  subroutine init_nhc()
!----------------------------------------------------------------------
    use extended_system
    use md_condition
    use md_multiplestep
    implicit none
    integer(4) :: i
    real(8) :: aa

    okt   =                  systemp
    gnd2kt=degree_of_freedom    *systemp

    aa=1d0
    wmts(1)=gnd2kt/(aa*fts0*fts0)
    do i = 2, nnhc
       wmts(i)=okt/(aa*fts1*fts1)
    enddo

    dnsc = dble(nsc)
    www=1.d0/(4.d0-4.d0**(1.d0/3.d0))

    !     ### set dt for thermostat (always XO) ###
    wdti2(1)=www*dt/(2.d0*dnsc)
    wdti2(2)=wdti2(1)
    wdti2(3)=(1.d0-4.d0*www)*dt/(2.d0*dnsc)
    wdti2(4)=wdti2(1)
    wdti2(5)=wdti2(1)
    wdti4(1)=www*dt/(2.d0*dnsc)/2.d0
    wdti4(2)=wdti4(1)
    wdti4(3)=(1.d0-4.d0*www)*dt/(2.d0*dnsc)/2.d0
    wdti4(4)=wdti4(1)
    wdti4(5)=wdti4(1)
    wdti8(1)=www*dt/(2.d0*dnsc)/2.d0/2.d0
    wdti8(2)=wdti8(1)
    wdti8(3)=(1.d0-4.d0*www)*dt/(2.d0*dnsc)/2.d0/2.d0
    wdti8(4)=wdti8(1)
    wdti8(5)=wdti8(1)
  end subroutine init_nhc
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to enter zero into variables/arrays
!! \author Yoshimichi Andoh
!<
  subroutine zero_thermostat
!----------------------------------------------------------------------
    implicit none

    vss=0d0
    rss=0d0
  end subroutine zero_thermostat
!----------------------------------------------------------------------
!>
!! \brief  Subroutine for thermostat
!! \author Yoshimichi Andoh
!! \details This code bases on M.E.Tuckerman et al.,
!!          J.Phys.A:Math.Gen.,39,5629(2006).
!<
!----------------------------------------------------------------------
  subroutine nhc_a_thermo()
!----------------------------------------------------------------------
    use trajectory_org
    use trajectory_mpi
    use extended_system
    use md_const
    use md_condition
    use md_multiplestep
    use unit_cell
    use kinetic_energy
    use atom_mass
    use param
    implicit none
    integer(4) :: i,j,k,i0,k0
    real(8) :: aa,totke
    real(8) :: totke_t,totke_r
    real(8) :: dkin_part,dkin_baro
    real(8) :: tbss, pbss, sttss, sppss, bsttss, bsppss
    real(8) :: scalev,scaleb

    call k_energy_scaler(totke)
    dkin_part=2d0*totke*rvkbolz

    call k_energy_baro(vboxg,baromass,dkin_baro)

    scalev=1d0
    scaleb=1d0

    !     ### main part of NHCP ###
    do i=1,nsc
       do j=1,nys

          gms(nnhc) =(wmts(nnhc-1)*vss(nnhc-1)**2-okt)/wmts(nnhc)
          vss(nnhc) =vss(nnhc) +gms(nnhc)*wdti4(j)
          bms(nnhc) =(wmbs(nnhc-1)*vssb(nnhc-1)**2-bkt)/wmbs(nnhc)
          vssb(nnhc)=vssb(nnhc)+bms(nnhc)*wdti4(j)

          do k=nnhc-1,2,-1
             aa     =dexp(-wdti8(j)*vss(k+1))
             gms(k) =(wmts(k-1)*vss(k-1)**2-okt)/wmts(k)
             vss(k) =vss(k)*aa*aa+gms(k)*wdti4(j)*aa
             aa     =dexp(-wdti8(j)*vssb(k+1))
             bms(k) =(wmbs(k-1)*vssb(k-1)**2-bkt)/wmbs(k)
             vssb(k)=vssb(k)*aa*aa+bms(k)*wdti4(j)*aa
          enddo

          aa     =dexp(-wdti8(j)*vss(2))
          gms(1) =(dkin_part-gnd2kt)/wmts(1)
          vss(1) =vss(1) *aa*aa+gms(1)*wdti4(j)*aa
          aa     =dexp(-wdti8(j)*vssb(2))
          bms(1) =(dkin_baro-dnd2kt)/wmbs(1)
          vssb(1)=vssb(1)*aa*aa+bms(1)*wdti4(j)*aa

          do k=1,nnhc
             rss(k) =rss(k) +vss(k) *wdti2(j)
             rssb(k)=rssb(k)+vssb(k)*wdti2(j)
          enddo

          !^^^ scale atomic velocities ^^^!
          aa=dexp(-vss(1) *wdti2(j))
          scalev=scalev*aa
          dkin_part=dkin_part*aa*aa

          !^^^ scale box velocities ^^^!
          aa=dexp(-vssb(1)*wdti2(j))
          scaleb=scaleb*aa
          dkin_baro=dkin_baro*aa*aa

          aa     =dexp(-wdti8(j)*vss(2))
          gms(1) =(dkin_part-gnd2kt)/wmts(1)
          vss(1) =vss(1) *aa*aa+gms(1)*wdti4(j)*aa
          aa     =dexp(-wdti8(j)*vssb(2))
          bms(1) =(dkin_baro-dnd2kt)/wmbs(1)
          vssb(1)=vssb(1)*aa*aa+bms(1)*wdti4(j)*aa

          do k=2,nnhc-1,+1
             aa     =dexp(-wdti8(j)*vss(k+1))
             gms(k) =(wmts(k-1)*vss(k-1)**2-okt)/wmts(k)
             vss(k) =vss(k)*aa*aa+gms(k)*wdti4(j)*aa
             aa     =dexp(-wdti8(j)*vssb(k+1))
             bms(k) =(wmbs(k-1)*vssb(k-1)**2-bkt)/wmbs(k)
             vssb(k)=vssb(k)*aa*aa+bms(k)*wdti4(j)*aa
          enddo

          gms(nnhc) =(wmts(nnhc-1)*vss(nnhc-1)**2-okt)/wmts(nnhc)
          vss(nnhc) =vss(nnhc) +gms(nnhc)*wdti4(j)
          bms(nnhc) =(wmbs(nnhc-1)*vssb(nnhc-1)**2-bkt)/wmbs(nnhc)
          vssb(nnhc)=vssb(nnhc)+bms(nnhc)*wdti4(j)

       enddo
    enddo

    !     ### scale atomic velocities ###
!$omp parallel do default(none) &
!$omp& private(k0,i0) &
!$omp& shared(nselfseg,lsegtop,lseg_natoms) &
!$omp& shared(wkv,scalev,mass,paranum,m2i)
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
#ifdef TIP4
       if(mass(paranum(m2i(i0))) .lt. 0) cycle
#endif
          wkv(:,i0)=wkv(:,i0)*scalev
       enddo
    enddo

    !     ### scale box belocities ###
    vboxg(1,1)=vboxg(1,1)*scaleb
    vboxg(2,2)=vboxg(2,2)*scaleb
    vboxg(3,3)=vboxg(3,3)*scaleb

    !     ### Hamiltonian calculation ###
    ttss(1) =0.5d0*wmts(1)*vss(1) *vss(1)
    bttss(1)=0.5d0*wmbs(1)*vssb(1)*vssb(1)
    ppss(1) =gnd2kt*rss(1)
    bppss(1)=dnd2kt*rssb(1)
    do i = 2,nnhc
       ttss(i) =0.5d0*wmts(i)*vss(i) *vss(i)
       bttss(i)=0.5d0*wmbs(i)*vssb(i)*vssb(i)
       ppss(i) =okt*rss(i)
       bppss(i)=bkt*rssb(i)
    enddo
    sttss =0.0d0
    sppss =0.0d0
    bsttss=0.0d0
    bsppss=0.0d0
    do i=1,nnhc
       sttss =sttss +ttss(i)
       sppss =sppss +ppss(i)
       bsttss=bsttss+bttss(i)
       bsppss=bsppss+bppss(i)
    enddo
    tbss=0.5d0*dkin_baro  ! dkin_baro [K]
    pbss=syspres*cellvol*rvkbolz  ! syspres*cellvol [J] -> pbss [K]
    st  =sttss+sppss+bsttss+bsppss+tbss+pbss

    st = st*md_BOLTZMANN  ! [K] -> [J]

  end subroutine nhc_a_thermo
!----------------------------------------------------------------------
!>
!! \brief  Subroutine for thermostat
!! \author Yoshimichi Andoh
!! \details This code bases on M.E.Tuckerman et al.,
!!          J.Phys.A:Math.Gen.,39,5629(2006).
!<
!----------------------------------------------------------------------
  subroutine nhc_pr_thermo()
!----------------------------------------------------------------------
    use trajectory_org
    use trajectory_mpi
    use extended_system
    use md_const
    use md_condition
    use md_multiplestep
    use unit_cell
    use kinetic_energy
    use atom_mass
    use param
    implicit none
    integer(4) :: i,j,k,i0,k0
    real(8) :: aa,totke
    real(8) :: dkin_part,dkin_baro
    real(8) :: tbss, pbss, sttss, sppss, bsttss, bsppss
    real(8) :: scalev,scaleb
    ! For NtT ensemble
    real(8) :: Trans_box(3,3), matrixG(3,3), matrixSG(3,3)

    call k_energy_scaler(totke)
    call k_energy_baro(vboxg,baromass,dkin_baro)

    scalev=1d0
    scaleb=1d0
    dkin_part=2d0*totke*rvkbolz

    !     ### main part of NHCP ###
    do i=1,nsc
       do j=1,nys

          gms(nnhc) =(wmts(nnhc-1)*vss(nnhc-1)**2-okt)/wmts(nnhc)
          vss(nnhc) =vss(nnhc) +gms(nnhc)*wdti4(j)
          bms(nnhc) =(wmbs(nnhc-1)*vssb(nnhc-1)**2-bkt)/wmbs(nnhc)
          vssb(nnhc)=vssb(nnhc)+bms(nnhc)*wdti4(j)

          do k=nnhc-1,2,-1
             aa     =dexp(-wdti8(j)*vss(k+1))
             gms(k) =(wmts(k-1)*vss(k-1)**2-okt)/wmts(k)
             vss(k) =vss(k)*aa*aa+gms(k)*wdti4(j)*aa
             aa     =dexp(-wdti8(j)*vssb(k+1))
             bms(k) =(wmbs(k-1)*vssb(k-1)**2-bkt)/wmbs(k)
             vssb(k)=vssb(k)*aa*aa+bms(k)*wdti4(j)*aa
          enddo

          aa     =dexp(-wdti8(j)*vss(2))
          gms(1) =(dkin_part-gnd2kt)/wmts(1)
          vss(1) =vss(1) *aa*aa+gms(1)*wdti4(j)*aa
          aa     =dexp(-wdti8(j)*vssb(2))
          bms(1) =(dkin_baro-dnd2kt)/wmbs(1)
          vssb(1)=vssb(1)*aa*aa+bms(1)*wdti4(j)*aa

          do k=1,nnhc
             rss(k) =rss(k) +vss(k) *wdti2(j)
             rssb(k)=rssb(k)+vssb(k)*wdti2(j)
          enddo

          !^^^ scale atomic velocities ^^^!
          aa=dexp(-vss(1) *wdti2(j))
          scalev=scalev*aa
          dkin_part=dkin_part*aa*aa

          !^^^ scale box velocities ^^^!
          aa=dexp(-vssb(1)*wdti2(j))
          scaleb=scaleb*aa
          dkin_baro=dkin_baro*aa*aa

          aa     =dexp(-wdti8(j)*vss(2))
          gms(1) =(dkin_part-gnd2kt)/wmts(1)
          vss(1) =vss(1) *aa*aa+gms(1)*wdti4(j)*aa
          aa     =dexp(-wdti8(j)*vssb(2))
          bms(1) =(dkin_baro-dnd2kt)/wmbs(1)
          vssb(1)=vssb(1)*aa*aa+bms(1)*wdti4(j)*aa

          do k=2,nnhc-1,+1
             aa     =dexp(-wdti8(j)*vss(k+1))
             gms(k) =(wmts(k-1)*vss(k-1)**2-okt)/wmts(k)
             vss(k) =vss(k)*aa*aa+gms(k)*wdti4(j)*aa
             aa     =dexp(-wdti8(j)*vssb(k+1))
             bms(k) =(wmbs(k-1)*vssb(k-1)**2-bkt)/wmbs(k)
             vssb(k)=vssb(k)*aa*aa+bms(k)*wdti4(j)*aa
          enddo

          gms(nnhc) =(wmts(nnhc-1)*vss(nnhc-1)**2-okt)/wmts(nnhc)
          vss(nnhc) =vss(nnhc) +gms(nnhc)*wdti4(j)
          bms(nnhc) =(wmbs(nnhc-1)*vssb(nnhc-1)**2-bkt)/wmbs(nnhc)
          vssb(nnhc)=vssb(nnhc)+bms(nnhc)*wdti4(j)

       enddo
    enddo

    !     ### scale atomic velocities ###
!$omp parallel do default(none) &
!$omp& private(k0,i0) &
!$omp& shared(nselfseg,lsegtop,lseg_natoms) &
!$omp& shared(wkv,scalev,mass,paranum,m2i)
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
#ifdef TIP4
       if(mass(paranum(m2i(i0))) .lt. 0) cycle
#endif
          wkv(:,i0)=wkv(:,i0)*scalev
       enddo
    enddo

    !     ### scale box belocities ###
    vboxg(1,1)=vboxg(1,1)*scaleb
    vboxg(2,2)=vboxg(2,2)*scaleb
    vboxg(3,3)=vboxg(3,3)*scaleb
    vboxg(1,2)=vboxg(1,2)*scaleb
    vboxg(1,3)=vboxg(1,3)*scaleb
    vboxg(2,3)=vboxg(2,3)*scaleb

    !     ### Hamiltonian calculation ###
    ttss(1) =0.5d0*wmts(1)*vss(1) *vss(1)
    bttss(1)=0.5d0*wmbs(1)*vssb(1)*vssb(1)
    ppss(1) =gnd2kt*rss(1)
    bppss(1)=dnd2kt*rssb(1)
    do i = 2,nnhc
       ttss(i) =0.5d0*wmts(i)*vss(i) *vss(i)
       bttss(i)=0.5d0*wmbs(i)*vssb(i)*vssb(i)
       ppss(i) =okt*rss(i)
       bppss(i)=bkt*rssb(i)
    enddo
    sttss =0.0d0
    sppss =0.0d0
    bsttss=0.0d0
    bsppss=0.0d0
    do i=1,nnhc
       sttss =sttss +ttss(i)
       sppss =sppss +ppss(i)
       bsttss=bsttss+bttss(i)
       bsppss=bsppss+bppss(i)
    enddo
    tbss=0.5d0*dkin_baro  ! dkin_baro [K]
    pbss=syspres*cellvol*rvkbolz  ! syspres*cellvol [J] -> pbss [K]
    ! For NtT ensemble
    Trans_box = TRANSPOSE( box )
    matrixG   = MATMUL( Trans_box, box )
    matrixSG  = MATMUL( sigmaS, matrixG )
    pbss      = pbss + 0.5d0*( matrixSG(1,1) + matrixSG(2,2) + matrixSG(3,3) )*rvkbolz
    st  =sttss+sppss+bsttss+bsppss+tbss+pbss

    st = st*md_BOLTZMANN  ! [K] -> [J]
  end subroutine nhc_pr_thermo
!----------------------------------------------------------------------
!>
!! \brief  Subroutine for thermostat
!! \author Yoshimichi Andoh
!! \details This code bases on M.E.Tuckerman et al.,
!!          J.Phys.A:Math.Gen.,39,5629(2006).
!<
!----------------------------------------------------------------------
  subroutine nhc_z_thermo()
!----------------------------------------------------------------------
    use trajectory_org
    use trajectory_mpi
    use extended_system
    use md_const
    use md_condition
    use md_multiplestep
    use unit_cell
    use kinetic_energy
    use atom_mass
    use param
    implicit none
    integer(4) :: i,j,k,i0,k0
    real(8) :: aa,totke
    real(8) :: dkin_part,dkin_baro
    real(8) :: tbss, pbss, sttss, sppss, bsttss, bsppss
    real(8) :: scalev,scaleb

    call k_energy_scaler(totke)
    call k_energy_baro(vboxg,baromass,dkin_baro)

    scalev=1d0
    scaleb=1d0
    dkin_part=2d0*totke*rvkbolz

    !     ### main part of NHCP ###
    do i=1,nsc
       do j=1,nys

          gms(nnhc) =(wmts(nnhc-1)*vss(nnhc-1)**2-okt)/wmts(nnhc)
          vss(nnhc) =vss(nnhc) +gms(nnhc)*wdti4(j)
          bms(nnhc) =(wmbs(nnhc-1)*vssb(nnhc-1)**2-bkt)/wmbs(nnhc)
          vssb(nnhc)=vssb(nnhc)+bms(nnhc)*wdti4(j)

          do k=nnhc-1,2,-1
             aa     =dexp(-wdti8(j)*vss(k+1))
             gms(k) =(wmts(k-1)*vss(k-1)**2-okt)/wmts(k)
             vss(k) =vss(k)*aa*aa+gms(k)*wdti4(j)*aa
             aa     =dexp(-wdti8(j)*vssb(k+1))
             bms(k) =(wmbs(k-1)*vssb(k-1)**2-bkt)/wmbs(k)
             vssb(k)=vssb(k)*aa*aa+bms(k)*wdti4(j)*aa
          enddo

          aa     =dexp(-wdti8(j)*vss(2))
          gms(1) =(dkin_part-gnd2kt)/wmts(1)
          vss(1) =vss(1) *aa*aa+gms(1)*wdti4(j)*aa
          aa     =dexp(-wdti8(j)*vssb(2))
          bms(1) =(dkin_baro-dnd2kt)/wmbs(1)
          vssb(1)=vssb(1)*aa*aa+bms(1)*wdti4(j)*aa

          do k=1,nnhc
             rss(k) =rss(k) +vss(k) *wdti2(j)
             rssb(k)=rssb(k)+vssb(k)*wdti2(j)
          enddo

          !^^^ scale atomic velocities ^^^!
          aa=dexp(-vss(1) *wdti2(j))
          scalev=scalev*aa
          dkin_part=dkin_part*aa*aa

          !^^^ scale box velocities ^^^!
          aa=dexp(-vssb(1)*wdti2(j))
          scaleb=scaleb*aa
          dkin_baro=dkin_baro*aa*aa

          aa     =dexp(-wdti8(j)*vss(2))
          gms(1) =(dkin_part-gnd2kt)/wmts(1)
          vss(1) =vss(1) *aa*aa+gms(1)*wdti4(j)*aa
          aa     =dexp(-wdti8(j)*vssb(2))
          bms(1) =(dkin_baro-dnd2kt)/wmbs(1)
          vssb(1)=vssb(1)*aa*aa+bms(1)*wdti4(j)*aa

          do k=2,nnhc-1,+1
             aa     =dexp(-wdti8(j)*vss(k+1))
             gms(k) =(wmts(k-1)*vss(k-1)**2-okt)/wmts(k)
             vss(k) =vss(k)*aa*aa+gms(k)*wdti4(j)*aa
             aa     =dexp(-wdti8(j)*vssb(k+1))
             bms(k) =(wmbs(k-1)*vssb(k-1)**2-bkt)/wmbs(k)
             vssb(k)=vssb(k)*aa*aa+bms(k)*wdti4(j)*aa
          enddo

          gms(nnhc) =(wmts(nnhc-1)*vss(nnhc-1)**2-okt)/wmts(nnhc)
          vss(nnhc) =vss(nnhc) +gms(nnhc)*wdti4(j)
          bms(nnhc) =(wmbs(nnhc-1)*vssb(nnhc-1)**2-bkt)/wmbs(nnhc)
          vssb(nnhc)=vssb(nnhc)+bms(nnhc)*wdti4(j)

       enddo
    enddo

    !     ### scale atomic velocities ###
!$omp parallel do default(none) &
!$omp& private(k0,i0) &
!$omp& shared(nselfseg,lsegtop,lseg_natoms) &
!$omp& shared(wkv,scalev,mass,paranum,m2i)
    do k0=1,nselfseg
       do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
#ifdef TIP4
       if(mass(paranum(m2i(i0))) .lt. 0) cycle
#endif
          wkv(:,i0)=wkv(:,i0)*scalev
       enddo
    enddo

    !     ### scale box belocities ###
    vboxg(1,1)=vboxg(1,1)*scaleb
    vboxg(2,2)=vboxg(2,2)*scaleb
    vboxg(3,3)=vboxg(3,3)*scaleb

    !     ### Hamiltonian calculation ###
    ttss(1) =0.5d0*wmts(1)*vss(1) *vss(1)
    bttss(1)=0.5d0*wmbs(1)*vssb(1)*vssb(1)
    ppss(1) =gnd2kt*rss(1)
    bppss(1)=dnd2kt*rssb(1)
    do i = 2,nnhc
       ttss(i) =0.5d0*wmts(i)*vss(i) *vss(i)
       bttss(i)=0.5d0*wmbs(i)*vssb(i)*vssb(i)
       ppss(i) =okt*rss(i)
       bppss(i)=bkt*rssb(i)
    enddo
    sttss =0.0d0
    sppss =0.0d0
    bsttss=0.0d0
    bsppss=0.0d0
    do i=1,nnhc
       sttss =sttss +ttss(i)
       sppss =sppss +ppss(i)
       bsttss=bsttss+bttss(i)
       bsppss=bsppss+bppss(i)
    enddo
    tbss=0.5d0*dkin_baro  ! dkin_baro [K]
    pbss=syspres*cellvol*rvkbolz  ! syspres*cellvol [J] -> pbss [K]
    st  =sttss+sppss+bsttss+bsppss+tbss+pbss

    st = st*md_BOLTZMANN  ! [K] -> [J]

  end subroutine nhc_z_thermo
!----------------------------------------------------------------------
  subroutine fmod_set_rss(i0, value, line_number)
    implicit none
    integer(4), intent(in) :: i0, line_number
    real(8), intent(in) :: value
    integer(4) :: i
    i = i0 + 1
    if (i > mnhc) then
       write(0,'(a,i0)')  'ERROR (line ', line_number, ') : '
       write(0,'(a,i0,a)')  &
            &    'The number of thermostat.positions is out of bounds.  '// &
            &    'It must be less than ', mnhc, '.'
       call modylas_abort()
    endif
    rss(i) = value
  end subroutine fmod_set_rss
!----------------------------------------------------------------------
  subroutine fmod_set_vss(i0, value, line_number)
    implicit none
    integer(4), intent(in) :: i0, line_number
    real(8), intent(in) :: value
    integer(4) :: i
    i = i0 + 1
    if (i > mnhc) then
       write(0,'(a,i0)')  'ERROR (line ', line_number, ') : '
       write(0,'(a,i0,a)')  &
            &    'The number of thermostat.velocities is out of bounds.  '// &
            &    'It must be less than ', mnhc, '.'
       call modylas_abort()
    endif
    vss(i) = value
  end subroutine fmod_set_vss
!----------------------------------------------------------------------
  subroutine fmod_set_rssb(i0, value, line_number)
    implicit none
    integer(4), intent(in) :: i0, line_number
    real(8), intent(in) :: value
    integer(4) :: i
    i = i0 + 1
    if (i > mnhc) then
       write(0,'(a,i0)')  'ERROR (line ', line_number, ') : '
       write(0,'(a,i0,a)')  &
            &    'The number of barostat.positions is out of bounds.  '// &
            &    'It must be less than ', mnhc, '.'
       call modylas_abort()
    endif
    rssb(i) = value
  end subroutine fmod_set_rssb
!----------------------------------------------------------------------
  subroutine fmod_set_vssb(i0, value, line_number)
    implicit none
    integer(4), intent(in) :: i0, line_number
    real(8), intent(in) :: value
    integer(4) :: i
    i = i0 + 1
    if (i > mnhc) then
       write(0,'(a,i0)')  'ERROR (line ', line_number, ') : '
       write(0,'(a,i0,a)')  &
            &    'The number of barostat.velocities is out of bounds.  '// &
            &    'It must be less than ', mnhc, '.'
       call modylas_abort()
    endif
    vssb(i) = value
  end subroutine fmod_set_vssb
!----------------------------------------------------------------------
  subroutine fmod_set_vboxg(v11,v12,v13,v21,v22,v23,v31,v32,v33)
    implicit none
    real(8), intent(in) :: v11,v12,v13,v21,v22,v23,v31,v32,v33
    vboxg(1,1) = v11
    vboxg(2,1) = v12
    vboxg(3,1) = v13
    vboxg(1,2) = v21
    vboxg(2,2) = v22
    vboxg(3,2) = v23
    vboxg(1,3) = v31
    vboxg(2,3) = v32
    vboxg(3,3) = v33
  end subroutine fmod_set_vboxg

! For NtT ensemble
! Get reference box (cell matrix)
  subroutine fmod_referece_box(value1,value2,value3,value4,value5,value6)
    use md_const
    implicit none
    real(8), intent(in) :: value1,value2,value3,value4,value5,value6
    real(8) :: frac
    frac=dcos(value4*deg2rad)-dcos(value5*deg2rad) &
         &     *dcos(value6*deg2rad)
    box0(1,1)=value1
    box0(1,2)=value2*dcos(value6*deg2rad)
    box0(1,3)=value3*dcos(value5*deg2rad)
    box0(2,1)=0.d0
    box0(2,2)=value2*dsin(value6*deg2rad)
    box0(2,3)=value3*frac/dsin(value6*deg2rad)
    box0(3,1)=0.d0
    box0(3,2)=0.d0
    box0(3,3)=dsqrt(value3**2-(box0(1,3)**2+box0(2,3)**2))
    box0=box0*1.d-10  !! Unit: A to m
  end subroutine fmod_referece_box

  ! Get matrix Sigma
  subroutine fmod_sigmaS(value1,value2,value3,value4,value5,value6)
    use extended_system
    implicit none
    real(8) :: value1,value2,value3,value4,value5,value6
    real(8) :: det_box0
    real(8) :: Inv_box0(3,3),InvTrans_box0(3,3)
    real(8) :: tmp_matrix(3,3)
    ! Determinant of box0
    det_box0 = box0(1,1)*box0(2,2)*box0(3,3)&
         &           + box0(2,1)*box0(3,2)*box0(1,3)&
         &           + box0(3,1)*box0(1,2)*box0(2,3)&
         &           - box0(1,1)*box0(3,2)*box0(2,3)&
         &           - box0(3,1)*box0(2,2)*box0(1,3)&
         &           - box0(2,1)*box0(1,2)*box0(3,3)
    ! Inverse matrix of box0
    Inv_box0(1,1) = (box0(2,2)*box0(3,3)-box0(2,3)*box0(3,2))/det_box0
    Inv_box0(1,2) = (box0(1,3)*box0(3,2)-box0(1,2)*box0(3,3))/det_box0
    Inv_box0(1,3) = (box0(1,2)*box0(2,3)-box0(1,3)*box0(2,2))/det_box0
    Inv_box0(2,1) = (box0(2,3)*box0(3,1)-box0(2,1)*box0(3,3))/det_box0
    Inv_box0(2,2) = (box0(1,1)*box0(3,3)-box0(1,3)*box0(3,1))/det_box0
    Inv_box0(2,3) = (box0(1,3)*box0(2,1)-box0(1,1)*box0(2,3))/det_box0
    Inv_box0(3,1) = (box0(2,1)*box0(3,2)-box0(2,2)*box0(3,1))/det_box0
    Inv_box0(3,2) = (box0(1,2)*box0(3,1)-box0(1,1)*box0(3,2))/det_box0
    Inv_box0(3,3) = (box0(1,1)*box0(2,2)-box0(1,2)*box0(2,1))/det_box0
    ! t - P_ext*1
    tmp_matrix(1,1) = -value1 - syspres
    tmp_matrix(1,2) = -value2
    tmp_matrix(1,3) = -value3
    tmp_matrix(2,1) = -tmp_matrix(1,2)
    tmp_matrix(2,2) = -value4 - syspres
    tmp_matrix(2,3) = -value5
    tmp_matrix(3,1) = -tmp_matrix(1,3)
    tmp_matrix(3,2) = -tmp_matrix(2,3)
    tmp_matrix(3,3) = -value6 - syspres
    ! Eqn.: sigmaS = det[box0]*Inv[box0]*( t - P_ext*1 )Inv[Trans[box0]]
    ! > Note that Inv[Tans[A]] = Trans[Inv[A]]
    InvTrans_box0 = TRANSPOSE( Inv_box0 )
    sigmaS = MATMUL( tmp_matrix, InvTrans_box0 )
    sigmaS = MATMUL( Inv_box0, sigmaS )
    sigmaS = det_box0*sigmaS

  end subroutine fmod_sigmaS

!----------------------------------------------------------------------
!>
!! \brief  Subroutine to initialize constants
!! \author Yoshimichi Andoh
!<
!----------------------------------------------------------------------
  subroutine init_nhca()
    use extended_system
    use md_const
    use ensemble_numbers
    use md_condition
    implicit none
    integer(4) :: i
    real(8) :: aa

    nbaro  = 3
    bkt    =systemp
    dnd2kt =dble(nbaro)*bkt     !! x,y,z (normal)
    if(      md_condition__ensemble == NPTLZ )then
       dnd2kt =dble(nbaro-1)*bkt   !! x,y
    elseif(  md_condition__ensemble == NLXLYPZT )then
       dnd2kt =dble(nbaro-2)*bkt   !! x,y
    elseif(  md_condition__ensemble == NPXPYPZT )then
       dnd2kt =dble(nbaro-0)*bkt   !! x,y
    elseif(  md_condition__ensemble == NLXLYLZT )then
       dnd2kt =dble(nbaro-3)*bkt   !! x,y
    elseif(  md_condition__ensemble == NPTLZxy )then
       dnd2kt =dble(nbaro-1)*bkt   !! x,y
    endif

    aa=1.d0
    wmbs(1)=dnd2kt/(aa*fbs0*fbs0)
    do i = 2, nnhc
       wmbs(i)=bkt/(aa*fbs1*fbs1)
    enddo

    baromass = dble(degree_of_freedom+nbaro)*systemp &
  &            / (dble(nbaro)*aa*fbs*fbs)

    rbaromass=1d0/baromass

    onethird=1d0/3d0
  end subroutine init_nhca

end module thermostat
