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
!! \brief  Modules to store global constants.
!<
!----------------------------------------------------------------------
!>
!! \brief  Modules to store mathmatical constants.
!! \author Tatsuya Sakashita
!<
module math_const
  real(8),parameter :: PI=3.1415926535897932384626433832795029d0
  real(8),parameter :: PI_sqrt=1.77245385090551602729816748334d0
  real(8),parameter :: r_PI_sqrt=1.0d0/PI_sqrt
end module math_const

!>
!! \brief  Modules to store physical constants.
!! \author Tatsuya Sakashita
!<
module physics_const
  use math_const, only : PI
  real(8),parameter :: md_AVOGADRO=6.02214199d+23
  real(8),parameter :: md_ELEMENTARY_CHARGE=1.602176462D-19
  real(8),parameter :: md_ATOMIC_MASS_UNIT=1.0d-3/md_AVOGADRO
  real(8),parameter :: md_VACUUM_DIELECTRIC_CONSTANT=8.854187817D-12
  real(8),parameter :: md_BOLTZMANN=1.3806503D-23
  real(8),parameter :: md_QQ_4PiE=md_ELEMENTARY_CHARGE **2 &
       &                                 * 0.25d0 / PI &
       &                                 /md_VACUUM_DIELECTRIC_CONSTANT
  real(8),parameter :: rvkbolz=1d0/md_BOLTZMANN
end module physics_const

!>
!! \brief  Modules to store unit conversion constants.
!! \author Tatsuya Sakashita
!<
module unit_conversion_const
  use math_const, only : PI
  use physics_const, only : md_AVOGADRO
  real(8),parameter :: md_DEGREE=PI/180.0d0
  real(8),parameter :: md_CALORIE=4.184d0
  real(8),parameter :: md_E_CONVERT=md_AVOGADRO*1.0d-3/md_CALORIE
  real(8),parameter :: md_F_CONVERT=md_AVOGADRO*1.0d-13/md_CALORIE
  real(8),parameter :: rad2deg=57.29577951308232088d0
  real(8),parameter :: deg2rad=1.0d0/rad2deg
end module unit_conversion_const

!>
!! \brief  Modules to wrap three modules.
!! \author Tatsuya Sakashita
!<
module md_const
  use math_const
  use physics_const
  use unit_conversion_const
end module md_const

!>
!! \brief  Modules to store I/O device number.
!! \author Tatsuya Sakashita
!<
module device_numbers
  integer(4),parameter :: f_rvf=10,f_trj=11,f_run=12,f_mntr=13
  integer(4),parameter :: f_force=14
  integer(4),parameter :: f_mdff=17, f_mdffbin=18
  integer(4),parameter :: f_mdbb=15
  integer(4),parameter :: f_pconst=19
  logical :: exist_mdffbin=.false.
  integer(4),parameter :: f_analef=50
  integer(4), parameter :: f_dcd=20
#ifdef DIVIDEDCD
  integer(4),parameter :: f_dcdxyz=3000
#endif
end module device_numbers

!>
!! \brief  Modules to store ensemble number.
!! \author Tatsuya Sakashita
!<
module ensemble_numbers
  integer(4),parameter :: NVE=10, NVT=20, NPT_A=30, NPT_PR=40
  integer(4),parameter :: NPT_Z=31,NPTLZ=32,NPTLZxy=33,NLXLYPZT=34,NLXLYLZT=35,NVLXLYLZT=36,NPXPYPZT=37 ! 37 is used below
  integer(4),parameter :: OPT=50
end module ensemble_numbers

!>
!! \brief  Modules to store force field number.
!! \author Tatsuya Sakashita
!<
module force_field_numbers
  integer(4),parameter :: CHARMM=100, OPLSAA=200, AMBER=300
  integer(4),parameter :: GAmbFF=310  !  General Amber FF (GAFF)
!  integer(4),parameter :: KREMER=400, PCFF = 500
  integer(4),parameter :: KREMER=400
end module force_field_numbers

!>
!! \brief  Modules to store margin parameters.
!! \author Tatsuya Sakashita
!<
module margin_sizes
  real(8), parameter :: na1cellmargin=2.00d0 ! ya  !!100% margin
  real(8), parameter :: margin_mnpc=5.0d0
  integer(4), parameter :: ndcellmargin=0
  integer(4), parameter :: bso_mgn=6  !! bsorder_margin to be transferred
end module margin_sizes
