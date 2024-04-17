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
!! \brief Module and subroutines to create table functions
!<
!----------------------------------------------------------------------
!>
!! \brief Module to create table functions
!! \author Yuta Asano, and Yoshimichi Andoh
!<
module table_functions
! References:
! J. Jung et al., J.Comput.Chem., 34, 2412 (2013).
! Steinbach & Brooks, J.Comput.Chem.,15,667(1994). [eq(10)-(13)]  #CHARMMFSW
implicit none

integer(4) table_density                             !table density
real(8),save:: table_density_cut2                    !table density * cutrad2
real(8),allocatable,save:: lj12_table            (:) !r^-12 table
real(8),allocatable,save:: lj12_table_delta      (:) !r^-12 table difference
real(8),allocatable,save:: lj12_force_table      (:) !r^-14 table, -{1/(12r)}(d/dr)r^-12
real(8),allocatable,save:: lj12_force_table_delta(:) !r^-14 table difference
real(8),allocatable,save:: lj6_table             (:) !r^-6  table
real(8),allocatable,save:: lj6_table_delta       (:) !r^-6  table difference
real(8),allocatable,save:: lj6_force_table       (:) !r^-8  table, -{1/(6r)}(d/dr)r^-6
real(8),allocatable,save:: lj6_force_table_delta (:) !r^-8  table difference
real(8),allocatable,save:: cl_table              (:) !r^-1  table
real(8),allocatable,save:: cl_table_delta        (:) !r^-1  table difference
real(8),allocatable,save:: cl_force_table        (:) !r^-3  table, -(1/r)(d/dr)r^-1
real(8),allocatable,save:: cl_force_table_delta  (:) !r^-3  table difference

contains
!-----------------------------------------------------------------------
  subroutine init_table
  use ewald_variables, only : ewald_alpha
  use math_const,      only : r_PI_sqrt
  use cutoff_radius,   only : cutrad2
  use lj_mod
  use md_periodic
  use unit_cell
  use domain, only : ncellx,ncelly,ncellz
  use mpi_tool, only : mpiend
  implicit none
  integer(4) dummy_integer                        !dummy for loop
  real(8) cutoff_distance                         !cutoff distance of potential function = 12 angstrom
  real(8) cutoff_distance2                        !square of above variable
  real(8) minimum_value_of_inter_atomic_distance  !minimum of inter atomic distance
  real(8) minimum_value_of_inter_atomic_distance2 !square of above variable
  integer(4) number_of_intersection               !number of intersection
  integer(4) number_of_table_point                !number of table point
  real(8) inter_atomic_distance2                  !inter atomic distance for creating table
  real(8) interpolate_distance                    !interpolate distance of tabulated function
  real(8) table_point_real                        !table point (real)
  integer(4) table_point                          !table point (int)
!
  real(8) :: a2,alpsp2,derfc,r2
  real(8) :: xcell5,ycell5,zcell5

  a2     = ewald_alpha * ewald_alpha
  alpsp2 = 2.0d0* ewald_alpha* r_PI_sqrt 

#ifdef CHARMMFSW
  Ron2=Ron**2; Roff2=Roff**2
  Ron3=Ron**3; Roff3=Roff**3
  Ron6=Ron**6; Roff6=Roff**6
  b1=Ron3*Roff3 ; b2=Roff3-Ron3
  b3=Roff2-Ron2 ; ib3=1d0/b3**3
  c6a=    1d0/b1
  c6b=  Roff3/b2
  c6c=   -2d0/b2
  c6d=    1d0/b2/Roff3
  iRon2 =1d0/Ron2
  iRoff3=1d0/Roff3
  iRoff6=1d0/Roff6
  c12a=    1d0/b1**2
  c12b=Roff6/(Roff6-Ron6)
#endif /**CHARMMFSW**/

  cutoff_distance2 = cutrad2                       !cutoff distance ^ 2

  !!!! define table function parameters ; should be input parameter
#ifdef TIP4
  table_density                          = 20      !table density; default is 20 (ref. jung)
  minimum_value_of_inter_atomic_distance = 0.1d-10 !O-M ditance is 0.15 A
#else
  table_density                          = 20      !table density; default is 20 (ref. jung)
  minimum_value_of_inter_atomic_distance = 0.5d-10 !lower bound of interaction distance ( jung's default is 1\aa)
#endif

  !!!! calculate table function parameters
  minimum_value_of_inter_atomic_distance2 = &
     &  minimum_value_of_inter_atomic_distance ** 2    !lower bound of interaction distance ^ 2
  number_of_intersection =                &
     &  int(  cutoff_distance2           &
     &       / minimum_value_of_inter_atomic_distance2) !number of intersection of table
  number_of_table_point = table_density   &
     & * number_of_intersection !number of table point; this is the array size of table function
  !!!! allocate arrays
  allocate( lj12_table            (0:number_of_table_point))
  allocate( lj12_table_delta      (0:number_of_table_point))
  allocate( lj12_force_table      (0:number_of_table_point))
  allocate( lj12_force_table_delta(0:number_of_table_point))
  allocate( lj6_table             (0:number_of_table_point))
  allocate( lj6_table_delta       (0:number_of_table_point))
  allocate( lj6_force_table       (0:number_of_table_point))
  allocate( lj6_force_table_delta (0:number_of_table_point))
  allocate( cl_table              (0:number_of_table_point))
  allocate( cl_table_delta        (0:number_of_table_point))
  allocate( cl_force_table        (0:number_of_table_point))
  allocate( cl_force_table_delta  (0:number_of_table_point))
  !!!!zero clear of allocated arrays
  lj12_table             = 0d0
  lj12_table_delta       = 0d0
  lj12_force_table       = 0d0
  lj12_force_table_delta = 0d0
  lj6_table              = 0d0
  lj6_table_delta        = 0d0
  lj6_force_table        = 0d0
  lj6_force_table_delta  = 0d0
  cl_table               = 0d0
  cl_table_delta         = 0d0
  cl_force_table         = 0d0
  cl_force_table_delta   = 0d0

  table_density_cut2 = table_density*cutoff_distance2

  !!!!input function values into arrays
  do dummy_integer = 1, number_of_table_point, 1  !! behaive as r^2

    !calculate inter atomic distance at table point; table_point = int(table_density * rc^2 / r^2)
    inter_atomic_distance2 =       &
     &    table_density* cutoff_distance2/ dble(dummy_integer)

    !Coulomb potential (Ewald real part)
    cl_table(dummy_integer) =  & ! erfc( ewald_alpha * r ) / r ; coulomb potential part
     &    derfc( ewald_alpha*dsqrt(inter_atomic_distance2) )  &      
     &   /dsqrt(inter_atomic_distance2)  
    !Coulomb force (Ewald real part)
    cl_force_table(dummy_integer)  &
     &  = ( alpsp2*dexp ( -a2*inter_atomic_distance2 )  & 
     &    + derfc( ewald_alpha                          &
     &            *dsqrt( inter_atomic_distance2 )      &
     &           )/dsqrt( inter_atomic_distance2 )      &
     &    )/ inter_atomic_distance2  
    ! [ alpsp2 * exp( -a2 * r^2 )+ derfc( ewald_alpha * r ) / r ] / r^2; coulomb force part

#ifdef CHARMMFSW
    !LJ potential
    if(    inter_atomic_distance2 .lt. Ron2) then
    lj12_table(dummy_integer) =  &
     &  (1d0 / (inter_atomic_distance2 ** 6)) - c12a   ! r ^ -12 - c12a  ; switching lennard-jones potential first half part
    lj6_table(dummy_integer) =   &
     &  (1d0 / (inter_atomic_distance2 ** 3)) - c6a    ! r ^ -6  - c6a   ; switching lennard-jones potential first half part
    elseif(inter_atomic_distance2 .ge. Ron2 .and.  &
     &     inter_atomic_distance2 .le. Roff2)then  
    lj12_table(dummy_integer) = c12b  &
     &  * (inter_atomic_distance2 **(-3   ) - iRoff6) ** 2 ! c12b * ( r ^ -6 - iRofff6 ) ^ 2 ; switching lennard-jones potential latter half part
    lj6_table(dummy_integer) = c6b    &
     &  * (inter_atomic_distance2 **(-1.5d0) - iRoff3) ** 2 ! c6b  * ( r ^ -3 - iRoff3  ) ^ 2 ; switching lennard-jones potential latter half part
    else
        ! default zero values
    endif
    !LJ force
    if(    inter_atomic_distance2 .lt. Ron2) then
    lj12_force_table(dummy_integer) =  &
     &  1d0 / (inter_atomic_distance2 ** 7) ! r ^ -14 ; switching lennard-jones force first half part
    lj6_force_table(dummy_integer) =   &
     &  1d0 / (inter_atomic_distance2 ** 4) ! r ^ -8  ; switching lennard-jones force first half part
    elseif(inter_atomic_distance2 .ge. Ron2 .and.  &
     &     inter_atomic_distance2 .le. Roff2)then
    r2 = inter_atomic_distance2                            ! define temporaly parameter r2; r ^ 2
    Sr = (Roff2-r2)**2*(Roff2+2d0*r2-3d0*Ron2)*ib3         ! switching parameter
    lj12_force_table(dummy_integer) = 1d0 / (r2 ** 7) * Sr ! switching lennard-jones force latter half part
    lj6_force_table(dummy_integer)  = 1d0 / (r2 ** 4) * Sr ! switching lennard-jones force latter half part
    else
        ! default zero values
    endif
#else
    !LJ potential
    lj12_table(dummy_integer)       &
     &  = 1d0 / (inter_atomic_distance2 ** 6) ! r ^ -12; lennard-jones potential part
    lj6_table(dummy_integer)        &
     &  = 1d0 / (inter_atomic_distance2 ** 3) ! r ^  -6; lennard-jones potential part
    !LJ force
    lj12_force_table(dummy_integer) &
     &  = 1d0 / (inter_atomic_distance2 ** 7) ! r ^ -14; lennard-jones force part
    lj6_force_table(dummy_integer)  &
     &  = 1d0 / (inter_atomic_distance2 ** 4) ! r ^  -8; lennard-jones force part
#endif /**CHARMMFSW**/

  end do ! dummy_integer

  !!!!store array difference
  do dummy_integer = 0, number_of_table_point - 1, 1
    lj12_table_delta(dummy_integer)           &
     &  = lj12_table(dummy_integer + 1)       &
     &  - lj12_table(dummy_integer    )
    lj12_force_table_delta(dummy_integer)     &
     &  = lj12_force_table(dummy_integer + 1) &
     &  - lj12_force_table(dummy_integer    )
    lj6_table_delta(dummy_integer)            &
     &  = lj6_table(dummy_integer + 1)        &
     &  - lj6_table(dummy_integer    )
    lj6_force_table_delta(dummy_integer)      &
     &  = lj6_force_table(dummy_integer + 1)  &
     &  - lj6_force_table(dummy_integer    )
    cl_table_delta(dummy_integer)             &
     &  = cl_table(dummy_integer + 1)         &
     &  - cl_table(dummy_integer    )
    cl_force_table_delta(dummy_integer)       &
     &  = cl_force_table(dummy_integer + 1)   &
     &  - cl_force_table(dummy_integer    )
  end do

  end subroutine init_table

end module table_functions
