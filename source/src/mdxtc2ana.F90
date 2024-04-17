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
!! \brief  Program to convert xtc file
!! \author Ryo Urano
!<
  program modylas_mdxtc2ana

    use session_name_mod
    use trajectory_org
    use unit_cell
    use cell_shape
    use mpi_tool
    use md_condition
    use device_numbers
    use file_xtc

    use atom_mass
    use thermostat
    use md_const
    use device_numbers
    use session_name_mod
    use param
    use mpi_tool
    use unit_cell
    use commandline_args
    use input_mod

#ifdef XTC
    
    implicit none
    real(8):: r2,r2_r
    integer(8):: i,j
  call mpistart()


  if (myrank==0) then

     call parse_args
     call parse_input(.false.)



    call xtcf % init(trim(session_name)// '.xtc')

    do while ( xtcf % STAT == 0 )
        call xtcf % read

        mdstep=    xtcf % STEP-xtc_start
        dt = xtcf % time*1e-12/mdstep
!!! assumu cubic or rectangular box
        cellx =     xtcf % box(1,1) *1e-9
        celly =     xtcf % box(2,2) *1e-9
        cellz =     xtcf % box(3,3) *1e-9
        xyz = xtcf%pos *1e-9
!!! do analysis 

        write(*,*) "mdstep,dt,cellx,celly,cellz"
        write(*,*) mdstep,dt,cellx,celly,cellz

           call mpiend
           stop
        end do
   ! 5. Close the file
    call xtcf % close

 endif

#else 
    write(*,*)"*******************************"
    write(*,*)"you need to compile with XTC option"
    write(*,*)"*******************************"
#endif 
    call mpiend
    stop
  end program modylas_mdxtc2ana
  
