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
!! \brief  Module and subroutines to control .mdbb file output.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module to control .mdbb file output.
!! \author Zhiye Tang
!<
module file_mdbb
    implicit none
  
  contains
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to open the mdbb file
!! \author Zhiye Tang
!<  
    subroutine open_mdbb
      use device_numbers
      use session_name_mod
      use mpi_tool
      use file_utility, only : create_file
      implicit none
  
      call create_file(f_mdbb, trim(session_name)// '.mdbb')
    end subroutine open_mdbb
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to write to the mdbb file
!! \author Zhiye Tang
!<  
    subroutine write_mdbb_entry(istep, atomID1, atomID2)
        use device_numbers
        implicit none
        integer(4) :: istep, atomID1, atomID2

        write(f_mdbb,'(A3,I10,A24,I8,A5,I8,A9)') "At ", istep, " bond breaking between ", atomID1, " and ", atomID2, " (1-base)"
        flush(f_mdbb)
    end subroutine write_mdbb_entry
  end module file_mdbb
  
