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
!! \brief  Module and subroutines for converting from GROMACS inputs
!!         to MODYLAS inputs.
!<
!----------------------------------------------------------------------
!>
!! \brief  Module for converting from GROMACS inputs to MODYLAS inputs.
!! \author Tetsuro Nagai
!<
module file_groext
  use iso_c_binding
  implicit none
  public parse_groext
  public convert_gmx_to_modylas
  public convert_modylas_to_gmx
  public record_gro
  public record_top
  public record_mdxyz
  public record_mdff
  private output_notice

  ! input parameters
  logical :: gmx_input_groext  = .false.
  logical :: gmx_output_groext = .false.
  character(LEN=100) :: gmx_include_path = '../top'
  logical :: gmx_read_top     = .false.
  logical :: gmx_read_gro     = .false.
  logical :: gmx_write_top    = .false.
  logical :: gmx_write_gro    = .false.
  logical :: gmx_write_mdff   = .false.
  logical :: gmx_write_mdxyz  = .false.
  logical :: gmx_write_mddef  = .false.
  logical :: gmx_convert_only = .false.
  character(LEN=1000) :: gmx_inc_keys = ""
  logical :: gmx_read_mdp     = .false.
  integer :: gmx_constraints  = 2

  character(:,c_char), pointer :: str_mdxyz
  character(:,c_char), pointer :: str_mdff
  character(:,c_char), pointer :: str_gro
  character(:,c_char), pointer :: str_top

  interface
     function c_get_mdxyz_ptr() bind(c)
       import c_ptr
       type(c_ptr) c_get_mdxyz_ptr
     end function c_get_mdxyz_ptr
     
     function c_get_mdff_ptr() bind(c)
       import c_ptr
       type(c_ptr) c_get_mdff_ptr
     end function c_get_mdff_ptr
     
     function c_get_gro_ptr() bind(c)
       import c_ptr
       type(c_ptr) c_get_gro_ptr
     end function c_get_gro_ptr
     
     function c_get_top_ptr() bind(c)
       import c_ptr
       type(c_ptr) c_get_top_ptr
     end function c_get_top_ptr

  end interface
  
contains
  
  subroutine parse_groext()
    use parse
    use mpi_tool
    use session_name_mod
    implicit none
    logical :: lerr
    character(LEN=100) :: chara
    character(LEN=1000) :: keys
    integer :: i, j, n, access

    call parse_get('/input/input_groext', chara, lerr)
    gmx_input_groext = (.not. lerr .and. trim(chara) == 'yes')
    
    call parse_get('/input/output_groext', chara, lerr)
    gmx_output_groext = (.not. lerr .and. trim(chara) == 'yes')

    if( .not. gmx_input_groext .and. .not. gmx_output_groext ) return

    if( gmx_input_groext .and. gmx_output_groext ) then
       if(myrank==0) write(0,*)"ERROR: input_groext and output_groext are exclusive"
       call modylas_abort
    end if

    call parse_get('/input/groext/include_path', chara, lerr)
    if( .not. lerr ) gmx_include_path = trim(chara)

    gmx_read_top = gmx_input_groext
    gmx_read_gro = gmx_input_groext

    call parse_get('/input/groext/read_top', chara, lerr)
    if( .not. lerr ) gmx_read_top = (trim(chara) == 'yes')

    if( gmx_input_groext .and. .not. (gmx_read_top .and. gmx_read_gro) ) then
       if(myrank==0) write(0,*)"ERROR: read_top and read_gro must be yes when input_groext is yes"
       call modylas_abort
    end if
    if( .not. gmx_input_groext .and. (gmx_read_top .or. gmx_read_gro) ) then
       if(myrank==0) write(0,*)"ERROR: read_top and read_gro must be no when input_groext is no"
       call modylas_abort
    end if

    call parse_get('/input/groext/read_gro', chara, lerr)
    if( .not. lerr ) gmx_read_gro = (trim(chara) == 'yes')

    gmx_write_top = gmx_output_groext
    gmx_write_gro = gmx_output_groext

    call parse_get('/input/groext/write_top', chara, lerr)
    if( .not. lerr ) gmx_write_top = (trim(chara) == 'yes')

    call parse_get('/input/groext/write_gro', chara, lerr)
    if( .not. lerr ) gmx_write_gro = (trim(chara) == 'yes')

    if( gmx_output_groext .and. gmx_write_gro ) then
       if( 0 < access(trim(session_name)//'.gro', "r") ) then
          if(myrank==0) write(0,*)"ERROR: gro file is required when output_groext and write_gro are yes"
          call modylas_abort
       end if
    end if

    gmx_write_mdff  = gmx_input_groext
    gmx_write_mdxyz = gmx_input_groext

    call parse_get('/input/groext/write_mdff', chara, lerr)
    if( .not. lerr ) gmx_write_mdff = (trim(chara) == 'yes')
    
    call parse_get('/input/groext/write_mdxyz', chara, lerr)
    if( .not. lerr ) gmx_write_mdxyz = (trim(chara) == 'yes')
    
    call parse_get('/input/groext/write_mddef', chara, lerr)
    if( .not. lerr ) gmx_write_mddef = (trim(chara) == 'yes')
    
    call parse_get('/input/groext/convert_only', chara, lerr)
    if( .not. lerr ) gmx_convert_only = (trim(chara) == 'yes')
    
    gmx_inc_keys = ""
    call parse_get('/input/groext/include_keys', chara, lerr)
    if( .not. lerr .and. trim(chara) == 'yes') then
       do
          call parse_get('/input/groext/inc_keys', chara, lerr)
          if( lerr ) exit
          gmx_inc_keys = trim(gmx_inc_keys) // " " // trim(chara)
       end do
    end if

    call parse_get('/input/groext/read_mdp', chara, lerr)
    if( .not. lerr ) gmx_read_mdp = (trim(chara) == 'yes')

    call parse_get('/input/groext/constraints', chara, lerr)
    if( .not. lerr ) then
       gmx_constraints = 0
       if(trim(chara) == 'none'     ) gmx_constraints = 1
       if(trim(chara) == 'hbonds'   ) gmx_constraints = 2
       if(trim(chara) == 'all-bonds') gmx_constraints = 3
       if( gmx_constraints == 0 ) then
          if(myrank==0) write(0,*)"ERROR: Invalid constraints value:", trim(chara), "Select from none, hbonds or all-bonds."
          call modylas_abort
       end if
    end if

  end subroutine parse_groext

!!>
!!! \brief  Subroutine to read gromacs gro and top file.
!!! \author 
!!<

  function f2c(c_str) result(f_str)
    character(*, c_char), intent(in)  :: c_str
    character(:, c_char), allocatable :: f_str
    f_str = c_str // c_null_char
  end function f2c

  function c2f(cptr) result(fptr)
    type(c_ptr), intent(in), value :: cptr
    character(:, c_char), pointer :: fptr
    interface
      function strlen(p) bind(C)
        import c_ptr,c_size_t
        type(c_ptr), value :: p
        integer(c_size_t) strlen
      end function
    end interface
    fptr => convert_c2f(cptr,strlen(cptr))
  contains
    function convert_c2f(p,len)
      type(c_ptr)            , intent(In) :: p
      integer(c_size_t)      , intent(In) :: len
      character(len, c_char) , pointer    :: convert_c2f
      call c_f_pointer(p,convert_c2f)
    end function
  end function

  subroutine convert_gmx_to_modylas()
    use session_name_mod
    use mpi_tool
    implicit none
    integer(c_int) :: i_read_mdp, i_write_mddef, i_constraints

    interface
       subroutine c_convert_gmx_modylas(c_session_name, c_include_path, &
            c_inc_keys, c_i_read_mdp, c_i_write_mddef, c_i_constraints) bind(c)
         import c_int
         character :: c_session_name(*), c_include_path(*), c_inc_keys(*)
         integer(c_int), value :: c_i_read_mdp, c_i_write_mddef, c_i_constraints
       end subroutine c_convert_gmx_modylas
    end interface

    if( gmx_read_mdp ) then
       i_read_mdp = 1
    else
       i_read_mdp = 0
    end if

    if( gmx_write_mddef ) then
       i_write_mddef = 1
    else
       i_write_mddef = 0
    end if       

    i_constraints = gmx_constraints

    call c_convert_gmx_modylas( &
         f2c(trim(session_name)), &
         f2c(trim(gmx_include_path)), &
         f2c(trim(gmx_inc_keys)), &
         i_read_mdp, i_write_mddef, i_constraints)
    
    str_mdxyz => c2f(c_get_mdxyz_ptr())
    str_mdff  => c2f(c_get_mdff_ptr())
    str_gro   => c2f(c_get_gro_ptr())
    str_top   => c2f(c_get_top_ptr())

    if( myrank .eq. mpiout ) then
       if( gmx_write_mdxyz .and. gmx_convert_only ) then
          call record_mdxyz
       end if
       if( gmx_write_mdff ) then
          call record_mdff
       end if
    end if

    call output_notice()

  end subroutine convert_gmx_to_modylas

  subroutine convert_modylas_to_gmx()
    use session_name_mod
    use mpi_tool
    implicit none

    interface
       subroutine c_convert_modylas_gmx(c_session_name) bind(c)
         character c_session_name(*)
       end subroutine c_convert_modylas_gmx
    end interface

    call c_convert_modylas_gmx(f2c(trim(session_name)))
    
    str_mdxyz => c2f(c_get_mdxyz_ptr())
    str_mdff  => c2f(c_get_mdff_ptr())
    str_gro   => c2f(c_get_gro_ptr())
    str_top   => c2f(c_get_top_ptr())

    if( myrank .eq. mpiout ) then
       if( gmx_write_gro ) then
          call update_gro_names
          call record_gro
       end if
       if( gmx_write_top ) then
          call record_top
       end if
    end if

    call output_notice()

  end subroutine convert_modylas_to_gmx

  subroutine record_gro
    use session_name_mod
    use file_utility, only : open_file
    implicit none
    integer(4), parameter :: f = 105

    call open_file(f, trim(session_name)//'.gro.rst')
    write(f, '(a$)') str_gro
    close(f)

  end subroutine record_gro

  subroutine record_top
    use session_name_mod
    use file_utility, only : open_file
    implicit none
    integer(4), parameter :: f = 106

    call open_file(f, trim(session_name)//'.top.rst')
    write(f, '(a$)') str_top
    close(f)

  end subroutine record_top

  subroutine update_gro_positions()
    use trajectory_org
    use unit_cell
    implicit none

    interface
       subroutine c_update_gro_positions(xyz, v, cell, n) bind(c)
         import
         integer(c_size_t), value :: n
         real(c_double), intent(in) :: xyz(n), v(n), cell(6)
       end subroutine c_update_gro_positions
    end interface

    real(8) :: cell(6)

    cell(1) = cellx
    cell(2) = celly
    cell(3) = cellz
    cell(4) = alpha
    cell(5) = beta
    cell(6) = gamma

    call c_update_gro_positions(xyz, v, cell, size(xyz,kind=c_size_t))
    
    str_gro => c2f(c_get_gro_ptr())

  end subroutine update_gro_positions

  subroutine update_gro_names()
    use session_name_mod
    implicit none

    interface
       subroutine c_update_gro_names(c_session_name) bind(c)
         character c_session_name(*)
       end subroutine c_update_gro_names
    end interface

    call c_update_gro_names(f2c(trim(session_name)))
    
    str_gro => c2f(c_get_gro_ptr())

  end subroutine update_gro_names

  subroutine record_mdff
    use session_name_mod
    use file_utility, only : open_file
    implicit none
    integer(4), parameter :: f = 105

    call open_file(f, trim(session_name)//'.mdff.rst')
    write(f, '(a$)') str_mdff
    close(f)

  end subroutine record_mdff

   subroutine record_mdxyz
    use session_name_mod
    use file_utility, only : open_file
    implicit none
    integer(4), parameter :: f = 106

    call open_file(f, trim(session_name)//'.mdxyz.rst')
    write(f, '(a$)') str_mdxyz
    close(f)

  end subroutine record_mdxyz

  subroutine output_notice()
    use session_name_mod
    use mpi_tool
    implicit none
    character(LEN=1000) :: chara

    if (myrank==0) then
       write(6,*) "GROEXT was executed. "

       chara = ""
       if( gmx_read_mdp ) chara = trim(chara) // " " // trim(session_name) // '.mdp'
       if( gmx_read_gro ) chara = trim(chara) // " " // trim(session_name) // '.gro'
       if( gmx_read_top ) chara = trim(chara) // " " // trim(session_name) // '.top'

       write(6,*) "Used input files: " // trim(chara)

       chara = ""
       if( gmx_write_gro   ) chara = trim(chara) // " " // trim(session_name) // '.gro'
       if( gmx_write_top   ) chara = trim(chara) // " " // trim(session_name) // '.top'
       if( gmx_write_mdff  ) chara = trim(chara) // " " // trim(session_name) // '.mdff.rst'
       if( gmx_write_mdxyz ) chara = trim(chara) // " " // trim(session_name) // '.mdxyz.rst'
       if( gmx_write_mddef ) chara = trim(chara) // " " // trim(session_name) // '.mddef.rst'

       write(6,*) "Generated files: " // trim(chara)

       chara = ""
       if( trim(gmx_include_path) /= '../top' ) chara = trim(chara) // " " // 'include_path=' // trim(gmx_include_path)
       if( trim(gmx_inc_keys    ) /= ''       ) chara = trim(chara) // " " // 'inc_keys='     // trim(gmx_inc_keys)
       if( gmx_convert_only                   ) chara = trim(chara) // " " // 'convert_only=yes'
       if( gmx_constraints == 1               ) chara = trim(chara) // " " // 'constraints=none'
       if( gmx_constraints == 2               ) chara = trim(chara) // " " // 'constraints=hbonds'
       if( gmx_constraints == 3               ) chara = trim(chara) // " " // 'constraints=all-bonds'

       write(6,*) "The option specified: " // trim(chara)
       write(6,*) ""



       write(6,*) "NOTES related to GROEXT: "
       write(6,*) "At this point, users are responsible for confirming that the input gro file ", &
                  "is wrapped at the level coarser than segments. The segments are generated ",   &
                  "from charge group of the s topology file. ",                                   &
                  " Atoms connected by shake should belong to the same segment, ",                &
                  "and segments should not be greater than functional groups ",                   &
                  "for numerical accuracy." 
       
       write(6,*) ""
       write(6,*) "It would be perhaps recommendable to perform " 
       write(6,*) "`gmx trajconv -f original.gro -o wrap_by_mol.gro -pbc mol` "
       write(6,*) "beforehand to make a gro file that is wrapped with respect to each molecule ", &
                  "(which is of course coarser than segments).  It should be also noted that " ,  &
                  "all coordinates should be in the original cell or the first image cells. "  ,  &
                  "Any particles should not be in the second or further neighboring image cells. "

       write(6,*) ""
       write(6,*) "For charmm force field, sgima_11^special and epsilon_11^special in <lj special pair> ",    &
                  "are generated ad hoc from sigma_14 and epsilon_14 and defaults sigma_44 and epsilon_44. " ,&
                  "Thus the parameters in this tag may not coincide with the CHARMM's original values. "

       write(6,*) ""
       write(6,*) "For oplsaa force field, all dihedral potential entries with function type = 1 are assumed to be ",       &
                  "improper dihedral, while entries with function type = 3 (and 9) are considered to be proper dihedral. ", &
                  "Should not these be the case, 1-4 interactions may be wrong. " ,                                         &
                  "You may need to manually switch the corresponding entries."

       write(6,*) ""
       write(6,*) "For GAFF force field (or Amber FF)" ,  & 
                  "All dihedral potential entries with function type = 1 will be converted to IMPROPER dihedral, ",      &
                  "as per gromacs's (somewhat weired) OLPS convention. ",        &
                  "Dihedral potential entries with function type = 4 will be correctly converted into improper dihedral.",&
                  "Entries with function type = 3 (and 9) are considered to be proper dihedral. "
    end if
  end subroutine output_notice


end module file_groext
