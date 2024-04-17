set (CMAKE_VERBOSE_MAKEFILE "On")

# timer
if(USE_MAPROF)
  include_directories(${PROJECT_BINARY_DIR}/benchmark/ma_prof/src) # for Fujitsu compiler
endif(USE_MAPROF)

# Fujitsu inline
if(INLINE_FUJITSU)
  set(INLINE_FUJITSU_FORTRAN_FLAGS "-x200 -xdir=${PROJECT_SOURCE_DIR}/src")
  message(STATUS "Use inline option for Fujitsu Fortran compiler: " ${INLINE_FUJITSU_FORTRAN_FLAGS})
  set_property(SOURCE regular_solid_harmonics_cartesian.f90 APPEND_STRING PROPERTY COMPILE_FLAGS " -xregular_solid_harmonics_cartesian.calculate_regular_harmonics_1dim_array")
  set_property(SOURCE fmm_l_m_index.f90 APPEND_STRING PROPERTY COMPILE_FLAGS " -xfmm_l_m_index.translate_l_m_to_1dim")
  set_property(SOURCE fmm_far.F90 APPEND_STRING PROPERTY COMPILE_FLAGS " ${INLINE_FUJITSU_FORTRAN_FLAGS} -xfmm_l_m_index.translate_l_m_to_1dim -xregular_solid_harmonics_cartesian.calculate_regular_harmonics_1dim_array")
endif(INLINE_FUJITSU)

if(XTC)
include_directories(../../external/libxdrfile/include)
link_directories(../../source/external/libxdrfile/build/lib/)
endif(XTC)


# Creating static library
set(LIB_SOURCES session_name.F90 random_seed.F90 commandline_args.F90 version.F90 global_constants.F90 global_variables.F90 param.F90 trajectory_org.F90 trajectory_mpi.F90 segments.F90 molecules.F90 mol_info.F90 atom_mass.F90 md_periodic.F90 md_condition.F90 ewald_variables.F90 unit_cell.F90 tip4p.F90 file_mdrun.F90 file_mdtrj.F90 file_mdmntr.F90 file_restart.F90 file_dcd.F90 file_dcd_divided.F90 app.F90 file_force.F90 file_mdxyz_bin.F90
  precision_mod.F90
  timer_mdls_mpi.F90
  force_wrap.F90 force_short.F90 force_middle.F90 force_long.F90 fmm_ewald.F90 fmm_far.F90
  coulomb.F90 lj.F90 table_funct.F90
  ewald.F90 pme_far.F90
  MAP_module.F90 CMAP.F90
  nve.F90 nvt.F90 kinetic_energy.F90
  external/mt19937ar.f external/pzfft3dv.f external/ffte.f
  mpi_tool.F90 domain.F90 dr_cntl.F90 subcell.F90 boundary.F90 dist_atom2cell.F90
  center_of_mass_variables.F90 shake_rattle_roll.F90
  npt_pr.F90 npt_z.F90 update.F90 npt_a.F90 pressure.F90 thermostat.F90 barostat.F90 cell_shape.F90
  array.F90 hash.F90 parse.F90 parse_input.F90 file_mdff.F90 file_xtc.F90
  external/diagonal.f matrix_inverse.F90
  comm_bound.F90 comm_fmm.F90 comm_pme.F90
  comm_direct2_dr.F90 comm_direct3_dr.F90
  fmm_near_dr.F90 pme_near_dr.F90
  opt.F90 position_constrain.F90 system_dipole.F90
  center_of_mass.F90 bond.F90 bond_morse.F90 bond_morse2.F90 angle.F90 angle_morse.F90 UB.F90 dihedral.F90 improper_torsion.F90 bond_breaking.F90 parse_special.F90 void123.F90 special14.F90 spherical_harmonics.F90 surface_term.F90 display_log.F90 md_oplsaa_special_divide.F90 file_utility.F90 file_mdbb.F90 file_application.F90 atom_virial.F90 mpi_3d_grid.F90 regular_singular_solid_harmonics.F90 math_functions.F90 fmm_l_m_index.f90 regular_solid_harmonics_cartesian.f90 fmm_subcell_index.f90 fmm_parameters.f90 openmp_tool.F90
regular_singular_rotation_real_1dim.F90 regular_singular_solid_harmonics_north.F90 profile_comm.F90
)

if(GROEXT)
  set(LIB_SOURCES ${LIB_SOURCES} groext/file_groext.F90 groext/ConvertGmxModylas.C groext/GromacsConverter.C groext/GromacsData.C 
    groext/ModylasData.C groext/ModylasTagData.C groext/ModylasUtil.C groext/String.C )
endif(GROEXT)

if(USE_MAPROF)
  set(LIB_SOURCES ${LIB_SOURCES} timer_open.F90 timer_close.F90)
endif(USE_MAPROF)

if(XTC)
  set(LIB_SOURCES ${LIB_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/../external/xdr.F90)
else ()
  set(LIB_SOURCES ${LIB_SOURCES} xdr_dummy.F90)
endif(XTC)

if(FJ_RDMA)
  set(LIB_SOURCES ${LIB_SOURCES} fj_rdma.F90 fj_rdma_lib.c)
endif(FJ_RDMA)

if(PROFILE_COMM)
  set(LIB_SOURCES ${LIB_SOURCES} profile_comm.F90)
endif(PROFILE_COMM)

add_library(modylas_lib ${LIB_SOURCES})
if(USE_MAPROF)
  list(APPEND LIST_PREPROCESSOR_DEF "PROF_MAPROF")
  target_link_libraries(modylas_lib ma_prof_lib)
endif(USE_MAPROF)

# if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
# set_source_files_properties (pme_far.F90 PROPERTIES COMPILE_FLAGS    "-fpe0") 
# set_source_files_properties (pme_far.F90 PROPERTIES COMPILE_FLAGS    "-O0") 
# endif()
## Lines for usage of PMlib & papi
#include_directories(/path/to/pmlib/include)
#include_directories(/path/to/papi/include)
#target_link_libraries(modylas_lib /path/to/pmlib/lib/libPMmpi.a)
#target_link_libraries(modylas_lib /path/to/pmlib/lib/libpapi_ext.a)
#target_link_libraries(modylas_lib /path/to/papi/lib/libpapi.a)
#target_link_libraries(modylas_lib /path/to/papi/lib/libpfm.a)
#target_link_libraries(modylas_lib stdc++)

if(GROEXT)
  target_link_libraries(modylas_lib stdc++)
endif(GROEXT)

target_compile_definitions(modylas_lib PUBLIC ${LIST_PREPROCESSOR_DEF})
get_target_property(TMP_COMPILE_DEFINITIONS modylas_lib COMPILE_DEFINITIONS)
message(STATUS "Detected #ifdef variables in CMAKE_Fortran_FLAGS: ${TMP_COMPILE_DEFINITIONS}")
set_target_properties(modylas_lib PROPERTIES OUTPUT_NAME modylas)

# Building executable files
add_executable(modylas main.F90)
target_link_libraries(modylas modylas_lib)

if(XTC)
message(STATUS "Set  XTC library option, Link libxdrfile.so or .a")
message(STATUS "export LD_LIBRARY_PATH=${PROJECT_SOURCE_DIR}/external/libxdrfile/build/lib/:\${LD_LIBRARY_PATH}")
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -lxdrfile")
# target_link_libraries(modylas_lib xdrfile)
else ()
message(STATUS "NO setting of  XTC library")
endif(XTC)



add_executable(modylas-text2bin modylas-text2bin.F90)
target_link_libraries(modylas-text2bin modylas_lib)

add_executable(modylas-mdtrj2xyz modylas-mdtrj2xyz.F90)
target_link_libraries(modylas-mdtrj2xyz modylas_lib)

add_executable(modylas-mdtrj2xtc mdtrj2xtc.F90)
add_executable(modylas-mdxtc2ana mdxtc2ana.F90)
target_link_libraries(modylas-mdtrj2xtc  modylas_lib)
target_link_libraries(modylas-mdxtc2ana modylas_lib)

#add_executable(modylas-readmdtrj_mbarout readmdtrj2mbarinput.F90)
#target_link_libraries(modylas-readmdtrj_mbarout  modylas_lib)

## Program for check_accuracy
#add_executable(check_accuracy check_accuracy.F90)
#target_link_libraries(check_accuracy modylas_lib)

if(GROEXT)
  if(${CMAKE_Fortran_COMPILER} MATCHES "frtpx")
    set_target_properties(modylas           PROPERTIES LINKER_LANGUAGE "CXX")
    set_target_properties(modylas-text2bin  PROPERTIES LINKER_LANGUAGE "CXX")
    set_target_properties(modylas-mdtrj2xyz PROPERTIES LINKER_LANGUAGE "CXX")
    set_target_properties(modylas-mdtrj2xtc PROPERTIES LINKER_LANGUAGE "CXX")
    set_target_properties(modylas-mdxtc2ana PROPERTIES LINKER_LANGUAGE "CXX")
#   set_target_properties(modylas-readmdtrj_mbarout PROPERTIES LINKER_LANGUAGE "CXX")
#   set_target_properties(check_accuracy    PROPERTIES LINKER_LANGUAGE "CXX")
  else(${CMAKE_Fortran_COMPILER} MATCHES "frtpx")
    set_target_properties(modylas PROPERTIES LINKER_LANGUAGE "Fortran")
    set_target_properties(modylas-text2bin PROPERTIES LINKER_LANGUAGE "Fortran")
    set_target_properties(modylas-mdtrj2xyz PROPERTIES LINKER_LANGUAGE "Fortran")
    set_target_properties(modylas-mdtrj2xtc PROPERTIES LINKER_LANGUAGE "Fortran")
    set_target_properties(modylas-mdxtc2ana PROPERTIES LINKER_LANGUAGE "Fortran")
#   set_target_properties(modylas-readmdtrj_mbarout PROPERTIES LINKER_LANGUAGE "Fortran")
#   set_target_properties(check_accuracy PROPERTIES LINKER_LANGUAGE "Fortran")
  endif(${CMAKE_Fortran_COMPILER} MATCHES "frtpx")
endif(GROEXT)

# Installing
# set(PROGS modylas modylas-text2bin modylas-mdtrj2xyz)
#set(PROGS modylas modylas-text2bin modylas-mdtrj2xyz  modylas-mdtrj2xtc modylas-mdxtc2ana modylas-readmdtrj_mbarout)
set(PROGS modylas modylas-text2bin modylas-mdtrj2xyz  modylas-mdtrj2xtc modylas-mdxtc2ana)
install(TARGETS ${PROGS} DESTINATION bin)


