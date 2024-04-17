# cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
# project(MODYLAS NONE)
# enable_language(Fortran)
set (CMAKE_VERBOSE_MAKEFILE "On")


# find_package(xdrfile REQUIRED)
# # xdrfileが見つかった場合のみ設定を行う
# if(xdrfile)
#   # インクルードパスを指定
#   include_directories(${xdrfile_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})
#   # リンクするライブラリ指定
#   target_link_libraries(modylas ${xdrfile_LIBS} -lxdrfile)
# endif(xdrfile)


# Make library for test
# set(MODYLAS_SOURCES  
# xdr.F90
# )
# add_library( ${MODYLAS_SOURCES})

if(XTC)
# set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
#   find_package(xdrfile)
#   if(xdrfile)
#     target_link_libraries(modylas ${xdrfile_LIBS} xdrfile)
#     set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -lxdrfile")
#   endif(xdrfile)
endif(XTC)
