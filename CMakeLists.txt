##############################################################
###    General                                             ###
##############################################################
# general settings
cmake_minimum_required (VERSION 2.8.11) # ???
project (H2Lib)

# The version number.
set(H2LIB_VERSION_MAJOR 1)
set(H2LIB_VERSION_MINOR 0)


##############################################################
###    Set Compiler Flags                                  ###
##############################################################
# Find Libraries
find_package (BLAS)
find_package (LAPACK)
if (BLAS_FOUND AND LAPACK_FOUND)
  set(CMAKE_C_FLAGS "-DUSEBLAS ${CMAKE_C_FLAGS}")
endif (BLAS_FOUND AND LAPACK_FOUND)


# Compiler flags
set(CMAKE_C_COMPILER gcc)
set(CMAKE_C_FLAGS "-O3 -Wall -march=native -funroll-loops -funswitch-loops ${CMAKE_C_FLAGS}")
set(CMAKE_C_FLAGS_DEBUG "-g -Wall ${CMAKE_C_FLAGS_DEBUG}")
set(CMAKE_C_FLAGS_RELEASE "-O3 -Wall ${CMAKE_C_FLAGS_RELEASE}")


##############################################################
###    Subdirectories                                      ###
##############################################################
add_subdirectory (Library)
add_subdirectory (Tests)


##############################################################
###    Include Directories                                 ###
##############################################################
include_directories( ${PROJECT_SOURCE_DIR}/Library )

