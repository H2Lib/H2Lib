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
###    Subdirectories                                      ###
##############################################################
add_subdirectory (Library)
add_subdirectory (Tests)


##############################################################
###    Set Compiler Flags                                  ###
##############################################################
# Compiler definitions
#add_definitions("-O3 -Wall")

# Compiler flags
set(CMAKE_C_FLAGS "-O3 -Wall -march=native -funroll-loops -funswitch-loops ${CMAKE_C_FLAGS}")
set(CMAKE_C_FLAGS_DEBUG "-g -Wall ${CMAKE_C_FLAGS_DEBUG}")
set(CMAKE_C_FLAGS_RELEASE "-O3 -Wall ${CMAKE_C_FLAGS_RELEASE}")


##############################################################
###    Include Directories                                 ###
##############################################################
include_directories( ${PROJECT_SOURCE_DIR}/Library )

