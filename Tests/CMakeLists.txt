##############################################################
###    Include Directories and External Libraries          ###
##############################################################
include_directories( ${PROJECT_SOURCE_DIR}/Library )

find_library(M_LIB m)


##############################################################
###    Build Executables                                   ###
##############################################################
set(SRC
  test_amatrix.c test_eigen.c test_h2compression.c 
  test_h2matrix.c test_hmatrix.c test_laplacebem2d.c 
  test_laplacebem3d.c
)

foreach( testsourcefile ${SRC} )
  string( REPLACE ".c" "" testname ${testsourcefile} )
  add_executable( ${testname} ${testsourcefile} )
  target_link_libraries( ${testname} LINK_PUBLIC H2Lib ${M_LIB})
endforeach( testsourcefile ${SRC} )


