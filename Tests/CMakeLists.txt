set(TESTS
   amatrix
   ddcluster
   eigen
   h2compression
   h2matrix
   helmholtzbem3d
   helmholtzbem3d_ocl
   hmatrix
   krylov
   krylovsolvers
   laplacebem2d
   laplacebem3d
   laplacebem3d_ocl
   tet3d
   tet3dp1
   tet3drt0
   tri2d
   tri2dp1
   tri2drt0
   )

include_directories(${PROJECT_SOURCE_DIR}/Library)

foreach(t ${TESTS})
  set(test_executable test_${t})
  add_executable(${test_executable} ${test_executable}.c)
  target_link_libraries(${test_executable} H2LIB ${gcov} -lm)
  add_test(${t} ${test_executable})
endforeach()


