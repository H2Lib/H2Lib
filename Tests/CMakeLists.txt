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
   laplacebem2d
   laplacebem3d
   laplacebem3d_ocl
   tet3d
   tri2d
   )

include_directories(${PROJECT_SOURCE_DIR}/Library)

foreach(t ${TESTS})
  set(test_executable test_${t})
  add_executable(${test_executable} ${test_executable}.c)
  target_link_libraries(${test_executable} H2LIB ${gcov} -lm)
  add_test(${t} ${test_executable})
endforeach()


