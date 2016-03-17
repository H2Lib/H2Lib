set(TESTS
   amatrix
   eigen
   h2compression
   h2matrix
   hmatrix
   laplacebem2d
   laplacebem3d
   )

include_directories(${PROJECT_SOURCE_DIR}/Library)

foreach(t ${TESTS})
  set(test_executable test_${t})
  add_executable(${test_executable} ${test_executable}.c)
  target_link_libraries(${test_executable} H2LIB ${gcov} -lm)
  add_test(${t} ${test_executable})
endforeach()


