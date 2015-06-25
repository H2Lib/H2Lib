set(SRC
  aca.c amatrix.c avector.c basic.c bem2d.c bem3d.c 
  block.c clusterbasis.c cluster.c clustergeometry.c 
  clusteroperator.c curve2d.c eigensolvers.c factorizations.c 
  gaussquad.c h2arith.c h2compression.c h2matrix.c h2update.c 
  harith.c hcoarsen.c hmatrix.c krylov.c laplacebem2d.c
  laplacebem3d.c macrosurface3d.c parameters.c rkmatrix.c 
  settings.c singquad1d.c singquad2d.c sparsematrix.c 
  sparsepattern.c surface3d.c truncation.c uniform.c
)

# Create a library called "Hello" which includes the source file "hello.cxx".
# The extension is already found. Any number of sources could be listed here.
add_library (H2Lib SHARED ${SRC})

# Make sure the compiler can find include files for our Hello library
# when other libraries or executables link to Hello
target_include_directories (H2Lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
