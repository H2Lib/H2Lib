
# ------------------------------------------------------------
# Components of the main library
# ------------------------------------------------------------

H2LIB_CORE0 = \
	Library/basic.c \
	Library/settings.c \
	Library/parameters.c \
	Library/opencl.c

H2LIB_CORE1 = \
	Library/avector.c \
	Library/realavector.c \
	Library/amatrix.c \
	Library/factorizations.c \
	Library/eigensolvers.c \
	Library/sparsematrix.c \
	Library/sparsepattern.c \
	Library/gaussquad.c \
	Library/krylov.c

H2LIB_CORE2 = \
	Library/cluster.c \
	Library/clustergeometry.c \
	Library/block.c \
	Library/clusterbasis.c \
	Library/clusteroperator.c \
	Library/uniform.c \
	Library/h2matrix.c \
	Library/rkmatrix.c \
	Library/hmatrix.c \
	Library/krylovsolvers.c

H2LIB_CORE3 = \
	Library/truncation.c \
	Library/harith.c \
	Library/harith2.c \
	Library/hcoarsen.c \
	Library/h2compression.c \
	Library/h2update.c \
	Library/h2arith.c \
	Library/aca.c \
	Library/visualize.c \
	Library/matrixnorms.c

H2LIB_DIRECTIONAL = \
	Library/dcluster.c \
	Library/dblock.c \
	Library/dclusterbasis.c \
	Library/duniform.c \
	Library/dh2matrix.c \
	Library/dh2compression.c \
	Library/dclusteroperator.c

H2LIB_SIMPLE =

H2LIB_FEM = \
	Library/tri2d.c \
	Library/tri2dp1.c \
	Library/tet3d.c \
	Library/tet3dp1.c\
	Library/ddcluster.c\
	Library/tri2drt0.c\
	Library/tet3drt0.c

H2LIB_BEM = \
	Library/curve2d.c \
	Library/singquad1d.c \
	Library/bem2d.c \
	Library/laplacebem2d.c \
	Library/surface3d.c \
	Library/macrosurface3d.c \
	Library/singquad2d.c \
	Library/bem3d.c \
	Library/oclbem3d.c \
	Library/laplacebem3d.c \
	Library/laplaceoclbem3d.c \
	Library/helmholtzbem3d.c \
	Library/helmholtzoclbem3d.c

SOURCES_libh2 := \
	$(H2LIB_CORE0) \
	$(H2LIB_CORE1) \
	$(H2LIB_CORE2) \
	$(H2LIB_CORE3) \
	$(H2LIB_DIRECTIONAL) \
	$(H2LIB_SIMPLE) \
	$(H2LIB_FEM) \
	$(H2LIB_BEM)

EXTRA_HEADERS = \
	Library/blas.h \
	Library/clsettings.h \
	Library/simd_avx.h \
	Library/simd_sse2.h \
	Library/simd.h

HEADERS_libh2 := $(SOURCES_libh2:.c=.h) $(EXTRA_HEADERS)

OBJECTS_libh2 := $(SOURCES_libh2:.c=.o)

DEPENDENCIES_libh2 := $(SOURCES_libh2:.c=.d)

# ------------------------------------------------------------
# Test programs
# ------------------------------------------------------------

SOURCES_tests := \
	Tests/test_amatrix.c \
	Tests/test_eigen.c \
	Tests/test_hmatrix.c \
	Tests/test_h2matrix.c \
	Tests/test_krylov.c \
	Tests/test_krylovsolvers.c \
	Tests/test_laplacebem2d.c \
	Tests/test_laplacebem3d.c \
	Tests/test_laplacebem3d_ocl.c \
	Tests/test_helmholtzbem3d.c \
	Tests/test_helmholtzbem3d_ocl.c \
	Tests/test_h2compression.c \
	Tests/test_tet3dp1.c \
	Tests/test_tri2dp1.c\
	Tests/test_ddcluster.c\
	Tests/test_tri2drt0.c\
	Tests/test_tet3drt0.c

OBJECTS_tests := \
	$(SOURCES_tests:.c=.o)

DEPENDENCIES_tests := \
	$(SOURCES_tests:.c=.d)

PROGRAMS_tests := \
	$(SOURCES_tests:.c=)

# ------------------------------------------------------------
# Example programs
# ------------------------------------------------------------

SOURCES_examples = \
	Examples/example_dh2matrix.c \
	Examples/example_visualize.c \
	Examples/example_amatrix_bem3d.c \
	Examples/example_hmatrix_bem3d.c \
	Examples/example_h2matrix_bem3d.c\
	Examples/example_fem_p1_harith.c\
	Examples/example_fem_rt0_harith.c

OBJECTS_examples = \
	$(SOURCES_examples:.c=.o)

DEPENDENCIES_examples = \
	$(SOURCES_examples:.c=.d)

PROGRAMS_examples = \
	$(SOURCES_examples:.c=)

# ------------------------------------------------------------
# All files
# ------------------------------------------------------------

SOURCES := \
	$(SOURCES_libh2) \
	$(SOURCES_tests) \
	$(SOURCES_examples)

HEADERS := \
	$(HEADER_libh2)

OBJECTS := \
	$(OBJECTS_libh2) \
	$(OBJECTS_tests) \
	$(OBJECTS_examples)

DEPENDENCIES := \
	$(DEPENDENCIES_libh2) \
	$(DEPENDENCIES_tests) \
	$(DEPENDENCIES_examples)

PROGRAMS := \
	$(PROGRAMS_tests) \
	$(PROGRAMS_examples)

# ------------------------------------------------------------
# Standard target
# ------------------------------------------------------------

all: $(PROGRAMS_tests) $(PROGRAMS_examples)

# ------------------------------------------------------------
# OpenCL sources
# ------------------------------------------------------------

Library/laplaceoclbem3d.c: Library/laplacebem3d.cl

Library/helmholtzoclbem3d.c: Library/helmholtzbem3d.cl

# ------------------------------------------------------------
# Build configuration
# ------------------------------------------------------------

ifeq ($(wildcard options.inc),)
$(OBJECTS): options.inc.default
include options.inc.default
else
$(OBJECTS): options.inc
include options.inc
endif

# ------------------------------------------------------------
# System-dependent parameters (e.g., name of compiler)
# ------------------------------------------------------------

ifeq ($(wildcard system.inc),)
$(OBJECTS): system.inc.linux
include system.inc.linux
else
$(OBJECTS): system.inc
include system.inc
endif

# ------------------------------------------------------------
# System-independent configuration (e.g., variants of algorithms)
# ------------------------------------------------------------

ifdef HARITH_RKMATRIX_QUICK_EXIT
CFLAGS += -DHARITH_RKMATRIX_QUICK_EXIT
endif

ifdef HARITH_AMATRIX_QUICK_EXIT
CFLAGS += -DHARITH_AMATRIX_QUICK_EXIT
endif

ifdef TRACE_MEMORY
CFLAGS += -DTRACE_MEMORY
endif

# ------------------------------------------------------------
# Rules for test programs
# ------------------------------------------------------------

$(PROGRAMS_tests): %: %.o
ifdef BRIEF_OUTPUT
	@echo Linking $@
	@$(CC) $(LDFLAGS) $< -o $@ -lh2 $(LIBS) 
else
	$(CC) $(LDFLAGS) $< -o $@ -lh2 $(LIBS) 
endif

$(PROGRAMS_tests) $(PROGRAMS_tools): libh2.a

$(OBJECTS_tests): %.o: %.c
ifdef BRIEF_OUTPUT
	@echo Compiling $<
	@$(GCC) -MT $@ -MM -I Library $< > $(<:%.c=%.d)
	@$(CC) $(CFLAGS) -I Library -c $< -o $@
else
	@$(GCC) -MT $@ -MM -I Library $< > $(<:%.c=%.d)
	$(CC) $(CFLAGS) -I Library -c $< -o $@
endif

-include $(DEPENDENCIES_tests) $(DEPENDENCIES_tools)
$(OBJECTS_tests): Makefile

# ------------------------------------------------------------
# Rules for example programs
# ------------------------------------------------------------

$(PROGRAMS_examples): %: %.o
ifdef BRIEF_OUTPUT
	@echo Linking $@
	@$(CC) $(LDFLAGS) $< -o $@ -lh2 $(LIBS) 
else
	$(CC) $(LDFLAGS) $< -o $@ -lh2 $(LIBS) 
endif

$(PROGRAMS_examples): libh2.a

$(OBJECTS_examples): %.o: %.c
ifdef BRIEF_OUTPUT
	@echo Compiling $<
	@$(GCC) -MT $@ -MM -I Library $< > $(<:%.c=%.d)
	@$(CC) $(CFLAGS) -I Library -c $< -o $@
else
	@$(GCC) -MT $@ -MM -I Library $< > $(<:%.c=%.d)
	$(CC) $(CFLAGS) -I Library -c $< -o $@
endif

-include $(DEPENDENCIES_examples)
$(OBJECTS_examples): Makefile

# ------------------------------------------------------------
# Rules for the Doxygen documentation
# ------------------------------------------------------------

doc:
	doxygen Doc/Doxyfile

# ------------------------------------------------------------
# Rules for the main library
# ------------------------------------------------------------

libh2.a: $(OBJECTS_libh2)
ifdef BRIEF_OUTPUT
	@echo Building $@
	@$(AR) $(ARFLAGS) $@ $(OBJECTS_libh2)
else
	$(AR) $(ARFLAGS) $@ $(OBJECTS_libh2)
endif

$(OBJECTS_libh2): %.o: %.c
ifdef BRIEF_OUTPUT
	@echo Compiling $<
	@$(GCC) -MT $@ -MM $< > $(<:%.c=%.d)
	@$(CC) $(CFLAGS) -c $< -o $@
else
	@$(GCC) -MT $@ -MM $< > $(<:%.c=%.d)
	$(CC) $(CFLAGS) -c $< -o $@
endif

-include $(DEPENDENCIES_libh2)
$(OBJECTS_libh2): Makefile

# ------------------------------------------------------------
# Useful additions
# ------------------------------------------------------------

.PHONY: clean cleandoc programs indent

clean:
	$(RM) -f $(OBJECTS) $(DEPENDENCIES) $(PROGRAMS) libh2.a

cleandoc:
	$(RM) -rf Doc/html Doc/latex

indent:
	indent -bap -br -nce -cdw -npcs \
	  -di10 -nbc -brs -blf -i2 -lp \
	  -T amatrix -T pamatrix -T pcamatrix \
	  -T avector -T pavector -T pcavector \
	  -T cluster -T pcluster -T pccluster \
	  -T block -T pblock -T pcblock \
	  -T rkmatrix -T prkmatrix -T pcrkmatrix \
	  -T hmatrix -T phmatrix -T pchmatrix \
	  -T uniform -T puniform -T pcuniform \
	  -T h2matrix -T ph2matrix -T pch2matrix \
	  -T dcluster -T pdcluster -T pcdcluster \
	  -T dblock -T pdblock -T pcdblock \
	  -T duniform -T pduniform -T pcduniform \
	  -T dh2matrix -T pdh2matrix -T pcdh2matrix \
	  $(SOURCES)

coverage:
	mkdir Coverage > /dev/null 2>&1; \
	lcov --base-directory . --directory . --capture \
	--output-file Coverage/coverage.info && \
	genhtml -o Coverage Coverage/coverage.info

cleangcov:
	$(RM) -rf Library/*.gcov Library/*.gcda Library/*.gcno \
	Tests/*.gcov Tests/*.gcda Tests/*.gcno;
	$(RM) -rf Coverage
