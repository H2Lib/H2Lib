#ifdef USE_OPENCL
#ifdef USE_OPENMP

#include "basic.h"
#include "krylov.h"
#include "laplaceoclbem3d.h"
#include "matrixnorms.h"

#ifdef USE_CAIRO
#include <cairo.h>
#endif

static uint problems = 0;

#define IS_IN_RANGE(a, b, c) (((a) <= (b)) && ((b) <= (c)))

/* Simple convenience wrapper for conjugate gradient solver */
static void
solve_cg_bem3d(matrixtype type, void *A, pavector b, pavector x,
	       real accuracy, uint steps)
{
  uint      n, iter;

  n = b->dim;
  assert(x->dim == n);

  random_real_avector(x);

  switch (type) {
  case AMATRIX:
    iter = solve_cg_amatrix_avector((pcamatrix) A, b, x, accuracy, steps);
    break;
  case HMATRIX:
    iter = solve_cg_hmatrix_avector((pchmatrix) A, b, x, accuracy, steps);
    break;
  case H2MATRIX:
    iter = solve_cg_h2matrix_avector((pch2matrix) A, b, x, accuracy, steps);
    break;
  default:
    iter = 0;
    printf("ERROR: unknown matrix type!\n");
    abort();
    break;
  }
  printf("CG iterations:\n");
  printf("  %d\n", iter);

}

static void
test_hmatrix_system(const char *apprxtype, pcamatrix Vfull,
		    pcamatrix KMfull, pblock block, pbem3d bem_slp,
		    phmatrix V, pbem3d bem_dlp, phmatrix KM, bool linear,
		    bool exterior, real low, real high)
{
  pavector  x, b;
  real      errorV, errorKM, error_solve, eps_solve;
  uint      steps;

  eps_solve = 1.0e-12;
  steps = 1000;

  printf("Testing: %s Hmatrix %s\n"
	 "====================================\n\n",
	 (exterior == true ? "exterior" : "interior"), apprxtype);

  SCHEDULE_OPENCL(0, 1, assemble_bem3d_hmatrix, bem_slp, block, V);
  SCHEDULE_OPENCL(0, 1, assemble_bem3d_hmatrix, bem_dlp, block, KM);

  errorV = norm2diff_amatrix_hmatrix(V, Vfull) / norm2_amatrix(Vfull);
  printf("rel. error V       : %.5e\n", errorV);
  errorKM = norm2diff_amatrix_hmatrix(KM, KMfull) / norm2_amatrix(KMfull);
  printf("rel. error K%c0.5*M : %.5e\n", (exterior == true ? '-' : '+'),
	 errorKM);

  x = new_avector(Vfull->rows);
  b = new_avector(KMfull->cols);

  printf("Solving Dirichlet problem:\n");

  if (linear == true) {
    projectL2_bem3d_l_avector(bem_dlp,
			      (exterior ==
			       true ? eval_dirichlet_fundamental2_laplacebem3d
			       : eval_dirichlet_fundamental_laplacebem3d), x,
			      NULL);
  }
  else {
    projectL2_bem3d_c_avector(bem_dlp,
			      (exterior ==
			       true ? eval_dirichlet_fundamental2_laplacebem3d
			       : eval_dirichlet_fundamental_laplacebem3d), x,
			      NULL);
  }
  clear_avector(b);
  addeval_hmatrix_avector(1.0, KM, x, b);

  solve_cg_bem3d(HMATRIX, V, b, x, eps_solve, steps);

  if (linear == true) {
    error_solve = normL2diff_l_bem3d(bem_slp, x,
				     (exterior ==
				      true ?
				      eval_neumann_fundamental2_laplacebem3d :
				      eval_neumann_fundamental_laplacebem3d),
				     NULL);
  }
  else {
    error_solve = normL2diff_c_bem3d(bem_slp, x,
				     (exterior ==
				      true ?
				      eval_neumann_fundamental2_laplacebem3d :
				      eval_neumann_fundamental_laplacebem3d),
				     NULL);
  }

  error_solve = error_solve
    / normL2_bem3d(bem_slp,
		   (exterior ==
		    true ? eval_neumann_fundamental2_laplacebem3d :
		    eval_neumann_fundamental_laplacebem3d), NULL);

  printf("rel. error neumann : %.5e       %s\n", error_solve,
	 (IS_IN_RANGE(low, error_solve, high) ? "    okay" : "NOT okay"));

  if (!IS_IN_RANGE(low, error_solve, high))
    problems++;

  printf("\n");

  del_avector(x);
  del_avector(b);
}

static void
test_h2matrix_system(const char *apprxtype, pcamatrix Vfull,
		     pcamatrix KMfull, pblock block, pbem3d bem_slp,
		     ph2matrix V, pbem3d bem_dlp, ph2matrix KM, bool linear,
		     bool exterior, real low, real high)
{
  pavector  x, b;
  real      errorV, errorKM, error_solve, eps_solve;
  uint      steps;

  eps_solve = 1.0e-12;
  steps = 1000;

  printf("Testing: %s H2matrix %s\n"
	 "====================================\n\n",
	 (exterior == true ? "exterior" : "interior"), apprxtype);

  assemble_bem3d_h2matrix_row_clusterbasis(bem_slp, V->rb);
  assemble_bem3d_h2matrix_col_clusterbasis(bem_slp, V->cb);
  SCHEDULE_OPENCL(0, 1, assemble_bem3d_h2matrix, bem_slp, V);

  assemble_bem3d_h2matrix_row_clusterbasis(bem_dlp, KM->rb);
  assemble_bem3d_h2matrix_col_clusterbasis(bem_dlp, KM->cb);
  SCHEDULE_OPENCL(0, 1, assemble_bem3d_h2matrix, bem_dlp, KM);

  errorV = norm2diff_amatrix_h2matrix(V, Vfull) / norm2_amatrix(Vfull);
  printf("rel. error V       : %.5e\n", errorV);
  errorKM = norm2diff_amatrix_h2matrix(KM, KMfull) / norm2_amatrix(KMfull);
  printf("rel. error K%c0.5*M : %.5e\n", (exterior == true ? '-' : '+'),
	 errorKM);

  x = new_avector(Vfull->rows);
  b = new_avector(KMfull->cols);

  printf("Solving Dirichlet problem:\n");

  if (linear == true) {
    projectL2_bem3d_l_avector(bem_dlp,
			      (exterior ==
			       true ? eval_dirichlet_fundamental2_laplacebem3d
			       : eval_dirichlet_fundamental_laplacebem3d), x,
			      NULL);
  }
  else {
    projectL2_bem3d_c_avector(bem_dlp,
			      (exterior ==
			       true ? eval_dirichlet_fundamental2_laplacebem3d
			       : eval_dirichlet_fundamental_laplacebem3d), x,
			      NULL);
  }
  clear_avector(b);
  addeval_h2matrix_avector(1.0, KM, x, b);

  solve_cg_bem3d(H2MATRIX, V, b, x, eps_solve, steps);

  if (linear == true) {
    error_solve = normL2diff_l_bem3d(bem_slp, x,
				     (exterior ==
				      true ?
				      eval_neumann_fundamental2_laplacebem3d :
				      eval_neumann_fundamental_laplacebem3d),
				     NULL);
  }
  else {
    error_solve = normL2diff_c_bem3d(bem_slp, x,
				     (exterior ==
				      true ?
				      eval_neumann_fundamental2_laplacebem3d :
				      eval_neumann_fundamental_laplacebem3d),
				     NULL);
  }

  error_solve = error_solve
    / normL2_bem3d(bem_slp,
		   (exterior ==
		    true ? eval_neumann_fundamental2_laplacebem3d :
		    eval_neumann_fundamental_laplacebem3d), NULL);

  printf("rel. error neumann : %.5e       %s\n", error_solve,
	 (IS_IN_RANGE(low, error_solve, high) ? "    okay" : "NOT okay"));

  if (!IS_IN_RANGE(low, error_solve, high))
    problems++;

  printf("\n");

  del_avector(x);
  del_avector(b);
}

int
main(int argc, char **argv)
{
  pmacrosurface3d mg;
  psurface3d gr;
  pamatrix  Vfull, KMfull;
  pbem3d    bem_slp, bem_dlp;
  pcluster  root;
  pblock block;
  phmatrix  V, KM;
  pclusterbasis Vrb, Vcb, KMrb, KMcb;
  ph2matrix V2, KM2;
  uint      n, q, clf, m, l;
  real      eta, delta, eps_aca;
  cl_device_id *devices;
  cl_uint   ndevices;
  uint     *idx;
  uint      i;

  init_h2lib(&argc, &argv);

  get_opencl_devices(&devices, &ndevices);
  ndevices = 1;
  set_opencl_devices(devices, ndevices, 2);

  n = 512;
  q = 2;
  clf = 16;
  eta = 1.0;

  mg = new_sphere_macrosurface3d();
  gr = build_from_macrosurface3d_surface3d(mg, REAL_SQRT(n * 0.125));
  n = gr->triangles;

  printf("Testing unit sphere with %d triangles\n", n);

  bem_slp =
    new_slp_laplace_ocl_bem3d(gr, q, q + 2, BASIS_CONSTANT_BEM3D,
			      BASIS_CONSTANT_BEM3D);
  bem_dlp =
    new_dlp_laplace_ocl_bem3d(gr, q, q + 2, BASIS_CONSTANT_BEM3D,
			      BASIS_CONSTANT_BEM3D, 0.5);
  root = build_bem3d_cluster(bem_slp, clf, BASIS_CONSTANT_BEM3D);
  block = build_nonstrict_block(root, root, &eta, admissible_max_cluster);

  max_pardepth = 0;

  Vfull = new_amatrix(n, n);
  KMfull = new_amatrix(n, n);
  idx = allocuint(n);
  for (i = 0; i < n; ++i) {
    idx[i] = i;
  }
  SCHEDULE_OPENCL(0, 1, bem_slp->nearfield, idx, idx, bem_slp, false, Vfull);
  SCHEDULE_OPENCL(0, 1, bem_dlp->nearfield, idx, idx, bem_dlp, false, KMfull);

  V = build_from_block_hmatrix(block, 0);
  KM = build_from_block_hmatrix(block, 0);

  printf("----------------------------------------\n");
  printf("Testing inner Boundary integral equations:\n");
  printf("----------------------------------------\n\n");

  /*
   * Test Interpolation
   */

  m = 4;

  setup_hmatrix_aprx_inter_row_bem3d(bem_slp, root, root, block, m);
  setup_hmatrix_aprx_inter_row_bem3d(bem_dlp, root, root, block, m);
  test_hmatrix_system("Interpolation row", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, false, false, 7.0e-2, 7.5e-2);

  setup_hmatrix_aprx_inter_col_bem3d(bem_slp, root, root, block, m);
  setup_hmatrix_aprx_inter_col_bem3d(bem_dlp, root, root, block, m);
  test_hmatrix_system("Interpolation column", Vfull, KMfull, block, bem_slp,
		      V, bem_dlp, KM, false, false, 7.0e-2, 7.5e-2);

  setup_hmatrix_aprx_inter_mixed_bem3d(bem_slp, root, root, block, m);
  setup_hmatrix_aprx_inter_mixed_bem3d(bem_dlp, root, root, block, m);
  test_hmatrix_system("Interpolation mixed", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, false, false, 7.0e-2, 7.5e-2);

  /*
   * Test Green
   */

  m = 5;
  l = 1;
  delta = 0.5;

  setup_hmatrix_aprx_green_row_bem3d(bem_slp, root, root, block, m, l, delta,
				     build_bem3d_cube_quadpoints);
  setup_hmatrix_aprx_green_row_bem3d(bem_dlp, root, root, block, m, l, delta,
				     build_bem3d_cube_quadpoints);
  test_hmatrix_system("Green row", Vfull, KMfull, block, bem_slp, V, bem_dlp,
		      KM, false, false, 7.0e-2, 7.5e-2);

  setup_hmatrix_aprx_green_col_bem3d(bem_slp, root, root, block, m, l, delta,
				     build_bem3d_cube_quadpoints);
  setup_hmatrix_aprx_green_col_bem3d(bem_dlp, root, root, block, m, l, delta,
				     build_bem3d_cube_quadpoints);
  test_hmatrix_system("Green column", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, false, false, 7.0e-2, 7.5e-2);

  setup_hmatrix_aprx_green_mixed_bem3d(bem_slp, root, root, block, m, l,
				       delta, build_bem3d_cube_quadpoints);
  setup_hmatrix_aprx_green_mixed_bem3d(bem_dlp, root, root, block, m, l,
				       delta, build_bem3d_cube_quadpoints);
  test_hmatrix_system("Green mixed", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, false, false, 7.0e-2, 7.5e-2);

  /*
   * Test Greenhybrid
   */

  m = 2;
  l = 1;
  delta = 1.0;
  eps_aca = 1.0e-2;

  setup_hmatrix_aprx_greenhybrid_row_bem3d(bem_slp, root, root, block, m, l,
					   delta, eps_aca,
					   build_bem3d_cube_quadpoints);
  setup_hmatrix_aprx_greenhybrid_row_bem3d(bem_dlp, root, root, block, m, l,
					   delta, eps_aca,
					   build_bem3d_cube_quadpoints);
  test_hmatrix_system("Greenhybrid row", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, false, false, 7.0e-2, 7.5e-2);

  setup_hmatrix_aprx_greenhybrid_col_bem3d(bem_slp, root, root, block, m, l,
					   delta, eps_aca,
					   build_bem3d_cube_quadpoints);
  setup_hmatrix_aprx_greenhybrid_col_bem3d(bem_dlp, root, root, block, m, l,
					   delta, eps_aca,
					   build_bem3d_cube_quadpoints);
  test_hmatrix_system("Greenhybrid column", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, false, false, 7.0e-2, 7.5e-2);

  setup_hmatrix_aprx_greenhybrid_mixed_bem3d(bem_slp, root, root, block, m, l,
					     delta, eps_aca,
					     build_bem3d_cube_quadpoints);
  setup_hmatrix_aprx_greenhybrid_mixed_bem3d(bem_dlp, root, root, block, m, l,
					     delta, eps_aca,
					     build_bem3d_cube_quadpoints);
  test_hmatrix_system("Greenhybrid mixed", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, false, false, 7.0e-2, 7.5e-2);

  /*
   * Test ACA / PACA / HCA
   */

  m = 3;
  eps_aca = 1.0e-2;

  /* Nearfield computation on GPU not applicable here yet! */

//  setup_hmatrix_aprx_aca_bem3d(bem_slp, root, root, block, eps_aca);
//  setup_hmatrix_aprx_aca_bem3d(bem_dlp, root, root, block, eps_aca);
//  test_hmatrix_system("ACA full pivoting", Vfull, KMfull, block, bem_slp, V,
//      bem_dlp, KM, false, false, 7.0e-2, 7.5e-2);
//
//  setup_hmatrix_aprx_paca_bem3d(bem_slp, root, root, block, eps_aca);
//  setup_hmatrix_aprx_paca_bem3d(bem_dlp, root, root, block, eps_aca);
//  test_hmatrix_system("ACA partial pivoting", Vfull, KMfull, block, bem_slp, V,
//      bem_dlp, KM, false, true, 7.0e-2, 7.5e-2);
  setup_hmatrix_aprx_hca_bem3d(bem_slp, root, root, block, m, eps_aca);
  setup_hmatrix_aprx_hca_bem3d(bem_dlp, root, root, block, m, eps_aca);
  test_hmatrix_system("HCA2", Vfull, KMfull, block, bem_slp, V, bem_dlp, KM,
		      false, false, 7.0e-2, 7.5e-2);

  del_hmatrix(V);
  del_hmatrix(KM);
  del_block(block);

  /*
   * H2-matrix
   */

  block = build_strict_block(root, root, &eta, admissible_max_cluster);

  Vrb = build_from_cluster_clusterbasis(root);
  Vcb = build_from_cluster_clusterbasis(root);
  KMrb = build_from_cluster_clusterbasis(root);
  KMcb = build_from_cluster_clusterbasis(root);

  V2 = build_from_block_h2matrix(block, Vrb, Vcb);
  KM2 = build_from_block_h2matrix(block, KMrb, KMcb);

  /*
   * Test Interpolation
   */

  m = 3;

  setup_h2matrix_aprx_inter_bem3d(bem_slp, Vrb, Vcb, block, m);
  setup_h2matrix_aprx_inter_bem3d(bem_dlp, KMrb, KMcb, block, m);
  test_h2matrix_system("Interpolation", Vfull, KMfull, block, bem_slp, V2,
		       bem_dlp, KM2, false, false, 7.0e-2, 7.5e-2);

  /*
   * Test Greenhybrid
   */

  m = 2;
  l = 1;
  delta = 1.0;
  eps_aca = 6.0e-3;

  setup_h2matrix_aprx_greenhybrid_bem3d(bem_slp, Vrb, Vcb, block, m, l, delta,
					eps_aca, build_bem3d_cube_quadpoints);
  setup_h2matrix_aprx_greenhybrid_bem3d(bem_dlp, KMrb, KMcb, block, m, l,
					delta, eps_aca,
					build_bem3d_cube_quadpoints);
  test_h2matrix_system("Greenhybrid", Vfull, KMfull, block, bem_slp, V2,
		       bem_dlp, KM2, false, false, 7.0e-2, 7.5e-2);

  /* Nearfield computation on GPU not applicable here yet! */

//  setup_h2matrix_aprx_greenhybrid_ortho_bem3d(bem_slp, Vrb, Vcb, block, m, l,
//      delta, eps_aca, build_bem3d_cube_quadpoints);
//  setup_h2matrix_aprx_greenhybrid_ortho_bem3d(bem_dlp, KMrb, KMcb, block, m, l,
//      delta, eps_aca, build_bem3d_cube_quadpoints);
//  test_h2matrix_system("Greenhybrid ortho", Vfull, KMfull, block, bem_slp, V2,
//      bem_dlp, KM2, false, false,7.0e-2, 7.5e-2);
  del_h2matrix(V2);
  del_h2matrix(KM2);
  del_block(block);
  freemem(root->idx);
  del_cluster(root);
  del_laplace_ocl_bem3d(bem_slp);
  del_laplace_ocl_bem3d(bem_dlp);

  del_amatrix(Vfull);
  del_amatrix(KMfull);
  del_surface3d(gr);
  del_macrosurface3d(mg);
  freemem(idx);

  (void) printf("----------------------------------------\n"
		"  %u matrices and\n"
		"  %u vectors still active\n"
		"  %u errors found\n", getactives_amatrix(),
		getactives_avector(), problems);

  uninit_h2lib();

  return problems;
}

#else
#include <stdio.h>

int
main(int argc, char **argv)
{
  (void) argc;
  (void) argv;

  fprintf(stderr,
	  "Need to activate OpenMP by \'-DUSE_OPENMP\' to use this test\n");

  return 0;
}
#endif

#else
#include <stdio.h>

int
main(int argc, char **argv)
{
  (void) argc;
  (void) argv;

  fprintf(stderr,
	  "Need to activate OpenCL by \'-DUSE_OPENCL\' to use this test\n");

  return 0;
}
#endif
