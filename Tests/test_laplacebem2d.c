#include "basic.h"
#include "krylov.h"
#include "krylovsolvers.h"
#include "laplacebem2d.h"
#include "matrixnorms.h"

static uint problems = 0;

#define IS_IN_RANGE(a, b, c) (((a) <= (b)) && ((b) <= (c)))

/* Compute the L_2-error for constant basis functions */
static    real
L2gamma_c_diff_norm2(pbem2d bem, pavector x, boundary_func2d rhs)
{
  pccurve2d gr = bem->gr;
  const     real(*gr_x)[2] = (const real(*)[2]) gr->x;
  const     uint(*gr_e)[2] = (const uint(*)[2]) gr->e;
  const     real(*gr_n)[2] = (const real(*)[2]) gr->n;
  const preal gr_g = (const preal) gr->g;
  uint      n = x->dim;
  uint      nq = bem->sq->n_single;
  real     *xx = bem->sq->x_single;
  real     *ww = bem->sq->w_single;

  const real *A, *B, *N;
  uint      e, q;
  real      norm, sum, X[2], tx, Ax, Bx;

  assert(n == gr->edges);

  norm = 0.0;
  for (e = 0; e < n; ++e) {
    A = gr_x[gr_e[e][0]];
    B = gr_x[gr_e[e][1]];
    N = gr_n[e];

    sum = 0.0;

    for (q = 0; q < nq; ++q) {
      tx = xx[q];
      Ax = 1.0 - tx;
      Bx = tx;

      X[0] = A[0] * Ax + B[0] * Bx;
      X[1] = A[1] * Ax + B[1] * Bx;

      sum += ww[q] * ABSSQR(rhs(X, N) - x->v[e]);
    }

    norm += sum * gr_g[e];
  }

  norm = REAL_SQRT(norm);

  return norm;
}

/* Simple convenience wrapper for conjugate gradient solver */
static void
solve_cg_bem2d(matrixtype type, void *A, pavector b, pavector x,
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
		    pcamatrix KMfull, pblock block, pbem2d bem_slp,
		    phmatrix V, pbem2d bem_dlp, phmatrix KM, bool exterior)
{
  pavector  x, b;
  real      errorV, errorKM, error_solve, eps_solve;
  uint      steps;

  eps_solve = 1.0e-12;
  steps = 1000;

  printf("Testing: %s Hmatrix %s\n"
	 "====================================\n\n",
	 (exterior == true ? "exterior" : "interior"), apprxtype);

  assemble_bem2d_hmatrix(bem_slp, block, V);
  assemble_bem2d_hmatrix(bem_dlp, block, KM);

  errorV = norm2diff_amatrix_hmatrix(V, Vfull) / norm2_amatrix(Vfull);
  printf("rel. error V       : %.5e\n", errorV);
  errorKM = norm2diff_amatrix_hmatrix(KM, KMfull) / norm2_amatrix(KMfull);
  printf("rel. error K%c0.5*M : %.5e\n", (exterior == true ? '-' : '+'),
	 errorKM);

  x = new_avector(Vfull->rows);
  b = new_avector(KMfull->cols);

  printf("Solving Dirichlet problem:\n");

  projectl2_bem2d_const_avector(bem_dlp,
				eval_dirichlet_quadratic_laplacebem2d, x);
  clear_avector(b);
  addeval_hmatrix_avector(1.0, KM, x, b);

  solve_cg_bem2d(HMATRIX, V, b, x, eps_solve, steps);
  if (exterior == true) {
    scale_avector(-1.0, x);
  }

  error_solve = L2gamma_c_diff_norm2(bem_slp, x,
				     eval_neumann_quadratic_laplacebem2d);

  clear_avector(x);
  error_solve = error_solve
    / L2gamma_c_diff_norm2(bem_slp, x, eval_neumann_quadratic_laplacebem2d);

  printf("rel. error neumann : %.5e       %s\n", error_solve,
	 (IS_IN_RANGE(3.0e-3, error_solve, 4.0e-3) ? "    okay" :
	  "NOT okay"));

  if (!IS_IN_RANGE(3.0e-3, error_solve, 4.0e-3))
    problems++;

  printf("\n");

  del_avector(x);
  del_avector(b);
}

static void
test_h2matrix_system(const char *apprxtype, pcamatrix Vfull,
		     pcamatrix KMfull, pblock block, pbem2d bem_slp,
		     ph2matrix V, pbem2d bem_dlp, ph2matrix KM, bool exterior)
{
  pavector  x, b;
  real      errorV, errorKM, error_solve, eps_solve;
  uint      steps;

  eps_solve = 1.0e-12;
  steps = 1000;

  printf("Testing: %s H2matrix %s\n"
	 "====================================\n\n",
	 (exterior == true ? "exterior" : "interior"), apprxtype);

  assemble_bem2d_h2matrix_row_clusterbasis(bem_slp, V->rb);
  assemble_bem2d_h2matrix_col_clusterbasis(bem_slp, V->cb);
  assemble_bem2d_h2matrix(bem_slp, block, V);

  assemble_bem2d_h2matrix_row_clusterbasis(bem_dlp, KM->rb);
  assemble_bem2d_h2matrix_col_clusterbasis(bem_dlp, KM->cb);
  assemble_bem2d_h2matrix(bem_dlp, block, KM);

  errorV = norm2diff_amatrix_h2matrix(V, Vfull) / norm2_amatrix(Vfull);
  printf("rel. error V       : %.5e\n", errorV);
  errorKM = norm2diff_amatrix_h2matrix(KM, KMfull) / norm2_amatrix(KMfull);
  printf("rel. error K%c0.5*M : %.5e\n", (exterior == true ? '-' : '+'),
	 errorKM);

  x = new_avector(Vfull->rows);
  b = new_avector(KMfull->cols);

  printf("Solving Dirichlet problem:\n");

  projectl2_bem2d_const_avector(bem_dlp,
				eval_dirichlet_quadratic_laplacebem2d, x);
  clear_avector(b);
  addeval_h2matrix_avector(1.0, KM, x, b);

  solve_cg_bem2d(H2MATRIX, V, b, x, eps_solve, steps);
  if (exterior == true) {
    scale_avector(-1.0, x);
  }

  error_solve = L2gamma_c_diff_norm2(bem_slp, x,
				     eval_neumann_quadratic_laplacebem2d);

  clear_avector(x);
  error_solve = error_solve
    / L2gamma_c_diff_norm2(bem_slp, x, eval_neumann_quadratic_laplacebem2d);

  printf("rel. error neumann : %.5e       %s\n", error_solve,
	 (IS_IN_RANGE(3.0e-3, error_solve, 4.0e-3) ? "    okay" :
	  "NOT okay"));

  if (!IS_IN_RANGE(3.0e-3, error_solve, 4.0e-3))
    problems++;

  printf("\n");

  del_avector(x);
  del_avector(b);
}

int
main(int argc, char **argv)
{
  pcurve2d  gr;
  pamatrix  Vfull, KMfull;
  pbem2d    bem_slp, bem_dlp;
  pcluster  root;
  pblock block;
  phmatrix  V, KM;
  pclusterbasis Vrb, Vcb, KMrb, KMcb;
  ph2matrix V2, KM2;
  uint      n, q, clf, m, l;
  real      eta, delta, eps_aca;

  init_h2lib(&argc, &argv);

  n = 579;
  q = 2;
  clf = 16;
  eta = 1.0;

  gr = new_circle_curve2d(n, 0.333);
  bem_slp = new_slp_laplace_bem2d(gr, q, BASIS_CONSTANT_BEM2D);
  bem_dlp = new_dlp_laplace_bem2d(gr, q, BASIS_CONSTANT_BEM2D,
				  BASIS_CONSTANT_BEM2D, 0.5);
  root = build_bem2d_cluster(bem_slp, clf, BASIS_CONSTANT_BEM2D);
  block = build_nonstrict_block(root, root, &eta, admissible_max_cluster);

  Vfull = new_amatrix(n, n);
  KMfull = new_amatrix(n, n);
  bem_slp->nearfield(NULL, NULL, bem_slp, false, Vfull);
  bem_dlp->nearfield(NULL, NULL, bem_dlp, false, KMfull);

  V = build_from_block_hmatrix(block, 0);
  KM = build_from_block_hmatrix(block, 0);

  printf("----------------------------------------\n");
  printf("Testing inner Boundary integral equations:\n");
  printf("----------------------------------------\n\n");

  /*
   * Test Interpolation
   */

  m = 9;

  setup_hmatrix_aprx_inter_row_bem2d(bem_slp, root, root, block, m);
  setup_hmatrix_aprx_inter_row_bem2d(bem_dlp, root, root, block, m);
  test_hmatrix_system("Interpolation row", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, false);

  setup_hmatrix_aprx_inter_col_bem2d(bem_slp, root, root, block, m);
  setup_hmatrix_aprx_inter_col_bem2d(bem_dlp, root, root, block, m);
  test_hmatrix_system("Interpolation column", Vfull, KMfull, block, bem_slp,
		      V, bem_dlp, KM, false);

  setup_hmatrix_aprx_inter_mixed_bem2d(bem_slp, root, root, block, m);
  setup_hmatrix_aprx_inter_mixed_bem2d(bem_dlp, root, root, block, m);
  test_hmatrix_system("Interpolation mixed", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, false);

  /*
   * Test Green
   */

  m = 10;
  l = 2;
  delta = 0.5;

  setup_hmatrix_aprx_green_row_bem2d(bem_slp, root, root, block, m, l, delta,
				     build_bem2d_rect_quadpoints);
  setup_hmatrix_aprx_green_row_bem2d(bem_dlp, root, root, block, m, l, delta,
				     build_bem2d_rect_quadpoints);
  test_hmatrix_system("Green row", Vfull, KMfull, block, bem_slp, V, bem_dlp,
		      KM, false);

  setup_hmatrix_aprx_green_col_bem2d(bem_slp, root, root, block, m, l, delta,
				     build_bem2d_rect_quadpoints);
  setup_hmatrix_aprx_green_col_bem2d(bem_dlp, root, root, block, m, l, delta,
				     build_bem2d_rect_quadpoints);
  test_hmatrix_system("Green column", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, false);

  setup_hmatrix_aprx_green_mixed_bem2d(bem_slp, root, root, block, m, l,
				       delta, build_bem2d_rect_quadpoints);
  setup_hmatrix_aprx_green_mixed_bem2d(bem_dlp, root, root, block, m, l,
				       delta, build_bem2d_rect_quadpoints);
  test_hmatrix_system("Green mixed", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, false);

  /*
   * Test Greenhybrid
   */

  m = 3;
  l = 1;
  delta = 1.0;
  eps_aca = 1.0e-7;

  setup_hmatrix_aprx_greenhybrid_row_bem2d(bem_slp, root, root, block, m, l,
					   delta, eps_aca,
					   build_bem2d_rect_quadpoints);
  setup_hmatrix_aprx_greenhybrid_row_bem2d(bem_dlp, root, root, block, m, l,
					   delta, eps_aca,
					   build_bem2d_rect_quadpoints);
  test_hmatrix_system("Greenhybrid row", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, false);

  setup_hmatrix_aprx_greenhybrid_col_bem2d(bem_slp, root, root, block, m, l,
					   delta, eps_aca,
					   build_bem2d_rect_quadpoints);
  setup_hmatrix_aprx_greenhybrid_col_bem2d(bem_dlp, root, root, block, m, l,
					   delta, eps_aca,
					   build_bem2d_rect_quadpoints);
  test_hmatrix_system("Greenhybrid column", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, false);

  setup_hmatrix_aprx_greenhybrid_mixed_bem2d(bem_slp, root, root, block, m, l,
					     delta, eps_aca,
					     build_bem2d_rect_quadpoints);
  setup_hmatrix_aprx_greenhybrid_mixed_bem2d(bem_dlp, root, root, block, m, l,
					     delta, eps_aca,
					     build_bem2d_rect_quadpoints);
  test_hmatrix_system("Greenhybrid mixed", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, false);

  /*
   * Test ACA / PACA / HCA
   */

  m = 4;
  eps_aca = 1.0e-6;

  setup_hmatrix_aprx_aca_bem2d(bem_slp, root, root, block, eps_aca);
  setup_hmatrix_aprx_aca_bem2d(bem_dlp, root, root, block, eps_aca);
  test_hmatrix_system("ACA full pivoting", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, false);

  setup_hmatrix_aprx_paca_bem2d(bem_slp, root, root, block, eps_aca);
  setup_hmatrix_aprx_paca_bem2d(bem_dlp, root, root, block, eps_aca);
  test_hmatrix_system("ACA partial pivoting", Vfull, KMfull, block, bem_slp,
		      V, bem_dlp, KM, false);

  setup_hmatrix_aprx_hca_bem2d(bem_slp, root, root, block, m, eps_aca);
  setup_hmatrix_aprx_hca_bem2d(bem_dlp, root, root, block, m, eps_aca);
  test_hmatrix_system("HCA2", Vfull, KMfull, block, bem_slp, V, bem_dlp, KM,
		      false);

  /*
   * H2-matrix
   */

  del_hmatrix(V);
  del_hmatrix(KM);
  del_block(block);

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

  m = 8;

  setup_h2matrix_aprx_inter_bem2d(bem_slp, Vrb, Vcb, block, m);
  setup_h2matrix_aprx_inter_bem2d(bem_dlp, KMrb, KMcb, block, m);
  test_h2matrix_system("Interpolation", Vfull, KMfull, block, bem_slp, V2,
		       bem_dlp, KM2, false);

  /*
   * Test Greenhybrid
   */

  m = 3;
  l = 1;
  delta = 1.0;
  eps_aca = 1.0e-7;

  setup_h2matrix_aprx_greenhybrid_bem2d(bem_slp, Vrb, Vcb, block, m, l, delta,
					eps_aca, build_bem2d_rect_quadpoints);
  setup_h2matrix_aprx_greenhybrid_bem2d(bem_dlp, KMrb, KMcb, block, m, l,
					delta, eps_aca,
					build_bem2d_rect_quadpoints);
  test_h2matrix_system("Greenhybrid", Vfull, KMfull, block, bem_slp, V2,
		       bem_dlp, KM2, false);

  setup_h2matrix_aprx_greenhybrid_ortho_bem2d(bem_slp, Vrb, Vcb, block, m, l,
					      delta, eps_aca,
					      build_bem2d_rect_quadpoints);
  setup_h2matrix_aprx_greenhybrid_ortho_bem2d(bem_dlp, KMrb, KMcb, block, m,
					      l, delta, eps_aca,
					      build_bem2d_rect_quadpoints);
  test_h2matrix_system("Greenhybrid ortho", Vfull, KMfull, block, bem_slp, V2,
		       bem_dlp, KM2, false);

  del_h2matrix(V2);
  del_h2matrix(KM2);
  del_block(block);
  del_bem2d(bem_dlp);

  printf("----------------------------------------\n");
  printf("Testing outer Boundary integral equations:\n");
  printf("----------------------------------------\n\n");

  bem_dlp = new_dlp_laplace_bem2d(gr, q, BASIS_CONSTANT_BEM2D,
				  BASIS_CONSTANT_BEM2D, -0.5);
  block = build_nonstrict_block(root, root, &eta, admissible_max_cluster);
  bem_dlp->nearfield(NULL, NULL, bem_dlp, false, KMfull);

  V = build_from_block_hmatrix(block, 0);
  KM = build_from_block_hmatrix(block, 0);

  /*
   * Test Interpolation
   */

  m = 9;

  setup_hmatrix_aprx_inter_row_bem2d(bem_slp, root, root, block, m);
  setup_hmatrix_aprx_inter_row_bem2d(bem_dlp, root, root, block, m);
  test_hmatrix_system("Interpolation row", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, true);

  setup_hmatrix_aprx_inter_col_bem2d(bem_slp, root, root, block, m);
  setup_hmatrix_aprx_inter_col_bem2d(bem_dlp, root, root, block, m);
  test_hmatrix_system("Interpolation column", Vfull, KMfull, block, bem_slp,
		      V, bem_dlp, KM, true);

  setup_hmatrix_aprx_inter_mixed_bem2d(bem_slp, root, root, block, m);
  setup_hmatrix_aprx_inter_mixed_bem2d(bem_dlp, root, root, block, m);
  test_hmatrix_system("Interpolation mixed", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, true);

  /*
   * Test Green
   */

  m = 10;
  l = 2;
  delta = 0.5;

  setup_hmatrix_aprx_green_row_bem2d(bem_slp, root, root, block, m, l, delta,
				     build_bem2d_rect_quadpoints);
  setup_hmatrix_aprx_green_row_bem2d(bem_dlp, root, root, block, m, l, delta,
				     build_bem2d_rect_quadpoints);
  test_hmatrix_system("Green row", Vfull, KMfull, block, bem_slp, V, bem_dlp,
		      KM, true);

  setup_hmatrix_aprx_green_col_bem2d(bem_slp, root, root, block, m, l, delta,
				     build_bem2d_rect_quadpoints);
  setup_hmatrix_aprx_green_col_bem2d(bem_dlp, root, root, block, m, l, delta,
				     build_bem2d_rect_quadpoints);
  test_hmatrix_system("Green column", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, true);

  setup_hmatrix_aprx_green_mixed_bem2d(bem_slp, root, root, block, m, l,
				       delta, build_bem2d_rect_quadpoints);
  setup_hmatrix_aprx_green_mixed_bem2d(bem_dlp, root, root, block, m, l,
				       delta, build_bem2d_rect_quadpoints);
  test_hmatrix_system("Green mixed", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, true);

  /*
   * Test Greenhybrid
   */

  m = 3;
  l = 1;
  delta = 1.0;
  eps_aca = 1.0e-7;

  setup_hmatrix_aprx_greenhybrid_row_bem2d(bem_slp, root, root, block, m, l,
					   delta, eps_aca,
					   build_bem2d_rect_quadpoints);
  setup_hmatrix_aprx_greenhybrid_row_bem2d(bem_dlp, root, root, block, m, l,
					   delta, eps_aca,
					   build_bem2d_rect_quadpoints);
  test_hmatrix_system("Greenhybrid row", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, true);

  setup_hmatrix_aprx_greenhybrid_col_bem2d(bem_slp, root, root, block, m, l,
					   delta, eps_aca,
					   build_bem2d_rect_quadpoints);
  setup_hmatrix_aprx_greenhybrid_col_bem2d(bem_dlp, root, root, block, m, l,
					   delta, eps_aca,
					   build_bem2d_rect_quadpoints);
  test_hmatrix_system("Greenhybrid column", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, true);

  setup_hmatrix_aprx_greenhybrid_mixed_bem2d(bem_slp, root, root, block, m, l,
					     delta, eps_aca,
					     build_bem2d_rect_quadpoints);
  setup_hmatrix_aprx_greenhybrid_mixed_bem2d(bem_dlp, root, root, block, m, l,
					     delta, eps_aca,
					     build_bem2d_rect_quadpoints);
  test_hmatrix_system("Greenhybrid mixed", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, true);

  /*
   * Test ACA / PACA / HCA
   */

  m = 4;
  eps_aca = 1.0e-6;

  setup_hmatrix_aprx_aca_bem2d(bem_slp, root, root, block, eps_aca);
  setup_hmatrix_aprx_aca_bem2d(bem_dlp, root, root, block, eps_aca);
  test_hmatrix_system("ACA full pivoting", Vfull, KMfull, block, bem_slp, V,
		      bem_dlp, KM, true);

  setup_hmatrix_aprx_paca_bem2d(bem_slp, root, root, block, eps_aca);
  setup_hmatrix_aprx_paca_bem2d(bem_dlp, root, root, block, eps_aca);
  test_hmatrix_system("ACA partial pivoting", Vfull, KMfull, block, bem_slp,
		      V, bem_dlp, KM, true);

  setup_hmatrix_aprx_hca_bem2d(bem_slp, root, root, block, m, eps_aca);
  setup_hmatrix_aprx_hca_bem2d(bem_dlp, root, root, block, m, eps_aca);
  test_hmatrix_system("HCA2", Vfull, KMfull, block, bem_slp, V, bem_dlp, KM,
		      true);

  /*
   * H2-matrix
   */

  del_hmatrix(V);
  del_hmatrix(KM);
  del_block(block);

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

  m = 8;

  setup_h2matrix_aprx_inter_bem2d(bem_slp, Vrb, Vcb, block, m);
  setup_h2matrix_aprx_inter_bem2d(bem_dlp, KMrb, KMcb, block, m);
  test_h2matrix_system("Interpolation", Vfull, KMfull, block, bem_slp, V2,
		       bem_dlp, KM2, true);

  /*
   * Test Greenhybrid
   */

  m = 3;
  l = 1;
  delta = 1.0;
  eps_aca = 1.0e-7;

  setup_h2matrix_aprx_greenhybrid_bem2d(bem_slp, Vrb, Vcb, block, m, l, delta,
					eps_aca, build_bem2d_rect_quadpoints);
  setup_h2matrix_aprx_greenhybrid_bem2d(bem_dlp, KMrb, KMcb, block, m, l,
					delta, eps_aca,
					build_bem2d_rect_quadpoints);
  test_h2matrix_system("Greenhybrid", Vfull, KMfull, block, bem_slp, V2,
		       bem_dlp, KM2, true);

  setup_h2matrix_aprx_greenhybrid_ortho_bem2d(bem_slp, Vrb, Vcb, block, m, l,
					      delta, eps_aca,
					      build_bem2d_rect_quadpoints);
  setup_h2matrix_aprx_greenhybrid_ortho_bem2d(bem_dlp, KMrb, KMcb, block, m,
					      l, delta, eps_aca,
					      build_bem2d_rect_quadpoints);
  test_h2matrix_system("Greenhybrid Ortho", Vfull, KMfull, block, bem_slp, V2,
		       bem_dlp, KM2, true);

  del_h2matrix(V2);
  del_h2matrix(KM2);
  del_block(block);
  freemem(root->idx);
  del_cluster(root);
  del_bem2d(bem_slp);
  del_bem2d(bem_dlp);
  del_amatrix(Vfull);
  del_amatrix(KMfull);
  del_curve2d(gr);

  (void) printf("----------------------------------------\n"
		"  %u matrices and\n"
		"  %u vectors still active\n"
		"  %u errors found\n", getactives_amatrix(),
		getactives_avector(), problems);

  uninit_h2lib();

  return problems;
}
