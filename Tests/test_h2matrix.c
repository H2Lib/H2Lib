#include <stdio.h>
#include "settings.h"
#include "hmatrix.h"
#include "harith.h"
#include "h2matrix.h"
#include "h2arith.h"
#include "truncation.h"

#include "laplacebem2d.h"

#ifdef USE_CAIRO
#include "cairo/cairo.h"
#endif

static uint problems = 0;

#define IS_IN_RANGE(a, b, c) (((a) < (b)) && ((b) < (c)))

#ifdef USE_COMPLEX
static field alpha = 1.0 + 1.0 * I;
#else
static field alpha = 1.0;
#endif

#ifdef USE_FLOAT
static real tolerance = 1.0e-6;
#else
static real tolerance = 1.0e-12;
#endif

int
main()
{
  ph2matrix h2, h2copy, L, R;
  pclusterbasis rb, cb, rbcopy, cbcopy, rblow, cblow, rbup, cbup;
  pclusteroperator rwf, cwf, rwflow, cwflow, rwfup, cwfup, rwfh2, cwfh2;
  ptruncmode tm;

  pavector  x, b;
  uint      n;
  real      error;
  pcurve2d  gr2;
  pbem2d    bem2;
  pcluster  root2;
  pblock    block2;
  uint      clf, m;
  real      tol, eta, delta, eps_aca;

  n = 579;
  tol = tolerance;

  clf = 16;
  eta = 1.0;
  m = 4;
  delta = 1.0;
  eps_aca = tolerance;

  gr2 = new_circle_curve2d(n, 0.333);
  bem2 = new_slp_laplace_bem2d(gr2, 2, BASIS_CONSTANT_BEM2D);
  root2 = build_bem2d_cluster(bem2, clf, BASIS_CONSTANT_BEM2D);
  block2 = build_strict_block(root2, root2, &eta, admissible_max_cluster);

  rb = build_from_cluster_clusterbasis(root2);
  cb = build_from_cluster_clusterbasis(root2);
  setup_h2matrix_aprx_greenhybrid_bem2d(bem2, rb, cb, block2, m, 1, delta,
					eps_aca, build_bem2d_rect_quadpoints);

  (void) printf("----------------------------------------\n"
		"Check %u x %u Cholesky factorization\n", n, n);

  (void) printf("Creating laplacebem2d SLP matrix\n");
  assemble_bem2d_h2matrix_row_clusterbasis(bem2, rb);
  assemble_bem2d_h2matrix_col_clusterbasis(bem2, cb);
  h2 = build_from_block_h2matrix(block2, rb, cb);
  assemble_bem2d_h2matrix(bem2, block2, h2);

  (void) printf("Creating random solution and right-hand side\n");
  x = new_avector(n);
  random_avector(x);
  b = new_avector(n);
  clear_avector(b);
  mvm_h2matrix_avector(alpha, false, h2, x, b);

  (void) printf("Copying matrix\n");

  rbcopy = clone_clusterbasis(h2->rb);
  cbcopy = clone_clusterbasis(h2->cb);
  h2copy = clone_h2matrix(h2, rbcopy, cbcopy);

  (void) printf("Computing Cholesky factorization\n");

  tm = new_releucl_truncmode();
  rblow = build_from_cluster_clusterbasis(root2);
  cblow = build_from_cluster_clusterbasis(root2);
  L = build_from_block_lower_h2matrix(block2, rblow, cblow);

  rwflow = prepare_row_clusteroperator(L->rb, L->cb, tm);
  cwflow = prepare_col_clusteroperator(L->rb, L->cb, tm);
  rwf = NULL;
  cwf = NULL;
  init_cholesky_h2matrix(h2, &rwf, &cwf, tm);
  choldecomp_h2matrix(h2, rwf, cwf, L, rwflow, cwflow, tm, tol);

  (void) printf("Solving\n");

  cholsolve_h2matrix_avector(L, b);

  add_avector(-alpha, x, b);
  error = norm2_avector(b) / norm2_avector(x);
  (void) printf("  Accuracy %g, %sokay\n", error,
		IS_IN_RANGE(0.0, error, 25.0 * tol) ? "" : "    NOT ");
  if (!IS_IN_RANGE(0.0, error, 25.0 * tol))
    problems++;

  rwfh2 = prepare_row_clusteroperator(h2copy->rb, h2copy->cb, tm);
  cwfh2 = prepare_col_clusteroperator(h2copy->rb, h2copy->cb, tm);

  (void) printf("Checking factorization\n");
  error = norm2_h2matrix(h2copy);
  addmul_h2matrix(-1.0, L, true, L, h2copy, rwfh2, cwfh2, tm, tol);
  error = norm2_h2matrix(h2copy) / error;
  (void) printf("  Accuracy %g, %sokay\n", error,
		IS_IN_RANGE(0.0, error, 25.0 * tol) ? "" : "    NOT ");
  if (!IS_IN_RANGE(0.0, error, 25.0 * tol))
    problems++;

  del_h2matrix(h2copy);
  del_h2matrix(h2);
  del_h2matrix(L);
  del_avector(b);
  del_avector(x);
  del_truncmode(tm);

  del_clusteroperator(rwflow);
  del_clusteroperator(cwflow);
  del_clusteroperator(rwf);
  del_clusteroperator(cwf);
  del_clusteroperator(rwfh2);
  del_clusteroperator(cwfh2);

  rb = build_from_cluster_clusterbasis(root2);
  cb = build_from_cluster_clusterbasis(root2);
  setup_h2matrix_aprx_greenhybrid_bem2d(bem2, rb, cb, block2, m, 1, delta,
					eps_aca, build_bem2d_rect_quadpoints);

  (void) printf("----------------------------------------\n"
		"Check %u x %u LR factorization\n", n, n);

  (void) printf("Creating laplacebem2d SLP matrix\n");
  assemble_bem2d_h2matrix_row_clusterbasis(bem2, rb);
  assemble_bem2d_h2matrix_col_clusterbasis(bem2, cb);
  h2 = build_from_block_h2matrix(block2, rb, cb);
  assemble_bem2d_h2matrix(bem2, block2, h2);

  (void) printf("Creating random solution and right-hand side\n");
  x = new_avector(n);
  random_avector(x);
  b = new_avector(n);
  clear_avector(b);
  mvm_h2matrix_avector(alpha, false, h2, x, b);

  (void) printf("Copying matrix\n");
  rbcopy = clone_clusterbasis(h2->rb);
  cbcopy = clone_clusterbasis(h2->cb);
  h2copy = clone_h2matrix(h2, rbcopy, cbcopy);

  (void) printf("Computing LR factorization\n");
  rblow = build_from_cluster_clusterbasis(root2);
  cblow = build_from_cluster_clusterbasis(root2);
  L = build_from_block_lower_h2matrix(block2, rblow, cblow);

  rbup = build_from_cluster_clusterbasis(root2);
  cbup = build_from_cluster_clusterbasis(root2);
  R = build_from_block_upper_h2matrix(block2, rbup, cbup);

  tm = new_releucl_truncmode();
  rwf = prepare_row_clusteroperator(h2->rb, h2->cb, tm);
  cwf = prepare_col_clusteroperator(h2->rb, h2->cb, tm);
  rwflow = prepare_row_clusteroperator(L->rb, L->cb, tm);
  cwflow = prepare_col_clusteroperator(L->rb, L->cb, tm);
  rwfup = prepare_row_clusteroperator(R->rb, R->cb, tm);
  cwfup = prepare_col_clusteroperator(R->rb, R->cb, tm);

  lrdecomp_h2matrix(h2, rwf, cwf, L, rwflow, cwflow, R, rwfup, cwfup, tm,
		    tol);

  (void) printf("Solving\n");
  lrsolve_h2matrix_avector(L, R, b);

  add_avector(-alpha, x, b);
  error = norm2_avector(b) / norm2_avector(x);
  (void) printf("  Accuracy %g, %sokay\n", error,
		IS_IN_RANGE(0.0, error, 25.0 * tol) ? "" : "    NOT ");
  if (!IS_IN_RANGE(0.0, error, 25.0 * tol))
    problems++;

  rwfh2 = prepare_row_clusteroperator(h2copy->rb, h2copy->cb, tm);
  cwfh2 = prepare_col_clusteroperator(h2copy->rb, h2copy->cb, tm);

  (void) printf("Checking factorization\n");
  error = norm2_h2matrix(h2copy);
  addmul_h2matrix(-1.0, L, false, R, h2copy, rwfh2, cwfh2, tm, tol);
  error = norm2_h2matrix(h2copy) / error;
  (void) printf("  Accuracy %g, %sokay\n", error,
		IS_IN_RANGE(0.0, error, 25.0 * tol) ? "" : "    NOT ");
  if (!IS_IN_RANGE(0.0, error, 25.0 * tol))
    problems++;

  /* Final clean-up */
  (void) printf("Cleaning up\n");

  del_h2matrix(h2);
  del_h2matrix(h2copy);
  del_h2matrix(L);
  del_h2matrix(R);
  del_clusteroperator(rwf);
  del_clusteroperator(cwf);
  del_clusteroperator(rwflow);
  del_clusteroperator(cwflow);
  del_clusteroperator(rwfup);
  del_clusteroperator(cwfup);
  del_clusteroperator(rwfh2);
  del_clusteroperator(cwfh2);

  del_avector(b);
  del_avector(x);
  del_truncmode(tm);

  freemem(root2->idx);
  del_bem2d(bem2);
  del_block(block2);
  del_cluster(root2);
  del_curve2d(gr2);

  (void) printf("----------------------------------------\n"
		"  %u matrices and\n"
		"  %u vectors still active\n"
		"  %u errors found\n", getactives_amatrix(),
		getactives_avector(), problems);

  return problems;
}
