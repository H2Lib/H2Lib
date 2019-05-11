#include <stdio.h>

#include "parameters.h"
#include "basic.h"
#include "hmatrix.h"
#include "clusterbasis.h"
#include "h2matrix.h"
#include "h2compression.h"
#include "laplacebem2d.h"
#include "matrixnorms.h"

static uint problems = 0;

#define IS_IN_RANGE(a, b, c) (((a) <= (b)) && ((b) <= (c)))

#ifdef USE_FLOAT
static real tolerance = 1.0e-6;
#else
static real tolerance = 5.0e-8;
#endif

int
main(int argc, char **argv)
{
  pamatrix  G;			/* Original matrix */
  pcluster  root;		/* Cluster tree */
  pblock    broot;		/* Block cluster tree */
  pcurve2d  gr;
  pbem2d    bem;		/* Auxiliary date for ie1d problem */
  pclusterbasis rb, cb, rb2, cb2;	/* Cluster basis */
  ph2matrix G2;			/* H^2-matrix approximation */
  pclusteroperator rbw, cbw;	/* Weight matrices for cb */
  pclusteroperator rw, cw, rw2, cw2;	/* Total weights for G2 */
  ptruncmode tm;		/* Truncation strategy */
  pclusterbasis rbnew, cbnew;	/* Truncated cluster basis */
  pclusteroperator rold2new, cold2new;	/* Transformation form old to new basis */
  ph2matrix G3;			/* H^2-matrix recompressed */
  phmatrix  Gh;			/* H-matrix */
  pclusterbasis rbh, cbh;	/* Adaptive cluster bases */
  ph2matrix G4;			/* H^2-matrix from H-matrix */
  pclusterbasis rbf, cbf;	/* Adaptive cluster bases */
  ph2matrix G5;			/* H^2-matrix from dense matrix */
  ph2matrix G6;			/* H^2-matrix from hierarchical compression */
  pavector  x, y;		/* Vectors for testing */
  pstopwatch sw;		/* Measure runtime */
  real      t_run;		/* Runtime */
  real      error, normG;	/* Norm and error estimate */
  size_t    sz;			/* Storage size */
  uint      q;			/* order of quadrature */
  uint      n;			/* Resolution */
  uint      m;			/* Approximation order */
  uint      l;			/* Number of quadrature subdivisions */
  uint      clf;		/* maximal leaf size */
  real      eta;		/* admissibility parameter */
  real      delta;		/* distance of bounding box for green approximation */
  real      eps_aca;		/* accuracy for ACA */
  real      eps;		/* Compression tolerance */

  /* Initialize the library */
  init_h2lib(&argc, &argv);

  sw = new_stopwatch();

  n = 1024;
  q = 2;

  clf = 16;
  eta = 1.0;
  m = 4;
  l = 1;
  delta = 1.0;
  eps_aca = tolerance * 0.1;

  eps = tolerance;

  (void) printf("========================================\n"
		"Creating curve2d circle\n");

  gr = new_circle_curve2d(n, 0.333);

  (void) printf("========================================\n"
		"Creating bem2d structure\n");
  bem = new_slp_laplace_bem2d(gr, q, BASIS_CONSTANT_BEM2D);

  G = 0;
  if (n <= 16384) {
    (void) printf("Creating dense matrix\n");
    G = new_amatrix(n, n);
    sz = getsize_amatrix(G);
    (void) printf("  %.2f KB (%.2f KB/DoF)\n", sz / 1024.0, sz / 1024.0 / n);

    start_stopwatch(sw);
    (void) printf("Filling dense matrix\n");
    bem->nearfield(NULL, NULL, bem, false, G);
    t_run = stop_stopwatch(sw);
    (void) printf("  %.2f seconds\n", t_run);
    normG = norm2_amatrix(G);
  }
  else
    (void) printf("Skipping dense matrix, it would require %.2f GB\n",
		  (real) sizeof(field) * n * n / 1073741824.0);

  (void) printf("========================================\n"
		"Creating cluster tree\n");
  root = build_bem2d_cluster(bem, clf, BASIS_CONSTANT_BEM2D);
  (void) printf("  %d clusters\n", root->desc);

  (void) printf("Creating block tree\n");
  broot = build_strict_block(root, root, &eta, admissible_max_cluster);
  (void) printf("  %d blocks\n", broot->desc);

  start_stopwatch(sw);
  (void) printf("========================================\n"
		"Creating cluster basis\n");
  rb = build_from_cluster_clusterbasis(root);
  cb = build_from_cluster_clusterbasis(root);
  setup_h2matrix_aprx_greenhybrid_bem2d(bem, rb, cb, broot, m, l, delta,
					eps_aca, build_bem2d_rect_quadpoints);
  assemble_bem2d_h2matrix_row_clusterbasis(bem, rb);
  assemble_bem2d_h2matrix_col_clusterbasis(bem, cb);
  sz = getsize_clusterbasis(rb);
  (void) printf("  %.2f KB (%.2f KB/DoF)\n"
		"  Rank sum %d (%.2f KB for coefficients)\n"
		"  Rank branch sum %d\n", sz / 1024.0, sz / 1024.0 / n,
		rb->ktree,
		((size_t) sizeof(field) * rb->ktree +
		 sizeof(avector)) / 1024.0, rb->kbranch);
  sz = getsize_clusterbasis(cb);
  (void) printf("  %.2f KB (%.2f KB/DoF)\n"
		"  Rank sum %d (%.2f KB for coefficients)\n"
		"  Rank branch sum %d\n", sz / 1024.0, sz / 1024.0 / n,
		cb->ktree,
		((size_t) sizeof(field) * cb->ktree +
		 sizeof(avector)) / 1024.0, cb->kbranch);

  (void) printf("Creating H^2-matrix structure\n");
  G2 = build_from_block_h2matrix(broot, rb, cb);
  sz = getsize_h2matrix(G2);
  (void) printf("  %.2f KB (%.2f KB/DoF)\n", sz / 1024.0, sz / 1024.0 / n);

  (void) printf("Filling H^2-matrix structure\n");
  assemble_bem2d_h2matrix(bem, broot, G2);
  t_run = stop_stopwatch(sw);
  (void) printf("  %.2f seconds\n", t_run);

  if (G) {
    x = new_avector(G->cols);
    y = new_avector(G->rows);

    (void) printf("Testing matrix-vector multiplication\n");
    fill_avector(x, 1.0 / G->rows);
    clear_avector(y);
    addeval_h2matrix_avector(1.0, G2, x, y);
    addeval_amatrix_avector(-1.0, G, x, y);
    error = norm2_avector(y);
    (void) printf("  Difference for constant vector: %.4g %s okay\n", error,
		  IS_IN_RANGE(0.0, error,
			      10.0 * tolerance) ? "       " : "   NOT ");
    if (!IS_IN_RANGE(0.0, error, 10.0 * tolerance))
      problems++;

    (void) printf("Rel. spectral error bound by power iteration\n");
    error = norm2diff_amatrix_h2matrix(G2, G) / normG;
    (void) printf("  %.4e                                %s okay\n", error,
		  IS_IN_RANGE(0.0, error,
			      10.0 * tolerance) ? "       " : "   NOT ");
    if (!IS_IN_RANGE(0.0, error, 10.0 * tolerance))
      problems++;

    del_avector(y);
    del_avector(x);
  }

  (void) printf("========================================\n"
		"Computing cluster basis weights\n");

  start_stopwatch(sw);

  rbw = build_from_clusterbasis_clusteroperator(rb);
  cbw = build_from_clusterbasis_clusteroperator(cb);

  weight_clusterbasis_clusteroperator(rb, rbw);
  weight_clusterbasis_clusteroperator(cb, cbw);

  t_run = stop_stopwatch(sw);

  sz = getsize_clusteroperator(rbw);
  (void) printf("  %.2f KB (%.2f KB/DoF)\n"
		"  %.2f seconds\n", sz / 1024.0, sz / 1024.0 / n, t_run);
  sz = getsize_clusteroperator(cbw);
  (void) printf("  %.2f KB (%.2f KB/DoF)\n"
		"  %.2f seconds\n", sz / 1024.0, sz / 1024.0 / n, t_run);

  (void) printf("----------------------------------------\n"
		"Computing local weights\n");

  start_stopwatch(sw);

  tm = new_blockreleucl_truncmode();
  rw = build_from_clusterbasis_clusteroperator(rb);
  cw = build_from_clusterbasis_clusteroperator(cb);

  localweights_h2matrix(G2, rbw, cbw, tm, rw, cw);

  t_run = stop_stopwatch(sw);

  sz = getsize_clusteroperator(rw);
  (void) printf("  %.2f KB (%.2f KB/DoF) for row weights\n", sz / 1024.0,
		sz / 1024.0 / n);

  sz = getsize_clusteroperator(cw);
  (void) printf("  %.2f KB (%.2f KB/DoF) for column weights\n"
		"  %.2f seconds\n", sz / 1024.0, sz / 1024.0 / n, t_run);

  (void) printf("Computing local weights by alternative method\n");

  rw2 = build_from_clusterbasis_clusteroperator(rb);
  cw2 = build_from_clusterbasis_clusteroperator(cb);

  rowweights_h2matrix(G2, rbw, cbw, tm, rw2);
  colweights_h2matrix(G2, rbw, cbw, tm, cw2);

  error = compareweights_clusteroperator(rw, rw2);
  (void) printf("  Relative row weight error %.4e      %s okay\n", error,
		IS_IN_RANGE(0.0, error,
			    10.0 * tolerance) ? "       " : "   NOT ");
  if (!IS_IN_RANGE(0.0, error, 10.0 * tolerance))
    problems++;

  error = compareweights_clusteroperator(cw, cw2);
  (void) printf("  Relative column weight error %.4e   %s okay\n", error,
		IS_IN_RANGE(0.0, error,
			    10.0 * tolerance) ? "       " : "   NOT ");
  if (!IS_IN_RANGE(0.0, error, 10.0 * tolerance))
    problems++;

  del_clusteroperator(cw2);
  del_clusteroperator(rw2);

  (void) printf("----------------------------------------\n"
		"Computing total weights\n");

  start_stopwatch(sw);

  accumulate_clusteroperator(rb, tm, rw);
  accumulate_clusteroperator(cb, tm, cw);

  t_run = stop_stopwatch(sw);

  sz = getsize_clusteroperator(rw);
  (void) printf("  %.2f KB (%.2f KB/DoF) for row weights\n", sz / 1024.0,
		sz / 1024.0 / n);

  sz = getsize_clusteroperator(cw);
  (void) printf("  %.2f KB (%.2f KB/DoF) for column weights\n"
		"  %.2f seconds\n", sz / 1024.0, sz / 1024.0 / n, t_run);

  (void) printf("----------------------------------------\n"
		"Computing truncated cluster basis, eps=%.2g\n", eps);

  start_stopwatch(sw);

  rbnew = clonestructure_clusterbasis(rb);
  rold2new = build_from_clusterbasis_clusteroperator(rb);
  truncate_clusterbasis(rb, rw, 0, tm, eps, rbnew, rold2new);

  cbnew = clonestructure_clusterbasis(cb);
  cold2new = build_from_clusterbasis_clusteroperator(cb);
  truncate_clusterbasis(cb, cw, 0, tm, eps, cbnew, cold2new);

  t_run = stop_stopwatch(sw);

  sz = getsize_clusterbasis(rbnew);
  (void) printf("  %.2f KB (%.2f KB/DoF) for new cluster basis\n"
		"  %.2f seconds\n"
		"  Rank sum %u\n"
		"  Rank branch sum %u\n", sz / 1024.0, sz / 1024.0 / n, t_run,
		rbnew->ktree, rbnew->kbranch);

  sz = getsize_clusterbasis(cbnew);
  (void) printf("  %.2f KB (%.2f KB/DoF) for new cluster basis\n"
		"  %.2f seconds\n"
		"  Rank sum %u\n"
		"  Rank branch sum %u\n", sz / 1024.0, sz / 1024.0 / n, t_run,
		cbnew->ktree, cbnew->kbranch);

  (void) printf("----------------------------------------\n"
		"Projecting to new basis\n");

  start_stopwatch(sw);
  G3 = build_projected_h2matrix(G2, rbnew, rold2new, cbnew, cold2new);
  t_run = stop_stopwatch(sw);

  sz = getsize_h2matrix(G3);
  (void) printf("  %.2f KB (%.2f KB/DoF)\n"
		"  %.2f seconds\n", sz / 1024.0, sz / 1024.0 / n, t_run);

  (void) printf("Rel. spectral error bound by power iteration\n");
  error = norm2diff_h2matrix(G2, G3) / normG;
  (void) printf("  %.4e                                %s okay\n", error,
		IS_IN_RANGE(0.0, error,
			    10.0 * tolerance) ? "       " : "   NOT ");
  if (!IS_IN_RANGE(0.0, error, 10.0 * tolerance))
    problems++;

#ifdef USE_NETCDF
  (void) printf("Writing to \"G3.nc\"\n");
  start_stopwatch(sw);
  write_cdfcomplete_h2matrix(G3, "G3.nc");
  t_run = stop_stopwatch(sw);

  (void) printf("Reading from \"G3.nc\"\n");
  start_stopwatch(sw);
  G4 = read_cdfcomplete_h2matrix("G3.nc");
  t_run = stop_stopwatch(sw);

  (void) printf("Relative spectral error bound by power iteration\n");
  error = norm2diff_h2matrix(G3, G4) / normG;
  (void) printf("  %.4e                                %s okay\n", error,
		IS_IN_RANGE(0.0, error,
			    10.0 * tolerance) ? "       " : "   NOT ");
  if (!IS_IN_RANGE(0.0, error, 10.0 * tolerance))
    problems++;

  del_h2matrix(G4);
#endif

  (void) printf("========================================\n"
		"Creating H-matrix structure\n");

  Gh = build_from_block_hmatrix(broot, m);
  sz = getsize_hmatrix(Gh);
  (void) printf("  %.2f KB (%.2f KB/DoF)\n", sz / 1024.0, sz / 1024.0 / n);

  sz = getfarsize_hmatrix(Gh);
  (void) printf("   far:  %.2f KB (%.2f KB/DoF)\n", sz / 1024.0,
		sz / 1024.0 / n);

  sz = getnearsize_hmatrix(Gh);
  (void) printf("   near: %.2f KB (%.2f KB/DoF)\n", sz / 1024.0,
		sz / 1024.0 / n);

  (void) printf("Filling H-matrix structure\n");
  setup_hmatrix_aprx_greenhybrid_row_bem2d(bem, root, root, broot, m, l,
					   delta, eps_aca,
					   build_bem2d_rect_quadpoints);
  assemble_bem2d_hmatrix(bem, broot, Gh);
  t_run = stop_stopwatch(sw);
  (void) printf("  %.2f seconds\n", t_run);

  (void) printf("----------------------------------------\n"
		"Building adaptive row cluster basis for H-matrix\n");

  start_stopwatch(sw);
  rbh = buildrowbasis_hmatrix(Gh, tm, eps);
  t_run = stop_stopwatch(sw);

  sz = getsize_clusterbasis(rbh);
  (void) printf("  %.2f KB (%.2f KB/DoF) for new cluster basis\n"
		"  %.2f seconds\n"
		"  Rank sum %u\n"
		"  Rank branch sum %u\n", sz / 1024.0, sz / 1024.0 / n, t_run,
		rbh->ktree, rbh->kbranch);

  (void) printf("----------------------------------------\n"
		"Building adaptive column cluster basis for H-matrix\n");

  start_stopwatch(sw);
  cbh = buildcolbasis_hmatrix(Gh, tm, eps);
  t_run = stop_stopwatch(sw);

  sz = getsize_clusterbasis(cbh);
  (void) printf("  %.2f KB (%.2f KB/DoF) for new H^2-matrix\n"
		"  %.2f seconds\n"
		"  Rank sum %u\n"
		"  Rank branch sum %u\n", sz / 1024.0, sz / 1024.0 / n, t_run,
		cbh->ktree, cbh->kbranch);

  (void) printf("----------------------------------------\n"
		"Approximating H-matrix by H^2-matrix\n");

  start_stopwatch(sw);
  G4 = build_projected_hmatrix_h2matrix(Gh, rbh, cbh);
  t_run = stop_stopwatch(sw);

  sz = getsize_h2matrix(G4);
  (void) printf("  %.2f KB (%.2f KB/DoF) for new H^2-matrix\n"
		"  %.2f seconds\n", sz / 1024.0, sz / 1024.0 / n, t_run);

  (void) printf("Rel. spectral error bound by power iteration\n");
  error = norm2diff_hmatrix_h2matrix(G4, Gh) / normG;
  (void) printf("  %.4e                                %s okay\n", error,
		IS_IN_RANGE(0.0, error,
			    10.0 * tolerance) ? "       " : "   NOT ");
  if (!IS_IN_RANGE(0.0, error, 10.0 * tolerance))
    problems++;

  (void) printf("========================================\n");
  (void) printf("Building H^2-matrix with hierarchical compression\n");
  start_stopwatch(sw);
  rb2 = build_from_cluster_clusterbasis(root);
  cb2 = build_from_cluster_clusterbasis(root);
  G6 = build_from_block_h2matrix(broot, rb2, cb2);
  setup_h2matrix_recomp_bem2d(bem, true, eps * 0.25);
  assemblehiercomp_bem2d_h2matrix(bem, broot, G6);
  t_run = stop_stopwatch(sw);

  sz = getsize_clusterbasis(G6->rb);
  (void) printf("  %.2f KB (%.2f KB/DoF) for row cluster basis\n"
		"  Rank sum %u\n"
		"  Rank branch sum %u\n", sz / 1024.0, sz / 1024.0 / n,
		G6->rb->ktree, G6->rb->kbranch);

  sz = getsize_clusterbasis(G6->cb);
  (void) printf("  %.2f KB (%.2f KB/DoF) for column cluster basis\n"
		"  Rank sum %u\n"
		"  Rank branch sum %u\n", sz / 1024.0, sz / 1024.0 / n,
		G6->cb->ktree, G6->cb->kbranch);

  sz = getsize_h2matrix(G6);
  (void) printf("  %.2f KB (%.2f KB/DoF) for new H^2-matrix\n"
		"  %.2f seconds\n", sz / 1024.0, sz / 1024.0 / n, t_run);

  (void) printf("Rel. spectral error bound by power iteration\n");
  error = norm2diff_hmatrix_h2matrix(G6, Gh) / normG;
  (void) printf("  %.4e                                %s okay\n", error,
		IS_IN_RANGE(0.0, error,
			    10.0 * tolerance) ? "       " : "   NOT ");
  if (!IS_IN_RANGE(0.0, error, 10.0 * tolerance))
    problems++;

  G5 = 0;
  if (G) {
    (void) printf("========================================\n"
		  "Building adaptive row cluster basis for dense matrix\n");

    start_stopwatch(sw);
    rbf = buildrowbasis_amatrix(G, broot, tm, eps);
    t_run = stop_stopwatch(sw);

    sz = getsize_clusterbasis(rbf);
    (void) printf("  %.2f KB (%.2f KB/DoF) for new cluster basis\n"
		  "  %.2f seconds\n"
		  "  Rank sum %u\n"
		  "  Rank branch sum %u\n", sz / 1024.0, sz / 1024.0 / n,
		  t_run, rbf->ktree, rbf->kbranch);

    (void) printf("----------------------------------------\n"
		  "Building adaptive column cluster basis for dense matrix\n");

    start_stopwatch(sw);
    cbf = buildcolbasis_amatrix(G, broot, tm, eps);
    t_run = stop_stopwatch(sw);

    sz = getsize_clusterbasis(cbf);
    (void) printf("  %.2f KB (%.2f KB/DoF) for new cluster basis\n"
		  "  %.2f seconds\n"
		  "  Rank sum %u\n"
		  "  Rank branch sum %u\n", sz / 1024.0, sz / 1024.0 / n,
		  t_run, cbf->ktree, cbf->kbranch);

    (void) printf("----------------------------------------\n"
		  "Approximating dense matrix by H^2-matrix\n");

    start_stopwatch(sw);
    G5 = build_projected_amatrix_h2matrix(G, broot, rbf, cbf);
    t_run = stop_stopwatch(sw);

    sz = getsize_h2matrix(G5);
    (void) printf("  %.2f KB (%.2f KB/DoF) for new H^2-matrix\n"
		  "  %.2f seconds\n", sz / 1024.0, sz / 1024.0 / n, t_run);

    (void) printf("Rel. spectral error bound by power iteration\n");
    error = norm2diff_amatrix_h2matrix(G5, G) / normG;
    (void) printf("  %.4e                                %s okay\n", error,
		  IS_IN_RANGE(0.0, error,
			      10.0 * tolerance) ? "       " : "   NOT ");
    if (!IS_IN_RANGE(0.0, error, 10.0 * tolerance))
      problems++;
  }

  (void) printf("========================================\n" "Cleaning up\n");

  del_h2matrix(G6);

  if (G5)
    del_h2matrix(G5);

  del_h2matrix(G4);
  del_hmatrix(Gh);
  del_h2matrix(G3);
  del_clusteroperator(rold2new);
  del_clusteroperator(cold2new);
  del_truncmode(tm);
  del_clusteroperator(cw);
  del_clusteroperator(rw);
  del_clusteroperator(cbw);
  del_clusteroperator(rbw);
  del_h2matrix(G2);
  del_block(broot);
  freemem(root->idx);
  del_cluster(root);
  del_bem2d(bem);
  del_curve2d(gr);
  del_stopwatch(sw);

  if (G)
    del_amatrix(G);

  (void) printf("  %u cluster bases,\n"
		"  %u matrices and\n"
		"  %u vectors still active\n"
		"  %u errors found\n", getactives_clusterbasis(),
		getactives_amatrix(), getactives_avector(), problems);

  return problems;
}
