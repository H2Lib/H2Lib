/*----------------------------------------------------------------------------*/
/* Working with linear p1 elements and hierarchical matrices for laplace                                            */
/*----------------------------------------------------------------------------*/
#include <stdio.h>

#include "parameters.h"		/* Read parameters interactively */

#include "tri2d.h"		/* 2-dimensional mesh */
#include "tri2dp1.h"		/* discretisation with linear p1 elements based on tri2d.h */

#include "tet3d.h"		/* 3-dimensional mesh */
#include "tet3dp1.h"		/* discretisation with linear p1 elements based on tet3d.h */

#include "ddcluster.h"		/* Domain decomposition clustering */
#include "hmatrix.h"		/* Hierarchical matrices */
#include "harith.h"		/* Arithmetic for hierarchical matrices */

#include "truncation.h"		/* Auxiliary function for truncation */
#include "hcoarsen.h"		/* Coarsening of hierarchical matrices */
#include "matrixnorms.h"	/* Computing difference norms for two matrices */

real
norm2lu_sparsematrix(pchmatrix LU, pcsparsematrix sp)
{
  avector   tmp1, tmp2;
  uint      rows = LU->rc->size;
  uint      cols = LU->cc->size;
  pavector  x, y;
  real      norm;
  uint      j;

  assert(sp->rows == rows);
  assert(sp->cols == cols);

  x = init_avector(&tmp1, cols);
  y = init_avector(&tmp2, rows);

  random_avector(x);
  norm = norm2_avector(x);
  scale_avector(1.0 / norm, x);

  for (j = 0; j < NORM_STEPS; j++) {
    // printf("norm = %g \n", sqrt( norm));
    clear_avector(y);
    mvm_sparsematrix_avector(1.0, false, sp, x, y);
    triangularsolve_hmatrix_avector(true, true, false, LU, y);
    triangularsolve_hmatrix_avector(false, false, false, LU, y);
    add_avector(-1.0, y, x);
    copy_avector(x, y);
    triangularsolve_hmatrix_avector(false, false, true, LU, y);
    triangularsolve_hmatrix_avector(true, true, true, LU, y);
    mvm_sparsematrix_avector(-1.0, true, sp, y, x);
    norm = norm2_avector(x);
    scale_avector(1.0 / norm, x);
  }

  uninit_avector(y);
  uninit_avector(x);

  return REAL_SQRT(norm);
}

real
norm2chol_sparsematrix(pchmatrix LU, pcsparsematrix sp)
{
  avector   tmp1, tmp2;
  uint      rows = LU->rc->size;
  uint      cols = LU->cc->size;
  pavector  x, y;
  real      norm;
  uint      j;

  assert(sp->rows == rows);
  assert(sp->cols == cols);

  x = init_avector(&tmp1, cols);
  y = init_avector(&tmp2, rows);

  random_avector(x);
  norm = norm2_avector(x);
  scale_avector(1.0 / norm, x);

  for (j = 0; j < NORM_STEPS; j++) {
    // printf("norm = %g \n", sqrt( norm));
    clear_avector(y);
    mvm_sparsematrix_avector(1.0, false, sp, x, y);
    triangularsolve_hmatrix_avector(true, false, false, LU, y);
    triangularsolve_hmatrix_avector(true, false, true, LU, y);
    add_avector(-1.0, y, x);
    copy_avector(x, y);
    triangularsolve_hmatrix_avector(true, false, false, LU, y);
    triangularsolve_hmatrix_avector(true, false, true, LU, y);
    mvm_sparsematrix_avector(-1.0, true, sp, y, x);
    norm = norm2_avector(x);
    scale_avector(1.0 / norm, x);
  }

  uninit_avector(y);
  uninit_avector(x);

  return REAL_SQRT(norm);
}




int
main(int argc, char **argv)
{

  uint      prob;		/* Used for selecting the problem: 2d or 3d? */
  uint      L;			/* Number of grid refinements */
  uint      clf;		/* Leafsize */
  uint      clustermode;	/* Used for selecting the cluster strategy */
  uint      decomp;		/* Used for selecting the type of the decomposition */
  uint      i, j;		/* Auxiliary variable for loops */
  ptri2d   *gr_2d;		/* 2d mesh hierarchy */
  ptri2dp1  p1_2d;		/* Linear p1 basis functions in 2d */
  ptet3d   *gr_3d;		/* 3d mesh hierarchy */
  ptet3dp1  p1_3d;		/* Linear p1 basis functions in 3d */
  psparsematrix sp;		/* Sparsematrix object */
  uint      ndof;		/* Number of degrees of freedom */
  uint     *idx;		/* Array of indices for the clustergeometry */
  pclustergeometry cg;		/* Clustergeometry object */
  uint      dim;		/* Dimension for domain decomposition clustering */
  uint     *flag;		/* Auxiliary array for domain decompostion clustering */
  pcluster  root;		/* Cluster tree object */
  pblock    broot;		/* Block cluster tree object */
  real      eta;		/* Admissibility parameter */
  phmatrix  hm;			/* Hierarchical matrices */
  real      tol_decomp;		/* Tolerance for decompositions */
  real      tol_coarsen;	/* Tolerance for coarsening of the hierarchical matrix */
  phmatrix  l;			/* Hierarchical matrices for storing the Cholesky factor */
  ptruncmode tm;		/* truncmode for truncations within the arithmetic */
  real      error;		/* Auxiliary variable for testing */
  pstopwatch sw;		/* Stopwatch for time measuring */
  real      time;		/* Variable for time measuring */
  size_t    sz, sz_rk, sz_full;	/* Variables for memory footprint */

  /* First initialise the library */
  init_h2lib(&argc, &argv);

  /*Initialise variables */
  eta = 2.0;
  tol_coarsen = 1e-16;
  tol_decomp = 1e-12;
  tm = new_releucl_truncmode();
  sw = new_stopwatch();

  /* Set problem parameters */
  printf("Select problem\n");
  prob = askforint("1 = fem2d,\t 2 = fem3d\t \n", "h2lib_fem", 1);
  L = askforint("Refinement?\n", "h2lib_L", 4);
  clf = askforint("Leafsize?\n", "h2lib_clf", 32);
  clustermode =
    askforint("1 = adaptive, \t 2 = ddcluster", "h2lib_clusterstrategy", 2);
  decomp =
    askforint("1 = LU-decomposition,\t, 2 = cholesky-decomposition\n",
	      "h2lib_decomp", 1);
  tol_decomp = askforreal("tol_decomp=?\n", "h2lib_tol_decomp", 1e-13);

  /* Build geometry and discretisation for laplace equation */
  if (prob == 1) {
    printf("========================================\n"
	   "  Create and fill fem2d sparsematrix\n");
    /* Mesh hierarchy */
    gr_2d = (ptri2d *) allocmem((size_t) sizeof(ptri2d) * (L + 1));
    gr_2d[0] = new_unitsquare_tri2d();	/* Set domain */
    //gr_2d[0]=new_unitcircle_tri2d();
    //gr_2d[0] = new_lshape_tri2d();
    for (i = 0; i < L; i++) {	/* Mesh refinements */
      gr_2d[i + 1] = refine_tri2d(gr_2d[i], NULL);
    }
    check_tri2d(gr_2d[L]);	/* Check mesh for inconsistencies */
    p1_2d = new_tri2dp1(gr_2d[L]);	/* Build discretisation */
    sp = build_tri2dp1_sparsematrix(p1_2d);	/* Build corresponding sparsematrix */
    assemble_tri2dp1_laplace_sparsematrix(p1_2d, sp, 0);	/* Fill the sparsematrix */
    ndof = p1_2d->ndof;
    printf("ndof = %u\n", ndof);
    /* Initialise index array for the cluster geometry */
    idx = allocuint(ndof);
    for (i = 0; i < ndof; i++)
      idx[i] = i;
    /* Build clustergeometry for the problem */
    cg = build_tri2dp1_clustergeometry(p1_2d, idx);
    del_tri2dp1(p1_2d);
    for (i = 0; i <= L; i++) {
      j = L - i;
      del_tri2d(gr_2d[j]);
    }
    freemem(gr_2d);
    dim = 2;			/* Set dimenison for domain decomposition clustering */
  }
  else {
    assert(prob == 2);
    printf("========================================\n"
	   "  Create and fill fem3d sparsematrix\n");
    /* Mesh hierarchy */
    gr_3d = (ptet3d *) allocmem((size_t) sizeof(ptri2d) * (L + 1));
    gr_3d[0] = new_unitcube_tet3d();	/* Set domain */
    for (i = 0; i < L; i++) {
      gr_3d[i + 1] = refine_tet3d(gr_3d[i], NULL);
    }
    check_tet3d(gr_3d[L]);	/* Check mesh for inconsistencies */
    p1_3d = new_tet3dp1(gr_3d[L]);	/* build discretisation */
    sp = build_tet3dp1_sparsematrix(p1_3d);	/* Build corresponding sparsematrix */
    assemble_tet3dp1_laplace_sparsematrix(p1_3d, sp, 0);	/* Fill the sparsematrix */
    ndof = p1_3d->ndof;
    printf("ndof = %u\n", ndof);
    /* initialise index array for the clustergeometry */
    idx = allocuint(ndof);
    for (i = 0; i < ndof; i++)
      idx[i] = i;
    /* Build clustergeometry for the problem */
    cg = build_tet3dp1_clustergeometry(p1_3d, idx);
    del_tet3dp1(p1_3d);
    for (i = 0; i <= L; i++) {
      j = L - i;
      del_tet3d(gr_3d[j]);
    }
    freemem(gr_3d);
    dim = 3;
  }

  /* Build and fill the corresponding hierarchical matrix */
  printf("========================================\n"
	 "  Create and fill hierarchical matrix\n");
  if (clustermode == 1) {
    /* Build cluster tree based on adaptive clustering */
    root = build_cluster(cg, ndof, idx, clf, H2_ADAPTIVE);
    /* Build block cluster tree using the euclidian (maximum) admissibilty condition */
    broot = build_nonstrict_block(root, root, &eta, admissible_2_cluster);
    /*Build hierarchical matrix from block */
    hm = build_from_block_hmatrix(broot, 0);
    /* Fill hierarchical matrix with sparsematrix entries */
    copy_sparsematrix_hmatrix(sp, hm);
    /* Compute error */
    error = norm2diff_sparsematrix_hmatrix(hm, sp);
    printf("error = %g\n", error);
    del_block(broot);
  }
  else {			/* clustermode == 2 */
    /* Initialise auxiliary array for domain decomposition clustering */
    flag = allocuint(ndof);
    for (i = 0; i < ndof; i++)
      flag[i] = 0;
    /* build cluster tree based adaptive doamin decomposition clustering */
    root = build_adaptive_dd_cluster(cg, ndof, idx, clf, sp, dim, flag);
    freemem(flag);
    /* Build block cluster tree using the domain decompostion admissibilty condition */
    broot = build_nonstrict_block(root, root, &eta, admissible_dd_cluster);
    /* Build hierarchical matrix from block cluster tree */
    hm = build_from_block_hmatrix(broot, 0);
    /* Fill hierarchical matrix with sparsematrix entries */
    copy_sparsematrix_hmatrix(sp, hm);
    /* Compute error */
    error = norm2diff_sparsematrix_hmatrix(hm, sp);
    printf("error = %g\n", error);
    del_block(broot);
  }

  /* Draw hierarchical matrix with cairo */
#if 0
#ifdef USE_CAIRO
  cairo_t  *cr;
  printf("Draw hmatrix to \"hm_p1.pdf\"\n");
  cr = new_cairopdf("../hm_p1.pdf", 1024.0, 1024.0);
  draw_cairo_hmatrix(cr, hm, false, 0);
  cairo_destroy(cr);
#endif
#endif

  /* Compute size of the hierarchical matrix */
  sz = getsize_hmatrix(hm);
  sz_rk = getfarsize_hmatrix(hm);
  sz_full = getnearsize_hmatrix(hm);
  printf("size original matrix:\t \t %.2f MB \t %.2f KB \t %.2f KB/DoF\n",
	 sz / 1024.0 / 1024.0, sz / 1024.0, sz / 1024.0 / ndof);
  printf("size rk: \t %.2f MB\t size full: \t %.2f MB \n",
	 sz_rk / 1024.0 / 1024.0, sz_full / 1024.0 / 1024.0);
  printf("Coarsen hmatrix\n");

  /* Coarsening of the hierarchical matrix to safe storage */
  coarsen_hmatrix(hm, tm, tol_coarsen, true);

  /* Compute size of the hierarchical matrix after coarsening */
  sz = getsize_hmatrix(hm);
  sz_rk = getfarsize_hmatrix(hm);
  sz_full = getnearsize_hmatrix(hm);
  printf("size original matrix:\t \t %.2f MB \t %.2f KB \t %.2f KB/DoF\n",
	 sz / 1024.0 / 1024.0, sz / 1024.0, sz / 1024.0 / ndof);
  printf("size rk: \t %.2f MB\t size full: \t %.2f MB \n",
	 sz_rk / 1024.0 / 1024.0, sz_full / 1024.0 / 1024.0);
  /* Compute error after coarseing */
  error = norm2diff_sparsematrix_hmatrix(hm, sp);
  printf("error = %g\n", error);

  /* Draw the hierarchical matrix after coarsening */
#if 0
#ifdef USE_CAIRO
  cairo_t  *cr2;
  printf("Draw hmatrix to \"hm_p1_coarsend.pdf\"\n");
  cr2 = new_cairopdf("../hm_p1_coarsend.pdf", 1024.0, 1024.0);
  draw_cairo_hmatrix(cr2, hm, false, 0);
  cairo_destroy(cr2);
#endif
#endif

  /* Compute decomposition of the hierarchical matrix */
  if (decomp == 1) {
    printf("========================================\n"
	   "  Test lu-decomposition\n");
    /* Compute size of the hierarchical matrix */
    sz = getsize_hmatrix(hm);
    printf("size original matrix:\t \t %.2f MB \t %.2f KB \t %.2f KB/DoF\n",
	   sz / 1024.0 / 1024.0, sz / 1024.0, sz / 1024.0 / ndof);
    /* Compute lu-decomposition of the hierarchical matrix */
    start_stopwatch(sw);
    lrdecomp_hmatrix(hm, 0, tol_decomp);
    time = stop_stopwatch(sw);
    /* Compute size of the lu-decomposition */
    sz = getsize_hmatrix(hm);
    sz_rk = getfarsize_hmatrix(hm);
    sz_full = getnearsize_hmatrix(hm);
    printf
      ("size lu-decomposition (lu):\t %.2f MB \t %.2f KB \t %.2f KB/DoF\n",
       sz / 1024.0 / 1024.0, sz / 1024.0, sz / 1024.0 / ndof);
    printf("size rk: \t %.2f MB\t size full: \t %.2f MB \n",
	   sz_rk / 1024.0 / 1024.0, sz_full / 1024.0 / 1024.0);
    printf("time = %f s\n", time);
#if 0
#ifdef USE_CAIRO
    cairo_t  *cr;
    printf("Draw lu-decomposition to \"hm_p1_lu.pdf\"\n");
    cr = new_cairopdf("../hm_p1_lu_staging.pdf", 512.0, 512.0);
    draw_cairo_hmatrix(cr, hm, false, 100);
    cairo_destroy(cr);
#endif
#endif
    /* Compute error */
    error = norm2lu_sparsematrix(hm, sp);
    printf("||I-(lu)^-1 sp||_2 ");
    printf("error = %g\n", error);
  }
  else {			/*decomp == 2 */
    printf("========================================\n"
	   "  Test cholesky-decomposition\n");
    /* Compute size of the hierarchical matrix */
    sz = getsize_hmatrix(hm);
    printf("size original matrix:\t %.2f MB \t %.2f KB \t %.2f KB/DoF\n",
	   sz / 1024.0 / 1024.0, sz / 1024.0, sz / 1024.0 / ndof);
    /* Compute cholesky-decomposition */
    start_stopwatch(sw);
    choldecomp_hmatrix(hm, 0, tol_decomp);
    time = stop_stopwatch(sw);
    /* Compute size of the cholesky-factor */
    l = clone_lower_hmatrix(false, hm);
    sz = getsize_hmatrix(l);
    sz_rk = getfarsize_hmatrix(l);
    sz_full = getnearsize_hmatrix(l);
    printf
      ("size chol-decomposition (l):\t %.2f MB \t %.2f KB \t %.2f KB/DoF\n",
       sz / 1024.0 / 1024.0, sz / 1024.0, sz / 1024.0 / ndof);
    printf("size rk: \t %.2f MB\t size full: \t %.2f MB \n",
	   sz_rk / 1024.0 / 1024.0, sz_full / 1024.0 / 1024.0);
    printf("time = %f s\n", time);
#if 0
#ifdef USE_CAIRO
    cairo_t  *cr;
    printf("Draw chol-factor to \"hm_chol_staging.pdf\"\n");
    cr = new_cairopdf("../hm_chol_staging.pdf", 512.0, 512.0);
    draw_cairo_hmatrix(cr, l, false, 1000);
    cairo_destroy(cr);
#endif
#endif
    /* Compute error */
    error = norm2chol_sparsematrix(hm, sp);
    printf("||I-(ll*)^-1 sp||_2 ");
    printf("error = %g\n", error);

    del_hmatrix(l);
  }

  /* Cleaning up */
  del_truncmode(tm);
  del_stopwatch(sw);
  del_hmatrix(hm);

  del_cluster(root);
  del_clustergeometry(cg);
  del_sparsematrix(sp);
  freemem(idx);
  (void) printf("  %u matrices and\n"
		"  %u vectors still active\n",
		getactives_amatrix(), getactives_avector());
  uninit_h2lib();
  printf("The end\n");
  return 0;
}
