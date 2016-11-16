

#include <stdio.h>

#ifdef USE_CAIRO
#include <cairo/cairo.h>
#endif

#include "basic.h"
#include "cluster.h"
#include "clustergeometry.h"
#include "block.h"
#include "tri2d.h"
#include "tri2dp1.h"
#include "ddcluster.h"
#include "hmatrix.h"
#include "matrixnorms.h"

static uint problems = 0;
#define IS_IN_RANGE(a, b, c) (((a) <= (b)) && ((b) <= (c)))

int
main(int argc, char **argv)
{
  uint      i, j;
  uint      L;			/* Number of grid refinements */
  uint      clf;		/* Maximum leaf size */
  real      eta;		/* Admissibility parameter */
  ptri2d   *gr;			/* Grids */
  ptri2dp1  p1;			/* P1 elements on grid gr[L] */
  pclustergeometry cg;		/* Clustergeometry structure for clustering */
  uint     *idx;		/* Index array for cluster */
  pcluster  root;		/* Cluster tree (for rows and columns) */
  pblock    broot;		/* Block tree */
  psparsematrix sp;		/* Sparsematrix */
  uint      dim;		/* Dimension for splitting ddcluster */
  uint     *flag;		/* Auxiliary array for dd-cluster */
  phmatrix  hm;			/* Hierarchical matrix */
  real      error;

  init_h2lib(&argc, &argv);

  L = 7;
  eta = 2.0;
  clf = 16;
  dim = 2;
  printf("========================================\n"
         "Testing ddcluster\n");
  printf("========================================\n"
	 "  Creating grid hierarchy\n");
  gr = (ptri2d *) allocmem((size_t) sizeof(ptri2d) * (L + 1));
  gr[0] = new_unitsquare_tri2d();
  printf("    Level 0: %u vertices, %u edges,"
	 " %u triangles\n", gr[0]->vertices, gr[0]->edges, gr[0]->triangles);
  for (i = 0; i < L; i++) {
    gr[i + 1] = refine_tri2d(gr[i], NULL);
    (void) printf("    Level %u: %u vertices, %u edges,"
		  " %u triangles\n",
		  i + 1,
		  gr[i + 1]->vertices, gr[i + 1]->edges,
		  gr[i + 1]->triangles);
  }
  printf("  Checking grid on level %d\n", L);
  check_tri2d(gr[L]);

  printf("========================================\n"
	 "  Creating tri2dp1 structure\n");
  p1 = new_tri2dp1(gr[L]);

  printf("========================================\n"
	 "  Building clustergeometry object\n");
  idx = allocuint(p1->ndof);
  for (i = 0; i < p1->ndof; i++)
    idx[i] = i;
  cg = build_tri2dp1_clustergeometry(p1, idx);

  printf("========================================\n"
	 "  Building sparsematrix\n");
  sp = build_tri2dp1_sparsematrix(p1);
  assemble_tri2dp1_laplace_sparsematrix(p1, sp, NULL);
  //print_sparsematrix(sp);

  printf("========================================\n"
	 "  Building domain decomposition cluster tree\n");
  //printf("f2->nidx = %u \n", f2->nidx);
  flag = allocuint(p1->ndof);
  for (i = 0; i < p1->ndof; i++)
    flag[i] = 0;

  //printf("0: hmin = %f, hmax = %f", cg->hmin[0], cg->hmax[0]);
  //printf("1: hmin = %f, hmax = %f \n", cg->hmin[1], cg->hmax[1]);
  root = build_adaptive_dd_cluster(cg, p1->ndof, idx, clf, sp, dim, flag);

  printf("========================================\n"
	 "  Building non strict block cluster tree\n");

  broot = build_nonstrict_block(root, root, &eta, admissible_dd_cluster);

#ifdef USE_CAIRO
  cairo_t  *cr;
  printf("  Drawing block tree to \"block.pdf\"\n");
  cr = new_cairopdf("../block.pdf", 512.0, 512.0);
  draw_cairo_block(cr, broot, 0);
  cairo_destroy(cr);
  printf("  Drawing grid to \"grid.pdf\"\n");
  draw_cairo_tri2d(gr[L], "../grid.pdf", false, -1);

#endif

  printf("========================================\n"
	 "  Building hierarchical matrix\n");

  hm = build_from_block_hmatrix(broot, 0);
  copy_sparsematrix_hmatrix(sp, hm);
  error = norm2diff_sparsematrix_hmatrix(hm, sp);
  printf("error = %g\n", error);

  if (!IS_IN_RANGE(0.0, error, 1.0e-16))
    problems++;

  printf("========================================\n" "  Cleaning up\n");
  for (i = 0; i <= L; i++) {
    j = L - i;
    del_tri2d(gr[j]);
  }
  freemem(gr);
  del_tri2dp1(p1);

  del_clustergeometry(cg);
  del_block(broot);
  del_cluster(root);
  del_hmatrix(hm);

  del_sparsematrix(sp);
  freemem(idx);
  freemem(flag);


  uninit_h2lib();

  return problems;
}
