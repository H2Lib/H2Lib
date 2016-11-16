
#include <math.h>
#include <stdio.h>

#include "basic.h"
#include "dblock.h"
#include "block.h"
#include "dclusterbasis.h"
#include "surface3d.h"
#include "macrosurface3d.h"
#include "clustergeometry.h"
#include "parameters.h"
#include "dh2compression.h"
#include "helmholtzbem3d.h"
#include "matrixnorms.h"

#ifndef USE_COMPLEX
int
main()
{
  return 0;
}
#else
int
main(int argc, char **argv)
{
  psurface3d gr;		/* BEM grid */
  pclustergeometry cgeo;	/* Description of cluster geometry */
  uint     *idx;		/* Index permutation */
  pcluster  root;		/* Root of cluster tree */
  pdcluster droot;		/* Root of directional cluster tree */
  pleveldir ld;			/* Levelwise directions */
  uint      res;		/* Cluster resolution */
  real      wave_k;		/* Wave number */
  real      wave_k_3;		/* Wave number / sqrt(3) */
  pfield    kvec;
  real      eta1;		/* Directional admissibility parameter */
  real      eta2;		/* Parabolic admissibility parameter */
  real      eta_eff;		/* Effective admissibility parameter */
  uint      nq;			/* Quadrature order */
  char      operator;		/* Operator */
  diradmdata dad;		/* Parameters for admissibility condition */
  pdblock   broot;		/* Root of block tree */
  pbem3d    bm;
  pamatrix  G;			/* Stiffness matrix */
  pdclusterbasis rb, cb;	/* Directional cluster bases */
  pdh2matrix Gh;		/* Directional H^2-matrix */
  real      norm, error;
  pstopwatch sw;
  ptruncmode tm;
#ifdef USE_CAIRO
  cairo_t  *cr;
#endif
  real      eps;
  real      hmin, hmax, anglemin, angleedge;
  real      t_setup, t_test;
  bool      gr_closed, gr_oriented;
  uint      n, dirs;
  size_t    sz, tsz;
  uint      i;

  init_h2lib(&argc, &argv);	/* Initialize the library */

  sw = new_stopwatch();

  tm = new_blockreleucl_truncmode();

  gr = build_interactive_surface3d();

  prepare_surface3d(gr);
  n = gr->triangles;

  getproperties_surface3d(gr, &hmin, &hmax, &anglemin, &angleedge);

  (void) printf("  %u vertices, %u edges, %u triangles\n"
		"  hmin %.2e, hmax %.2e, anglemin %g, angleedge %g\n",
		gr->vertices, gr->edges, gr->triangles,
		hmin, hmax, anglemin, angleedge);

  gr_closed = isclosed_surface3d(gr);
  gr_oriented = isoriented_surface3d(gr);

  (void) printf("  Surface is%s closed\n"
		"  Surface is%s oriented\n"
		"----------------------------------------\n",
		(gr_closed ? "" : " NOT\a"), (gr_oriented ? "" : " NOT\a"));

  res = askforint("Cluster resolution?", "h2lib_leafsize", 16);

  operator = askforchar("S)ingle or d)ouble layer operator?",
			"h2lib_operator", "sd", 's');

  wave_k = askforreal("Wave number?", "h2lib_wavek", 1.0);

  eta1 = askforreal("Directional admissibility parameter?", "h2lib_eta1",
		    10.0);

  eta2 = askforreal("Parabolic admissibility parameter?", "h2lib_eta2", 2.0);

  nq = askforint("Quadrature order?", "h2lib_quadorder", 3);

  eps = askforreal("Compression tolerance?", "h2lib_compeps", 1e-4);

  (void) printf("========================================\n"
		"Creating BEM object\n");
  kvec = allocfield(3);
  wave_k_3 = wave_k / REAL_SQRT(3.0);
  kvec[0] = wave_k_3;
  kvec[1] = wave_k_3;
  kvec[2] = wave_k_3;
  bm = (operator == 's' ?
	new_slp_helmholtz_bem3d(kvec, gr, nq, nq,
				BASIS_CONSTANT_BEM3D) :
	new_dlp_helmholtz_bem3d(kvec, gr, nq, nq,
				BASIS_CONSTANT_BEM3D, BASIS_CONSTANT_BEM3D,
				0.0));

  (void) printf("Creating cluster tree\n");
  cgeo = build_bem3d_const_clustergeometry(bm, &idx);
  root = build_adaptive_cluster(cgeo, n, idx, res);
  (void) printf("  %u clusters\n", root->desc);
  del_clustergeometry(cgeo);

  (void) printf("Creating directional cluster tree\n");
  droot = buildfromcluster_dcluster(root);
  ld = builddirections_box_dcluster(droot, eta1 / wave_k);
  dirs = getalldirections_dcluster(droot);
  for (i = 0; i <= ld->depth; i++)
    (void) printf("  Level %2u: maxdiam %g, %u directions\n",
		  i, ld->maxdiam[i], ld->directions[i]);
  (void) printf("  %u directions (%.1f per cluster)\n",
		dirs, (real) dirs / root->desc);

  (void) printf("Creating directional block tree\n");
  dad.eta1 = eta1;
  dad.eta2 = eta2;
  dad.wave_k = wave_k;
  dad.xy = allocreal(3);
  dad.ld = ld;
  broot = build_dblock(droot, droot, 0, parabolic_admissibility, &dad);
  eta_eff = getmaxeta_dblock(broot);
  (void) printf("  %u blocks (%.1f per cluster)\n"
		"  %.1f effective eta\n",
		broot->desc, (real) broot->desc / root->desc, eta_eff);

#ifdef USE_CAIRO
  (void) printf("Drawing block tree to \"dblock.png\"\n");
  cr = new_cairopng(1024, 1024);
  cairodraw_dblock(cr, broot, 0);
  write_cairopng(cr, "dblock.png");
  cairo_destroy(cr);
#endif

#ifdef USE_OPENMP
  (void) printf("Maximal parallelization depth %u\n", max_pardepth);
#endif

  (void) printf("----------------------------------------\n"
		"Creating dense matrix\n");
  G = new_amatrix(n, n);
  sz = getsize_amatrix(G);
  (void) printf("  %.1f MB (%.1f KB/DoF)\n", sz / 1048576.0, sz / 1024.0 / n);

  (void) printf("Filling dense matrix\n");
  start_stopwatch(sw);

  bm->nearfield(0, 0, bm, false, G);

  t_setup = stop_stopwatch(sw);
  sz = getsize_amatrix(G);
  (void) printf("  %.1f seconds\n", t_setup);
  (void) fflush(stdout);

  (void) printf("Computing spectral norm\n");
  start_stopwatch(sw);
  norm = norm2_amatrix(G);
  t_test = stop_stopwatch(sw);
  (void) printf("  Spectral norm of G: %.4e\n"
		"  %.2f seconds\n", norm, t_test);

  (void) printf("----------------------------------------\n"
		"Computing adaptive directional row basis, eps=%g\n", eps);
  start_stopwatch(sw);
  rb = buildrow_amatrix_dclusterbasis(G, broot, tm, eps);
  t_setup = stop_stopwatch(sw);
  sz = getsize_dclusterbasis(rb);
  (void) printf("  %.2f MB (%.2f KB/DoF)\n"
		"  %.2f seconds\n"
		"  Maximal rank %u\n"
		"  Rank sum %u\n"
		"  Rank branch sum %u\n",
		sz / 1048576.0, sz / 1024.0 / n,
		t_setup,
		getmaxrank_dclusterbasis(rb), rb->ktree, rb->kbranch);

  (void) printf("Computing adaptive directional column basis, eps=%g\n", eps);
  start_stopwatch(sw);
  cb = buildcol_amatrix_dclusterbasis(G, broot, tm, eps);
  t_setup = stop_stopwatch(sw);
  sz = getsize_dclusterbasis(cb);
  (void) printf("  %.2f MB (%.2f KB/DoF)\n"
		"  %.2f seconds\n"
		"  Maximal rank %u\n"
		"  Rank sum %u\n"
		"  Rank branch sum %u\n",
		sz / 1048576.0, sz / 1024.0 / n,
		t_setup,
		getmaxrank_dclusterbasis(cb), cb->ktree, cb->kbranch);

  (void) printf("Creating directional H^2-matrix\n");
  start_stopwatch(sw);
  Gh = buildfromblock_dh2matrix(broot, rb, cb);
  t_setup = stop_stopwatch(sw);
  sz = getsize_dh2matrix(Gh);
  tsz = gettotalsize_dh2matrix(Gh);
  (void) printf("  %.2f MB (%.2f KB/DoF)\n"
		"  Total %.2f MB (%.2f KB/DoF)\n"
		"  %.2f seconds\n",
		sz / 1048576.0, sz / 1024.0 / n,
		tsz / 1048576.0, tsz / 1024.0 / n, t_setup);

  (void) printf("Filling directional H^2-matrix\n");
  start_stopwatch(sw);
  projectdense_dh2matrix(G, Gh);
  t_setup = stop_stopwatch(sw);
  (void) printf("  %.2f seconds\n", t_setup);
  (void) fflush(stdout);

  (void) printf("Computing error norm\n");
  start_stopwatch(sw);
  error = norm2diff_amatrix_dh2matrix(Gh, G);
  t_test = stop_stopwatch(sw);
  (void) printf("  Spectral error: %.4e (%.4e)\n"
		"  %.2f seconds\n", error, error / norm, t_test);
  (void) fflush(stdout);

#ifdef USE_CAIRO
  (void) printf("Drawing directional DH2-matrix to \"dh2.png\"\n");
  cr = new_cairopng(1024, 1024);
  draw_cairo_dh2matrix(cr, Gh, true, true, 0);
  write_cairopng(cr, "dh2.png");
  cairo_destroy(cr);
#endif

  del_dh2matrix(Gh);
  del_dclusterbasis(cb);
  del_dclusterbasis(rb);

  (void) printf("========================================\n" "Cleaning up\n");

  del_amatrix(G);
  del_dblock(broot);
  del_dcluster(droot);
  del_leveldir(ld);
  del_cluster(root);
  freemem(idx);
  del_surface3d(gr);
  del_truncmode(tm);
  del_stopwatch(sw);

  (void) printf("  %u directional cluster bases,\n"
		"  %u matrices and\n"
		"  %u vectors still active\n",
		getactives_dclusterbasis(),
		getactives_amatrix(), getactives_avector());

  uninit_h2lib();

  return 0;
}
#endif
