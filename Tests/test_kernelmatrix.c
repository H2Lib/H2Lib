
#include "kernelmatrix.h"

#include <stdio.h>

#include "basic.h"
#include "h2compression.h"
#include "h2matrix.h"
#include "parameters.h"
#include "matrixnorms.h"

static field
kernel_newton(const real *xx, const real *yy, void *data)
{
  real norm2;

  (void) data;

  norm2 = REAL_SQR(xx[0] - yy[0]) + REAL_SQR(xx[1] - yy[1]);

  return (norm2 == 0.0 ? 0.0 : 1.0 / REAL_SQRT(norm2));
}

static field
kernel_exp(const real *xx, const real *yy, void *data)
{
  real norm2;

  (void) data;

  norm2 = REAL_SQR(xx[0] - yy[0]) + REAL_SQR(xx[1] - yy[1]);

  return REAL_EXP(-norm2);
}

static field
kernel_log(const real *xx, const real *yy, void *data)
{
  real norm2;

  (void) data;

  norm2 = REAL_SQR(xx[0] - yy[0]) + REAL_SQR(xx[1] - yy[1]);

  return (norm2 == 0.0 ? 0.0 : -0.5*REAL_LOG(norm2));
}

int
main(int argc, char **argv)
{
  pkernelmatrix km;
  pclustergeometry cg;
  pcluster root;
  pblock broot;
  pclusterbasis cb;
  ph2matrix Gh1, Gh2;
  pamatrix G;
  pstopwatch sw;
  char kernel;
  uint points;
  uint m, leafsize;
  uint *idx;
  size_t sz;
  real eps, eta;
  real t_setup, norm, error;
  uint i;

  init_h2lib(&argc, &argv);

  sw = new_stopwatch();

  kernel = askforchar("Kernel function? N)ewton, L)ogarithmic, or E)xponential?", "h2lib_kernelfunc", "nle", 'l');

  points = askforint("Number of points?", "h2lib_kernelpoints", 8192);

  m = askforint("Interpolation order?", "h2lib_interorder", 5);

  leafsize = askforint("Cluster resolution?", "h2lib_leafsize", 2*m*m);

  eps = askforreal("Recompression tolerance?", "h2lib_comptol", 1e-4);

  eta = 2.0;
  
  (void) printf("Creating kernelmatrix object for %u points, order %u\n",
		points, m);
  km = new_kernelmatrix(2, points, m);
  switch(kernel) {
  case 'e':
    (void) printf("  Exponential kernel function\n");
    km->kernel = kernel_exp;
    break;
  case 'n':
    (void) printf("  Newton kernel function\n");
    km->kernel = kernel_newton;
    break;
  default:
    (void) printf("  Logarithmic kernel function\n");
    km->kernel = kernel_log;
  }

  /* Random points in [-1,1]^2 */
  for(i=0; i<points; i++) {
    km->x[i][0] = 2.0 * rand() / RAND_MAX;
    km->x[i][1] = 2.0 * rand() / RAND_MAX;
  }

  (void) printf("Creating clustergeometry object\n");
  cg = creategeometry_kernelmatrix(km);

  (void) printf("Creating cluster tree\n");
  idx = (uint *) allocmem(sizeof(uint) * points);
  for(i=0; i<points; i++)
    idx[i] = i;
  root = build_adaptive_cluster(cg, points, idx, leafsize);
  (void) printf("  %u clusters, depth %u\n",
		root->desc, getdepth_cluster(root));

  (void) printf("Creating block tree\n");
  broot = build_strict_block(root, root, &eta, admissible_2_cluster);
  (void) printf("  %u blocks, depth %u\n",
		broot->desc, getdepth_block(broot));

  (void) printf("Creating cluster basis\n");
  cb = build_from_cluster_clusterbasis(root);

  (void) printf("Filling cluster basis\n");
  start_stopwatch(sw);
  fill_clusterbasis_kernelmatrix(km, cb);
  t_setup = stop_stopwatch(sw);
  sz = getsize_clusterbasis(cb);
  (void) printf("  %.2f seconds\n"
		"  %.1f MB\n"
		"  %.1f KB/DoF\n",
		t_setup, sz / 1048576.0, sz / 1024.0 / points);

  (void) printf("Creating H^2-matrix\n");
  Gh1 = build_from_block_h2matrix(broot, cb, cb);
  sz = getsize_h2matrix(Gh1);
  (void) printf("  %.1f MB\n"
		"  %.1f KB/DoF\n",
		sz / 1048576.0, sz / 1024.0 / points);

  (void) printf("Filling H^2-matrix\n");
  start_stopwatch(sw);
  fill_h2matrix_kernelmatrix(km, Gh1);
  t_setup = stop_stopwatch(sw);
  (void) printf("  %.2f seconds\n",
		t_setup);

  (void) printf("Filling reference matrix\n");
  G = new_amatrix(points, points);
  start_stopwatch(sw);
  fillN_kernelmatrix(0, 0, km, G);
  t_setup = stop_stopwatch(sw);
  sz = getsize_amatrix(G);
  (void) printf("  %.2f seconds\n"
		"  %.1f MB\n"
		"  %.1f KB/DoF\n",
		t_setup, sz / 1048576.0, sz / 1024.0 / points);

  (void) printf("Computing norm\n");
  norm = norm2_amatrix(G);
  (void) printf("  Spectral norm %.3e\n",
		norm);

  (void) printf("Computing approximation error\n");
  error = norm2diff_amatrix_h2matrix(Gh1, G);
  (void) printf("  Spectral error %.3e (%.3e)\n",
		error, error/norm);
  
  (void) printf("Recompressing H^2-matrix, eps=%g\n",
		eps);
  start_stopwatch(sw);
  Gh2 = compress_h2matrix_h2matrix(Gh1, false, false, 0, eps);
  t_setup = stop_stopwatch(sw);
  sz = getsize_h2matrix(Gh2);
  (void) printf("  %.2f seconds\n"
		"  %.1f MB\n"
		"  %.1f KB/DoF\n",
		t_setup, sz / 1048576.0, sz / 1024.0 / points);

  (void) printf("Computing approximation error\n");
  error = norm2diff_amatrix_h2matrix(Gh2, G);
  (void) printf("  Spectral error %.3e (%.3e)\n",
		error, error/norm);

  uninit_h2lib();

  return 0;
}
