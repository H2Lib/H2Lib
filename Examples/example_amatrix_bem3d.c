#include <stdio.h>

#include "basic.h"
#include "krylovsolvers.h"
#include "laplacebem3d.h"

/****************************************************
 * This examples sets up single layer potential operator(SLP),
 * double layer potential operator(DLP) as well as the mass matrix M
 * as dense matrices in order to solve the interior Dirichlet
 * problem for the Laplace equation.
 ****************************************************/

int
main(int argc, char **argv)
{
  pstopwatch sw;
  pmacrosurface3d mg;
  psurface3d gr;
  pbem3d    bem_slp, bem_dlp;
  uint      q_reg, q_sing;
  basisfunctionbem3d basis;
  pamatrix  V, KM;
  pavector  gd, b, x;
  real      eps_solve;
  uint      maxiter;
  real      t, size, norm;

  /* Init the H2Lib, should be called before any other function. */
  init_h2lib(&argc, &argv);

  /****************************************************
   * Set up basic parameters
   ****************************************************/

  /* Number of quadrature points for regular integrals. */
  q_reg = 2;
  /* Number of quadrature points for singular integrals. */
  q_sing = q_reg + 2;
  /* Basis functions that should be used. */
  basis = BASIS_CONSTANT_BEM3D;

  /* absolute norm of the residuum for CG-method */
  eps_solve = 1.0e-10;

  /* maximum number of CG-steps that should be performed. */
  maxiter = 500;

  /* Stopwatch for measuring the time. */
  sw = new_stopwatch();

  /****************************************************
   * Create geometry
   ****************************************************/

  /* Create abstract geometry of a sphere. */
  mg = new_sphere_macrosurface3d();
  /* Mesh the abstract geometry with 32 levels of refinement. */
  gr = build_from_macrosurface3d_surface3d(mg, 32);
  printf("Created geometry with %d vertices, %d edges and %d triangles\n",
	 gr->vertices, gr->edges, gr->triangles);

  /****************************************************
   * Set up bem objects
   ****************************************************/

  /* Create a new BEM-object, that can compute entries of SLP operator. */
  bem_slp = new_slp_laplace_bem3d(gr, q_reg, q_sing, basis, basis);
  /* Create a new BEM-object, that can compute entries of DLP operator
   * and 0.5*I. */
  bem_dlp = new_dlp_laplace_bem3d(gr, q_reg, q_sing, basis, basis, 0.5);

  /****************************************************
   * Assemble H-matrix SLP
   ****************************************************/

  printf("Assemble dense matrix V:\n");

  /* Create amatrix structure. */
  V = new_amatrix(gr->triangles, gr->triangles);

  start_stopwatch(sw);
  /* Assemble entries of V. */
  assemble_bem3d_amatrix(bem_slp, V);
  t = stop_stopwatch(sw);
  /* Get the total memory footprint for V. */
  size = getsize_amatrix(V) / 1024.0 / 1024.0;

  printf("  %.2f s\n", t);
  printf("  %.3f MB\n", size);

  /****************************************************
   * Assemble H-matrix 0.5I + K
   ****************************************************/

  printf("Assemble dense matrix 0.5M + K:\n");

  /* Create amatrix structure. */
  KM = new_amatrix(gr->triangles, gr->triangles);

  start_stopwatch(sw);
  /* Assemble entries of KM. */
  assemble_bem3d_amatrix(bem_dlp, KM);
  t = stop_stopwatch(sw);
  /* Get the total memory footprint for KM. */
  size = getsize_amatrix(KM) / 1024.0 / 1024.0;

  printf("  %.2f s\n", t);
  printf("  %.3f MB\n", size);

  /****************************************************
   * Create Dirichlet data
   ****************************************************/

  /* Create new vector to store L2-projection of the Dirichlet data. */
  gd = new_avector(gr->triangles);
  printf("Compute L2-projection of Dirichlet data:\n");
  start_stopwatch(sw);
  /* L2-projection onto the space of piecewise constant function
   * on the boundary. */
  projectL2_bem3d_c_avector(bem_dlp, eval_dirichlet_fundamental_laplacebem3d,
			    gd, (void *) bem_dlp);
  t = stop_stopwatch(sw);
  size = getsize_avector(gd) / 1024.0 / 1024.0;
  printf("  %.2f s\n", t);
  printf("  %.3f MB\n", size);

  /****************************************************
   * Compute right-hand-side b = (0.5M + K)*gd
   ****************************************************/
  /* Create new vector to store right-hand-side. */
  b = new_avector(gr->triangles);
  printf("Compute right-hand-side:\n");
  start_stopwatch(sw);
  clear_avector(b);
  /* H-matrix vector product. */
  addeval_amatrix_avector(1.0, KM, gd, b);
  t = stop_stopwatch(sw);
  size = getsize_avector(b) / 1024.0 / 1024.0;
  printf("  %.2f s\n", t);
  printf("  %.3f MB\n", size);

  /****************************************************
   * Solve linear system V x = b using CG-method.
   ****************************************************/

  /* Create new vector to store the solution coefficients. */
  x = new_avector(gr->triangles);
  printf("Solve linear system:\n");
  start_stopwatch(sw);
  /* Call the CG-solver for H-matrices. */
  solve_cg_amatrix_avector(V, b, x, eps_solve, maxiter);
  t = stop_stopwatch(sw);
  size = getsize_avector(x) / 1024.0 / 1024.0;
  printf("  %.2f s\n", t);
  printf("  %.3f MB\n", size);

  /****************************************************
   * Compute L2-error compared to analytical solution
   * of the Neumann data.
   ****************************************************/

  printf
    ("Compute L2-error against analytical solution of the Neumann data:\n");
  start_stopwatch(sw);
  /* L2-norm of gn(x) - \sum_i x_i \varphi_(x)  */
  norm = normL2diff_c_bem3d(bem_slp, x, eval_neumann_fundamental_laplacebem3d,
			    NULL);
  t = stop_stopwatch(sw);
  printf("Abs. L2-error:\n");
  printf("  %.2f s\n", t);
  printf("  %.5e\n", norm);

  printf("Rel. L2-error:\n");
  start_stopwatch(sw);
  /* Compute L2-norm of gn and update norm to relative L2-norm.  */
  norm /= normL2_bem3d(bem_slp, eval_neumann_fundamental_laplacebem3d, NULL);
  t = stop_stopwatch(sw);
  printf("  %.2f s\n", t);
  printf("  %.5e\n", norm);

  /****************************************************
   * cleanup
   ****************************************************/

  del_avector(x);
  del_avector(b);
  del_avector(gd);
  del_amatrix(V);
  del_amatrix(KM);
  del_bem3d(bem_slp);
  del_bem3d(bem_dlp);
  del_macrosurface3d(mg);
  del_surface3d(gr);
  del_stopwatch(sw);

  /* Uninit the H2Lib. */
  uninit_h2lib();

  return 0;
}
