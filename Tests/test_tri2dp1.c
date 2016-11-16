
#include "tri2dp1.h"
#include "parameters.h"
#include "krylov.h"
#include "krylovsolvers.h"
#include "basic.h"

#include <stdio.h>

static uint problems = 0;

#define IS_IN_RANGE(a, b, c) (((a) <= (b)) && ((b) <= (c)))


static    field
sin_solution(const real * x, void *data)
{
  (void) data;

  return sin(M_PI * x[0]) * sin(M_PI * x[1]);
}

static    field
sin_rhs(const real * x, void *data)
{
  (void) data;

  return 2.0 * M_PI * M_PI * sin(M_PI * x[0]) * sin(M_PI * x[1]);
}

int
main(int argc, char **argv)
{
  ptri2d   *gr;
  ptri2dp1 *dc;
  ptri2dref *rf;
  psparsematrix A, Af;
  pavector  xd, b, x;
  uint      L;
  real      error;
  pstopwatch sw;
  real      runtime;
  uint      i;
  real      eps;
  uint      steps;
  uint      iter;

  init_h2lib(&argc, &argv);

  sw = new_stopwatch();

  L = 8;
  eps = 1e-12;
  steps = 2500;

  (void) printf("----------------------------------------\n"
		"Creating mesh hierarchy\n");
  gr = (ptri2d *) allocmem(sizeof(ptri2d) * (L + 1));
  rf = (ptri2dref *) allocmem(sizeof(ptri2dref) * L);

  gr[0] = new_unitsquare_tri2d();
  (void) printf("  Level %2u: %u vertices, %u edges, %u triangles\n",
		0, gr[0]->vertices, gr[0]->edges, gr[0]->triangles);
  for (i = 0; i < L; i++) {
    gr[i + 1] = refine_tri2d(gr[i], rf + i);
    (void) printf("  Level %2u: %u vertices, %u edges, %u triangles\n",
		  i + 1, gr[i + 1]->vertices, gr[i + 1]->edges,
		  gr[i + 1]->triangles);
  }

  printf("Draw grid Level %u\n", 3);
  draw_cairo_tri2d(gr[3], "mesh", 0, 0);

  (void) printf("Creating discretizations\n");
  dc = (ptri2dp1 *) allocmem(sizeof(ptri2dp1) * (L + 1));
  for (i = 0; i <= L; i++) {
    dc[i] = new_tri2dp1(gr[i]);
    (void) printf("  Level %2u: %u degrees of freedom, %u fixed\n",
		  i, dc[i]->ndof, dc[i]->nfix);
  }

  for (i = 6; i <= L; i++) {
    (void) printf("Testing level %u\n" "  Setting up matrix\n", i);
    start_stopwatch(sw);
    A = build_tri2dp1_sparsematrix(dc[i]);
    Af = build_tri2dp1_interaction_sparsematrix(dc[i]);
    assemble_tri2dp1_laplace_sparsematrix(dc[i], A, Af);
    runtime = stop_stopwatch(sw);
    (void) printf("  %.1f seconds\n"
		  "  %.1f MB (interaction %.1f MB)\n"
		  "  %.1f KB/DoF (interaction %.1f KB/DoF)\n"
		  "  %u non-zero entries (interaction %u)\n",
		  runtime,
		  getsize_sparsematrix(A) / 1048576.0,
		  getsize_sparsematrix(Af) / 1048576.0,
		  getsize_sparsematrix(A) / 1024.0 / dc[i]->ndof,
		  getsize_sparsematrix(Af) / 1024.0 / dc[i]->ndof,
		  A->nz, Af->nz);

    (void) printf("  Setting up Dirichlet data\n");
    xd = new_avector(dc[i]->nfix);
    assemble_tri2dp1_dirichlet_avector(dc[i], sin_solution, 0, xd);

    (void) printf("  Setting up right-hand side\n");
    b = new_avector(dc[i]->ndof);
    assemble_tri2dp1_functional_avector(dc[i], sin_rhs, 0, b);

    (void) printf("  Starting iteration\n");
    x = new_avector(dc[i]->ndof);
    random_real_avector(x);

    start_stopwatch(sw);
    iter = solve_cg_sparsematrix_avector(A, b, x, eps, steps);
    runtime = stop_stopwatch(sw);
    (void) printf("  %u CG iterations\n", iter);

    (void) printf("\n"
		  "  %.1f seconds\n"
		  "  %.1f seconds per step\n", runtime, runtime / iter);

    error = norml2_tri2dp1(dc[i], sin_solution, 0, x, xd);
    (void) printf("  rel. L^2 error %.4e     %s\n", error,
		  (IS_IN_RANGE(1.0e-4, error, 4.0e-3) ? "    okay" :
		   "NOT okay"));
    if (!IS_IN_RANGE(1.0e-4, error, 4.0e-3))
      problems++;



    del_avector(x);
    del_avector(b);
    del_avector(xd);
    del_sparsematrix(Af);
    del_sparsematrix(A);
  }

  (void) printf("----------------------------------------\n" "Cleaning up\n");
  for (i = 0; i <= L; i++)
    del_tri2dp1(dc[i]);
  freemem(dc);

  for (i = 0; i < L; i++)
    del_tri2dref(rf[i]);
  freemem(rf);

  for (i = 0; i <= L; i++)
    del_tri2d(gr[i]);
  freemem(gr);

  del_stopwatch(sw);

  uninit_h2lib();

  return problems;
}
