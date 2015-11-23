
#include <stdio.h>
#include "amatrix.h"
#include "krylov.h"

static void
jacobi(void *pdata, pavector r)
{
  pamatrix  A = (pamatrix) pdata;

  diagsolve_amatrix_avector(false, A, r);
}

static void
gauss_seidel(void *pdata, pavector r)
{
  pamatrix  A = (pamatrix) pdata;

  triangularsolve_amatrix_avector(true, false, false, A, r);
}

int
main()
{
  pamatrix  A;
  pamatrix  qr;
  pavector  b, x;
  pavector  r, p, a, q, rhat, tau;
  real      error;
  uint      n, kmax;
  uint      k, steps;
  uint      problems;

  problems = 0;

  (void) printf("Testing conjugate gradient method\n");
  n = 47;
  A = new_amatrix(n, n);
  b = new_avector(n);
  x = new_avector(n);

  r = new_avector(n);
  p = new_avector(n);
  a = new_avector(n);

  random_spd_amatrix(A, 1.0);
  random_avector(x);
  clear_avector(b);

  init_cg((addeval_t) addeval_amatrix_avector, A, b, x, r, p, a);
  steps = 0;
  error = norm2_avector(r);
  while (steps < 2 * n && error > 1e-8) {
    step_cg((addeval_t) addeval_amatrix_avector, A, b, x, r, p, a);
    steps++;
    error = norm2_avector(r);

    printf("  Step %u: residual %.2e\n", steps, error);
  }
  (void) printf("  %u steps, residual norm %.2e:", steps, error);
  if (steps <= n && error <= 1e-8)
    printf("    Okay\n");
  else {
    printf("    NOT Okay\n");
    problems++;
  }

  (void) printf("Testing preconditioned conjugate gradient method\n");
  q = new_avector(n);
  random_avector(x);

  init_pcg((addeval_t) addeval_amatrix_avector, A,
	   jacobi, A, b, x, r, q, p, a);
  steps = 0;
  error = norm2_avector(r);
  while (steps < 2 * n && error > 1e-8) {
    step_pcg((addeval_t) addeval_amatrix_avector, A,
	     jacobi, A, b, x, r, q, p, a);
    steps++;
    error = norm2_avector(r);

    printf("  Step %u: residual %.2e\n", steps, error);
  }
  (void) printf("  %u steps, residual norm %.2e:", steps, error);
  if (steps <= n && error <= 1e-8)
    printf("    Okay\n");
  else {
    printf("    NOT Okay\n");
    problems++;
  }

  del_avector(q);
  del_avector(a);
  del_avector(p);
  del_avector(r);

  (void) printf("Testing GMRES method\n");
  kmax = 10;
  rhat = new_avector(n);
  r = new_avector(n);
  q = new_avector(n);
  qr = new_amatrix(n, kmax);
  tau = new_avector(kmax);

  random_invertible_amatrix(A, 1.0);
  random_avector(x);
  clear_avector(b);

  init_gmres((addeval_t) addeval_amatrix_avector, A,
	     b, x, rhat, q, &k, qr, tau);
  steps = 0;
  while (steps < 2 * n && residualnorm_gmres(rhat, k) > 1e-8) {
    if (k + 1 >= kmax)
      finish_gmres((addeval_t) addeval_amatrix_avector, A,
		   b, x, rhat, q, &k, qr, tau);

    step_gmres((addeval_t) addeval_amatrix_avector, A,
	       b, x, rhat, q, &k, qr, tau);
    steps++;

    printf("  Step %u: dimension %u, residual %.2e\n",
	   steps, k, ABS(rhat->v[k]));
  }
  finish_gmres((addeval_t) addeval_amatrix_avector, A,
	       b, x, rhat, q, &k, qr, tau);
  copy_avector(b, r);
  addeval_amatrix_avector(-1.0, A, x, r);
  error = norm2_avector(r);
  (void) printf("  %u steps, residual norm %.2e (%.2e):",
		steps, error, residualnorm_gmres(rhat, k));
  if (steps <= n && error <= 1e-8 && residualnorm_gmres(rhat, k) <= 1e-8)
    printf("    Okay\n");
  else {
    printf("    NOT Okay\n");
    problems++;
  }

  (void) printf("Testing preconditioned GMRES method\n");
  random_avector(x);

  init_pgmres((addeval_t) addeval_amatrix_avector, A,
	      gauss_seidel, A, b, x, rhat, q, &k, qr, tau);
  steps = 0;
  while (steps < 2 * n && residualnorm_pgmres(rhat, k) > 1e-8) {
    if (k + 1 >= kmax)
      finish_pgmres((addeval_t) addeval_amatrix_avector, A,
		    gauss_seidel, A, b, x, rhat, q, &k, qr, tau);

    step_pgmres((addeval_t) addeval_amatrix_avector, A,
		gauss_seidel, A, b, x, rhat, q, &k, qr, tau);
    steps++;

    printf("  Step %u: dimension %u, preconditioned residual %.2e\n",
	   steps, k, ABS(rhat->v[k]));
  }
  finish_pgmres((addeval_t) addeval_amatrix_avector, A,
		gauss_seidel, A, b, x, rhat, q, &k, qr, tau);
  copy_avector(b, r);
  addeval_amatrix_avector(-1.0, A, x, r);
  gauss_seidel(A, r);
  error = norm2_avector(r);
  (void) printf("  %u steps, preconditioned residual norm %.2e (%.2e):",
		steps, error, residualnorm_pgmres(rhat, k));
  if (steps <= n && error <= 1e-8 && residualnorm_pgmres(rhat, k) < 1e-8)
    printf("    Okay\n");
  else {
    printf("    NOT Okay\n");
    problems++;
  }

  del_avector(tau);
  del_amatrix(qr);
  del_avector(q);
  del_avector(r);
  del_avector(rhat);
  del_avector(x);
  del_avector(b);
  del_amatrix(A);

  printf("----------------------------------------\n"
	 "  %u matrices and\n"
	 "  %u vectors still active\n"
	 "  %u errors found\n",
	 getactives_amatrix(), getactives_avector(), problems);

  return (int) problems;
}
