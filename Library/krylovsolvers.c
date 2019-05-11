/* ------------------------------------------------------------
 * This is the file "krylovsolvers.c" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2016
 * ------------------------------------------------------------ */

/** @file krylovsolvers.c
 *  @author Steffen B&ouml;rm */

#include "krylovsolvers.h"
#include "basic.h"
#include "krylov.h"

/* ------------------------------------------------------------
 * Conjugated gradients method
 * ------------------------------------------------------------ */

uint
solve_cg_avector(void *A, addeval_t addeval_A, pcavector b, pavector x,
		 real eps, uint maxiter)
{
  pavector  r, p, a;
  real      norm, error;
  uint      n, iter;

  n = x->dim;

  assert(b->dim == n);

  r = new_avector(n);
  p = new_avector(n);
  a = new_avector(n);

  norm = norm2_avector(b);

  init_cg(addeval_A, A, b, x, r, p, a);
  error = norm2_avector(r);

  iter = 0;
  while (error > eps * norm && iter + 1 != maxiter) {
    step_cg(addeval_A, A, b, x, r, p, a);
    error = norm2_avector(r);

    iter++;
  }

  del_avector(a);
  del_avector(p);
  del_avector(r);

  return iter;
}

uint
solve_cg_amatrix_avector(pcamatrix A, pcavector b, pavector x, real eps,
			 uint maxiter)
{
  return solve_cg_avector((void *) A, (addeval_t) addeval_amatrix_avector, b,
			  x, eps, maxiter);
}

uint
solve_cg_sparsematrix_avector(pcsparsematrix A, pcavector b, pavector x,
			      real eps, uint maxiter)
{
  return solve_cg_avector((void *) A,
			  (addeval_t) addeval_sparsematrix_avector, b, x, eps,
			  maxiter);
}

uint
solve_cg_hmatrix_avector(pchmatrix A, pcavector b, pavector x, real eps,
			 uint maxiter)
{
  return solve_cg_avector((void *) A, (addeval_t) addeval_hmatrix_avector, b,
			  x, eps, maxiter);
}

uint
solve_cg_h2matrix_avector(pch2matrix A, pcavector b, pavector x, real eps,
			  uint maxiter)
{
  return solve_cg_avector((void *) A, (addeval_t) addeval_h2matrix_avector, b,
			  x, eps, maxiter);
}

uint
solve_cg_dh2matrix_avector(pcdh2matrix A, pcavector b, pavector x,
			   real eps, uint maxiter)
{
  return solve_cg_avector((void *) A, (addeval_t) addeval_dh2matrix_avector,
			  b, x, eps, maxiter);
}

/* ------------------------------------------------------------
 * Preconditioned conjugated gradients method
 * ------------------------------------------------------------ */

uint
solve_pcg_avector(void *A, addeval_t addeval_A, prcd_t prcd, void *pdata,
		  pcavector b, pavector x, real eps, uint maxiter)
{
  pavector  r, q, p, a;
  real      norm, error;
  uint      n, iter;

  n = x->dim;

  assert(b->dim == n);

  r = new_avector(n);
  q = new_avector(n);
  p = new_avector(n);
  a = new_avector(n);

  norm = norm2_avector(b);

  init_pcg(addeval_A, A, prcd, pdata, b, x, r, q, p, a);
  error = norm2_avector(r);

  iter = 0;
  while (error > eps * norm && iter + 1 != maxiter) {
    step_pcg(addeval_A, A, prcd, pdata, b, x, r, q, p, a);
    error = norm2_avector(r);

    iter++;
  }

  del_avector(a);
  del_avector(p);
  del_avector(q);
  del_avector(r);

  return iter;
}

uint
solve_pcg_amatrix_avector(pcamatrix A, prcd_t prcd, void *pdata,
			  pcavector b, pavector x, real eps, uint maxiter)
{
  return solve_pcg_avector((void *) A, (addeval_t) addeval_amatrix_avector,
			   prcd, pdata, b, x, eps, maxiter);
}

uint
solve_pcg_sparsematrix_avector(pcsparsematrix A, prcd_t prcd, void *pdata,
			       pcavector b, pavector x, real eps,
			       uint maxiter)
{
  return solve_pcg_avector((void *) A,
			   (addeval_t) addeval_sparsematrix_avector, prcd,
			   pdata, b, x, eps, maxiter);
}

uint
solve_pcg_hmatrix_avector(pchmatrix A, prcd_t prcd, void *pdata,
			  pcavector b, pavector x, real eps, uint maxiter)
{
  return solve_pcg_avector((void *) A, (addeval_t) addeval_hmatrix_avector,
			   prcd, pdata, b, x, eps, maxiter);
}

uint
solve_pcg_h2matrix_avector(pch2matrix A, prcd_t prcd, void *pdata,
			   pcavector b, pavector x, real eps, uint maxiter)
{
  return solve_pcg_avector((void *) A, (addeval_t) addeval_h2matrix_avector,
			   prcd, pdata, b, x, eps, maxiter);
}

uint
solve_pcg_dh2matrix_avector(pcdh2matrix A, prcd_t prcd, void *pdata,
			    pcavector b, pavector x, real eps, uint maxiter)
{
  return solve_pcg_avector((void *) A, (addeval_t) addeval_dh2matrix_avector,
			   prcd, pdata, b, x, eps, maxiter);
}

/* ------------------------------------------------------------
 * Generalized minimal residual method
 * ------------------------------------------------------------ */

uint
solve_gmres_avector(void *A, addeval_t addeval_A, pcavector b, pavector x,
		    real eps, uint maxiter, uint kmax)
{
  pavector  rhat, r, q, tau;
  pamatrix  qr;
  real      norm, error;
  uint      n, iter, k;

  n = x->dim;

  assert(b->dim == n);

  rhat = new_avector(n);
  r = new_avector(n);
  q = new_avector(n);
  qr = new_amatrix(n, kmax);
  tau = new_avector(kmax);

  norm = norm2_avector(b);

  init_gmres(addeval_A, A, b, x, rhat, q, &k, qr, tau);
  error = residualnorm_gmres(rhat, k);

  iter = 0;
  while (error > eps * norm && iter + 1 != maxiter) {
    if (k + 1 >= kmax) {
      finish_gmres(addeval_A, A, b, x, rhat, q, &k, qr, tau);
    }

    step_gmres(addeval_A, A, b, x, rhat, q, &k, qr, tau);
    error = residualnorm_gmres(rhat, k);

    iter++;
  }
  finish_gmres(addeval_A, A, b, x, rhat, q, &k, qr, tau);

  del_avector(tau);
  del_amatrix(qr);
  del_avector(q);
  del_avector(r);
  del_avector(rhat);

  return iter;
}

uint
solve_gmres_amatrix_avector(pcamatrix A, pcavector b, pavector x, real eps,
			    uint maxiter, uint kmax)
{
  return solve_gmres_avector((void *) A, (addeval_t) addeval_amatrix_avector,
			     b, x, eps, maxiter, kmax);
}

uint
solve_gmres_sparsematrix_avector(pcsparsematrix A, pcavector b, pavector x,
				 real eps, uint maxiter, uint kmax)
{
  return solve_gmres_avector((void *) A,
			     (addeval_t) addeval_sparsematrix_avector, b, x,
			     eps, maxiter, kmax);
}

uint
solve_gmres_hmatrix_avector(pchmatrix A, pcavector b, pavector x, real eps,
			    uint maxiter, uint kmax)
{
  return solve_gmres_avector((void *) A, (addeval_t) addeval_hmatrix_avector,
			     b, x, eps, maxiter, kmax);
}

uint
solve_gmres_h2matrix_avector(pch2matrix A, pcavector b, pavector x,
			     real eps, uint maxiter, uint kmax)
{
  return solve_gmres_avector((void *) A, (addeval_t) addeval_h2matrix_avector,
			     b, x, eps, maxiter, kmax);
}

uint
solve_gmres_dh2matrix_avector(pcdh2matrix A, pcavector b, pavector x,
			      real eps, uint maxiter, uint kmax)
{
  return solve_gmres_avector((void *) A,
			     (addeval_t) addeval_dh2matrix_avector, b, x, eps,
			     maxiter, kmax);
}

/* ------------------------------------------------------------
 * Preconditioned generalized minimal residual method
 * ------------------------------------------------------------ */

uint
solve_pgmres_avector(void *A, addeval_t addeval_A, prcd_t prcd,
		     void *pdata, pcavector b, pavector x, real eps,
		     uint maxiter, uint kmax)
{
  pavector  rhat, r, q, tau;
  pamatrix  qr;
  real      norm, error;
  uint      n, iter, k;

  n = x->dim;

  assert(b->dim == n);

  rhat = new_avector(n);
  r = new_avector(n);
  q = new_avector(n);
  qr = new_amatrix(n, kmax);
  tau = new_avector(kmax);

  copy_avector(b, r);
  prcd(pdata, r);
  norm = norm2_avector(r);

  init_pgmres(addeval_A, A, prcd, pdata, b, x, rhat, q, &k, qr, tau);
  error = residualnorm_pgmres(rhat, k);

  iter = 0;
  while (error > eps * norm && iter + 1 != maxiter) {
    if (k + 1 >= kmax) {
      finish_pgmres(addeval_A, A, prcd, pdata, b, x, rhat, q, &k, qr, tau);
    }

    step_pgmres(addeval_A, A, prcd, pdata, b, x, rhat, q, &k, qr, tau);
    error = residualnorm_pgmres(rhat, k);

    iter++;
  }

  if (k > 0)
    finish_pgmres(addeval_A, A, prcd, pdata, b, x, rhat, q, &k, qr, tau);

  del_avector(tau);
  del_amatrix(qr);
  del_avector(q);
  del_avector(r);
  del_avector(rhat);

  return iter;
}

uint
solve_pgmres_amatrix_avector(pcamatrix A, prcd_t prcd, void *pdata,
			     pcavector b, pavector x, real eps, uint maxiter,
			     uint kmax)
{
  return solve_pgmres_avector((void *) A, (addeval_t) addeval_amatrix_avector,
			      prcd, pdata, b, x, eps, maxiter, kmax);
}

uint
solve_pgmres_sparsematrix_avector(pcsparsematrix A, prcd_t prcd,
				  void *pdata, pcavector b, pavector x,
				  real eps, uint maxiter, uint kmax)
{
  return solve_pgmres_avector((void *) A,
			      (addeval_t) addeval_sparsematrix_avector, prcd,
			      pdata, b, x, eps, maxiter, kmax);
}

uint
solve_pgmres_hmatrix_avector(pchmatrix A, prcd_t prcd, void *pdata,
			     pcavector b, pavector x, real eps, uint maxiter,
			     uint kmax)
{
  return solve_pgmres_avector((void *) A, (addeval_t) addeval_hmatrix_avector,
			      prcd, pdata, b, x, eps, maxiter, kmax);
}

uint
solve_pgmres_h2matrix_avector(pch2matrix A, prcd_t prcd, void *pdata,
			      pcavector b, pavector x, real eps, uint maxiter,
			      uint kmax)
{
  return solve_pgmres_avector((void *) A,
			      (addeval_t) addeval_h2matrix_avector, prcd,
			      pdata, b, x, eps, maxiter, kmax);
}

uint
solve_pgmres_dh2matrix_avector(pcdh2matrix A, prcd_t prcd, void *pdata,
			       pcavector b, pavector x, real eps,
			       uint maxiter, uint kmax)
{
  return solve_pgmres_avector((void *) A,
			      (addeval_t) addeval_dh2matrix_avector, prcd,
			      pdata, b, x, eps, maxiter, kmax);
}
