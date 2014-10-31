
/* ------------------------------------------------------------
   This is the file "krylov.c" of the H2Lib package.
   All rights reserved, Steffen Boerm 2010
   ------------------------------------------------------------ */

#include "krylov.h"
#include "factorizations.h"

/* ------------------------------------------------------------
   Standard conjugate gradient method
   ------------------------------------------------------------ */

void
init_cg(addeval_t addeval,
	void *matrix,
	pcavector b, pavector x, pavector r, pavector p, pavector a)
{
  (void) a;

  copy_avector(b, r);		/* r = b - A x */
  addeval(-1.0, matrix, x, r);

  copy_avector(r, p);		/* p = r */
}

void
step_cg(addeval_t addeval,
	void *matrix,
	pcavector b, pavector x, pavector r, pavector p, pavector a)
{
  field     gamma, lambda, mu;

  (void) b;

  clear_avector(a);		/* a = A p */
  addeval(1.0, matrix, p, a);

  gamma = dotprod_avector(p, a);	/* lambda = <p, r> / <p, a> */
  lambda = dotprod_avector(p, r) / gamma;

  add_avector(lambda, p, x);	/* x = x + lambda p */

  add_avector(-lambda, a, r);	/* r = r - lambda a */

  mu = dotprod_avector(r, a) / gamma;	/* p = r - mu p */
  scale_avector(-mu, p);
  add_avector(1.0, r, p);
}

real
evalfunctional_cg(addeval_t addeval,
		  void *matrix, pcavector b, pcavector x, pcavector r)
{
  (void) addeval;
  (void) matrix;

  return -0.5 * (dotprod_avector(r, x) + dotprod_avector(b, x));
}

/* ------------------------------------------------------------
   Preconditioned conjugate gradient method
   ------------------------------------------------------------ */

void
init_pcg(addeval_t addeval, void *matrix, prcd_t prcd, void *pdata, pcavector b,	/* Right-hand side */
	 pavector x,		/* Approximate solution */
	 pavector r,		/* Residual b-Ax */
	 pavector q,		/* Preconditioned residual */
	 pavector p,		/* Search direction */
	 pavector a)
{
  (void) p;
  (void) a;

  copy_avector(b, r);		/* r = b - A x */
  addeval(-1.0, matrix, x, r);

  copy_avector(r, q);		/* q = N r */
  if (prcd)
    prcd(pdata, q);

  copy_avector(q, p);		/* p = q */
}

void
step_pcg(addeval_t addeval, void *matrix, prcd_t prcd, void *pdata, pcavector b,	/* Right-hand side */
	 pavector x,		/* Approximate solution */
	 pavector r,		/* Residual b-Ax */
	 pavector q,		/* Preconditioned residual */
	 pavector p,		/* Search direction */
	 pavector a)
{
  field     gamma, lambda, mu;

  (void) b;

  clear_avector(a);		/* a = A p */
  addeval(1.0, matrix, p, a);

  gamma = dotprod_avector(p, a);	/* lambda = <p, r> / <p, a> */
  lambda = dotprod_avector(p, r) / gamma;

  add_avector(lambda, p, x);	/* x = x + lambda p */

  add_avector(-lambda, a, r);	/* r = r - lambda a */

  copy_avector(r, q);		/* q = N r */
  if (prcd)
    prcd(pdata, q);

  mu = dotprod_avector(q, a) / gamma;	/* p = r - mu p */
  scale_avector(-mu, p);
  add_avector(1.0, q, p);
}

/* ------------------------------------------------------------
   Biconjugate gradient method
   ------------------------------------------------------------ */

/* cf. Yousef Saad, Iterative Methods for Sparse Linear Systems,
   Section 7.3 */

void
init_bicg(addeval_t addeval, addeval_t addevaltrans, void *matrix, pcavector b,	/* Right-hand side */
	  pavector x,		/* Approximate solution */
	  pavector r,		/* Residual b-Ax */
	  pavector rt,		/* Adjoint residual */
	  pavector p,		/* Search direction */
	  pavector pt,		/* Adjoint search direction */
	  pavector a, pavector at)
{
  (void) a;
  (void) at;
  (void) addevaltrans;

  copy_avector(b, r);		/* r = b - A x */
  addeval(-1.0, matrix, x, r);

  copy_avector(r, rt);		/* r^* = r */

  copy_avector(r, p);		/* p = r */

  copy_avector(rt, pt);		/* p^* = r^* */
}

void
step_bicg(addeval_t addeval, addeval_t addevaltrans, void *matrix, pcavector b,	/* Right-hand side */
	  pavector x,		/* Approximate solution */
	  pavector r,		/* Residual b-Ax */
	  pavector rt,		/* Adjoint residual */
	  pavector p,		/* Search direction */
	  pavector pt,		/* Adjoint search direction */
	  pavector a, pavector at)
{
  field     alpha, beta, gamma, mu;

  (void) b;

  clear_avector(a);		/* a = A p */
  addeval(1.0, matrix, p, a);

  clear_avector(at);		/* a^* = A^* p^* */
  addevaltrans(1.0, matrix, pt, at);

  gamma = dotprod_avector(a, pt);
  mu = dotprod_avector(r, rt);
  alpha = mu / gamma;

  add_avector(alpha, p, x);	/* x = x + alpha p */

  add_avector(-alpha, a, r);	/* r = r - alpha a */

  add_avector(-alpha, at, rt);	/* r^* = r^* - alpha a^* */

  beta = dotprod_avector(r, rt) / mu;

  scale_avector(beta, p);	/* p = r + beta p */
  add_avector(1.0, r, p);

  scale_avector(beta, pt);	/* p^* = r^* + beta p^* */
  add_avector(1.0, rt, pt);
}

/* ------------------------------------------------------------
   Stabilized biconjugate gradient method
   ------------------------------------------------------------ */

/* cf. Yousef Saad, Iterative Methods for Sparse Linear Systems,
   Section 7.4.2
   Slight modification: the intermediate vector s is stored in r. */

void
init_bicgstab(addeval_t addeval, void *matrix, pcavector b,	/* Right-hand side */
	      pavector x,	/* Approximate solution */
	      pavector r,	/* Residual b-Ax */
	      pavector rt,	/* Adjoint residual */
	      pavector p,	/* Search direction */
	      pavector a, pavector as)
{
  (void) a;
  (void) as;

  copy_avector(b, r);		/* r = b - A x */
  addeval(-1.0, matrix, x, r);

  copy_avector(r, rt);		/* r^* = r */

  copy_avector(r, p);		/* p = r */
}

void
step_bicgstab(addeval_t addeval, void *matrix, pcavector b,	/* Right-hand side */
	      pavector x,	/* Approximate solution */
	      pavector r,	/* Residual b-Ax */
	      pavector rt,	/* Adjoint residual */
	      pavector p,	/* Search direction */
	      pavector a, pavector as)
{
  field     alpha, beta, omega, mu;

  (void) b;

  clear_avector(a);		/* a = A p */
  addeval(1.0, matrix, p, a);

  mu = dotprod_avector(r, rt);
  alpha = mu / dotprod_avector(a, rt);

  add_avector(-alpha, a, r);	/* r = r - alpha a */

  clear_avector(as);		/* as = A r */
  addeval(1.0, matrix, r, as);

  omega = dotprod_avector(as, r) / dotprod_avector(as, as);

  add_avector(alpha, p, x);	/* x = x + alpha p + omega s */
  add_avector(omega, r, x);

  add_avector(-omega, as, r);	/* r = r - omega as */

  beta = dotprod_avector(r, rt) / mu * alpha / omega;

  scale_avector(beta, p);	/* p = r + beta (p - omega a) */
  add_avector(1.0, r, p);
  add_avector(-beta * omega, a, p);
}

/* ------------------------------------------------------------
   GMRES method with Householder QR
   ------------------------------------------------------------ */

static    field
findapply_givens(pfield a, pfield b)
{
  field     aa = *a;
  field     bb = *b;
  field     s, c, t;
  field     rho;

  if (ABSSQR(bb) < ABSSQR(aa)) {
    t = bb / aa;
    c = 1.0 / REAL_SQRT(ABSSQR(t) + 1.0);
    s = t * c;

    rho = s;
  }
  else if (ABSSQR(aa) < ABSSQR(bb)) {
    t = aa / bb;
    s = 1.0 / REAL_SQRT(ABSSQR(t) + 1.0);
    c = t * s;

    rho = 1.0 / c;
  }
  else {
    c = 1.0;
    s = 0.0;

    rho = 1.0;
  }

  *a = c * aa + s * bb;
  *b = 0.0;

  return rho;
}

static void
apply_givens(field rho, pfield a, pfield b)
{
  field     aa = *a;
  field     bb = *b;
  field     s, c;

  if (ABSSQR(rho) < 1.0) {
    s = rho;
    c = REAL_SQRT(1.0 - ABSSQR(s));
  }
  else if (ABSSQR(rho) > 1.0) {
    c = 1.0 / rho;
    s = REAL_SQRT(1.0 - ABSSQR(c));
  }
  else {
    c = 1.0;
    s = 0.0;
  }

  *a = c * aa + s * bb;
  *b = -s * aa + c * bb;
}

void
init_gmres(addeval_t addeval, void *matrix, pcavector b,	/* Right-hand side */
	   pavector x,		/* Approximate solution */
	   pavector rhat,	/* Transformed residual */
	   pavector q,		/* Next search direction */
	   uint * kk,		/* Dimension of Krylov space */
	   pamatrix qr,		/* QR factorization of Krylov matrix */
	   pavector tau)
{				/* Scaling factors for elementary reflectors */
  avector   tmp1;
  amatrix   tmp2;
  pavector  r;
  pamatrix  qr_k;
  uint      kmax = qr->cols;

  assert(b->dim == x->dim);
  assert(b->dim == q->dim);
  assert(b->dim == qr->rows);
  assert(kmax <= tau->dim);

  if (kmax < 1)
    return;

  /* Residual r in the first column of qr */
  r = init_column_avector(&tmp1, qr, 0);
  copy_avector(b, r);
  addeval(-1.0, matrix, x, r);
  uninit_avector(r);

  /* Compute factorization */
  qr_k = init_sub_amatrix(&tmp2, qr, qr->rows, 0, 1, 0);

  qrdecomp_amatrix(qr_k, tau);

  /* Construct first orthogonal direction */
  clear_avector(q);
  q->v[0] = 1.0;
  qreval_amatrix_avector(false, qr_k, tau, q);
  uninit_amatrix(qr_k);

  /* Set up transformed residual */
  clear_avector(rhat);
  rhat->v[0] = qr->a[0];

  /* Set dimension */
  *kk = 0;
}

void
step_gmres(addeval_t addeval, void *matrix, pcavector b,	/* Right-hand side */
	   pavector x,		/* Approximate solution */
	   pavector rhat,	/* Transformed residual */
	   pavector q,		/* Next search direction */
	   uint * kk,		/* Dimension of Krylov space */
	   pamatrix qr,		/* QR factorization of Krylov matrix */
	   pavector tau)
{				/* Scaling factors for elementary reflectors */
  avector   tmp1, tmp2;
  amatrix   tmp3;
  pavector  a, tau_k;
  pamatrix  qr_k;
  field     rho;
  uint      k = *kk;
  uint      kmax = qr->cols;
  uint      i;

  (void) b;
  (void) x;

  if (k + 1 >= kmax)
    return;

  /* (k+1)-th Krylov vector A q in the (k+1)-th column of qr */
  a = init_column_avector(&tmp1, qr, k + 1);
  clear_avector(a);
  addeval(1.0, matrix, q, a);
  uninit_avector(a);

  /* Apply previous reflections */
  qr_k = init_sub_amatrix(&tmp3, qr, qr->rows, 0, k + 1, 0);
  qreval_amatrix_avector(true, qr_k, tau, a);
  uninit_amatrix(qr_k);

  /* Compute next reflection */
  qr_k = init_sub_amatrix(&tmp3, qr, qr->rows - (k + 1), k + 1, 1, k + 1);
  tau_k = init_sub_avector(&tmp2, tau, 1, k + 1);
  qrdecomp_amatrix(qr_k, tau_k);
  uninit_avector(tau_k);
  uninit_amatrix(qr_k);

  /* Construct next orthogonal direction */
  qr_k = init_sub_amatrix(&tmp3, qr, qr->rows, 0, k + 2, 0);
  clear_avector(q);
  q->v[k + 1] = 1.0;
  qreval_amatrix_avector(false, qr_k, tau, q);
  uninit_amatrix(qr_k);

  /* Apply preceding Givens rotations */
  for (i = 0; i < k; i++) {
    rho = qr->a[(i + 1) + (i + 1) * qr->ld];
    apply_givens(rho, qr->a + i + (k + 1) * qr->ld,
		 qr->a + (i + 1) + (k + 1) * qr->ld);
  }

  /* Eliminate subdiagonal */
  rho =
    findapply_givens(qr->a + k + (k + 1) * qr->ld,
		     qr->a + (k + 1) + (k + 1) * qr->ld);
  qr->a[(k + 1) + (k + 1) * qr->ld] = rho;
  apply_givens(rho, rhat->v + k, rhat->v + (k + 1));

  /* Increase dimension */
  *kk = k + 1;
}

void
finish_gmres(addeval_t addeval, void *matrix, pcavector b,	/* Right-hand side */
	     pavector x,	/* Approximate solution */
	     pavector rhat,	/* Transformed residual */
	     pavector q,	/* Next search direction */
	     uint * kk,		/* Dimension of Krylov space */
	     pamatrix qr,	/* QR factorization of Krylov matrix */
	     pavector tau)
{				/* Scaling factors for elementary reflectors */
  avector   tmp1;
  amatrix   tmp2;
  pamatrix  qr_k;
  pavector  rhat_k;
  uint      k = *kk;

  (void) addeval;
  (void) matrix;
  (void) b;
  (void) q;

  rhat_k = init_sub_avector(&tmp1, rhat, k, 0);
  qr_k = init_sub_amatrix(&tmp2, qr, k, 0, k, 1);

  triangularsolve_amatrix_avector(false, false, false, qr_k, rhat_k);

  uninit_amatrix(qr_k);
  uninit_avector(rhat_k);

  rhat_k = init_sub_avector(&tmp1, rhat, rhat->dim - k, k);
  clear_avector(rhat_k);
  uninit_avector(rhat_k);

  qr_k = init_sub_amatrix(&tmp2, qr, qr->rows, 0, k, 0);
  qreval_amatrix_avector(false, qr_k, tau, rhat);
  uninit_amatrix(qr_k);

  add_avector(1.0, rhat, x);

  init_gmres(addeval, matrix, b, x, rhat, q, kk, qr, tau);
}

/* ------------------------------------------------------------
   Preconditioned GMRES method with Householder QR
   ------------------------------------------------------------ */

void
init_pgmres(addeval_t addeval, void *matrix, prcd_t prcd, void *pdata, pcavector b,	/* Right-hand side */
	    pavector x,		/* Approximate solution */
	    pavector rhat,	/* Transformed residual */
	    pavector q,		/* Next search direction */
	    uint * kk,		/* Dimension of Krylov space */
	    pamatrix qr,	/* QR factorization of Krylov matrix */
	    pavector tau)
{				/* Scaling factors for elementary reflectors */
  avector   tmp1;
  amatrix   tmp2;
  pavector  r;
  pamatrix  qr_k;
  uint      kmax = qr->cols;

  assert(b->dim == x->dim);
  assert(b->dim == q->dim);
  assert(b->dim == qr->rows);
  assert(kmax <= tau->dim);

  if (kmax < 1)
    return;

  /* Residual r in the first column of qr */
  r = init_column_avector(&tmp1, qr, 0);
  copy_avector(b, r);
  addeval(-1.0, matrix, x, r);
  prcd(pdata, r);
  uninit_avector(r);

  /* Compute factorization */
  qr_k = init_sub_amatrix(&tmp2, qr, qr->rows, 0, 1, 0);

  qrdecomp_amatrix(qr_k, tau);

  /* Construct first orthogonal direction */
  clear_avector(q);
  q->v[0] = 1.0;
  qreval_amatrix_avector(false, qr_k, tau, q);
  uninit_amatrix(qr_k);

  /* Set up transformed residual */
  clear_avector(rhat);
  rhat->v[0] = qr->a[0];

  /* Set dimension */
  *kk = 0;
}

void
step_pgmres(addeval_t addeval, void *matrix, prcd_t prcd, void *pdata, pcavector b,	/* Right-hand side */
	    pavector x,		/* Approximate solution */
	    pavector rhat,	/* Transformed residual */
	    pavector q,		/* Next search direction */
	    uint * kk,		/* Dimension of Krylov space */
	    pamatrix qr,	/* QR factorization of Krylov matrix */
	    pavector tau)
{				/* Scaling factors for elementary reflectors */
  avector   tmp1, tmp2;
  amatrix   tmp3;
  pavector  a, tau_k;
  pamatrix  qr_k;
  field     rho;
  uint      k = *kk;
  uint      kmax = qr->cols;
  uint      i;

  (void) b;
  (void) x;

  if (k + 1 >= kmax)
    return;

  /* (k+1)-th Krylov vector A q in the (k+1)-th column of qr */
  a = init_column_avector(&tmp1, qr, k + 1);
  clear_avector(a);
  addeval(1.0, matrix, q, a);
  prcd(pdata, a);
  uninit_avector(a);

  /* Apply previous reflections */
  qr_k = init_sub_amatrix(&tmp3, qr, qr->rows, 0, k + 1, 0);
  qreval_amatrix_avector(true, qr_k, tau, a);
  uninit_amatrix(qr_k);

  /* Compute next reflection */
  qr_k = init_sub_amatrix(&tmp3, qr, qr->rows - (k + 1), k + 1, 1, k + 1);
  tau_k = init_sub_avector(&tmp2, tau, 1, k + 1);
  qrdecomp_amatrix(qr_k, tau_k);
  uninit_avector(tau_k);
  uninit_amatrix(qr_k);

  /* Construct next orthogonal direction */
  qr_k = init_sub_amatrix(&tmp3, qr, qr->rows, 0, k + 2, 0);
  clear_avector(q);
  q->v[k + 1] = 1.0;
  qreval_amatrix_avector(false, qr_k, tau, q);
  uninit_amatrix(qr_k);

  /* Apply preceding Givens rotations */
  for (i = 0; i < k; i++) {
    rho = qr->a[(i + 1) + (i + 1) * qr->ld];
    apply_givens(rho, qr->a + i + (k + 1) * qr->ld,
		 qr->a + (i + 1) + (k + 1) * qr->ld);
  }

  /* Eliminate subdiagonal */
  rho =
    findapply_givens(qr->a + k + (k + 1) * qr->ld,
		     qr->a + (k + 1) + (k + 1) * qr->ld);
  qr->a[(k + 1) + (k + 1) * qr->ld] = rho;
  apply_givens(rho, rhat->v + k, rhat->v + (k + 1));

  /* Increase dimension */
  *kk = k + 1;
}

void
finish_pgmres(addeval_t addeval, void *matrix, prcd_t prcd, void *pdata, pcavector b,	/* Right-hand side */
	      pavector x,	/* Approximate solution */
	      pavector rhat,	/* Transformed residual */
	      pavector q,	/* Next search direction */
	      uint * kk,	/* Dimension of Krylov space */
	      pamatrix qr,	/* QR factorization of Krylov matrix */
	      pavector tau)
{				/* Scaling factors for elementary reflectors */
  avector   tmp1;
  amatrix   tmp2;
  pamatrix  qr_k;
  pavector  rhat_k;
  uint      k = *kk;

  (void) addeval;
  (void) matrix;
  (void) b;
  (void) q;

  rhat_k = init_sub_avector(&tmp1, rhat, k, 0);
  qr_k = init_sub_amatrix(&tmp2, qr, k, 0, k, 1);

  triangularsolve_amatrix_avector(false, false, false, qr_k, rhat_k);

  uninit_amatrix(qr_k);
  uninit_avector(rhat_k);

  rhat_k = init_sub_avector(&tmp1, rhat, rhat->dim - k, k);
  clear_avector(rhat_k);
  uninit_avector(rhat_k);

  qr_k = init_sub_amatrix(&tmp2, qr, qr->rows, 0, k, 0);
  qreval_amatrix_avector(false, qr_k, tau, rhat);
  uninit_amatrix(qr_k);

  add_avector(1.0, rhat, x);

  init_pgmres(addeval, matrix, prcd, pdata, b, x, rhat, q, kk, qr, tau);
}
