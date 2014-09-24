
/* ------------------------------------------------------------
   This is the file "eigensolvers.c" of the H2Lib package.
   All rights reserved, Steffen Boerm 2009
   ------------------------------------------------------------ */

#include "eigensolvers.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "basic.h"
#include "factorizations.h"

/* ------------------------------------------------------------
   Constructors and destructors
   ------------------------------------------------------------ */

ptridiag
init_tridiag(ptridiag T, uint dim)
{
  T->d = allocfield(dim);
  T->u = (dim > 1 ? allocfield(dim-1) : NULL);
  T->l = (dim > 1 ? allocfield(dim-1) : NULL);
  T->dim = dim;
  T->owner = NULL;

  return T;
}

ptridiag
init_sub_tridiag(ptridiag T, ptridiag src, uint dim, uint off)
{
  assert(off + dim <= src->dim);

  T->d = src->d + off;
  T->u = (dim > 1 ? src->u + off : NULL);
  T->l = (dim > 1 ? src->l + off : NULL);
  T->dim = dim;
  T->owner = src;

  return T;
}

ptridiag
init_vec_tridiag(ptridiag T, pavector src, uint dim)
{
  assert(3*dim-2 <= src->dim);

  T->d = src->v;
  T->u = (dim > 1 ? src->v+dim : NULL);
  T->l = (dim > 1 ? src->v+2*dim-1 : NULL);
  T->dim = dim;
  T->owner = (ptridiag) src;

  return T;
}

void
uninit_tridiag(ptridiag T)
{
  if(T->owner == NULL) {
    freemem(T->d);
    if(T->dim > 1) {
      freemem(T->l);
      freemem(T->u);
    }
    else {
      assert(T->l == NULL);
      assert(T->u == NULL);
    }
  }
}

ptridiag
new_tridiag(uint dim)
{
  ptridiag T;

  T = (ptridiag) allocmem(sizeof(tridiag));

  init_tridiag(T, dim);

  return T;
}

void
del_tridiag(ptridiag T)
{
  uninit_tridiag(T);

  freemem(T);
}

void
copy_tridiag(pctridiag T, ptridiag Tcopy)
{
  uint dim = T->dim;
  pfield d = T->d;
  pfield u = T->u;
  pfield l = T->l;
  uint i;

  assert(Tcopy->dim >= T->dim);

  for(i=0; i<dim; i++)
    Tcopy->d[i] = d[i];
  for(i=0; i<dim-1; i++) {
    Tcopy->u[i] = u[i];
    Tcopy->l[i] = l[i];
  }
}

real
check_tridiag(pctridiag T, pcamatrix Ts)
{
  uint n = T->dim;
  pcfield d = T->d;
  pcfield l = T->l;
  pcfield u = T->u;
  pcfield a = Ts->a;
  uint lda = Ts->ld;
  uint i, j;
  field val;
  real norm;

  assert(Ts->rows >= n);
  assert(Ts->cols >= n);

  norm = 0.0;
  for(j=0; j<n; j++)
    for(i=0; i<n; i++) {
      val = a[i+j*lda];
      if(i == j)
	val -= d[i];
      else if(i == j+1)
	val -= l[j];
      else if(i+1 == j)
	val -= u[i];
      norm += ABSSQR(val);
    }

  return REAL_SQRT(norm);
}

real
check_lower_tridiag(pctridiag T, pcamatrix Ts)
{
  uint n = T->dim;
  pcfield d = T->d;
  pcfield l = T->l;
  pcfield a = Ts->a;
  uint lda = Ts->ld;
  uint i, j;
  field val;
  real norm;

  assert(Ts->rows >= n);
  assert(Ts->cols >= n);

  norm = 0.0;
  for(j=0; j<n; j++)
    for(i=0; i<n; i++) {
      val = a[i+j*lda];
      if(i == j)
	val -= d[i];
      else if(i == j+1)
	val -= l[j];

      if(ABS(val) > 1e-10)
	printf("%u %u %g\n", i, j, ABS(val));
      norm += ABSSQR(val);
    }

  return REAL_SQRT(norm);
}

/* ------------------------------------------------------------
   Givens rotation
   ------------------------------------------------------------ */

static void
givens(field a, field b, pfield c, pfield s)
{
  real norm, norm2;
  
  norm2 = ABSSQR(a) + ABSSQR(b);
  if(norm2 > 0.0) {
    norm = REAL_SQRT(norm2);
    *c = CONJ(a) / norm;
    *s = CONJ(b) / norm;
  }
  else {
    *c = 1.0;
    *s = 0.0;
  }
}

/* ------------------------------------------------------------
   One QR step for the symmetric tridiagonal matrix
   ------------------------------------------------------------ */

void
qrstep_tridiag(ptridiag T, field shift, pamatrix Q)
{
  field a0, a1, b0, b0t, x, y, carry;
  field c, s;
  pfield a = T->d;
  pfield b = T->l;
  pfield qa;
  uint ldq;
  uint n = T->dim;
  uint i, j;

  assert(Q == NULL || Q->cols == n);

  qa = (Q ? Q->a : NULL);
  ldq = (Q ? Q->ld : 0);

  /* Quick exit if trivial matrix */
  if(n < 2)
    return;

  /* Determine Givens rotation */
  givens(a[0]-shift, b[0], &c, &s);

  /* Apply to first two rows */
  a0   = c*a[0] + s*b[0];
  b0   = -CONJ(s)*a[0] + CONJ(c)*b[0];
  b0t  = c*CONJ(b[0]) + s*a[1];
  a1   = -CONJ(s)*CONJ(b[0]) + CONJ(c)*a[1];

  /* Apply to first two columns */
  a[0] = CONJ(c)*a0 + CONJ(s)*b0t;
  b[0] = CONJ(c)*b0 + CONJ(s)*a1;
  a[1] = -s*b0 + c*a1;

  /* Apply to first two columns of Q */
  if(Q)
    for(j=0; j<Q->rows; j++) {
      x = qa[j];
      y = qa[j+ldq];
      qa[j] =  CONJ(c)*x + CONJ(s)*y;
      qa[j+ldq] = -s*x + c*y;
    }

  /* Chase away coefficient at (2,0) */
  for(i=1; i<n-1; i++) {
    /* Apply Givens rotation to rows i, i+1 */
    carry = CONJ(s)*b[i];
    b[i] = c*b[i];

    /* Determine Givens rotation to eliminate carry... */
    givens(b[i-1], carry, &c, &s);
    
    /* ... and apply it */
    b[i-1] = c*b[i-1] + s*carry;

    /* Apply it also to remainder of rows i, i+1, */
    a0     = c*a[i] + s*b[i];
    b0     = -CONJ(s)*a[i] + CONJ(c)*b[i];
    b0t    = c*CONJ(b[i]) + s*a[i+1];
    a1     = -CONJ(s)*CONJ(b[i]) + CONJ(c)*a[i+1];

    /* apply it to the columns i, i+1, */
    a[i]   = CONJ(c)*a0 + CONJ(s)*b0t;
    b[i]   = CONJ(c)*b0 + CONJ(s)*a1;
    a[i+1] = -s*b0 + c*a1;

    /* and if necessary to the columns i, i+1 of Q */
    if(Q)
      for(j=0; j<Q->rows; j++) {
	x = qa[j+i*ldq];
	y = qa[j+(i+1)*ldq];
	qa[j+i*ldq] = CONJ(c)*x + CONJ(s)*y;
	qa[j+(i+1)*ldq] = -s*x + c*y;
      }
  }
}

/* ------------------------------------------------------------
   QR iteration for self-adjoint tridiagonal matrices
   ------------------------------------------------------------ */

struct _evp_sort_data {
  ptridiag T;
  pamatrix U, Vt;
};

static uint
evp_leq(uint i, uint j, void *data)
{
  struct _evp_sort_data *esd = (struct _evp_sort_data *) data;
  ptridiag T = esd->T;

  assert(i < T->dim);
  assert(j < T->dim);

  return (REAL(T->d[i]) <= REAL(T->d[j]));
}

static uint
evp_geq(uint i, uint j, void *data)
{
  struct _evp_sort_data *esd = (struct _evp_sort_data *) data;
  ptridiag T = esd->T;

  assert(i < T->dim);
  assert(j < T->dim);

  return (REAL(T->d[i]) >= REAL(T->d[j]));
}

static void
evp_swap(uint i, uint j, void *data)
{
  struct _evp_sort_data *esd = (struct _evp_sort_data *) data;
  ptridiag T = esd->T;
  pamatrix U = esd->U;
  pamatrix Vt = esd->Vt;
  uint ldu, ldv;
  field h;
  uint k;

  assert(i < T->dim);
  assert(j < T->dim);
  assert(U == NULL || U->cols >= T->dim);
  assert(Vt == NULL || Vt->rows >= T->dim);

  h = T->d[i];
  T->d[i] = T->d[j];
  T->d[j] = h;

  if(U) {
    ldu = U->ld;
    for(k=0; k<U->rows; k++) {
      h = U->a[k+i*ldu];
      U->a[k+i*ldu] = U->a[k+j*ldu];
      U->a[k+j*ldu] = h;
    }
  }

  if(Vt) {
    ldv = Vt->ld;
    for(k=0; k<Vt->cols; k++) {
      h = Vt->a[i+k*ldv];
      Vt->a[i+k*ldv] = Vt->a[j+k*ldv];
      Vt->a[j+k*ldv] = h;
    }
  }
}

uint
sb_muleig_tridiag(ptridiag T, pamatrix Q, uint maxiter)
{
  tridiag tmp;
  amatrix tmp2;
  struct _evp_sort_data esd;
  ptridiag Tsub;
  pamatrix Qsub;
  pfield a = T->d;
  pfield b = T->l;
  uint n = T->dim;
  field lambda;
  real d, m;
  uint off, dim, iter;

  /* Find first non-zero subdiagonal coefficient */
  off = 0;
  while(off < n-1 &&
	ABS(b[off]) < H2_QR_EPS * (ABS(a[off]) + ABS(a[off+1])))
    off++;

  iter = 0;
  while(off < n-1 && (maxiter == 0 || iter < maxiter)) {
    iter++;

    /* Find dimension of unreduced submatrix */
    dim = 2;
    while(dim < n-off &&
	  ABS(b[off+dim-1]) >= H2_QR_EPS * (ABS(a[off+dim-1]) + ABS(a[off+dim])))
      dim++;

    /* Determine Wilkinson shift */
    m = REAL(a[off+dim-2] + a[off+dim-1]) * 0.5;
    d = REAL(a[off+dim-2] - a[off+dim-1]) * 0.5;
    lambda = (d < 0.0 ?
	      m + REAL_SQRT(REAL_SQR(d) + ABSSQR(b[off+dim-2])) :
	      m - REAL_SQRT(REAL_SQR(d) + ABSSQR(b[off+dim-2])));

    /* Set up submatrices */
    Tsub = init_sub_tridiag(&tmp, T, dim, off);
    Qsub = (Q ? init_sub_amatrix(&tmp2, Q, Q->rows, 0, dim, off) : NULL);

    /* Perform QR step on submatrices */
    qrstep_tridiag(Tsub, lambda, Qsub);

    /* Clean up submatrices */
    if(Qsub)
      uninit_amatrix(Qsub);
    uninit_tridiag(Tsub);

    /* Find non-zero subdiagonal coefficient */
    while(off < n-1 &&
	  ABS(b[off]) < H2_QR_EPS*(ABS(a[off]) + ABS(a[off+1])))
      off++;
  }

  /* Sort eigenvalues */
  esd.T = T;
  esd.U = Q;
  esd.Vt = 0;
  heapsort(n, evp_leq, evp_swap, &esd);

  return iter;
}

#ifdef USE_BLAS
IMPORT_PREFIX void
dsteqr_(const char *compz,
	const unsigned *n,
	double *d,
	double *e,
	double *z,
	const unsigned *ldz,
	double *work,
	int *info);

IMPORT_PREFIX void
dstev_(const char *jobz,
       const unsigned *n,
       double *d,
       double *e,
       double *z,
       const unsigned *ldz,
       double *work,
       int *info);

uint
muleig_tridiag(ptridiag T, pamatrix Q)
{
  double *work;
  int info = 0;

  if(T->dim > 1) {
    if(Q) {
      work = allocfield(2 * T->dim + 2);
      dsteqr_("Vectors", &T->dim, T->d, T->l, Q->a, &Q->ld, work, &info);
      assert(info >= 0);
      freemem(work);
    }
    else {
      dstev_("No Vectors", &T->dim, T->d, T->l, NULL, &u_one, NULL, &info);
      assert(info >= 0);
    }
  }
  return (info != 0);
}
#else
uint
muleig_tridiag(ptridiag T, pamatrix Q)
{
  uint iter, maxiter;

  maxiter = 32 * T->dim;
  iter = sb_muleig_tridiag(T, Q, maxiter);

  return !(iter < maxiter);
}
#endif

uint
eig_tridiag(ptridiag T, pamatrix Q)
{
  if(Q)
    identity_amatrix(Q);

  return muleig_tridiag(T, Q);
}

/* ------------------------------------------------------------
   Hessenberg tridiagonalization
   ------------------------------------------------------------ */

void
sb_tridiagonalize_amatrix(pamatrix A,
			  ptridiag T, pamatrix Q)
{
  pfield aa;
  uint lda;
  pfield qa;
  uint ldq;
  pfield d, l, u;
  uint n;
  uint i, j, k;
  real norm2, norm;
  field first, alpha, beta, gamma;

  n = A->rows;

  assert(A->cols == n);
  assert(Q->rows == n);
  assert(Q->cols == n);
  assert(T->dim == n);

  aa = A->a;
  lda = A->ld;

  d = T->d;
  l = T->l;
  u = T->u;

  qa = (Q ? Q->a : 0);
  ldq = (Q ? Q->ld : 0);

  if(Q)
    identity_amatrix(Q);

  /* Quick exit if trivial matrix */
  if(n < 1)
    return;
  else if(n < 2) {
    d[0] = aa[0];
    return;
  }

  for(k=0; k<n-1; k++) {
    /* Compute norm of k-th column */
    norm2 = ABSSQR(aa[(k+1)+k*lda]);
    for(i=k+2; i<n; i++)
      norm2 += ABSSQR(aa[i+k*lda]);
    norm = REAL_SQRT(norm2);

    if(norm2 == 0.0) {
      d[k] = aa[k+k*lda];
      l[k] = 0.0;
      u[k] = 0.0;
    }
    else {
      /* Determine Householder reflection vector */
      first = aa[(k+1)+k*lda];
      alpha = -SIGN(first) * norm;
      gamma = first - alpha;

      /* Compute 2 / |v|^2 */
      beta = 1.0 / (norm2 - REAL(CONJ(alpha) * first));

      /* Normalize reflection vector */
      for(i=k+2; i<n; i++)
	aa[i+k*lda] /= gamma;
      beta *= ABSSQR(gamma);

      /* Compute k-th column */
      d[k] = aa[k+k*lda];
      l[k] = alpha;
      u[k] = CONJ(alpha);

      /* Update remaining columns */
      for(j=k+1; j<n; j++) {
	gamma = aa[(k+1)+j*lda];
	for(i=k+2; i<n; i++)
	  gamma += CONJ(aa[i+k*lda]) * aa[i+j*lda];

	gamma *= beta;

	aa[(k+1)+j*lda] -= gamma;
	for(i=k+2; i<n; i++)
	  aa[i+j*lda] -= gamma * aa[i+k*lda];
      }

      /* Update remaining rows */
      for(j=k+1; j<n; j++) {
	gamma = aa[j+(k+1)*lda];
	for(i=k+2; i<n; i++)
	  gamma += aa[i+k*lda] * aa[j+i*lda];

	gamma *= beta;

	aa[j+(k+1)*lda] -= gamma;
	for(i=k+2; i<n; i++)
	  aa[j+i*lda] -= gamma * CONJ(aa[i+k*lda]);
      }

      if(Q) {
	/* Update Q */
	for(j=0; j<n; j++) {
	  gamma = qa[j+(k+1)*ldq];
	  for(i=k+2; i<n; i++)
	    gamma += aa[i+k*ldq] * qa[j+i*ldq];

	  gamma *= beta;

	  qa[j+(k+1)*ldq] -= gamma;
	  for(i=k+2; i<n; i++)
	    qa[j+i*ldq] -= gamma * CONJ(aa[i+k*ldq]);
	}
      }
    }
  }
  T->d[k] = aa[k+k*lda];
}

#ifdef USE_BLAS
IMPORT_PREFIX void
dlarf_(const char *side,
       const unsigned *m,
       const unsigned *n,
       const double *v,
       const unsigned *incv,
       const double *tau,
       double *c,
       const unsigned *ldc,
       double *work);

void
tridiagonalize_amatrix(pamatrix A, pavector work,
		       ptridiag T, pamatrix Q)
{
  pfield aa;
  uint lda;
  pfield qa;
  uint ldq;
  pfield d, l, u;
  uint n, n1;
  uint i, k;
  real norm2, norm;
  field first, alpha, beta, gamma;

  n = A->rows;

  assert(A->cols == n);
  assert(Q->rows == n);
  assert(Q->cols == n);
  assert(T->dim == n);

  aa = A->a;
  lda = A->ld;

  d = T->d;
  l = T->l;
  u = T->u;

  qa = (Q ? Q->a : 0);
  ldq = (Q ? Q->ld : 0);

  if(Q)
    identity_amatrix(Q);

  /* Quick exit if trivial matrix */
  if(n < 1)
    return;
  else if(n < 2) {
    d[0] = aa[0];
    return;
  }

  for(k=0; k<n-1; k++) {
    /* Compute norm of k-th column */
    norm2 = ABSSQR(aa[(k+1)+k*lda]);
    for(i=k+2; i<n; i++)
      norm2 += ABSSQR(aa[i+k*lda]);
    norm = REAL_SQRT(norm2);

    if(norm2 == 0.0) {
      d[k] = aa[k+k*lda];
      l[k] = 0.0;
      u[k] = 0.0;
    }
    else {
      /* Determine Householder reflection vector */
      first = aa[(k+1)+k*lda];
      alpha = -SIGN(first) * norm;
      gamma = first - alpha;

      /* Compute 2 / |v|^2 */
      beta = 1.0 / (norm2 - REAL(CONJ(alpha) * first));

      /* Normalize reflection vector */
      for(i=k+2; i<n; i++)
	aa[i+k*lda] /= gamma;
      beta *= ABSSQR(gamma);

      /* Compute k-th column */
      d[k] = aa[k+k*lda];
      l[k] = alpha;
      u[k] = CONJ(alpha);

      /* Update remaining columns */
      first = aa[(k+1)+k*lda];
      aa[(k+1)+k*lda] = 1.0;
      n1 = n-k-1;
      dlarf_("Left",
	     &n1, &n1,
	     aa+(k+1)+k*lda, &u_one,
	     &beta,
	     aa+(k+1)+(k+1)*lda, &lda,
	     work->v);

      /* Update remaining rows */
      dlarf_("Right",
	     &n1, &n1,
	     aa+(k+1)+k*lda, &u_one,
	     &beta,
	     aa+(k+1)+(k+1)*lda, &lda,
	     work->v);

      if(Q) {
	/* Update Q */
	dlarf_("Right",
	       &n, &n1,
	       aa+(k+1)+k*lda, &u_one,
	       &beta,
	       qa+(k+1)*ldq, &ldq,
	       work->v);
      }
      aa[(k+1)+k*lda] = first;
    }
  }
  T->d[k] = aa[k+k*lda];
}
#else
void
tridiagonalize_amatrix(pamatrix A, pavector work,
		       ptridiag T, pamatrix Q)
{
  (void) work;

  sb_tridiagonalize_amatrix(A, T, Q);
}
#endif

/* ------------------------------------------------------------
   Self-adjoint eigenvalue problem
   ------------------------------------------------------------ */

uint
sb_eig_amatrix(pamatrix A, pavector lambda, pamatrix Q, uint maxiter)
{
  tridiag tmp;
  ptridiag T;
  uint n, iter;
  uint i;

  n = A->rows;

  assert(A->cols == n);
  assert(Q == 0 || Q->rows == n);
  assert(Q == 0 || Q->cols == n);

  /* Set up auxiliary matrix */
  T = init_tridiag(&tmp, n);

  /* Tridiagonalize A */
  sb_tridiagonalize_amatrix(A, T, Q);

  /* Solve tridiagonal eigenproblem */
  iter = sb_muleig_tridiag(T, Q, maxiter);

  /* Copy eigenvalues */
  if(lambda->dim < n)
    resize_avector(lambda, n);
  for(i=0; i<n; i++)
    lambda->v[i] = T->d[i];

  /* Clean up */
  uninit_tridiag(T);

  return iter;
}

#ifdef USE_BLAS
#if defined(THREADSAFE_LAPACK) || !defined(USE_OPENMP)
IMPORT_PREFIX void
dsyev_(const char *jobz,
       const char *uplo,
       const unsigned *n,
       double *a,
       const unsigned *lda,
       const double *w,
       double *work,
       const unsigned *lwork,
       int *info);

uint
eig_amatrix(pamatrix A, pavector lambda, pamatrix Q)
{
  double *work;
  unsigned lwork;
  uint n = A->rows;
  int info = 0;

  /* Quick exit if trivial matrix */
  if(n < 2) {
    if(Q)
      identity_amatrix(Q);
    return 0;
  }

  lwork = 12*n;
  work = allocfield(lwork);

  if(Q) {
    dsyev_("Vectors", "Lower triangle", &n, A->a, &A->ld, lambda->v,
	   work, &lwork, &info);
    copy_amatrix(false, A, Q);
  }
  else
    dsyev_("No vectors", "Lower triangle", &n, A->a, &A->ld, lambda->v,
	   work, &lwork, &info);

  freemem(work);

  return (info != 0);
}
#else
static uint
workaround_eig_amatrix(pamatrix A, pavector work,
		       pavector lambda, pamatrix Q)
{
  tridiag tmp;
  ptridiag T;
  uint n, info;
  uint i;

  n = A->rows;

  assert(A->cols == n);
  assert(Q == 0 || Q->rows == n);
  assert(Q == 0 || Q->cols == n);

  /* Set up auxiliary matrix */
  T = init_tridiag(&tmp, n);

  /* Tridiagonalize A */
  tridiagonalize_amatrix(A, work, T, Q);

  /* Solve tridiagonal eigenproblem */
  info = muleig_tridiag(T, Q);

  /* Copy eigenvalues */
  if(lambda->dim < n)
    resize_avector(lambda, n);
  for(i=0; i<n; i++)
    lambda->v[i] = T->d[i];

  /* Clean up */
  uninit_tridiag(T);

  return info;
}

uint
eig_amatrix(pamatrix A, pavector lambda, pamatrix Q)
{
  pavector work;
  avector worktmp;
  uint lwork;
  uint info;

  assert(A->rows == A->cols);

  lwork = A->rows;
  work = init_avector(&worktmp, lwork);
  info = workaround_eig_amatrix(A, work, lambda, Q);
  uninit_avector(work);

  return (info != 0);
}
#endif
#else
uint
eig_amatrix(pamatrix A, pavector lambda, pamatrix Q)
{
  uint iter, maxiter;

  assert(A->rows == A->cols);

  maxiter = 32 * A->rows;

  iter = sb_eig_amatrix(A, lambda, Q, maxiter);

  return !(iter < maxiter);
}
#endif

uint
geig_amatrix(pamatrix A, pamatrix M, pavector lambda, pamatrix Q)
{
  uint info;

  /* Compute Cholesky factorization M = L L^* */
  choldecomp_amatrix(M);

  /* Compute B = L^{-1} A L^{-*} */
  lowersolve_amatrix_amatrix(false, false, M, true, A);
  lowersolve_amatrix_amatrix(false, false, M, false, A);

  /* Solve eigenvalue problem for B */
  info = eig_amatrix(A, lambda, Q);

  /* Compute eigenvector matrix L^{-*} Q */
  uppersolve_amatrix_amatrix(false, true, M, false, Q);

  return info;
}

/* ------------------------------------------------------------
   One SVD step for a sub-bidiagonal matrix
   ------------------------------------------------------------ */

void
svdstep_tridiag(ptridiag T, field shift, pamatrix U, pamatrix Vt)
{
  field a0, a1, b0, b0t, x, y, carry;
  field c, s;
  pfield ua, va;
  uint ldu, ldv;
  pfield a = T->d;
  pfield b = T->l;
  uint n = T->dim;
  uint i, j;

  assert(U == NULL || U->cols == n);
  assert(Vt == NULL || Vt->rows == n);

  ua = (U ? U->a : NULL);
  ldu = (U ? U->ld : 0);

  va = (Vt ? Vt->a : NULL);
  ldv = (Vt ? Vt->ld : 0);

  /* Quick exit if trivial matrix */
  if(n < 2)
    return;

  /* Determine Givens rotation */
  givens(ABSSQR(a[0])-shift, b[0]*CONJ(a[0]), &c, &s);

  /* Apply to first two rows */
  a0  = c*a[0] + s*b[0];
  b0  = -CONJ(s)*a[0] + CONJ(c)*b[0];
  b0t = s*a[1];
  a1  = CONJ(c)*a[1];

  /* Apply to first two columns of U */
  if(U)
    for(j=0; j<U->rows; j++) {
      x = ua[j];
      y = ua[j+ldu];
      ua[j] = CONJ(c)*x + CONJ(s)*y;
      ua[j+ldu] = -s*x + c*y;
    }

  /* Determine Givens rotation to eliminate (0,1) */
  givens(CONJ(a0), CONJ(b0t), &c, &s);

  /* Apply to first two columns of rows 0 and 1 */
  a[0] = CONJ(c)*a0 + CONJ(s)*b0t;
  b[0] = CONJ(c)*b0 + CONJ(s)*a1;
  a[1] = -s*b0 + c*a1;

  /* Apply to first two rows of Vt */
  if(Vt)
    for(j=0; j<Vt->cols; j++) {
      x = va[j*ldv];
      y = va[1+j*ldv];
      va[j*ldv] = c*x + s*y;
      va[1+j*ldv] = -CONJ(s)*x + CONJ(c)*y;
    }
  
  /* Chase away coefficient at (2,0) */
  for(i=1; i<n-1; i++) {
    /* Apply Givens rotation to columns i-1, i of row i+1 */
    carry = CONJ(s)*b[i];
    b[i] = c*b[i];

    /* Determine Givens rotation to eliminate carry... */
    givens(b[i-1], carry, &c, &s);

    /* ... and apply it to rows i,i+1 of column i-1, ... */
    b[i-1] = c*b[i-1] + s*carry;

    /* ... rows i,i+1 of column i, ... */
    x = a[i];
    y = b[i];
    a[i] = c*x + s*y;
    b[i] = -CONJ(s)*x + CONJ(c)*y;

    /* ... and rows i,i+1 of column i+1 */
    carry = s*a[i+1];
    a[i+1] = CONJ(c)*a[i+1];

    /* Apply to columns i,i+1 of U */
    if(U)
      for(j=0; j<U->rows; j++) {
	x = ua[j+i*ldu];
	y = ua[j+(i+1)*ldu];
	ua[j+i*ldu] = CONJ(c)*x + CONJ(s)*y;
	ua[j+(i+1)*ldu] = -s*x + c*y;
      }

    /* Determine Givens rotation to eliminate carry... */
    givens(CONJ(a[i]), CONJ(carry), &c, &s);

    /* ... and apply it to columns i,i+1 of row i, ... */
    a[i] = CONJ(c)*a[i] + CONJ(s)*carry;

    /* ... and columns i,i+1 of row i+1 */
    x = b[i];
    y = a[i+1];
    b[i] = CONJ(c)*x + CONJ(s)*y;
    a[i+1] = -s*x + c*y;

    /* Apply to rows i,i+1 of Vt */
    if(Vt)
      for(j=0; j<Vt->cols; j++) {
	x = va[i+j*ldv];
	y = va[(i+1)+j*ldv];
	va[i+j*ldv] = c*x + s*y;
	va[(i+1)+j*ldv] = -CONJ(s)*x + CONJ(c)*y;
      }
  }
}

/* ------------------------------------------------------------
   Eliminate subdiagonal element for zero diagonal element
   ------------------------------------------------------------ */

static void
elimsubdiag(ptridiag T, pamatrix Vt)
{
  pfield a = T->d;
  pfield b = T->l;
  uint n = T->dim;
  pfield va;
  uint ldv;
  field c, s;
  field x, y, carry;
  uint i, j;

  assert(Vt == NULL || Vt->rows == n);

  va = (Vt ? Vt->a : NULL);
  ldv = (Vt ? Vt->ld : 0);

  /* Quick exit if trivial matrix */
  if(n < 2)
    return;

  /* Determine Givens rotation to eliminate subdiagonal... */
  givens(CONJ(a[1]), CONJ(b[0]), &c, &s);

  /* ... and apply it column-wise to the second row */
  a[1] = CONJ(c)*a[1] + CONJ(s)*b[0];
  b[0] = 0.0;

  /* ... and row-wise to Vt */
  if(Vt) {
    for(j=0; j<Vt->cols; j++) {
      x = va[j*ldv];
      y = va[1+j*ldv];
      va[j*ldv] = c*x + s*y;
      va[1+j*ldv] = -CONJ(s)*x + CONJ(c)*y;
    }
  }

  /* Chase away coefficient at (2,0) */
  for(i=1; i<n-1; i++) {
    /* Apply Givens rotation to row i+1 */
    y = b[i];
    b[i] = c*y;
    carry = CONJ(s)*y;

    /* Determine Givens rotation to eliminate (i+1,0)... */
    givens(CONJ(a[i+1]), CONJ(carry), &c, &s);

    /* ... and apply column-wise it to the (i+1)th row */
    a[i+1] = CONJ(c)*a[i+1] + CONJ(s)*carry;

    /* ... and row-wise to Vt */
    if(Vt) {
      for(j=0; j<Vt->cols; j++) {
	x = va[j*ldv];
	y = va[(i+1)+j*ldv];
	va[j*ldv] = c*x + s*y;
	va[(i+1)+j*ldv] = -CONJ(s)*x + CONJ(c)*y;
      }
    }
  }
}

/* ------------------------------------------------------------
   Singular value decomposition of a sub-bidiagonal matrix
   ------------------------------------------------------------ */

uint
sb_mulsvd_tridiag(ptridiag T, pamatrix U, pamatrix Vt, uint maxiter)
{
  tridiag tmp;
  amatrix tmp2, tmp3;
  struct _evp_sort_data esd;
  ptridiag Tsub;
  pamatrix Usub, Vsub;
  pfield a = T->d;
  pfield b = T->l;
  uint n = T->dim;
  pfield ua;
  uint ldu;
  real a0, a1, d, m, Tnorm;
  field b0, lambda;
  uint i, j, off, dim, iter;

  /* Quick exit if trivial matrix */
  if(T->dim < 2)
    return 0;

  ua = (U ? U->a : NULL);
  ldu = (U ? U->ld : 0);

  /* Compute maximum norm of T */
  Tnorm = ABS(a[0]);
  for(i=1; i<n; i++) {
    d = ABS(b[i-1]) + ABS(a[i]);
    if(d > Tnorm)
      Tnorm = d;
  }

  /* Find first non-zero subdiagonal coefficient */
  off = 0;
  while(off < n-1 &&
	ABS(b[off]) < H2_QR_EPS * (ABS(a[off]) + ABS(a[off+1])))
    off++;

  iter = 0;
  while(off < n-1 && (maxiter == 0 || iter < maxiter)) {
    iter++;

    /* Find dimension of unreduced submatrix */
    dim = 2;
    while(dim < n-off &&
	  ABS(b[off+dim-1]) >= H2_QR_EPS * (ABS(a[off+dim-1]) + ABS(a[off+dim])))
      dim++;

    /* Check for zero diagonal entries */
    i = 0;
    while(i < dim &&
	  ABS(a[off+i]) >= H2_QR_EPS * Tnorm)
      i++;

    if(i+1 < dim) {
      /* Eliminate subdiagonal entry */
      Tsub = init_sub_tridiag(&tmp, T, dim-i, off+i);
      Vsub = (Vt ? init_sub_amatrix(&tmp3, Vt, dim-i, off+i, Vt->cols, 0) : NULL);

      elimsubdiag(Tsub, Vsub);

      if(Vsub)
	uninit_amatrix(Vsub);
      uninit_tridiag(Tsub);
    }
    else {
      /* Set up submatrices */
      Tsub = init_sub_tridiag(&tmp, T, dim, off);
      Usub = (U ? init_sub_amatrix(&tmp2, U, U->rows, 0, dim, off) : NULL);
      Vsub = (Vt ? init_sub_amatrix(&tmp3, Vt, dim, off, Vt->cols, 0) : NULL);
      
      /* Determine Wilkinson shift for T T^* */
      a0 = ABSSQR(a[off+dim-2]);
      if(off+dim > 2)
	a0 += ABSSQR(b[off+dim-3]);
      b0 = b[off+dim-2] * CONJ(a[off+dim-2]);
      a1 = ABSSQR(b[off+dim-2]) + ABSSQR(a[off+dim-1]);

      m = (a0 + a1) * 0.5;
      d = (a0 - a1) * 0.5;
      lambda = (d < 0.0 ?
		m + REAL_SQRT(REAL_SQR(d) + ABSSQR(b0)) :
		m - REAL_SQRT(REAL_SQR(d) + ABSSQR(b0)));
      
      /* Perform SVD step on submatrices */
      svdstep_tridiag(Tsub, lambda, Usub, Vsub);

      /* Clean up submatrices */
      if(Vsub)
	uninit_amatrix(Vsub);
      if(Usub)
	uninit_amatrix(Usub);
      uninit_tridiag(Tsub);
    }

    /* Find non-zero subdiagonal coefficient */
    while(off < n-1 &&
	  ABS(b[off]) < H2_QR_EPS*(ABS(a[off]) + ABS(a[off+1])))
      off++;
  }

  /* Ensure positive signs */
  for(i=0; i<n; i++)
    if(ABS(a[i]) != a[i]) {
      assert(a[i] != 0.0);

      lambda = a[i] / ABS(a[i]);
      a[i] /= lambda;

      if(U)
	for(j=0; j<U->rows; j++)
	  ua[j+i*ldu] *= lambda;
    }

  /* Sort singular values */
  esd.T = T;
  esd.U = U;
  esd.Vt = Vt;
  heapsort(n, evp_geq, evp_swap, &esd);

  return iter;
}

#ifdef USE_BLAS
IMPORT_PREFIX void
dbdsqr_(const char *uplo,
	const unsigned *n,
	const unsigned *ncvt,
	const unsigned *nru,
	const unsigned *ncc,
	double *d,
	double *e,
	double *vt,
	const unsigned *ldvt,
	double *u,
	const unsigned *ldu,
	double *c,
	const unsigned *ldc,
	double *work,
	int *info);

uint
mulsvd_tridiag(ptridiag T, pamatrix U, pamatrix Vt)
{
  uint lwork;
  double *work;
  int info;

  if(T->dim < 1)
    return 0;

  lwork = 4 * T->dim;
  work = (double *) allocmem((size_t) sizeof(double) * lwork);
  dbdsqr_("Lower bidiagonal",
	  &T->dim,
	  (Vt ? &Vt->cols : &u_zero),
	  (U ? &U->rows : &u_zero),
	  &u_zero,
	  T->d, T->l,
	  (Vt ? Vt->a : 0), (Vt ? &Vt->ld : &u_one),
	  (U ? U->a : 0), (U ? &U->ld : &u_one),
	  0, &u_one,
	  work, &info);

  freemem(work);

  return (info != 0);
}
#else
uint
mulsvd_tridiag(ptridiag T, pamatrix U, pamatrix Vt)
{
  uint iter, maxiter;

  if(T->dim < 1)
    return 0;

  maxiter = 32 * T->dim;

  iter = sb_mulsvd_tridiag(T, U, Vt, maxiter);

  return !(iter < maxiter);
}
#endif

uint
svd_tridiag(ptridiag T, pamatrix U, pamatrix Vt)
{
  if(U)
    identity_amatrix(U);

  if(Vt)
    identity_amatrix(Vt);

  return mulsvd_tridiag(T, U, Vt);
}

/* ------------------------------------------------------------
   Bidiagonalize a matrix
   ------------------------------------------------------------ */

void
sb_bidiagonalize_amatrix(pamatrix A,
			 ptridiag T, pamatrix U, pamatrix Vt)
{
  pfield a, ua, va, d, l, tau;
  field alpha, beta, gamma, diag;
  real norm, norm2;
  uint rows, cols, lda, ldu, ldv;
  uint dim;
  uint i, j, k;

  rows = A->rows;
  cols = A->cols;
  dim = UINT_MIN(rows, cols);

  assert(T->dim == dim);
  assert(U == NULL || U->rows >= A->rows);
  assert(U == NULL || U->cols >= UINT_MIN(A->rows, A->cols));
  assert(Vt == NULL || Vt->cols >= A->cols);
  assert(Vt == NULL || Vt->rows >= UINT_MIN(A->rows, A->cols));

  a = A->a;
  lda = A->ld;

  d = T->d;
  l = T->l;
  tau = d;

  ua = (U ? U->a : NULL);
  ldu = (U ? U->ld : 0);

  va = (Vt ? Vt->a : NULL);
  ldv = (Vt ? Vt->ld : 0);

  if(U)
    identity_amatrix(U);
  if(Vt)
    identity_amatrix(Vt);


  if(rows > cols) {		/* Make A upper triangular */
    for(k=0; k<cols; k++) {
      norm2 = ABSSQR(a[k+k*lda]);
      for(i=k+1; i<rows; i++)
	norm2 += ABSSQR(a[i+k*lda]);
      norm = REAL_SQRT(norm2);

      beta = 0.0;
      if(norm2 > 0.0) {
	/* Determine Householder reflection vector v */
	diag = a[k+k*lda];
	alpha = -SIGN(diag) * norm;
	gamma = diag - alpha;

	/* Compute norm of v */
	beta = 1.0 / (norm2 - CONJ(alpha) * diag);

	/* Normalize reflection vector */
	for(i=k+1; i<rows; i++)
	  a[i+k*lda] /= gamma;
	beta *= ABSSQR(gamma);
	a[k+k*lda] = alpha;

	/* Update columns k+1,...,cols */
	for(j=k+1; j<cols; j++) {
	  gamma = a[k+j*lda];
	  for(i=k+1; i<rows; i++)
	    gamma += CONJ(a[i+k*lda]) * a[i+j*lda];

	  gamma *= beta;

	  a[k+j*lda] -= gamma;
	  for(i=k+1; i<rows; i++)
	    a[i+j*lda] -= gamma * a[i+k*lda];
	}
      }
      tau[k] = beta;
    }
    
    /* Apply reflections in reversed order to U */
    if(U)
      for(k=cols; k-->0; ) {
	beta = tau[k];
	if(beta != 0.0) {
	  for(j=0; j<U->cols; j++) {
	    gamma = ua[k+j*ldu];
	    for(i=k+1; i<rows; i++)
	      gamma += CONJ(a[i+k*lda]) * ua[i+j*ldu];
	    
	    gamma *= beta;
	    
	    ua[k+j*ldu] -= gamma;
	    for(i=k+1; i<rows; i++)
	      ua[i+j*ldu] -= gamma * a[i+k*lda];
	  }
	}
      }

    /* Clear subdiagonal entries */
    for(k=0; k<cols; k++)
      for(i=k+1; i<cols; i++)
	a[i+k*lda] = 0.0;

    /* Now A is quadratic */
    rows = cols;
  }
  else if(cols > rows) {	/* Make A lower triangular */
    for(k=0; k<rows; k++) {
      norm2 = ABSSQR(a[k+k*lda]);
      for(i=k+1; i<cols; i++)
	norm2 += ABSSQR(a[k+i*lda]);
      norm = REAL_SQRT(norm2);

      beta = 0.0;
      if(norm2 > 0.0) {
	/* Determine Householder reflection vector v */
	diag = a[k+k*lda];
	alpha = -SIGN(diag) * norm;
	gamma = diag - alpha;

	/* Compute norm of v */
	beta = 1.0 / (norm2 - CONJ(alpha) * diag);

	/* Normalize reflection vector */
	for(i=k+1; i<cols; i++)
	  a[k+i*lda] /= gamma;
	beta *= ABSSQR(gamma);
	a[k+k*lda] = alpha;

	/* Update rows k+1,...,rows */
	for(j=k+1; j<rows; j++) {
	  gamma = a[j+k*lda];
	  for(i=k+1; i<cols; i++)
	    gamma += CONJ(a[k+i*lda]) * a[j+i*lda];

	  gamma *= beta;

	  a[j+k*lda] -= gamma;
	  for(i=k+1; i<cols; i++)
	    a[j+i*lda] -= gamma * a[k+i*lda];
	}
      }
      tau[k] = beta;
    }
    
    /* Apply reflections in reversed order to Vt */
    if(Vt)
      for(k=rows; k-->0; ) {
	beta = tau[k];
	if(beta != 0.0) {
	  for(j=0; j<Vt->rows; j++) {
	    gamma = va[j+k*ldv];
	    for(i=k+1; i<cols; i++)
	      gamma += CONJ(a[k+i*lda]) * va[j+i*ldv];
	    
	    gamma *= beta;
	
	    va[j+k*ldv] -= gamma;
	    for(i=k+1; i<cols; i++)
	      va[j+i*ldv] -= gamma * a[k+i*lda];
	  }
	}
      }
    
    /* Clear superdiagonal entries */
    for(k=0; k<rows; k++)
      for(i=k+1; i<rows; i++)
	a[k+i*lda] = 0.0;

    /* Now A is quadratic */
    cols = rows;
  }

  assert(dim == rows);
  assert(dim == cols);

  /* Golub-Kahan bidiagonalization */
  for(k=0; k<dim; k++) {
    /* Eliminate (k,k+1) to (k,cols) by column reflections */
    norm2 = ABSSQR(a[k+k*lda]);
    for(i=k+1; i<cols; i++)
      norm2 += ABSSQR(a[k+i*lda]);
    norm = REAL_SQRT(norm2);

    if(norm2 > 0.0) {
      /* Determine Householder reflection vector v */
      diag = a[k+k*lda];
      alpha = -SIGN(diag) * norm;
      gamma = diag - alpha;

      /* Compute 2 / |v|^2 */
      beta = 1.0 / (norm2 - REAL(CONJ(alpha) * diag));

      /* Normalize reflection vector */
      for(i=k+1; i<cols; i++)
	a[k+i*lda] /= gamma;
      beta *= ABSSQR(gamma);

      /* Update rows k+1,...,rows */
      for(j=k+1; j<rows; j++) {
	gamma = a[j+k*lda];
	for(i=k+1; i<cols; i++)
	  gamma += CONJ(a[k+i*lda]) * a[j+i*lda];

	gamma *= beta;

	a[j+k*lda] -= gamma;
	for(i=k+1; i<cols; i++)
	  a[j+i*lda] -= gamma * a[k+i*lda];
      }

      /* Update columns of Vt */
      if(Vt)
	for(j=0; j<Vt->cols; j++) {
	  gamma = va[k+j*ldv];
	  for(i=k+1; i<cols; i++)
	    gamma += a[k+i*lda] * va[i+j*ldv];
	  
	  gamma *= beta;
	  
	  va[k+j*ldv] -= gamma;
	  for(i=k+1; i<cols; i++)
	    va[i+j*ldv] -= gamma * CONJ(a[k+i*lda]);
	}
      
      /* Store diagonal element */
      d[k] = alpha;
    }
    else
      d[k] = 0.0;

    /* Eliminate (k+2,k) to (rows,k) by row reflections */
    if(k+1 < rows) {
      norm2 = ABSSQR(a[(k+1)+k*lda]);
      for(i=k+2; i<rows; i++)
	norm2 += ABSSQR(a[i+k*lda]);
      norm = REAL_SQRT(norm2);

      if(norm2 > 0.0) {
	/* Determine Householder reflection vector v */
	diag = a[(k+1)+k*lda];
	alpha = -SIGN(diag) * norm;
	gamma = diag - alpha;
	
	/* Compute 2 / |v|^2 */
	beta = 1.0 / (norm2 - REAL(CONJ(alpha) * diag));
	
	/* Normalize reflection vector */
	for(i=k+2; i<rows; i++)
	  a[i+k*lda] /= gamma;
	beta *= ABSSQR(gamma);

	/* Update columns k+1,...,cols */
	for(j=k+1; j<cols; j++) {
	  gamma = a[(k+1)+j*lda];
	  for(i=k+2; i<rows; i++)
	    gamma += CONJ(a[i+k*lda]) * a[i+j*lda];
	  
	  gamma *= beta;
	  
	  a[(k+1)+j*lda] -= gamma;
	  for(i=k+2; i<rows; i++)
	    a[i+j*lda] -= gamma * a[i+k*lda];
	}
	
	/* Update rows of U */
	if(U)
	  for(j=0; j<U->rows; j++) {
	    gamma = ua[j+(k+1)*ldu];
	    for(i=k+2; i<rows; i++)
	      gamma += a[i+k*lda] * ua[j+i*ldu];
	    
	    gamma *= beta;
	    
	    ua[j+(k+1)*ldu] -= gamma;
	    for(i=k+2; i<rows; i++)
	      ua[j+i*ldu] -= gamma * CONJ(a[i+k*lda]);
	  }
	
	/* Store subdiagonal element */
	l[k] = alpha;
      }
      else
	l[k] = 0.0;
    }
  }
}

#ifdef USE_BLAS
IMPORT_PREFIX void
dlarf_(const char *side,
       const unsigned *m,
       const unsigned *n,
       const double *v,
       const unsigned *incv,
       const double *tau,
       double *c,
       const unsigned *ldc,
       double *work);

void
bidiagonalize_amatrix(pamatrix A, pavector work,
		      ptridiag T, pamatrix U, pamatrix Vt)
{
  pfield a, ua, va, d, l, tau;
  field alpha, beta, gamma, diag;
  real norm, norm2;
  uint rows, cols, lda, ldu, ldv;
  uint rows1, cols1;
  uint dim;
  uint i, k;

  rows = A->rows;
  cols = A->cols;
  dim = UINT_MIN(rows, cols);

  assert(T->dim == dim);
  assert(U == NULL || U->rows >= A->rows);
  assert(U == NULL || U->cols >= UINT_MIN(A->rows, A->cols));
  assert(Vt == NULL || Vt->cols >= A->cols);
  assert(Vt == NULL || Vt->rows >= UINT_MIN(A->rows, A->cols));
  assert(work->dim >= UINT_MAX(A->rows, A->cols));
  
  a = A->a;
  lda = A->ld;

  d = T->d;
  l = T->l;
  tau = d;

  ua = (U ? U->a : NULL);
  ldu = (U ? U->ld : 0);

  va = (Vt ? Vt->a : NULL);
  ldv = (Vt ? Vt->ld : 0);

  if(U)
    identity_amatrix(U);
  if(Vt)
    identity_amatrix(Vt);


  if(rows > cols) {		/* Make A upper triangular */
    assert(work->dim >= cols);

    for(k=0; k<cols; k++) {
      norm2 = ABSSQR(a[k+k*lda]);
      for(i=k+1; i<rows; i++)
	norm2 += ABSSQR(a[i+k*lda]);
      norm = REAL_SQRT(norm2);

      beta = 0.0;
      if(norm2 > 0.0) {
	/* Determine Householder reflection vector v */
	diag = a[k+k*lda];
	alpha = -SIGN(diag) * norm;
	gamma = diag - alpha;

	/* Compute norm of v */
	beta = 1.0 / (norm2 - CONJ(alpha) * diag);

	/* Normalize reflection vector */
	for(i=k+1; i<rows; i++)
	  a[i+k*lda] /= gamma;
	beta *= ABSSQR(gamma);
	a[k+k*lda] = 1.0;

	/* Update columns k+1,...,cols */
	rows1 = rows-k;
	cols1 = cols-k-1;
	dlarf_("Left",
	       &rows1, &cols1,
	       a+k+k*lda, &u_one,
	       &beta,
	       a+k+(k+1)*lda, &lda,
	       work->v);

	/* Set new diagonal entry */
	a[k+k*lda] = alpha;
      }
      tau[k] = beta;
    }
    
    /* Apply reflections in reversed order to U */
    if(U)
      for(k=cols; k-->0; ) {
	beta = tau[k];
	if(beta != 0.0) {
	  alpha = a[k+k*lda];
	  a[k+k*lda] = 1.0;

	  rows1 = rows-k;
	  cols1 = U->cols;
	  dlarf_("Left",
		 &rows1, &cols1,
		 a+k+k*lda, &u_one,
		 &beta,
		 ua+k, &ldu,
		 work->v);

	  a[k+k*lda] = alpha;
	}
      }

    /* Clear subdiagonal entries */
    for(k=0; k<cols; k++)
      for(i=k+1; i<cols; i++)
	a[i+k*lda] = 0.0;

    /* Now A is quadratic */
    rows = cols;
  }
  else if(cols > rows) {	/* Make A lower triangular */
    assert(work->dim >= rows);

    for(k=0; k<rows; k++) {
      norm2 = ABSSQR(a[k+k*lda]);
      for(i=k+1; i<cols; i++)
	norm2 += ABSSQR(a[k+i*lda]);
      norm = REAL_SQRT(norm2);

      beta = 0.0;
      if(norm2 > 0.0) {
	/* Determine Householder reflection vector v */
	diag = a[k+k*lda];
	alpha = -SIGN(diag) * norm;
	gamma = diag - alpha;

	/* Compute norm of v */
	beta = 1.0 / (norm2 - CONJ(alpha) * diag);

	/* Normalize reflection vector */
	for(i=k+1; i<cols; i++)
	  a[k+i*lda] /= gamma;
	beta *= ABSSQR(gamma);
	a[k+k*lda] = 1.0;

	/* Update rows k+1,...,rows */
	rows1 = rows-k-1;
	cols1 = cols-k;
	dlarf_("Right",
	       &rows1, &cols1,
	       a+k+k*lda, &lda,
	       &beta,
	       a+(k+1)+k*lda, &lda,
	       work->v);

	/* Set new diagonal entry */
	a[k+k*lda] = alpha;
      }
      tau[k] = beta;
    }
    
    /* Apply reflections in reversed order to Vt */
    if(Vt)
      for(k=rows; k-->0; ) {
	beta = tau[k];
	if(beta != 0.0) {
	  alpha = a[k+k*lda];
	  a[k+k*lda] = 1.0;

	  rows1 = Vt->rows;
	  cols1 = cols-k;
	  dlarf_("Right",
		 &rows1, &cols1,
		 a+k+k*lda, &lda,
		 &beta,
		 va+k*ldv, &ldv,
		 work->v);

	  a[k+k*lda] = alpha;
	}
      }
    
    /* Clear superdiagonal entries */
    for(k=0; k<rows; k++)
      for(i=k+1; i<rows; i++)
	a[k+i*lda] = 0.0;

    /* Now A is quadratic */
    cols = rows;
  }

  assert(dim == rows);
  assert(dim == cols);

  /* Golub-Kahan bidiagonalization */
  for(k=0; k<dim; k++) {
    /* Eliminate (k,k+1) to (k,cols) by column reflections */
    norm2 = ABSSQR(a[k+k*lda]);
    for(i=k+1; i<cols; i++)
      norm2 += ABSSQR(a[k+i*lda]);
    norm = REAL_SQRT(norm2);

    if(norm2 > 0.0) {
      /* Determine Householder reflection vector v */
      diag = a[k+k*lda];
      alpha = -SIGN(diag) * norm;
      gamma = diag - alpha;

      /* Compute 2 / |v|^2 */
      beta = 1.0 / (norm2 - REAL(CONJ(alpha) * diag));

      /* Normalize reflection vector */
      for(i=k+1; i<cols; i++)
	a[k+i*lda] /= gamma;
      beta *= ABSSQR(gamma);
      a[k+k*lda] = 1.0;

      /* Update rows k+1,...,rows */
      rows1 = rows-k-1;
      cols1 = cols-k;
      dlarf_("Right",
	     &rows1, &cols1,
	     a+k+k*lda, &lda,
	     &beta,
	     a+(k+1)+k*lda, &lda,
	     work->v);

      /* Update columns of Vt */
      if(Vt) {
	rows1 = cols-k;
	cols1 = Vt->cols;
	dlarf_("Left",
	       &rows1, &cols1,
	       a+k+k*lda, &lda,
	       &beta,
	       va+k, &ldv,
	       work->v);
      }
      
      /* Store diagonal element */
      d[k] = alpha;
    }
    else
      d[k] = 0.0;

    /* Eliminate (k+2,k) to (rows,k) by row reflections */
    if(k+1 < rows) {
      norm2 = ABSSQR(a[(k+1)+k*lda]);
      for(i=k+2; i<rows; i++)
	norm2 += ABSSQR(a[i+k*lda]);
      norm = REAL_SQRT(norm2);

      if(norm2 > 0.0) {
	/* Determine Householder reflection vector v */
	diag = a[(k+1)+k*lda];
	alpha = -SIGN(diag) * norm;
	gamma = diag - alpha;
	
	/* Compute 2 / |v|^2 */
	beta = 1.0 / (norm2 - REAL(CONJ(alpha) * diag));
	
	/* Normalize reflection vector */
	for(i=k+2; i<rows; i++)
	  a[i+k*lda] /= gamma;
	beta *= ABSSQR(gamma);
	a[(k+1)+k*lda] = 1.0;

	/* Update columns k+1,...,cols */
	rows1 = rows-k-1;
	cols1 = cols-k-1;
	dlarf_("Left",
	       &rows1, &cols1,
	       a+(k+1)+k*lda, &u_one,
	       &beta,
	       a+(k+1)+(k+1)*lda, &lda,
	       work->v);

	/* Update rows of U */
	if(U) {
	  rows1 = U->rows;
	  cols1 = rows-k-1;
	  dlarf_("Right",
		 &rows1, &cols1,
		 a+(k+1)+k*lda, &u_one,
		 &beta,
		 ua+(k+1)*ldu, &ldu,
		 work->v);
	}
	
	/* Store subdiagonal element */
	l[k] = alpha;
      }
      else
	l[k] = 0.0;
    }
  }
}
#else
void
bidiagonalize_amatrix(pamatrix A, pavector work,
		      ptridiag T, pamatrix U, pamatrix Vt)
{
  (void) work;

  sb_bidiagonalize_amatrix(A, T, U, Vt);
}
#endif

void
bidiagonalize_verified_amatrix(pamatrix A, pavector work,
			       ptridiag T, pamatrix U, pamatrix Vt)
{
  pamatrix Acopy, Utmp, Vttmp;
  amatrix tmp1, tmp2, tmp3;
  pavector Td, Tl;
  avector tmp4, tmp5;
  uint k;
  real error;

  k = UINT_MIN(A->rows, A->cols);

  Acopy = init_amatrix(&tmp1, A->rows, A->cols);
  copy_amatrix(false, A, Acopy);

  Utmp = init_amatrix(&tmp2, A->rows, k);
  Vttmp = init_amatrix(&tmp3, k, A->cols);

  bidiagonalize_amatrix(A, work, T, Utmp, Vttmp);

  if(U)
    copy_amatrix(false, Utmp, U);

  if(Vt)
    copy_amatrix(false, Vttmp, Vt);

  if(T->dim > 0) {
    Td = init_pointer_avector(&tmp4, T->d, T->dim);
    Tl = init_pointer_avector(&tmp5, T->l, T->dim-1);
    bidiagmul_amatrix(1.0, false, Utmp, Td, Tl);
    addmul_amatrix(-1.0, false, Utmp, false, Vttmp, Acopy);
    error = normfrob_amatrix(Acopy);
    if(error > 1e-12) {
      printf("  bidiag: %.4g\n", error);
      abort();
    }
    uninit_avector(Tl);
    uninit_avector(Td);
  }

  uninit_amatrix(Vttmp);
  uninit_amatrix(Utmp);
  uninit_amatrix(Acopy);
}

/* ------------------------------------------------------------
   Singular value decomposition of an arbitrary matrix
   ------------------------------------------------------------ */

uint
sb_svd_amatrix(pamatrix A,
	       pavector sigma, pamatrix U, pamatrix Vt, uint maxiter)
{
  tridiag tmp;
  ptridiag T;
  uint dim;
  uint iter;
  uint i;

  dim = UINT_MIN(A->rows, A->cols);

  /* Set up auxiliary matrix */
  T = init_tridiag(&tmp, dim);

  /* Bidiagonalize A */
  sb_bidiagonalize_amatrix(A, T, U, Vt);

  /* Compute SVD of bidiagonal matrix */
  iter = sb_mulsvd_tridiag(T, U, Vt, maxiter);

  /* Copy singular values */
  for(i=0; i<dim; i++)
    sigma->v[i] = T->d[i];

  /* Clean up */
  uninit_tridiag(T);

  return iter;
}

#ifdef USE_BLAS
#if defined(THREADSAFE_LAPACK) || !defined(USE_OPENMP)
IMPORT_PREFIX void
dgesvd_(const char *jobu,
	const char *jobvt,
	const unsigned *m,
	const unsigned *n,
	double *a,
	const unsigned *lda,
	double *s,
	double *u,
	const unsigned *ldu,
	double *vt,
	const unsigned *ldvt,
	double *work,
	const unsigned *lwork,
	int *info);

uint
svd_amatrix(pamatrix A, pavector sigma, pamatrix U, pamatrix Vt)
{
  double *work;
  unsigned lwork;
  int info = 0;

  if(A->rows > 0 && A->cols > 0) {
    lwork = 10 * UINT_MAX(A->rows, A->cols);
    work = allocfield(lwork);

    dgesvd_((U ? "Skinny left vectors" : "No left vectors"),
	    (Vt ? "Skinny right vectors" : "No right vectors"),
	    &A->rows, &A->cols,
	    A->a, &A->ld,
	    sigma->v,
	    (U ? U->a : NULL), (U ? &U->ld : &u_one),
	    (Vt ? Vt->a : NULL), (Vt ? &Vt->ld : &u_one),
	    work, &lwork,
	    &info);

    freemem(work);
  }

  return (info != 0);
}
#else
static uint
workaround_svd_amatrix(pamatrix A, pavector work,
		       pavector sigma, pamatrix U, pamatrix Vt)
{
  tridiag tmp;
  ptridiag T;
  uint dim;
  uint info;
  uint i;

  dim = UINT_MIN(A->rows, A->cols);

  /* Set up auxiliary matrix */
  T = init_tridiag(&tmp, dim);

  /* Bidiagonalize A */
  bidiagonalize_amatrix(A, work, T, U, Vt);

  /* Compute SVD of bidiagonal matrix */
  info = mulsvd_tridiag(T, U, Vt);

  /* Copy singular values */
  for(i=0; i<dim; i++)
    sigma->v[i] = T->d[i];

  /* Clean up */
  uninit_tridiag(T);

  return info;
}

uint
svd_amatrix(pamatrix A, pavector sigma, pamatrix U, pamatrix Vt)
{
  pavector work;
  avector worktmp;
  uint lwork;
  uint info;

  lwork = UINT_MAX(A->rows, A->cols) + 3*UINT_MIN(A->rows, A->cols);
  work = init_avector(&worktmp, lwork);
  info = workaround_svd_amatrix(A, work, sigma, U, Vt);
  uninit_avector(work);

  return (info != 0);
}
#endif
#else
uint
svd_amatrix(pamatrix A, pavector sigma, pamatrix U, pamatrix Vt)
{
  pavector work;
  avector worktmp;
  uint lwork;
  uint iter, maxiter;

  lwork = 3*UINT_MIN(A->rows, A->cols);
  work = init_avector(&worktmp, lwork);
  maxiter = 32 * A->rows;
  iter = sb_svd_amatrix(A, work, sigma, U, Vt, maxiter);
  uninit_avector(work);

  return !(iter < maxiter);
}
#endif

uint
svd_verified_amatrix(pamatrix A, pavector sigma, pamatrix U, pamatrix Vt)
{
  pamatrix Acopy;
  pamatrix Utmp, Vttmp, E;
  amatrix tmp1, tmp2, tmp3, tmp4;
  real errorU, errorVt, errorA;
  uint res;
  uint k;

  Acopy = init_amatrix(&tmp1, A->rows, A->cols);
  copy_amatrix(false, A, Acopy);

  k = UINT_MIN(A->rows, A->cols);
  Utmp = init_amatrix(&tmp2, A->rows, k);
  Vttmp = init_amatrix(&tmp3, k, A->cols);
  E = init_amatrix(&tmp4, k, k);

  res = svd_amatrix(A, sigma, Utmp, Vttmp);

  if(U)
    copy_amatrix(false, Utmp, U);

  if(Vt)
    copy_amatrix(false, Vttmp, Vt);

  errorU = check_ortho_amatrix(false, Utmp);
  errorVt = check_ortho_amatrix(true, Vttmp);

  diagmul_amatrix(1.0, false, Utmp, sigma);
  addmul_amatrix(-1.0, false, Utmp, false, Vttmp, Acopy);
  errorA = normfrob_amatrix(Acopy);

  if(errorU > 1e-12 || errorVt > 1e-12 || errorA > 1e-12) {
    (void) printf("  errorU %.4g, errorVt %.4g, errorA %.4g\n",
		  errorU, errorVt, errorA);
    abort();
  }

  uninit_amatrix(E);
  uninit_amatrix(Vttmp);
  uninit_amatrix(Utmp);
  uninit_amatrix(Acopy);

  return res;
}
