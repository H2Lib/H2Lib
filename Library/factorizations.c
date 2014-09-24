
/* ------------------------------------------------------------
   This is the file "factorizations.c" of the H2Lib package.
   All rights reserved, Steffen Boerm 2010
   ------------------------------------------------------------ */

#include "factorizations.h"

#include "settings.h"
#include "basic.h"

#include <math.h>
#include <stdio.h>

/* ------------------------------------------------------------
   Diagonal matrices
   ------------------------------------------------------------ */

void
diagsolve_amatrix_avector(bool atrans, pcamatrix a, pavector x)
{
  pcfield aa = a->a;
  uint lda = a->ld;
  pfield xv = x->v;
  uint n = UINT_MIN(a->rows, a->cols);
  uint i;

  if(atrans) {
    for(i=0; i<n; i++)
      xv[i] /= CONJ(aa[i+i*lda]);
  }
  else {
    for(i=0; i<n; i++)
      xv[i] /= aa[i+i*lda];
  }
}

#ifdef USE_BLAS
IMPORT_PREFIX void
dscal_(const unsigned *n,
       const double *alpha,
       double *x,
       const unsigned *incx);

void
diagsolve_amatrix_amatrix(bool atrans, pcamatrix a,
			  bool xtrans, pamatrix x)
{
  uint n = UINT_MIN(a->rows, a->cols);
  uint lda = a->ld;
  uint ldx = x->ld;
  double *aa = a->a;
  double *xa = x->a;
  double alpha;
  uint i;

  if(xtrans) {
    for(i=0; i<n; i++) {
      alpha = (atrans ?
	       1.0 / CONJ(aa[i+i*lda]) :
	       1.0 / aa[i+i*lda]);
      dscal_(&x->rows, &alpha, xa+i*ldx, &u_one);
    }
  }
  else {
    for(i=0; i<n; i++) {
      alpha = (atrans ?
	       1.0 / CONJ(aa[i+i*lda]) :
	       1.0 / aa[i+i*lda]);
      dscal_(&x->cols, &alpha, xa+i, &ldx);
    }
  }
}
#else
void
diagsolve_amatrix_amatrix(bool atrans, pcamatrix a,
			  bool xtrans, pamatrix x)
{
  uint n = UINT_MIN(a->rows, a->cols);
  uint lda = a->ld;
  uint ldx = x->ld;
  double *aa = a->a;
  double *xa = x->a;
  double alpha;
  uint i, j;

  if(xtrans) {
    for(i=0; i<n; i++) {
      alpha = (atrans ?
	       1.0 / CONJ(aa[i+i*lda]) :
	       1.0 / aa[i+i*lda]);
      for(j=0; j<x->rows; j++)
	xa[j+i*ldx] *= alpha;
    }
  }
  else {
    for(i=0; i<n; i++) {
      alpha = (atrans ?
	       1.0 / CONJ(aa[i+i*lda]) :
	       1.0 / aa[i+i*lda]);
      for(j=0; j<x->cols; j++)
	xa[i+j*ldx] *= alpha;
    }
  }
}
#endif

#ifdef USE_BLAS
IMPORT_PREFIX void
dscal_(const unsigned *n,
       const double *alpha,
       double *x,
       const unsigned *incx);

void
diageval_amatrix_amatrix(bool atrans, pcamatrix a,
			 bool xtrans, pamatrix x)
{
  uint n = UINT_MIN(a->rows, a->cols);
  uint lda = a->ld;
  uint ldx = x->ld;
  double *aa = a->a;
  double *xa = x->a;
  double alpha;
  uint i;

  if(xtrans) {
    for(i=0; i<n; i++) {
      alpha = (atrans ?
	       CONJ(aa[i+i*lda]) :
	       aa[i+i*lda]);
      dscal_(&x->rows, &alpha, xa+i*ldx, &u_one);
    }
  }
  else {
    for(i=0; i<n; i++) {
      alpha = (atrans ?
	       CONJ(aa[i+i*lda]) :
	       aa[i+i*lda]);
      dscal_(&x->cols, &alpha, xa+i, &ldx);
    }
  }
}
#else
void
diageval_amatrix_amatrix(bool atrans, pcamatrix a,
			 bool xtrans, pamatrix x)
{
  uint n = UINT_MIN(a->rows, a->cols);
  uint lda = a->ld;
  uint ldx = x->ld;
  double *aa = a->a;
  double *xa = x->a;
  double alpha;
  uint i, j;

  if(xtrans) {
    for(i=0; i<n; i++) {
      alpha = (atrans ?
	       CONJ(aa[i+i*lda]) :
	       aa[i+i*lda]);
      for(j=0; j<x->rows; j++)
	xa[j+i*ldx] *= alpha;
    }
  }
  else {
    for(i=0; i<n; i++) {
      alpha = (atrans ?
	       CONJ(aa[i+i*lda]) :
	       aa[i+i*lda]);
      for(j=0; j<x->cols; j++)
	xa[i+j*ldx] *= alpha;
    }
  }
}
#endif

/* ------------------------------------------------------------
   Triangular matrices
   ------------------------------------------------------------ */

#ifdef USE_BLAS
IMPORT_PREFIX void
dtrtrs_(const char *uplo,
	const char *trans,
	const char *diag,
	const unsigned *n,
	const unsigned *nrhs,
	const double *a,
	const unsigned *lda,
	double *b,
	const unsigned *ldb,
	int *info);

void
lowersolve_amatrix_avector(bool aunit, bool atrans, pcamatrix a, pavector x)
{
  uint n = UINT_MIN(a->rows, a->cols);
  int info;

  if(atrans) {
    if(aunit) {
      dtrtrs_("Upper", "Transposed", "Unit triangular",
	      &n, &u_one,
	      a->a, &a->ld, x->v, &x->dim,
	      &info);
      assert(info == 0);
    }
    else {
      dtrtrs_("Upper", "Transposed", "Non-unit triangular",
	      &n, &u_one,
	      a->a, &a->ld, x->v, &x->dim,
	      &info);
      assert(info == 0);
    }
  }
  else {
    if(aunit) {
      dtrtrs_("Lower", "Not transposed", "Unit triangular",
	      &n, &u_one,
	      a->a, &a->ld, x->v, &x->dim,
	      &info);
      assert(info == 0);
    }
    else {
      dtrtrs_("Lower", "Not transposed", "Non-unit triangular",
	      &n, &u_one,
	      a->a, &a->ld, x->v, &x->dim,
	      &info);
      assert(info == 0);
    }
  }
}

void
uppersolve_amatrix_avector(bool aunit, bool atrans, pcamatrix a, pavector x)
{
  uint n = UINT_MIN(a->rows, a->cols);
  int info;

  if(atrans) {
    if(aunit) {
      dtrtrs_("Lower", "Transposed", "Unit triangular",
	      &n, &u_one,
	      a->a, &a->ld, x->v, &x->dim,
	      &info);
      assert(info == 0);
    }
    else {
      dtrtrs_("Lower", "Transposed", "Non-unit triangular",
	      &n, &u_one,
	      a->a, &a->ld, x->v, &x->dim,
	      &info);
      assert(info == 0);
    }
  }
  else {
    if(aunit) {
      dtrtrs_("Upper", "Not transposed", "Unit triangular",
	      &n, &u_one,
	      a->a, &a->ld, x->v, &x->dim,
	      &info);
      assert(info == 0);
    }
    else {
      dtrtrs_("Upper", "Not transposed", "Non-unit triangular",
	      &n, &u_one,
	      a->a, &a->ld, x->v, &x->dim,
	      &info);
      assert(info == 0);
    }
  }
}
#else
void
lowersolve_amatrix_avector(bool aunit, bool atrans, pcamatrix a, pavector x)
{
  pcfield aa = a->a;
  uint lda = a->ld;
  pfield xv = x->v;
  uint n = UINT_MIN(a->rows, a->cols);
  field newval;
  uint i, j;

  if(atrans) {
    for(j=0; j<n; j++) {
      newval = (aunit ? xv[j] : (xv[j] /= CONJ(aa[j+j*lda])));
      for(i=j+1; i<n; i++)
	xv[i] -= CONJ(aa[j+i*lda]) * newval;
    }
  }
  else {
    for(j=0; j<n; j++) {
      newval = (aunit ? xv[j] : (xv[j] /= aa[j+j*lda]));
      for(i=j+1; i<n; i++)
	xv[i] -= aa[i+j*lda] * newval;
    }
  }
}

void
uppersolve_amatrix_avector(bool aunit, bool atrans, pcamatrix a, pavector x)
{
  pcfield aa = a->a;
  uint lda = a->ld;
  pfield xv = x->v;
  uint n = UINT_MIN(a->rows, a->cols);
  field newval;
  uint i, j;

  if(atrans) {
    for(j=n; j-->0; ) {
      newval = (aunit ? xv[j] : (xv[j] /= CONJ(aa[j+j*lda])));
      for(i=0; i<j; i++)
	xv[i] -= CONJ(aa[j+i*lda]) * newval;
    }
  }
  else {
    for(j=n; j-->0; ) {
      newval = (aunit ? xv[j] : (xv[j] /= aa[j+j*lda]));
      for(i=0; i<j; i++)
	xv[i] -= aa[i+j*lda] * newval;
    }
  }
}
#endif

void
triangularsolve_amatrix_avector(bool alower, bool atrans, bool aunit,
			pcamatrix a, pavector x)
{
  if(alower) {
    if(atrans)
      uppersolve_amatrix_avector(aunit, atrans, a, x);
    else
      lowersolve_amatrix_avector(aunit, atrans, a, x);
  }
  else {
    if(atrans)
      lowersolve_amatrix_avector(aunit, atrans, a, x);
    else
      uppersolve_amatrix_avector(aunit, atrans, a, x);
  }
}

#ifdef USE_BLAS
IMPORT_PREFIX void
dscal_(const unsigned *n,
       const double *alpha,
       double *x,
       const unsigned *incx);

IMPORT_PREFIX void
dger_(const unsigned *m,
      const unsigned *n,
      const double *alpha,
      const double *x,
      const unsigned *incx,
      const double *y,
      const unsigned *incy,
      double *a,
      const unsigned *lda);

void
lowersolve_amatrix_amatrix(bool aunit, bool atrans, pcamatrix a,
			   bool xtrans, pamatrix x)
{
  uint n = UINT_MIN(a->rows, a->cols);
  uint lda = a->ld;
  uint ldx = x->ld;
  double *aa = a->a;
  double *xa = x->a;
  uint k, nk;
  double alpha;

  if(atrans) {
    if(xtrans) {
      assert(x->cols >= n);

      for(k=0; k<n; k++) {
	nk = n-k-1;
	if(!aunit) {
	  alpha = 1.0 / aa[k+k*lda];
	  dscal_(&x->rows, &alpha, xa+k*ldx, &u_one);
	}
	dger_(&x->rows, &nk, &f_minusone,
	      xa+k*ldx, &u_one, aa+k+(k+1)*lda, &lda,
	      xa+(k+1)*ldx, &ldx);
      }
    }
    else {
      assert(x->rows >= n);

      for(k=0; k<n; k++) {
	nk = n-k-1;
	if(!aunit) {
	  alpha = 1.0 / aa[k+k*lda];
	  dscal_(&x->cols, &alpha, xa+k, &ldx);
	}
	dger_(&nk, &x->cols, &f_minusone,
	      aa+k+(k+1)*lda, &lda, xa+k, &ldx,
	      xa+(k+1), &ldx);
      }
    }
  }
  else {
    if(xtrans) {
      assert(x->cols >= n);

      for(k=0; k<n; k++) {
	nk = n-k-1;
	if(!aunit) {
	  alpha = 1.0 / aa[k+k*lda];
	  dscal_(&x->rows, &alpha, xa+k*ldx, &u_one);
	}
	dger_(&x->rows, &nk, &f_minusone,
	      xa+k*ldx, &u_one, aa+(k+1)+k*lda, &u_one,
	      xa+(k+1)*ldx, &ldx);
      }
    }
    else {
      assert(x->rows >= n);

      for(k=0; k<n; k++) {
	nk = n-k-1;
	if(!aunit) {
	  alpha = 1.0 / aa[k+k*lda];
	  dscal_(&x->cols, &alpha, xa+k, &ldx);
	}
	dger_(&nk, &x->cols, &f_minusone,
	      aa+(k+1)+k*lda, &u_one, xa+k, &ldx,
	      xa+(k+1), &ldx);
      }
    }
  }
}

void
uppersolve_amatrix_amatrix(bool aunit, bool atrans, pcamatrix a,
			   bool xtrans, pamatrix x)
{
  uint n = UINT_MIN(a->rows, a->cols);
  uint lda = a->ld;
  uint ldx = x->ld;
  double *aa = a->a;
  double *xa = x->a;
  uint k;
  double alpha;

  if(atrans) {
    if(xtrans) {
      assert(x->cols >= n);

      for(k=n; k-->0; ) {
	if(!aunit) {
	  alpha = 1.0 / aa[k+k*lda];
	  dscal_(&x->rows, &alpha, xa+k*ldx, &u_one);
	}
	dger_(&x->rows, &k, &f_minusone,
	      xa+k*ldx, &u_one, aa+k, &lda,
	      xa, &ldx);
      }
    }
    else {
      assert(x->rows >= n);

      for(k=n; k-->0; ) {
	if(!aunit) {
	  alpha = 1.0 / aa[k+k*lda];
	  dscal_(&x->cols, &alpha, xa+k, &ldx);
	}
	dger_(&k, &x->cols, &f_minusone,
	      aa+k, &lda, xa+k, &ldx,
	      xa, &ldx);
      }
    }
  }
  else {
    if(xtrans) {
      assert(x->cols >= n);

      for(k=n; k-->0; ) {
	if(!aunit) {
	  alpha = 1.0 / aa[k+k*lda];
	  dscal_(&x->rows, &alpha, xa+k*ldx, &u_one);
	}
	dger_(&x->rows, &k, &f_minusone,
	      xa+k*ldx, &u_one, aa+k*lda, &u_one,
	      xa, &ldx);
      }
    }
    else {
      assert(x->rows >= n);

      for(k=n; k-->0; ) {
	if(!aunit) {
	  alpha = 1.0 / aa[k+k*lda];
	  dscal_(&x->cols, &alpha, xa+k, &ldx);
	}
	dger_(&k, &x->cols, &f_minusone,
	      aa+k*lda, &u_one, xa+k, &ldx,
	      xa, &ldx);
      }
    }
  }
}
#else
void
lowersolve_amatrix_amatrix(bool aunit, bool atrans, pcamatrix a,
			   bool xtrans, pamatrix x)
{
  uint n = UINT_MIN(a->rows, a->cols);
  uint lda = a->ld;
  uint ldx = x->ld;
  pfield aa = a->a;
  pfield xa = x->a;
  uint i, j, k;
  field alpha;

  if(atrans) {
    if(xtrans) {
      assert(x->cols >= n);

      for(k=0; k<n; k++) {
	if(!aunit) {
	  alpha = 1.0 / aa[k+k*lda];
	  for(i=0; i<x->rows; i++)
	    xa[i+k*ldx] *= alpha;
	}
	for(i=0; i<x->rows; i++)
	  for(j=k+1; j<n; j++)
	    xa[i+j*ldx] -= xa[i+k*ldx] * aa[k+j*lda];
      }
    }
    else {
      assert(x->rows >= n);

      for(k=0; k<n; k++) {
	if(!aunit) {
	  alpha = 1.0 / CONJ(aa[k+k*lda]);
	  for(j=0; j<x->cols; j++)
	    xa[k+j*ldx] *= alpha;
	}
	for(i=k+1; i<n; i++)
	  for(j=0; j<x->cols; j++)
	    xa[i+j*ldx] -= CONJ(aa[k+i*lda]) * xa[k+j*ldx];
      }
    }
  }
  else {
    if(xtrans) {
      assert(x->cols >= n);

      for(k=0; k<n; k++) {
	if(!aunit) {
	  alpha = 1.0 / CONJ(aa[k+k*lda]);
	  for(i=0; i<x->rows; i++)
	    xa[i+k*ldx] *= alpha;
	}
	for(i=0; i<x->rows; i++)
	  for(j=k+1; j<n; j++)
	    xa[i+j*ldx] -= xa[i+k*ldx] * CONJ(aa[j+k*lda]);
      }
    }
    else {
      assert(x->rows >= n);

      for(k=0; k<n; k++) {
	if(!aunit) {
	  alpha = 1.0 / aa[k+k*lda];
	  for(j=0; j<x->cols; j++)
	    xa[k+j*ldx] *= alpha;
	}
	for(i=k+1; i<n; i++)
	  for(j=0; j<x->cols; j++)
	    xa[i+j*ldx] -= aa[i+k*lda] * xa[k+j*ldx];
      }
    }
  }
}

void
uppersolve_amatrix_amatrix(bool aunit, bool atrans, pcamatrix a,
			   bool xtrans, pamatrix x)
{
  uint n = UINT_MIN(a->rows, a->cols);
  uint lda = a->ld;
  uint ldx = x->ld;
  pfield aa = a->a;
  pfield xa = x->a;
  uint i, j, k;
  field alpha;

  if(atrans) {
    if(xtrans) {
      assert(x->cols >= n);

      for(k=n; k-->0; ) {
	if(!aunit) {
	  alpha = 1.0 / aa[k+k*lda];
	  for(i=0; i<x->rows; i++)
	    xa[i+k*ldx] *= alpha;
	}
	for(i=0; i<x->rows; i++)
	  for(j=0; j<k; j++)
	    xa[i+j*ldx] -= xa[i+k*ldx] * aa[k+j*lda];
      }
    }
    else {
      assert(x->rows >= n);

      for(k=n; k-->0; ) {
	if(!aunit) {
	  alpha = 1.0 / CONJ(aa[k+k*lda]);
	  for(j=0; j<x->cols; j++)
	    xa[k+j*ldx] *= alpha;
	}
	for(i=0; i<k; i++)
	  for(j=0; j<x->cols; j++)
	    xa[i+j*ldx] -= CONJ(aa[k+i*lda]) * xa[k+j*ldx];
      }
    }
  }
  else {
    if(xtrans) {
      assert(x->cols >= n);

      for(k=n; k-->0; ) {
	if(!aunit) {
	  alpha = 1.0 / CONJ(aa[k+k*lda]);
	  for(i=0; i<x->rows; i++)
	    xa[i+k*ldx] *= alpha;
	}
	for(i=0; i<x->rows; i++)
	  for(j=0; j<k; j++)
	    xa[i+j*ldx] -= xa[i+k*ldx] * CONJ(aa[j+k*lda]);
      }
    }
    else {
      assert(x->rows >= n);

      for(k=n; k-->0; ) {
	if(!aunit) {
	  alpha = 1.0 / aa[k+k*lda];
	  for(j=0; j<x->cols; j++)
	    xa[k+j*ldx] *= alpha;
	}
	for(i=0; i<k; i++)
	  for(j=0; j<x->cols; j++)
	    xa[i+j*ldx] -= aa[i+k*lda] * xa[k+j*ldx];
      }
    }
  }
}
#endif

void
triangularsolve_amatrix_amatrix(bool alower, bool atrans, bool aunit,
				pcamatrix a, bool xtrans, pamatrix x)
{
  if(alower) {
    if(atrans)
      uppersolve_amatrix_amatrix(aunit, atrans, a, xtrans, x);
    else
      lowersolve_amatrix_amatrix(aunit, atrans, a, xtrans, x);
  }
  else {
    if(atrans)
      lowersolve_amatrix_amatrix(aunit, atrans, a, xtrans, x);
    else
      uppersolve_amatrix_amatrix(aunit, atrans, a, xtrans, x);
  }
}

#ifdef USE_BLAS
IMPORT_PREFIX void
dtrmv_(const char *uplo,
       const char *trans,
       const char *diag,
       const unsigned *n,
       const double *a,
       const unsigned *lda,
       double *x,
       const unsigned *incx);

void
lowereval_amatrix_avector(pcamatrix a, bool aunit, pavector x)
{
  uint n = UINT_MIN(a->rows, a->cols);

  if(n == 0) 			/* Quick exit */
    return;

  dtrmv_("Lower", "Not transposed",
	 (aunit ? "Unit triangular" : "Non-unit triangular"),
	 &n, a->a, &a->ld, x->v, &u_one);
}

void
uppereval_amatrix_avector(pcamatrix a, bool aunit, pavector x)
{
  uint n = UINT_MIN(a->rows, a->cols);

  if(n == 0)			/* Quick exit */
    return;

  dtrmv_("Upper", "Not transposed",
	 (aunit ? "Unit triangular" : "Non-unit triangular"),
	 &n, a->a, &a->ld, x->v, &u_one);
}
#else
void
lowereval_amatrix_avector(pcamatrix a, bool aunit, pavector x)
{
  pcfield aa = a->a;
  uint lda = a->ld;
  pfield xv = x->v;
  uint n = UINT_MIN(a->rows, a->cols);
  field newval;
  uint i, j;

  assert(x->dim >= a->rows);
  assert(x->dim >= a->cols);

  for(i=a->rows; i-->n; ) {
    newval = 0.0;
    for(j=0; j<a->cols; j++)
      newval += aa[i+j*lda] * xv[j];
    xv[i] = newval;
  }

  for(; i-->0; ) {
    newval = (aunit ? xv[i] : xv[i] * aa[i+i*lda]);
    for(j=0; j<i; j++)
      newval += aa[i+j*lda] * xv[j];
    xv[i] = newval;
  }
}

void
uppereval_amatrix_avector(pcamatrix a, bool aunit, pavector x)
{
  pcfield aa = a->a;
  uint lda = a->ld;
  pfield xv = x->v;
  uint n = UINT_MIN(a->rows, a->cols);
  field newval;
  uint i, j;

  assert(x->dim >= a->rows);
  assert(x->dim >= a->cols);

  for(i=0; i<n; i++) {
    newval = (aunit ? xv[i] : xv[i] * aa[i+i*lda]);
    for(j=i+1; j<a->cols; j++)
      newval += aa[i+j*lda] * xv[j];
    xv[i] = newval;
  }

  for(; i<a->rows; i++)
    xv[i] = 0.0;
}
#endif

void
triangulareval_amatrix_avector(bool alower, bool aunit, pcamatrix a, pavector x)
{
  if(alower)
    lowereval_amatrix_avector(a, aunit, x);
  else
    uppereval_amatrix_avector(a, aunit, x);
}

#ifdef USE_BLAS
IMPORT_PREFIX void
dscal_(const unsigned *n,
       const double *alpha,
       double *x,
       const unsigned *incx);

void
lowereval_amatrix_amatrix(bool aunit, bool atrans, pcamatrix a,
			  bool xtrans, pamatrix x)
{
  uint n = UINT_MIN(a->rows, a->cols);
  uint lda = a->ld;
  uint ldx = x->ld;
  double *aa = a->a;
  double *xa = x->a;
  uint k, nk, l;
  double alpha;

  if(x->rows == 0 || x->cols == 0)
    return;

  if(atrans) {
    if(xtrans) {
      assert(x->cols >= a->cols);

      for(k=n; k<a->cols; k++)
	for(l=0; l<x->rows; l++)
	  xa[l+k*ldx] = 0.0;

      for(k=n; k-->0; ) {
	nk = a->cols-k-1;
	dger_(&x->rows, &nk, &f_one,
	      xa+k*ldx, &u_one, aa+k+(k+1)*lda, &lda,
	      xa+(k+1)*ldx, &ldx);
	if(!aunit) {
	  alpha = aa[k+k*lda];
	  dscal_(&x->rows, &alpha, xa+k*ldx, &u_one);
	}
      }
    }
    else {
      assert(x->rows >= a->cols);

      for(k=n; k<a->cols; k++)
	for(l=0; l<x->cols; l++)
	  xa[k+l*ldx] = 0.0;

      for(k=n; k-->0; ) {
	nk = a->cols-k-1;
	dger_(&nk, &x->cols, &f_one,
	      aa+k+(k+1)*lda, &lda, xa+k, &ldx,
	      xa+(k+1), &ldx);
	if(!aunit) {
	  alpha = aa[k+k*lda];
	  dscal_(&x->cols, &alpha, xa+k, &ldx);
	}
      }
    }
  }
  else {
    if(xtrans) {
      assert(x->cols >= a->rows);

      for(k=n; k<a->rows; k++)
	for(l=0; l<x->rows; l++)
	  xa[l+k*ldx] = 0.0;

      for(k=n; k-->0; ) {
	nk = a->rows-k-1;
	dger_(&x->rows, &nk, &f_one,
	      xa+k*ldx, &u_one, aa+(k+1)+k*lda, &u_one,
	      xa+(k+1)*ldx, &ldx);
	if(!aunit) {
	  alpha = aa[k+k*lda];
	  dscal_(&x->rows, &alpha, xa+k*ldx, &u_one);
	}
      }
    }
    else {
      assert(x->rows >= a->rows);

      for(k=n; k<a->rows; k++)
	for(l=0; l<x->cols; l++)
	  xa[k+l*ldx] = 0.0;

      for(k=n; k-->0; ) {
	nk = a->rows-k-1;
	dger_(&nk, &x->cols, &f_one,
	      aa+(k+1)+k*lda, &u_one, xa+k, &ldx,
	      xa+(k+1), &ldx);
	if(!aunit) {
	  alpha = aa[k+k*lda];
	  dscal_(&x->cols, &alpha, xa+k, &ldx);
	}
      }
    }
  }
}

void
uppereval_amatrix_amatrix(bool aunit, bool atrans, pcamatrix a,
			  bool xtrans, pamatrix x)
{
  uint n = UINT_MIN(a->rows, a->cols);
  uint lda = a->ld;
  uint ldx = x->ld;
  double *aa = a->a;
  double *xa = x->a;
  uint k, k2, l;
  double alpha;

  if(x->rows == 0 || x->cols == 0)
    return;

  if(atrans) {
    if(xtrans) {
      assert(x->cols >= a->cols);

      for(k=n; k<a->cols; k++)
	for(l=0; l<x->rows; l++)
	  xa[l+k*ldx] = 0.0;

      for(k=0; k<n; k++) {
	k2 = UINT_MIN(k, a->cols);
	dger_(&x->rows, &k2, &f_one,
	      xa+k*ldx, &u_one, aa+k, &lda,
	      xa, &ldx);
	if(!aunit) {
	  alpha = aa[k+k*lda];
	  dscal_(&x->rows, &alpha, xa+k*ldx, &u_one);
	}
      }

      for(; k<a->rows; k++)
	dger_(&x->rows, &a->cols, &f_one,
	      xa+k*ldx, &u_one, aa+k, &lda,
	      xa, &ldx);
    }
    else {
      assert(x->rows >= a->cols);

      for(k=n; k<a->cols; k++)
	for(l=0; l<x->cols; l++)
	  xa[k+l*ldx] = 0.0;

      for(k=0; k<n; k++) {
	k2 = UINT_MIN(k, a->cols);
	dger_(&k2, &x->cols, &f_one,
	      aa+k, &lda, xa+k, &ldx,
	      xa, &ldx);
	if(!aunit) {
	  alpha = aa[k+k*lda];
	  dscal_(&x->cols, &alpha, xa+k, &ldx);
	}
      }

      for(; k<a->rows; k++)
	dger_(&a->cols, &x->cols, &f_one,
	      aa+k, &lda, xa+k, &ldx,
	      xa, &ldx);
    }
  }
  else {
    if(xtrans) {
      assert(x->cols >= a->rows);

      for(k=n; k<a->rows; k++)
	for(l=0; l<x->rows; l++)
	  xa[l+k*ldx] = 0.0;

      for(k=0; k<n; k++) {
	k2 = UINT_MIN(k, a->rows);
	dger_(&x->rows, &k2, &f_one,
	      xa+k*ldx, &u_one, aa+k*lda, &u_one,
	      xa, &ldx);
	if(!aunit) {
	  alpha = aa[k+k*lda];
	  dscal_(&x->rows, &alpha, xa+k*ldx, &u_one);
	}
      }

      for(; k<a->cols; k++)
	dger_(&x->rows, &a->rows, &f_one,
	      xa+k*ldx, &u_one, aa+k*lda, &u_one,
	      xa, &ldx);
    }
    else {
      assert(x->rows >= a->rows);

      for(k=n; k<a->rows; k++)
	for(l=0; l<x->cols; l++)
	  xa[k+l*ldx] = 0.0;

      for(k=0; k<n; k++) {
	k2 = UINT_MIN(k, a->rows);
	dger_(&k2, &x->cols, &f_one,
	      aa+k*lda, &u_one, xa+k, &ldx,
	      xa, &ldx);
	if(!aunit) {
	  alpha = aa[k+k*lda];
	  dscal_(&x->cols, &alpha, xa+k, &ldx);
	}
      }

      for(; k<a->cols; k++)
	dger_(&a->rows, &x->cols, &f_one,
	      aa+k*lda, &u_one, xa+k, &ldx,
	      xa, &ldx);
    }
  }
}
#else
void
lowereval_amatrix_amatrix(bool aunit, bool atrans, pcamatrix a,
			  bool xtrans, pamatrix x)
{
  uint n = UINT_MIN(a->rows, a->cols);
  uint lda = a->ld;
  uint ldx = x->ld;
  pfield aa = a->a;
  pfield xa = x->a;
  uint i, j, k;
  field alpha;

  if(atrans) {
    if(xtrans) {
      assert(x->cols >= a->cols);

      for(k=n; k<a->cols; k++)
	for(i=0; i<x->rows; i++)
	  xa[i+k*ldx] = 0.0;

      for(k=n; k-->0; ) {
	for(i=0; i<x->rows; i++)
	  for(j=k+1; j<a->cols; j++)
	    xa[i+j*ldx] += xa[i+k*ldx] * aa[k+j*lda];
	if(!aunit) {
	  alpha = aa[k+k*lda];
	  for(i=0; i<x->rows; i++)
	    xa[i+k*ldx] *= alpha;
	}
      }
    }
    else {
      assert(x->rows >= a->cols);

      for(k=n; k<a->cols; k++)
	for(j=0; j<x->cols; j++)
	  xa[k+j*ldx] = 0.0;

      for(k=n; k-->0; ) {
	for(i=k+1; i<a->cols; i++)
	  for(j=0; j<x->cols; j++)
	    xa[i+j*ldx] += CONJ(aa[k+i*lda]) * xa[k+j*ldx];
	if(!aunit) {
	  alpha = CONJ(aa[k+k*lda]);
	  for(j=0; j<x->cols; j++)
	    xa[k+j*ldx] *= alpha;
	}
      }
    }
  }
  else {
    if(xtrans) {
      assert(x->cols >= a->rows);

      for(k=n; k<a->rows; k++)
	for(i=0; i<x->rows; i++)
	  xa[i+k*ldx] = 0.0;

      for(k=n; k-->0; ) {
	for(i=0; i<x->rows; i++)
	  for(j=k+1; j<a->rows; j++)
	    xa[i+j*ldx] += xa[i+k*ldx] * CONJ(aa[j+k*lda]);
	if(!aunit) {
	  alpha = CONJ(aa[k+k*lda]);
	  for(i=0; i<x->rows; i++)
	    xa[i+k*ldx] *= alpha;
	}
      }
    }
    else {
      assert(x->rows >= n);

      for(k=n; k<a->rows; k++)
	for(j=0; j<x->cols; j++)
	  xa[k+j*ldx] = 0.0;

      for(k=n; k-->0; ) {
	for(j=0; j<x->cols; j++)
	  for(i=k+1; i<a->rows; i++)
	    xa[i+j*ldx] += aa[i+k*lda] * xa[k+j*ldx];
	if(!aunit) {
	  alpha = aa[k+k*lda];
	  for(j=0; j<x->cols; j++)
	    xa[k+j*ldx] *= alpha;
	}
      }
    }
  }
}

void
uppereval_amatrix_amatrix(bool aunit, bool atrans, pcamatrix a,
			  bool xtrans, pamatrix x)
{
  uint n = UINT_MIN(a->rows, a->cols);
  uint lda = a->ld;
  uint ldx = x->ld;
  pfield aa = a->a;
  pfield xa = x->a;
  uint i, j, k, k2;
  field alpha;

  if(atrans) {
    if(xtrans) {
      assert(x->cols >= n);

      for(k=n; k<a->cols; k++)
	for(i=0; i<x->rows; i++)
	  xa[i+k*ldx] = 0.0;

      for(k=0; k<n; k++) {
	k2 = UINT_MIN(k, a->cols);
	for(j=0; j<k2; j++)
	  for(i=0; i<x->rows; i++)
	    xa[i+j*ldx] += xa[i+k*lda] * aa[k+j*lda];
	if(!aunit) {
	  alpha = aa[k+k*lda];
	  for(i=0; i<x->rows; i++)
	    xa[i+k*ldx] *= alpha;
	}
      }

      for(; k<a->rows; k++)
	for(j=0; j<a->cols; j++)
	  for(i=0; i<x->rows; i++)
	    xa[i+j*ldx] += xa[i+k*ldx] * aa[k+j*lda];
    }
    else {
      assert(x->rows >= a->rows);

      for(k=n; k<a->cols; k++)
	for(j=0; j<x->cols; j++)
	  xa[k+j*ldx] = 0.0;

      for(k=0; k<n; k++) {
	k2 = UINT_MIN(k, a->cols);
	for(i=0; i<k2; i++)
	  for(j=0; j<x->cols; j++)
	    xa[i+j*ldx] += CONJ(aa[k+i*lda]) * xa[k+j*ldx];
	if(!aunit) {
	  alpha = CONJ(aa[k+k*lda]);
	  for(j=0; j<x->cols; j++)
	    xa[k+j*ldx] *= alpha;
	}
      }

      for(; k<a->rows; k++)
	for(i=0; i<a->cols; i++)
	  for(j=0; j<x->cols; j++)
	    xa[i+j*ldx] += aa[k+i*lda] * xa[k+j*ldx];
    }
  }
  else {
    if(xtrans) {
      assert(x->cols >= a->cols);

      for(k=n; k<a->rows; k++)
	for(i=0; i<x->rows; i++)
	  xa[i+k*ldx] = 0.0;

      for(k=0; k<n; k++) {
	k2 = UINT_MIN(k, a->rows);
	for(i=0; i<x->rows; i++)
	  for(j=0; j<k2; j++)
	    xa[i+j*ldx] += xa[i+k*ldx] * CONJ(aa[j+k*lda]);
	if(!aunit) {
	  alpha = CONJ(aa[k+k*lda]);
	  for(i=0; i<x->rows; i++)
	    xa[i+k*ldx] *= alpha;
	}
      }
      
      for(; k<a->cols; k++)
	for(i=0; i<x->rows; i++)
	  for(j=0; j<a->rows; j++)
	    xa[i+j*ldx] += xa[i+k*ldx] * aa[j+k*lda];
    }
    else {
      assert(x->rows >= a->cols);

      for(k=n; k<a->rows; k++)
	for(j=0; j<x->cols; j++)
	  xa[k+j*ldx] = 0.0;

      for(k=0; k<n; k++) {
	k2 = UINT_MIN(k, a->rows);
	for(j=0; j<x->cols; j++)
	  for(i=0; i<k2; i++)
	    xa[i+j*ldx] += aa[i+k*lda] * xa[k+j*ldx];
	if(!aunit) {
	  alpha = aa[k+k*lda];
	  for(j=0; j<x->cols; j++)
	    xa[k+j*ldx] *= alpha;
	}
      }

      for(; k<a->cols; k++)
	for(j=0; j<x->cols; j++)
	  for(i=0; i<a->rows; i++)
	    xa[i+j*ldx] += aa[i+k*lda] * xa[k+j*ldx];
    }
  }
}
#endif

void
triangulareval_amatrix_amatrix(bool alower, bool atrans, bool aunit,
			       pcamatrix a, bool xtrans, pamatrix x)
{
  if(alower) {
    if(atrans)
      uppereval_amatrix_amatrix(aunit, atrans, a, xtrans, x);
    else
      lowereval_amatrix_amatrix(aunit, atrans, a, xtrans, x);
  }
  else {
    if(atrans)
      lowereval_amatrix_amatrix(aunit, atrans, a, xtrans, x);
    else
      uppereval_amatrix_amatrix(aunit, atrans, a, xtrans, x);
  }
}

void
triangularaddmul_amatrix_amatrix(field alpha,
				 bool alower, bool atrans, pcamatrix a,
				 bool blower, bool btrans, pcamatrix b,
				 pamatrix c)
{
  pcfield aa = a->a;
  uint lda = a->ld;
  pcfield ba = b->a;
  uint ldb = b->ld;
  pfield ca = c->a;
  uint ldc = c->ld;
  uint aoff, adim, ainc, boff, bdim, binc;
  uint j;
#ifndef USE_BLAS
  uint i, k;
#endif

  if(atrans) {
    assert(c->rows == a->cols);

    ainc = lda;
    lda = 1;

    if(btrans) {
      assert(a->rows == b->cols);
      assert(c->cols == b->rows);

      binc = 1;

      for(j=0; j<a->rows; j++) {
	if(alower) {		/* A^* upper triangular */
	  aoff = 0;
	  adim = UINT_MIN(j + 1, a->cols);
	}
	else {			/* A^* lower triangular */
	  aoff = j;
	  adim = a->cols - UINT_MIN(j, a->cols);
	}

	if(blower) {		/* B^* upper triangular */
	  boff = j;
	  bdim = b->rows - UINT_MIN(j, b->rows);
	}
	else {			/* B^* lower triangular */
	  boff = 0;
	  bdim = UINT_MIN(j + 1, b->rows);
	}

#ifdef USE_BLAS
	dger_(&adim, &bdim, &alpha,
	      aa + aoff * ainc + j*lda, &ainc,
	      ba + boff * binc + j*ldb, &binc,
	      ca + aoff + boff * ldc, &ldc);
#else
	for(k=0; k<bdim; k++)
	  for(i=0; i<adim; i++)
	    ca[(aoff+i) + (boff+k)*ldc] +=
	      alpha * CONJ(aa[(aoff+i)*ainc + j*lda]) * CONJ(ba[(boff+k)*binc + j*ldb]);
#endif
      }
    }
    else {
      assert(a->rows == b->rows);
      assert(c->cols == b->cols);

      binc = ldb;
      ldb = 1;

      for(j=0; j<a->rows; j++) {
	if(alower) {		/* A^* upper triangular */
	  aoff = 0;
	  adim = UINT_MIN(j + 1, a->cols);
	}
	else {			/* A^* lower triangular */
	  aoff = j;
	  adim = a->cols - UINT_MIN(j, a->cols);
	}

	if(blower) {		/* B lower triangular */
	  boff = 0;
	  bdim = UINT_MIN(j + 1, b->cols);
	}
	else {			/* B upper triangular */
	  boff = j;
	  bdim = b->cols - UINT_MIN(j, b->cols);
	}

#ifdef USE_BLAS
	dger_(&adim, &bdim, &alpha,
	      aa + aoff * ainc + j*lda, &ainc,
	      ba + boff * binc + j*ldb, &binc,
	      ca + aoff + boff * ldc, &ldc);
#else
	for(k=0; k<bdim; k++)
	  for(i=0; i<adim; i++)
	    ca[(aoff+i) + (boff+k)*ldc] +=
	      alpha * CONJ(aa[(aoff+i)*ainc + j*lda]) * ba[(boff+k)*binc + j*ldb];
#endif
      }
    }
  }
  else {
    assert(c->rows == a->rows);

    ainc = 1;

    if(btrans) {
      assert(a->cols == b->cols);
      assert(c->cols == b->rows);

      binc = 1;

      for(j=0; j<a->cols; j++) {
	if(alower) {		/* A lower triangular */
	  aoff = j;
	  adim = a->rows - UINT_MIN(j, a->rows);
	}
	else {			/* A upper triangular */
	  aoff = 0;
	  adim = UINT_MIN(j + 1, a->rows);
	}

	if(blower) {		/* B^* upper triangular */
	  boff = j;
	  bdim = b->rows - UINT_MIN(j, b->rows);
	}
	else {			/* B^* lower triangular */
	  boff = 0;
	  bdim = UINT_MIN(j + 1, b->rows);
	}

#ifdef USE_BLAS
	dger_(&adim, &bdim, &alpha,
	      aa + aoff * ainc + j*lda, &ainc,
	      ba + boff * binc + j*ldb, &binc,
	      ca + aoff + boff * ldc, &ldc);
#else
	for(k=0; k<bdim; k++)
	  for(i=0; i<adim; i++)
	    ca[(aoff+i) + (boff+k)*ldc] +=
	      alpha * aa[(aoff+i)*ainc + j*lda] * CONJ(ba[(boff+k)*binc + j*ldb]);
#endif
      }
    }
    else {
      assert(a->cols == b->rows);
      assert(c->cols == b->cols);

      binc = ldb;
      ldb = 1;

      for(j=0; j<a->cols; j++) {
	if(alower) {		/* A lower triangular */
	  aoff = j;
	  adim = a->rows - UINT_MIN(j, a->rows);
	}
	else {			/* A upper triangular */
	  aoff = 0;
	  adim = UINT_MIN(j + 1, a->rows);
	}

	if(blower) {		/* B lower triangular */
	  boff = 0;
	  bdim = UINT_MIN(j + 1, b->cols);
	}
	else {			/* B upper triangular */
	  boff = j;
	  bdim = b->cols - UINT_MIN(j, b->cols);
	}

#ifdef USE_BLAS
	dger_(&adim, &bdim, &alpha,
	      aa + aoff * ainc + j*lda, &ainc,
	      ba + boff * binc + j*ldb, &binc,
	      ca + aoff + boff * ldc, &ldc);
#else

	for(k=0; k<bdim; k++)
	  for(i=0; i<adim; i++)
	    ca[(aoff+i) + (boff+k)*ldc] +=
	      alpha * aa[(aoff+i)*ainc + j*lda] * ba[(boff+k)*binc + j*ldb];
#endif
      }
    }
  }
}

void
copy_lower_amatrix(pcamatrix a, bool aunit, pamatrix b)
{
  pfield aa = a->a;
  uint lda = a->ld;
  pfield ba = b->a;
  uint ldb = b->ld;
  uint rows;
  uint cols;
  uint i, j;

  rows = UINT_MIN(a->rows, b->rows);
  cols = UINT_MIN(a->cols, b->cols);

  for(j=0; j<cols; j++) {
    if(aunit) {
      for(i=0; i<rows && i<j; i++)
	ba[i+j*ldb] = 0.0;
      
      if(i==j) {
	ba[i+j*ldb] = 1.0;
	i++;
      }
    }
    else 
      for(i=0; i<rows && i<j; i++)
	ba[i+j*ldb] = 0.0;

    for(; i<rows; i++)
      ba[i+j*ldb] = aa[i+j*lda];
  }
  for(; j<b->cols; j++)
    for(i=0; i<rows; i++)
      ba[i+j*ldb] = 0.0;
}

void
copy_upper_amatrix(pcamatrix a, bool aunit, pamatrix b)
{
  pfield aa = a->a;
  uint lda = a->ld;
  pfield ba = b->a;
  uint ldb = b->ld;
  uint rows;
  uint cols;
  uint i, j;

  rows = UINT_MIN(a->rows, b->rows);
  cols = UINT_MIN(a->cols, b->cols);

  for(j=0; j<cols; j++) {
    if(aunit) {
      for(i=0; i<rows && i<j; i++)
	ba[i+j*ldb] = aa[i+j*lda];

      ba[i+j*ldb] = 1.0;
      i++;
    }
    else {
      for(i=0; i<rows && i<=j; i++)
	ba[i+j*ldb] = aa[i+j*lda];
    }

    for(; i<b->rows; i++)
      ba[i+j*ldb] = 0.0;
  }
}

/* ------------------------------------------------------------
   LR decomposition
   ------------------------------------------------------------ */

#ifdef USE_BLAS
IMPORT_PREFIX void
dscal_(const unsigned *n,
       const double *alpha,
       double *x,
       const unsigned *incx);

IMPORT_PREFIX void
dger_(const unsigned *m,
      const unsigned *n,
      const double *alpha,
      const double *x,
      const unsigned *incx,
      const double *y,
      const unsigned *incy,
      double *a,
      const unsigned *lda);

uint
lrdecomp_amatrix(pamatrix a)
{
  double *aa = a->a;
  uint lda = a->ld;
  uint n = a->rows;
  double alpha;
  uint i, n1;

  assert(n == a->cols);

  for(i=0; i<n-1; i++) {
    if(aa[i+i*lda] == 0.0)
      return i+1;

    alpha = 1.0 / aa[i+i*lda];

    n1 = n-i-1;
    dscal_(&n1, &alpha, aa+(i+1)+i*lda, &u_one);
    dger_(&n1, &n1,
	  &f_minusone,
	  aa+(i+1)+i*lda, &u_one,
	  aa+i+(i+1)*lda, &lda,
	  aa+(i+1)+(i+1)*lda, &lda);
  }

  if(aa[i+i*lda] == 0.0)
    return i+1;

  return 0;
}
#else
uint
lrdecomp_amatrix(pamatrix a)
{
  pfield aa = a->a;
  uint lda = a->ld;
  uint n = a->rows;
  field alpha;
  uint i, j, k;

  assert(n == a->cols);

  for(i=0; i<n-1; i++) {
    if(aa[i+i*lda] == 0.0)
      return i+1;

    alpha = 1.0 / aa[i+i*lda];

    for(j=i+1; j<n; j++) {
      aa[j+i*lda] *= alpha;
      for(k=i+1; k<n; k++)
	aa[j+k*lda] -= aa[j+i*lda] * aa[i+k*lda];
    }
  }

  if(aa[i+i*lda] == 0.0)
    return i+1;

  return 0;
}
#endif

void
lrsolve_amatrix_avector(pcamatrix a, pavector x)
{
  lowersolve_amatrix_avector(true, false, a, x);
  uppersolve_amatrix_avector(false, false, a, x);
}

/* ------------------------------------------------------------
   Cholesky decomposition
   ------------------------------------------------------------ */

#ifdef USE_BLAS
IMPORT_PREFIX void
dscal_(const unsigned *n,
       const double *alpha,
       double *x,
       const unsigned *incx);

IMPORT_PREFIX void
dsyr_(const char *uplo,
      const unsigned *n,
      const double *alpha,
      const double *x,
      const unsigned *incx,
      double *a,
      const unsigned *lda);

uint
choldecomp_amatrix(pamatrix a)
{
  double *aa = a->a;
  uint lda = a->ld;
  uint n = a->rows;
  double diag, alpha;
  uint i, n1;

  assert(n == a->cols);

  for(i=0; i<n-1; i++) {
    diag = REAL(aa[i+i*lda]);

    if(ABS(aa[i+i*lda] - diag) > 1e-12 ||
       diag <= 0.0)
      return i+1;

    aa[i+i*lda] = REAL_SQRT(aa[i+i*lda]);
    alpha = 1.0 / aa[i+i*lda];
    n1 = n-i-1;
    dscal_(&n1, &alpha, aa+(i+1)+i*lda, &u_one);
    dsyr_("Lower part", &n1,
	  &f_minusone,
	  aa+(i+1)+i*lda, &u_one,
	  aa+(i+1)+(i+1)*lda, &lda);
  }

  diag = REAL(aa[i+i*lda]);
  if(ABS(aa[i+i*lda] - diag) > 1e-12 ||
     diag <= 0.0)
    return i+1;
  
  aa[i+i*lda] = REAL_SQRT(diag);

  return 0;
}
#else
uint
choldecomp_amatrix(pamatrix a)
{
  pfield aa = a->a;
  uint lda = a->ld;
  uint n = a->rows;
  real diag, alpha;
  uint i, j, k;

  assert(n == a->cols);

  for(i=0; i<n-1; i++) {
    diag = REAL(aa[i+i*lda]);

    if(ABS(aa[i+i*lda] - diag) > 1e-12 ||
       diag <= 0.0)
      return i+1;

    aa[i+i*lda] = REAL_SQRT(diag);
    alpha = 1.0 / aa[i+i*lda];
    for(j=i+1; j<n; j++)
      aa[j+i*lda] *= alpha;

    for(j=i+1; j<n; j++)
      for(k=i+1; k<=j; k++)
	aa[j+k*lda] -= aa[j+i*lda] * CONJ(aa[k+i*lda]);
  }

  diag = REAL(aa[i+i*lda]);
  if(ABS(aa[i+i*lda] - diag) > 1e-12 ||
     diag <= 0.0)
    return i+1;
  
  aa[i+i*lda] = REAL_SQRT(diag);

  return 0;
}
#endif

void
cholsolve_amatrix_avector(pcamatrix a, pavector x)
{
  lowersolve_amatrix_avector(false, false, a, x);
  uppersolve_amatrix_avector(false, true, a, x);
}

/* ------------------------------------------------------------
   LDL^T decomposition
   ------------------------------------------------------------ */

#ifdef USE_BLAS
IMPORT_PREFIX void
dscal_(const unsigned *n,
       const double *alpha,
       double *x,
       const unsigned *incx);

IMPORT_PREFIX void
dsyr_(const char *uplo,
      const unsigned *n,
      const double *alpha,
      const double *x,
      const unsigned *incx,
      double *a,
      const unsigned *lda);

uint
ldltdecomp_amatrix(pamatrix a)
{
  double *aa = a->a;
  uint lda = a->ld;
  uint n = a->rows;
  double diag, alpha;
  uint i, n1;

  assert(n == a->cols);

  for(i=0; i<n-1; i++) {
    diag = REAL(aa[i+i*lda]);

    if(ABS(aa[i+i*lda] - diag) > 1e-12 ||
       diag == 0.0)
      return i+1;

    alpha = 1.0 / diag;
    n1 = n-i-1;
    dscal_(&n1, &alpha, aa+(i+1)+i*lda, &u_one);

    alpha = -diag;
    dsyr_("Lower part", &n1,
	  &alpha,
	  aa+(i+1)+i*lda, &u_one,
	  aa+(i+1)+(i+1)*lda, &lda);
  }

  diag = REAL(aa[i+i*lda]);
  if(ABS(aa[i+i*lda] - diag) > 1e-12 ||
     diag == 0.0)
    return i+1;

  return 0;
}
#else
uint
ldltdecomp_amatrix(pamatrix a)
{
  pfield aa = a->a;
  uint lda = a->ld;
  uint n = a->rows;
  real diag, alpha;
  uint i, j, k;

  assert(n == a->cols);

  for(i=0; i<n-1; i++) {
    diag = REAL(aa[i+i*lda]);

    if(ABS(aa[i+i*lda] - diag) > 1e-12 ||
       diag == 0.0)
      return i+1;

    alpha = 1.0 / diag;
    for(j=i+1; j<n; j++)
      aa[j+i*lda] *= alpha;

    for(j=i+1; j<n; j++)
      for(k=i+1; k<=j; k++)
	aa[j+k*lda] -= diag * aa[j+i*lda] * CONJ(aa[k+i*lda]);
  }

  diag = REAL(aa[i+i*lda]);
  if(ABS(aa[i+i*lda] - diag) > 1e-12 ||
     diag == 0.0)
    return i+1;

  return 0;
}
#endif

void
ldltsolve_amatrix_avector(pcamatrix a, pavector x)
{
  lowersolve_amatrix_avector(true, false, a, x);
  diagsolve_amatrix_avector(false, a, x);
  uppersolve_amatrix_avector(true, true, a, x);
}

/* ------------------------------------------------------------
   Orthogonal decompositions
   ------------------------------------------------------------ */

#ifdef USE_BLAS
IMPORT_PREFIX void
dgeqrf_(const unsigned *m,
	const unsigned *n,
	double *a,
	const unsigned *lda,
	double *tau,
	double *work,
	const int *lwork,
	int *info);

void
qrdecomp_amatrix(pamatrix a, pavector tau)
{
  uint rows = a->rows;
  uint cols = a->cols;
  uint refl = UINT_MIN(rows, cols);
  double *work;
  int lwork, info;

  assert(a->ld >= rows);
  /* Quick exit if no reflections used */
  if(refl == 0)
    return;

  lwork = 4 * cols;
  work = allocfield(lwork);

  if(tau->dim < refl)
    resize_avector(tau, refl);

  dgeqrf_(&rows, &cols, a->a, &a->ld, tau->v, work, &lwork, &info);
  assert(info == 0);

  freemem(work);
}
#else
void
qrdecomp_amatrix(pamatrix a, pavector tau)
{
  pfield aa = a->a;
  uint lda = a->ld;
  pfield tauv;
  uint rows = a->rows;
  uint cols = a->cols;
  uint refl = UINT_MIN(rows, cols);
  field alpha, beta, gamma, diag;
  real norm2, norm;
  uint i, j, k;

  /* Provide enough storage for scaling factors */
  if(tau->dim < refl)
    resize_avector(tau, refl);
  tauv = tau->v;

  for(k=0; k<refl; k++) {
    /* Compute norm of k-th column */
    norm2 = 0.0;
    for(i=k; i<rows; i++)
      norm2 += ABSSQR(aa[i+k*lda]);
    norm = REAL_SQRT(norm2);

    if(norm2 == 0.0)
      tauv[k] = 0.0;
    else {
      /* Determine reflection vector v */
      diag = aa[k+k*lda];
      alpha = -SIGN(diag) * norm;

      /* Compute norm of v */
      beta = 1.0 / (norm2 - CONJ(alpha) * diag);

      /* Rescale to ensure v_1 = 1 */
      beta *= ABSSQR(diag - alpha);
      gamma = 1.0 / (diag - alpha);
      for(i=k+1; i<rows; i++)
	aa[i+k*lda] *= gamma;
      tauv[k] = beta;

      /* Compute k-th column */
      aa[k+k*lda] = alpha;

      /* Update remaining columns */
      for(j=k+1; j<cols; j++) {
	gamma = aa[k+j*lda];
	for(i=k+1; i<rows; i++)
	  gamma += CONJ(aa[i+k*lda]) * aa[i+j*lda];
	
	gamma *= beta;

	aa[k+j*lda] -= gamma;
	for(i=k+1; i<rows; i++)
	  aa[i+j*lda] -= gamma * aa[i+k*lda];
      }
    }
  }
}
#endif

#ifdef USE_BLAS
IMPORT_PREFIX void
dormqr_(const char *side,
	const char *trans,
	const unsigned *m,
	const unsigned *n,
	const unsigned *k,
	const double *a,
	const unsigned *lda,
	const double *tau,
	double *c,
	const unsigned *ldc,
	double *work,
	const int *lwork,
	int *info);

void
qrsubst_amatrix_avector(bool qtrans, pcamatrix a, pcavector tau, pavector x)
{
  uint rows = a->rows;
  uint cols = a->cols;
  uint refl;
  double work[4];
  int lwork, info;

  refl = UINT_MIN(rows, cols);

  if(refl < 1)
    return;

  assert(tau->dim >= refl);
  assert(x->dim >= rows);

  lwork = 4;

  if(qtrans) {
    dormqr_("Left", "Transposed",
	    &rows, &u_one, &refl,
	    a->a, &a->ld, tau->v,
	    x->v, &x->dim,
	    work, &lwork, &info);
    assert(info == 0);
  }
  else {
    dormqr_("Left", "Not Transposed",
	    &rows, &u_one, &refl,
	    a->a, &a->ld, tau->v,
	    x->v, &x->dim,
	    work, &lwork, &info);
    assert(info == 0);
  }
}

void
qrsubst_amatrix_amatrix(bool qtrans, pcamatrix a, pcavector tau,
			pamatrix x)
{
  uint rows = a->rows;
  uint cols = a->cols;
  uint refl;
  double *work;
  int lwork, info;

  refl = UINT_MIN(rows, cols);

  if(refl < 1 || x->cols < 1)
    return;

  assert(tau->dim >= refl);
  assert(x->rows >= rows);

  lwork = 4 * x->cols;
  work = (double *) allocmem((size_t) sizeof(double) * lwork);

  if(qtrans) {
    dormqr_("Left", "Transposed",
	    &rows, &x->cols, &refl,
	    a->a, &a->ld, tau->v,
	    x->a, &x->ld,
	    work, &lwork, &info);
    assert(info == 0);
  }
  else {
    dormqr_("Left", "Not Transposed",
	    &rows, &x->cols, &refl,
	    a->a, &a->ld, tau->v,
	    x->a, &x->ld,
	    work, &lwork, &info);
    assert(info == 0);
  }

  freemem(work);
}
#else
void
qrsubst_amatrix_avector(bool qtrans, pcamatrix a, pcavector tau, pavector x)
{
  pcfield aa = a->a;
  uint lda = a->ld;
  pcfield tauv = tau->v;
  pfield xv = x->v;
  uint rows = a->rows;
  uint cols = a->cols;
  uint refl;
  field beta, gamma;
  uint i, k;

  refl = UINT_MIN(rows, cols);

  assert(tau->dim >= refl);
  assert(x->dim >= rows);

  if(qtrans) {
    for(k=0; k<refl; k++) {
      beta = tauv[k];
      
      if(beta != 0.0) {
	gamma = xv[k];
	for(i=k+1; i<rows; i++)
	  gamma += CONJ(aa[i+k*lda]) * xv[i];
	
	gamma *= beta;
	
	xv[k] -= gamma;
	for(i=k+1; i<rows; i++)
	  xv[i] -= gamma * aa[i+k*lda];
      }
    }
  }
  else {
    for(k=refl; k-->0; ) {
      beta = tauv[k];
      
      if(beta != 0.0) {
	gamma = xv[k];
	for(i=k+1; i<rows; i++)
	  gamma += CONJ(aa[i+k*lda]) * xv[i];
	
	gamma *= beta;
	
	xv[k] -= gamma;
	for(i=k+1; i<rows; i++)
	  xv[i] -= gamma * aa[i+k*lda];
      }
    }
  }
}

void
qrsubst_amatrix_amatrix(bool qtrans, pcamatrix a, pcavector tau,
			pamatrix x)
{
  pcfield aa = a->a;
  uint lda = a->ld;
  pcfield tauv = tau->v;
  pfield xa = x->a;
  uint ldx = x->ld;
  uint rows = a->rows;
  uint cols = a->cols;
  uint refl;
  field beta, gamma;
  uint i, j, k;

  refl = UINT_MIN(rows, cols);

  assert(tau->dim >= refl);
  assert(x->rows >= rows);

  if(qtrans) {
    for(k=0; k<refl; k++) {
      beta = tauv[k];
      
      if(beta != 0.0)
	for(j=0; j<x->cols; j++) {
	  gamma = xa[k+j*ldx];
	  for(i=k+1; i<rows; i++)
	    gamma += CONJ(aa[i+k*lda]) * xa[i+j*ldx];
	
	  gamma *= beta;
	
	  xa[k+j*ldx] -= gamma;
	  for(i=k+1; i<rows; i++)
	    xa[i+j*ldx] -= gamma * aa[i+k*lda];
	}
    }
  }
  else {
    for(k=refl; k-->0; ) {
      beta = tauv[k];
      
      if(beta != 0.0)
	for(j=0; j<x->cols; j++) {
	  gamma = xa[k+j*ldx];
	  for(i=k+1; i<rows; i++)
	    gamma += CONJ(aa[i+k*lda]) * xa[i+j*ldx];
	  
	  gamma *= beta;
	  
	  xa[k+j*ldx] -= gamma;
	  for(i=k+1; i<rows; i++)
	    xa[i+j*ldx] -= gamma * aa[i+k*lda];
	}
    }
  }
}
#endif

void
qrsolve_amatrix_avector(pcamatrix a, pcavector tau, pavector x)
{
  qrsubst_amatrix_avector(true, a, tau, x);
  uppersolve_amatrix_avector(false, false, a, x);
}

void
qrinvert_amatrix(pamatrix a)
{
  pamatrix acopy;
  pavector tau, b;
  amatrix atmp;
  avector ttmp, btmp;
  uint n = a->rows;
  uint i;

  assert(n == a->cols);

  acopy = init_amatrix(&atmp, n, n);
  tau = init_avector(&ttmp, n);

  copy_amatrix(false, a, acopy);

  qrdecomp_amatrix(acopy, tau);

  identity_amatrix(a);

  for(i=0; i<n; i++) {
    b = init_column_avector(&btmp, a, i);
    qrsolve_amatrix_avector(acopy, tau, b);
    uninit_avector(b);
  }

  uninit_avector(tau);
  uninit_amatrix(acopy);
}

#ifdef USE_BLAS
IMPORT_PREFIX void
dorgqr_(const unsigned *m,
	const unsigned *n,
	const unsigned *k,
	double *a,
	const unsigned *lda,
	double *tau,
	double *work,
	const unsigned *lwork,
	unsigned *info);

void
qrexpand_amatrix(pcamatrix a, pcavector tau, pamatrix q)
{
  double *work;
  uint refl, lwork, info;

  refl = UINT_MIN3(q->cols, a->rows, a->cols);

  /* Quick exit if no reflections used */
  if(refl == 0) {
    identity_amatrix(q);
    return;
  }

  copy_sub_amatrix(false, a, q);

  lwork = 4 * a->rows;
  work = allocfield(lwork);

  dorgqr_(&q->rows, &q->cols, &refl,
	  q->a, &q->ld, tau->v, work, &lwork, &info);
  assert(info == 0);

  freemem(work);
}
#else
void
qrexpand_amatrix(pcamatrix a, pcavector tau, pamatrix q)
{
  pcfield aa = a->a;
  uint lda = a->ld;
  pcfield tauv = tau->v;
  pfield qa = q->a;
  uint ldq = q->ld;
  uint rows = a->rows;
  uint cols = a->cols;
  uint refl;
  field beta, gamma;
  uint i, j, k;

  /* Determine number of relevant elementary reflections */
  refl = UINT_MIN3(q->cols, rows, cols);

  assert(tau->dim >= refl);
  assert(q->rows >= rows);
  assert(q->cols <= rows);

  /* Create identity matrix */
  for(j=0; j<q->cols; j++) {
    for(i=0; i<rows; i++)
      qa[i+j*ldq] = 0.0;
    qa[j+j*ldq] = 1.0;
  }

  /* Apply reflections in reversed order */
  for(k=refl; k-->0; ) {
    beta = tauv[k];
    if(beta != 0.0)
      for(j=k; j<q->cols; j++) {
	gamma = qa[k+j*ldq];
	for(i=k+1; i<rows; i++)
	  gamma += CONJ(aa[i+k*lda]) * qa[i+j*ldq];

	gamma *= beta;

	qa[k+j*ldq] -= gamma;
	for(i=k+1; i<rows; i++)
	  qa[i+j*ldq] -= gamma * aa[i+k*lda];
      }
  }
}
#endif
