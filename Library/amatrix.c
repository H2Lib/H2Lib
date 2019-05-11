/* ------------------------------------------------------------
 * This is the file "amatrix.c" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

#include "amatrix.h"

#include "settings.h"
#include "basic.h"

#include <math.h>
#include <stdio.h>

static uint active_amatrix = 0;

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

pamatrix
init_amatrix(pamatrix a, uint rows, uint cols)
{
  assert(a != NULL);

  a->a = (rows > 0 && cols > 0 ? allocmatrix(rows, cols) : NULL);
  a->ld = rows;
  a->rows = rows;
  a->cols = cols;
  a->owner = NULL;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_amatrix++;

  return a;
}

pamatrix
init_sub_amatrix(pamatrix a, pamatrix src, uint rows, uint roff,
		 uint cols, uint coff)
{
  longindex ld = src->ld;

  assert(a != NULL);
  assert(src != NULL);
  assert(roff + rows <= src->rows);
  assert(coff + cols <= src->cols);

  a->a = src->a + roff + ld * coff;
  a->ld = src->ld;
  a->rows = rows;
  a->cols = cols;
  a->owner = src;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_amatrix++;

  return a;
}

pamatrix
init_vec_amatrix(pamatrix a, pavector src, uint rows, uint cols)
{
  assert(a != NULL);
  assert(src != NULL);
  assert(rows <= src->dim / cols);

  a->a = src->v;
  a->ld = rows;
  a->rows = rows;
  a->cols = cols;
  a->owner = (pamatrix) src;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_amatrix++;

  return a;
}

pamatrix
init_pointer_amatrix(pamatrix a, pfield src, uint rows, uint cols)
{
  assert(a != NULL);
  assert(rows == 0 || cols == 0 || src != NULL);

  a->a = src;
  a->rows = rows;
  a->ld = rows;
  a->cols = rows;
  a->owner = src;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_amatrix++;

  return a;
}

pamatrix
init_zero_amatrix(pamatrix a, uint rows, uint cols)
{
  longindex lda;
  uint      i, j;

  init_amatrix(a, rows, cols);
  lda = a->ld;

  for (j = 0; j < cols; j++)
    for (i = 0; i < rows; i++)
      a->a[i + j * lda] = 0.0;

  return a;
}

pamatrix
init_identity_amatrix(pamatrix a, uint rows, uint cols)
{
  longindex lda;
  uint      i, j;

  init_amatrix(a, rows, cols);
  lda = a->ld;

  for (j = 0; j < cols; j++)
    for (i = 0; i < rows; i++)
      a->a[i + j * lda] = 0.0;

  for (i = 0; i < rows && i < cols; i++)
    a->a[i + i * lda] = 1.0;

  return a;
}

void
uninit_amatrix(pamatrix a)
{
  assert(a != NULL);

  if (!a->owner && a->a != NULL) {
    freemem(a->a);
  }

  assert(active_amatrix > 0);

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_amatrix--;
}

pamatrix
new_amatrix(uint rows, uint cols)
{
  pamatrix  a;

  a = (pamatrix) allocmem(sizeof(amatrix));

  init_amatrix(a, rows, cols);

  return a;
}

pamatrix
new_sub_amatrix(pamatrix src, uint rows, uint roff, uint cols, uint coff)
{
  pamatrix  a;

  a = (pamatrix) allocmem(sizeof(amatrix));

  init_sub_amatrix(a, src, rows, roff, cols, coff);

  return a;
}

pamatrix
new_pointer_amatrix(field * src, uint rows, uint cols)
{
  pamatrix  a;

  a = (pamatrix) allocmem(sizeof(amatrix));

  init_pointer_amatrix(a, src, rows, cols);

  return a;
}

pamatrix
new_zero_amatrix(uint rows, uint cols)
{
  pamatrix  a;

  a = (pamatrix) allocmem(sizeof(amatrix));

  init_zero_amatrix(a, rows, cols);

  return a;
}

pamatrix
new_identity_amatrix(uint rows, uint cols)
{
  pamatrix  a;

  a = (pamatrix) allocmem(sizeof(amatrix));

  init_identity_amatrix(a, rows, cols);

  return a;
}

void
del_amatrix(pamatrix a)
{
  assert(a != NULL);

  uninit_amatrix(a);
  freemem(a);
}

void
resize_amatrix(pamatrix a, uint rows, uint cols)
{
  assert(a != NULL);
  assert(a->owner == NULL);

  if (rows != a->rows || cols != a->cols) {
    if (a->a != NULL) {
      freemem(a->a);
    }
    a->a = allocmatrix(rows, cols);
    a->rows = rows;
    a->cols = cols;
    a->ld = rows;
  }
}

void
resizecopy_amatrix(pamatrix a, uint rows, uint cols)
{
  pfield    new_a;
  longindex lda, ldn;
  uint      i, j;

  assert(a != NULL);
  assert(a->owner == NULL);

  if (rows != a->rows || cols != a->cols) {
    new_a = allocmatrix(rows, cols);
    lda = a->ld;
    ldn = rows;

    for (j = 0; j < cols && j < a->cols; j++) {
      for (i = 0; i < rows && i < a->rows; i++)
	new_a[i + j * ldn] = a->a[i + j * lda];
      for (; i < rows; i++)
	new_a[i + j * ldn] = 0.0;
    }
    for (; j < cols; j++)
      for (i = 0; i < rows; i++)
	new_a[i + j * ldn] = 0.0;

    if (a->a != NULL) {
      freemem(a->a);
    }

    a->a = new_a;
    a->rows = rows;
    a->cols = cols;
    a->ld = rows;
  }
}

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

uint
getactives_amatrix()
{
  return active_amatrix;
}

size_t
getsize_amatrix(pcamatrix a)
{
  size_t    sz;

  sz = sizeof(amatrix);
  if (a->owner == NULL)
    sz += (size_t) sizeof(field) * a->rows * a->cols;

  return sz;
}

size_t
getsize_heap_amatrix(pcamatrix a)
{
  size_t    sz;

  sz = 0;
  if (a->owner == NULL)
    sz += (size_t) sizeof(field) * a->rows * a->cols;

  return sz;
}

/* ------------------------------------------------------------
 * Simple utility functions
 * ------------------------------------------------------------ */

void
clear_amatrix(pamatrix a)
{
  uint      rows = a->rows;
  longindex lda = a->ld;
  uint      i, j;

  for (j = 0; j < a->cols; j++)
    for (i = 0; i < rows; i++)
      a->a[i + j * lda] = 0.0;
}

void
clear_lower_amatrix(pamatrix a, bool strict)
{
  uint      rows = a->rows;
  uint      cols = a->cols;
  longindex lda = a->ld;
  uint      i, j;

  if (strict) {
    for (j = 0; j < cols; j++) {
      for (i = j + 1; i < rows; i++) {
	a->a[i + j * lda] = 0.0;
      }
    }
  }
  else {
    for (j = 0; j < cols; j++) {
      for (i = j; i < rows; i++) {
	a->a[i + j * lda] = 0.0;
      }
    }
  }
}

void
clear_upper_amatrix(pamatrix a, bool strict)
{
  uint      rows = a->rows;
  uint      cols = a->cols;
  longindex lda = a->ld;
  uint      i, j;

  if (strict) {
    for (j = 0; j < cols; j++) {
      for (i = 0; i < UINT_MIN(j, rows); i++) {
	a->a[i + j * lda] = 0.0;
      }
    }
  }
  else {
    for (j = 0; j < cols; j++) {
      for (i = 0; i <= UINT_MIN(j, rows - 1); i++) {
	a->a[i + j * lda] = 0.0;
      }
    }
  }
}

void
identity_amatrix(pamatrix a)
{
  uint      rows = a->rows;
  uint      cols = a->cols;
  longindex lda = a->ld;
  uint      i, j;

  for (j = 0; j < cols; j++) {
    for (i = 0; i < rows; i++) {
      a->a[i + j * lda] = 0.0;
    }
  }

  for (i = 0; i < cols && i < rows; i++) {
    a->a[i + i * lda] = 1.0;
  }
}

void
random_amatrix(pamatrix a)
{
  uint      rows = a->rows;
  uint      cols = a->cols;
  longindex lda = a->ld;

  uint      i, j;

  for (j = 0; j < cols; j++) {
    for (i = 0; i < rows; i++) {
      a->a[i + j * lda] = FIELD_RAND();
    }
  }
}

void
random_invertible_amatrix(pamatrix a, real alpha)
{
  pfield    aa = a->a;
  uint      rows = a->rows;
  longindex lda = a->ld;
  real      sum;
  uint      i, j;

  assert(rows == a->cols);

  for (j = 0; j < rows; j++) {
    for (i = 0; i < rows; i++) {
      aa[i + j * lda] = FIELD_RAND();
    }
    aa[j + j * lda] = 0.0;
  }

  for (j = 0; j < rows; j++) {
    sum = 0.0;
    for (i = 0; i < rows; i++)
      sum += ABS(aa[i + j * lda]);
    aa[j + j * lda] = sum + alpha;
  }
}

void
random_selfadjoint_amatrix(pamatrix a)
{
  pfield    aa = a->a;
  uint      rows = a->rows;
  longindex lda = a->ld;
  uint      i, j;

  assert(rows == a->cols);

  for (j = 0; j < rows; j++) {
    for (i = 0; i < j; i++) {
      aa[i + j * lda] = FIELD_RAND();
      aa[j + i * lda] = CONJ(aa[i + j * lda]);
    }
    aa[j + j * lda] = REAL_RAND();
  }
}

void
random_spd_amatrix(pamatrix a, real alpha)
{
  pfield    aa = a->a;
  uint      rows = a->rows;
  longindex lda = a->ld;
  real      sum;
  uint      i, j;

  assert(rows == a->cols);

  for (j = 0; j < rows; j++) {
    for (i = 0; i < j; i++) {
      aa[i + j * lda] = FIELD_RAND();
      aa[j + i * lda] = CONJ(aa[i + j * lda]);
    }
    aa[j + j * lda] = 0.0;
  }

  for (j = 0; j < rows; j++) {
    sum = 0.0;
    for (i = 0; i < rows; i++)
      sum += ABS(aa[i + j * lda]);
    aa[j + j * lda] = sum + alpha;
  }
}

void
copy_amatrix(bool atrans, pcamatrix a, pamatrix b)
{
  if (atrans) {
    assert(a->rows == b->cols);
    assert(a->cols == b->rows);

    copy_sub_amatrix(true, a, b);
  }
  else {
    assert(a->rows == b->rows);
    assert(a->cols == b->cols);

    copy_sub_amatrix(false, a, b);
  }
}

void
copy_colpiv_amatrix(bool atrans, pcamatrix a, const uint * colpiv, pamatrix b)
{
  pcfield   aa = a->a;
  longindex lda = a->ld;
  pfield    ba = b->a;
  longindex ldb = b->ld;
  uint      rows = a->rows;
  uint      cols = a->cols;
  uint      i, j, j1;

  if (atrans) {
    assert(b->rows == cols);
    assert(b->cols == rows);

    for (j = 0; j < cols; j++) {
      j1 = colpiv[j];
      assert(j1 < cols);

      for (i = 0; i < rows; i++)
	ba[j + i * ldb] = CONJ(aa[i + j1 * lda]);
    }
  }
  else {
    assert(b->rows == rows);
    assert(b->cols == cols);

    for (j = 0; j < cols; j++) {
      j1 = colpiv[j];
      assert(j1 < cols);

      for (i = 0; i < rows; i++)
	ba[i + j * ldb] = aa[i + j1 * lda];
    }
  }
}

pamatrix
clone_amatrix(pcamatrix src)
{
  pamatrix  trg;

  trg = new_amatrix(src->rows, src->cols);
  copy_amatrix(false, src, trg);

  return trg;
}

void
copy_sub_amatrix(bool atrans, pcamatrix a, pamatrix b)
{
  longindex lda = a->ld;
  longindex ldb = b->ld;
  uint      i, j, rows, cols;

  if (atrans) {
    rows = UINT_MIN(a->cols, b->rows);
    cols = UINT_MIN(a->rows, b->cols);

    for (i = 0; i < rows; i++) {
      for (j = 0; j < cols; j++) {
	b->a[i + j * ldb] = CONJ(a->a[j + i * lda]);
      }
    }
  }
  else {
    rows = UINT_MIN(a->rows, b->rows);
    cols = UINT_MIN(a->cols, b->cols);

    for (j = 0; j < cols; j++) {
      for (i = 0; i < rows; i++) {
	b->a[i + j * ldb] = a->a[i + j * lda];
      }
    }
  }
}

void
print_amatrix(pcamatrix a)
{
  uint      rows = a->rows;
  uint      cols = a->cols;
  longindex lda = a->ld;
  uint      i, j;

  (void) printf("amatrix(%u,%u,%u)\n", rows, cols, a->ld);
  if (rows == 0 || cols == 0)
    return;

  for (i = 0; i < rows; i++) {
    (void) printf("  (" FIELD_CS(+.5, e), FIELD_ARG(a->a[i]));
    for (j = 1; j < cols; j++)
      (void) printf(" | " FIELD_CS(+.5, e), FIELD_ARG(a->a[i + j * lda]));
    (void) printf(")\n");
  }
}

void
print_matlab_amatrix(pcamatrix a)
{
  uint      rows = a->rows;
  uint      cols = a->cols;
  longindex lda = a->ld;
  uint      i, j;

  if (rows == 0)
    (void) printf("  [ ]\n");
  else if (rows == 1) {
    if (cols == 0)
      (void) printf("  [ ]\n");
    else {
      (void) printf("  [" FIELD_CS(.5, e), FIELD_ARG(a->a[0]));
      for (j = 1; j < cols; j++)
	(void) printf(" " FIELD_CS(.5, e), FIELD_ARG(a->a[j * lda]));
      (void) printf("]\n");
    }
  }
  else {
    if (cols == 0) {
      (void) printf("  [");
      for (i = 1; i < rows; i++)
	(void) printf(" ;");
      (void) printf(" ]\n");
    }
    else {
      (void) printf("  [" FIELD_CS(.5, e), FIELD_ARG(a->a[0]));
      for (j = 1; j < cols; j++)
	(void) printf(" " FIELD_CS(.5, e), FIELD_ARG(a->a[j * lda]));
      (void) printf(" ;\n");
      for (i = 1; i < rows - 1; i++) {
	(void) printf("  " FIELD_CS(.5, e), FIELD_ARG(a->a[i]));
	for (j = 1; j < cols; j++)
	  (void) printf(" " FIELD_CS(.5, e), FIELD_ARG(a->a[i + j * lda]));
	(void) printf(" ;\n");
      }
      (void) printf("  " FIELD_CS(.5, e), FIELD_ARG(a->a[i]));
      for (j = 1; j < cols; j++)
	(void) printf(" " FIELD_CS(.5, e), FIELD_ARG(a->a[i + j * lda]));
      (void) printf("]\n");
    }
  }
}

real
check_ortho_amatrix(bool atrans, pcamatrix a)
{
  pamatrix  a2;
  uint      rows = a->rows;
  uint      cols = a->cols;
  real      norm, error;

  error = 0.0;

  if (atrans) {
    a2 = new_identity_amatrix(rows, rows);
    addmul_amatrix(-1.0, false, a, true, a, a2);
    norm = normfrob_amatrix(a2);
    del_amatrix(a2);

    error = REAL_MAX(error, norm);
  }
  else {
    a2 = new_identity_amatrix(cols, cols);
    addmul_amatrix(-1.0, true, a, false, a, a2);
    norm = normfrob_amatrix(a2);
    del_amatrix(a2);

    error = REAL_MAX(error, norm);
  }

  return error;
}

/* ------------------------------------------------------------
 * Basic linear algebra
 * ------------------------------------------------------------ */

#ifdef USE_BLAS
void
scale_amatrix(field alpha, pamatrix a)
{
  uint      rows = a->rows;
  uint      cols = a->cols;
  longindex lda = a->ld;
  uint      j;

  for (j = 0; j < cols; j++) {
    h2_scal(&rows, &alpha, a->a + j * lda, &u_one);
  }

}
#else
void
scale_amatrix(field alpha, pamatrix a)
{
  uint      rows = a->rows;
  uint      cols = a->cols;
  longindex lda = a->ld;
  uint      i, j;

  for (j = 0; j < cols; j++) {
    for (i = 0; i < rows; i++) {
      a->a[i + j * lda] *= alpha;
    }
  }
}
#endif

void
conjugate_amatrix(pamatrix a)
{
#ifdef USE_COMPLEX
  const uint rows = a->rows;
  const uint cols = a->cols;
  const longindex ld = a->ld;
  field    *aa = a->a;

  uint      i, j;

  for (j = 0; j < cols; ++j) {
    for (i = 0; i < rows; ++i) {
      aa[i + j * ld] = CONJ(aa[i + j * ld]);
    }
  }
#endif
}

#ifdef USE_BLAS
field
dotprod_amatrix(pcamatrix a, pcamatrix b)
{
  field     sum;
  uint      rows, cols;
  longindex lda = a->ld;
  longindex ldb = b->ld;
  uint      j;

  rows = a->rows;
  cols = a->cols;

  assert(rows == b->rows);
  assert(cols == b->cols);

  sum = 0.0;
  for (j = 0; j < cols; j++) {
    sum += h2_dot(&rows, a->a + j * lda, &u_one, b->a + j * ldb, &u_one);
  }

  return sum;
}
#else
field
dotprod_amatrix(pcamatrix a, pcamatrix b)
{
  field     sum;
  uint      rows, cols;
  longindex lda = a->ld;
  longindex ldb = b->ld;
  uint      i, j;

  rows = a->rows;
  cols = a->cols;

  assert(rows == b->rows);
  assert(cols == b->cols);

  sum = 0.0;
  for (j = 0; j < cols; j++) {
    for (i = 0; i < rows; i++) {
      sum += CONJ(a->a[i + j * lda]) * b->a[i + j * ldb];
    }
  }

  return sum;
}
#endif

real
norm2_amatrix(pcamatrix A)
{
  return norm2_matrix((mvm_t) mvm_amatrix_avector, (void *) A, A->rows,
		      A->cols);
}

#ifdef USE_BLAS
real
normfrob_amatrix(pcamatrix a)
{
  real      sum;
  longindex lda = a->ld;
  uint      j;

  sum = 0.0;
  for (j = 0; j < a->cols; j++) {
    sum += REAL_SQR(h2_nrm2(&a->rows, a->a + j * lda, &u_one));
  }

  return REAL_SQRT(sum);
}
#else
real
normfrob_amatrix(pcamatrix a)
{
  real      sum;
  longindex lda = a->ld;
  uint      i, j;

  sum = 0.0;
  for (j = 0; j < a->cols; j++) {
    for (i = 0; i < a->rows; i++) {
      sum += ABSSQR(a->a[i + j * lda]);
    }
  }

  return REAL_SQRT(sum);
}
#endif

#ifdef USE_BLAS
real
normfrob2_amatrix(pcamatrix a)
{
  real      sum;
  longindex lda = a->ld;
  uint      j;

  sum = 0.0;
  for (j = 0; j < a->cols; j++) {
    sum += h2_dot(&a->rows, a->a + j * lda, &u_one, a->a + j * lda, &u_one);
  }

  return sum;
}
#else
real
normfrob2_amatrix(pcamatrix a)
{
  real      sum;
  longindex lda = a->ld;
  uint      i, j;

  sum = 0.0;
  for (j = 0; j < a->cols; j++) {
    for (i = 0; i < a->rows; i++) {
      sum += ABSSQR(a->a[i + j * lda]);
    }
  }

  return sum;
}
#endif

real
norm2diff_amatrix(pcamatrix a, pcamatrix b)
{
  return norm2diff_matrix((mvm_t) mvm_amatrix_avector, (void *) a,
			  (mvm_t) mvm_amatrix_avector, (void *) b, a->rows,
			  a->cols);
}

#ifdef USE_BLAS
void
addeval_amatrix_avector(field alpha, pcamatrix a, pcavector src, pavector trg)
{
  assert(src->dim >= a->cols);
  assert(trg->dim >= a->rows);

  if (a->rows > 0 && a->cols > 0) {
    h2_gemv(_h2_ntrans, &a->rows, &a->cols, &alpha, a->a, &a->ld, src->v,
	    &u_one, &f_one, trg->v, &u_one);
  }
}

void
addevaltrans_amatrix_avector(field alpha, pcamatrix a, pcavector src,
			     pavector trg)
{
  assert(src->dim >= a->rows);
  assert(trg->dim >= a->cols);

  if (a->rows > 0 && a->cols > 0) {
    h2_gemv(_h2_adj, &a->rows, &a->cols, &alpha, a->a, &a->ld, src->v, &u_one,
	    &f_one, trg->v, &u_one);
  }
}
#else
void
addeval_amatrix_avector(field alpha, pcamatrix a, pcavector src, pavector trg)
{
  field     sum;
  longindex lda = a->ld;
  uint      i, j;

  assert(src->dim >= a->cols);
  assert(trg->dim >= a->rows);

  for (i = 0; i < a->rows; i++) {
    sum = f_zero;
    for (j = 0; j < a->cols; j++) {
      sum += a->a[i + j * lda] * src->v[j];
    }
    trg->v[i] += alpha * sum;
  }
}

void
addevaltrans_amatrix_avector(field alpha, pcamatrix a, pcavector src,
			     pavector trg)
{
  field     sum;
  longindex lda = a->ld;
  uint      i, j;

  assert(src->dim >= a->rows);
  assert(trg->dim >= a->cols);

  for (j = 0; j < a->cols; j++) {
    sum = f_zero;
    for (i = 0; i < a->rows; i++) {
      sum += CONJ(a->a[i + j * lda]) * src->v[i];
    }
    trg->v[j] += alpha * sum;
  }
}
#endif

void
mvm_amatrix_avector(field alpha, bool atrans, pcamatrix a, pcavector src,
		    pavector trg)
{
  if (atrans)
    addevaltrans_amatrix_avector(alpha, a, src, trg);
  else
    addeval_amatrix_avector(alpha, a, src, trg);
}

#ifdef USE_BLAS
void
add_amatrix(field alpha, bool atrans, pcamatrix a, pamatrix b)
{
  longindex lda = a->ld;
  longindex ldb = b->ld;
  uint      rows = a->rows;
  uint      cols = a->cols;
  uint      i;

  if (atrans) {
    assert(rows <= b->cols);
    assert(cols <= b->rows);

    for (i = 0; i < cols; ++i) {
      h2_gerc(&u_one, &rows, &alpha, &f_one, &u_one, a->a + i * lda, &u_one,
	      b->a + i, &b->ld);
    }
  }
  else {
    assert(rows <= b->rows);
    assert(cols <= b->cols);

    for (i = 0; i < cols; i++) {
      h2_axpy(&rows, &alpha, a->a + i * lda, &u_one, b->a + i * ldb, &u_one);
    }
  }
}
#else
void
add_amatrix(field alpha, bool atrans, pcamatrix a, pamatrix b)
{
  longindex lda = a->ld;
  longindex ldb = b->ld;
  uint      rows = a->rows;
  uint      cols = a->cols;
  uint      i, j;

  if (atrans) {
    assert(rows <= b->cols);
    assert(cols <= b->rows);

    for (i = 0; i < rows; i++) {
      for (j = 0; j < cols; j++) {
	b->a[j + i * ldb] += alpha * CONJ(a->a[i + j * lda]);
      }
    }
  }
  else {
    assert(rows <= b->rows);
    assert(cols <= b->cols);

    for (j = 0; j < cols; j++) {
      for (i = 0; i < rows; i++) {
	b->a[i + j * ldb] += alpha * a->a[i + j * lda];
      }
    }
  }
}
#endif

#ifdef USE_BLAS
void
addmul_amatrix(field alpha, bool atrans, pcamatrix a, bool btrans,
	       pcamatrix b, pamatrix c)
{

  if (a->rows == 0 || a->cols == 0 || b->rows == 0 || b->cols == 0) {
    return;
  }

  if (atrans) {
    if (btrans) {
      assert(a->cols <= c->rows);
      assert(b->rows <= c->cols);
      assert(a->rows == b->cols);

      h2_gemm(_h2_adj, _h2_adj, &a->cols, &b->rows, &a->rows, &alpha, a->a,
	      &a->ld, b->a, &b->ld, &f_one, c->a, &c->ld);
    }
    else {
      assert(a->cols <= c->rows);
      assert(b->cols <= c->cols);
      assert(a->rows == b->rows);

      h2_gemm(_h2_adj, _h2_ntrans, &a->cols, &b->cols, &a->rows, &alpha, a->a,
	      &a->ld, b->a, &b->ld, &f_one, c->a, &c->ld);
    }
  }
  else {
    if (btrans) {
      assert(a->rows <= c->rows);
      assert(b->rows <= c->cols);
      assert(a->cols == b->cols);

      h2_gemm(_h2_ntrans, _h2_adj, &a->rows, &b->rows, &a->cols, &alpha, a->a,
	      &a->ld, b->a, &b->ld, &f_one, c->a, &c->ld);
    }
    else {
      assert(a->rows <= c->rows);
      assert(b->cols <= c->cols);
      assert(a->cols == b->rows);

      h2_gemm(_h2_ntrans, _h2_ntrans, &a->rows, &b->cols, &a->cols, &alpha,
	      a->a, &a->ld, b->a, &b->ld, &f_one, c->a, &c->ld);
    }
  }
}
#else
void
addmul_amatrix(field alpha, bool atrans, pcamatrix a, bool btrans,
	       pcamatrix b, pamatrix c)
{
  uint      rows, cols, mid;
  pcfield   aa = a->a;
  pcfield   ba = b->a;
  pfield    ca = c->a;
  longindex lda = a->ld;
  longindex ldb = b->ld;
  longindex ldc = c->ld;
  register uint i, j, k;

  if (a->rows == 0 || a->cols == 0 || b->rows == 0 || b->cols == 0) {
    return;
  }

  if (atrans) {
    if (btrans) {
      assert(a->cols <= c->rows);
      assert(b->rows <= c->cols);
      assert(a->rows == b->cols);

      rows = a->cols;
      mid = a->rows;
      cols = b->rows;

      for (i = 0; i < rows; i++) {
	for (k = 0; k < cols; k++) {
	  for (j = 0; j < mid; j++) {
	    ca[i + k * ldc] += alpha
	      * CONJ(aa[j + i * lda]) * CONJ(ba[k + j * ldb]);
	  }
	}
      }
    }
    else {
      assert(a->cols <= c->rows);
      assert(b->cols <= c->cols);
      assert(a->rows == b->rows);

      rows = a->cols;
      mid = a->rows;
      cols = b->cols;

      for (k = 0; k < cols; k++) {
	for (i = 0; i < rows; i++) {
	  for (j = 0; j < mid; j++) {
	    ca[i + k * ldc] +=
	      alpha * CONJ(aa[j + i * lda]) * ba[j + k * ldb];
	  }
	}
      }
    }
  }
  else {
    if (btrans) {
      assert(a->rows <= c->rows);
      assert(b->rows <= c->cols);
      assert(a->cols == b->cols);

      rows = a->rows;
      mid = a->cols;
      cols = b->rows;

      for (j = 0; j < mid; j++) {
	for (k = 0; k < cols; k++) {
	  for (i = 0; i < rows; i++) {
	    ca[i + k * ldc] +=
	      alpha * aa[i + j * lda] * CONJ(ba[k + j * ldb]);
	  }
	}
      }
    }
    else {
      assert(a->rows <= c->rows);
      assert(b->cols <= c->cols);
      assert(a->cols == b->rows);

      rows = a->rows;
      mid = a->cols;
      cols = b->cols;

      for (k = 0; k < cols; k++) {
	for (j = 0; j < mid; j++) {
	  for (i = 0; i < rows; i++) {
	    ca[i + k * ldc] += alpha * aa[i + j * lda] * ba[j + k * ldb];
	  }
	}
      }
    }
  }
}
#endif

#ifdef USE_BLAS
void
bidiagmul_amatrix(field alpha, bool atrans, pamatrix a, pcavector d,
		  pcavector l)
{
  field     beta, gamma;
  longindex lda = a->ld;
  uint      j;

  if (atrans) {
    assert(a->rows == d->dim);
    assert(l->dim + 1 == d->dim);

    if (a->rows < 1 || a->cols < 1) {
      return;
    }

    for (j = 0; j + 1 < a->rows; j++) {
      gamma = CONJ(alpha * d->v[j]);
      h2_scal(&a->cols, &gamma, a->a + j, &a->ld);
      beta = CONJ(alpha * l->v[j]);
      h2_axpy(&a->cols, &beta, a->a + (j + 1), &a->ld, a->a + j, &a->ld);
    }
    gamma = CONJ(alpha * d->v[j]);
    h2_scal(&a->cols, &gamma, a->a + j, &a->ld);
  }
  else {
    assert(a->cols == d->dim);
    assert(l->dim + 1 == d->dim);

    if (a->rows < 1 || a->cols < 1) {
      return;
    }

    for (j = 0; j + 1 < a->cols; j++) {
      gamma = alpha * d->v[j];
      h2_scal(&a->rows, &gamma, a->a + j * lda, &u_one);
      beta = alpha * l->v[j];
      h2_axpy(&a->rows, &beta, a->a + (j + 1) * lda, &u_one, a->a + j * lda,
	      &u_one);
    }
    gamma = alpha * d->v[j];
    h2_scal(&a->rows, &gamma, a->a + j * lda, &u_one);
  }
}
#else
void
bidiagmul_amatrix(field alpha, bool atrans, pamatrix a, pcavector d,
		  pcavector l)
{
  field     beta, gamma;
  longindex lda = a->ld;
  uint      i, j;

  if (atrans) {
    assert(a->rows == d->dim);
    assert(l->dim + 1 == d->dim);

    if (a->rows < 1 || a->cols < 1) {
      return;
    }

    for (j = 0; j + 1 < a->rows; j++) {
      gamma = CONJ(alpha * d->v[j]);
      beta = CONJ(alpha * l->v[j]);
      for (i = 0; i < a->cols; i++)
	a->a[j + i * lda] = gamma * a->a[j + i * lda]
	  + beta * a->a[(j + 1) + i * lda];
    }
    gamma = CONJ(alpha * d->v[j]);
    for (i = 0; i < a->cols; i++)
      a->a[j + i * lda] *= gamma;
  }
  else {
    assert(a->cols == d->dim);
    assert(l->dim + 1 == d->dim);

    if (a->rows < 1 || a->cols < 1) {
      return;
    }

    for (j = 0; j + 1 < a->cols; j++) {
      gamma = alpha * d->v[j];
      beta = alpha * l->v[j];
      for (i = 0; i < a->rows; i++)
	a->a[i + j * lda] = gamma * a->a[i + j * lda]
	  + beta * a->a[i + (j + 1) * lda];
    }
    gamma = alpha * d->v[j];
    for (i = 0; i < a->rows; i++)
      a->a[i + j * lda] *= gamma;
  }
}
#endif
