
/* ------------------------------------------------------------
 * This is the file "sparsematrix.c" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2012
 * ------------------------------------------------------------ */

#include "sparsematrix.h"

#include <assert.h>
#include <stdio.h>

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

psparsematrix
new_raw_sparsematrix(uint rows, uint cols, uint nz)
{
  uint      i;

  psparsematrix sp;

  sp = (psparsematrix) allocmem(sizeof(sparsematrix));
  sp->row = allocuint(rows + 1);
  sp->col = allocuint(nz);
  sp->coeff = allocfield(nz);

  sp->rows = rows;
  sp->cols = cols;
  sp->nz = nz;

  for (i = 0; i < rows + 1; ++i) {
    sp->row[i] = 0;
  }

  return sp;
}

psparsematrix
new_identity_sparsematrix(uint rows, uint cols)
{
  uint      i, n;

  psparsematrix sp;

  n = UINT_MIN(rows, cols);

  sp = (psparsematrix) allocmem(sizeof(sparsematrix));
  sp->row = allocuint(rows + 1);
  sp->col = allocuint(n);
  sp->coeff = allocfield(n);

  sp->rows = rows;
  sp->cols = cols;
  sp->nz = n;

  for (i = 0; i < n; ++i) {
    sp->row[i] = i;
    sp->col[i] = i;
    sp->coeff[i] = 1.0;
  }
  for (; i < rows + 1; ++i) {
    sp->row[i] = n;
  }

  return sp;
}

psparsematrix
new_zero_sparsematrix(psparsepattern sp)
{
  psparsematrix A;
  ppatentry e;
  uint     *row;
  uint     *col;
  pfield    coeff;
  uint      rows, cols;
  uint      nz;
  uint      i, j;

  assert(sp != NULL);

  rows = sp->rows;
  cols = sp->cols;

  /* Find out total number of non-zero entries */
  nz = 0;
  for (i = 0; i < sp->rows; i++)
    for (e = sp->row[i]; e != NULL; e = e->next)
      nz++;

  /* Create empty sparse matrix */
  A = new_raw_sparsematrix(rows, cols, nz);
  row = A->row;
  col = A->col;
  coeff = A->coeff;

  /* Initialize pattern */
  j = 0;
  for (i = 0; i < rows; i++) {
    row[i] = j;
    for (e = sp->row[i]; e != NULL; e = e->next) {
      assert(e->row == i);
      col[j] = e->col;
      coeff[j] = 0.0;
      j++;
    }
  }
  assert(j == nz);
  row[i] = j;

  /* Ensure that diagonal entries come first */
  sort_sparsematrix(A);

  return A;
}

void
del_sparsematrix(psparsematrix a)
{
  assert(a != NULL);

  if (a->coeff != NULL)
    freemem(a->coeff);
  if (a->col != NULL)
    freemem(a->col);
  if (a->row != NULL)
    freemem(a->row);
  freemem(a);
}

/* ------------------------------------------------------------
 * Access methods
 * ------------------------------------------------------------ */

field
addentry_sparsematrix(psparsematrix a, uint row, uint col, field x)
{
  uint      i;

  assert(a != NULL);
  assert(row < a->rows);
  assert(col < a->cols);

  for (i = a->row[row]; i < a->row[row + 1] && a->col[i] != col; i++);

  assert(i < a->row[row + 1]);
  assert(a->col[i] == col);

  return a->coeff[i] += x;
}

void
setentry_sparsematrix(psparsematrix a, uint row, uint col, field x)
{
  uint      i;

  assert(a != NULL);
  assert(row < a->rows);
  assert(col < a->cols);

  for (i = a->row[row]; i < a->row[row + 1] && a->col[i] != col; i++);
  assert(i < a->row[row + 1]);
  assert(a->col[i] == col);

  a->coeff[i] = x;
}

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

size_t
getsize_sparsematrix(pcsparsematrix a)
{
  size_t    sz;

  sz = (size_t) sizeof(sparsematrix);
  sz += (size_t) sizeof(uint) * (a->rows + 1);
  sz += (size_t) sizeof(uint) * a->nz;
  sz += (size_t) sizeof(field) * a->nz;

  return sz;
}

/* ------------------------------------------------------------
 * Simple utility functions
 * ------------------------------------------------------------ */

static void
swap(uint i1, uint i2, uint * col, pfield coeff)
{
  uint      hcol;
  field     hcoeff;

  hcol = col[i1];
  col[i1] = col[i2];
  col[i2] = hcol;

  hcoeff = coeff[i1];
  coeff[i1] = coeff[i2];
  coeff[i2] = hcoeff;
}

void
sort_sparsematrix(psparsematrix a)
{
  uint     *row = a->row;
  uint     *col = a->col;
  pfield    coeff = a->coeff;
  uint      rows = a->rows;
  uint      i, j, k;

  for (i = 0; i < rows; i++) {
    k = row[i];
    for (j = row[i]; j < row[i + 1] && i != col[j]; j++);
    if (j < row[i + 1])
      swap(k, j, col, coeff);
  }
}

void
clear_sparsematrix(psparsematrix a)
{
  uint      i;

  for (i = 0; i < a->nz; i++)
    a->coeff[i] = 0.0;
}

void
print_sparsematrix(pcsparsematrix a)
{
  uint     *row;
  uint     *col;
  pfield    coeff;
  uint      rows;
  uint      i, j;

  assert(a != NULL);

  row = a->row;
  col = a->col;
  coeff = a->coeff;
  rows = a->rows;

  (void) printf("sparsematrix(%u,%u,%u)\n", rows, a->cols, a->nz);
  if (a->nz > 0) {
    for (i = 0; i < rows; i++) {
      (void) printf("  %u:", i);

      for (j = row[i]; j < row[i + 1]; j++)
	(void) printf(" (%u " FIELD_CS(.5, e) ")", col[j],
		      FIELD_ARG(coeff[j]));

      (void) printf("\n");
    }
  }
}

void
print_eps_sparsematrix(pcsparsematrix a, const char *filename, uint offset)
{
  const uint *row;
  const uint *col;
  const field *coeff;
  uint      rows, cols, nz;
  FILE     *out;
  real      val, maxval, scale;
  uint      i, j, r;

  assert(a != NULL);

  rows = a->rows;
  cols = a->cols;
  nz = a->nz;
  row = a->row;
  col = a->col;
  coeff = a->coeff;

  maxval = 0.0;
  if (nz > 0) {
    maxval = ABS(coeff[0]);
    for (i = 1; i < nz; i++) {
      val = ABS(coeff[i]);
      if (maxval < val)
	maxval = val;
    }
  }

  scale = 500.0 / UINT_MAX(rows, cols);

  out = fopen(filename, "w");
  assert(out != NULL);

  (void) fprintf(out,
		 "%%!PS-adobe-2.0 EPSF-2.0\n"
		 "%%%%BoundingBox: %d %d %d %d\n"
		 "%d dup translate\n"
		 "%f dup scale\n"
		 "0.1 dup translate\n"
		 "0 setgray\n"
		 "/bx {dup 0 exch 1 exch sub setrgbcolor moveto 0.8 0.0 rlineto 0.0 0.8 rlineto -0.8 0.0 rlineto closepath fill} def\n",
		 offset, offset, offset + (int) (scale * cols + 0.5),
		 offset + (int) (scale * rows + 0.5), offset, scale);

  for (i = 0; i < rows; i++)
    for (r = row[i]; r < row[i + 1]; r++) {
      j = col[r];

      (void) fprintf(out, "%u %u %f bx\n", j, rows - 1 - i,
		     ABS(coeff[r]) / maxval);
    }

  (void) fprintf(out, "showpage\n");
  (void) fclose(out);
}

/* ------------------------------------------------------------
 * Basic linear algebra
 * ------------------------------------------------------------ */

void
addeval_sparsematrix_avector(field alpha, pcsparsematrix a, pcavector x,
			     pavector y)
{
  uint     *row;
  uint     *col;
  field    *coeff;
  uint      rows;
  pcfield   xv;
  pfield    yv;
  register field sum;
  uint      i, j;

  assert(a != NULL);
  assert(x != NULL);
  assert(y != NULL);
  assert(a->rows == y->dim);
  assert(a->cols == x->dim);

  row = a->row;
  col = a->col;
  coeff = a->coeff;
  rows = a->rows;
  xv = x->v;
  yv = y->v;

  for (i = 0; i < rows; i++) {
    sum = 0.0;
    for (j = row[i]; j < row[i + 1]; j++) {
      sum += coeff[j] * xv[col[j]];
    }
    yv[i] += alpha * sum;
  }
}

void
addevaltrans_sparsematrix_avector(field alpha, pcsparsematrix a,
				  pcavector x, pavector y)
{
  uint     *row;
  uint     *col;
  field    *coeff;
  uint      rows;
  pcfield   xv;
  pfield    yv;
  register field val;
  uint      i, j;

  assert(a != NULL);
  assert(x != NULL);
  assert(y != NULL);
  assert(a->rows == x->dim);
  assert(a->cols == y->dim);

  row = a->row;
  col = a->col;
  coeff = a->coeff;
  rows = a->rows;
  xv = x->v;
  yv = y->v;

  for (i = 0; i < rows; i++) {
    val = alpha * xv[i];
    for (j = row[i]; j < row[i + 1]; j++) {
      yv[col[j]] += coeff[j] * val;
    }
  }
}

void
mvm_sparsematrix_avector(field alpha, bool trans, pcsparsematrix a,
			 pcavector x, pavector y)
{
  if (trans) {
    addevaltrans_sparsematrix_avector(alpha, a, x, y);
  }
  else {
    addeval_sparsematrix_avector(alpha, a, x, y);
  }
}

real
norm2_sparsematrix(pcsparsematrix S)
{
  return norm2_matrix((mvm_t) mvm_sparsematrix_avector, (void *) S, S->rows,
		      S->rows);
}

real
norm2diff_sparsematrix(pcsparsematrix a, pcsparsematrix b)
{
  return norm2diff_matrix((mvm_t) mvm_sparsematrix_avector, (void *) a,
			  (mvm_t) mvm_sparsematrix_avector, (void *) b,
			  a->rows, a->rows);
}

void
add_sparsematrix_amatrix(field alpha, bool atrans, pcsparsematrix a,
			 pamatrix b)
{
  const uint *row = a->row;
  const uint *col = a->col;
  pcfield   coeff = a->coeff;
  uint      rows = a->rows;
  uint      cols = a->cols;
  uint      ldb = b->ld;
  uint      i, j, k;

#ifdef HARITH_SPARSEMATRIX_QUICK_EXIT
  if (a->nz == 0)
    return;
#endif
  if (a->nz > 0) {
    if (atrans) {
      assert(b->rows == cols);
      assert(b->cols == rows);

      for (i = 0; i < rows; i++) {
	for (k = row[i]; k < row[i + 1]; k++) {
	  j = col[k];
	  b->a[j + i * ldb] += alpha * coeff[k];
	}
      }
    }
    else {
      assert(b->rows == rows);
      assert(b->cols == cols);

      for (i = 0; i < rows; i++) {
	for (k = row[i]; k < row[i + 1]; k++) {
	  j = col[k];
	  b->a[i + j * ldb] += alpha * coeff[k];
	}
      }
    }
  }
}
