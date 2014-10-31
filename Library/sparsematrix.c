
/* ------------------------------------------------------------
   This is the file "sparsematrix.c" of the H2Lib package.
   All rights reserved, Steffen Boerm 2012
   ------------------------------------------------------------ */

#include "sparsematrix.h"

#include <assert.h>
#include <stdio.h>

/* ------------------------------------------------------------
   Constructors and destructors
   ------------------------------------------------------------ */

psparsematrix
new_raw_sparsematrix(uint rows, uint cols, uint nz)
{
  psparsematrix sp;

  sp = (psparsematrix) allocmem(sizeof(sparsematrix));
  sp->row = allocuint(rows + 1);
  sp->col = allocuint(nz);
  sp->coeff = allocfield(nz);

  sp->rows = rows;
  sp->cols = cols;
  sp->nz = nz;

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

  freemem(a->coeff);
  freemem(a->col);
  freemem(a->row);
  freemem(a);
}

/* ------------------------------------------------------------
   Access methods
   ------------------------------------------------------------ */

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
   Statistics
   ------------------------------------------------------------ */

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
   Simple utility functions
   ------------------------------------------------------------ */

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

  for (i = 0; i < rows; i++) {
    (void) printf("  %u:", i);

    for (j = row[i]; j < row[i + 1]; j++)
      (void) printf(" (%u %g)", col[j], coeff[j]);

    (void) printf("\n");
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
    maxval = fabs(coeff[0]);
    for (i = 1; i < nz; i++) {
      val = fabs(coeff[i]);
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
		 offset, offset,
		 offset + (int) (scale * cols + 0.5),
		 offset + (int) (scale * rows + 0.5), offset, scale);

  for (i = 0; i < rows; i++)
    for (r = row[i]; r < row[i + 1]; r++) {
      j = col[r];

      (void) fprintf(out,
		     "%u %u %f bx\n",
		     j, rows - 1 - i, fabs(coeff[r]) / maxval);
    }

  (void) fprintf(out, "showpage\n");
  (void) fclose(out);
}

/* ------------------------------------------------------------
   Basic linear algebra
   ------------------------------------------------------------ */

void
addeval_sparsematrix_avector(field alpha, pcsparsematrix a,
			     pcavector x, pavector y)
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
    for (j = row[i]; j < row[i + 1]; j++)
      sum += coeff[j] * xv[col[j]];
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
    for (j = row[i]; j < row[i + 1]; j++)
      yv[col[j]] += coeff[j] * val;
  }
}

void
mvm_sparsematrix_avector(field alpha, bool trans, pcsparsematrix a,
			 pcavector x, pavector y)
{
  if (trans)
    addevaltrans_sparsematrix_avector(alpha, a, x, y);
  else
    addeval_sparsematrix_avector(alpha, a, x, y);
}

real
norm2_sparsematrix(pcsparsematrix a)
{
  avector   tmp1, tmp2;
  pavector  x, y;
  real      norm;
  uint      i;

  x = init_avector(&tmp1, a->cols);
  y = init_avector(&tmp2, a->rows);

  random_avector(x);
  norm = norm2_avector(x);
  i = 0;
  while (i < NORM_STEPS && norm > 0.0) {
    scale_avector(1.0 / norm, x);

    clear_avector(y);
    addeval_sparsematrix_avector(1.0, a, x, y);

    clear_avector(x);
    addevaltrans_sparsematrix_avector(1.0, a, y, x);

    norm = norm2_avector(x);
    i++;
  }

  uninit_avector(y);
  uninit_avector(x);

  return REAL_SQRT(norm);
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

  if (atrans) {
    assert(b->rows == cols);
    assert(b->cols == rows);

    for (i = 0; i < rows; i++)
      for (k = row[i]; k < row[i + 1]; k++) {
	j = col[k];
	b->a[j + i * ldb] += alpha * coeff[k];
      }
  }
  else {
    assert(b->rows == rows);
    assert(b->cols == cols);

    for (i = 0; i < rows; i++)
      for (k = row[i]; k < row[i + 1]; k++) {
	j = col[k];
	b->a[i + k * ldb] += alpha * coeff[k];
      }
  }
}
