/* ------------------------------------------------------------
 * This is the file "avector.c" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

#include "avector.h"

#include "settings.h"
#include "basic.h"

#include <math.h>
#include <stdio.h>

static uint active_avector = 0;

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

pavector
init_avector(pavector v, uint dim)
{
  assert(v != NULL);

  v->v = (dim > 0 ? allocfield(dim) : NULL);
  v->dim = dim;
  v->owner = NULL;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_avector++;

  return v;
}

pavector
init_sub_avector(pavector v, pavector src, uint dim, uint off)
{
  assert(v != NULL);
  assert(src != NULL);
  assert(off + dim <= src->dim);

  v->v = src->v + off;
  v->dim = dim;
  v->owner = src;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_avector++;

  return v;
}

pavector
init_zero_avector(pavector v, uint dim)
{
  assert(v != NULL);

  init_avector(v, dim);
  clear_avector(v);

  return v;
}

pavector
init_column_avector(pavector v, pamatrix src, uint col)
{
  longindex lda;

  assert(v != NULL);
  assert(src != NULL);
  assert(col < src->cols);

  lda = src->ld;

  v->v = src->a + col * lda;
  v->dim = src->rows;
  v->owner = src;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_avector++;

  return v;
}

pavector
init_pointer_avector(pavector v, pfield src, uint dim)
{
  assert(v != NULL);
  assert(dim == 0 || src != NULL);

  v->v = src;
  v->dim = dim;
  v->owner = src;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_avector++;

  return v;
}

void
uninit_avector(pavector v)
{
  if (!v->owner)
    freemem(v->v);

  assert(active_avector > 0);

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_avector--;
}

pavector
new_avector(uint dim)
{
  pavector  v;

  v = (pavector) allocmem(sizeof(avector));

  init_avector(v, dim);

  return v;
}

pavector
new_sub_avector(pavector src, uint dim, uint off)
{
  pavector  v;

  v = (pavector) allocmem(sizeof(avector));

  init_sub_avector(v, src, dim, off);

  return v;
}

pavector
new_zero_avector(uint dim)
{
  pavector  v;

  v = (pavector) allocmem(sizeof(amatrix));

  init_zero_avector(v, dim);

  return v;
}

pavector
new_pointer_avector(pfield src, uint dim)
{
  pavector  v;

  v = (pavector) allocmem(sizeof(avector));

  init_pointer_avector(v, src, dim);

  return v;
}

void
del_avector(pavector v)
{
  uninit_avector(v);
  freemem(v);
}

void
resize_avector(pavector v, uint dim)
{
  assert(v->owner == NULL);

  if (dim != v->dim) {
    freemem(v->v);
    v->v = allocfield(dim);
    v->dim = dim;
  }
}

void
shrink_avector(pavector v, uint dim)
{
  assert(dim <= v->dim);

  v->dim = dim;
}

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

uint
getactives_avector()
{
  return active_avector;
}

size_t
getsize_avector(pcavector v)
{
  size_t    sz;

  sz = sizeof(avector);
  if (v->owner == NULL)
    sz += (size_t) sizeof(field) * v->dim;

  return sz;
}

size_t
getsize_heap_avector(pcavector v)
{
  size_t    sz;

  sz = 0;
  if (v->owner == NULL)
    sz += (size_t) sizeof(field) * v->dim;

  return sz;
}

/* ------------------------------------------------------------
 * Simple utility functions
 * ------------------------------------------------------------ */

void
clear_avector(pavector v)
{
  uint      i;

  for (i = 0; i < v->dim; i++)
    v->v[i] = 0.0;
}

void
fill_avector(pavector v, field x)
{
  uint      i;

  for (i = 0; i < v->dim; i++)
    v->v[i] = x;
}

void
random_avector(pavector v)
{
  uint      i;

  for (i = 0; i < v->dim; i++) {
    v->v[i] = FIELD_RAND();
  }
}

void
random_real_avector(pavector v)
{
  uint      i;

  for (i = 0; i < v->dim; i++) {
#ifdef USE_COMPLEX
    v->v[i] = REAL_RAND() + 0.0 * I;
#else
    v->v[i] = REAL_RAND();
#endif
  }
}

void
copy_avector(pcavector v, pavector w)
{
  uint      i;

  assert(v->dim == w->dim);

  for (i = 0; i < v->dim; i++)
    w->v[i] = v->v[i];
}

void
copy_sub_avector(pcavector v, pavector w)
{
  uint      i, n;

  n = UINT_MIN(v->dim, w->dim);

  for (i = 0; i < n; i++)
    w->v[i] = v->v[i];
}

void
print_avector(pcavector v)
{
  uint      dim = v->dim;
  uint      i;

  (void) printf("avector(%u)\n", dim);
  if (dim == 0)
    return;

  (void) printf("  (" FIELD_CS(.5, e), FIELD_ARG(v->v[0]));
  for (i = 1; i < dim; i++)
    (void) printf(" " FIELD_CS(.5, e), FIELD_ARG(v->v[i]));
  (void) printf(")\n");
}

/* ------------------------------------------------------------
 * Very basic linear algebra
 * ------------------------------------------------------------ */

#ifdef USE_BLAS
void
scale_avector(field alpha, pavector v)
{
  h2_scal(&v->dim, &alpha, v->v, &u_one);
}
#else
void
scale_avector(field alpha, pavector v)
{
  uint      i;

  for (i = 0; i < v->dim; i++)
    v->v[i] *= alpha;
}
#endif

#ifdef USE_BLAS
real
norm2_avector(pcavector v)
{
  return h2_nrm2(&v->dim, v->v, &u_one);
}
#else
real
norm2_avector(pcavector v)
{
  real      sum;
  uint      i;

  sum = 0.0;
  for (i = 0; i < v->dim; i++)
    sum += ABSSQR(v->v[i]);

  return REAL_SQRT(sum);
}
#endif

#ifdef USE_BLAS
field
dotprod_avector(pcavector x, pcavector y)
{
  assert(x->dim == y->dim);

  return h2_dot(&x->dim, x->v, &u_one, y->v, &u_one);
}
#else
field
dotprod_avector(pcavector x, pcavector y)
{
  register field alpha;
  uint      i;

  assert(x->dim == y->dim);

  alpha = 0.0;
  for (i = 0; i < x->dim; i++)
    alpha += CONJ(x->v[i]) * y->v[i];

  return alpha;
}
#endif

#ifdef USE_BLAS
void
add_avector(field alpha, pcavector x, pavector y)
{
  assert(y->dim >= x->dim);

  h2_axpy(&x->dim, &alpha, x->v, &u_one, y->v, &u_one);
}
#else
void
add_avector(field alpha, pcavector x, pavector y)
{
  uint      i;

  assert(y->dim >= x->dim);

  for (i = 0; i < x->dim; i++) {
    y->v[i] += alpha * x->v[i];
  }
}
#endif
