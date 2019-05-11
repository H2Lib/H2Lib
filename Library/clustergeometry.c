/* ------------------------------------------------------------
 * This is the file "clustergeometry.c" of the H2Lib package.
 * All rights reserved, Knut Reimer 2009
 * ------------------------------------------------------------ */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "clustergeometry.h"

#include "basic.h"
#include "settings.h"

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

pclustergeometry
new_clustergeometry(uint dim, uint nidx)
{
  pclustergeometry cf;

  uint      i;
  real     *buf;

  cf = (pclustergeometry) allocmem((size_t) sizeof(clustergeometry));
  cf->dim = dim;
  cf->nidx = nidx;
  cf->x = (real **) allocmem((size_t) nidx * sizeof(real *));
  cf->hmin = allocreal(dim);
  cf->hmax = allocreal(dim);
  cf->smin = (real **) allocmem((size_t) nidx * sizeof(real *));
  cf->smax = (real **) allocmem((size_t) nidx * sizeof(real *));
  cf->w = allocreal(nidx);
  cf->buf = buf = allocreal(3 * nidx * dim);

  for (i = 0; i < nidx; i++) {
    cf->x[i] = buf;
    buf += dim;

    cf->smin[i] = buf;
    buf += dim;

    cf->smax[i] = buf;
    buf += dim;

    cf->w[i] = 1.0;
  }

  return cf;
}

void
del_clustergeometry(pclustergeometry cf)
{
  freemem(cf->buf);
  freemem(cf->smax);
  freemem(cf->smin);
  freemem(cf->w);
  freemem(cf->hmax);
  freemem(cf->hmin);
  freemem(cf->x);
  freemem(cf);
}

/* ------------------------------------------------------------
 * Auxiliary routines
 * ------------------------------------------------------------ */

void
update_point_bbox_clustergeometry(pclustergeometry cf, uint size, uint * idx)
{
  uint      i, j;

  assert(size > 0);

  for (j = 0; j < cf->dim; j++) {
    cf->hmin[j] = cf->x[idx[0]][j];
    cf->hmax[j] = cf->x[idx[0]][j];
  }

  for (i = 1; i < size; i++) {
    for (j = 0; j < cf->dim; j++) {
      if (cf->x[idx[i]][j] < cf->hmin[j]) {
	cf->hmin[j] = cf->x[idx[i]][j];
      }
      if (cf->x[idx[i]][j] > cf->hmax[j]) {
	cf->hmax[j] = cf->x[idx[i]][j];
      }
    }
  }
}

void
update_support_bbox_cluster(pclustergeometry cf, pcluster t)
{
  uint      i, j;

  assert(t->size > 0);

  for (j = 0; j < cf->dim; j++) {
    t->bmin[j] = cf->smin[t->idx[0]][j];
    t->bmax[j] = cf->smax[t->idx[0]][j];
  }

  for (i = 1; i < t->size; i++) {
    for (j = 0; j < cf->dim; j++) {
      if (cf->smin[t->idx[i]][j] < t->bmin[j]) {
	t->bmin[j] = cf->smin[t->idx[i]][j];
      }
      if (cf->smax[t->idx[i]][j] > t->bmax[j]) {
	t->bmax[j] = cf->smax[t->idx[i]][j];
      }
    }
  }
}

void
update_bbox_cluster(pcluster t)
{
  uint      i, j;

  assert(t->son > 0);

  for (j = 0; j < t->dim; j++) {
    t->bmin[j] = t->son[0]->bmin[j];
    t->bmax[j] = t->son[0]->bmax[j];
  }

  for (i = 1; i < t->sons; i++) {
    for (j = 0; j < t->dim; j++) {
      t->bmin[j] = REAL_MIN(t->bmin[j], t->son[i]->bmin[j]);
      t->bmax[j] = REAL_MAX(t->bmax[j], t->son[i]->bmax[j]);
    }
  }
}
