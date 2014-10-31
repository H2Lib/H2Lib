/* ------------------------------------------------------------
 This is the file "clustergeometry.c" of the H2Lib package.
 All rights reserved, Knut Reimer 2009
 ------------------------------------------------------------ */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "clustergeometry.h"

#include "basic.h"
#include "settings.h"

/* ------------------------------------------------------------
 Constructors and destructors
 ------------------------------------------------------------ */

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
 Clustering strategies
 ------------------------------------------------------------ */

pcluster
build_adaptive_cluster(pclustergeometry cf, uint size, uint * idx, uint clf)
{
  pcluster  t;

  uint      direction;
  uint      size0, size1;
  uint      i, j;
  real      a, m;

  if (size > clf) {

    update_point_bbox_clustergeometry(cf, size, idx);

    /* compute the direction of partition */
    direction = 0;
    a = cf->hmax[0] - cf->hmin[0];

    for (j = 1; j < cf->dim; j++) {
      m = cf->hmax[j] - cf->hmin[j];
      if (a < m) {
	a = m;
	direction = j;
      }
    }

    /* build sons */
    if (a > 0.0) {
      m = (cf->hmax[direction] + cf->hmin[direction]) / 2.0;
      size0 = 0;
      size1 = 0;

      for (i = 0; i < size; i++) {
	if (cf->x[idx[i]][direction] < m) {
	  j = idx[i];
	  idx[i] = idx[size0];
	  idx[size0] = j;
	  size0++;
	}
	else {
	  size1++;
	}
      }
      t = new_cluster(size, idx, 2, cf->dim);

      t->son[0] = build_adaptive_cluster(cf, size0, idx, clf);
      t->son[1] = build_adaptive_cluster(cf, size1, idx + size0, clf);

      update_bbox_cluster(t);
    }
    else {
      assert(a == 0.0);
      t = new_cluster(size, idx, 0, cf->dim);
      update_support_bbox_cluster(cf, t);
    }
  }
  else {
    t = new_cluster(size, idx, 0, cf->dim);
    update_support_bbox_cluster(cf, t);
  }

  update_cluster(t);

  return t;
}

pcluster
build_regular_cluster(pclustergeometry cf, uint size, uint * idx,
		      uint clf, uint direction)
{
  pcluster  t;

  uint      newd;
  uint      size0, size1;
  uint      i, j;
  real      a, b, m;

  if (size > clf) {
    size0 = 0;
    size1 = 0;

    if (direction < cf->dim - 1) {
      newd = direction + 1;
    }
    else {
      newd = 0;
    }

    m = cf->hmax[direction] - cf->hmin[direction];

    if (m > 0.0) {
      m = (cf->hmax[direction] + cf->hmin[direction]) / 2.0;

      for (i = 0; i < size; i++) {
	if (cf->x[idx[i]][direction] < m) {
	  j = idx[i];
	  idx[i] = idx[size0];
	  idx[size0] = j;
	  size0++;
	}
	else {
	  size1++;
	}
      }

      /* build sons */
      if (size0 > 0) {
	if (size1 > 0) {
	  /* both sons are not empty */
	  t = new_cluster(size, idx, 2, cf->dim);

	  a = cf->hmin[direction];
	  b = cf->hmax[direction];
	  cf->hmax[direction] = m;

	  t->son[0] = build_regular_cluster(cf, size0, idx, clf, newd);

	  cf->hmax[direction] = b;
	  cf->hmin[direction] = m;

	  t->son[1] =
	    build_regular_cluster(cf, size1, idx + size0, clf, newd);

	  cf->hmin[direction] = a;

	  update_bbox_cluster(t);
	}
	else {
	  /* only the first son is not empty */
	  t = new_cluster(size, idx, 1, cf->dim);

	  b = cf->hmax[direction];
	  cf->hmax[direction] = m;

	  t->son[0] = build_regular_cluster(cf, size, idx, clf, newd);

	  cf->hmax[direction] = b;

	  update_bbox_cluster(t);
	}
      }
      else {
	/* only the second son is not empty */ assert(size1 > 0);

	t = new_cluster(size, idx, 1, cf->dim);

	a = cf->hmin[direction];
	cf->hmin[direction] = m;

	t->son[0] = build_regular_cluster(cf, size, idx, clf, newd);

	cf->hmin[direction] = a;

	update_bbox_cluster(t);
      }
    }
    else {
      assert(m == 0.0);
      t = new_cluster(size, idx, 1, cf->dim);

      t->son[0] = build_regular_cluster(cf, size, idx, clf, newd);

      update_bbox_cluster(t);
    }
  }
  else {
    t = new_cluster(size, idx, 0, cf->dim);
    update_support_bbox_cluster(cf, t);
  }

  update_cluster(t);

  return t;
}

/* auxiliary routine for build_simsub_cluster */
static pcluster
build_help_cluster(pclustergeometry cf, uint * idx, uint size,
		   uint clf, uint direction, uint * leaves)
{
  pcluster  s;

  uint      size0, size1;
  uint      i, j;
  real      a, b, m;

  if (direction < cf->dim) {
    size0 = 0;
    size1 = 0;

    m = cf->hmax[direction] - cf->hmin[direction];

    if (m > 0.0) {
      m = (cf->hmax[direction] + cf->hmin[direction]) / 2.0;

      for (i = 0; i < size; i++) {
	if (cf->x[idx[i]][direction] < m) {
	  j = idx[i];
	  idx[i] = idx[size0];
	  idx[size0] = j;
	  size0++;
	}
	else {
	  size1++;
	}
      }

      /* build sons */
      if (size0 > 0) {
	if (size1 > 0) {
	  /* both sons are not empty */
	  s = new_cluster(size, idx, 2, cf->dim);
	  a = cf->hmin[direction];
	  b = cf->hmax[direction];

	  cf->hmax[direction] = m;
	  s->son[0] = build_help_cluster(cf, idx, size0, clf, direction + 1,
					 leaves);
	  cf->hmax[direction] = b;

	  cf->hmin[direction] = m;
	  s->son[1] = build_help_cluster(cf, idx + size0, size1, clf,
					 direction + 1, leaves);
	  cf->hmin[direction] = a;
	}
	else {
	  /* only the first son is not empty */
	  s = new_cluster(size, idx, 1, cf->dim);
	  b = cf->hmax[direction];

	  cf->hmax[direction] = m;
	  s->son[0] = build_help_cluster(cf, idx, size, clf, direction + 1,
					 leaves);
	  cf->hmax[direction] = b;
	}
      }
      else {
	/* only the second son is not empty */ assert(size1 > 0);

	s = new_cluster(size, idx, 1, cf->dim);
	a = cf->hmin[direction];

	cf->hmin[direction] = m;
	s->son[0] = build_help_cluster(cf, idx, size, clf, direction + 1,
				       leaves);
	cf->hmin[direction] = a;
      }
    }
    else {
      assert(m == 0.0);
      s = new_cluster(size, idx, 1, cf->dim);
      s->son[0] =
	build_help_cluster(cf, idx, size, clf, direction + 1, leaves);
    }
  }
  else {
    s = new_cluster(size, idx, 0, cf->dim);
    for (i = 0; i < cf->dim; i++) {
      s->bmin[i] = cf->hmin[i];
      s->bmax[i] = cf->hmax[i];
    }
    leaves[0]++;
  }

  return s;
}

/* auxiliary routine for build_simsub_cluster */
static void
leaves_into_sons(pclustergeometry cf, uint clf, pcluster s,
		 pcluster t, uint * leaves)
{
  uint      i;

  if (s->sons > 0) {
    for (i = 0; i < s->sons; i++) {
      leaves_into_sons(cf, clf, s->son[i], t, leaves);
    }
    free(s->bmin);
    free(s->bmax);
    free(s->son);
  }
  else {
    for (i = 0; i < cf->dim; i++) {
      cf->hmin[i] = s->bmin[i];
      cf->hmax[i] = s->bmax[i];
    }
    t->son[*leaves] = build_simsub_cluster(cf, s->size, s->idx, clf);
    (*leaves)++;
  }
}

pcluster
build_simsub_cluster(pclustergeometry cf, uint size, uint * idx, uint clf)
{
  pcluster  t;

  uint      leaves;
  pcluster  s;

  leaves = 0;

  if (size > clf) {
    s = build_help_cluster(cf, idx, size, clf, 0, &leaves);
    t = new_cluster(size, idx, leaves, cf->dim);

    leaves = 0;
    leaves_into_sons(cf, clf, s, t, &leaves);

    update_bbox_cluster(t);
  }
  else {
    t = new_cluster(size, idx, 0, cf->dim);
    update_support_bbox_cluster(cf, t);
  }

  update_cluster(t);

  return t;
}

pcluster
build_pca_cluster(pclustergeometry cf, uint size, uint * idx, uint clf)
{
  const uint dim = cf->dim;

  pamatrix  C, Q;
  pavector  lambda, v;
  real     *x, *y;
  real      w;
  uint      i, j, k, size0, size1;

  pcluster  t;

  size0 = 0;
  size1 = 0;

  if (size > clf) {
    x = allocreal(dim);
    y = allocreal(dim);

    /* determine weight of current cluster */
    w = 0.0;
    for (i = 0; i < size; ++i) {
      w += cf->w[idx[i]];
    }
    w = 1.0 / w;

    for (j = 0; j < dim; ++j) {
      x[j] = 0.0;
    }

    /* determine center of mass */
    for (i = 0; i < size; ++i) {
      for (j = 0; j < dim; ++j) {
	x[j] += cf->w[idx[i]] * cf->x[idx[i]][j];
      }
    }
    for (j = 0; j < dim; ++j) {
      x[j] *= w;
    }

    C = new_zero_amatrix(dim, dim);
    Q = new_zero_amatrix(dim, dim);
    lambda = new_avector(dim);
    clear_avector(lambda);

    /* setup covariance matrix */
    for (i = 0; i < size; ++i) {

      for (j = 0; j < dim; ++j) {
	y[j] = cf->x[idx[i]][j] - x[j];
      }

      for (j = 0; j < dim; ++j) {
	for (k = 0; k < dim; ++k) {
	  C->a[j + k * C->ld] += cf->w[idx[i]] * y[j] * y[k];
	}
      }
    }

    /* get eigenvalues and eigenvectors of covariance matrix */
    eig_amatrix(C, lambda, Q);

    /* get eigenvector from largest eigenvalue */
    v = new_avector(0);
    init_column_avector(v, Q, dim - 1);

    /* separate cluster with v as separation-plane */
    for (i = 0; i < size; ++i) {
      /* x_i - X */
      for (j = 0; j < dim; ++j) {
	y[j] = cf->x[idx[i]][j] - x[j];
      }

      /* <y,v> */
      w = 0.0;
      for (j = 0; j < dim; ++j) {
	w += y[j] * v->v[j];
      }

      if (w >= 0.0) {
	j = idx[i];
	idx[i] = idx[size0];
	idx[size0] = j;
	size0++;
      }
      else {
	size1++;
      }
    }

    assert(size0 + size1 == size);

    del_amatrix(Q);
    del_amatrix(C);
    del_avector(lambda);
    del_avector(v);
    freemem(x);
    freemem(y);

    /* recursion */
    if (size0 > 0) {
      if (size1 > 0) {
	t = new_cluster(size, idx, 2, cf->dim);

	t->son[0] = build_pca_cluster(cf, size0, idx, clf);
	t->son[1] = build_pca_cluster(cf, size1, idx + size0, clf);

	update_bbox_cluster(t);
      }
      else {
	t = new_cluster(size, idx, 1, cf->dim);
	t->son[0] = build_pca_cluster(cf, size0, idx, clf);

	update_bbox_cluster(t);
      }
    }
    else {
      assert(size1 > 0);
      t = new_cluster(size, idx, 1, cf->dim);
      t->son[0] = build_pca_cluster(cf, size1, idx, clf);

      update_bbox_cluster(t);
    }

  }
  else {
    t = new_cluster(size, idx, 0, cf->dim);
    update_support_bbox_cluster(cf, t);
  }

  update_cluster(t);

  return t;
}

pcluster
build_cluster(pclustergeometry cf, uint size, uint * idx, uint clf,
	      clustermode mode)
{
  pcluster  t;

  if (mode == H2_ADAPTIVE) {
    t = build_adaptive_cluster(cf, size, idx, clf);
  }
  else if (mode == H2_REGULAR) {
    update_point_bbox_clustergeometry(cf, size, idx);
    t = build_regular_cluster(cf, size, idx, clf, 0);
  }
  else if (mode == H2_PCA) {
    t = build_pca_cluster(cf, size, idx, clf);
  }
  else {
    assert(mode == H2_SIMSUB);
    update_point_bbox_clustergeometry(cf, size, idx);
    t = build_simsub_cluster(cf, size, idx, clf);
  }

  return t;
}

/* ------------------------------------------------------------
 Auxiliary routines
 ------------------------------------------------------------ */

void
update_point_bbox_clustergeometry(pclustergeometry cf, uint size, uint * idx)
{
  uint      i, j;

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
