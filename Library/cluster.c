/* ------------------------------------------------------------
 * This is the file "cluster.c" of the H2Lib package.
 * All rights reserved, Knut Reimer 2009
 * ------------------------------------------------------------ */

#include <assert.h>
#include <math.h>

#include "basic.h"
#include "cluster.h"

#ifdef USE_NETCDF
#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#endif

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

pcluster
new_cluster(uint size, uint * idx, uint sons, uint dim)
{
  pcluster  t;

  uint      i;

  t = (pcluster) allocmem((size_t) sizeof(cluster));
  t->type = 0;
  t->size = size;
  t->dim = dim;
  t->idx = idx;
  t->bmin = allocreal(dim);
  t->bmax = allocreal(dim);
  t->sons = sons;
  if (sons > 0) {
    t->son = (pcluster *) allocmem((size_t) sons * sizeof(pcluster));
    for (i = 0; i < sons; i++)
      t->son[i] = 0;
  }
  else
    t->son = 0;

  return t;
}

void
del_cluster(pcluster t)
{
  uint      i;

  if (t->sons > 0) {
    for (i = 0; i < t->sons; i++)
      del_cluster(t->son[i]);
    freemem(t->son);
  }
  freemem(t->bmax);
  freemem(t->bmin);
  freemem(t);
}

void
update_cluster(pcluster t)
{
  uint      i, desc;

  desc = 1;
  for (i = 0; i < t->sons; i++)
    desc += t->son[i]->desc;

  t->desc = desc;
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

  assert(size > 0);

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

      /* build sons */
      if (size0 > 0) {
	if (size1 > 0) {
	  /* both sons are not empty */
	  t = new_cluster(size, idx, 2, cf->dim);

	  t->son[0] = build_adaptive_cluster(cf, size0, idx, clf);
	  t->son[1] = build_adaptive_cluster(cf, size1, idx + size0, clf);

	  update_bbox_cluster(t);
	}
	else {
	  /* only the first son is not empty */
	  assert(size0 == size);

	  t = new_cluster(size, idx, 0, cf->dim);
	  update_bbox_cluster(t);
	}
      }
      else {
	/* only the second son is not empty */
	assert(size1 > 0);
	assert(size1 == size);

	t = new_cluster(size, idx, 0, cf->dim);
	update_bbox_cluster(t);
      }
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

  assert(size > 0);

  if (size > clf) {
    size0 = 0;
    size1 = 0;

    update_point_bbox_clustergeometry(cf, size, idx);

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
	/* only the second son is not empty */
	assert(size1 > 0);

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

  assert(size > 0);

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
  pavector  v;
  prealavector lambda;
  real     *x, *y;
  real      w;
  uint      i, j, k, size0, size1;

  pcluster  t;

  assert(size > 0);

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
    lambda = new_realavector(dim);

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
    del_realavector(lambda);
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
 * Statistics
 * ------------------------------------------------------------ */

uint
getdepth_cluster(pccluster t)
{
  uint      l, lt;
  uint      i;

  l = 0;
  if (t->son) {
    for (i = 0; i < t->sons; i++) {
      lt = getdepth_cluster(t->son[i]);
      l = UINT_MAX(lt, l);
    }
    l++;
  }

  return l;
}

uint
getmindepth_cluster(pccluster t)
{
  uint      l, lt;
  uint      i;

  l = 0;
  if (t->son) {
    l = getmindepth_cluster(t->son[0]);
    for (i = 1; i < t->sons; i++) {
      lt = getmindepth_cluster(t->son[i]);
      l = UINT_MIN(lt, l);
    }
    l++;
  }

  return l;
}

/* ------------------------------------------------------------
 * Structure Adaptation
 * ------------------------------------------------------------ */

void
extend_cluster(pcluster t, uint depth)
{
  uint      i;

  if ((t->son) && (depth > 0))
    for (i = 0; i < t->sons; i++)
      extend_cluster(t->son[i], depth - 1);

  else if ((!t->son) && (depth > 0)) {
    t->sons = 1;
    t->son = (pcluster *) allocmem(sizeof(pcluster));
    t->son[0] = new_cluster(t->size, t->idx, 0, t->dim);
    for (i = 0; i < t->dim; i++) {
      t->son[0]->bmin[i] = t->bmin[i];
      t->son[0]->bmax[i] = t->bmax[i];
    }
    extend_cluster(t->son[0], depth - 1);
  }
  update_cluster(t);
}

void
cut_cluster(pcluster t, uint depth)
{
  uint      i;

  if ((t->son) && (depth > 0))
    for (i = 0; i < t->sons; i++)
      cut_cluster(t->son[i], depth - 1);

  else if ((t->son) && (depth == 0)) {
    for (i = 0; i < t->sons; i++)
      del_cluster(t->son[i]);
    t->sons = 0;
    freemem(t->son);
    t->son = NULL;
  }
  update_cluster(t);
}

void
balance_cluster(pcluster t, uint depth)
{
  uint      i;

  if ((t->son) && (depth > 0))
    for (i = 0; i < t->sons; i++)
      balance_cluster(t->son[i], depth - 1);

  else if ((t->son) && (depth == 0))
    cut_cluster(t, depth);

  else if ((!t->son) && (depth > 0))
    extend_cluster(t, depth);

  update_cluster(t);
}

void
coarsen_cluster(pcluster t, uint minsize)
{
  bool      a = false;
  uint      i;

  for (i = 0; i < t->sons; i++) {
    if (t->son[i]->size < minsize) {
      a = true;
      break;
    }
  }

  if (a) {
    for (i = 0; i < t->sons; i++)
      del_cluster(t->son[i]);
    t->sons = 0;
    freemem(t->son);
    t->son = NULL;
  }
  else
    for (i = 0; i < t->sons; i++)
      coarsen_cluster(t->son[i], minsize);

  update_cluster(t);
}

void
setsons_cluster(pcluster t, uint sons)
{
  uint      i;

  assert(t->sons == 0);

  t->sons = sons;
  t->son = (pcluster *) allocmem(sizeof(pcluster) * sons);

  for (i = 0; i < sons; i++)
    t->son[i] = 0;
}

/* ------------------------------------------------------------
 * Bounding box
 * ------------------------------------------------------------ */

real
getdiam_2_cluster(pccluster t)
{
  real      diam2;
  uint      i;

  diam2 = 0.0;
  for (i = 0; i < t->dim; i++)
    diam2 += REAL_SQR(t->bmax[i] - t->bmin[i]);

  return REAL_SQRT(diam2);
}

real
getdist_2_cluster(pccluster t, pccluster s)
{
  real      dist2;
  uint      i;

  assert(t->dim == s->dim);

  dist2 = 0.0;
  for (i = 0; i < t->dim; i++)
    if (t->bmax[i] < s->bmin[i])	/* t to the "left" of s */
      dist2 += REAL_SQR(s->bmin[i] - t->bmax[i]);
    else if (s->bmax[i] < t->bmin[i])	/* t to the "right" of s */
      dist2 += REAL_SQR(t->bmin[i] - s->bmax[i]);

  return REAL_SQRT(dist2);
}

real
getdiam_max_cluster(pccluster t)
{

  real      diam_max, diam;
  uint      i;

  diam_max = 0.0;
  for (i = 0; i < t->dim; i++) {
    diam = t->bmax[i] - t->bmin[i];
    if (diam > diam_max) {
      diam_max = diam;
    }
  }

  return diam_max;
}

real
getdist_max_cluster(pccluster t, pccluster s)
{

  real      dist_max, dist;
  uint      i;

  assert(t->dim == s->dim);

  dist_max = 0.0;
  for (i = 0; i < t->dim; i++) {
    if (t->bmax[i] < s->bmin[i]) {
      dist = s->bmin[i] - t->bmax[i];
      if (dist > dist_max) {
	dist_max = dist;
      }
    }
    else if (s->bmax[i] < t->bmin[i]) {
      dist = t->bmin[i] - s->bmax[i];
      if (dist > dist_max) {
	dist_max = dist;
      }
    }
  }

  return dist_max;
}

/* ------------------------------------------------------------
 * Hierarchical iterator
 * ------------------------------------------------------------ */

void
iterate_cluster(pccluster t, uint tname,
		void (*pre) (pccluster t, uint tname, void *data),
		void (*post) (pccluster t, uint tname, void *data),
		void *data)
{
  uint      tname1;
  uint      i;

  if (pre)
    pre(t, tname, data);

  tname1 = tname + 1;
  for (i = 0; i < t->sons; i++) {
    iterate_cluster(t->son[i], tname1, pre, post, data);

    tname1 += t->son[i]->desc;
  }
  assert(tname1 == tname + t->desc);

  if (post)
    post(t, tname, data);
}

void
iterate_parallel_cluster(pccluster t, uint tname, uint pardepth,
			 void (*pre) (pccluster t, uint tname, void *data),
			 void (*post) (pccluster t, uint tname, void *data),
			 void *data)
{
  uint     *tname1, tname2;
#ifdef USE_OPENMP
  uint      nthreads;		/* HACK: Solaris workaround */
#endif
  uint      i;

  if (pre)
    pre(t, tname, data);

  if (t->sons > 0) {
    tname1 = (uint *) allocmem((size_t) sizeof(uint) * t->sons);

    tname2 = tname + 1;
    for (i = 0; i < t->sons; i++) {
      tname1[i] = tname2;
      tname2 += t->son[i]->desc;
    }
    assert(tname2 == tname + t->desc);

#ifdef USE_OPENMP
    nthreads = t->sons;
    (void) nthreads;
#pragma omp parallel for if(pardepth > 0), num_threads(nthreads)
#endif
    for (i = 0; i < t->sons; i++)
      iterate_parallel_cluster(t->son[i], tname1[i],
			       (pardepth > 0 ? pardepth - 1 : 0), pre, post,
			       data);

    freemem(tname1);
  }

  if (post)
    post(t, tname, data);
}

/* ------------------------------------------------------------
 * Enumeration
 * ------------------------------------------------------------ */

static void
enumerate(pcluster t, uint tname, pcluster *tn)
{
  uint      tname1;
  uint      i;

  tn[tname] = t;

  tname1 = tname + 1;
  for (i = 0; i < t->sons; i++) {
    enumerate(t->son[i], tname1, tn);

    tname1 += t->son[i]->desc;
  }
  assert(tname1 == tname + t->desc);
}

pcluster *
enumerate_cluster(pcluster t)
{
  pcluster *tn;

  tn = (pcluster *) allocmem((size_t) sizeof(pcluster) * t->desc);

  enumerate(t, 0, tn);

  return tn;
}

/* ------------------------------------------------------------
 * File I/O
 * ------------------------------------------------------------ */

#ifdef USE_NETCDF
static void
write_count(pccluster t, size_t * clusters, size_t * coeffs)
{
  uint      i;

  /* Increase cluster counter */
  (*clusters)++;

  /* Add number of coefficients of transfer matrix */
  (*coeffs) += 2 * t->dim;

  /* Handle sons */
  for (i = 0; i < t->sons; i++)
    write_count(t->son[i], clusters, coeffs);
}

static void
write_cdf(pccluster t,
	  size_t clusters, size_t coeffs,
	  size_t * clusteridx, size_t * coeffidx,
	  int nc_file, int nc_sons, int nc_size, int nc_coeff)
{
  size_t    start, count;
  ptrdiff_t stride;
  int       val, result;
  uint      i;

  assert(*clusteridx <= clusters);

  /* Write number of sons to nc_sons[*clusteridx] */
  start = *clusteridx;
  count = 1;
  stride = 1;
  val = t->sons;
  result = nc_put_vars(nc_file, nc_sons, &start, &count, &stride, &val);
  assert(result == NC_NOERR);

  /* Write size of cluster to nc_size[*clusteridx] */
  val = t->size;
  result = nc_put_vars(nc_file, nc_size, &start, &count, &stride, &val);
  assert(result == NC_NOERR);

  /* Increase cluster index */
  (*clusteridx)++;

  /* Handle sons */
  for (i = 0; i < t->sons; i++)
    write_cdf(t->son[i], clusters, coeffs, clusteridx, coeffidx,
	      nc_file, nc_sons, nc_size, nc_coeff);

  /* Write bounding box */
  start = *coeffidx;
  assert(start + 2 * t->dim <= coeffs);
  count = t->dim;
  result = nc_put_vars(nc_file, nc_coeff, &start, &count, &stride, t->bmin);
  assert(result == NC_NOERR);
  start += t->dim;

  result = nc_put_vars(nc_file, nc_coeff, &start, &count, &stride, t->bmax);
  assert(result == NC_NOERR);
  start += t->dim;
  (*coeffidx) = start;
}

void
write_cdf_cluster(pccluster t, const char *name)
{
  size_t    clusters, clusteridx;
  size_t    coeffs, coeffidx;
  int       nc_file, nc_sons, nc_size, nc_idx, nc_coeff;
  int       nc_clusters, nc_coeffs, nc_totalsize, nc_dim;
  int       result;

  /* Count number of clusters and coefficients */
  clusters = 0;
  coeffs = 0;
  write_count(t, &clusters, &coeffs);

  /* Create NetCDF file */
  result = nc_create(name, NC_64BIT_OFFSET, &nc_file);
  assert(result == NC_NOERR);

  /* Define "clusters" dimension */
  result = nc_def_dim(nc_file, "clusters", clusters, &nc_clusters);
  assert(result == NC_NOERR);

  /* Define "coeffs" dimension */
  result = nc_def_dim(nc_file, "coeffs", coeffs, &nc_coeffs);
  assert(result == NC_NOERR);

  /* Define "totalsize" dimension */
  result = nc_def_dim(nc_file, "totalsize", t->size, &nc_totalsize);
  assert(result == NC_NOERR);

  /* Define "dim" dimension */
  result = nc_def_dim(nc_file, "dim", t->dim, &nc_dim);
  assert(result == NC_NOERR);

  /* Define "sons" variable */
  result = nc_def_var(nc_file, "sons", NC_INT, 1, &nc_clusters, &nc_sons);
  assert(result == NC_NOERR);

  /* Define "size" variable */
  result = nc_def_var(nc_file, "size", NC_INT, 1, &nc_clusters, &nc_size);
  assert(result == NC_NOERR);

  /* Define "idx" variable */
  result = nc_def_var(nc_file, "idx", NC_INT, 1, &nc_totalsize, &nc_idx);
  assert(result == NC_NOERR);

  /* Define "coeff" variable */
  result = nc_def_var(nc_file, "coeff", NC_DOUBLE, 1, &nc_coeffs, &nc_coeff);
  assert(result == NC_NOERR);

  /* Finish NetCDF define mode */
  result = nc_enddef(nc_file);
  assert(result == NC_NOERR);

  /* Write index to NetCDF variable */
  result = nc_put_var(nc_file, nc_idx, t->idx);

  /* Write coefficiens to NetCDF variables */
  clusteridx = 0;
  coeffidx = 0;
  write_cdf(t, clusters, coeffs, &clusteridx, &coeffidx,
	    nc_file, nc_sons, nc_size, nc_coeff);
  assert(clusteridx == clusters);
  assert(coeffidx == coeffs);

  /* Close file */
  result = nc_close(nc_file);
  assert(result == NC_NOERR);
}

static void
prefix_name(char *buf, int bufsize, const char *prefix, const char *name)
{
  if (prefix)
    snprintf(buf, bufsize, "%s_%s", prefix, name);
  else
    strncpy(buf, name, bufsize);
}

void
write_cdfpart_cluster(pccluster t, int nc_file, const char *prefix)
{
  size_t    clusters, clusteridx;
  size_t    coeffs, coeffidx;
  char     *buf;
  int       bufsize;
  int       nc_sons, nc_size, nc_idx, nc_coeff;
  int       nc_clusters, nc_coeffs, nc_totalsize, nc_dim;
  int       result;

  /* Prepare buffer for prefixed names */
  bufsize = strlen(prefix) + 16;
  buf = (char *) allocmem(sizeof(char) * bufsize);

  /* Count number of clusters and coefficients */
  clusters = 0;
  coeffs = 0;
  write_count(t, &clusters, &coeffs);

  /* Switch NetCDF file to define mode */
  result = nc_redef(nc_file);
  assert(result == NC_NOERR || result == NC_EINDEFINE);

  /* Define "clusters" dimension */
  prefix_name(buf, bufsize, prefix, "clusters");
  result = nc_def_dim(nc_file, buf, clusters, &nc_clusters);
  assert(result == NC_NOERR);

  /* Define "coeffs" dimension */
  prefix_name(buf, bufsize, prefix, "coeffs");
  result = nc_def_dim(nc_file, buf, coeffs, &nc_coeffs);
  assert(result == NC_NOERR);

  /* Define "totalsize" dimension */
  prefix_name(buf, bufsize, prefix, "totalsize");
  result = nc_def_dim(nc_file, buf, t->size, &nc_totalsize);
  assert(result == NC_NOERR);

  /* Define "dim" dimension */
  prefix_name(buf, bufsize, prefix, "dim");
  result = nc_def_dim(nc_file, buf, t->dim, &nc_dim);
  assert(result == NC_NOERR);

  /* Define "sons" variable */
  prefix_name(buf, bufsize, prefix, "sons");
  result = nc_def_var(nc_file, buf, NC_INT, 1, &nc_clusters, &nc_sons);
  assert(result == NC_NOERR);

  /* Define "size" variable */
  prefix_name(buf, bufsize, prefix, "size");
  result = nc_def_var(nc_file, buf, NC_INT, 1, &nc_clusters, &nc_size);
  assert(result == NC_NOERR);

  /* Define "idx" variable */
  prefix_name(buf, bufsize, prefix, "idx");
  result = nc_def_var(nc_file, buf, NC_INT, 1, &nc_totalsize, &nc_idx);
  assert(result == NC_NOERR);

  /* Define "coeff" variable */
  prefix_name(buf, bufsize, prefix, "coeff");
  result = nc_def_var(nc_file, buf, NC_DOUBLE, 1, &nc_coeffs, &nc_coeff);
  assert(result == NC_NOERR);

  /* Finish NetCDF define mode */
  result = nc_enddef(nc_file);
  assert(result == NC_NOERR);

  /* Write index to NetCDF variable */
  result = nc_put_var(nc_file, nc_idx, t->idx);

  /* Write coefficiencs to NetCDF variables */
  clusteridx = 0;
  coeffidx = 0;
  write_cdf(t, clusters, coeffs, &clusteridx, &coeffidx,
	    nc_file, nc_sons, nc_size, nc_coeff);
  assert(clusteridx == clusters);
  assert(coeffidx == coeffs);

  /* Clean up */
  nc_sync(nc_file);
  freemem(buf);
}

static pcluster
read_cdf_part(int nc_file, size_t clusters, size_t coeffs,
	      int nc_sons, int nc_size, int nc_coeff,
	      uint * idx, int dim, size_t * clusteridx, size_t * coeffidx)
{
  pcluster  t, t1;
  uint     *idx1;
  uint      size;
  uint      sons;
  uint      i;
  size_t    start, count;
  ptrdiff_t stride;
  int       val, result;

  /* Get number of sons */
  start = *clusteridx;
  count = 1;
  stride = 1;
  result = nc_get_vars(nc_file, nc_sons, &start, &count, &stride, &val);
  assert(result == NC_NOERR);
  sons = val;

  /* Get size of cluster */
  result = nc_get_vars(nc_file, nc_size, &start, &count, &stride, &val);
  assert(result == NC_NOERR);
  size = val;

  /* Create new cluster */
  t = new_cluster(size, idx, sons, dim);

  /* Increase cluster index */
  (*clusteridx)++;

  /* Handle sons */
  if (sons > 0) {
    idx1 = idx;
    for (i = 0; i < sons; i++) {
      t1 = read_cdf_part(nc_file, clusters, coeffs,
			 nc_sons, nc_size, nc_coeff,
			 idx1, dim, clusteridx, coeffidx);
      t->son[i] = t1;

      idx1 += t1->size;
    }
    assert(idx1 == idx + size);
  }

  /* Get bounding box */
  start = (*coeffidx);
  count = dim;
  result = nc_get_vars(nc_file, nc_coeff, &start, &count, &stride, t->bmin);
  start += dim;

  result = nc_get_vars(nc_file, nc_coeff, &start, &count, &stride, t->bmax);
  start += dim;
  (*coeffidx) = start;

  /* Finish initialization */
  update_cluster(t);

  return t;
}

pcluster
read_cdf_cluster(const char *name)
{
  pcluster  t;
  uint     *idx;
  size_t    clusters, clusteridx;
  size_t    coeffs, coeffidx;
  char      dimname[NC_MAX_NAME + 1];
  int       nc_file, nc_sons, nc_size, nc_idx, nc_coeff;
  int       nc_clusters, nc_coeffs, nc_totalsize, nc_dim;
  size_t    dim, totalsize;
  int       result;

  /* Open NetCDF file */
  result = nc_open(name, NC_NOWRITE, &nc_file);
  assert(result == NC_NOERR);

  /* Get "clusters" dimension */
  result = nc_inq_dimid(nc_file, "clusters", &nc_clusters);
  assert(result == NC_NOERR);
  result = nc_inq_dim(nc_file, nc_clusters, dimname, &clusters);
  assert(result == NC_NOERR);

  /* Get "coeffs" dimension */
  result = nc_inq_dimid(nc_file, "coeffs", &nc_coeffs);
  assert(result == NC_NOERR);
  result = nc_inq_dim(nc_file, nc_coeffs, dimname, &coeffs);
  assert(result == NC_NOERR);

  /* Get "totalsize" dimension */
  result = nc_inq_dimid(nc_file, "totalsize", &nc_totalsize);
  assert(result == NC_NOERR);
  result = nc_inq_dim(nc_file, nc_totalsize, dimname, &totalsize);
  assert(result == NC_NOERR);

  /* Get "dim" dimension */
  result = nc_inq_dimid(nc_file, "dim", &nc_dim);
  assert(result == NC_NOERR);
  result = nc_inq_dim(nc_file, nc_dim, dimname, &dim);
  assert(result == NC_NOERR);

  /* Get "sons" variable */
  result = nc_inq_varid(nc_file, "sons", &nc_sons);
  assert(result == NC_NOERR);

  /* Get "size" variable */
  result = nc_inq_varid(nc_file, "size", &nc_size);
  assert(result == NC_NOERR);

  /* Get "idx" variable */
  result = nc_inq_varid(nc_file, "idx", &nc_idx);
  assert(result == NC_NOERR);

  /* Get "coeff" variable */
  result = nc_inq_varid(nc_file, "coeff", &nc_coeff);
  assert(result == NC_NOERR);

  /* Read index */
  idx = (uint *) allocmem(sizeof(uint) * totalsize);
  nc_get_var(nc_file, nc_idx, idx);

  /* Read coefficients from NetCDF variables */
  clusteridx = 0;
  coeffidx = 0;
  t = read_cdf_part(nc_file, clusters, coeffs,
		    nc_sons, nc_size, nc_coeff,
		    idx, dim, &clusteridx, &coeffidx);
  assert(clusteridx == clusters);
  assert(coeffidx == coeffs);

  /* Close NetCDF file */
  nc_close(nc_file);

  return t;
}

pcluster
read_cdfpart_cluster(int nc_file, const char *prefix)
{
  pcluster  t;
  uint     *idx;
  size_t    clusters, clusteridx;
  size_t    coeffs, coeffidx;
  char      dimname[NC_MAX_NAME + 1];
  char     *buf;
  int       bufsize;
  int       nc_sons, nc_size, nc_idx, nc_coeff;
  int       nc_clusters, nc_coeffs, nc_totalsize, nc_dim;
  size_t    dim, totalsize;
  int       result;

  /* Prepare buffer for prefixed names */
  bufsize = strlen(prefix) + 16;
  buf = (char *) allocmem(sizeof(char) * bufsize);

  /* Get "clusters" dimension */
  prefix_name(buf, bufsize, prefix, "clusters");
  result = nc_inq_dimid(nc_file, buf, &nc_clusters);
  assert(result == NC_NOERR);
  result = nc_inq_dim(nc_file, nc_clusters, dimname, &clusters);
  assert(result == NC_NOERR);

  /* Get "coeffs" dimension */
  prefix_name(buf, bufsize, prefix, "coeffs");
  result = nc_inq_dimid(nc_file, buf, &nc_coeffs);
  assert(result == NC_NOERR);
  result = nc_inq_dim(nc_file, nc_coeffs, dimname, &coeffs);
  assert(result == NC_NOERR);

  /* Get "totalsize" dimension */
  prefix_name(buf, bufsize, prefix, "totalsize");
  result = nc_inq_dimid(nc_file, buf, &nc_totalsize);
  assert(result == NC_NOERR);
  result = nc_inq_dim(nc_file, nc_totalsize, dimname, &totalsize);
  assert(result == NC_NOERR);

  /* Get "dim" dimension */
  prefix_name(buf, bufsize, prefix, "dim");
  result = nc_inq_dimid(nc_file, buf, &nc_dim);
  assert(result == NC_NOERR);
  result = nc_inq_dim(nc_file, nc_dim, dimname, &dim);
  assert(result == NC_NOERR);

  /* Get "sons" variable */
  prefix_name(buf, bufsize, prefix, "sons");
  result = nc_inq_varid(nc_file, buf, &nc_sons);
  assert(result == NC_NOERR);

  /* Get "size" variable */
  prefix_name(buf, bufsize, prefix, "size");
  result = nc_inq_varid(nc_file, buf, &nc_size);
  assert(result == NC_NOERR);

  /* Get "idx" variable */
  prefix_name(buf, bufsize, prefix, "idx");
  result = nc_inq_varid(nc_file, buf, &nc_idx);
  assert(result == NC_NOERR);

  /* Get "coeff" variable */
  prefix_name(buf, bufsize, prefix, "coeff");
  result = nc_inq_varid(nc_file, buf, &nc_coeff);
  assert(result == NC_NOERR);

  /* Read index */
  idx = (uint *) allocmem(sizeof(uint) * totalsize);
  nc_get_var(nc_file, nc_idx, idx);

  /* Read coefficients from NetCDF variables */
  clusteridx = 0;
  coeffidx = 0;
  t = read_cdf_part(nc_file, clusters, coeffs,
		    nc_sons, nc_size, nc_coeff,
		    idx, dim, &clusteridx, &coeffidx);
  assert(clusteridx == clusters);
  assert(coeffidx == coeffs);

  /* Clean up */
  freemem(buf);

  return t;
}
#endif
