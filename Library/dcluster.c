
/* ------------------------------------------------------------
 * This is the file "dcluster.c" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2015
 * ------------------------------------------------------------ */

#include "dcluster.h"

#include "basic.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static uint active_dcluster = 0;

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

pdcluster
new_dcluster(uint size, uint * idx, uint sons, uint dim)
{
  pdcluster t;
  uint      i;

  t = (pdcluster) allocmem(sizeof(dcluster));
  t->size = size;
  t->idx = idx;
  t->sons = sons;
  t->dim = dim;
  t->desc = 0;

  t->bmin = allocreal(dim);
  t->bmax = allocreal(dim);

  t->directions = 0;
  t->dir = 0;
  t->dirson = 0;

  t->sons = sons;
  t->son = 0;
  if (sons > 0) {
    t->son = (pdcluster *) allocmem(sizeof(pdcluster) * sons);

    for (i = 0; i < sons; i++)
      t->son[i] = 0;
  }

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_dcluster++;

  return t;
}

void
update_dcluster(pdcluster t)
{
  uint      desc;
  uint      i;

  desc = 1;
  for (i = 0; i < t->sons; i++)
    desc += t->son[i]->desc;

  t->desc = desc;
}

void
del_dcluster(pdcluster t)
{
  uint      i;

  freemem(t->bmin);
  freemem(t->bmax);

  if (t->sons > 0) {
    for (i = 0; i < t->sons; i++)
      del_dcluster(t->son[i]);
    freemem(t->son);

    if (t->directions > 0) {
      for (i = 0; i < t->sons; i++)
	freemem(t->dirson[i]);
      freemem(t->dirson);
    }
  }

  freemem(t);

  assert(active_dcluster > 0);

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_dcluster--;
}

pleveldir
new_leveldir(uint depth, uint dim)
{
  pleveldir ld;
  uint      l;

  ld = (pleveldir) allocmem(sizeof(leveldir));
  ld->depth = depth;
  ld->dim = dim;
  ld->maxdiam = allocreal(depth + 1);
  ld->splits = (uint *) allocmem(sizeof(uint) * (depth + 1));
  ld->directions = (uint *) allocmem(sizeof(uint) * (depth + 1));
  ld->dir = (preal **) allocmem(sizeof(preal *) * (depth + 1));

  for (l = 0; l <= depth; l++) {
    ld->maxdiam[l] = 0.0;
    ld->splits[l] = 0;
    ld->directions[l] = 0;
    ld->dir[l] = 0;
  }
  ld->dirmem = 0;

  return ld;
}

void
del_leveldir(pleveldir ld)
{
  uint      depth;
  uint      l;

  depth = ld->depth;

  for (l = 0; l <= depth; l++)
    if (ld->directions[l]) {
      freemem(ld->dir[l]);
      ld->dir[l] = 0;
    }

  if (ld->dirmem) {
    freemem(ld->dirmem);
    ld->dirmem = 0;
  }

  freemem(ld->dir);
  freemem(ld->directions);
  freemem(ld->splits), freemem(ld->maxdiam);
  freemem(ld);
}

/* ------------------------------------------------------------
 * Construction of dcluster trees
 * ------------------------------------------------------------ */

pdcluster
buildfromcluster_dcluster(pccluster t)
{
  pdcluster s;
  uint      i;

  s = new_dcluster(t->size, t->idx, t->sons, t->dim);

  for (i = 0; i < t->dim; i++) {
    s->bmin[i] = t->bmin[i];
    s->bmax[i] = t->bmax[i];
  }

  for (i = 0; i < t->sons; i++)
    s->son[i] = buildfromcluster_dcluster(t->son[i]);

  update_dcluster(s);

  return s;
}

/* ------------------------------------------------------------
 * Bounding box
 * ------------------------------------------------------------ */

real
diam_dcluster(pcdcluster t)
{
  real      diam2;
  uint      i;

  diam2 = 0.0;
  for (i = 0; i < t->dim; i++)
    diam2 += REAL_SQR(t->bmax[i] - t->bmin[i]);

  return REAL_SQRT(diam2);
}

real
dist_dcluster(pcdcluster t, pcdcluster s)
{
  real      dist2;
  uint      i;

  assert(t->dim == s->dim);

  dist2 = 0.0;
  for (i = 0; i < t->dim; i++)
    if (t->bmax[i] < s->bmin[i])
      dist2 += REAL_SQR(s->bmin[i] - t->bmax[i]);
    else if (s->bmax[i] < t->bmin[i])
      dist2 += REAL_SQR(t->bmin[i] - s->bmax[i]);

  return REAL_SQRT(dist2);
}

real
middist_dcluster(pcdcluster t, pcdcluster s)
{
  real      dist2;
  uint      i;

  assert(t->dim == s->dim);

  dist2 = 0.0;
  for (i = 0; i < t->dim; i++)
    dist2 += REAL_SQR(0.5 * (t->bmax[i] + t->bmin[i]
			     - s->bmax[i] - s->bmin[i]));

  return REAL_SQRT(dist2);
}

/* ------------------------------------------------------------
 * Match directions
 * ------------------------------------------------------------ */

uint
finddirection_dcluster(pcdcluster t, real alpha, pcreal d)
{
  real      norm, normmin;
  uint      dim = t->dim;
  uint      i, k, kmin;

  normmin = 4.0;
  kmin = 0;

  for (k = 0; k < t->directions; k++) {
    norm = REAL_SQR(alpha * d[0] - t->dir[k][0]);
    for (i = 1; i < dim; i++)
      norm += REAL_SQR(alpha * d[i] - t->dir[k][i]);

    if (norm < normmin) {
      normmin = norm;
      kmin = k;
    }
  }

  return kmin;
}

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

uint
getactives_dcluster()
{
  return active_dcluster;
}

size_t
getsize_dcluster(pcdcluster t)
{
  size_t    sz;
  uint      i;

  sz = (size_t) sizeof(dcluster) + 2 * t->dim * sizeof(real);
  if (t->sons > 0) {
    sz += (size_t) sizeof(pdcluster) * t->sons;
    for (i = 0; i < t->sons; i++)
      sz += getsize_dcluster(t->son[i]);
  }

  return sz;
}

uint
getdepth_dcluster(pcdcluster t)
{
  uint      d, d1;
  uint      i;

  d = 0;

  if (t->sons > 0) {
    for (i = 0; i < t->sons; i++) {
      d1 = getdepth_dcluster(t->son[i]);
      if (d1 > d)
	d = d1;
    }
    d++;
  }

  return d;
}

uint
getalldirections_dcluster(pcdcluster t)
{
  uint      dirs;
  uint      i;

  dirs = t->directions;

  for (i = 0; i < t->sons; i++)
    dirs += getalldirections_dcluster(t->son[i]);

  return dirs;
}

/* ------------------------------------------------------------
 * Directions via the surface of boxes
 * ------------------------------------------------------------ */

static void
compute_maxdiam(pdcluster t, uint l, pleveldir ld)
{
  real      diam;
  uint      i;

  assert(l <= ld->depth);

  diam = diam_dcluster(t);

  if (diam > ld->maxdiam[l])
    ld->maxdiam[l] = diam;

  for (i = 0; i < t->sons; i++)
    compute_maxdiam(t->son[i], l + 1, ld);
}

static void
compute_splits(pleveldir ld, real eta1)
{
  uint      depth = ld->depth;
  real      r, invr;
  uint      l, splits;

  switch (ld->dim) {
  case 1:
    for (l = 0; l <= depth; l++)
      if (ld->maxdiam[l] <= 0.5 * eta1)
	ld->directions[l] = ld->splits[l] = 0;
      else {
	ld->splits[l] = 1;
	ld->directions[l] = 2;
      }
    break;

  case 2:
    for (l = 0; l <= depth; l++)
      if (ld->maxdiam[l] <= 0.5 * eta1)
	ld->directions[l] = ld->splits[l] = 0;
      else {
	r = eta1 / ld->maxdiam[l];

	invr = 1.0 / r;
	splits = (uint) invr;
	if (splits < invr)
	  splits++;

	ld->splits[l] = splits;
	ld->directions[l] = 4 * splits;
      }
    break;

  case 3:
    for (l = 0; l <= depth; l++)
      if (ld->maxdiam[l] <= 0.5 * eta1)
	ld->directions[l] = ld->splits[l] = 0;
      else {
	r = eta1 * M_SQRT1_2 / ld->maxdiam[l];

	invr = 1.0 / r;
	splits = (uint) invr;
	if (splits < invr)
	  splits++;

	ld->splits[l] = splits;
	ld->directions[l] = 6 * splits * splits;
      }
    break;
  }
}

static void
compute_directions(pleveldir ld)
{
  uint      depth = ld->depth;
  preal    *dir;
  real      inorm;
  uint      l, splits;
  uint      i, j, k;

  switch (ld->dim) {
  case 1:
    for (l = 0; l <= depth; l++) {
      ld->dir[l][0][0] = 1.0;
      ld->dir[l][1][0] = -1.0;
    }
    break;

  case 2:
    for (l = 0; l <= depth; l++) {
      splits = ld->splits[l];
      dir = ld->dir[l];

      /* Create directions on the boundary of the square [-1,1]^2 */
      k = 0;

      /* Right side 0 */
      for (i = 0; i < splits; i++) {
	dir[k][0] = 1.0;
	dir[k][1] = -1.0 + 2.0 * (i + 0.5) / splits;
	k++;
      }

      /* Top side 1 */
      for (i = 0; i < splits; i++) {
	dir[k][0] = -1.0 + 2.0 * (i + 0.5) / splits;
	dir[k][1] = 1.0;
	k++;
      }

      /* Left side 2 */
      for (i = 0; i < splits; i++) {
	dir[k][0] = -1.0;
	dir[k][1] = -1.0 + 2.0 * (i + 0.5) / splits;
	k++;
      }

      /* Bottom side 3 */
      for (i = 0; i < splits; i++) {
	dir[k][0] = -1.0 + 2.0 * (i + 0.5) / splits;
	dir[k][1] = -1.0;
	k++;
      }

      assert(k == ld->directions[l]);

      /* Project to unit circle */
      for (k = 0; k < ld->directions[l]; k++) {
	inorm = 1.0 / REAL_SQRT(REAL_SQR(dir[k][0]) +
				REAL_SQR(dir[k][1]) + REAL_SQR(dir[k][2]));
	dir[k][0] *= inorm;
	dir[k][1] *= inorm;
      }
    }
    break;

  case 3:
    for (l = 0; l <= depth; l++) {
      splits = ld->splits[l];
      dir = ld->dir[l];

      /* Create directions on the boundary of the cube [-1,1]^3 */
      k = 0;

      /* Front face 0 */
      for (i = 0; i < splits; i++)
	for (j = 0; j < splits; j++) {
	  dir[k][2] = 1.0;
	  dir[k][0] = -1.0 + 2.0 * (i + 0.5) / splits;
	  dir[k][1] = -1.0 + 2.0 * (j + 0.5) / splits;
	  k++;
	}

      /* Back face 1 */
      for (i = 0; i < splits; i++)
	for (j = 0; j < splits; j++) {
	  dir[k][2] = -1.0;
	  dir[k][0] = -1.0 + 2.0 * (i + 0.5) / splits;
	  dir[k][1] = -1.0 + 2.0 * (j + 0.5) / splits;
	  k++;
	}

      /* Right face 2 */
      for (i = 0; i < splits; i++)
	for (j = 0; j < splits; j++) {
	  dir[k][0] = 1.0;
	  dir[k][1] = -1.0 + 2.0 * (i + 0.5) / splits;
	  dir[k][2] = -1.0 + 2.0 * (j + 0.5) / splits;
	  k++;
	}

      /* Left face 3 */
      for (i = 0; i < splits; i++)
	for (j = 0; j < splits; j++) {
	  dir[k][0] = -1.0;
	  dir[k][1] = -1.0 + 2.0 * (i + 0.5) / splits;
	  dir[k][2] = -1.0 + 2.0 * (j + 0.5) / splits;
	  k++;
	}

      /* Top face 4 */
      for (i = 0; i < splits; i++)
	for (j = 0; j < splits; j++) {
	  dir[k][1] = 1.0;
	  dir[k][0] = -1.0 + 2.0 * (i + 0.5) / splits;
	  dir[k][2] = -1.0 + 2.0 * (j + 0.5) / splits;
	  k++;
	}

      /* Bottom face 5 */
      for (i = 0; i < splits; i++)
	for (j = 0; j < splits; j++) {
	  dir[k][1] = -1.0;
	  dir[k][0] = -1.0 + 2.0 * (i + 0.5) / splits;
	  dir[k][2] = -1.0 + 2.0 * (j + 0.5) / splits;
	  k++;
	}

      assert(k == ld->directions[l]);

      /* Project to unit sphere */
      for (k = 0; k < ld->directions[l]; k++) {
	inorm = 1.0 / REAL_SQRT(REAL_SQR(dir[k][0]) +
				REAL_SQR(dir[k][1]) + REAL_SQR(dir[k][2]));
	dir[k][0] *= inorm;
	dir[k][1] *= inorm;
	dir[k][2] *= inorm;
      }
    }
    break;
  }
}

static void
set_directions(pcleveldir ld, uint l, pdcluster t)
{
  uint      iota, j;

  assert(l <= ld->depth);

  t->directions = ld->directions[l];
  t->dir = (pcreal *) ld->dir[l];
  t->dirson = 0;

  if (t->sons > 0 && t->directions > 0) {
    t->dirson = (uint **) allocmem(sizeof(uint *) * t->sons);

    for (j = 0; j < t->sons; j++) {
      set_directions(ld, l + 1, t->son[j]);

      t->dirson[j] = (uint *) allocmem(sizeof(uint) * t->directions);
      for (iota = 0; iota < t->directions; iota++)
	/*
	   t->dirson[j][iota] = (t->son[j]->directions > 0 ?
	   finddirection_dcluster(t->son[j], 1.0, t->dir[iota]) :
	   0);
	 */
	t->dirson[j][iota] = (t->son[j]->directions > 0 ?
			      finddirection_leveldir(ld, l + 1, 1.0,
						     t->dir[iota]) : 0);
    }
  }
}

pleveldir
builddirections_box_dcluster(pdcluster t, real eta1)
{
  pleveldir ld;
  preal     dirmem;
  uint      depth;
  uint      l, iota;
  size_t    dirsz;

  depth = getdepth_dcluster(t);

  ld = new_leveldir(depth, t->dim);

  compute_maxdiam(t, 0, ld);

  compute_splits(ld, eta1);

  dirsz = 0;
  for (l = 0; l <= depth; l++)
    dirsz += ld->directions[l] * ld->dim;

  ld->dirmem = dirmem = allocreal(dirsz);

  for (l = 0; l <= depth; l++) {
    ld->dir[l] = (preal *) allocmem(sizeof(preal) * ld->directions[l]);
    for (iota = 0; iota < ld->directions[l]; iota++) {
      ld->dir[l][iota] = dirmem;
      dirmem += ld->dim;
    }
  }
  assert(dirmem == ld->dirmem + dirsz);

  compute_directions(ld);

  set_directions(ld, 0, t);

  return ld;
}

uint
finddirection_leveldir(pcleveldir ld, uint l, real alpha, pcreal d)
{
  real      ad[3];
#ifndef NDEBUG
  real      inorm;
#endif
  uint      splits = ld->splits[l];
  uint      face;
  uint      i, j, k;

  k = 0;

  /* Quick exit if there are no directions */
  if (splits == 0)
    return 0;

  switch (ld->dim) {
  case 1:
    if (alpha * d[0] < 0.0)
      k = 1;
    else
      k = 0;
    break;

  case 2:
    ad[0] = REAL_ABS(alpha * d[0]);
    ad[1] = REAL_ABS(alpha * d[1]);

    if (ad[0] >= ad[1]) {
      /* ad[0] is maximal, use left or right side */
      ad[1] = alpha * d[1] / ad[0];
      if (alpha * d[0] < 0) {
	/* Left side 2 */
	face = 2;
	i = (uint) (0.5 * (ad[1] + 1.0) * splits);
	ad[0] = -1.0;
      }
      else {
	/* Right side 0 */
	face = 0;
	i = (uint) (0.5 * (ad[1] + 1.0) * splits);
	ad[0] = 1.0;
      }
    }
    else {
      /* ad[1] is maximal, use bottom or top side */
      ad[0] = alpha * d[0] / ad[1];
      if (alpha * d[1] < 0) {
	/* Bottom side 3 */
	face = 3;
	i = (uint) (0.5 * (ad[0] + 1.0) * splits);
	ad[1] = -1.0;
      }
      else {
	/* Top side 1 */
	face = 1;
	i = (uint) (0.5 * (ad[0] + 1.0) * splits);
	ad[1] = 1.0;
      }
    }

    /* Handle special case */
    if (i == splits)
      i--;

    assert(i < splits);

    k = i + splits * face;

    /* Check if the index matches the correct direction */
#ifndef NDEBUG
    inorm = 1.0 / REAL_SQRT(REAL_SQR(d[0]) + REAL_SQR(d[1]));
    ad[0] = alpha * d[0] * inorm;
    ad[1] = alpha * d[1] * inorm;

    assert(REAL_ABS(ld->dir[l][k][0] - ad[0]) <= 1.0 / splits);
    assert(REAL_ABS(ld->dir[l][k][1] - ad[1]) <= 1.0 / splits);
#endif
    break;

  case 3:
    ad[0] = REAL_ABS(alpha * d[0]);
    ad[1] = REAL_ABS(alpha * d[1]);
    ad[2] = REAL_ABS(alpha * d[2]);

    if (ad[0] >= ad[1])
      if (ad[0] >= ad[2]) {
	/* ad[0] is maximal, use left or right face */
	ad[1] = alpha * d[1] / ad[0];
	ad[2] = alpha * d[2] / ad[0];
	if (alpha * d[0] < 0) {
	  /* Left face 3 */
	  face = 3;
	  i = (uint) (0.5 * (ad[1] + 1.0) * splits);
	  j = (uint) (0.5 * (ad[2] + 1.0) * splits);
	  ad[0] = -1.0;
	}
	else {
	  /* Right face 2 */
	  face = 2;
	  i = (uint) (0.5 * (ad[1] + 1.0) * splits);
	  j = (uint) (0.5 * (ad[2] + 1.0) * splits);
	  ad[0] = 1.0;
	}
      }
      else {
	/* ad[2] is maximal, use back or front face */
	ad[0] = alpha * d[0] / ad[2];
	ad[1] = alpha * d[1] / ad[2];
	if (alpha * d[2] < 0) {
	  /* Back face 1 */
	  face = 1;
	  i = (uint) (0.5 * (ad[0] + 1.0) * splits);
	  j = (uint) (0.5 * (ad[1] + 1.0) * splits);
	  ad[2] = -1.0;
	}
	else {
	  /* Front face 0 */
	  face = 0;
	  i = (uint) (0.5 * (ad[0] + 1.0) * splits);
	  j = (uint) (0.5 * (ad[1] + 1.0) * splits);
	  ad[2] = 1.0;
	}
      }
    else if (ad[1] >= ad[2]) {
      /* ad[1] is maximal, use bottom or top face */
      ad[0] = alpha * d[0] / ad[1];
      ad[2] = alpha * d[2] / ad[1];
      if (alpha * d[1] < 0) {
	/* Bottom face 5 */
	face = 5;
	i = (uint) (0.5 * (ad[0] + 1.0) * splits);
	j = (uint) (0.5 * (ad[2] + 1.0) * splits);
	ad[1] = -1.0;
      }
      else {
	/* Top face 4 */
	face = 4;
	i = (uint) (0.5 * (ad[0] + 1.0) * splits);
	j = (uint) (0.5 * (ad[2] + 1.0) * splits);
	ad[1] = 1.0;
      }
    }
    else {
      /* ad[2] is maximal, use back or front face */
      ad[0] = alpha * d[0] / ad[2];
      ad[1] = alpha * d[1] / ad[2];
      if (alpha * d[2] < 0) {
	/* Back face 1 */
	face = 1;
	i = (uint) (0.5 * (ad[0] + 1.0) * splits);
	j = (uint) (0.5 * (ad[1] + 1.0) * splits);
	ad[2] = -1.0;
      }
      else {
	/* Front face 0 */
	face = 0;
	i = (uint) (0.5 * (ad[0] + 1.0) * splits);
	j = (uint) (0.5 * (ad[1] + 1.0) * splits);
	ad[2] = 1.0;
      }
    }

    /* Handle special cases */
    if (i == splits)
      i--;
    if (j == splits)
      j--;

    assert(i < splits);
    assert(j < splits);

    k = j + splits * (i + splits * face);

#ifndef NDEBUG
    /* Check if the index matches the correct direction */
    inorm = 1.0 / REAL_SQRT(REAL_SQR(d[0]) + REAL_SQR(d[1]) + REAL_SQR(d[2]));
    ad[0] = alpha * d[0] * inorm;
    ad[1] = alpha * d[1] * inorm;
    ad[2] = alpha * d[2] * inorm;

    assert(REAL_ABS(ld->dir[l][k][0] - ad[0]) <= 1.0 / splits);
    assert(REAL_ABS(ld->dir[l][k][1] - ad[1]) <= 1.0 / splits);
    assert(REAL_ABS(ld->dir[l][k][2] - ad[2]) <= 1.0 / splits);
#endif
    break;
  }

  return k;
}
