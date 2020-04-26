

/* ------------------------------------------------------------
 * This is the file "kernelmatrix.c" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2018
 * ------------------------------------------------------------ */

#include "kernelmatrix.h"

#include "basic.h"

/* ------------------------------------------------------------
 * Create an empty @ref kernelmatrix object
 * ------------------------------------------------------------ */

pkernelmatrix
new_kernelmatrix(uint dim, uint points, uint m)
{
  pkernelmatrix km;
  real *x0;
  uint i;

  /* Initialize basic structure */
  km = (pkernelmatrix) allocmem(sizeof(kernelmatrix));
  km->dim = dim;
  km->points = points;
  km->m = m;

  /* Empty kernel callback function */
  km->kernel = 0;
  km->data = 0;

  /* Initialize arrays for point coordinates */
  km->x = (real **) allocmem(sizeof(real *) * points);
  km->x[0] = x0 = allocreal(points * dim);
  for(i=1; i<points; i++) {
    x0 += dim;
    km->x[i] = x0;
  }

  /* Initialize Chebyshev points */
  km->xi_ref = allocreal(m);
  for(i=0; i<m; i++)
    km->xi_ref[i] = REAL_COS(M_PI * (m - i - 0.5) / m);

  return km;
}

/* ------------------------------------------------------------
 * Delete a @ref kernelmatrix object
 * ------------------------------------------------------------ */

void
del_kernelmatrix(pkernelmatrix km)
{
  freemem(km->xi_ref);
  freemem(km->x[0]);
  freemem(km->x);
  freemem(km);
}

/* ------------------------------------------------------------
 * Create a clustergeometry object
 * ------------------------------------------------------------ */

pclustergeometry
creategeometry_kernelmatrix(pckernelmatrix km)
{
  uint dim = km->dim;
  uint points = km->points;
  const real **x = (const real **) km->x;
  pclustergeometry cg;
  uint i, j;

  cg = new_clustergeometry(dim, points);

  for(i=0; i<points; i++)
    for(j=0; j<dim; j++)
      cg->x[i][j] = cg->smin[i][j] = cg->smax[i][j] = x[i][j];

  return cg;
}

/* ------------------------------------------------------------
 * Fill a dense matrix
 * ------------------------------------------------------------ */

void
fillN_kernelmatrix(const uint *ridx, const uint *cidx, pckernelmatrix km,
		   pamatrix N)
{
  const real **x = (const real **) km->x;
  uint rows = N->rows;
  uint cols = N->cols;
  pfield Na = N->a;
  longindex ldN = N->ld;
  uint i, j, ii, jj;

  if(ridx) {
    if(cidx) {
      for(j=0; j<cols; j++) {
	jj = cidx[j];

	for(i=0; i<rows; i++) {
	  ii = ridx[i];

	  Na[i+j*ldN] = km->kernel(x[ii], x[jj], km->data);
	}
      }
    }
    else {
      assert(cols <= km->points);

      for(j=0; j<cols; j++) {
	for(i=0; i<rows; i++) {
	  ii = ridx[i];

	  Na[i+j*ldN] = km->kernel(x[ii], x[j], km->data);
	}
      }
    }
  }
  else {
    assert(rows <= km->points);

    if(cidx) {
      for(j=0; j<cols; j++) {
	jj = cidx[j];

	for(i=0; i<rows; i++)
	  Na[i+j*ldN] = km->kernel(x[i], x[jj], km->data);
      }
    }
    else {
      assert(cols <= km->points);

      for(j=0; j<cols; j++)
	for(i=0; i<rows; i++)
	  Na[i+j*ldN] = km->kernel(x[i], x[j], km->data);
    }
  }
}

/* ------------------------------------------------------------
 * Create transformed interpolation points for a bounding box
 * ------------------------------------------------------------ */

static real **
transform_points(uint dim, real *bmin, real *bmax,
		 uint m, const real *xi_ref)
{
  real **xi;
  real *xi0;
  real mid, rad;
  uint i, j;

  xi = (real **) allocmem(sizeof(real *) * dim);

  xi[0] = xi0 = (real *) allocmem(sizeof(real) * m * dim);
  for(i=1; i<dim; i++) {
    xi0 += m;
    xi[i] = xi0;
  }

  for(i=0; i<dim; i++) {
    mid = 0.5 * (bmax[i] + bmin[i]);
    rad = 0.5 * (bmax[i] - bmin[i]);

    for(j=0; j<m; j++)
      xi[i][j] = mid + rad * xi_ref[j];
  }

  return xi;
}

/* ------------------------------------------------------------
 * Fill a coupling matrix
 * ------------------------------------------------------------ */

static void
fillS_1d(uint dim, uint i0, uint j0,
	 const real **rxi, const real **cxi,
	 pckernelmatrix km,
	 real *xx, real *yy, pamatrix S)
{
  uint m = km->m;
  pfield Sa;
  longindex ldS;
  uint i, j;

  if(dim > 1) {
    dim--;

    for(i=0; i<m; i++) {
      xx[dim] = rxi[dim][i];

      for(j=0; j<m; j++) {
	yy[dim] = cxi[dim][j];

	fillS_1d(dim, i+i0*m, j+j0*m, rxi, cxi, km,
		 xx, yy, S);
      }
    }
  }
  else {
    Sa = S->a;
    ldS = S->ld;

    assert((m-1) + i0*m < S->rows);
    assert((m-1) + j0*m < S->cols);

    dim--;

    for(i=0; i<m; i++) {
      xx[dim] = rxi[dim][i];

      for(j=0; j<m; j++) {
	yy[dim] = cxi[dim][j];

	Sa[(i+i0*m)+(j+j0*m)*ldS] = km->kernel(xx, yy, km->data);
      }
    }
  }
}

static void
fillS_2d(uint dim, uint i0, uint j0,
	 const real **rxi, const real **cxi,
	 pckernelmatrix km,
	 real *xx, real *yy, pamatrix S)
{
  uint m = km->m;
  pfield Sa;
  longindex ldS;
  uint i1, j1, i2, j2;

  if(dim > 2) {
    dim -= 2;

    for(i1=0; i1<m; i1++) {
      xx[dim] = rxi[dim][i1];

      for(i2=0; i2<m; i2++) {
	xx[dim+1] = rxi[dim+1][i2];

	for(j1=0; j1<m; j1++) {
	  yy[dim] = cxi[dim][j1];

	  for(j2=0; j2<m; j2++) {
	    yy[dim+1] = cxi[dim+1][j2];

	    fillS_2d(dim, i1+m*(i2+m*i0), j1+m*(j2+m*j0), rxi, cxi, km,
		     xx, yy, S);
	  }
	}
      }
    }
  }
  else if(dim == 2) {
    Sa = S->a;
    ldS = S->ld;

    assert((m*m-1) + i0*m*m < S->rows);
    assert((m*m-1) + j0*m*m < S->cols);

    dim -= 2;

    for(i1=0; i1<m; i1++) {
      xx[dim] = rxi[dim][i1];

      for(i2=0; i2<m; i2++) {
	xx[dim+1] = rxi[dim+1][i2];

	for(j1=0; j1<m; j1++) {
	  yy[dim] = cxi[dim][j1];

	  for(j2=0; j2<m; j2++) {
	    yy[dim+1] = cxi[dim+1][j2];

	    Sa[(i1+m*(i2+m*i0)) + (j1+m*(j2+m*j0))*ldS] = km->kernel(xx, yy, km->data);
	  }
	}
      }
    }
  }
  else {
    assert(dim == 1);

    fillS_1d(dim, i0, j0, rxi, cxi, km, xx, yy, S);
  }
}

void
fillS_kernelmatrix(pccluster rc, pccluster cc,
		   pckernelmatrix km, pamatrix S)
{
  uint dim = km->dim;
  uint m = km->m;
  const real *xi_ref = km->xi_ref;
  real **rxi, **cxi;
  real *xx, *yy;

  assert(rc->dim == dim);
  assert(cc->dim == dim);

  /* Compute transformed interpolation points for the row cluster */
  rxi = transform_points(dim, rc->bmin, rc->bmax, m, xi_ref);

  /* Compute transformed interpolation points for the column cluster */
  cxi = transform_points(dim, cc->bmin, cc->bmax, m, xi_ref);

  /* Allocate storage for coordinate vectors */
  xx = (real *) allocmem(sizeof(real) * dim);
  yy = (real *) allocmem(sizeof(real) * dim);

  /* Fill S recursively by dimension */
  if(dim >= 2)
    fillS_2d(dim, 0, 0, (const real **) rxi, (const real **) cxi,
	     km, xx, yy, S);
  else
    fillS_1d(dim, 0, 0, (const real **) rxi, (const real **) cxi,
	     km, xx, yy, S);

  /* Clean up */
  freemem(yy);
  freemem(xx);
  freemem(cxi[0]);
  freemem(cxi);
  freemem(rxi[0]);
  freemem(rxi);
}

/* ------------------------------------------------------------
 * Evaluate a Lagrange polynomial
 * ------------------------------------------------------------ */

static real
eval_lagrange(uint m, const real *xi, uint i, real t)
{
  real d, n;
  uint j;

  d = n = 1.0;

  for(j=0; j<i; j++) {
    d *= t - xi[j];
    n *= xi[i] - xi[j];
  }

  for(j=i+1; j<m; j++) {
    d *= t - xi[j];
    n *= xi[i] - xi[j];
  }

  return d / n;
}

/* ------------------------------------------------------------
 * Fill a leaf matrix
 * ------------------------------------------------------------ */

static void
fillV_1d(uint dim, uint i, uint j0,
	 const real **txi,
	 pckernelmatrix km,
	 const real *xx, field alpha0, pamatrix V)
{
  uint m = km->m;
  pfield Va;
  longindex ldV;
  field alpha;
  uint j;

  if(dim > 0) {
    dim--;
    
    for(j=0; j<m; j++) {
      alpha = alpha0 * eval_lagrange(m, txi[dim], j, xx[dim]);

      fillV_1d(dim, i, j+j0*m, txi, km, xx, alpha, V);
    }
  }
  else {
    Va = V->a;
    ldV = V->ld;

    assert(i < V->rows);
    assert(j0 < V->cols);

    Va[i+j0*ldV] = alpha0;
  }
}

void
fillV_kernelmatrix(pccluster tc,
		   pckernelmatrix km, pamatrix V)
{
  uint dim = km->dim;
  uint m = km->m;
  const real *xi_ref = km->xi_ref;
  const real **x = (const real **) km->x;
  const uint *idx = tc->idx;
  real **txi;
  uint i;

  assert(tc->dim == dim);

  /* Compute transformed interpolation points for the cluster */
  txi = transform_points(dim, tc->bmin, tc->bmax, m, xi_ref);

  /* Fill V recursively by dimension, each row individually */
  for(i=0; i<tc->size; i++)
    fillV_1d(dim, i, 0, (const real **) txi, km, x[idx[i]], 1.0, V);

  /* Clean up */
  freemem(txi[0]);
  freemem(txi);
}

/* ------------------------------------------------------------
 * Fill a transfer matrix
 * ------------------------------------------------------------ */

static void
fillE_1d(uint dim, uint i0, uint j0,
	 const real **sxi, const real **fxi,
	 pckernelmatrix km,
	 field alpha0, pamatrix E)
{
  uint m = km->m;
  pfield Ea;
  longindex ldE;
  field alpha;
  uint i, j;

  if(dim > 0) {
    dim--;
    
    for(i=0; i<m; i++)
      for(j=0; j<m; j++) {
	alpha = alpha0 * eval_lagrange(m, fxi[dim], j, sxi[dim][i]);

	fillE_1d(dim, i+i0*m, j+j0*m, sxi, fxi, km, alpha, E);
      }
  }
  else {
    Ea = E->a;
    ldE = E->ld;

    assert(i0 < E->rows);
    assert(j0 < E->cols);

    Ea[i0+j0*ldE] = alpha0;
  }
}

void
fillE_kernelmatrix(pccluster sc, pccluster fc,
		   pckernelmatrix km, pamatrix E)
{
  uint dim = km->dim;
  uint m = km->m;
  const real *xi_ref = km->xi_ref;
  real **sxi, **fxi;

  assert(sc->dim == dim);
  assert(fc->dim == dim);

  /* Compute transformed interpolation points for the son cluster */
  sxi = transform_points(dim, sc->bmin, sc->bmax, m, xi_ref);

  /* Compute transformed interpolation points for the father cluster */
  fxi = transform_points(dim, fc->bmin, fc->bmax, m, xi_ref);

  /* Fill E recursively by dimension */
  fillE_1d(dim, 0, 0, (const real **) sxi, (const real **) fxi, km, 1.0, E);

  /* Clean up */
  freemem(fxi[0]);
  freemem(fxi);
  freemem(sxi[0]);
  freemem(sxi);
}

/* ------------------------------------------------------------
 * Fill a cluster basis
 * ------------------------------------------------------------ */

static void
fill_clusterbasis(pckernelmatrix km, uint k, pclusterbasis cb)
{
  uint i;

  for(i=0; i<cb->sons; i++)
    fill_clusterbasis(km, k, cb->son[i]);

  setrank_clusterbasis(cb, k);

  if(cb->sons > 0) {
    for(i=0; i<cb->sons; i++)
      fillE_kernelmatrix(cb->son[i]->t, cb->t, km, &cb->son[i]->E);
  }
  else
    fillV_kernelmatrix(cb->t, km, &cb->V);

  update_clusterbasis(cb);
}

void
fill_clusterbasis_kernelmatrix(pckernelmatrix km, pclusterbasis cb)
{
  uint dim = km->dim;
  uint m = km->m;
  uint k;
  uint i;

  k = 1;
  for(i=0; i<dim; i++)
    k *= m;

  fill_clusterbasis(km, k, cb);
}

/* ------------------------------------------------------------
 * Fill an H^2-matrix
 * ------------------------------------------------------------ */

void
fill_h2matrix_kernelmatrix(pckernelmatrix km, ph2matrix G)
{
  uint rsons, csons;
  uint i, j;

  if(G->son) {
    rsons = G->rsons;
    csons = G->csons;

    for(j=0; j<csons; j++)
      for(i=0; i<rsons; i++)
	fill_h2matrix_kernelmatrix(km, G->son[i+j*rsons]);
  }
  else if(G->u)
    fillS_kernelmatrix(G->rb->t, G->cb->t, km, &G->u->S);
  else {
    assert(G->f);
    fillN_kernelmatrix(G->rb->t->idx, G->cb->t->idx, km, G->f);
  }
}
