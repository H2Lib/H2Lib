
/* ------------------------------------------------------------
 * This is the file "ie1d.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2016
 * ------------------------------------------------------------ */

#include "ie1d.h"
#include "basic.h"
#include "cluster.h"
#include "gaussquad.h"

/* ------------------------------------------------------------
 * Constructor and destructor
 * ------------------------------------------------------------ */

pie1d
new_ie1d(uint n, uint q)
{
  pie1d ie;

  /* Allocate storage */
  ie = (pie1d) allocmem(sizeof(ie1d));

  /* Initialize entries */
  ie->n = n;
  ie->mmax = 0;
  ie->xi = 0;
  ie->m = 0;

  /* Initialize callback functions */
  ie->nearfield = nearfield_ie1d;
  ie->fillrk = 0;
  ie->fillV = 0;
  ie->fillE = 0;
  ie->fillS = 0;

  /* Initialize Gauss quadrature */
  ie->q = q;
  ie->xq = allocreal(q);
  ie->wq = allocreal(q);
  assemble_gauss(q, ie->xq, ie->wq);

  return ie;
}

void
del_ie1d(pie1d ie)
{
  uint m;

  /* Release quadrature points */
  freemem(ie->wq);
  freemem(ie->xq);

  /* Release interpolation points */
  if(ie->xi) {
    for(m=1; m<=ie->mmax; m++)
      freemem(ie->xi[m]);
    freemem(ie->xi);

    /* Make sure the released storage cannot be used by mistake */
    ie->xi = 0;
  }

  /* Release object */
  freemem(ie);
}

/* ------------------------------------------------------------
 * Construct a cluster tree
 * ------------------------------------------------------------ */

static pcluster
buildcluster(uint size, uint off, uint *idx, real h, uint leafsize)
{
  pcluster t;
  uint size0, size1;

  t = 0;
  if(size <= leafsize)
    t = new_cluster(size, idx+off, 0, 1);
  else {
    t = new_cluster(size, idx+off, 2, 1);

    size0 = size / 2;
    t->son[0] = buildcluster(size0, off, idx, h, leafsize);

    size1 = size - size0;
    t->son[1] = buildcluster(size1, off+size0, idx, h, leafsize);
  }

  t->bmin[0] = off * h;
  t->bmax[0] = (off+size) * h;

  update_cluster(t);

  return t;
}

pcluster
build_ie1d_cluster(pcie1d ie, uint leafsize)
{
  pcluster root;
  uint n = ie->n;
  uint *idx;
  uint i;

  idx = (uint *) allocmem(sizeof(uint) * n);
  for(i=0; i<n; i++)
    idx[i] = i;

  root = buildcluster(n, 0, idx, 1.0/n, leafsize);

  return root;
}

/* ------------------------------------------------------------
 * Determine approximation orders
 * ------------------------------------------------------------ */

static real
maxorder(real m0, real m1, pccluster t)
{
  real mt;
  uint i;

  mt = m0;

  if(t->sons > 0) {
    mt = 0.0;
    for(i=0; i<t->sons; i++)
      mt = REAL_MAX(mt, maxorder(m0, m1, t->son[i]));

    mt += m1;
  }

  return mt;
}

static void
setorder(real mt, real m0, real m1, uint tname, pccluster t, uint *m)
{
  uint tname1;
  uint i;

  m[tname] = (uint) mt;

  tname1 = tname+1;
  for(i=0; i<t->sons; i++) {
    setorder(mt-m1, m0, m1, tname1, t->son[i], m);
    
    tname1 += t->son[i]->desc;
  }
  assert(tname1 == tname+t->desc);
}

void
prepare_orders_ie1d(pie1d ie, pccluster root, real m0, real m1)
{
  real mroot;
  uint mmax;
  uint i, m;

  assert(ie->xi == 0);

  mroot = maxorder(m0, m1, root);

  ie->mmax = mmax = (uint) mroot;
  ie->xi = (preal *) allocmem(sizeof(preal) * (mmax+1));

  ie->xi[0] = 0;
  for(m=1; m<=mmax; m++) {
    ie->xi[m] = allocreal(m);
    for(i=0; i<m; i++)
      ie->xi[m][i] = REAL_COS(M_PI * (m-i-0.5) / m);
  }

  ie->m = (uint *) allocmem(sizeof(uint) * root->desc);
  setorder(mroot, m0, m1, 0, root, ie->m);
}

/* ------------------------------------------------------------
 * Nearfield integration
 * ------------------------------------------------------------ */

static real
antiderivative(real x, real y)
{
  if(x == y)
    return 0.0;
  else
    return 0.25 * REAL_SQR(x-y) * (3.0 - 2.0 * REAL_LOG(REAL_ABS(x-y)));
}

void
nearfield_ie1d(const uint *ridx, const uint *cidx, pcie1d ie,
	       pamatrix N)
{
  uint rows = N->rows;
  uint cols = N->cols;
  real h = 1.0 / ie->n;
  uint i, j, ii, jj;
  real x0, x1, y0, y1, val;

  for(j=0; j<cols; j++) {
    jj = (cidx ? cidx[j] : j);
    y0 = jj * h;
    y1 = (jj+1) * h;

    for(i=0; i<rows; i++) {
      ii = (ridx ? ridx[i] : i);
      x0 = ii * h;
      x1 = (ii+1) * h;

      val = antiderivative(x1, y0) + antiderivative(x0, y1) - 2.0 * antiderivative(x0, y0);

      setentry_amatrix(N, i, j, val);
    }
  }
}

/* ------------------------------------------------------------
 * Approximation by Taylor expansion
 * ------------------------------------------------------------ */

void
fillrk_taylor_ie1d(pccluster rc, uint rname,
		   pccluster cc, uint cname,
		   pcie1d ie, prkmatrix r)
{
  uint rows = getrows_rkmatrix(r);
  uint cols = getcols_rkmatrix(r);
  pamatrix A, B;
  field val;
  real xt, x0, x1;
  real yt, y0, y1;
  real lnu, rnu;
  real h = 1.0 / ie->n;
  uint i, ii, j, jj, m, nu;
  
  assert(rc->size == rows);
  assert(cc->size == cols);

  if(getdiam_2_cluster(rc) <= getdiam_2_cluster(cc)) {
    /* Determine order */
    m = ie->m[rname];
    setrank_rkmatrix(r, m);

    /* Compute center of expansion */
    xt = 0.5 * (rc->bmax[0] + rc->bmin[0]);

    /* Fill A */
    A = getA_rkmatrix(r);

    for(i=0; i<rows; i++) {
      /* Get row index */
      ii = rc->idx[i];

      /* Get support of row basis function */
      x0 = ii*h;
      x1 = (ii+1)*h;

      /* Compute values of first Taylor monomial's antiderivative */
      lnu = x0 - xt;
      rnu = x1 - xt;

      for(nu=0; nu<m; nu++) {
	/* Compute integral of Taylor monomial */
	val = rnu - lnu;

	/* Store coefficient */
	setentry_amatrix(A, i, nu, val);

	/* Update Taylor monomial */
	lnu *= (x0 - xt) / (nu + 2.0);
	rnu *= (x1 - xt) / (nu + 2.0);
      }
    }

    /* Fill B */
    B = getB_rkmatrix(r);

    for(j=0; j<cols; j++) {
      /* Get column index */
      jj = cc->idx[j];

      /* Get support of column basis function */
      y0 = jj*h;
      y1 = (jj+1)*h;

      /* First column, antiderivative of log |y-xt| */
      if(m > 0) {
	lnu = (xt - y0) * REAL_LOG(REAL_ABS(y0 - xt)) - (xt - y0);
	rnu = (xt - y1) * REAL_LOG(REAL_ABS(y1 - xt)) - (xt - y1);

	val = rnu - lnu;

	setentry_amatrix(B, j, 0, CONJ(val));
      }

      /* Second column, antiderivative of 1/(y-xt) */
      if(m > 1) {
	lnu = REAL_LOG(REAL_ABS(xt - y0));
	rnu = REAL_LOG(REAL_ABS(xt - y1));

	val = rnu - lnu;

	setentry_amatrix(B, j, 1, CONJ(val));
      }

      /* Remaining columns, antiderivatives of (-1)^(nu-1) (nu-1)!/(y-xt)^nu */
      lnu = 1.0 / (xt - y0);
      rnu = 1.0 / (xt - y1);

      for(nu=2; nu<m; nu++) {
	val = rnu - lnu;

	setentry_amatrix(B, j, nu, CONJ(val));

	lnu *= -(nu-1.0) / (xt - y0);
	rnu *= -(nu-1.0) / (xt - y1);
      }
    }
  }
  else {
    /* Determine order */
    m = ie->m[cname];
    setrank_rkmatrix(r, m);

    /* Compute center of expansion */
    yt = 0.5 * (cc->bmax[0] + cc->bmin[0]);

    /* Fill A */
    A = getA_rkmatrix(r);

    for(i=0; i<rows; i++) {
      /* Get row index */
      ii = rc->idx[i];

      /* Get support of row basis function */
      x0 = ii*h;
      x1 = (ii+1)*h;

      /* First column, antiderivative of -log |x-yt| */
      if(m > 0) {
	lnu = (yt - x0) * REAL_LOG(REAL_ABS(yt - x0)) - (yt - x0);
	rnu = (yt - x1) * REAL_LOG(REAL_ABS(yt - x1)) - (yt - x1);

	val = rnu - lnu;

	setentry_amatrix(A, i, 0, val);
      }

      /* Second column, antiderivative of 1/(x-yt) */
      if(m > 1) {
	lnu = REAL_LOG(REAL_ABS(x0 - yt));
	rnu = REAL_LOG(REAL_ABS(x1 - yt));

	val = rnu - lnu;

	setentry_amatrix(A, i, 1, val);
      }

      /* Remaining columns, antiderivatives of (nu-1)!/(x-yt)^nu */
      lnu = -1.0 / (x0 - yt);
      rnu = -1.0 / (x1 - yt);

      for(nu=2; nu<m; nu++) {
	val = rnu - lnu;

	setentry_amatrix(A, i, nu, val);

	lnu *= (nu-1.0) / (x0 - yt);
	rnu *= (nu-1.0) / (x1 - yt);
      }
    }

    /* Fill B */
    B = getB_rkmatrix(r);

    for(j=0; j<cols; j++) {
      /* Get column index */
      jj = cc->idx[j];

      /* Get support of column basis function */
      y0 = jj*h;
      y1 = (jj+1)*h;

      /* Compute values of first Taylor monomial's antiderivative */
      lnu = y0 - yt;
      rnu = y1 - yt;

      for(nu=0; nu<m; nu++) {
	/* Compute integral of Taylor monomial */
	val = rnu - lnu;

	/* Store coefficient */
	setentry_amatrix(B, j, nu, CONJ(val));

	/* Update Taylor monomial */
	lnu *= (y0 - yt) / (nu + 2.0);
	rnu *= (y1 - yt) / (nu + 2.0);
      }
    }
  }
}

void
fillV_taylor_ie1d(pccluster t, uint tname,
		  pcie1d ie, pamatrix V)
{
  uint size = t->size;
  field val;
  real xt, x0, x1;
  real lnu, rnu;
  real h = 1.0 / ie->n;
  uint i, ii, m, nu;
  
  /* Determine order */
  m = ie->m[tname];
  
  assert(V->rows == t->size);
  assert(V->cols == m);

  /* Compute center of expansion */
  xt = 0.5 * (t->bmax[0] + t->bmin[0]);

  for(i=0; i<size; i++) {
    /* Get row index */
    ii = t->idx[i];

    /* Get support of row basis function */
    x0 = ii*h;
    x1 = (ii+1)*h;
    
    /* Compute values of first Taylor monomial's antiderivative */
    lnu = x0 - xt;
    rnu = x1 - xt;
    
    for(nu=0; nu<m; nu++) {
      /* Compute integral of Taylor monomial */
      val = rnu - lnu;
      
      /* Store coefficient */
      setentry_amatrix(V, i, nu, val);
      
      /* Update Taylor monomial */
      lnu *= (x0 - xt) / (nu + 2.0);
      rnu *= (x1 - xt) / (nu + 2.0);
    }
  }
}

void
fillE_taylor_ie1d(pccluster s, uint sname,
		  pccluster f, uint fname,
		  pcie1d ie, pamatrix E)
{
  field val;
  real xs, xf;
  uint sm, fm;
  uint nu, mu;

  /* Determine son's and father's order */
  sm = ie->m[sname];
  fm = ie->m[fname];

  assert(E->rows == sm);
  assert(E->cols == fm);

  /* Compute centers of expansion */
  xs = 0.5 * (s->bmax[0] + s->bmin[0]);
  xf = 0.5 * (f->bmax[0] + f->bmin[0]);

  /* Derivatives vanish if nu > mu */
  for(mu=0; mu<fm; mu++)
    for(nu=mu+1; nu<sm; nu++)
      setentry_amatrix(E, nu, mu, 0.0);

  /* Evaluate derivatives of father's Taylor monomial in son's center */
  for(mu=0; mu<fm; mu++) {
    /* Higher derivatives are ignored if sm < fm */
    val = 1.0;
    for(nu=mu; nu>=sm; nu--)
      val *= (xs - xf) / (mu - nu + 1.0);

    for(; nu>0; nu--) {
      setentry_amatrix(E, nu, mu, val);

      val *= (xs - xf) / (mu - nu + 1.0);
    }
    setentry_amatrix(E, nu, mu, val);
  }
}

void
fillS_taylor_ie1d(pccluster rc, uint rname,
		  pccluster cc, uint cname,
		  pcie1d ie, pamatrix S)
{
  field val, sval;
  real xr, xc;
  uint rm, cm, m;
  uint nu, mu;

  /* Determine row and column cluster's order */
  rm = ie->m[rname];
  cm = ie->m[cname];
  m = UINT_MIN(rm, cm);

  assert(S->rows == rm);
  assert(S->cols == cm);

  /* Compute centers of expansion */
  xr = 0.5 * (rc->bmax[0] + rc->bmin[0]);
  xc = 0.5 * (cc->bmax[0] + cc->bmin[0]);

  /* Derivatives of order beyond m are ignored */
  for(mu=0; mu<m; mu++)
    for(nu=m-mu; nu<rm; nu++)
      setentry_amatrix(S, nu, mu, 0.0);
  for(; mu<cm; mu++)
    for(nu=0; nu<rm; nu++)
      setentry_amatrix(S, nu, mu, 0.0);

  /* Evaluate kernel function in centers of expansion */
  val = -REAL_LOG(REAL_ABS(xr - xc));
  setentry_amatrix(S, 0, 0, val);

  /* Evaluate derivatives */
  val = -1.0 / (xr - xc);
  for(nu=1; nu<m; nu++) {
    sval = val;

    for(mu=0; mu<=nu; mu++) {
      setentry_amatrix(S, nu-mu, mu, sval);

      sval *= -1.0;
    }

    val *= -(nu / (xr - xc));
  }
}

void
setup_aprx_taylor_ie1d(pie1d ie, pccluster root, real m0, real m1)
{
  ie->nearfield = nearfield_ie1d;
  ie->fillrk = fillrk_taylor_ie1d;
  ie->fillV = fillV_taylor_ie1d;
  ie->fillE = fillE_taylor_ie1d;
  ie->fillS = fillS_taylor_ie1d;

  prepare_orders_ie1d(ie, root, m0, m1);
}

/* ------------------------------------------------------------
 * Approximation by interpolation
 * ------------------------------------------------------------ */

static field
eval_kernel(real x, real y)
{
  return -REAL_LOG(REAL_ABS(x-y));
}

static real
eval_lagrange(uint nu, uint m, pcreal xi, real x)
{
  real val;
  uint i;

  val = 1.0;

  for(i=0; i<nu; i++)
    val *= (x-xi[i]) / (xi[nu]-xi[i]);

  for(i=nu+1; i<m; i++)
    val *= (x-xi[i]) / (xi[nu]-xi[i]);

  return val;
}

void
fillrk_interpolation_ie1d(pccluster rc, uint rname,
			  pccluster cc, uint cname,
			  pcie1d ie, prkmatrix r)
{
  uint rows = getrows_rkmatrix(r);
  uint cols = getcols_rkmatrix(r);
  pamatrix A, B;
  preal xi;
  field val;
  real x, x0, x1;
  real y, y0, y1;
  real mid, rad;
  real h = 1.0 / ie->n;
  uint i, ii, j, jj, k, m, nu;
  
  assert(rc->size == rows);
  assert(cc->size == cols);

  if(getdiam_2_cluster(rc) <= getdiam_2_cluster(cc)) {
    /* Determine order */
    m = ie->m[rname];
    setrank_rkmatrix(r, m);

    /* Prepare transformed interpolation points */
    xi = allocreal(m);
    mid = (rc->bmax[0] + rc->bmin[0]) * 0.5;
    rad = (rc->bmax[0] - rc->bmin[0]) * 0.5;
    for(nu=0; nu<m; nu++)
      xi[nu] = mid + ie->xi[m][nu] * rad;

    /* Fill A */
    A = getA_rkmatrix(r);

    for(i=0; i<rows; i++) {
      /* Get row index */
      ii = rc->idx[i];

      /* Get support of row basis function */
      x0 = ii*h;
      x1 = (ii+1)*h;
      mid = (x1 + x0) * 0.5;
      rad = (x1 - x0) * 0.5;

      for(nu=0; nu<m; nu++) {
	/* Evaluate integral by quadrature */
	val = 0.0;
	for(k=0; k<ie->q; k++) {
	  x = mid + ie->xq[k] * rad;
	  val += ie->wq[k] * eval_lagrange(nu, m, xi, x);
	}
	val *= rad;

	/* Store result */
	setentry_amatrix(A, i, nu, val);
      }
    }

    /* Fill B */
    B = getB_rkmatrix(r);

    for(j=0; j<cols; j++) {
      /* Get column index */
      jj = cc->idx[j];

      /* Get support of column basis function */
      y0 = jj*h;
      y1 = (jj+1)*h;
      mid = (y1 + y0) * 0.5;
      rad = (y1 - y0) * 0.5;

      for(nu=0; nu<m; nu++) {
	/* Evaluate integral by quadrature */
	val = 0.0;
	for(k=0; k<ie->q; k++) {
	  y = mid + ie->xq[k] * rad;
	  val += ie->wq[k] * eval_kernel(xi[nu], y);
	}
	val *= rad;

	/* Store result */
	setentry_amatrix(B, j, nu, CONJ(val));
      }
    }

    /* Clean up */
    freemem(xi);
  }
  else {
    /* Determine order */
    m = ie->m[cname];
    setrank_rkmatrix(r, m);

    /* Prepare transformed interpolation points */
    xi = allocreal(m);
    mid = (cc->bmax[0] + cc->bmin[0]) * 0.5;
    rad = (cc->bmax[0] - cc->bmin[0]) * 0.5;
    for(nu=0; nu<m; nu++)
      xi[nu] = mid + ie->xi[m][nu] * rad;

    /* Fill A */
    A = getA_rkmatrix(r);

    for(i=0; i<rows; i++) {
      /* Get row index */
      ii = rc->idx[i];

      /* Get support of row basis function */
      x0 = ii*h;
      x1 = (ii+1)*h;
      mid = (x1 + x0) * 0.5;
      rad = (x1 - x0) * 0.5;

      for(nu=0; nu<m; nu++) {
	/* Evaluate integral by quadrature */
	val = 0.0;
	for(k=0; k<ie->q; k++) {
	  x = mid + ie->xq[k] * rad;
	  val += ie->wq[k] * eval_kernel(x, xi[nu]);
	}
	val *= rad;

	/* Store result */
	setentry_amatrix(A, i, nu, val);
      }
    }

    /* Fill B */
    B = getB_rkmatrix(r);

    for(j=0; j<cols; j++) {
      /* Get column index */
      jj = cc->idx[j];

      /* Get support of column basis function */
      y0 = jj*h;
      y1 = (jj+1)*h;
      mid = (y1 + y0) * 0.5;
      rad = (y1 - y0) * 0.5;

      for(nu=0; nu<m; nu++) {
	/* Evaluate integral by quadrature */
	val = 0.0;
	for(k=0; k<ie->q; k++) {
	  y = mid + ie->xq[k] * rad;
	  val += ie->wq[k] * eval_lagrange(nu, m, xi, y);
	}
	val *= rad;

	/* Store result */
	setentry_amatrix(B, j, nu, CONJ(val));
      }
    }

    /* Clean up */
    freemem(xi);
  }
}

void
fillV_interpolation_ie1d(pccluster t, uint tname,
			 pcie1d ie, pamatrix V)
{
  uint size = t->size;
  preal xi;
  field val;
  real x, x0, x1;
  real mid, rad;
  real h = 1.0 / ie->n;
  uint i, ii, k, m, nu;

  /* Determine order */
  m = ie->m[tname];

  assert(V->rows == t->size);
  assert(V->cols == m);

  /* Prepare transformed interpolation points */
  xi = allocreal(m);
  mid = (t->bmax[0] + t->bmin[0]) * 0.5;
  rad = (t->bmax[0] - t->bmin[0]) * 0.5;
  for(nu=0; nu<m; nu++)
    xi[nu] = mid + ie->xi[m][nu] * rad;

  for(i=0; i<size; i++) {
    /* Get index */
    ii = t->idx[i];

    /* Get support of row basis function */
    x0 = ii*h;
    x1 = (ii+1)*h;
    mid = (x1 + x0) * 0.5;
    rad = (x1 - x0) * 0.5;

    for(nu=0; nu<m; nu++) {
      /* Evaluate integral by quadrature */
      val = 0.0;
      for(k=0; k<ie->q; k++) {
	x = mid + ie->xq[k] * rad;
	val += ie->wq[k] * eval_lagrange(nu, m, xi, x);
      }
      val *= rad;

      /* Store result */
      setentry_amatrix(V, i, nu, val);
    }
  }
}

void
fillE_interpolation_ie1d(pccluster s, uint sname,
			 pccluster f, uint fname,
			 pcie1d ie, pamatrix E)
{
  preal sxi, fxi;
  field val;
  real mid, rad;
  uint sm, fm;
  uint nu, mu;

  /* Determine son's and father's order */
  sm = ie->m[sname];
  fm = ie->m[fname];

  assert(E->rows == sm);
  assert(E->cols == fm);

  /* Prepare son's transformed interpolation points */
  sxi = allocreal(sm);
  mid = (s->bmax[0] + s->bmin[0]) * 0.5;
  rad = (s->bmax[0] - s->bmin[0]) * 0.5;
  for(nu=0; nu<sm; nu++)
    sxi[nu] = mid + ie->xi[sm][nu] * rad;

  /* Prepare father's transformed interpolation points */
  fxi = allocreal(fm);
  mid = (f->bmax[0] + f->bmin[0]) * 0.5;
  rad = (f->bmax[0] - f->bmin[0]) * 0.5;
  for(mu=0; mu<fm; mu++)
    fxi[mu] = mid + ie->xi[fm][mu] * rad;

  /* Evaluate father's Lagrange polynomial in son's interpolation points */
  for(mu=0; mu<fm; mu++)
    for(nu=0; nu<sm; nu++) {
      val = eval_lagrange(mu, fm, fxi, sxi[nu]);

      setentry_amatrix(E, nu, mu, val);
    }
}

void
fillS_interpolation_ie1d(pccluster rc, uint rname,
			 pccluster cc, uint cname,
			 pcie1d ie, pamatrix S)
{
  preal rxi, cxi;
  field val;
  real mid, rad;
  uint rm, cm;
  uint nu, mu;

  /* Determine row and column cluster's order */
  rm = ie->m[rname];
  cm = ie->m[cname];

  assert(S->rows == rm);
  assert(S->cols == cm);

  /* Prepare row cluster's transformed interpolation points */
  rxi = allocreal(rm);
  mid = (rc->bmax[0] + rc->bmin[0]) * 0.5;
  rad = (rc->bmax[0] - rc->bmin[0]) * 0.5;
  for(nu=0; nu<rm; nu++)
    rxi[nu] = mid + ie->xi[rm][nu] * rad;

  /* Prepare column cluster's transformed interpolation points */
  cxi = allocreal(cm);
  mid = (cc->bmax[0] + cc->bmin[0]) * 0.5;
  rad = (cc->bmax[0] - cc->bmin[0]) * 0.5;
  for(mu=0; mu<cm; mu++)
    cxi[mu] = mid + ie->xi[cm][mu] * rad;

  /* Evaluate kernel function in interpolation points */
  for(mu=0; mu<cm; mu++)
    for(nu=0; nu<rm; nu++) {
      val = eval_kernel(rxi[nu], cxi[mu]);

      setentry_amatrix(S, nu, mu, val);
    }
}

void
setup_aprx_interpolation_ie1d(pie1d ie, pccluster root, real m0, real m1)
{
  ie->nearfield = nearfield_ie1d;
  ie->fillrk = fillrk_interpolation_ie1d;
  ie->fillV = fillV_interpolation_ie1d;
  ie->fillE = fillE_interpolation_ie1d;
  ie->fillS = fillS_interpolation_ie1d;

  prepare_orders_ie1d(ie, root, m0, m1);
}

/* ------------------------------------------------------------
 * Fill a hierarchical matrix
 * ------------------------------------------------------------ */

static void
fill_hmatrix(phmatrix G, uint rname, uint cname, pcie1d ie)
{
  uint rsons, csons;
  uint rname1, cname1;
  uint i, j;

  if(G->son) {
    /* Handle submatrices */
    rsons = G->rsons;
    csons = G->csons;

    /* Compute name of first column son */
    cname1 = (G->son[0]->cc == G->cc ? cname : cname+1);

    for(j=0; j<csons; j++) {
      /* Compute name of first row son */
      rname1 = (G->son[0]->rc == G->rc ? rname : rname+1);

      for(i=0; i<rsons; i++) {
	/* Fill submatrices */
	fill_hmatrix(G->son[i+j*rsons], rname1, cname1, ie);

	/* Switch to next row son */
	rname1 += G->son[i]->rc->desc;
      }
      assert(rname1 == rname+G->rc->desc);

      /* Switch to next column son */
      cname1 += G->son[j*rsons]->cc->desc;
    }
    assert(cname1 == cname+G->cc->desc);
  }
  else if(G->f) {
    /* Fill inadmissible submatrix */
    ie->nearfield(G->rc->idx, G->cc->idx, ie, G->f);
  }
  else if(G->r) {
    /* Fill admissible submatrix */
    ie->fillrk(G->rc, rname, G->cc, cname, ie, G->r);
  }
}

void
fill_hmatrix_ie1d(pcie1d ie, phmatrix G)
{
  fill_hmatrix(G, 0, 0, ie);
}

/* ------------------------------------------------------------
 * Fill a cluster basis
 * ------------------------------------------------------------ */

static void
fill_clusterbasis(pclusterbasis cb, uint tname, pcie1d ie)
{
  uint i, tname1;

  if(cb->son) {
    /* Compute name of first son */
    tname1 = tname+1;
    for(i=0; i<cb->sons; i++) {
      /* Prepare son's cluster basis */
      fill_clusterbasis(cb->son[i], tname1, ie);

      /* Switch to next son */
      tname1 += cb->son[i]->t->desc;
    }
    assert(tname1 == tname+cb->t->desc);

    /* Set the rank */
    setrank_clusterbasis(cb, ie->m[tname]);

    /* Compute name of first son */
    tname1 = tname+1;
    for(i=0; i<cb->sons; i++) {
      /* Fill transfer matrix */
      ie->fillE(cb->son[i]->t, tname1, cb->t, tname,
		ie, getE_clusterbasis(cb->son[i]));

      /* Switch to next son */
      tname1 += cb->son[i]->t->desc;
    }
    assert(tname1 == tname+cb->t->desc);
  }
  else {
    /* Set the rank */
    setrank_clusterbasis(cb, ie->m[tname]);

    /* Fill leaf matrix */
    ie->fillV(cb->t, tname, ie, getV_clusterbasis(cb));
  }
}

void
fill_clusterbasis_ie1d(pcie1d ie, pclusterbasis cb)
{
  fill_clusterbasis(cb, 0, ie);
}

/* ------------------------------------------------------------
 * Fill an H2-matrix
 * ------------------------------------------------------------ */

static void
fill_h2matrix(ph2matrix G, uint rname, uint cname, pcie1d ie)
{
  uint rsons, csons;
  uint rname1, cname1;
  uint i, j;

  if(G->son) {
    /* Handle submatrices */
    rsons = G->rsons;
    csons = G->csons;

    /* Compute name of first column son */
    cname1 = (G->son[0]->cb == G->cb ? cname : cname+1);

    for(j=0; j<csons; j++) {
      /* Compute name of first row son */
      rname1 = (G->son[0]->rb == G->rb ? rname : rname+1);

      for(i=0; i<rsons; i++) {
	/* Fill submatrices */
	fill_h2matrix(G->son[i+j*rsons], rname1, cname1, ie);

	/* Switch to next row son */
	rname1 += G->son[i]->rb->t->desc;
      }
      assert(rname1 == rname+G->rb->t->desc);

      /* Switch to next column son */
      cname1 += G->son[j*rsons]->cb->t->desc;
    }
    assert(cname1 == cname+G->cb->t->desc);
  }
  else if(G->f) {
    /* Fill inadmissible submatrix */
    ie->nearfield(G->rb->t->idx, G->cb->t->idx, ie, G->f);
  }
  else if(G->u) {
    /* Fill admissible submatrix */
    ie->fillS(G->rb->t, rname, G->cb->t, cname, ie,
	      getS_uniform(G->u));
  }
}

void
fill_h2matrix_ie1d(pcie1d ie, ph2matrix G)
{
  fill_h2matrix(G, 0, 0, ie);
}
