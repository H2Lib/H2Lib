/* ------------------------------------------------------------
 This is the file "laplacebem3d.c" of the H2Lib package.
 All rights reserved, Sven Christophersen 2011
 ------------------------------------------------------------ */

/**
 * @file laplacebem3d.c
 * @author Sven Christophersen
 * @date 2011
 */

/* C STD LIBRARY */
/* CORE 0 */
#include "basic.h"
#include "parameters.h"
/* CORE 1 */
/* CORE 2 */
/* CORE 3 */
/* SIMPLE */
/* PARTICLES */
/* BEM */
#include "laplacebem3d.h"

/* This constant comes from the fundamental solution of the Laplace-equation,
 * which is defined as: @f$ g(x,y) := \frac{1}{4 \pi} \frac{1}{ \left\lVert x - y \right\rVert}  @f$ .
 * Therefore the constant takes the value @f$ \frac{1}{4 \pi}@f$ .
 * */
#define KERNEL_CONST_BEM3D 0.0795774715459476679

static void
fill_slp_cc_laplacebem3d(const uint * ridx, const uint * cidx,
			 pcbem3d bem, bool ntrans, pamatrix N)
{
  const pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const real *gr_g = (const real *) gr->g;
  field    *aa = N->a;
  uint      rows = N->rows;
  uint      cols = N->cols;
  longindex ld = N->ld;

  const real *A_t, *B_t, *C_t, *A_s, *B_s, *C_s;
  const uint *tri_t, *tri_s;
  real     *xq, *yq, *wq;
  uint      tp[3], sp[3];
  real      Ax, Bx, Cx, Ay, By, Cy, tx, sx, ty, sy, dx, dy, dz, factor,
    factor2;
  field     sum;
  uint      q, nq, ss, tt, s, t;

  if (ntrans == true) {
    for (t = 0; t < cols; ++t) {
      tt = (ridx == NULL ? t : ridx[t]);
      tri_t = gr_t[tt];
      factor = gr_g[tt] * KERNEL_CONST_BEM3D;
      for (s = 0; s < rows; ++s) {
	ss = (cidx == NULL ? s : cidx[s]);
	tri_s = gr_t[ss];
	factor2 = factor * gr_g[ss];

	select_quadrature_singquad2d(bem->sq, tri_t, tri_s, tp, sp, &xq, &yq,
				     &wq, &nq, &sum);
	wq += 9 * nq;

	A_t = gr_x[tri_t[tp[0]]];
	B_t = gr_x[tri_t[tp[1]]];
	C_t = gr_x[tri_t[tp[2]]];
	A_s = gr_x[tri_s[sp[0]]];
	B_s = gr_x[tri_s[sp[1]]];
	C_s = gr_x[tri_s[sp[2]]];

	for (q = 0; q < nq; ++q) {
	  tx = xq[q];
	  sx = xq[q + nq];
	  ty = yq[q];
	  sy = yq[q + nq];
	  Ax = 1.0 - tx;
	  Bx = tx - sx;
	  Cx = sx;
	  Ay = 1.0 - ty;
	  By = ty - sy;
	  Cy = sy;

	  dx = A_t[0] * Ax + B_t[0] * Bx + C_t[0] * Cx
	    - (A_s[0] * Ay + B_s[0] * By + C_s[0] * Cy);
	  dy = A_t[1] * Ax + B_t[1] * Bx + C_t[1] * Cx
	    - (A_s[1] * Ay + B_s[1] * By + C_s[1] * Cy);
	  dz = A_t[2] * Ax + B_t[2] * Bx + C_t[2] * Cx
	    - (A_s[2] * Ay + B_s[2] * By + C_s[2] * Cy);

	  sum += wq[q] / REAL_SQRT(dx * dx + dy * dy + dz * dz);
	}

	aa[s + t * ld] = sum * factor2;
      }
    }
  }
  else {
    for (s = 0; s < cols; ++s) {
      ss = (cidx == NULL ? s : cidx[s]);
      tri_s = gr_t[ss];
      factor = gr_g[ss] * KERNEL_CONST_BEM3D;

      for (t = 0; t < rows; ++t) {
	tt = (ridx == NULL ? t : ridx[t]);
	tri_t = gr_t[tt];
	factor2 = factor * gr_g[tt];

	select_quadrature_singquad2d(bem->sq, tri_t, tri_s, tp, sp, &xq, &yq,
				     &wq, &nq, &sum);
	wq += 9 * nq;

	A_t = gr_x[tri_t[tp[0]]];
	B_t = gr_x[tri_t[tp[1]]];
	C_t = gr_x[tri_t[tp[2]]];
	A_s = gr_x[tri_s[sp[0]]];
	B_s = gr_x[tri_s[sp[1]]];
	C_s = gr_x[tri_s[sp[2]]];

	for (q = 0; q < nq; ++q) {
	  tx = xq[q];
	  sx = xq[q + nq];
	  ty = yq[q];
	  sy = yq[q + nq];
	  Ax = 1.0 - tx;
	  Bx = tx - sx;
	  Cx = sx;
	  Ay = 1.0 - ty;
	  By = ty - sy;
	  Cy = sy;

	  dx = A_t[0] * Ax + B_t[0] * Bx + C_t[0] * Cx
	    - (A_s[0] * Ay + B_s[0] * By + C_s[0] * Cy);
	  dy = A_t[1] * Ax + B_t[1] * Bx + C_t[1] * Cx
	    - (A_s[1] * Ay + B_s[1] * By + C_s[1] * Cy);
	  dz = A_t[2] * Ax + B_t[2] * Bx + C_t[2] * Cx
	    - (A_s[2] * Ay + B_s[2] * By + C_s[2] * Cy);

	  sum += wq[q] / REAL_SQRT(dx * dx + dy * dy + dz * dz);
	}

	aa[t + s * ld] = sum * factor2;
      }
    }
  }
}

static void
fill_dlp_cc_laplacebem3d(const uint * ridx, const uint * cidx,
			 pcbem3d bem, bool ntrans, pamatrix N)
{
  const pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const real *gr_g = (const real *) gr->g;
  field    *aa = N->a;
  uint      rows = N->rows;
  uint      cols = N->cols;
  longindex ld = N->ld;

  const real *A_t, *B_t, *C_t, *A_s, *B_s, *C_s, *ns;
  const uint *tri_t, *tri_s;
  real     *xq, *yq, *wq;
  uint      tp[3], sp[3];
  real      norm, dx, dy, dz, Ax, Bx, Cx, Ay, By, Cy, tx, sx, ty, sy, factor,
    factor2;
  field     res;
  uint      q, tt, ss, nq, t, s;

  if (ntrans == true) {
    for (t = 0; t < cols; ++t) {
      tt = (ridx == NULL ? t : ridx[t]);
      tri_t = gr_t[tt];
      factor = gr_g[tt] * KERNEL_CONST_BEM3D;

      for (s = 0; s < rows; ++s) {
	ss = (cidx == NULL ? s : cidx[s]);
	tri_s = gr_t[ss];
	ns = gr_n[ss];
	factor2 = factor * gr_g[ss];

	if (tt == ss) {
	  res = 0.5 * bem->alpha * gr_g[tt];

	  aa[t * (ld + 1)] = res;
	}
	else {

	  (void) select_quadrature_singquad2d(bem->sq, tri_t, tri_s, tp, sp,
					      &xq, &yq, &wq, &nq, &res);
	  wq += 9 * nq;

	  A_t = gr_x[tri_t[tp[0]]];
	  B_t = gr_x[tri_t[tp[1]]];
	  C_t = gr_x[tri_t[tp[2]]];
	  A_s = gr_x[tri_s[sp[0]]];
	  B_s = gr_x[tri_s[sp[1]]];
	  C_s = gr_x[tri_s[sp[2]]];

	  for (q = 0; q < nq; ++q) {
	    tx = xq[q];
	    sx = xq[q + nq];
	    ty = yq[q];
	    sy = yq[q + nq];
	    Ax = 1.0 - tx;
	    Bx = tx - sx;
	    Cx = sx;
	    Ay = 1.0 - ty;
	    By = ty - sy;
	    Cy = sy;

	    dx = A_t[0] * Ax + B_t[0] * Bx + C_t[0] * Cx
	      - (A_s[0] * Ay + B_s[0] * By + C_s[0] * Cy);
	    dy = A_t[1] * Ax + B_t[1] * Bx + C_t[1] * Cx
	      - (A_s[1] * Ay + B_s[1] * By + C_s[1] * Cy);
	    dz = A_t[2] * Ax + B_t[2] * Bx + C_t[2] * Cx
	      - (A_s[2] * Ay + B_s[2] * By + C_s[2] * Cy);

	    norm = dx * dx + dy * dy + dz * dz;

	    res += wq[q] * (dx * ns[0] + dy * ns[1] + dz * ns[2])
	      / (norm * REAL_SQRT(norm));
	  }

	  aa[s + t * ld] = res * factor2;
	}
      }
    }
  }
  else {
    for (s = 0; s < cols; ++s) {
      ss = (cidx == NULL ? s : cidx[s]);
      tri_s = gr_t[ss];
      ns = gr_n[ss];
      factor = gr_g[ss] * KERNEL_CONST_BEM3D;

      for (t = 0; t < rows; ++t) {
	tt = (ridx == NULL ? t : ridx[t]);
	tri_t = gr_t[tt];
	factor2 = factor * gr_g[tt];

	if (tt == ss) {
	  res = 0.5 * bem->alpha * gr_g[tt];

	  aa[t * (ld + 1)] = res;
	}
	else {

	  (void) select_quadrature_singquad2d(bem->sq, tri_t, tri_s, tp, sp,
					      &xq, &yq, &wq, &nq, &res);
	  wq += 9 * nq;

	  A_t = gr_x[tri_t[tp[0]]];
	  B_t = gr_x[tri_t[tp[1]]];
	  C_t = gr_x[tri_t[tp[2]]];
	  A_s = gr_x[tri_s[sp[0]]];
	  B_s = gr_x[tri_s[sp[1]]];
	  C_s = gr_x[tri_s[sp[2]]];

	  for (q = 0; q < nq; ++q) {
	    tx = xq[q];
	    sx = xq[q + nq];
	    ty = yq[q];
	    sy = yq[q + nq];
	    Ax = 1.0 - tx;
	    Bx = tx - sx;
	    Cx = sx;
	    Ay = 1.0 - ty;
	    By = ty - sy;
	    Cy = sy;

	    dx = A_t[0] * Ax + B_t[0] * Bx + C_t[0] * Cx
	      - (A_s[0] * Ay + B_s[0] * By + C_s[0] * Cy);
	    dy = A_t[1] * Ax + B_t[1] * Bx + C_t[1] * Cx
	      - (A_s[1] * Ay + B_s[1] * By + C_s[1] * Cy);
	    dz = A_t[2] * Ax + B_t[2] * Bx + C_t[2] * Cx
	      - (A_s[2] * Ay + B_s[2] * By + C_s[2] * Cy);

	    norm = dx * dx + dy * dy + dz * dz;

	    res += wq[q] * (dx * ns[0] + dy * ns[1] + dz * ns[2])
	      / (norm * REAL_SQRT(norm));
	  }

	  aa[t + s * ld] = res * factor2;
	}
      }
    }
  }
}

static void
fill_dlp_cl_laplacebem3d(const uint * ridx, const uint * cidx,
			 pcbem3d bem, bool ntrans, pamatrix N)
{
  const pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const preal gr_g = (const preal) gr->g;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const uint triangles = gr->triangles;
  plistnode *v2t = bem->v2t;
  field    *quad;
  field    *aa = N->a;
  uint      rows = N->rows;
  uint      cols = N->cols;
  longindex ld = N->ld;

  ptri_list tl, tl1;
  pvert_list vl;
  const real *A_t, *B_t, *C_t, *A_s, *B_s, *C_s, *ns;
  const uint *tri_t, *tri_s;
  plistnode v;
  real     *xq, *xq2, *yq, *yq2, *wq, *mass;
  uint      tp[3], sp[3], tri_tp[3], tri_sp[3];
  real      norm, Ax, Bx, Cx, Ay, By, Cy, tx, sx, ty, sy, dx, dy, dz, factor,
    factor2;
  field     res, base;
  uint      i, j, t, s, q, nq, cj;
  uint      ii, jj, tt, ss, vv;

  clear_amatrix(N);

  quad = allocfield(bem->sq->nmax);

  if (ntrans == true) {

    tl = NULL;

    cj = 0;
    for (i = 0; i < rows; ++i) {
      ii = (cidx == NULL ? i : cidx[i]);
      for (v = v2t[ii], vv = v->data; v->next != NULL;
	   v = v->next, vv = v->data) {

	tl1 = tl;
	while (tl1 && tl1->t != vv) {
	  tl1 = tl1->next;
	}

	if (tl1 == NULL) {
	  tl1 = tl = new_tri_list(tl);
	  tl->t = vv;
	  cj++;
	}

	tl1->vl = new_vert_list(tl1->vl);
	tl1->vl->v = i;
      }
    }

    for (s = 0, tl1 = tl; s < cj; s++, tl1 = tl1->next) {
      ss = tl1->t;
      assert(ss < triangles);
      factor = gr_g[ss] * KERNEL_CONST_BEM3D;
      tri_s = gr_t[ss];
      ns = gr_n[ss];
      for (t = 0; t < cols; t++) {
	tt = (ridx == NULL ? t : ridx[t]);
	assert(tt < triangles);
	tri_t = gr_t[tt];

	if (tt != ss) {
	  factor2 = factor * gr_g[tt];

	  select_quadrature_singquad2d(bem->sq, tri_t, tri_s, tp, sp, &xq,
				       &yq, &wq, &nq, &base);
	  xq2 = xq + nq;
	  yq2 = yq + nq;

	  for (i = 0; i < 3; ++i) {
	    tri_tp[i] = tri_t[tp[i]];
	    tri_sp[i] = tri_s[sp[i]];
	  }

	  A_t = gr_x[tri_tp[0]];
	  B_t = gr_x[tri_tp[1]];
	  C_t = gr_x[tri_tp[2]];
	  A_s = gr_x[tri_sp[0]];
	  B_s = gr_x[tri_sp[1]];
	  C_s = gr_x[tri_sp[2]];

	  for (q = 0; q < nq; ++q) {
	    tx = xq[q];
	    sx = xq2[q];
	    ty = yq[q];
	    sy = yq2[q];
	    Ax = 1.0 - tx;
	    Bx = tx - sx;
	    Cx = sx;
	    Ay = 1.0 - ty;
	    By = ty - sy;
	    Cy = sy;

	    dx = A_t[0] * Ax + B_t[0] * Bx + C_t[0] * Cx
	      - (A_s[0] * Ay + B_s[0] * By + C_s[0] * Cy);
	    dy = A_t[1] * Ax + B_t[1] * Bx + C_t[1] * Cx
	      - (A_s[1] * Ay + B_s[1] * By + C_s[1] * Cy);
	    dz = A_t[2] * Ax + B_t[2] * Bx + C_t[2] * Cx
	      - (A_s[2] * Ay + B_s[2] * By + C_s[2] * Cy);

	    norm = dx * dx + dy * dy + dz * dz;

	    quad[q] = (dx * ns[0] + dy * ns[1] + dz * ns[2])
	      / (norm * REAL_SQRT(norm));
	  }

	  vl = tl1->vl;
	  while (vl) {
	    j = vl->v;
	    if (j < rows) {
	      jj = cidx == NULL ? j : cidx[j];
	      for (i = 0; i < 3; ++i) {
		if (jj == tri_sp[i]) {
		  res = base;

		  for (q = 0; q < nq; ++q) {
		    res += wq[q] * quad[q];
		  }

		  aa[j + t * ld] += res * factor2;
		}
		wq += nq;
	      }
	      wq -= 3 * nq;
	    }
	    vl = vl->next;
	  }

	}
	else {

	  for (i = 0; i < 3; ++i) {
	    tri_sp[i] = tri_s[i];
	  }

	  mass = bem->mass;
	  factor2 = bem->alpha * gr_g[tt];

	  if (cidx != NULL) {

	    vl = tl1->vl;
	    while (vl) {
	      j = vl->v;
	      if (j < rows) {
		jj = cidx == NULL ? j : cidx[j];
		for (i = 0; i < 3; ++i) {
		  if (jj == tri_sp[i]) {
		    aa[j + t * ld] += factor2 * *mass;
		  }
		  mass++;
		}
		mass = bem->mass;
	      }
	      vl = vl->next;
	    }
	  }
	}
      }
    }
  }
  else {

    tl = NULL;

    cj = 0;
    for (i = 0; i < cols; ++i) {
      ii = (cidx == NULL ? i : cidx[i]);
      for (v = v2t[ii], vv = v->data; v->next != NULL;
	   v = v->next, vv = v->data) {

	tl1 = tl;
	while (tl1 && tl1->t != vv) {
	  tl1 = tl1->next;
	}

	if (tl1 == NULL) {
	  tl1 = tl = new_tri_list(tl);
	  tl->t = vv;
	  cj++;
	}

	tl1->vl = new_vert_list(tl1->vl);
	tl1->vl->v = i;
      }
    }

    for (s = 0, tl1 = tl; s < cj; s++, tl1 = tl1->next) {
      ss = tl1->t;
      assert(ss < triangles);
      factor = gr_g[ss] * KERNEL_CONST_BEM3D;
      tri_s = gr_t[ss];
      ns = gr_n[ss];
      for (t = 0; t < rows; t++) {
	tt = (ridx == NULL ? t : ridx[t]);
	assert(tt < triangles);
	tri_t = gr_t[tt];

	if (tt != ss) {
	  factor2 = factor * gr_g[tt];

	  select_quadrature_singquad2d(bem->sq, tri_t, tri_s, tp, sp, &xq,
				       &yq, &wq, &nq, &base);

	  for (i = 0; i < 3; ++i) {
	    tri_tp[i] = tri_t[tp[i]];
	    tri_sp[i] = tri_s[sp[i]];
	  }

	  A_t = gr_x[tri_tp[0]];
	  B_t = gr_x[tri_tp[1]];
	  C_t = gr_x[tri_tp[2]];
	  A_s = gr_x[tri_sp[0]];
	  B_s = gr_x[tri_sp[1]];
	  C_s = gr_x[tri_sp[2]];

	  for (q = 0; q < nq; ++q) {
	    tx = xq[q];
	    sx = xq[q + nq];
	    ty = yq[q];
	    sy = yq[q + nq];
	    Ax = 1.0 - tx;
	    Bx = tx - sx;
	    Cx = sx;
	    Ay = 1.0 - ty;
	    By = ty - sy;
	    Cy = sy;

	    dx = A_t[0] * Ax + B_t[0] * Bx + C_t[0] * Cx
	      - (A_s[0] * Ay + B_s[0] * By + C_s[0] * Cy);
	    dy = A_t[1] * Ax + B_t[1] * Bx + C_t[1] * Cx
	      - (A_s[1] * Ay + B_s[1] * By + C_s[1] * Cy);
	    dz = A_t[2] * Ax + B_t[2] * Bx + C_t[2] * Cx
	      - (A_s[2] * Ay + B_s[2] * By + C_s[2] * Cy);

	    norm = dx * dx + dy * dy + dz * dz;

	    quad[q] = (dx * ns[0] + dy * ns[1] + dz * ns[2])
	      / (norm * REAL_SQRT(norm));
	  }

	  vl = tl1->vl;
	  while (vl) {
	    j = vl->v;
	    if (j < cols) {
	      jj = cidx == NULL ? j : cidx[j];
	      for (i = 0; i < 3; ++i) {
		if (jj == tri_sp[i]) {
		  res = base;

		  for (q = 0; q < nq; ++q) {
		    res += wq[q] * quad[q];
		  }

		  aa[t + j * ld] += res * factor2;
		}
		wq += nq;
	      }
	      wq -= 3 * nq;
	    }
	    vl = vl->next;
	  }

	}
	else {

	  for (i = 0; i < 3; ++i) {
	    tri_sp[i] = tri_s[i];
	  }

	  mass = bem->mass;
	  factor2 = bem->alpha * gr_g[tt];

	  vl = tl1->vl;
	  while (vl) {
	    j = vl->v;
	    if (j < cols) {
	      jj = cidx == NULL ? j : cidx[j];
	      for (i = 0; i < 3; ++i) {
		if (jj == tri_sp[i]) {
		  aa[t + j * ld] += factor2 * *mass;
		}
		mass++;
	      }
	      mass = bem->mass;
	    }
	    vl = vl->next;
	  }
	}
      }
    }
  }

  del_tri_list(tl);
  freemem(quad);
}

static void
fill_slp_ll_laplacebem3d(const uint * ridx, const uint * cidx,
			 pcbem3d bem, bool ntrans, pamatrix N)
{
  const pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const preal gr_g = (const preal) gr->g;
  const uint triangles = gr->triangles;
  plistnode *v2t = bem->v2t;
  field    *quad;
  field    *aa = N->a;
  uint      rows = N->rows;
  uint      cols = N->cols;
  longindex ld = N->ld;

  ptri_list tl_r, tl1_r, tl_c, tl1_c;
  pvert_list vl_r, vl_c;
  const real *A_t, *B_t, *C_t, *A_s, *B_s, *C_s;
  const uint *tri_t, *tri_s;
  plistnode v;
  real     *xq, *yq, *wq, *ww;
  uint      tp[3], sp[3], tri_tp[3], tri_sp[3];
  real      norm, Ax, Bx, Cx, Ay, By, Cy, tx, sx, ty, sy, dx, dy, dz, factor,
    factor2, res, base;
  uint      i, j, t, s, k, l, rj, cj, tt, ss, q, nq, ii, jj, vv;

  quad = allocfield(bem->sq->nmax);

  clear_amatrix(N);

  if (ntrans == true) {

    tl_r = NULL;
    tl_c = NULL;

    rj = 0;
    for (i = 0; i < cols; ++i) {
      ii = (ridx == NULL ? i : ridx[i]);
      for (v = v2t[ii], vv = v->data; v->next != NULL;
	   v = v->next, vv = v->data) {

	tl1_r = tl_r;
	while (tl1_r && tl1_r->t != vv) {
	  tl1_r = tl1_r->next;
	}

	if (tl1_r == NULL) {
	  tl1_r = tl_r = new_tri_list(tl_r);
	  tl_r->t = vv;
	  rj++;
	}

	tl1_r->vl = new_vert_list(tl1_r->vl);
	tl1_r->vl->v = i;
      }
    }

    cj = 0;
    for (i = 0; i < rows; ++i) {
      ii = (cidx == NULL ? i : cidx[i]);
      for (v = v2t[ii], vv = v->data; v->next != NULL;
	   v = v->next, vv = v->data) {

	tl1_c = tl_c;
	while (tl1_c && tl1_c->t != vv) {
	  tl1_c = tl1_c->next;
	}

	if (tl1_c == NULL) {
	  tl1_c = tl_c = new_tri_list(tl_c);
	  tl_c->t = vv;
	  cj++;
	}

	tl1_c->vl = new_vert_list(tl1_c->vl);
	tl1_c->vl->v = i;
      }
    }

    for (s = 0, tl1_c = tl_c; s < cj; s++, tl1_c = tl1_c->next) {
      ss = tl1_c->t;
      assert(ss < triangles);
      factor = gr_g[ss] * KERNEL_CONST_BEM3D;
      tri_s = gr_t[ss];
      for (t = 0, tl1_r = tl_r; t < rj; t++, tl1_r = tl1_r->next) {
	tt = tl1_r->t;
	assert(tt < triangles);
	factor2 = factor * gr_g[tt];
	tri_t = gr_t[tt];

	select_quadrature_singquad2d(bem->sq, tri_t, tri_s, tp, sp, &xq, &yq,
				     &wq, &nq, &base);

	for (i = 0; i < 3; ++i) {
	  tri_tp[i] = tri_t[tp[i]];
	  tri_sp[i] = tri_s[sp[i]];
	}

	A_t = gr_x[tri_tp[0]];
	B_t = gr_x[tri_tp[1]];
	C_t = gr_x[tri_tp[2]];
	A_s = gr_x[tri_sp[0]];
	B_s = gr_x[tri_sp[1]];
	C_s = gr_x[tri_sp[2]];

	for (q = 0; q < nq; ++q) {
	  tx = xq[q];
	  sx = xq[q + nq];
	  ty = yq[q];
	  sy = yq[q + nq];
	  Ax = 1.0 - tx;
	  Bx = tx - sx;
	  Cx = sx;
	  Ay = 1.0 - ty;
	  By = ty - sy;
	  Cy = sy;

	  dx = A_t[0] * Ax + B_t[0] * Bx + C_t[0] * Cx
	    - (A_s[0] * Ay + B_s[0] * By + C_s[0] * Cy);
	  dy = A_t[1] * Ax + B_t[1] * Bx + C_t[1] * Cx
	    - (A_s[1] * Ay + B_s[1] * By + C_s[1] * Cy);
	  dz = A_t[2] * Ax + B_t[2] * Bx + C_t[2] * Cx
	    - (A_s[2] * Ay + B_s[2] * By + C_s[2] * Cy);

	  norm = REAL_SQRT(dx * dx + dy * dy + dz * dz);

	  quad[q] = 1.0 / norm;
	}

	vl_c = tl1_c->vl;
	while (vl_c) {
	  j = vl_c->v;
	  assert(j < rows);
	  jj = ((cidx == NULL) ? j : cidx[j]);
	  for (k = 0; k < 3; ++k) {
	    if (jj == tri_sp[k]) {
	      vl_r = tl1_r->vl;
	      while (vl_r) {
		i = vl_r->v;
		assert(i < cols);
		ii = ((ridx == NULL) ? i : ridx[i]);
		for (l = 0; l < 3; ++l) {
		  if (ii == tri_tp[l]) {
		    res = base;

		    ww = wq + (l + k * 3) * nq;
		    for (q = 0; q < nq; ++q) {
		      res += ww[q] * quad[q];
		    }

		    aa[j + i * ld] += res * factor2;
		  }
		}
		vl_r = vl_r->next;
	      }
	    }
	  }
	  vl_c = vl_c->next;
	}
      }
    }
  }
  else {

    tl_r = NULL;
    tl_c = NULL;

    rj = 0;
    for (i = 0; i < rows; ++i) {
      ii = (ridx == NULL ? i : ridx[i]);
      for (v = v2t[ii], vv = v->data; v->next != NULL;
	   v = v->next, vv = v->data) {

	tl1_r = tl_r;
	while (tl1_r && tl1_r->t != vv) {
	  tl1_r = tl1_r->next;
	}

	if (tl1_r == NULL) {
	  tl1_r = tl_r = new_tri_list(tl_r);
	  tl_r->t = vv;
	  rj++;
	}

	tl1_r->vl = new_vert_list(tl1_r->vl);
	tl1_r->vl->v = i;
      }
    }

    cj = 0;
    for (i = 0; i < cols; ++i) {
      ii = (cidx == NULL ? i : cidx[i]);
      for (v = v2t[ii], vv = v->data; v->next != NULL;
	   v = v->next, vv = v->data) {

	tl1_c = tl_c;
	while (tl1_c && tl1_c->t != vv) {
	  tl1_c = tl1_c->next;
	}

	if (tl1_c == NULL) {
	  tl1_c = tl_c = new_tri_list(tl_c);
	  tl_c->t = vv;
	  cj++;
	}

	tl1_c->vl = new_vert_list(tl1_c->vl);
	tl1_c->vl->v = i;
      }
    }

    for (s = 0, tl1_c = tl_c; s < cj; s++, tl1_c = tl1_c->next) {
      ss = tl1_c->t;
      assert(ss < triangles);
      factor = gr_g[ss] * KERNEL_CONST_BEM3D;
      tri_s = gr_t[ss];
      for (t = 0, tl1_r = tl_r; t < rj; t++, tl1_r = tl1_r->next) {
	tt = tl1_r->t;
	assert(tt < triangles);
	factor2 = factor * gr_g[tt];
	tri_t = gr_t[tt];

	select_quadrature_singquad2d(bem->sq, tri_t, tri_s, tp, sp, &xq, &yq,
				     &wq, &nq, &base);

	for (i = 0; i < 3; ++i) {
	  tri_tp[i] = tri_t[tp[i]];
	  tri_sp[i] = tri_s[sp[i]];
	}

	A_t = gr_x[tri_tp[0]];
	B_t = gr_x[tri_tp[1]];
	C_t = gr_x[tri_tp[2]];
	A_s = gr_x[tri_sp[0]];
	B_s = gr_x[tri_sp[1]];
	C_s = gr_x[tri_sp[2]];

	for (q = 0; q < nq; ++q) {
	  tx = xq[q];
	  sx = xq[q + nq];
	  ty = yq[q];
	  sy = yq[q + nq];
	  Ax = 1.0 - tx;
	  Bx = tx - sx;
	  Cx = sx;
	  Ay = 1.0 - ty;
	  By = ty - sy;
	  Cy = sy;

	  dx = A_t[0] * Ax + B_t[0] * Bx + C_t[0] * Cx
	    - (A_s[0] * Ay + B_s[0] * By + C_s[0] * Cy);
	  dy = A_t[1] * Ax + B_t[1] * Bx + C_t[1] * Cx
	    - (A_s[1] * Ay + B_s[1] * By + C_s[1] * Cy);
	  dz = A_t[2] * Ax + B_t[2] * Bx + C_t[2] * Cx
	    - (A_s[2] * Ay + B_s[2] * By + C_s[2] * Cy);

	  norm = REAL_SQRT(dx * dx + dy * dy + dz * dz);

	  quad[q] = 1.0 / norm;
	}

	vl_c = tl1_c->vl;
	while (vl_c) {
	  j = vl_c->v;
	  assert(j < cols);
	  jj = ((cidx == NULL) ? j : cidx[j]);
	  for (k = 0; k < 3; ++k) {
	    if (jj == tri_sp[k]) {
	      vl_r = tl1_r->vl;
	      while (vl_r) {
		i = vl_r->v;
		assert(i < rows);
		ii = ((ridx == NULL) ? i : ridx[i]);
		for (l = 0; l < 3; ++l) {
		  if (ii == tri_tp[l]) {
		    res = base;

		    ww = wq + (l + k * 3) * nq;
		    for (q = 0; q < nq; ++q) {
		      res += ww[q] * quad[q];
		    }

		    aa[i + j * ld] += res * factor2;
		  }
		}
		vl_r = vl_r->next;
	      }
	    }
	  }
	  vl_c = vl_c->next;
	}
      }
    }
  }

  del_tri_list(tl_r);
  del_tri_list(tl_c);

  freemem(quad);
}

static void
fill_dlp_ll_laplacebem3d(const uint * ridx, const uint * cidx,
			 pcbem3d bem, bool ntrans, pamatrix N)
{
  const pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const preal gr_g = (const preal) gr->g;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const uint triangles = gr->triangles;
  plistnode *v2t = bem->v2t;
  field    *quad;
  field    *aa = N->a;
  uint      rows = N->rows;
  uint      cols = N->cols;
  longindex ld = N->ld;

  ptri_list tl_r, tl1_r, tl_c, tl1_c;
  pvert_list vl_r, vl_c;
  const real *A_t, *B_t, *C_t, *A_s, *B_s, *C_s, *ns;
  const uint *tri_t, *tri_s;
  plistnode v;
  real     *xq, *yq, *wq, *ww, *mass;
  uint      tp[3], sp[3], tri_tp[3], tri_sp[3];
  real      res, norm, base, Ax, Bx, Cx, Ay, By, Cy, tx, sx, ty, sy, dx, dy,
    dz, factor, factor2;
  uint      i, j, k, l, t, s, tt, ss, q, nq, rj, cj, ii, jj, vv;

  quad = allocfield(bem->sq->nmax);

  clear_amatrix(N);

  if (ntrans == true) {

    tl_r = NULL;
    tl_c = NULL;

    rj = 0;
    for (i = 0; i < cols; ++i) {
      ii = (ridx == NULL ? i : ridx[i]);
      for (v = v2t[ii], vv = v->data; v->next != NULL;
	   v = v->next, vv = v->data) {

	tl1_r = tl_r;
	while (tl1_r && tl1_r->t != vv) {
	  tl1_r = tl1_r->next;
	}

	if (tl1_r == NULL) {
	  tl1_r = tl_r = new_tri_list(tl_r);
	  tl_r->t = vv;
	  rj++;
	}

	tl1_r->vl = new_vert_list(tl1_r->vl);
	tl1_r->vl->v = i;
      }
    }

    cj = 0;
    for (i = 0; i < rows; ++i) {
      ii = (cidx == NULL ? i : cidx[i]);
      for (v = v2t[ii], vv = v->data; v->next != NULL;
	   v = v->next, vv = v->data) {

	tl1_c = tl_c;
	while (tl1_c && tl1_c->t != vv) {
	  tl1_c = tl1_c->next;
	}

	if (tl1_c == NULL) {
	  tl1_c = tl_c = new_tri_list(tl_c);
	  tl_c->t = vv;
	  cj++;
	}

	tl1_c->vl = new_vert_list(tl1_c->vl);
	tl1_c->vl->v = i;
      }
    }

    for (s = 0, tl1_c = tl_c; s < cj; s++, tl1_c = tl1_c->next) {
      ss = tl1_c->t;
      assert(ss < triangles);
      factor = gr_g[ss] * KERNEL_CONST_BEM3D;
      tri_s = gr_t[ss];
      ns = gr_n[ss];
      for (t = 0, tl1_r = tl_r; t < rj; t++, tl1_r = tl1_r->next) {
	tt = tl1_r->t;
	assert(tt < triangles);
	tri_t = gr_t[tt];

	if (tt != ss) {
	  factor2 = factor * gr_g[tt];

	  select_quadrature_singquad2d(bem->sq, tri_t, tri_s, tp, sp, &xq,
				       &yq, &wq, &nq, &base);

	  for (i = 0; i < 3; ++i) {
	    tri_tp[i] = tri_t[tp[i]];
	    tri_sp[i] = tri_s[sp[i]];
	  }

	  A_t = gr_x[tri_tp[0]];
	  B_t = gr_x[tri_tp[1]];
	  C_t = gr_x[tri_tp[2]];
	  A_s = gr_x[tri_sp[0]];
	  B_s = gr_x[tri_sp[1]];
	  C_s = gr_x[tri_sp[2]];

	  for (q = 0; q < nq; ++q) {
	    tx = xq[q];
	    sx = xq[q + nq];
	    ty = yq[q];
	    sy = yq[q + nq];
	    Ax = 1.0 - tx;
	    Bx = tx - sx;
	    Cx = sx;
	    Ay = 1.0 - ty;
	    By = ty - sy;
	    Cy = sy;

	    dx = A_t[0] * Ax + B_t[0] * Bx + C_t[0] * Cx
	      - (A_s[0] * Ay + B_s[0] * By + C_s[0] * Cy);
	    dy = A_t[1] * Ax + B_t[1] * Bx + C_t[1] * Cx
	      - (A_s[1] * Ay + B_s[1] * By + C_s[1] * Cy);
	    dz = A_t[2] * Ax + B_t[2] * Bx + C_t[2] * Cx
	      - (A_s[2] * Ay + B_s[2] * By + C_s[2] * Cy);

	    norm = dx * dx + dy * dy + dz * dz;

	    quad[q] = (dx * ns[0] + dy * ns[1] + dz * ns[2])
	      / (norm * REAL_SQRT(norm));
	  }

	  vl_c = tl1_c->vl;
	  while (vl_c) {
	    j = vl_c->v;
	    assert(j < rows);
	    jj = ((cidx == NULL) ? j : cidx[j]);
	    for (k = 0; k < 3; ++k) {
	      if (jj == tri_sp[k]) {
		vl_r = tl1_r->vl;
		while (vl_r) {
		  i = vl_r->v;
		  assert(i < cols);
		  ii = ((ridx == NULL) ? i : ridx[i]);
		  for (l = 0; l < 3; ++l) {
		    if (ii == tri_tp[l]) {
		      res = base;

		      ww = wq + (l + k * 3) * nq;
		      for (q = 0; q < nq; ++q) {
			res += ww[q] * quad[q];
		      }

		      aa[j + i * ld] += res * factor2;
		    }
		  }
		  vl_r = vl_r->next;
		}
	      }
	    }
	    vl_c = vl_c->next;
	  }
	}
	else {

	  for (i = 0; i < 3; ++i) {
	    tri_tp[i] = tri_t[i];
	    tri_sp[i] = tri_s[i];
	  }

	  mass = bem->mass;
	  factor2 = bem->alpha * gr_g[tt];

	  vl_c = tl1_c->vl;
	  while (vl_c) {
	    j = vl_c->v;
	    assert(j < rows);
	    jj = ((cidx == NULL) ? j : cidx[j]);
	    for (k = 0; k < 3; ++k) {
	      if (jj == tri_sp[k]) {
		vl_r = tl1_r->vl;
		while (vl_r) {
		  i = vl_r->v;
		  assert(i < cols);
		  ii = ((ridx == NULL) ? i : ridx[i]);
		  for (l = 0; l < 3; ++l) {
		    if (ii == tri_tp[l]) {
		      aa[j + i * ld] += mass[l + k * 3] * factor2;
		    }
		  }
		  vl_r = vl_r->next;
		}
	      }
	    }
	    vl_c = vl_c->next;
	  }

	}

      }
    }
  }
  else {

    tl_r = NULL;
    tl_c = NULL;

    rj = 0;
    for (i = 0; i < rows; ++i) {
      ii = (ridx == NULL ? i : ridx[i]);
      for (v = v2t[ii], vv = v->data; v->next != NULL;
	   v = v->next, vv = v->data) {

	tl1_r = tl_r;
	while (tl1_r && tl1_r->t != vv) {
	  tl1_r = tl1_r->next;
	}

	if (tl1_r == NULL) {
	  tl1_r = tl_r = new_tri_list(tl_r);
	  tl_r->t = vv;
	  rj++;
	}

	tl1_r->vl = new_vert_list(tl1_r->vl);
	tl1_r->vl->v = i;
      }
    }

    cj = 0;
    for (i = 0; i < cols; ++i) {
      ii = (cidx == NULL ? i : cidx[i]);
      for (v = v2t[ii], vv = v->data; v->next != NULL;
	   v = v->next, vv = v->data) {

	tl1_c = tl_c;
	while (tl1_c && tl1_c->t != vv) {
	  tl1_c = tl1_c->next;
	}

	if (tl1_c == NULL) {
	  tl1_c = tl_c = new_tri_list(tl_c);
	  tl_c->t = vv;
	  cj++;
	}

	tl1_c->vl = new_vert_list(tl1_c->vl);
	tl1_c->vl->v = i;
      }
    }

    for (s = 0, tl1_c = tl_c; s < cj; s++, tl1_c = tl1_c->next) {
      ss = tl1_c->t;
      assert(ss < triangles);
      factor = gr_g[ss] * KERNEL_CONST_BEM3D;
      tri_s = gr_t[ss];
      ns = gr_n[ss];
      for (t = 0, tl1_r = tl_r; t < rj; t++, tl1_r = tl1_r->next) {
	tt = tl1_r->t;
	assert(tt < triangles);
	tri_t = gr_t[tt];

	if (tt != ss) {
	  factor2 = factor * gr_g[tt];

	  select_quadrature_singquad2d(bem->sq, tri_t, tri_s, tp, sp, &xq,
				       &yq, &wq, &nq, &base);

	  for (i = 0; i < 3; ++i) {
	    tri_tp[i] = tri_t[tp[i]];
	    tri_sp[i] = tri_s[sp[i]];
	  }

	  A_t = gr_x[tri_tp[0]];
	  B_t = gr_x[tri_tp[1]];
	  C_t = gr_x[tri_tp[2]];
	  A_s = gr_x[tri_sp[0]];
	  B_s = gr_x[tri_sp[1]];
	  C_s = gr_x[tri_sp[2]];

	  for (q = 0; q < nq; ++q) {
	    tx = xq[q];
	    sx = xq[q + nq];
	    ty = yq[q];
	    sy = yq[q + nq];
	    Ax = 1.0 - tx;
	    Bx = tx - sx;
	    Cx = sx;
	    Ay = 1.0 - ty;
	    By = ty - sy;
	    Cy = sy;

	    dx = A_t[0] * Ax + B_t[0] * Bx + C_t[0] * Cx
	      - (A_s[0] * Ay + B_s[0] * By + C_s[0] * Cy);
	    dy = A_t[1] * Ax + B_t[1] * Bx + C_t[1] * Cx
	      - (A_s[1] * Ay + B_s[1] * By + C_s[1] * Cy);
	    dz = A_t[2] * Ax + B_t[2] * Bx + C_t[2] * Cx
	      - (A_s[2] * Ay + B_s[2] * By + C_s[2] * Cy);

	    norm = dx * dx + dy * dy + dz * dz;

	    quad[q] = (dx * ns[0] + dy * ns[1] + dz * ns[2])
	      / (norm * REAL_SQRT(norm));
	  }

	  vl_c = tl1_c->vl;
	  while (vl_c) {
	    j = vl_c->v;
	    assert(j < cols);
	    jj = ((cidx == NULL) ? j : cidx[j]);
	    for (k = 0; k < 3; ++k) {
	      if (jj == tri_sp[k]) {
		vl_r = tl1_r->vl;
		while (vl_r) {
		  i = vl_r->v;
		  assert(i < rows);
		  ii = ((ridx == NULL) ? i : ridx[i]);
		  for (l = 0; l < 3; ++l) {
		    if (ii == tri_tp[l]) {
		      res = base;

		      ww = wq + (l + k * 3) * nq;
		      for (q = 0; q < nq; ++q) {
			res += ww[q] * quad[q];
		      }

		      aa[i + j * ld] += res * factor2;
		    }
		  }
		  vl_r = vl_r->next;
		}
	      }
	    }
	    vl_c = vl_c->next;
	  }
	}
	else {

	  for (i = 0; i < 3; ++i) {
	    tri_tp[i] = tri_t[i];
	    tri_sp[i] = tri_s[i];
	  }

	  mass = bem->mass;
	  factor2 = bem->alpha * gr_g[tt];

	  vl_c = tl1_c->vl;
	  while (vl_c) {
	    j = vl_c->v;
	    assert(j < cols);
	    jj = ((cidx == NULL) ? j : cidx[j]);
	    for (k = 0; k < 3; ++k) {
	      if (jj == tri_sp[k]) {
		vl_r = tl1_r->vl;
		while (vl_r) {
		  i = vl_r->v;
		  assert(i < rows);
		  ii = ((ridx == NULL) ? i : ridx[i]);
		  for (l = 0; l < 3; ++l) {
		    if (ii == tri_tp[l]) {
		      aa[i + j * ld] += mass[l + k * 3] * factor2;
		    }
		  }
		  vl_r = vl_r->next;
		}
	      }
	    }
	    vl_c = vl_c->next;
	  }
	}
      }
    }
  }

  del_tri_list(tl_r);
  del_tri_list(tl_c);

  freemem(quad);
}

static void
fill_kernel_laplacebem3d(const real(*X)[3], const real(*Y)[3], pamatrix V)
{
  uint      rows = V->rows;
  uint      cols = V->cols;
  uint      ld = V->ld;

  uint      i, j;
  real      dx, dy, dz;

  for (j = 0; j < cols; ++j) {
    for (i = 0; i < rows; ++i) {
      dx = X[i][0] - Y[j][0];
      dy = X[i][1] - Y[j][1];
      dz = X[i][2] - Y[j][2];

      V->a[i + j * ld] = KERNEL_CONST_BEM3D
	/ REAL_SQRT(dx * dx + dy * dy + dz * dz);
    }
  }
}

static void
fill_dny_kernel_laplacebem3d(const real(*X)[3], const real(*Y)[3],
			     const real(*NY)[3], pamatrix V)
{
  uint      rows = V->rows;
  uint      cols = V->cols;
  uint      ld = V->ld;

  uint      i, j;
  real      dx, dy, dz, norm2;

  for (j = 0; j < cols; ++j) {
    for (i = 0; i < rows; ++i) {
      dx = X[i][0] - Y[j][0];
      dy = X[i][1] - Y[j][1];
      dz = X[i][2] - Y[j][2];

      norm2 = 1.0 / (dx * dx + dy * dy + dz * dz);

      V->a[i + j * ld] = KERNEL_CONST_BEM3D
	* (dx * NY[j][0] + dy * NY[j][1] + dz * NY[j][2]) * norm2
	* REAL_SQRT(norm2);
    }
  }
}

static void
fill_dnx_dny_kernel_laplacebem3d(const real(*X)[3],
				 const real(*NX)[3], const real(*Y)[3],
				 const real(*NY)[3], pamatrix V)
{
  uint      rows = V->rows;
  uint      cols = V->cols;
  uint      ld = V->ld;

  uint      i, j;
  real      h[3];
  real      dx, dy, dz, norm2, dot1, kernel;

  for (j = 0; j < cols; ++j) {
    for (i = 0; i < rows; ++i) {
      dx = X[i][0] - Y[j][0];
      dy = X[i][1] - Y[j][1];
      dz = X[i][2] - Y[j][2];

      norm2 = 1.0 / (dx * dx + dy * dy + dz * dz);

      dot1 = -3.0 * norm2 * (NY[j][0] * dx + NY[j][1] * dy + NY[j][2] * dz);
      h[0] = dx * dot1 + NY[j][0];
      h[1] = dy * dot1 + NY[j][1];
      h[2] = dz * dot1 + NY[j][2];
      kernel = REAL_SQRT(norm2) * norm2
	* (h[0] * NX[i][0] + h[1] * NX[i][1] + h[2] * NX[i][2]);

      V->a[i + j * ld] = KERNEL_CONST_BEM3D * kernel;
    }
  }
}

static void
fill_kernel_c_laplacebem3d(const uint * idx, const real(*Z)[3],
			   pcbem3d bem, pamatrix V)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const preal gr_g = (const preal) gr->g;
  uint      rows = V->rows;
  uint      cols = V->cols;
  uint      ld = V->ld;

  uint      nq = bem->sq->n_single;
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single + 3 * nq;

  const real *A, *B, *C;
  uint      s, ss, i, q;
  real      gs_fac, sum, kernel, x, y, z, tx, sx, Ax, Bx, Cx;

  /*
   *  integrate kernel function over first variable with constant basisfunctions
   */

  for (s = 0; s < rows; ++s) {
    ss = (idx == NULL ? s : idx[s]);
    gs_fac = gr_g[ss] * KERNEL_CONST_BEM3D;
    A = gr_x[gr_t[ss][0]];
    B = gr_x[gr_t[ss][1]];
    C = gr_x[gr_t[ss][2]];

    for (i = 0; i < cols; ++i) {

      sum = 0.0;

      for (q = 0; q < nq; ++q) {
	tx = xx[q];
	sx = yy[q];
	Ax = 1.0 - tx;
	Bx = tx - sx;
	Cx = sx;

	x = Z[i][0] - (A[0] * Ax + B[0] * Bx + C[0] * Cx);
	y = Z[i][1] - (A[1] * Ax + B[1] * Bx + C[1] * Cx);
	z = Z[i][2] - (A[2] * Ax + B[2] * Bx + C[2] * Cx);

	kernel = REAL_SQRT(x * x + y * y + z * z);

	sum += ww[q] / kernel;
      }

      V->a[s + i * (longindex) ld] = sum * gs_fac;
    }
  }
}

static void
fill_kernel_l_laplacebem3d(const uint * idx, const real(*Z)[3],
			   pcbem3d bem, pamatrix V)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const preal gr_g = (const preal) gr->g;
  plistnode *v2t = bem->v2t;
  uint      rows = V->rows;
  uint      cols = V->cols;
  field    *aa = V->a;
  uint      ld = V->ld;
  uint      nq = bem->sq->n_single;
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single;
  real      base = bem->sq->base_single;
  field    *quad;

  ptri_list tl, tl1;
  pvert_list vl;
  const real *A, *B, *C;
  uint      tri_tp[3];
  plistnode v;
  uint      t, i, j, k, q, rj;
  real      gt_fac, sum, kernel, x, y, z, tx, sx, Ax, Bx, Cx;
  longindex ii, tt, vv;

  quad = allocfield(nq);

  clear_amatrix(V);

  tl = NULL;

  rj = 0;
  for (i = 0; i < rows; ++i) {
    ii = (idx == NULL ? i : idx[i]);
    for (v = v2t[ii], vv = v->data; v->next != NULL;
	 v = v->next, vv = v->data) {

      tl1 = tl;
      while (tl1 && tl1->t != vv) {
	tl1 = tl1->next;
      }

      if (tl1 == NULL) {
	tl1 = tl = new_tri_list(tl);
	tl->t = vv;
	rj++;
      }

      tl1->vl = new_vert_list(tl1->vl);
      tl1->vl->v = i;
    }
  }

  for (t = 0, tl1 = tl; t < rj; t++, tl1 = tl1->next) {
    tt = tl1->t;
    gt_fac = gr_g[tt] * KERNEL_CONST_BEM3D;
    A = gr_x[gr_t[tt][0]];
    B = gr_x[gr_t[tt][1]];
    C = gr_x[gr_t[tt][2]];

    for (i = 0; i < 3; ++i) {
      tri_tp[i] = gr_t[tt][i];
    }

    for (j = 0; j < cols; ++j) {

      for (q = 0; q < nq; ++q) {
	tx = xx[q];
	sx = yy[q];
	Ax = 1.0 - tx;
	Bx = tx - sx;
	Cx = sx;

	x = A[0] * Ax + B[0] * Bx + C[0] * Cx - Z[j][0];
	y = A[1] * Ax + B[1] * Bx + C[1] * Cx - Z[j][1];
	z = A[2] * Ax + B[2] * Bx + C[2] * Cx - Z[j][2];

	kernel = REAL_SQRT(x * x + y * y + z * z);

	quad[q] = 1.0 / kernel;
      }

      ww = bem->sq->w_single;
      vl = tl1->vl;
      while (vl) {
	k = vl->v;
	if (k < rows) {
	  ii = idx == NULL ? k : idx[k];
	  for (i = 0; i < 3; ++i) {
	    if (ii == tri_tp[i]) {
	      sum = base;

	      for (q = 0; q < nq; ++q) {
		sum += ww[q] * quad[q];
	      }

	      aa[k + j * ld] += sum * gt_fac;
	    }
	    ww += nq;
	  }
	  ww = bem->sq->w_single;
	}
	vl = vl->next;
      }

    }
  }

  del_tri_list(tl);

  freemem(quad);
}

static void
fill_dnz_kernel_c_laplacebem3d(const uint * idx, const real(*Z)[3],
			       const real(*N)[3], pcbem3d bem, pamatrix V)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const preal gr_g = (const preal) gr->g;
  uint      rows = V->rows;
  uint      cols = V->cols;
  uint      ld = V->ld;

  uint      nq = bem->sq->n_single;
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single + 3 * nq;

  const real *A, *B, *C;
  uint      t, tt, i, q;
  real      gt_fac, sum, kernel, dx, dy, dz, tx, sx, Ax, Bx, Cx;

  for (t = 0; t < rows; ++t) {
    tt = (idx == NULL ? t : idx[t]);
    gt_fac = gr_g[tt] * KERNEL_CONST_BEM3D;
    A = gr_x[gr_t[tt][0]];
    B = gr_x[gr_t[tt][1]];
    C = gr_x[gr_t[tt][2]];

    for (i = 0; i < cols; ++i) {

      sum = 0.0;

      for (q = 0; q < nq; ++q) {
	tx = xx[q];
	sx = yy[q];
	Ax = 1.0 - tx;
	Bx = tx - sx;
	Cx = sx;

	dx = A[0] * Ax + B[0] * Bx + C[0] * Cx - Z[i][0];
	dy = A[1] * Ax + B[1] * Bx + C[1] * Cx - Z[i][1];
	dz = A[2] * Ax + B[2] * Bx + C[2] * Cx - Z[i][2];

	kernel = dx * dx + dy * dy + dz * dz;

	sum += ww[q] * (dx * N[i][0] + dy * N[i][1] + dz * N[i][2])
	  / (kernel * REAL_SQRT(kernel));

      }

      V->a[t + i * ld] = sum * gt_fac;
    }
  }
}

static void
fill_dnz_kernel_l_laplacebem3d(const uint * idx, const real(*Z)[3],
			       const real(*N)[3], pcbem3d bem, pamatrix V)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const preal gr_g = (const preal) gr->g;
  plistnode *v2t = bem->v2t;
  uint      rows = V->rows;
  uint      cols = V->cols;
  field    *aa = V->a;
  uint      ld = V->ld;
  uint      nq = bem->sq->n_single;
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single;
  real      base = bem->sq->base_single;
  field    *quad;

  ptri_list tl, tl1;
  pvert_list vl;
  const real *A, *B, *C;
  uint      tri_tp[3];
  plistnode v;
  uint      t, i, j, k, q, rj;
  real      gt_fac, sum, kernel, x, y, z, tx, sx, Ax, Bx, Cx;
  longindex ii, tt, vv;

  quad = allocfield(nq);

  clear_amatrix(V);

  tl = NULL;

  rj = 0;
  for (i = 0; i < rows; ++i) {
    ii = (idx == NULL ? i : idx[i]);
    for (v = v2t[ii], vv = v->data; v->next != NULL;
	 v = v->next, vv = v->data) {

      tl1 = tl;
      while (tl1 && tl1->t != vv) {
	tl1 = tl1->next;
      }

      if (tl1 == NULL) {
	tl1 = tl = new_tri_list(tl);
	tl->t = vv;
	rj++;
      }

      tl1->vl = new_vert_list(tl1->vl);
      tl1->vl->v = i;
    }
  }

  for (t = 0, tl1 = tl; t < rj; t++, tl1 = tl1->next) {
    tt = tl1->t;
    gt_fac = gr_g[tt] * KERNEL_CONST_BEM3D;
    A = gr_x[gr_t[tt][0]];
    B = gr_x[gr_t[tt][1]];
    C = gr_x[gr_t[tt][2]];

    for (i = 0; i < 3; ++i) {
      tri_tp[i] = gr_t[tt][i];
    }

    for (j = 0; j < cols; ++j) {

      for (q = 0; q < nq; ++q) {
	tx = xx[q];
	sx = yy[q];
	Ax = 1.0 - tx;
	Bx = tx - sx;
	Cx = sx;

	x = A[0] * Ax + B[0] * Bx + C[0] * Cx - Z[j][0];
	y = A[1] * Ax + B[1] * Bx + C[1] * Cx - Z[j][1];
	z = A[2] * Ax + B[2] * Bx + C[2] * Cx - Z[j][2];

	kernel = x * x + y * y + z * z;

	quad[q] = (x * N[j][0] + y * N[j][1] + z * N[j][2])
	  / (kernel * REAL_SQRT(kernel));
      }

      ww = bem->sq->w_single;
      vl = tl1->vl;
      while (vl) {
	k = vl->v;
	if (k < rows) {
	  ii = idx == NULL ? k : idx[k];
	  for (i = 0; i < 3; ++i) {
	    if (ii == tri_tp[i]) {
	      sum = base;

	      for (q = 0; q < nq; ++q) {
		sum += ww[q] * quad[q];
	      }

	      aa[k + j * ld] += sum * gt_fac;
	    }
	    ww += nq;
	  }
	  ww = bem->sq->w_single;
	}
	vl = vl->next;
      }
    }
  }

  del_tri_list(tl);

  freemem(quad);
}

static void
fill_dnzdcol_kernel_c_laplacebem3d(const uint * idx,
				   const real(*Z)[3], const real(*N)[3],
				   pcbem3d bem, pamatrix V)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const preal gr_g = (const preal) gr->g;
  uint      rows = V->rows;
  uint      cols = V->cols;
  uint      ld = V->ld;

  uint      nq = bem->sq->n_single;
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single + 3 * nq;

  const real *A, *B, *C, *ns;
  uint      s, ss, i, q;
  real      gs_fac, sum, kernel, norm2, dot1, tx, sx, Ax, Bx, Cx;
  real      dxy[3], h[3];

  for (s = 0; s < rows; ++s) {
    ss = (idx == NULL ? s : idx[s]);
    gs_fac = gr_g[ss] * KERNEL_CONST_BEM3D;
    A = gr_x[gr_t[ss][0]];
    B = gr_x[gr_t[ss][1]];
    C = gr_x[gr_t[ss][2]];
    ns = gr_n[ss];

    for (i = 0; i < cols; ++i) {

      sum = 0.0;

      for (q = 0; q < nq; ++q) {
	tx = xx[q];
	sx = yy[q];
	Ax = 1.0 - tx;
	Bx = tx - sx;
	Cx = sx;

	dxy[0] = Z[i][0] - (A[0] * Ax + B[0] * Bx + C[0] * Cx);
	dxy[1] = Z[i][1] - (A[1] * Ax + B[1] * Bx + C[1] * Cx);
	dxy[2] = Z[i][2] - (A[2] * Ax + B[2] * Bx + C[2] * Cx);

	norm2 = 1.0 / (dxy[0] * dxy[0] + dxy[1] * dxy[1] + dxy[2] * dxy[2]);

	dot1 = -3.0 * norm2
	  * (ns[0] * dxy[0] + ns[1] * dxy[1] + ns[2] * dxy[2]);
	h[0] = dxy[0] * dot1 + ns[0];
	h[1] = dxy[1] * dot1 + ns[1];
	h[2] = dxy[2] * dot1 + ns[2];
	kernel = REAL_SQRT(norm2) * norm2
	  * (h[0] * N[i][0] + h[1] * N[i][1] + h[2] * N[i][2]);

	sum += ww[q] * kernel;

      }

      V->a[s + i * ld] = sum * gs_fac;
    }
  }
}

static void
fill_dnzdcol_kernel_l_laplacebem3d(const uint * idx,
				   const real(*Z)[3], const real(*N)[3],
				   pcbem3d bem, pamatrix V)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  plistnode *v2t = bem->v2t;
  uint      rows = V->rows;
  uint      cols = V->cols;
  field    *aa = V->a;
  uint      ld = V->ld;
  uint      nq = bem->sq->n_single;
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single;
  real      base = bem->sq->base_single;
  field    *quad;

  ptri_list tl, tl1;
  pvert_list vl;
  const real *A, *B, *C, *nt;
  uint      tri_tp[3];
  plistnode v;
  uint      t, i, j, k, q, rj;
  real      gt_fac, sum, kernel, norm2, dot1, tx, sx, Ax, Bx, Cx;
  real      dxy[3], h[3];
  longindex ii, tt, vv;

  quad = allocfield(nq);

  clear_amatrix(V);

  tl = NULL;

  rj = 0;
  for (i = 0; i < rows; ++i) {
    ii = (idx == NULL ? i : idx[i]);
    for (v = v2t[ii], vv = v->data; v->next != NULL;
	 v = v->next, vv = v->data) {

      tl1 = tl;
      while (tl1 && tl1->t != vv) {
	tl1 = tl1->next;
      }

      if (tl1 == NULL) {
	tl1 = tl = new_tri_list(tl);
	tl->t = vv;
	rj++;
      }

      tl1->vl = new_vert_list(tl1->vl);
      tl1->vl->v = i;
    }
  }

  for (t = 0, tl1 = tl; t < rj; t++, tl1 = tl1->next) {
    tt = tl1->t;
    gt_fac = gr_g[tt] * KERNEL_CONST_BEM3D;
    nt = gr_n[tt];
    A = gr_x[gr_t[tt][0]];
    B = gr_x[gr_t[tt][1]];
    C = gr_x[gr_t[tt][2]];

    for (i = 0; i < 3; ++i) {
      tri_tp[i] = gr_t[tt][i];
    }

    for (j = 0; j < cols; ++j) {

      for (q = 0; q < nq; ++q) {
	tx = xx[q];
	sx = yy[q];
	Ax = 1.0 - tx;
	Bx = tx - sx;
	Cx = sx;

	dxy[0] = Z[j][0] - (A[0] * Ax + B[0] * Bx + C[0] * Cx);
	dxy[1] = Z[j][1] - (A[1] * Ax + B[1] * Bx + C[1] * Cx);
	dxy[2] = Z[j][2] - (A[2] * Ax + B[2] * Bx + C[2] * Cx);

	norm2 = 1.0 / (dxy[0] * dxy[0] + dxy[1] * dxy[1] + dxy[2] * dxy[2]);

	dot1 = -3.0 * norm2
	  * (nt[0] * dxy[0] + nt[1] * dxy[1] + nt[2] * dxy[2]);
	h[0] = dxy[0] * dot1 + nt[0];
	h[1] = dxy[1] * dot1 + nt[1];
	h[2] = dxy[2] * dot1 + nt[2];
	kernel = REAL_SQRT(norm2) * norm2
	  * (h[0] * N[j][0] + h[1] * N[j][1] + h[2] * N[j][2]);

	quad[q] = kernel;
      }

      ww = bem->sq->w_single;
      vl = tl1->vl;
      while (vl) {
	k = vl->v;
	if (k < rows) {
	  ii = idx == NULL ? k : idx[k];
	  for (i = 0; i < 3; ++i) {
	    if (ii == tri_tp[i]) {
	      sum = base;

	      for (q = 0; q < nq; ++q) {
		sum += ww[q] * quad[q];
	      }

	      aa[k + j * ld] += sum * gt_fac;
	    }
	    ww += nq;
	  }
	  ww = bem->sq->w_single;
	}
	vl = vl->next;
      }
    }
  }

  del_tri_list(tl);

  freemem(quad);
}

static void
fill_dcol_kernel_col_c_laplacebem3d(const uint * idx,
				    const real(*Z)[3], pcbem3d bem,
				    pamatrix V)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const preal gr_g = (const preal) gr->g;
  uint      rows = V->rows;
  uint      cols = V->cols;
  uint      ld = V->ld;

  uint      nq = bem->sq->n_single;
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single + 3 * nq;

  const real *A, *B, *C, *ns;
  uint      s, ss, i, q;
  real      gs_fac, sum, kernel, dx, dy, dz, tx, sx, Ax, Bx, Cx;

  /*
   *  integrate kernel function over first variable with constant basisfunctions
   */

  for (s = 0; s < rows; ++s) {
    ss = (idx == NULL ? s : idx[s]);
    gs_fac = gr_g[ss] * KERNEL_CONST_BEM3D;
    ns = gr_n[ss];
    A = gr_x[gr_t[ss][0]];
    B = gr_x[gr_t[ss][1]];
    C = gr_x[gr_t[ss][2]];

    for (i = 0; i < cols; ++i) {

      sum = 0.0;

      for (q = 0; q < nq; ++q) {
	tx = xx[q];
	sx = yy[q];
	Ax = 1.0 - tx;
	Bx = tx - sx;
	Cx = sx;

	dx = Z[i][0] - (A[0] * Ax + B[0] * Bx + C[0] * Cx);
	dy = Z[i][1] - (A[1] * Ax + B[1] * Bx + C[1] * Cx);
	dz = Z[i][2] - (A[2] * Ax + B[2] * Bx + C[2] * Cx);

	kernel = dx * dx + dy * dy + dz * dz;

	sum += ww[q] * (dx * ns[0] + dy * ns[1] + dz * ns[2])
	  / (kernel * REAL_SQRT(kernel));

      }

      V->a[s + i * ld] = sum * gs_fac;
    }
  }
}

static void
fill_dcol_kernel_col_l_laplacebem3d(const uint * idx,
				    const real(*Z)[3], pcbem3d bem,
				    pamatrix V)
{
  pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const preal gr_g = (const preal) gr->g;
  plistnode *v2t = bem->v2t;
  uint      rows = V->rows;
  uint      cols = V->cols;
  field    *aa = V->a;
  uint      ld = V->ld;
  uint      nq = bem->sq->n_single;
  real     *xx = bem->sq->x_single;
  real     *yy = bem->sq->y_single;
  real     *ww = bem->sq->w_single;
  real      base = bem->sq->base_single;
  field    *quad;

  ptri_list tl, tl1;
  pvert_list vl;
  const real *A, *B, *C, *nt;
  uint      tri_tp[3];
  plistnode v;
  uint      t, i, j, k, q, rj;
  real      gt_fac, sum, kernel, dx, dy, dz, tx, sx, Ax, Bx, Cx;
  longindex ii, tt, vv;

  quad = allocfield(nq);

  clear_amatrix(V);

  tl = NULL;

  rj = 0;
  for (i = 0; i < rows; ++i) {
    ii = (idx == NULL ? i : idx[i]);
    for (v = v2t[ii], vv = v->data; v->next != NULL;
	 v = v->next, vv = v->data) {

      tl1 = tl;
      while (tl1 && tl1->t != vv) {
	tl1 = tl1->next;
      }

      if (tl1 == NULL) {
	tl1 = tl = new_tri_list(tl);
	tl->t = vv;
	rj++;
      }

      tl1->vl = new_vert_list(tl1->vl);
      tl1->vl->v = i;
    }
  }

  for (t = 0, tl1 = tl; t < rj; t++, tl1 = tl1->next) {
    tt = tl1->t;
    gt_fac = gr_g[tt] * KERNEL_CONST_BEM3D;
    nt = gr_n[tt];
    A = gr_x[gr_t[tt][0]];
    B = gr_x[gr_t[tt][1]];
    C = gr_x[gr_t[tt][2]];

    for (i = 0; i < 3; ++i) {
      tri_tp[i] = gr_t[tt][i];
    }

    for (j = 0; j < cols; ++j) {

      for (q = 0; q < nq; ++q) {
	tx = xx[q];
	sx = yy[q];
	Ax = 1.0 - tx;
	Bx = tx - sx;
	Cx = sx;

	dx = Z[j][0] - (A[0] * Ax + B[0] * Bx + C[0] * Cx);
	dy = Z[j][1] - (A[1] * Ax + B[1] * Bx + C[1] * Cx);
	dz = Z[j][2] - (A[2] * Ax + B[2] * Bx + C[2] * Cx);

	kernel = dx * dx + dy * dy + dz * dz;

	quad[q] = (dx * nt[0] + dy * nt[1] + dz * nt[2])
	  / (kernel * REAL_SQRT(kernel));
      }

      ww = bem->sq->w_single;
      vl = tl1->vl;
      while (vl) {
	k = vl->v;
	if (k < rows) {
	  ii = idx == NULL ? k : idx[k];
	  for (i = 0; i < 3; ++i) {
	    if (ii == tri_tp[i]) {
	      sum = base;

	      for (q = 0; q < nq; ++q) {
		sum += ww[q] * quad[q];
	      }

	      aa[k + j * ld] += sum * gt_fac;
	    }
	    ww += nq;
	  }
	  ww = bem->sq->w_single;
	}
	vl = vl->next;
      }
    }
  }

  del_tri_list(tl);

  freemem(quad);
}

pbem3d
new_slp_laplace_bem3d(pcsurface3d gr, uint q, basisfunctionbem3d basis)
{
  pkernelbem3d kernels;

  pbem3d    bem;

  bem = new_bem3d(gr);
  kernels = bem->kernels;

  bem->sq = build_singquad2d(q, q + 2);

  if (basis == BASIS_LINEAR_BEM3D) {
    setup_vertex_to_triangle_map_bem3d(bem);
  }

  bem->basis_neumann = basis;

  kernels->fundamental = fill_kernel_laplacebem3d;
  kernels->dny_fundamental = fill_dny_kernel_laplacebem3d;
  kernels->dnx_dny_fundamental = fill_dnx_dny_kernel_laplacebem3d;

  if (basis == BASIS_CONSTANT_BEM3D) {
    bem->N_neumann = gr->triangles;

    bem->nearfield = fill_slp_cc_laplacebem3d;

    kernels->lagrange_row = assemble_bem3d_lagrange_const_amatrix;
    kernels->lagrange_col = assemble_bem3d_lagrange_const_amatrix;

    kernels->fundamental_row = fill_kernel_c_laplacebem3d;
    kernels->fundamental_col = fill_kernel_c_laplacebem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_c_laplacebem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_c_laplacebem3d;

    kernels->kernel_row = fill_kernel_c_laplacebem3d;
    kernels->kernel_col = fill_kernel_c_laplacebem3d;
    kernels->dnz_kernel_row = fill_dnz_kernel_c_laplacebem3d;
    kernels->dnz_kernel_col = fill_dnz_kernel_c_laplacebem3d;
  }
  else {
    assert(basis == BASIS_LINEAR_BEM3D);

    weight_basisfunc_ll_singquad2d(bem->sq->x_id, bem->sq->y_id,
				   bem->sq->w_id, bem->sq->n_id);
    weight_basisfunc_ll_singquad2d(bem->sq->x_edge, bem->sq->y_edge,
				   bem->sq->w_edge, bem->sq->n_edge);
    weight_basisfunc_ll_singquad2d(bem->sq->x_vert, bem->sq->y_vert,
				   bem->sq->w_vert, bem->sq->n_vert);
    weight_basisfunc_ll_singquad2d(bem->sq->x_dist, bem->sq->y_dist,
				   bem->sq->w_dist, bem->sq->n_dist);
    weight_basisfunc_l_singquad2d(bem->sq->x_single, bem->sq->y_single,
				  bem->sq->w_single, bem->sq->n_single);

    bem->N_neumann = gr->vertices;

    bem->nearfield = fill_slp_ll_laplacebem3d;

    kernels->lagrange_row = assemble_bem3d_lagrange_linear_amatrix;
    kernels->lagrange_col = assemble_bem3d_lagrange_linear_amatrix;

    kernels->fundamental_row = fill_kernel_l_laplacebem3d;
    kernels->fundamental_col = fill_kernel_l_laplacebem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_l_laplacebem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_l_laplacebem3d;

    kernels->kernel_row = fill_kernel_l_laplacebem3d;
    kernels->kernel_col = fill_kernel_l_laplacebem3d;
    kernels->dnz_kernel_row = fill_dnz_kernel_l_laplacebem3d;
    kernels->dnz_kernel_col = fill_dnz_kernel_l_laplacebem3d;
  }

  return bem;
}

pbem3d
new_dlp_laplace_bem3d(pcsurface3d gr, uint q,
		      basisfunctionbem3d basis_neumann,
		      basisfunctionbem3d basis_dirichlet, field alpha)
{
  pkernelbem3d kernels;

  pbem3d    bem;

  bem = new_bem3d(gr);
  kernels = bem->kernels;

  bem->sq = build_singquad2d(q, q + 2);

  if (basis_neumann == BASIS_LINEAR_BEM3D || basis_dirichlet
      == BASIS_LINEAR_BEM3D) {
    setup_vertex_to_triangle_map_bem3d(bem);
  }

  bem->basis_neumann = basis_neumann;
  bem->basis_dirichlet = basis_dirichlet;

  bem->alpha = alpha;

  kernels->fundamental = fill_kernel_laplacebem3d;
  kernels->dny_fundamental = fill_dny_kernel_laplacebem3d;
  kernels->dnx_dny_fundamental = fill_dnx_dny_kernel_laplacebem3d;

  if (basis_neumann == BASIS_CONSTANT_BEM3D) {
    bem->N_neumann = gr->triangles;
  }
  else {
    assert(basis_neumann == BASIS_LINEAR_BEM3D);
    bem->N_neumann = gr->vertices;
  }
  if (basis_dirichlet == BASIS_CONSTANT_BEM3D) {
    bem->N_dirichlet = gr->triangles;
  }
  else {
    assert(basis_dirichlet == BASIS_LINEAR_BEM3D);
    bem->N_dirichlet = gr->vertices;
  }

  if (basis_neumann == BASIS_CONSTANT_BEM3D && basis_dirichlet
      == BASIS_CONSTANT_BEM3D) {
    bem->nearfield = fill_dlp_cc_laplacebem3d;

    kernels->lagrange_row = assemble_bem3d_lagrange_const_amatrix;
    kernels->lagrange_col = assemble_bem3d_dn_lagrange_const_amatrix;

    kernels->fundamental_row = fill_kernel_c_laplacebem3d;
    kernels->fundamental_col = fill_kernel_c_laplacebem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_c_laplacebem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_c_laplacebem3d;

    kernels->kernel_row = fill_kernel_c_laplacebem3d;
    kernels->kernel_col = fill_dcol_kernel_col_c_laplacebem3d;
    kernels->dnz_kernel_row = fill_dnz_kernel_c_laplacebem3d;
    kernels->dnz_kernel_col = fill_dnzdcol_kernel_c_laplacebem3d;
  }
  else if (basis_neumann == BASIS_LINEAR_BEM3D
	   && basis_dirichlet == BASIS_CONSTANT_BEM3D) {
    bem->mass = allocfield(3);
    /* TODO MASS-MATRIX */

  }
  else if (basis_neumann == BASIS_CONSTANT_BEM3D
	   && basis_dirichlet == BASIS_LINEAR_BEM3D) {

    weight_basisfunc_cl_singquad2d(bem->sq->x_id, bem->sq->y_id,
				   bem->sq->w_id, bem->sq->n_id);
    weight_basisfunc_cl_singquad2d(bem->sq->x_edge, bem->sq->y_edge,
				   bem->sq->w_edge, bem->sq->n_edge);
    weight_basisfunc_cl_singquad2d(bem->sq->x_vert, bem->sq->y_vert,
				   bem->sq->w_vert, bem->sq->n_vert);
    weight_basisfunc_cl_singquad2d(bem->sq->x_dist, bem->sq->y_dist,
				   bem->sq->w_dist, bem->sq->n_dist);
    weight_basisfunc_l_singquad2d(bem->sq->x_single, bem->sq->y_single,
				  bem->sq->w_single, bem->sq->n_single);

    bem->nearfield = fill_dlp_cl_laplacebem3d;

    kernels->lagrange_row = assemble_bem3d_lagrange_const_amatrix;
    kernels->lagrange_col = assemble_bem3d_dn_lagrange_linear_amatrix;

    kernels->fundamental_row = fill_kernel_c_laplacebem3d;
    kernels->fundamental_col = fill_kernel_l_laplacebem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_c_laplacebem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_l_laplacebem3d;

    kernels->kernel_row = fill_kernel_c_laplacebem3d;
    kernels->kernel_col = fill_dcol_kernel_col_l_laplacebem3d;
    kernels->dnz_kernel_row = fill_dnz_kernel_c_laplacebem3d;
    kernels->dnz_kernel_col = fill_dnzdcol_kernel_l_laplacebem3d;

    bem->mass = allocfield(3);
    bem->mass[0] = 1.0 / 6.0;
    bem->mass[1] = 1.0 / 6.0;
    bem->mass[2] = 1.0 / 6.0;

  }
  else {
    assert(basis_neumann == BASIS_LINEAR_BEM3D && basis_dirichlet
	   == BASIS_LINEAR_BEM3D);

    weight_basisfunc_ll_singquad2d(bem->sq->x_id, bem->sq->y_id,
				   bem->sq->w_id, bem->sq->n_id);
    weight_basisfunc_ll_singquad2d(bem->sq->x_edge, bem->sq->y_edge,
				   bem->sq->w_edge, bem->sq->n_edge);
    weight_basisfunc_ll_singquad2d(bem->sq->x_vert, bem->sq->y_vert,
				   bem->sq->w_vert, bem->sq->n_vert);
    weight_basisfunc_ll_singquad2d(bem->sq->x_dist, bem->sq->y_dist,
				   bem->sq->w_dist, bem->sq->n_dist);
    weight_basisfunc_l_singquad2d(bem->sq->x_single, bem->sq->y_single,
				  bem->sq->w_single, bem->sq->n_single);

    bem->nearfield = fill_dlp_ll_laplacebem3d;

    kernels->lagrange_row = assemble_bem3d_lagrange_linear_amatrix;
    kernels->lagrange_col = assemble_bem3d_dn_lagrange_linear_amatrix;

    kernels->fundamental_row = fill_kernel_l_laplacebem3d;
    kernels->fundamental_col = fill_kernel_l_laplacebem3d;
    kernels->dnz_fundamental_row = fill_dnz_kernel_l_laplacebem3d;
    kernels->dnz_fundamental_col = fill_dnz_kernel_l_laplacebem3d;

    kernels->kernel_row = fill_kernel_l_laplacebem3d;
    kernels->kernel_col = fill_dcol_kernel_col_l_laplacebem3d;
    kernels->dnz_kernel_row = NULL;
    kernels->dnz_kernel_col = fill_dnzdcol_kernel_l_laplacebem3d;

    bem->mass = allocfield(9);
    bem->mass[0] = 1.0 / 12.0;
    bem->mass[1] = 1.0 / 24.0;
    bem->mass[2] = 1.0 / 24.0;
    bem->mass[3] = 1.0 / 24.0;
    bem->mass[4] = 1.0 / 12.0;
    bem->mass[5] = 1.0 / 24.0;
    bem->mass[6] = 1.0 / 24.0;
    bem->mass[7] = 1.0 / 24.0;
    bem->mass[8] = 1.0 / 12.0;
  }

  return bem;
}

field
eval_dirichlet_linear_laplacebem3d(const real * x, const real * n)
{
  (void) n;
  return (x[0] + x[1] + x[2]);
}

field
eval_neumann_linear_laplacebem3d(const real * x, const real * n)
{
  (void) x;
  return (n[0] + n[1] + n[2]);
}

field
eval_dirichlet_quadratic_laplacebem3d(const real * x, const real * n)
{
  (void) n;
  return (x[0] * x[0] - x[2] * x[2]);
}

field
eval_neumann_quadratic_laplacebem3d(const real * x, const real * n)
{
  return 2.0 * (n[0] * x[0] - n[2] * x[2]);
}

field
eval_dirichlet_fundamental_laplacebem3d(const real * x, const real * n)
{
  real      d[3];

  field     kernel;

  (void) n;

  d[0] = x[0] - (1.2);
  d[1] = x[1] - (1.2);
  d[2] = x[2] - (1.2);

  kernel = REAL_SQR(d[0]) + REAL_SQR(d[1]) + REAL_SQR(d[2]);
  kernel = 1.0 / REAL_SQRT(kernel);

  return kernel;
}

field
eval_neumann_fundamental_laplacebem3d(const real * x, const real * n)
{
  real      d[3];

  field     kernel;

  d[0] = x[0] - (1.2);
  d[1] = x[1] - (1.2);
  d[2] = x[2] - (1.2);

  kernel = REAL_SQR(d[0]) + REAL_SQR(d[1]) + REAL_SQR(d[2]);
  kernel = -(d[0] * n[0] + d[1] * n[1] + d[2] * n[2])
    / (kernel * REAL_SQRT(kernel));

  return kernel;
}

field
eval_dirichlet_fundamental2_laplacebem3d(const real * x, const real * n)
{
  real      d[3];

  field     kernel;

  (void) n;

  d[0] = x[0] - (0.25);
  d[1] = x[1] - (0.25);
  d[2] = x[2] - (0.25);

  kernel = REAL_SQR(d[0]) + REAL_SQR(d[1]) + REAL_SQR(d[2]);
  kernel = 1.0 / REAL_SQRT(kernel);

  return kernel;
}

field
eval_neumann_fundamental2_laplacebem3d(const real * x, const real * n)
{
  real      d[3];

  field     kernel;

  d[0] = x[0] - (0.25);
  d[1] = x[1] - (0.25);
  d[2] = x[2] - (0.25);

  kernel = REAL_SQR(d[0]) + REAL_SQR(d[1]) + REAL_SQR(d[2]);
  kernel = -(d[0] * n[0] + d[1] * n[1] + d[2] * n[2])
    / (kernel * REAL_SQRT(kernel));

  return kernel;
}

uint
build_interactive_laplacebem3d(pcsurface3d gr, char op,
			       basisfunctionbem3d basis_neumann,
			       basisfunctionbem3d basis_dirichlet, uint q,
			       void **G, real * time, char *filename)
{
  pbem3d    bem = NULL;
  admissible admiss = NULL;
  quadpoints3d quadpoints = NULL;
  pstopwatch sw = NULL;
  pclustergeometry neumann_cg = NULL, dirichlet_cg = NULL;
  uint     *neumann_idx = NULL, *dirichlet_idx = NULL;
  pcluster  neumann = NULL, dirichlet = NULL;
  pclusterbasis neumann_cb = NULL, dirichlet_cb = NULL;
  pblock    tree = NULL;
  FILE     *file;
  real      eta = 0.0, eps = 0.0, delta = 0.0, accur_recomp =
    0.0, accur_coarsen = 0.0, accur_hiercomp = 0.0;
  uint      N = 0, M, clf = 0, huhh2 = 0, method = 0, mi = 0, m = 0, l =
    0, type = 0;
  char      recomptech = 'n';
  bool      recomp = false, coarsen = false, hiercomp = false;
  field     alpha = 0.5;

  *G = NULL;

  sw = new_stopwatch();

  huhh2 = 1;
  method = 1;
  eta = 1.0;
  mi = 2;
  m = 2;
  l = 1;
  delta = 0.5;
  eps = 1.0e-4;
  quadpoints = build_bem3d_cube_quadpoints;
  clf = 0;
  recomp = false;
  accur_recomp = 1.0e-4;
  recomptech = 'n';
  coarsen = false;
  accur_coarsen = 1.0e-4;
  hiercomp = false;
  accur_hiercomp = 1.0e-4;
  admiss = admissible_max_cluster;

  type = 0;

  if (op == 'd') {
    alpha = askforreal("alpha = ?", "h2lib_alpha_bem3d", alpha);
  }

  huhh2 = askforint("Type of initial matrix-approximation:\n"
		    "  [1] h-matrix,\n"
		    "  [2] h2-matrix?\n", "h2lib_approxtype_bem3d", huhh2);

  switch (huhh2) {
  case 1:
    method = askforint("Method used for h-matrix-approximation:\n"
		       "  [ 1] Interpolation row-cluster\n"
		       "  [ 2] Interpolation col-cluster\n"
		       "  [ 3] Interpolation mixed\n"
		       "  [ 4] Green row-cluster\n"
		       "  [ 5] Green col-cluster\n"
		       "  [ 6] Green mixed\n"
		       "  [ 7] Greenhybrid row-cluster\n"
		       "  [ 8] Greenhybrid col-cluster\n"
		       "  [ 9] Greenhybrid mixed\n"
		       "  [10] ACA with full pivoting\n"
		       "  [11] ACA with partial pivoting\n"
		       "  [12] HCA\n", "h2lib_method_bem3d", method);

    switch (method) {
    case 1:			/* Interpolation row */
    case 2:			/* Interpolation col */
    case 3:			/* Interpolation mixed */
      mi = askforint("Interpolation order m?", "h2lib_mi_bem3d", mi);
      eta = 2.0;
      admiss = admissible_2_cluster;
      clf = 2 * mi * mi * mi;
      clf = askforint("Minimal leafsize?", "h2lib_clf_bem3d", clf);
      if (filename != NULL) {
	if (!(file = fopen(filename, "r"))) {
	  file = fopen(filename, "a+");
	  fprintf(file, "#m\tclf\tsize\ttime\tabs error\trel error\n");
	}
	fclose(file);
	file = fopen(filename, "a+");
	fprintf(file, "%d\t%d\t", mi, clf);
	fclose(file);
      }
      break;
    case 4:			/* Green row */
    case 5:			/* Green col */
    case 6:			/* Green mixed */
      delta = 0.5;
      quadpoints =
	askforchar("[c]ube parameterization?", "h2lib_greenparam_bem3d", "c",
		   'c')
	== 'c' ? build_bem3d_cube_quadpoints : NULL;
      admiss =
	(quadpoints == build_bem3d_cube_quadpoints) ?
	admissible_max_cluster : admissible_sphere_cluster;
      m = askforint("Quadrature order m?", "h2lib_m_bem3d", m);
      l = askforint("Number subdivisions l?", "h2lib_l_bem3d", l);
      delta = askforreal("delta = x * diam(t)?", "h2lib_delta_bem3d", delta);
      clf = 2 * 12 * m * m * l * l;
      clf = askforint("Minimal leafsize?", "h2lib_clf_bem3d", clf);
      if (filename != NULL) {
	if (!(file = fopen(filename, "r"))) {
	  file = fopen(filename, "a+");
	  fprintf(file,
		  "#m\tl\tdelta\tparam\tclf\tsize\ttime\tabs error\trel error\n");
	}
	fclose(file);
	file = fopen(filename, "a+");
	fprintf(file, "%d\t%d\t%.2f\t%c\t%d\t", m, l, delta,
		quadpoints == build_bem3d_cube_quadpoints ? 'c' : 's', clf);
	fclose(file);
      }
      break;
    case 7:			/* Greenhybrid row */
    case 8:			/* Greenhybrid col */
    case 9:			/* Greenhybrid mixed */
      delta = 1.0;
      quadpoints =
	askforchar("[c]ube parameterization?", "h2lib_greenparam_bem3d", "c",
		   'c')
	== 'c' ? build_bem3d_cube_quadpoints : NULL;
      admiss =
	(quadpoints == build_bem3d_cube_quadpoints) ?
	admissible_max_cluster : admissible_sphere_cluster;
      m = askforint("Quadrature order m?", "h2lib_m_bem3d", m);
      l = askforint("Number subdivisions l?", "h2lib_l_bem3d", l);
      eps = askforreal("ACA accuracy?", "h2lib_eps_bem3d", eps);
      delta = askforreal("delta = x * diam(t)?", "h2lib_delta_bem3d", delta);
      clf = REAL_LOG(gr->triangles);
      clf = askforint("Minimal leafsize?", "h2lib_clf_bem3d", clf);
      if (filename != NULL) {
	if (!(file = fopen(filename, "r"))) {
	  file = fopen(filename, "a+");
	  fprintf(file,
		  "#m\tl\tdelta\teps\tparam\tclf\tsize\ttime\tabs error\trel error\n");
	}
	fclose(file);
	file = fopen(filename, "a+");
	fprintf(file, "%d\t%d\t%.2f\t%.2e\t%c\t%d\t", m, l, delta, eps,
		quadpoints == build_bem3d_cube_quadpoints ? 'c' : 's', clf);
	fclose(file);
      }
      break;
    case 10:			/* ACA full pivoting */
    case 11:			/* ACA partial pivoting */
      admiss = admissible_max_cluster;
      eps = askforreal("ACA accuracy?", "h2lib_eps_bem3d", eps);
      clf = REAL_LOG(gr->triangles);
      clf = askforint("Minimal leafsize?", "h2lib_clf_bem3d", clf);
      if (filename != NULL) {
	if (!(file = fopen(filename, "r"))) {
	  file = fopen(filename, "a+");
	  fprintf(file, "#eps\tclf\tsize\ttime\tabs error\trel error\n");
	}
	fclose(file);
	file = fopen(filename, "a+");
	fprintf(file, "%.2e\t%d\t", eps, clf);
	fclose(file);
      }
      break;
    case 12:			/* HCA with interpolation */
      mi = askforint("Interpolation order m?", "h2lib_mi_bem3d", mi);
      eps = askforreal("ACA accuracy?", "h2lib_eps_bem3d", eps);
      eta = 2.0;
      admiss = admissible_2_cluster;
      clf = REAL_LOG(gr->triangles);
      clf = askforint("Minimal leafsize?", "h2lib_clf_bem3d", clf);
      if (filename != NULL) {
	if (!(file = fopen(filename, "r"))) {
	  file = fopen(filename, "a+");
	  fprintf(file, "#m\teps\tclf\tsize\ttime\tabs error\trel error\n");
	}
	fclose(file);
	file = fopen(filename, "a+");
	fprintf(file, "%d\t%.2e\t%d\t", mi, eps, clf);
	fclose(file);
      }
      break;
    default:
      break;
    }

    eta = askforreal("eta?", "h2lib_eta_bem3d", eta);

    break;
  case 2:
    method = askforint("Method used for h2-matrix-approximation:\n"
		       "  [1] Interpolation\n"
		       "  [2] Greenhybrid\n"
		       "  [3] Greenhybrid ortho\n", "h2lib_method_bem3d",
		       method);

    switch (method) {
    case 1:			/* Interpolation */
      admiss = admissible_2_cluster;
      eta = 2.0;
      m = askforint("Interpolation order m?", "h2lib_m_bem3d", mi);
      clf = 2 * mi * mi * mi;
      clf = askforint("Minimal leafsize?", "h2lib_clf_bem3d", clf);
      if (filename != NULL) {
	if (!(file = fopen(filename, "r"))) {
	  file = fopen(filename, "a+");
	  fprintf(file, "#m\tclf\tsize\ttime\tabs error\trel error\n");
	}
	fclose(file);
	file = fopen(filename, "a+");
	fprintf(file, "%d\t%d\t", mi, clf);
	fclose(file);
      }
      break;
    case 2:			/* Greenhybrid */
      delta = 1.0;
      quadpoints =
	askforchar("[c]ube parameterization?", "h2lib_greenparam_bem3d", "c",
		   'c')
	== 'c' ? build_bem3d_cube_quadpoints : NULL;
      admiss =
	(quadpoints == build_bem3d_cube_quadpoints) ?
	admissible_max_cluster : admissible_sphere_cluster;
      m = askforint("Quadrature order m?", "h2lib_m_bem3d", m);
      l = askforint("Number subdivisions l?", "h2lib_l_bem3d", l);
      eps = askforreal("ACA accuracy?", "h2lib_eps_bem3d", eps);
      delta = askforreal("delta = x * diam(t)?", "h2lib_delta_bem3d", delta);
      clf = REAL_LOG(gr->triangles);
      clf = askforint("Minimal leafsize?", "h2lib_clf_bem3d", clf);
      if (filename != NULL) {
	if (!(file = fopen(filename, "r"))) {
	  file = fopen(filename, "a+");
	  fprintf(file,
		  "#m\tl\tdelta\teps\tparam\tclf\tsize\ttime\tabs error\trel error\n");
	}
	fclose(file);
	file = fopen(filename, "a+");
	fprintf(file, "%d\t%d\t%.2f\t%.2e\t%c\t%d\t", m, l, delta, eps,
		quadpoints == build_bem3d_cube_quadpoints ? 'c' : 's', clf);
	fclose(file);
      }
      break;
    case 3:			/* Greenhybrid ortho */
      delta = 1.0;
      quadpoints =
	askforchar("[c]ube parameterization?", "h2lib_greenparam_bem3d", "c",
		   'c')
	== 'c' ? build_bem3d_cube_quadpoints : NULL;
      admiss =
	(quadpoints == build_bem3d_cube_quadpoints) ?
	admissible_max_cluster : admissible_sphere_cluster;
      m = askforint("Quadrature order m?", "h2lib_m_bem3d", m);
      l = askforint("Number subdivisions l?", "h2lib_l_bem3d", l);
      eps = askforreal("ACA accuracy?", "h2lib_eps_bem3d", eps);
      delta = askforreal("delta = x * diam(t)?", "h2lib_delta_bem3d", delta);
      clf = REAL_LOG(gr->triangles);
      clf = askforint("Minimal leafsize?", "h2lib_clf_bem3d", clf);
      if (filename != NULL) {
	if (!(file = fopen(filename, "r"))) {
	  file = fopen(filename, "a+");
	  fprintf(file,
		  "#m\tl\tdelta\teps\tparam\tclf\tsize\ttime\tabs error\trel error\n");
	}
	fclose(file);
	file = fopen(filename, "a+");
	fprintf(file, "%d\t%d\t%.2f\t%.2e\t%c\t%d\t", m, l, delta, eps,
		quadpoints == build_bem3d_cube_quadpoints ? 'c' : 's', clf);
	fclose(file);
      }
      break;
    default:
      break;
    }

    eta = askforreal("eta?", "h2lib_eta_bem3d", eta);

    break;
  default:
    printf("ERROR: unknown Matrix-type!\n");
    abort();
    break;
  }

  if (askforchar("Recompression?", "h2lib_recomp_bem3d", "yn",
		 recomp ? 'y' : 'n')
      == 'y') {
    accur_recomp = askforreal("Recompression accuracy?",
			      "h2lib_accur_recomp_bem3d", accur_recomp);
    recomp = true;
  }
  else {
    recomp = false;
  }

  if (huhh2 == 1) {
    recomptech = askforchar("Further recompression techniques ?\n"
			    "  [c]oarsen,\n"
			    "  [h]ierarchical compression,\n"
			    "  [n]one?\n", "h2lib_recomptech_bem3d", "chn",
			    recomptech);
    if (recomptech == 'c') {
      coarsen = true;
      accur_coarsen = askforreal("Coarsen accuracy?",
				 "h2lib_accur_coarsen_bem3d", accur_coarsen);
      hiercomp = false;
      accur_hiercomp = 0.0;
    }
    else if (recomptech == 'h') {
      coarsen = false;
      accur_coarsen = 0.0;
      hiercomp = true;
      accur_hiercomp = askforreal("Hierarchical compression accuracy?",
				  "h2lib_accur_hiercomp_bem3d",
				  accur_hiercomp);
    }
    else {
      coarsen = false;
      accur_coarsen = 0.0;
      hiercomp = false;
      accur_hiercomp = 0.0;
    }
  }

  if (op == 's') {
    bem = new_slp_laplace_bem3d(gr, q, basis_neumann);
  }
  else {
    assert(op == 'd');
    bem = new_dlp_laplace_bem3d(gr, q, basis_neumann, basis_dirichlet, alpha);
  }

  N = bem->N_neumann;
  M = bem->N_dirichlet;

  start_stopwatch(sw);

  neumann_cg = build_bem3d_clustergeometry(bem, &neumann_idx,
					   bem->basis_neumann);
  neumann = build_adaptive_cluster(neumann_cg, N, neumann_idx, clf);
  dirichlet = neumann;
  del_clustergeometry(neumann_cg);

  if (op == 'd') {
    dirichlet_cg = build_bem3d_clustergeometry(bem, &dirichlet_idx,
					       bem->basis_dirichlet);
    dirichlet = build_adaptive_cluster(dirichlet_cg, M, dirichlet_idx, clf);
    del_clustergeometry(dirichlet_cg);
  }

  switch (huhh2) {
  case 1:
    tree = build_nonstrict_block(neumann, dirichlet, (void *) &eta, admiss);

    /* h-matrix approximations */
    switch (method) {
    case 1:			/* Interpolation row */
      setup_hmatrix_aprx_inter_row_bem3d(bem, neumann, dirichlet, tree, mi);
      setup_hmatrix_recomp_bem3d(bem, recomp, accur_recomp, coarsen,
				 accur_coarsen);
      break;
    case 2:			/* Interpolation col */
      setup_hmatrix_aprx_inter_col_bem3d(bem, neumann, dirichlet, tree, mi);
      setup_hmatrix_recomp_bem3d(bem, recomp, accur_recomp, coarsen,
				 accur_coarsen);
      break;
    case 3:			/* Interpolation mixed */
      setup_hmatrix_aprx_inter_mixed_bem3d(bem, neumann, dirichlet, tree, mi);
      setup_hmatrix_recomp_bem3d(bem, recomp, accur_recomp, coarsen,
				 accur_coarsen);
      break;
    case 4:			/* Green row */
      setup_hmatrix_aprx_green_row_bem3d(bem, neumann, dirichlet, tree, m, l,
					 delta, quadpoints);
      setup_hmatrix_recomp_bem3d(bem, recomp, accur_recomp, coarsen,
				 accur_coarsen);
      break;
    case 5:			/* Green col */
      setup_hmatrix_aprx_green_col_bem3d(bem, neumann, dirichlet, tree, m, l,
					 delta, quadpoints);
      setup_hmatrix_recomp_bem3d(bem, recomp, accur_recomp, coarsen,
				 accur_coarsen);
      break;
    case 6:			/* Green mixed */
      setup_hmatrix_aprx_green_mixed_bem3d(bem, neumann, dirichlet, tree, m,
					   l, delta, quadpoints);
      setup_hmatrix_recomp_bem3d(bem, recomp, accur_recomp, coarsen,
				 accur_coarsen);
      break;
    case 7:			/* Greenhybrid row */
      setup_hmatrix_aprx_greenhybrid_row_bem3d(bem, neumann, dirichlet, tree,
					       m, l, delta, eps, quadpoints);
      setup_hmatrix_recomp_bem3d(bem, recomp, accur_recomp, coarsen,
				 accur_coarsen);
      break;
    case 8:			/* Greenhybrid col */
      setup_hmatrix_aprx_greenhybrid_col_bem3d(bem, neumann, dirichlet, tree,
					       m, l, delta, eps, quadpoints);
      setup_hmatrix_recomp_bem3d(bem, recomp, accur_recomp, coarsen,
				 accur_coarsen);
      break;
    case 9:			/* Greenhybrid mixed */
      setup_hmatrix_aprx_greenhybrid_mixed_bem3d(bem, neumann, dirichlet,
						 tree, m, l, delta, eps,
						 quadpoints);
      setup_hmatrix_recomp_bem3d(bem, recomp, accur_recomp, coarsen,
				 accur_coarsen);
      break;
    case 10:			/* ACA full pivoting */
      setup_hmatrix_aprx_aca_bem3d(bem, neumann, dirichlet, tree, eps);
      setup_hmatrix_recomp_bem3d(bem, recomp, accur_recomp, coarsen,
				 accur_coarsen);
      break;
    case 11:			/* ACA partial pivoting */
      setup_hmatrix_aprx_paca_bem3d(bem, neumann, dirichlet, tree, eps);
      setup_hmatrix_recomp_bem3d(bem, recomp, accur_recomp, coarsen,
				 accur_coarsen);
      break;
    case 12:			/* HCA with interpolation */
      setup_hmatrix_aprx_hca_bem3d(bem, neumann, dirichlet, tree, mi, eps);
      setup_hmatrix_recomp_bem3d(bem, recomp, accur_recomp, coarsen,
				 accur_coarsen);
      break;
    default:
      break;
    }

    if (coarsen) {
      *G = build_from_block_hmatrix(tree, 0);
      type = 1;
      assemblecoarsen_bem3d_hmatrix(bem, tree, (phmatrix) *G);
    }
    else if (hiercomp) {
      del_block(tree);
      tree = build_strict_block(neumann, dirichlet, (void *) &eta, admiss);
      neumann_cb = build_from_cluster_clusterbasis(neumann);
      dirichlet_cb = neumann_cb;
      if (op == 'd') {
	dirichlet_cb = build_from_cluster_clusterbasis(dirichlet);
      }
      *G = build_from_block_h2matrix(tree, neumann_cb, dirichlet_cb);
      type = 3;
      setup_h2matrix_recomp_bem3d(bem, hiercomp, accur_hiercomp);
      assemblehiercomp_bem3d_h2matrix(bem, tree, (ph2matrix) *G);
    }
    else {
      *G = build_from_block_hmatrix(tree, 0);
      type = 1;
      assemble_bem3d_hmatrix(bem, tree, (phmatrix) *G);
    }

    break;
  case 2:
    type = 2;

    neumann_cb = build_from_cluster_clusterbasis(neumann);
    dirichlet_cb = build_from_cluster_clusterbasis(dirichlet);

    tree = build_strict_block(neumann, dirichlet, (void *) &eta, admiss);
    *G = build_from_block_h2matrix(tree, neumann_cb, dirichlet_cb);

    /* H2-matrix approximations */
    switch (method) {
    case 1:			/* Interpolation H2 */
      setup_h2matrix_aprx_inter_bem3d(bem, neumann_cb, dirichlet_cb, tree,
				      mi);
      break;
    case 2:			/* Greenhybrid H2 */
      setup_h2matrix_aprx_greenhybrid_bem3d(bem, neumann_cb, dirichlet_cb,
					    tree, m, l, delta, eps,
					    quadpoints);
      break;
    case 3:			/* Greenhybrid ortho H2 */
      setup_h2matrix_aprx_greenhybrid_ortho_bem3d(bem, neumann_cb,
						  dirichlet_cb, tree, m, l,
						  delta, eps, quadpoints);
      break;
    default:
      break;
    }

    assemble_bem3d_h2matrix_row_clusterbasis(bem, neumann_cb);
    assemble_bem3d_h2matrix_col_clusterbasis(bem, dirichlet_cb);

    assemble_bem3d_h2matrix(bem, tree, (ph2matrix) *G);

    if (recomp == true) {
      recompress_inplace_h2matrix((ph2matrix) *G, 0, accur_recomp);
    }
    break;
  }

  *time = stop_stopwatch(sw);

  del_stopwatch(sw);
  del_bem3d(bem);
  del_block(tree);

  return type;
}
