/* ------------------------------------------------------------
 This is the file "singquad2d.c" of the H2Lib package.
 All rights reserved, Sven Christophersen 2011-2013
 ------------------------------------------------------------ */

/* C STD LIBRARY */
#include <stdio.h>

/* CORE 0 */
#include "basic.h"

/* CORE 1 */

/* CORE 2 */

/* CORE 3 */

/* SIMPLE */

/* PARTICLES */

/* BEM */
#include "singquad2d.h"

static void build_triangle_singquad2d(psingquad2d sq, real * xq, real * wq) {
  uint nq = sq->n_single;
  uint vnq = ROUNDUP(nq, VREAL);
  real *xx = sq->x_single;
  real *yy = sq->y_single;
  real *ww = sq->w_single + 3 * vnq;
  uint q = sq->q;

  uint i, j, p;

  p = 0;

  for (i = 0; i < q; ++i) {
    for (j = 0; j < q; ++j) {
      ww[p] = wq[i] * wq[j] * xq[i];
      xx[p] = xq[i];
      yy[p] = xq[i] * xq[j];
      p++;
    }
  }

  assert(p == nq);

}

static void build_pow_dist_singquad2d(psingquad2d sq, real * xq, real * wq) {
  uint nq = sq->n_dist;
  uint vnq = ROUNDUP(nq, VREAL);
  real *xx = sq->x_dist;
  real *yy = sq->y_dist;
  real *ww = sq->w_dist + 9 * vnq;
  uint q = sq->q;

  uint i, j, k, l, p;

  p = 0;
  
   for (i = 0; i < q; ++i) {
     for (j = 0; j < q; ++j) {
       for (k = 0; k < q; ++k) {
         for (l = 0; l < q; ++l) {
           ww[p] = wq[i] * wq[j] * wq[k] * wq[l] * xq[i] * xq[k];
           xx[p] = xq[i];
           xx[p + vnq] = xq[i] * xq[j];
           yy[p] = xq[k];
           yy[p + vnq] = xq[k] * xq[l];
           p++;
         }
       }
     }
   }

  assert(p == nq);
}

static void build_pow_vert_singquad2d(psingquad2d sq, real * xq, real * wq) {
  uint nq = sq->n_vert;
  uint vnq = ROUNDUP(nq, VREAL);
  real *xx = sq->x_vert;
  real *yy = sq->y_vert;
  real *ww = sq->w_vert + 9 * vnq;
  uint q = sq->q2;

  uint i, j, k, l, p;
  real eta1, eta2, eta3, xi, wi, wj, wk, wl;

  p = 0;

  for (i = 0; i < q; ++i) {
    eta1 = xq[i];
    wi = wq[i];
    for (j = 0; j < q; ++j) {
      eta2 = xq[j];
      wj = wq[j] * wi;
      for (k = 0; k < q; ++k) {
        eta3 = xq[k];
        wk = wq[k] * wj;
        for (l = 0; l < q; ++l) {
          xi = xq[l];
          wl = wq[l] * wk;

          ww[p] = wl * xi * xi * xi * eta1;
          xx[p] = xi;
          xx[p + vnq] = eta2 * xi;
          yy[p] = -(-eta1 * xi);
          yy[p + vnq] = -(-eta1 * eta3 * xi);
          p++;

          ww[p] = wl * xi * xi * xi * eta1;
          xx[p] = eta1 * xi;
          xx[p + vnq] = eta1 * eta3 * xi;
          yy[p] = -(-xi);
          yy[p + vnq] = -(-eta2 * xi);
          p++;

        }
      }
    }
  }

  assert(p == nq);

}

static void build_pow_edge_singquad2d(psingquad2d sq, real * xq, real * wq) {
  uint nq = sq->n_edge;
  uint vnq = ROUNDUP(nq, VREAL);
  real *xx = sq->x_edge;
  real *yy = sq->y_edge;
  real *ww = sq->w_edge + 9 * vnq;
  uint q = sq->q2;

  uint i, j, k, l, p;
  real eta1, eta2, xi1, xi2, wi, wj, wk, wl;

  p = 0;

  for (i = 0; i < q; ++i) {
    eta1 = xq[i];
    wi = wq[i];
    for (j = 0; j < q; ++j) {
      eta2 = xq[j];
      wj = wq[j] * wi;
      for (k = 0; k < q; ++k) {
        xi1 = xq[k];
        wk = wq[k] * wj;
        for (l = 0; l < q; ++l) {
          xi2 = xq[l] * (1.0 - xi1);
          wl = wq[l] * wk * (1.0 - xi1);

          ww[p] = wl * xi1 * xi1 * eta1;
          xx[p] = xi1 + xi2;
          xx[p + vnq] = xi1 * eta1 * eta2;
          yy[p] = xi1 * (1.0 - eta1) + xi2;
          yy[p + vnq] = -(xi1 * (eta1 - 1.0));
          p++;

          ww[p] = wl * xi1 * xi1 * eta1;
          xx[p] = xi1 + xi2;
          xx[p + vnq] = xi1 * eta1;
          yy[p] = xi1 * (1.0 - eta1 * eta2) + xi2;
          yy[p + vnq] = -(xi1 * (eta1 * eta2 - 1.0));
          p++;

          ww[p] = wl * xi1 * xi1 * eta1;
          xx[p] = xi1 + xi2;
          xx[p + vnq] = xi1;
          yy[p] = xi1 * eta1 + xi2;
          yy[p + vnq] = xi1 * eta1 * eta2;
          p++;

          ww[p] = wl * xi1 * xi1 * eta1;
          xx[p] = xi1 * (1.0 - eta1 + eta1 * eta2) + xi2;
          xx[p + vnq] = xi1 * eta1 * eta2;
          yy[p] = xi1 + xi2;
          yy[p + vnq] = xi1;
          p++;

          ww[p] = wl * xi1 * xi1 * eta1;
          xx[p] = xi1 * eta1 * eta2 + xi2;
          xx[p + vnq] = xi1 * eta1 * eta2;
          yy[p] = xi1 + (xi2);
          yy[p + vnq] = xi1 * eta1;
          p++;

          ww[p] = wl * xi1 * xi1 * eta1;
          xx[p] = xi1 * eta1 + xi2;
          xx[p + vnq] = xi1 * eta1;
          yy[p] = xi1 + xi2;
          yy[p + vnq] = xi1 * eta1 * eta2;
          p++;

        }
      }
    }
  }

  assert(p == nq);

}

static void build_pow_id_singquad2d(psingquad2d sq, real * xq, real * wq) {
  uint nq = sq->n_id;
  uint vnq = ROUNDUP(nq, VREAL);
  real *xx = sq->x_id;
  real *yy = sq->y_id;
  real *ww = sq->w_id + 9 * vnq;
  uint q = sq->q2;

  uint i, j, k, l, p;
  real eta, xi1, xi2, xi3, wi, wj, wk, wl;

  p = 0;

  for (i = 0; i < q; ++i) {
    eta = xq[i];
    wi = wq[i];
    for (j = 0; j < q; ++j) {
      xi1 = xq[j];
      wj = wq[j] * wi;
      for (k = 0; k < q; ++k) {
        xi2 = xq[k] * (1.0 - xi1);
        wk = wq[k] * wj * (1.0 - xi1);
        for (l = 0; l < q; ++l) {
          xi3 = xq[l] * (1.0 - xi1 - xi2);
          wl = wq[l] * wk * (1.0 - xi1 - xi2);

          ww[p] = wl * xi1;
          xx[p] = xi1 + xi2 + xi3;
          xx[p + vnq] = xi1 + xi2;
          yy[p] = xi1 * (1.0 - eta) + xi2 + xi3;
          yy[p + vnq] = xi2;
          p++;

          ww[p] = wl * xi1;
          xx[p] = xi1 * (1.0 - eta) + xi2 + xi3;
          xx[p + vnq] = xi1 * (1.0 - eta) + xi2;
          yy[p] = xi1 + xi2 + xi3;
          yy[p + vnq] = xi2;
          p++;

          ww[p] = wl * xi1;
          xx[p] = xi1 + xi2 + xi3;
          xx[p + vnq] = xi1 * eta + xi2;
          yy[p] = xi2 + xi3;
          yy[p + vnq] = xi2;
          p++;

          ww[p] = wl * xi1;
          xx[p] = xi1 * (1.0 - eta) + xi2 + xi3;
          xx[p + vnq] = xi2;
          yy[p] = xi1 + xi2 + xi3;
          yy[p + vnq] = xi1 + xi2;
          p++;

          ww[p] = wl * xi1;
          xx[p] = xi1 + xi2 + xi3;
          xx[p + vnq] = xi2;
          yy[p] = xi1 * (1.0 - eta) + xi2 + xi3;
          yy[p + vnq] = xi1 * (1.0 - eta) + xi2;
          p++;

          ww[p] = wl * xi1;
          xx[p] = xi2 + xi3;
          xx[p + vnq] = xi2;
          yy[p] = xi1 + xi2 + xi3;
          yy[p + vnq] = xi1 * eta + xi2;
          p++;

        }
      }
    }
  }

  assert(p == nq);

}

#ifdef USE_TRIQUADPOINTS
static void init_triquadpoints(pcsurface3d gr, psingquad2d sq) {
  uint n = gr->triangles;
  uint(*gr_t)[3] = gr->t;
  real(*gr_x)[3] = gr->x;
  uint q2 = sq->q * sq->q;
  uint q2ld = ROUNDUP(q2, VREAL);
  real *xx = sq->x_single;
  real *yy = sq->y_single;

  real *tri_x, *tri_y, *tri_z, *A, *B, *C;
  real tx, sx, Ax, Bx, Cx;
  uint q, t;

  sq->tri_x = tri_x = allocreal(n * q2ld);
  sq->tri_y = tri_y = allocreal(n * q2ld);
  sq->tri_z = tri_z = allocreal(n * q2ld);

  for (t = 0; t < n; ++t) {
    A = gr_x[gr_t[t][0]];
    B = gr_x[gr_t[t][1]];
    C = gr_x[gr_t[t][2]];
    for (q = 0; q < q2; ++q) {
      tx = xx[q];
      sx = yy[q];
      Ax = 1.0 - tx;
      Bx = tx - sx;
      Cx = sx;

      tri_x[q + t * q2ld] = A[0] * Ax + B[0] * Bx + C[0] * Cx;
      tri_y[q + t * q2ld] = A[1] * Ax + B[1] * Bx + C[1] * Cx;
      tri_z[q + t * q2ld] = A[2] * Ax + B[2] * Bx + C[2] * Cx;
    }
  }
}
#endif

psingquad2d build_singquad2d(pcsurface3d gr, uint q, uint q2) {
  uint i, nq, nq2;
  real *x, *w, *x2, *w2;

  psingquad2d sq;

  x = allocreal(q);
  w = allocreal(q);
  x2 = allocreal(q2);
  w2 = allocreal(q2);

  assemble_gauss(q, x, w);
  for (i = 0; i < q; ++i) {
    x[i] = 0.5 + 0.5 * x[i];
    w[i] = w[i] * 0.5;
  }
  nq = q * q * q * q;

  assemble_gauss(q2, x2, w2);
  for (i = 0; i < q2; ++i) {
    x2[i] = 0.5 + 0.5 * x2[i];
    w2[i] = w2[i] * 0.5;
  }
  nq2 = q2 * q2 * q2 * q2;

  sq = allocmem((size_t ) sizeof(singquad2d));

  sq->q2 = q2;

  sq->n_id = 6 * nq2;
  sq->x_id = allocreal((size_t) 2 * ROUNDUP(sq->n_id, VREAL));
  sq->y_id = allocreal((size_t) 2 * ROUNDUP(sq->n_id, VREAL));
  sq->w_id = allocreal((size_t) 10 * ROUNDUP(sq->n_id, VREAL));
  sq->base_id = 0.0;

  build_pow_id_singquad2d(sq, x2, w2);

  sq->n_edge = 6 * nq2;
  sq->x_edge = allocreal((size_t) 2 * ROUNDUP(sq->n_edge, VREAL));
  sq->y_edge = allocreal((size_t) 2 * ROUNDUP(sq->n_edge, VREAL));
  sq->w_edge = allocreal((size_t) 10 * ROUNDUP(sq->n_edge, VREAL));
  sq->base_edge = 0.0;

  build_pow_edge_singquad2d(sq, x2, w2);

  sq->n_vert = 2 * nq2;
  sq->x_vert = allocreal((size_t) 2 * ROUNDUP(sq->n_vert, VREAL));
  sq->y_vert = allocreal((size_t) 2 * ROUNDUP(sq->n_vert, VREAL));
  sq->w_vert = allocreal((size_t) 10 * ROUNDUP(sq->n_vert, VREAL));
  sq->base_vert = 0.0;

  build_pow_vert_singquad2d(sq, x2, w2);

  sq->q = q;

  sq->n_dist = nq;
  sq->x_dist = allocreal((size_t) 2 * ROUNDUP(sq->n_dist, VREAL));
  sq->y_dist = allocreal((size_t) 2 * ROUNDUP(sq->n_dist, VREAL));
  sq->w_dist = allocreal((size_t) 10 * ROUNDUP(sq->n_dist, VREAL));
  sq->base_dist = 0.0;

  build_pow_dist_singquad2d(sq, x, w);

  sq->n_single = q * q;
  sq->x_single = allocreal((size_t) ROUNDUP(sq->n_single, VREAL));
  sq->y_single = allocreal((size_t) ROUNDUP(sq->n_single, VREAL));
  sq->w_single = allocreal((size_t) 4 * ROUNDUP(sq->n_single, VREAL));
  sq->base_single = 0.0;

  build_triangle_singquad2d(sq, x, w);

  sq->nmax = 6 * nq2;

#ifdef USE_TRIQUADPOINTS
  init_triquadpoints(gr, sq);
#endif

  freemem(x);
  freemem(w);
  freemem(x2);
  freemem(w2);

  return sq;
}

void del_singquad2d(psingquad2d sq) {
  assert(sq != NULL);

  if (sq->w_dist != NULL)
    freemem(sq->w_dist);
  if (sq->w_vert != NULL)
    freemem(sq->w_vert);
  if (sq->w_edge != NULL)
    freemem(sq->w_edge);
  if (sq->w_id != NULL)
    freemem(sq->w_id);
  if (sq->w_single != NULL)
    freemem(sq->w_single);

  if (sq->x_dist != NULL)
    freemem(sq->x_dist);
  if (sq->x_vert != NULL)
    freemem(sq->x_vert);
  if (sq->x_edge != NULL)
    freemem(sq->x_edge);
  if (sq->x_id != NULL)
    freemem(sq->x_id);
  if (sq->x_single != NULL)
    freemem(sq->x_single);

  if (sq->y_dist != NULL)
    freemem(sq->y_dist);
  if (sq->y_vert != NULL)
    freemem(sq->y_vert);
  if (sq->y_edge != NULL)
    freemem(sq->y_edge);
  if (sq->y_id != NULL)
    freemem(sq->y_id);
  if (sq->y_single != NULL)
    freemem(sq->y_single);

#ifdef USE_TRIQUADPOINTS
  if (sq->tri_x != NULL) {
    freemem(sq->tri_x);
  }
  if (sq->tri_y != NULL) {
    freemem(sq->tri_y);
  }
  if (sq->tri_z != NULL) {
    freemem(sq->tri_z);
  }
#endif

  freemem(sq);
}

void weight_basisfunc_ll_singquad2d(real * x, real * y, real * w, uint nq) {
  uint vnq = ROUNDUP(nq, VREAL);

  uint i, j, k, idx;
  real bx, by;

  idx = 0;
  for (j = 0; j < 3; ++j) {
    for (i = 0; i < 3; ++i) {
      for (k = 0; k < nq; ++k) {
        bx =
            (i == 0 ?
                1.0 - x[k + 0 * vnq] :
                (i == 1 ? x[k + 0 * vnq] - x[k + 1 * vnq] : x[k + 1 * vnq]));
        by =
            (j == 0 ?
                1.0 - y[k + 0 * vnq] :
                (j == 1 ? y[k + 0 * vnq] - y[k + 1 * vnq] : y[k + 1 * vnq]));
        assert(idx < 9);
        w[k + idx * vnq] = w[k + 9 * vnq] * bx * by;
      }
      idx++;
    }
  }
}

void weight_basisfunc_cl_singquad2d(real * x, real * y, real * w, uint nq) {
  uint vnq = ROUNDUP(nq, VREAL);

  uint i, k;
  real by;

  (void) x;

  for (i = 0; i < 3; ++i) {
    for (k = 0; k < nq; ++k) {
      by = (
          i == 0 ? 1.0 - y[k + 0 * vnq] :
                   (i == 1 ? y[k + 0 * vnq] - y[k + 1 * vnq] : y[k + 1 * vnq]));
      w[k + i * vnq] = w[k + 9 * vnq] * by;
    }
  }
}

void weight_basisfunc_lc_singquad2d(real * x, real * y, real * w, uint nq) {
  uint vnq = ROUNDUP(nq, VREAL);

  uint i, k;
  real by;

  (void) y;

  for (i = 0; i < 3; ++i) {
    for (k = 0; k < nq; ++k) {
      by = (
          i == 0 ? 1.0 - x[k + 0 * vnq] :
                   (i == 1 ? x[k + 0 * vnq] - x[k + 1 * vnq] : x[k + 1 * vnq]));
      w[k + i * vnq] = w[k + 9 * vnq] * by;
    }
  }
}

void weight_basisfunc_l_singquad2d(real * x, real * y, real * w, uint nq) {
  uint vnq = ROUNDUP(nq, VREAL);

  uint i, k;
  real bx;

  for (i = 0; i < 3; ++i) {
    for (k = 0; k < nq; ++k) {
      bx = (i == 0 ? 1.0 - x[k] : (i == 1 ? x[k] - y[k] : y[k]));
      w[k + i * vnq] = w[k + 3 * vnq] * bx;
    }
  }
}

uint fast_select_quadrature(uint (*geo_t)[3], uint t, uint s) {
  uint p;

  p = 0;
  p += (geo_t[t][0] == geo_t[s][0]);
  p += (geo_t[t][0] == geo_t[s][1]);
  p += (geo_t[t][0] == geo_t[s][2]);
  p += (geo_t[t][1] == geo_t[s][0]);
  p += (geo_t[t][1] == geo_t[s][1]);
  p += (geo_t[t][1] == geo_t[s][2]);
  p += (geo_t[t][2] == geo_t[s][0]);
  p += (geo_t[t][2] == geo_t[s][1]);
  p += (geo_t[t][2] == geo_t[s][2]);

  return p;
}

uint select_quadrature_singquad2d(pcsingquad2d sq, const uint * tv,
    const uint * sv, uint * tp, uint * sp, real ** x, real ** y, real ** w,
    uint * n, real * base) {

  uint p, q, i, j;

  p = (tv[0] == sv[0]) + (tv[0] == sv[1]) + (tv[0] == sv[2]) + (tv[1] == sv[0])
      + (tv[1] == sv[1]) + (tv[1] == sv[2]) + (tv[2] == sv[0])
      + (tv[2] == sv[1]) + (tv[2] == sv[2]);

  tp[0] = 0, tp[1] = 1, tp[2] = 2;
  sp[0] = 0, sp[1] = 1, sp[2] = 2;

  switch (p) {
  case 0: /* DISTANT */
    *x = sq->x_dist;
    *y = sq->y_dist;
    *w = sq->w_dist;
    *n = sq->n_dist;
    *base = sq->base_dist;
    return p;
    break;
  case 1: /* VERTEX */
    *x = sq->x_vert;
    *y = sq->y_vert;
    *w = sq->w_vert;
    *n = sq->n_vert;
    *base = sq->base_vert;
    break;
  case 2: /* EDGE */
    *x = sq->x_edge;
    *y = sq->y_edge;
    *w = sq->w_edge;
    *n = sq->n_edge;
    *base = sq->base_edge;
    break;
  case 3: /* IDENTICAL */
    *x = sq->x_id;
    *y = sq->y_id;
    *w = sq->w_id;
    *n = sq->n_id;
    *base = sq->base_id;
    return p;
    break;
  default:
    printf("ERROR: Unknown quadrature situation!\n");
    abort();
    break;
  }

  p = 0;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      if (tv[i] == sv[j]) {
        tp[p] = i;
        sp[p] = j;
        p++;
        break;
      }
    }
  }

  q = p;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < q && tv[i] != tv[tp[j]]; j++)
      ;
    if (j == q)
      tp[q++] = i;
  }
  assert(q == 3);

  q = p;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < q && sv[i] != sv[sp[j]]; j++)
      ;
    if (j == q)
      sp[q++] = i;
  }
  assert(q == 3);

  assert(p <= 3);

  return p;

}
