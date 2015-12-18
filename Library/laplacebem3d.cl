#include "Library/clsettings.h"

#ifdef USE_FLOAT
#define KERNEL_CONST_LAPLACEBEM3D 7.95774715459477e-02f
#define KERNEL_CONST_2 7.95774715459477e-02f
#define KERNEL_CONST_8 7.95774715459477e-02f
#else
#define KERNEL_CONST_LAPLACEBEM3D 7.95774715459477e-02
#define KERNEL_CONST_2 3.978873577297385e-02
#define KERNEL_CONST_8 9.9471839432434625e-03
#endif

inline field slp_eval(real dx, real dy, real dz) {
  real norm, norm2;

  norm2 = dx * dx + dy * dy + dz * dz;

#ifndef USE_FLOAT
  norm = convert_double(native_rsqrt(convert_float(norm2)));

  norm = norm * (3.0 - norm2 * norm * norm);
#else
  norm = native_rsqrt(norm2);
#endif

#ifdef USE_COMPLEX
  return (field) (norm, 0.0);
#else
  return norm;
#endif
}

inline field dlp_eval(real dx, real dy, real dz, real *N) {
  real norm, norm2;

  norm2 = dx * dx + dy * dy + dz * dz;

#ifndef USE_FLOAT
  norm = convert_double(native_rsqrt(convert_float(norm2)));
  norm = norm * (3.0 - norm2 * norm * norm);
#else
  norm = native_rsqrt(norm2);
#endif

#ifdef USE_COMPLEX
  return (field) ((N[0] * dx + N[1] * dy + N[2] * dz) * norm * norm * norm, 0.0);
#else
  return (N[0] * dx + N[1] * dy + N[2] * dz) * norm * norm * norm;
#endif
}

uint fast_select_quadrature( __global uint *geo_t, uint t, uint s, uint n) {
  uint p;

  p = 0;
  p += (geo_t[t] == geo_t[s]);
  p += (geo_t[t] == geo_t[s + n]);
  p += (geo_t[t] == geo_t[s + 2 * n]);
  p += (geo_t[t + n] == geo_t[s]);
  p += (geo_t[t + n] == geo_t[s + n]);
  p += (geo_t[t + n] == geo_t[s + 2 * n]);
  p += (geo_t[t + 2 * n] == geo_t[s]);
  p += (geo_t[t + 2 * n] == geo_t[s + n]);
  p += (geo_t[t + 2 * n] == geo_t[s + 2 * n]);

  return p;
}

uint select_quadrature( __global uint *geo_t, uint t, uint s, uint *tp,
    uint *sp, uint n) {

  uint p, q, i, j;

  p = 0;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      if (geo_t[t + i * n] == geo_t[s + j * n]) {
        tp[p] = i;
        sp[p] = j;
        p++;
      }
    }
  }

  q = p;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < q && geo_t[t + i * n] != geo_t[t + tp[j] * n]; j++)
    ;
    if (j == q)
    tp[q++] = i;
  }

  q = p;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < q && geo_t[s + i * n] != geo_t[s + sp[j] * n]; j++)
    ;
    if (j == q)
    sp[q++] = i;
  }

  return p;
}

real REAL_SQR(real x) {
  return x * x;
}

real REAL_SQRT(real x) {
  return sqrt(x);
}

real gramdet2(__global real *A, __global real *B, __global real *C) {
  real dx1, dx2, dx3, dy1, dy2, dy3, n[3];
  dx1 = B[0] - A[0];
  dx2 = B[1] - A[1];
  dx3 = B[2] - A[2];

  dy1 = C[0] - A[0];
  dy2 = C[1] - A[1];
  dy3 = C[2] - A[2];

  n[0] = dx2 * dy3 - dx3 * dy2;
  n[1] = dx3 * dy1 - dx1 * dy3;
  n[2] = dx1 * dy2 - dx2 * dy1;

  return REAL_SQR(n[0]) + REAL_SQR(n[1]) + REAL_SQR(n[2]);
}

void normal(__global real *A, __global real *B, __global real *C, real *N) {
  real dx1, dx2, dx3, dy1, dy2, dy3;
  dx1 = B[0] - A[0];
  dx2 = B[1] - A[1];
  dx3 = B[2] - A[2];

  dy1 = C[0] - A[0];
  dy2 = C[1] - A[1];
  dy3 = C[2] - A[2];

  N[0] = dx2 * dy3 - dx3 * dy2;
  N[1] = dx3 * dy1 - dx1 * dy3;
  N[2] = dx1 * dy2 - dx2 * dy1;
}

field slp_cc_dist(__constant real *quad, uint nq2, __global real *A_t,
    __global real *B_t, __global real *C_t, __global real *A_s,
    __global real *B_s, __global real *C_s) {

  uint i, j, k, l;
  real tx, sx, ty, sy, wi, wj, wk, w, dx, dy, dz, dx1, dy1, dz1, dx2, dy2, dz2,
  dx3, dy3, dz3;
  field sum;
  real A[3], BA_t[3], CB_t[3], BA_s[3], CB_s[3];
  __local real quad_x[20], quad_w[20];

  A[0] = A_t[0] - A_s[0];
  A[1] = A_t[1] - A_s[1];
  A[2] = A_t[2] - A_s[2];

  BA_t[0] = (B_t[0] - A_t[0]);
  BA_t[1] = (B_t[1] - A_t[1]);
  BA_t[2] = (B_t[2] - A_t[2]);

  CB_t[0] = (C_t[0] - B_t[0]);
  CB_t[1] = (C_t[1] - B_t[1]);
  CB_t[2] = (C_t[2] - B_t[2]);

  BA_s[0] = (B_s[0] - A_s[0]);
  BA_s[1] = (B_s[1] - A_s[1]);
  BA_s[2] = (B_s[2] - A_s[2]);

  CB_s[0] = (C_s[0] - B_s[0]);
  CB_s[1] = (C_s[1] - B_s[1]);
  CB_s[2] = (C_s[2] - B_s[2]);

  nq2 /= 2;

  for (i = 0; i < nq2; ++i) {
    quad_x[i] = quad[2 * i];
    quad_w[i] = quad[2 * i + 1];
  }

  sum = f_zero;

  for (i = 0; i < nq2; i += 1) {
    tx = quad_x[i];
    wi = tx * quad_w[i];

    dx1 = A[0] + BA_t[0] * tx;
    dy1 = A[1] + BA_t[1] * tx;
    dz1 = A[2] + BA_t[2] * tx;
    for (j = 0; j < nq2; j += 1) {
      sx = tx * quad_x[j];
      wj = wi * quad_w[j];

      dx2 = dx1 + CB_t[0] * sx;
      dy2 = dy1 + CB_t[1] * sx;
      dz2 = dz1 + CB_t[2] * sx;

      for (k = 0; k < nq2; k += 1) {
        ty = quad_x[k];
        wk = wj * ty * quad_w[k];

        dx3 = dx2 - BA_s[0] * ty;
        dy3 = dy2 - BA_s[1] * ty;
        dz3 = dz2 - BA_s[2] * ty;
        for (l = 0; l < nq2; l += 1) {
          sy = ty * quad_x[l];

          dx = dx3 - CB_s[0] * sy;
          dy = dy3 - CB_s[1] * sy;
          dz = dz3 - CB_s[2] * sy;

          sum += wk * quad_w[l] * slp_eval(dx, dy, dz);
        }
      }
    }
  }

  return sum * REAL_SQRT(gramdet2(A_t, B_t, C_t) * gramdet2(A_s, B_s, C_s));
}

field slp_cc_vert(__constant real *quad, uint nq2, __global real *A,
    __global real *B_t, __global real *C_t, __global real *B_s,
    __global real *C_s) {
  uint i, j, k, l;
  real tx, sx, ty, sy, wi, wj, wk, dx, dy, dz;
  field sum, sum2;
  real BA_t[3], CB_t[3], BA_s[3], CB_s[3];
  __local real quad_x[22], quad_w[22];

  BA_t[0] = (B_t[0] - A[0]);
  BA_t[1] = (B_t[1] - A[1]);
  BA_t[2] = (B_t[2] - A[2]);

  CB_t[0] = (C_t[0] - B_t[0]);
  CB_t[1] = (C_t[1] - B_t[1]);
  CB_t[2] = (C_t[2] - B_t[2]);

  BA_s[0] = (B_s[0] - A[0]);
  BA_s[1] = (B_s[1] - A[1]);
  BA_s[2] = (B_s[2] - A[2]);

  CB_s[0] = (C_s[0] - B_s[0]);
  CB_s[1] = (C_s[1] - B_s[1]);
  CB_s[2] = (C_s[2] - B_s[2]);

  nq2 /= 2;

  for (i = 0; i < nq2; ++i) {
    quad_x[i] = quad[2 * i];
    quad_w[i] = quad[2 * i + 1];
  }

  sum = f_zero;

  for (i = 0; i < nq2; i += 1) {
    wi = quad_w[i] * quad_x[i];
    for (j = 0; j < nq2; j += 1) {
      wj = wi * quad_w[j];
      for (k = 0; k < nq2; k += 1) {
        wk = wj * quad_w[k];
        for (l = 0; l < nq2; l += 1) {
          tx = quad_x[l];
          sx = quad_x[j] * quad_x[l];
          ty = quad_x[i] * quad_x[l];
          sy = quad_x[i] * quad_x[k] * quad_x[l];

          dx = BA_t[0] * tx + CB_t[0] * sx - BA_s[0] * ty - CB_s[0] * sy;
          dy = BA_t[1] * tx + CB_t[1] * sx - BA_s[1] * ty - CB_s[1] * sy;
          dz = BA_t[2] * tx + CB_t[2] * sx - BA_s[2] * ty - CB_s[2] * sy;

          sum2 = slp_eval(dx, dy, dz);

          dx = BA_t[0] * ty + CB_t[0] * sy - BA_s[0] * tx - CB_s[0] * sx;
          dy = BA_t[1] * ty + CB_t[1] * sy - BA_s[1] * tx - CB_s[1] * sx;
          dz = BA_t[2] * ty + CB_t[2] * sy - BA_s[2] * tx - CB_s[2] * sx;

          sum2 += slp_eval(dx, dy, dz);

          sum += sum2 * wk * quad_w[l] * quad_x[l] * quad_x[l] * quad_x[l];
        }
      }
    }
  }

  return sum * REAL_SQRT(gramdet2(A, B_t, C_t) * gramdet2(A, B_s, C_s));
}

field slp_cc_edge(__constant real *quad, uint nq2, __global real *A,
    __global real *B, __global real *C_t, __global real *C_s) {
  uint i, j, k, l;
  real eta1, eta2, xi1, tx, sx, sy, wi, wj, dx, dy, dz;
  field sum, sum2;
  real BA[3], CB_t[3], CB_s[3];

  BA[0] = B[0] - A[0];
  BA[1] = B[1] - A[1];
  BA[2] = B[2] - A[2];

  CB_t[0] = (C_t[0] - B[0]);
  CB_t[1] = (C_t[1] - B[1]);
  CB_t[2] = (C_t[2] - B[2]);

  CB_s[0] = (C_s[0] - B[0]);
  CB_s[1] = (C_s[1] - B[1]);
  CB_s[2] = (C_s[2] - B[2]);

  sum = f_zero;

  for (i = 0; i < nq2; i += 2) {
    eta1 = quad[i];
    wi = quad[i + 1];
    for (j = 0; j < nq2; j += 2) {
      eta2 = quad[j];
      wj = wi * quad[j + 1];
      for (k = 0; k < nq2; k += 2) {
        xi1 = quad[k];

        tx = xi1 * eta1;
        sx = xi1 * eta1 * eta2;
        sy = -xi1 * (eta1 - r_one);

        dx = BA[0] * tx + CB_t[0] * sx - CB_s[0] * sy;
        dy = BA[1] * tx + CB_t[1] * sx - CB_s[1] * sy;
        dz = BA[2] * tx + CB_t[2] * sx - CB_s[2] * sy;

        sum2 = slp_eval(dx, dy, dz);

        tx = xi1 * eta1 * eta2;
        sx = xi1 * eta1;
        sy = -xi1 * (eta1 * eta2 - r_one);

        dx = BA[0] * tx + CB_t[0] * sx - CB_s[0] * sy;
        dy = BA[1] * tx + CB_t[1] * sx - CB_s[1] * sy;
        dz = BA[2] * tx + CB_t[2] * sx - CB_s[2] * sy;
        sum2 += slp_eval(dx, dy, dz);

        tx = xi1 * (r_one - eta1);
        sx = xi1;
        sy = xi1 * eta1 * eta2;

        dx = BA[0] * tx + CB_t[0] * sx - CB_s[0] * sy;
        dy = BA[1] * tx + CB_t[1] * sx - CB_s[1] * sy;
        dz = BA[2] * tx + CB_t[2] * sx - CB_s[2] * sy;

        sum2 += slp_eval(dx, dy, dz);

        tx = xi1 * eta1 * (eta2 - r_one);
        sx = xi1 * eta1 * eta2;
        sy = xi1;

        dx = BA[0] * tx + CB_t[0] * sx - CB_s[0] * sy;
        dy = BA[1] * tx + CB_t[1] * sx - CB_s[1] * sy;
        dz = BA[2] * tx + CB_t[2] * sx - CB_s[2] * sy;

        sum2 += slp_eval(dx, dy, dz);

        tx = xi1 * (eta1 * eta2 - r_one);
        sx = xi1 * eta1 * eta2;
        sy = xi1 * eta1;

        dx = BA[0] * tx + CB_t[0] * sx - CB_s[0] * sy;
        dy = BA[1] * tx + CB_t[1] * sx - CB_s[1] * sy;
        dz = BA[2] * tx + CB_t[2] * sx - CB_s[2] * sy;

        sum2 += slp_eval(dx, dy, dz);

        tx = xi1 * (eta1 - r_one);
        sx = xi1 * eta1;
        sy = xi1 * eta1 * eta2;

        dx = BA[0] * tx + CB_t[0] * sx - CB_s[0] * sy;
        dy = BA[1] * tx + CB_t[1] * sx - CB_s[1] * sy;
        dz = BA[2] * tx + CB_t[2] * sx - CB_s[2] * sy;

        sum2 += slp_eval(dx, dy, dz);

        sum += sum2 * wj * quad[k + 1] * (r_one - xi1) * xi1 * xi1 * eta1;
      }
    }
  }

  return sum * REAL_SQRT(gramdet2(A, B, C_t) * gramdet2(A, B, C_s));
}

field slp_cc_iden(__constant real *quad, uint nq2, __global real *A,
    __global real *B, __global real *C) {
  uint i, j, k, l;
  real eta, xi1, xi2, tx, sx, wi, wj, wk, dx, dy, dz, dx1, dz1, dx2, dz2, dx3,
  dz3, sum3;
  field sum, sum2;

  dx1 = B[0] - A[0];
  dx2 = B[1] - A[1];
  dx3 = B[2] - A[2];

  dz1 = C[0] - B[0];
  dz2 = C[1] - B[1];
  dz3 = C[2] - B[2];

  sum = f_zero;

  for (i = 0; i < nq2; i += 2) {
    eta = quad[i];
    wi = quad[i + 1];
    for (j = 0; j < nq2; j += 2) {
      xi1 = quad[j];
      wj = wi * quad[j + 1];

      tx = xi1 * eta;
      sx = xi1;

      dx = dx1 * tx + dz1 * sx;
      dy = dx2 * tx + dz2 * sx;
      dz = dx3 * tx + dz3 * sx;

      sum2 = slp_eval(dx, dy, dz);

      //////////////////////////////////////////////////

      tx = xi1;
      sx = xi1 * eta;

      dx = dx1 * tx + dz1 * sx;
      dy = dx2 * tx + dz2 * sx;
      dz = dx3 * tx + dz3 * sx;

      sum2 += slp_eval(dx, dy, dz);

      //////////////////////////////////////////////////

      tx = xi1 * eta;
      sx = xi1 * (eta - r_one);

      dx = dx1 * tx + dz1 * sx;
      dy = dx2 * tx + dz2 * sx;
      dz = dx3 * tx + dz3 * sx;

      sum2 += slp_eval(dx, dy, dz);

      sum3 = 0.0;
      for (k = 0; k < nq2; k += 2) {
        xi2 = quad[k] * (r_one - xi1);
        wk = wj * quad[k + 1] * (r_one - xi1);
        for (l = 0; l < nq2; l += 2) {
          sum3 += wk * quad[l + 1] * (r_one - xi1 - xi2) * xi1;
        }
      }

      sum += r_two * sum2 * sum3;
    }
  }

  return sum * gramdet2(A, B, C);
}

__kernel void assemble_slp_cc_list_0(__constant real *xwq, uint nq2,
    __global uint *geo_t, __global real *geo_x, uint triangles,
    __global uint *ridx, __global uint *cidx, __global field *N, uint workitems) {

  uint index, tt, ss;
  field sum;
  __global real *A_t, *B_t, *C_t, *A_s, *B_s, *C_s;

  index = get_global_id(0);

  if (index < workitems) {
    tt = ridx[index];
    ss = cidx[index];

    A_t = geo_x + 3 * geo_t[tt + 0 * triangles];
    B_t = geo_x + 3 * geo_t[tt + 1 * triangles];
    C_t = geo_x + 3 * geo_t[tt + 2 * triangles];

    A_s = geo_x + 3 * geo_t[ss + 0 * triangles];
    B_s = geo_x + 3 * geo_t[ss + 1 * triangles];
    C_s = geo_x + 3 * geo_t[ss + 2 * triangles];

    sum = slp_cc_dist(xwq, nq2, A_t, B_t, C_t, A_s, B_s, C_s);

    N[index] = sum * KERNEL_CONST_2;
  }
}

__kernel void assemble_slp_cc_list_1(__constant real *xwq, uint nq2,
    __global uint *geo_t, __global real *geo_x, uint triangles,
    __global uint *ridx, __global uint *cidx, __global field *N, uint workitems) {

  uint index, tt, ss;
  uint tp[3], sp[3];
  field sum;
  __global real *A, *B_t, *C_t, *B_s, *C_s;

  index = get_global_id(0);

  if (index < workitems) {
    tt = ridx[index];
    ss = cidx[index];

    select_quadrature(geo_t, tt, ss, tp, sp, triangles);

    A = geo_x + 3 * geo_t[tt + tp[0] * triangles];
    B_t = geo_x + 3 * geo_t[tt + tp[1] * triangles];
    C_t = geo_x + 3 * geo_t[tt + tp[2] * triangles];

    B_s = geo_x + 3 * geo_t[ss + sp[1] * triangles];
    C_s = geo_x + 3 * geo_t[ss + sp[2] * triangles];

    sum = slp_cc_vert(xwq, nq2, A, B_t, C_t, B_s, C_s);

    N[index] = sum * KERNEL_CONST_2;
  }
}

__kernel void assemble_slp_cc_list_2(__constant real *xwq, uint nq2,
    __global uint *geo_t, __global real *geo_x, uint triangles,
    __global uint *ridx, __global uint *cidx, __global field *N, uint workitems) {

  uint index, tt, ss;
  uint tp[3], sp[3];
  field sum;
  __global real *A, *B, *C_t, *C_s;

  index = get_global_id(0);

  if (index < workitems) {
    tt = ridx[index];
    ss = cidx[index];

    select_quadrature(geo_t, tt, ss, tp, sp, triangles);

    A = geo_x + 3 * geo_t[tt + tp[0] * triangles];
    B = geo_x + 3 * geo_t[tt + tp[1] * triangles];
    C_t = geo_x + 3 * geo_t[tt + tp[2] * triangles];

    C_s = geo_x + 3 * geo_t[ss + sp[2] * triangles];

    sum = slp_cc_edge(xwq, nq2, A, B, C_t, C_s);

    N[index] = sum * KERNEL_CONST_2;
  }
}

__kernel void assemble_slp_cc_list_3(__constant real *xwq, uint nq2,
    __global uint *geo_t, __global real *geo_x, uint triangles,
    __global uint *ridx, __global uint *cidx, __global field *N, uint workitems) {

  uint index, tt, ss;
  field sum;
  __global real *A, *B, *C;

  index = get_global_id(0);

  if (index < workitems) {
    tt = ridx[index];
    ss = cidx[index];

    A = geo_x + 3 * geo_t[tt + 0 * triangles];
    B = geo_x + 3 * geo_t[tt + 1 * triangles];
    C = geo_x + 3 * geo_t[tt + 2 * triangles];

    sum = slp_cc_iden(xwq, nq2, A, B, C);

    N[index] = sum * KERNEL_CONST_2;
  }
}

//__kernel void assemble_kernel_c(__constant real *xwq, uint nq2,
//    __global uint *geo_t, __global real *geo_x, uint triangles,
//    __global real *Z, __global real *V, uint workitems) {
//
//  __global real *A, *B, *C;
//  real Z2[3];
//  real sum, tx, sx, w, x2, y2, z2, x, y, z;
//  uint index, tt, q1, q2;
//
//  index = get_global_id(0);
//
//  if (index < workitems) {
//
//    tt = *((__global uint*) (V + index) + 0);
//
//    A = geo_x + 3 * geo_t[tt + 0 * triangles];
//    B = geo_x + 3 * geo_t[tt + 1 * triangles];
//    C = geo_x + 3 * geo_t[tt + 2 * triangles];
//    Z2[0] = Z[index + 0 * workitems];
//    Z2[1] = Z[index + 1 * workitems];
//    Z2[2] = Z[index + 2 * workitems];
//
//    sum = 0.0;
//
//    for (q1 = 0; q1 < nq2; q1 += 2) {
//      tx = xwq[q1];
//      w = tx * xwq[q1 + 1];
//      x2 = A[0] - Z2[0] + (B[0] - A[0]) * tx;
//      y2 = A[1] - Z2[1] + (B[1] - A[1]) * tx;
//      z2 = A[2] - Z2[2] + (B[2] - A[2]) * tx;
//      for (q2 = 0; q2 < nq2; q2 += 2) {
//        sx = tx * xwq[q2];
//
//        x = x2 + (C[0] - B[0]) * sx;
//        y = y2 + (C[1] - B[1]) * sx;
//        z = z2 + (C[2] - B[2]) * sx;
//
//        sum += w * xwq[q2 + 1] * slp_eval(x, y, z);
//      }
//    }
//
//#if REALLY_FAST_DP_RSQRT
//    V[index] = KERNEL_CONST_2 * sum * REAL_SQRT(gramdet2(A, B, C));
//#else
//    V[index] = KERNEL_CONST_LAPLACEBEM3D * sum * REAL_SQRT(gramdet2(A, B, C));
//#endif
//  }
//}
//
//__kernel void assemble_lagrange_c(__constant real *xwq, uint nq2,
//    __global uint *geo_t, __global real *geo_x, uint triangles,
//    __constant real *xint, uint m, __global real *ab, __global real *V,
//    uint workitems) {
//
//  __global real *A, *B, *C;
//  real sum, lagr, tx, sx, w, x2, y2, z2, x, y, z, ax, bx, ay, by, az, bz;
//  uint index, j, jx, jy, jz, l, tt, q1, q2;
//
//  index = get_global_id(0);
//
//  if (index < workitems) {
//
//    tt = *((__global uint*) (V + index) + 0);
//    j = *((__global uint*) (V + index) + 1);
//    jx = j / (m * m);
//    jy = (j - jx * m * m) / m;
//    jz = (j - jx * m * m) % m;
//
//    A = geo_x + 3 * geo_t[tt + 0 * triangles];
//    B = geo_x + 3 * geo_t[tt + 1 * triangles];
//    C = geo_x + 3 * geo_t[tt + 2 * triangles];
//    ax = ab[index + 0 * workitems];
//    bx = ab[index + 1 * workitems];
//    ay = ab[index + 2 * workitems];
//    by = ab[index + 3 * workitems];
//    az = ab[index + 4 * workitems];
//    bz = ab[index + 5 * workitems];
//
//    sum = 0.0;
//    for (q1 = 0; q1 < nq2; q1 += 2) {
//      tx = xwq[q1];
//      w = tx * xwq[q1 + 1];
//      x2 = A[0] - ax + (B[0] - A[0]) * tx;
//      y2 = A[1] - ay + (B[1] - A[1]) * tx;
//      z2 = A[2] - az + (B[2] - A[2]) * tx;
//      for (q2 = 0; q2 < nq2; q2 += 2) {
//        sx = tx * xwq[q2];
//
//        x = x2 + (C[0] - B[0]) * sx;
//        y = y2 + (C[1] - B[1]) * sx;
//        z = z2 + (C[2] - B[2]) * sx;
//
//        lagr = 1.0;
//        for (l = 0; l < m; ++l) {
//          lagr = (l == jx) ? lagr : lagr * (x - bx * xint[l]);
//          lagr = (l == jy) ? lagr : lagr * (y - by * xint[l]);
//          lagr = (l == jz) ? lagr : lagr * (z - bz * xint[l]);
//        }
//
//        sum += w * xwq[q2 + 1] * lagr;
//      }
//    }
//
//    lagr = 1.0;
//    for (l = 0; l < m; ++l) {
//      lagr = (l == jx) ? lagr : lagr * bx * (xint[jx] - xint[l]);
//      lagr = (l == jy) ? lagr : lagr * by * (xint[jy] - xint[l]);
//      lagr = (l == jz) ? lagr : lagr * bz * (xint[jz] - xint[l]);
//    }
//
//    V[index] = sum / lagr * REAL_SQRT(gramdet2(A, B, C));
//  }
//}

field dlp_cc_dist(__constant real *quad, uint nq2, __global real *A_t,
    __global real *B_t, __global real *C_t, __global real *A_s,
    __global real *B_s, __global real *C_s, real *n_s) {

  uint i, j, k, l;
  real tx, sx, ty, sy, wi, wj, wk, w, dx, dy, dz, dx1, dy1, dz1, dx2, dy2, dz2,
  dx3, dy3, dz3;
  field sum;
  real A[3], BA_t[3], CB_t[3], BA_s[3], CB_s[3];
  __local real quad_x[20], quad_w[20];

  A[0] = A_t[0] - A_s[0];
  A[1] = A_t[1] - A_s[1];
  A[2] = A_t[2] - A_s[2];

  BA_t[0] = (B_t[0] - A_t[0]);
  BA_t[1] = (B_t[1] - A_t[1]);
  BA_t[2] = (B_t[2] - A_t[2]);

  CB_t[0] = (C_t[0] - B_t[0]);
  CB_t[1] = (C_t[1] - B_t[1]);
  CB_t[2] = (C_t[2] - B_t[2]);

  BA_s[0] = (B_s[0] - A_s[0]);
  BA_s[1] = (B_s[1] - A_s[1]);
  BA_s[2] = (B_s[2] - A_s[2]);

  CB_s[0] = (C_s[0] - B_s[0]);
  CB_s[1] = (C_s[1] - B_s[1]);
  CB_s[2] = (C_s[2] - B_s[2]);

  nq2 /= 2;

  for (i = 0; i < nq2; ++i) {
    quad_x[i] = quad[2 * i];
    quad_w[i] = quad[2 * i + 1];
  }

  sum = f_zero;

  for (i = 0; i < nq2; i += 1) {
    tx = quad_x[i];
    wi = tx * quad_w[i];

    dx1 = A[0] + BA_t[0] * tx;
    dy1 = A[1] + BA_t[1] * tx;
    dz1 = A[2] + BA_t[2] * tx;
    for (j = 0; j < nq2; j += 1) {
      sx = tx * quad_x[j];
      wj = wi * quad_w[j];

      dx2 = dx1 + CB_t[0] * sx;
      dy2 = dy1 + CB_t[1] * sx;
      dz2 = dz1 + CB_t[2] * sx;

      for (k = 0; k < nq2; k += 1) {
        ty = quad_x[k];
        wk = wj * ty * quad_w[k];

        dx3 = dx2 - BA_s[0] * ty;
        dy3 = dy2 - BA_s[1] * ty;
        dz3 = dz2 - BA_s[2] * ty;
        for (l = 0; l < nq2; l += 1) {
          sy = ty * quad_x[l];

          dx = dx3 - CB_s[0] * sy;
          dy = dy3 - CB_s[1] * sy;
          dz = dz3 - CB_s[2] * sy;

          sum += wk * quad_w[l] * dlp_eval(dx, dy, dz, n_s);
        }
      }
    }
  }

  return sum * REAL_SQRT(gramdet2(A_t, B_t, C_t));
}

field dlp_cc_vert(__constant real *quad, uint nq2, __global real *A,
    __global real *B_t, __global real *C_t, __global real *B_s,
    __global real *C_s, real *n_s) {
  uint i, j, k, l;
  real tx, sx, ty, sy, wi, wj, wk, dx, dy, dz;
  field sum, sum2;
  real BA_t[3], CB_t[3], BA_s[3], CB_s[3];
  __local real quad_x[22], quad_w[22];

  BA_t[0] = (B_t[0] - A[0]);
  BA_t[1] = (B_t[1] - A[1]);
  BA_t[2] = (B_t[2] - A[2]);

  CB_t[0] = (C_t[0] - B_t[0]);
  CB_t[1] = (C_t[1] - B_t[1]);
  CB_t[2] = (C_t[2] - B_t[2]);

  BA_s[0] = (B_s[0] - A[0]);
  BA_s[1] = (B_s[1] - A[1]);
  BA_s[2] = (B_s[2] - A[2]);

  CB_s[0] = (C_s[0] - B_s[0]);
  CB_s[1] = (C_s[1] - B_s[1]);
  CB_s[2] = (C_s[2] - B_s[2]);

  nq2 /= 2;

  for (i = 0; i < nq2; ++i) {
    quad_x[i] = quad[2 * i];
    quad_w[i] = quad[2 * i + 1];
  }

  sum = f_zero;

  for (i = 0; i < nq2; i += 1) {
    wi = quad_w[i] * quad_x[i];
    for (j = 0; j < nq2; j += 1) {
      wj = wi * quad_w[j];
      for (k = 0; k < nq2; k += 1) {
        wk = wj * quad_w[k];
        for (l = 0; l < nq2; l += 1) {
          tx = quad_x[l];
          sx = quad_x[j] * quad_x[l];
          ty = quad_x[i] * quad_x[l];
          sy = quad_x[i] * quad_x[k] * quad_x[l];

          dx = BA_t[0] * tx + CB_t[0] * sx - BA_s[0] * ty - CB_s[0] * sy;
          dy = BA_t[1] * tx + CB_t[1] * sx - BA_s[1] * ty - CB_s[1] * sy;
          dz = BA_t[2] * tx + CB_t[2] * sx - BA_s[2] * ty - CB_s[2] * sy;

          sum2 = dlp_eval(dx, dy, dz, n_s);

          dx = BA_t[0] * ty + CB_t[0] * sy - BA_s[0] * tx - CB_s[0] * sx;
          dy = BA_t[1] * ty + CB_t[1] * sy - BA_s[1] * tx - CB_s[1] * sx;
          dz = BA_t[2] * ty + CB_t[2] * sy - BA_s[2] * tx - CB_s[2] * sx;

          sum2 += dlp_eval(dx, dy, dz, n_s);

          sum += sum2 * wk * quad_w[l] * quad_x[l] * quad_x[l] * quad_x[l];
        }
      }
    }
  }

  return sum * REAL_SQRT(gramdet2(A, B_t, C_t));
}

field dlp_cc_edge(__constant real *quad, uint nq2, __global real *A,
    __global real *B, __global real *C_t, __global real *C_s, real *n_s) {
  uint i, j, k, l;
  real eta1, eta2, xi1, tx, sx, sy, wi, wj, dx, dy, dz;
  field sum, sum2;
  real BA[3], CB_t[3], CB_s[3];

  BA[0] = B[0] - A[0];
  BA[1] = B[1] - A[1];
  BA[2] = B[2] - A[2];

  CB_t[0] = (C_t[0] - B[0]);
  CB_t[1] = (C_t[1] - B[1]);
  CB_t[2] = (C_t[2] - B[2]);

  CB_s[0] = (C_s[0] - B[0]);
  CB_s[1] = (C_s[1] - B[1]);
  CB_s[2] = (C_s[2] - B[2]);

  sum = f_zero;

  for (i = 0; i < nq2; i += 2) {
    eta1 = quad[i];
    wi = quad[i + 1];
    for (j = 0; j < nq2; j += 2) {
      eta2 = quad[j];
      wj = wi * quad[j + 1];
      for (k = 0; k < nq2; k += 2) {
        xi1 = quad[k];

        tx = xi1 * eta1;
        sx = xi1 * eta1 * eta2;
        sy = -xi1 * (eta1 - r_one);

        dx = BA[0] * tx + CB_t[0] * sx - CB_s[0] * sy;
        dy = BA[1] * tx + CB_t[1] * sx - CB_s[1] * sy;
        dz = BA[2] * tx + CB_t[2] * sx - CB_s[2] * sy;

        sum2 = dlp_eval(dx, dy, dz, n_s);

        tx = xi1 * eta1 * eta2;
        sx = xi1 * eta1;
        sy = -xi1 * (eta1 * eta2 - r_one);

        dx = BA[0] * tx + CB_t[0] * sx - CB_s[0] * sy;
        dy = BA[1] * tx + CB_t[1] * sx - CB_s[1] * sy;
        dz = BA[2] * tx + CB_t[2] * sx - CB_s[2] * sy;
        sum2 += dlp_eval(dx, dy, dz, n_s);

        tx = xi1 * (r_one - eta1);
        sx = xi1;
        sy = xi1 * eta1 * eta2;

        dx = BA[0] * tx + CB_t[0] * sx - CB_s[0] * sy;
        dy = BA[1] * tx + CB_t[1] * sx - CB_s[1] * sy;
        dz = BA[2] * tx + CB_t[2] * sx - CB_s[2] * sy;

        sum2 += dlp_eval(dx, dy, dz, n_s);

        tx = xi1 * eta1 * (eta2 - r_one);
        sx = xi1 * eta1 * eta2;
        sy = xi1;

        dx = BA[0] * tx + CB_t[0] * sx - CB_s[0] * sy;
        dy = BA[1] * tx + CB_t[1] * sx - CB_s[1] * sy;
        dz = BA[2] * tx + CB_t[2] * sx - CB_s[2] * sy;

        sum2 += dlp_eval(dx, dy, dz, n_s);

        tx = xi1 * (eta1 * eta2 - r_one);
        sx = xi1 * eta1 * eta2;
        sy = xi1 * eta1;

        dx = BA[0] * tx + CB_t[0] * sx - CB_s[0] * sy;
        dy = BA[1] * tx + CB_t[1] * sx - CB_s[1] * sy;
        dz = BA[2] * tx + CB_t[2] * sx - CB_s[2] * sy;

        sum2 += dlp_eval(dx, dy, dz, n_s);

        tx = xi1 * (eta1 - r_one);
        sx = xi1 * eta1;
        sy = xi1 * eta1 * eta2;

        dx = BA[0] * tx + CB_t[0] * sx - CB_s[0] * sy;
        dy = BA[1] * tx + CB_t[1] * sx - CB_s[1] * sy;
        dz = BA[2] * tx + CB_t[2] * sx - CB_s[2] * sy;

        sum2 += dlp_eval(dx, dy, dz, n_s);

        sum += sum2 * wj * quad[k + 1] * (r_one - xi1) * xi1 * xi1 * eta1;
      }
    }
  }

  return sum * REAL_SQRT(gramdet2(A, B, C_t));
}

field dlp_cc_iden(__constant real *quad, uint nq2, __global real *A,
    __global real *B, __global real *C, real *n) {
  real alpha = 0.5;

#ifdef USE_COMPLEX
  return (field) (0.5 * alpha * REAL_SQRT(gramdet2(A, B, C)), 0.0);
#else
  return 0.5 * alpha * REAL_SQRT(gramdet2(A, B, C));
#endif
}

__kernel void assemble_dlp_cc_list_0(__constant real *xwq, uint nq2,
    __global uint *geo_t, __global real *geo_x, uint triangles,
    __global uint *ridx, __global uint *cidx, __global field *N, uint workitems) {

  uint index, tt, ss;
  field sum;
  real N_s[3];
  __global real *A_t, *B_t, *C_t, *A_s, *B_s, *C_s;

  index = get_global_id(0);

  if (index < workitems) {
    tt = ridx[index];
    ss = cidx[index];

    A_t = geo_x + 3 * geo_t[tt + 0 * triangles];
    B_t = geo_x + 3 * geo_t[tt + 1 * triangles];
    C_t = geo_x + 3 * geo_t[tt + 2 * triangles];

    A_s = geo_x + 3 * geo_t[ss + 0 * triangles];
    B_s = geo_x + 3 * geo_t[ss + 1 * triangles];
    C_s = geo_x + 3 * geo_t[ss + 2 * triangles];

    normal(A_s, B_s, C_s, N_s);

    sum = dlp_cc_dist(xwq, nq2, A_t, B_t, C_t, A_s, B_s, C_s, N_s);

    N[index] = sum * KERNEL_CONST_8;
  }
}

__kernel void assemble_dlp_cc_list_1(__constant real *xwq, uint nq2,
    __global uint *geo_t, __global real *geo_x, uint triangles,
    __global uint *ridx, __global uint *cidx, __global field *N, uint workitems) {

  uint index, tt, ss;
  uint tp[3], sp[3];
  field sum;
  real N_s[3];
  __global real *A, *B_t, *C_t, *B_s, *C_s;

  index = get_global_id(0);

  if (index < workitems) {
    tt = ridx[index];
    ss = cidx[index];

    select_quadrature(geo_t, tt, ss, tp, sp, triangles);

    A = geo_x + 3 * geo_t[tt + tp[0] * triangles];
    B_t = geo_x + 3 * geo_t[tt + tp[1] * triangles];
    C_t = geo_x + 3 * geo_t[tt + tp[2] * triangles];

    B_s = geo_x + 3 * geo_t[ss + sp[1] * triangles];
    C_s = geo_x + 3 * geo_t[ss + sp[2] * triangles];

    normal(geo_x + 3 * geo_t[ss + 0 * triangles],
        geo_x + 3 * geo_t[ss + 1 * triangles],
        geo_x + 3 * geo_t[ss + 2 * triangles], N_s);

    sum = dlp_cc_vert(xwq, nq2, A, B_t, C_t, B_s, C_s, N_s);

    N[index] = sum * KERNEL_CONST_8;
  }
}

__kernel void assemble_dlp_cc_list_2(__constant real *xwq, uint nq2,
    __global uint *geo_t, __global real *geo_x, uint triangles,
    __global uint *ridx, __global uint *cidx, __global field *N, uint workitems) {

  uint index, tt, ss;
  uint tp[3], sp[3];
  field sum;
  real N_s[3];
  __global real *A, *B, *C_t, *C_s;

  index = get_global_id(0);

  if (index < workitems) {
    tt = ridx[index];
    ss = cidx[index];

    select_quadrature(geo_t, tt, ss, tp, sp, triangles);

    A = geo_x + 3 * geo_t[tt + tp[0] * triangles];
    B = geo_x + 3 * geo_t[tt + tp[1] * triangles];
    C_t = geo_x + 3 * geo_t[tt + tp[2] * triangles];

    C_s = geo_x + 3 * geo_t[ss + sp[2] * triangles];

    normal(geo_x + 3 * geo_t[ss + 0 * triangles],
        geo_x + 3 * geo_t[ss + 1 * triangles],
        geo_x + 3 * geo_t[ss + 2 * triangles], N_s);

    sum = dlp_cc_edge(xwq, nq2, A, B, C_t, C_s, N_s);

    N[index] = sum * KERNEL_CONST_8;
  }
}

__kernel void assemble_dlp_cc_list_3(__constant real *xwq, uint nq2,
    __global uint *geo_t, __global real *geo_x, uint triangles,
    __global uint *ridx, __global uint *cidx, __global field *N, uint workitems) {

  uint index, tt, ss;
  field sum;
  real N_s[3];
  __global real *A, *B, *C;

  index = get_global_id(0);

  if (index < workitems) {
    tt = ridx[index];
    ss = cidx[index];

    A = geo_x + 3 * geo_t[tt + 0 * triangles];
    B = geo_x + 3 * geo_t[tt + 1 * triangles];
    C = geo_x + 3 * geo_t[tt + 2 * triangles];

    normal(A, B, C, N_s);

    sum = dlp_cc_iden(xwq, nq2, A, B, C, N_s);

    N[index] = sum;
  }
}

