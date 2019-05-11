/* ------------------------------------------------------------
 This is the file "helmholtzoclbem3d.c" of the H2Lib package.
 All rights reserved, Sven Christophersen 2015
 ------------------------------------------------------------ */

/**
 * @file helmholtzoclbem3d.c
 * @author Sven Christophersen
 * @date 2015
 */

#ifndef LIBRARY_HELMHOLTZOCLBEM3D_C_
#define LIBRARY_HELMHOLTZOCLBEM3D_C_

/* C STD LIBRARY */
/* CORE 0 */
/* CORE 1 */
/* CORE 2 */
/* CORE 3 */
/* SIMPLE */
/* PARTICLES */
/* BEM */
#include "oclbem3d.h"
#include "helmholtzbem3d.h"
#include "helmholtzbem3d.cl"

#ifdef USE_OPENMP
#ifdef USE_OPENCL
#ifdef USE_COMPLEX

static void
fill_cpu_wrapper_helmholtzbem3d(void *data)
{
  nearfield_args *nf_args = (nearfield_args *) data;
  pcbem3d   bem = nf_args->bem;

  if (nf_args->dist) {
    bem->nearfield_far(nf_args->ridx, nf_args->cidx, nf_args->bem,
		       nf_args->ntrans, nf_args->N);
  }
  else {
    bem->nearfield(nf_args->ridx, nf_args->cidx, nf_args->bem,
		   nf_args->ntrans, nf_args->N);
  }
}

static void
fill_cc_ocl_gpu_wrapper_helmholtzbem3d(void *data, uint kernel)
{
  merge_data_nf *mdata = (merge_data_nf *) data;

  cl_int    res;
  cl_uint   nq2;
  cl_uint   triangles;
  uint      current_device;
  uint      current_queue;
  uint      current_kernel_id;
  uint      current_thread;
  cl_uint   num_devices = ocl_system.num_devices;
  cl_uint   num_queues = ocl_system.queues_per_device;
  cl_event  h2d[2], calc;
  cl_kernel oclkernel;
  cl_command_queue queue;
  size_t    global_off[] = {
    0
  };
  size_t    local_size[] = {
    128
  };
  size_t    global_size[] = {
    ((mdata->pos + local_size[0] - 1) / local_size[0]) * local_size[0]
  };
  field     k = mdata->bem->k;
  field     alpha = mdata->bem->alpha;

  /****************************************************
   * Determine device, queue and kernel.
   ****************************************************/

  current_device = omp_get_thread_num() / num_queues;
  current_queue = omp_get_thread_num() % num_queues;
  current_thread = current_queue + current_device * num_queues;
  current_kernel_id = current_queue + current_device * num_queues
    + kernel * num_devices * num_queues;
  oclkernel = ocl_bem3d.kernels[current_kernel_id];
  queue = ocl_system.queues[current_thread];

  /****************************************************
   * Transfer input data.
   ****************************************************/

  res = clEnqueueWriteBuffer(queue, ocl_bem3d.mem_ridx[current_thread],
			     CL_FALSE, 0, mdata->pos * sizeof(uint),
			     mdata->ridx, 0, NULL, h2d + 0);
  CL_CHECK(res)

    res = clEnqueueWriteBuffer(queue, ocl_bem3d.mem_cidx[current_thread],
			       CL_FALSE, 0, mdata->pos * sizeof(uint),
			       mdata->cidx, 0, NULL, h2d + 1);
  CL_CHECK(res)

  /****************************************************
   * Setup kernel arguments for 'assemble_xxx_cc_list_z'
   ****************************************************/
    triangles = ocl_bem3d.triangles;

  if (kernel == 0 || kernel == 4) {
    nq2 = ocl_bem3d.nq;
    res = clSetKernelArg(oclkernel, 0, sizeof(cl_mem),
			 &ocl_bem3d.mem_q_xw[current_device]);
    CL_CHECK(res)
      res = clSetKernelArg(oclkernel, 1, sizeof(cl_uint), &nq2);
  }
  else {
    nq2 = ocl_bem3d.nq2;
    res = clSetKernelArg(oclkernel, 0, sizeof(cl_mem),
			 &ocl_bem3d.mem_q2_xw[current_device]);
    CL_CHECK(res)
      res = clSetKernelArg(oclkernel, 1, sizeof(cl_uint), &nq2);
  }
  CL_CHECK(res)
    res = clSetKernelArg(oclkernel, 2, sizeof(cl_mem),
			 &ocl_bem3d.mem_gr_t[current_device]);
  CL_CHECK(res)
    res = clSetKernelArg(oclkernel, 3, sizeof(cl_mem),
			 &ocl_bem3d.mem_gr_x[current_device]);
  CL_CHECK(res)
    res = clSetKernelArg(oclkernel, 4, sizeof(cl_uint), &triangles);
  CL_CHECK(res)
    res = clSetKernelArg(oclkernel, 5, sizeof(cl_mem),
			 &ocl_bem3d.mem_ridx[current_thread]);
  CL_CHECK(res)
    res = clSetKernelArg(oclkernel, 6, sizeof(cl_mem),
			 &ocl_bem3d.mem_cidx[current_thread]);
  CL_CHECK(res)
    res = clSetKernelArg(oclkernel, 7, sizeof(cl_mem),
			 &ocl_bem3d.mem_N[current_thread]);
  CL_CHECK(res)
    res = clSetKernelArg(oclkernel, 8, sizeof(cl_uint), &mdata->pos);
  CL_CHECK(res)
    res = clSetKernelArg(oclkernel, 9, sizeof(field), &k);
  CL_CHECK(res)
    res = clSetKernelArg(oclkernel, 10, sizeof(field), &alpha);
  CL_CHECK(res)

  /****************************************************
   * Invoke the kernel 'assemble_xxx_cc_list_z' on the GPU
   ****************************************************/
    res = clEnqueueNDRangeKernel(queue, oclkernel, 1, global_off, global_size,
				 local_size, 2, h2d, &calc);
  CL_CHECK(res)

  /****************************************************
   * Copy results back.
   ****************************************************/
    res =
    clEnqueueReadBuffer(queue, ocl_bem3d.mem_N[current_thread], CL_TRUE, 0,
			mdata->pos * sizeof(field), mdata->N, 1, &calc, NULL);
  CL_CHECK(res)

    clReleaseEvent(h2d[0]);
  clReleaseEvent(h2d[1]);
  clReleaseEvent(calc);
}

static void
fill_slp_cc_dist_gpu_wrapper_helmholtzbem3d(void *data)
{
  fill_cc_ocl_gpu_wrapper_helmholtzbem3d(data, 0);
}

static void
fill_slp_cc_vert_gpu_wrapper_helmholtzbem3d(void *data)
{
  fill_cc_ocl_gpu_wrapper_helmholtzbem3d(data, 1);
}

static void
fill_slp_cc_edge_gpu_wrapper_helmholtzbem3d(void *data)
{
  fill_cc_ocl_gpu_wrapper_helmholtzbem3d(data, 2);
}

static void
fill_slp_cc_iden_gpu_wrapper_helmholtzbem3d(void *data)
{
  fill_cc_ocl_gpu_wrapper_helmholtzbem3d(data, 3);
}

static void
fill_slp_cc_sing_cpu_wrapper_helmholtzbem3d(void *data)
{
  merge_data_nf *mdata = (merge_data_nf *) data;
  pcbem3d   bem = mdata->bem;
  const pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const real *gr_g = (const real *) gr->g;
  field    *aa = mdata->N;
  uint      pos = mdata->pos;
  field   **addr = mdata->addr;
  field     k = bem->k;

  const real *A_t, *B_t, *C_t, *A_s, *B_s, *C_s;
  const uint *tri_t, *tri_s;
  real     *xq, *yq, *wq;
  uint      tp[3], sp[3];
  real      Ax, Bx, Cx, Ay, By, Cy, tx, sx, ty, sy, dx, dy, dz, factor, norm,
    norm2, rnorm, base;
  field     sum;
  uint      q, nq, ss, tt, t;

  for (t = 0; t < pos; ++t) {
    tt = ((uint *) aa)[0];
    tri_t = gr_t[tt];
    ss = ((uint *) aa)[1];
    tri_s = gr_t[ss];
    factor = gr_g[ss] * gr_g[tt] * KERNEL_CONST_HELMHOLTZBEM3D;
    select_quadrature_singquad2d(bem->sq, tri_t, tri_s, tp, sp, &xq, &yq, &wq,
				 &nq, &base);
    wq += 9 * nq;

    A_t = gr_x[tri_t[tp[0]]];
    B_t = gr_x[tri_t[tp[1]]];
    C_t = gr_x[tri_t[tp[2]]];
    A_s = gr_x[tri_s[sp[0]]];
    B_s = gr_x[tri_s[sp[1]]];
    C_s = gr_x[tri_s[sp[2]]];

    sum = base;

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

      norm2 = dx * dx + dy * dy + dz * dz;
      rnorm = REAL_RSQRT(norm2);

      if (IMAG(k) == 0.0) {
	norm = REAL(k) * norm2 * rnorm;
	sum += wq[q] * (cos(norm) + I * sin(norm)) * rnorm;
      }
      else {
	norm = norm2 * rnorm;
	sum += wq[q] * exp(-IMAG(k) * norm)
	  * (cos(REAL(k) * norm) + I * sin(REAL(k) * norm)) * rnorm;
      }

    }
    *(addr[t]) = sum * factor;
  }
}

void
fill_slp_cc_near_task_helmholtzbem3d(const uint * ridx, const uint * cidx,
				     pcbem3d bem, bool ntrans, pamatrix N)
{
  nearfield_args *nf_args;

  nf_args = (nearfield_args *) allocmem(sizeof(nearfield_args));
  nf_args->ridx = ridx;
  nf_args->cidx = cidx;
  nf_args->dist = false;
  nf_args->bem = bem;
  nf_args->ntrans = ntrans;
  nf_args->N = N;

  if (nf == NULL) {
    op_cb.dist = fill_slp_cc_dist_gpu_wrapper_helmholtzbem3d;
    op_cb.vert = fill_slp_cc_vert_gpu_wrapper_helmholtzbem3d;
    op_cb.edge = fill_slp_cc_edge_gpu_wrapper_helmholtzbem3d;
    op_cb.iden = fill_slp_cc_iden_gpu_wrapper_helmholtzbem3d;
    op_cb.sing_cpu = fill_slp_cc_sing_cpu_wrapper_helmholtzbem3d;

    nf = new_ocltaskgroup(GPU_FIRST, NULL, merge_nf, cleanup_nf_merge,
			  distribute_nf, close_nf,
			  fill_cpu_wrapper_helmholtzbem3d,
			  fill_slp_cc_dist_gpu_wrapper_helmholtzbem3d,
			  getsize_nf, split_nf_task, cleanup_nf_task,
			  (void *) &op_cb);
  }

  add_task_taskgroup(&nf, (void *) nf_args);

}

void
fill_slp_cc_far_task_helmholtzbem3d(const uint * ridx, const uint * cidx,
				    pcbem3d bem, bool ntrans, pamatrix N)
{
  nearfield_args *nf_args;

  nf_args = (nearfield_args *) allocmem(sizeof(nearfield_args));
  nf_args->ridx = ridx;
  nf_args->cidx = cidx;
  nf_args->dist = true;
  nf_args->bem = bem;
  nf_args->ntrans = ntrans;
  nf_args->N = N;

  if (nf == NULL) {
    op_cb.dist = fill_slp_cc_dist_gpu_wrapper_helmholtzbem3d;
    op_cb.vert = fill_slp_cc_vert_gpu_wrapper_helmholtzbem3d;
    op_cb.edge = fill_slp_cc_edge_gpu_wrapper_helmholtzbem3d;
    op_cb.iden = fill_slp_cc_iden_gpu_wrapper_helmholtzbem3d;
    op_cb.sing_cpu = fill_slp_cc_sing_cpu_wrapper_helmholtzbem3d;

    nf = new_ocltaskgroup(GPU_FIRST, NULL, merge_nf, cleanup_nf_merge,
			  distribute_nf, close_nf,
			  fill_cpu_wrapper_helmholtzbem3d,
			  fill_slp_cc_dist_gpu_wrapper_helmholtzbem3d,
			  getsize_nf, split_nf_task, cleanup_nf_task,
			  (void *) &op_cb);
  }

  add_task_taskgroup(&nf, (void *) nf_args);

}

static void
fill_dlp_cc_dist_gpu_wrapper_helmholtzbem3d(void *data)
{
  fill_cc_ocl_gpu_wrapper_helmholtzbem3d(data, 4);
}

static void
fill_dlp_cc_vert_gpu_wrapper_helmholtzbem3d(void *data)
{
  fill_cc_ocl_gpu_wrapper_helmholtzbem3d(data, 5);
}

static void
fill_dlp_cc_edge_gpu_wrapper_helmholtzbem3d(void *data)
{
  fill_cc_ocl_gpu_wrapper_helmholtzbem3d(data, 6);
}

static void
fill_dlp_cc_iden_gpu_wrapper_helmholtzbem3d(void *data)
{
  fill_cc_ocl_gpu_wrapper_helmholtzbem3d(data, 7);
}

static void
fill_dlp_cc_sing_cpu_wrapper_helmholtzbem3d(void *data)
{
  merge_data_nf *mdata = (merge_data_nf *) data;
  pcbem3d   bem = mdata->bem;
  const pcsurface3d gr = bem->gr;
  const     real(*gr_x)[3] = (const real(*)[3]) gr->x;
  const     uint(*gr_t)[3] = (const uint(*)[3]) gr->t;
  const     real(*gr_n)[3] = (const real(*)[3]) gr->n;
  const real *gr_g = (const real *) gr->g;
  field    *aa = mdata->N;
  uint      pos = mdata->pos;
  field   **addr = mdata->addr;
  field     k = bem->k;

  const real *A_t, *B_t, *C_t, *A_s, *B_s, *C_s, *ns;
  const uint *tri_t, *tri_s;
  real     *xq, *yq, *wq;
  uint      tp[3], sp[3];
  real      Ax, Bx, Cx, Ay, By, Cy, tx, sx, ty, sy, dx, dy, dz, factor, norm,
    norm2, rnorm, base;
  field     sum;
  uint      q, nq, ss, tt, t;

  for (t = 0; t < pos; ++t) {
    tt = ((uint *) aa)[0];
    tri_t = gr_t[tt];
    ss = ((uint *) aa)[1];
    tri_s = gr_t[ss];
    ns = gr_n[ss];
    factor = gr_g[ss] * gr_g[tt] * KERNEL_CONST_HELMHOLTZBEM3D;

    if (tt == ss) {
      sum = 0.5 * bem->alpha * gr_g[tt];

      *(addr[t]) = sum;
    }
    else {

      (void) select_quadrature_singquad2d(bem->sq, tri_t, tri_s, tp, sp, &xq,
					  &yq, &wq, &nq, &base);
      wq += 9 * nq;

      A_t = gr_x[tri_t[tp[0]]];
      B_t = gr_x[tri_t[tp[1]]];
      C_t = gr_x[tri_t[tp[2]]];
      A_s = gr_x[tri_s[sp[0]]];
      B_s = gr_x[tri_s[sp[1]]];
      C_s = gr_x[tri_s[sp[2]]];

      sum = base;

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

	norm2 = dx * dx + dy * dy + dz * dz;
	rnorm = REAL_RSQRT(norm2);

	if (IMAG(k) == 0.0) {
	  norm = REAL(k) * norm2 * rnorm;
	  sum += wq[q] * (cos(norm) + I * sin(norm)) * rnorm * rnorm * rnorm
	    * (1.0 - I * norm) * (dx * ns[0] + dy * ns[1] + dz * ns[2]);
	}
	else {
	  norm = norm2 * rnorm;
	  sum += wq[q] * exp(-IMAG(k) * norm)
	    * (cos(REAL(k) * norm) + I * sin(REAL(k) * norm)) * rnorm
	    * rnorm * rnorm * (1.0 - I * k * norm)
	    * (dx * ns[0] + dy * ns[1] + dz * ns[2]);
	}
      }
      *(addr[t]) = sum * factor;
    }
  }
}

void
fill_dlp_cc_near_task_helmholtzbem3d(const uint * ridx, const uint * cidx,
				     pcbem3d bem, bool ntrans, pamatrix N)
{
  nearfield_args *nf_args;

  nf_args = (nearfield_args *) allocmem(sizeof(nearfield_args));
  nf_args->ridx = ridx;
  nf_args->cidx = cidx;
  nf_args->dist = false;
  nf_args->bem = bem;
  nf_args->ntrans = ntrans;
  nf_args->N = N;

  if (nf == NULL) {
    op_cb.dist = fill_dlp_cc_dist_gpu_wrapper_helmholtzbem3d;
    op_cb.vert = fill_dlp_cc_vert_gpu_wrapper_helmholtzbem3d;
    op_cb.edge = fill_dlp_cc_edge_gpu_wrapper_helmholtzbem3d;
    op_cb.iden = fill_dlp_cc_iden_gpu_wrapper_helmholtzbem3d;
    op_cb.sing_cpu = fill_dlp_cc_sing_cpu_wrapper_helmholtzbem3d;

    nf = new_ocltaskgroup(GPU_FIRST, NULL, merge_nf, cleanup_nf_merge,
			  distribute_nf, close_nf,
			  fill_cpu_wrapper_helmholtzbem3d,
			  fill_dlp_cc_dist_gpu_wrapper_helmholtzbem3d,
			  getsize_nf, split_nf_task, cleanup_nf_task,
			  (void *) &op_cb);
  }

  add_task_taskgroup(&nf, (void *) nf_args);

}

void
fill_dlp_cc_far_task_helmholtzbem3d(const uint * ridx, const uint * cidx,
				    pcbem3d bem, bool ntrans, pamatrix N)
{
  nearfield_args *nf_args;

  nf_args = (nearfield_args *) allocmem(sizeof(nearfield_args));
  nf_args->ridx = ridx;
  nf_args->cidx = cidx;
  nf_args->dist = true;
  nf_args->bem = bem;
  nf_args->ntrans = ntrans;
  nf_args->N = N;

  if (nf == NULL) {
    op_cb.dist = fill_dlp_cc_dist_gpu_wrapper_helmholtzbem3d;
    op_cb.vert = fill_dlp_cc_vert_gpu_wrapper_helmholtzbem3d;
    op_cb.edge = fill_dlp_cc_edge_gpu_wrapper_helmholtzbem3d;
    op_cb.iden = fill_dlp_cc_iden_gpu_wrapper_helmholtzbem3d;
    op_cb.sing_cpu = fill_dlp_cc_sing_cpu_wrapper_helmholtzbem3d;

    nf = new_ocltaskgroup(GPU_FIRST, NULL, merge_nf, cleanup_nf_merge,
			  distribute_nf, close_nf,
			  fill_cpu_wrapper_helmholtzbem3d,
			  fill_dlp_cc_dist_gpu_wrapper_helmholtzbem3d,
			  getsize_nf, split_nf_task, cleanup_nf_task,
			  (void *) &op_cb);
  }

  add_task_taskgroup(&nf, (void *) nf_args);

}

static void
init_helmholtzbem3d_opencl(pcbem3d bem)
{
  uint      num_kernels;
  const char *kernel_names[] = {
    "assemble_slp_cc_list_0", "assemble_slp_cc_list_1",
    "assemble_slp_cc_list_2", "assemble_slp_cc_list_3",
    "assemble_dlp_cc_list_0", "assemble_dlp_cc_list_1",
    "assemble_dlp_cc_list_2", "assemble_dlp_cc_list_3"
  };
  cl_uint   num_devices = ocl_system.num_devices;
  cl_uint   num_queues = ocl_system.queues_per_device;
  cl_uint   nthreads = num_devices * num_queues;
  cl_int    res;
  cl_mem_flags mem_rflags, mem_wflags;
  uint      i, j;
  real     *gr_x;
  uint     *gr_t;
  real     *gr_n;
  real     *q_xw, *x, *w;
  uint      q, v, t;

  num_kernels = sizeof(kernel_names) / sizeof(kernel_names[0]);
  mem_rflags = CL_MEM_READ_ONLY;
  mem_wflags = CL_MEM_WRITE_ONLY;

  if (ocl_bem3d.num_kernels == 0) {

    /****************************************************
     * Setup all necessary kernels
     ****************************************************/

    setup_kernels(helmholtzbem3d_ocl_src, num_kernels, kernel_names,
		  &ocl_bem3d.kernels);
    ocl_bem3d.num_kernels = num_kernels;

    /****************************************************
     * Create buffers for matrix chunks
     ****************************************************/

    ocl_bem3d.mem_N = (cl_mem *) allocmem(nthreads * sizeof(cl_mem));
    ocl_bem3d.mem_ridx = (cl_mem *) allocmem(nthreads * sizeof(cl_mem));
    ocl_bem3d.mem_cidx = (cl_mem *) allocmem(nthreads * sizeof(cl_mem));

    for (i = 0; i < num_devices; ++i) {
      for (j = 0; j < num_queues; ++j) {
	ocl_bem3d.mem_N[j + i * num_queues] =
	  clCreateBuffer(ocl_system.contexts[i], mem_wflags,
			 ocl_system.max_package_size, NULL, &res);
	CL_CHECK(res)

	  ocl_bem3d.mem_ridx[j + i * num_queues] =
	  clCreateBuffer(ocl_system.contexts[i], mem_rflags,
			 ocl_system.max_package_size / sizeof(field) *
			 sizeof(uint), NULL, &res);
	CL_CHECK(res)

	  ocl_bem3d.mem_cidx[j + i * num_queues] =
	  clCreateBuffer(ocl_system.contexts[i], mem_rflags,
			 ocl_system.max_package_size / sizeof(field) *
			 sizeof(uint), NULL, &res);
	CL_CHECK(res)
      }
    }

    /****************************************************
     * Create buffer for non singular quadrature rules
     * and copy the contents
     ****************************************************/

    q = bem->sq->q;

    x = allocreal(q);
    w = allocreal(q);
    q_xw = allocreal(2 * q);

    assemble_gauss(q, x, w);

    for (i = 0; i < q; ++i) {
      q_xw[2 * i] = 0.5 + 0.5 * x[i];
      q_xw[2 * i + 1] = 0.5 * w[i];
    }

    ocl_bem3d.mem_q_xw = (cl_mem *) allocmem(num_devices * sizeof(cl_mem));
    for (i = 0; i < num_devices; ++i) {
      ocl_bem3d.mem_q_xw[i] =
	clCreateBuffer(ocl_system.contexts[i], mem_rflags,
		       2 * q * sizeof(real), NULL, &res);
      CL_CHECK(res)
	res = clEnqueueWriteBuffer(ocl_system.queues[i * num_queues],
				   ocl_bem3d.mem_q_xw[i], CL_TRUE, 0,
				   2 * q * sizeof(real), q_xw, 0, NULL, NULL);
      CL_CHECK(res);
    }

    ocl_bem3d.nq = 2 * q;

    freemem(x);
    freemem(w);
    freemem(q_xw);

    /****************************************************
     * Create buffer for singular quadrature rules
     * and copy the contents
     ****************************************************/

    q = bem->sq->q2;

    x = allocreal(q);
    w = allocreal(q);
    q_xw = allocreal(2 * q);

    assemble_gauss(q, x, w);

    for (i = 0; i < q; ++i) {
      q_xw[2 * i] = 0.5 + 0.5 * x[i];
      q_xw[2 * i + 1] = 0.5 * w[i];
    }

    ocl_bem3d.mem_q2_xw = (cl_mem *) allocmem(num_devices * sizeof(cl_mem));
    for (i = 0; i < num_devices; ++i) {
      ocl_bem3d.mem_q2_xw[i] = clCreateBuffer(ocl_system.contexts[i],
					      mem_rflags,
					      2 * q * sizeof(real), NULL,
					      &res);
      CL_CHECK(res)
	res = clEnqueueWriteBuffer(ocl_system.queues[i * num_queues],
				   ocl_bem3d.mem_q2_xw[i], CL_TRUE, 0,
				   2 * q * sizeof(real), q_xw, 0, NULL, NULL);
      CL_CHECK(res);
    }

    ocl_bem3d.nq2 = 2 * q;

    freemem(x);
    freemem(w);
    freemem(q_xw);

    /****************************************************
     * Create buffers for geometry data and copy the contents
     ****************************************************/

    v = bem->gr->vertices;
    t = bem->gr->triangles;

    gr_x = allocreal(3 * v);
    gr_t = allocuint(3 * t);
    gr_n = allocreal(3 * t);

    for (i = 0; i < v; ++i) {
      gr_x[3 * i + 0] = bem->gr->x[i][0];
      gr_x[3 * i + 1] = bem->gr->x[i][1];
      gr_x[3 * i + 2] = bem->gr->x[i][2];
    }

    for (i = 0; i < t; ++i) {
      gr_t[i + 0 * t] = bem->gr->t[i][0];
      gr_t[i + 1 * t] = bem->gr->t[i][1];
      gr_t[i + 2 * t] = bem->gr->t[i][2];
    }

    for (i = 0; i < t; ++i) {
      gr_n[3 * i + 0] = bem->gr->n[i][0];
      gr_n[3 * i + 1] = bem->gr->n[i][1];
      gr_n[3 * i + 2] = bem->gr->n[i][2];
    }

    ocl_bem3d.mem_gr_t = (cl_mem *) allocmem(num_devices * sizeof(cl_mem));
    ocl_bem3d.mem_gr_x = (cl_mem *) allocmem(num_devices * sizeof(cl_mem));

    for (i = 0; i < num_devices; ++i) {
      ocl_bem3d.mem_gr_x[i] =
	clCreateBuffer(ocl_system.contexts[i], mem_rflags,
		       3 * v * sizeof(real), NULL, &res);
      CL_CHECK(res)

	ocl_bem3d.mem_gr_t[i] =
	clCreateBuffer(ocl_system.contexts[i], mem_rflags,
		       3 * t * sizeof(uint), NULL, &res);
      CL_CHECK(res)

	res = clEnqueueWriteBuffer(ocl_system.queues[i * num_queues],
				   ocl_bem3d.mem_gr_x[i], CL_TRUE, 0,
				   3 * v * sizeof(real), gr_x, 0, NULL, NULL);
      CL_CHECK(res);

      res = clEnqueueWriteBuffer(ocl_system.queues[i * num_queues],
				 ocl_bem3d.mem_gr_t[i], CL_TRUE, 0,
				 3 * t * sizeof(uint), gr_t, 0, NULL, NULL);
      CL_CHECK(res);
    }

    ocl_bem3d.triangles = t;

    freemem(gr_x);
    freemem(gr_t);
    freemem(gr_n);

    omp_init_lock(&nf_lock);
    omp_init_lock(&nf_dist_lock);
    omp_init_lock(&nf_vert_lock);
    omp_init_lock(&nf_edge_lock);
    omp_init_lock(&nf_iden_lock);
  }
}

static void
uninit_helmholtzbem3d_opencl()
{
  cl_uint   num_devices = ocl_system.num_devices;
  cl_uint   queues_per_device = ocl_system.queues_per_device;
  cl_uint   num_kernels = ocl_bem3d.num_kernels;

  uint      i, j, k;

  /****************************************************
   * Cleanup OpenCL objects belonging to 'ocl_bem3d'
   ****************************************************/

  for (i = 0; i < num_devices; ++i) {
    for (j = 0; j < queues_per_device; ++j) {
      clReleaseMemObject(ocl_bem3d.mem_N[j + i * queues_per_device]);
      clReleaseMemObject(ocl_bem3d.mem_ridx[j + i * queues_per_device]);
      clReleaseMemObject(ocl_bem3d.mem_cidx[j + i * queues_per_device]);
    }
    clReleaseMemObject(ocl_bem3d.mem_gr_t[i]);
    clReleaseMemObject(ocl_bem3d.mem_gr_x[i]);
    clReleaseMemObject(ocl_bem3d.mem_q_xw[i]);
    clReleaseMemObject(ocl_bem3d.mem_q2_xw[i]);
  }

  freemem(ocl_bem3d.mem_N);
  ocl_bem3d.mem_N = NULL;
  freemem(ocl_bem3d.mem_ridx);
  ocl_bem3d.mem_ridx = NULL;
  freemem(ocl_bem3d.mem_cidx);
  ocl_bem3d.mem_cidx = NULL;
  freemem(ocl_bem3d.mem_gr_t);
  ocl_bem3d.mem_gr_t = NULL;
  freemem(ocl_bem3d.mem_gr_x);
  ocl_bem3d.mem_gr_x = NULL;
  freemem(ocl_bem3d.mem_q_xw);
  ocl_bem3d.mem_q_xw = NULL;
  freemem(ocl_bem3d.mem_q2_xw);
  ocl_bem3d.mem_q2_xw = NULL;

  for (k = 0; k < num_kernels; ++k) {
    for (i = 0; i < num_devices; ++i) {
      for (j = 0; j < queues_per_device; ++j) {
	clReleaseKernel(ocl_bem3d.kernels[j + i * queues_per_device
					  +
					  k * num_devices *
					  queues_per_device]);
      }
    }
  }
  ocl_bem3d.num_kernels = 0;

  freemem(ocl_bem3d.kernels);
  ocl_bem3d.kernels = NULL;

  /****************************************************
   * Cleanup OpenCL objects belonging to 'ocl_system'
   ****************************************************/

  for (i = 0; i < num_devices; ++i) {
    for (j = 0; j < queues_per_device; ++j) {
      clReleaseCommandQueue(ocl_system.queues[j + i * queues_per_device]);
    }
    clReleaseContext(ocl_system.contexts[i]);
  }

  freemem(ocl_system.contexts);
  ocl_system.contexts = NULL;
  freemem(ocl_system.queues);
  ocl_system.queues = NULL;
  freemem(ocl_system.platforms);
  ocl_system.platforms = NULL;
  freemem(ocl_system.devices);
  ocl_system.devices = NULL;

  omp_destroy_lock(&nf_lock);
  omp_destroy_lock(&nf_dist_lock);
  omp_destroy_lock(&nf_vert_lock);
  omp_destroy_lock(&nf_edge_lock);
  omp_destroy_lock(&nf_iden_lock);
}

pbem3d
new_slp_helmholtz_ocl_bem3d(field k, pcsurface3d gr, uint q_regular,
			    uint q_singular, basisfunctionbem3d row_basis,
			    basisfunctionbem3d col_basis)
{

  pbem3d    bem;

  bem = new_slp_helmholtz_bem3d(k, gr, q_regular, q_singular, row_basis,
				col_basis);

  if (row_basis == BASIS_CONSTANT_BEM3D && col_basis == BASIS_CONSTANT_BEM3D) {
    bem->nearfield = fill_slp_cc_near_task_helmholtzbem3d;
    bem->nearfield_far = fill_slp_cc_far_task_helmholtzbem3d;
  }

  init_helmholtzbem3d_opencl(bem);

  return bem;
}

pbem3d
new_dlp_helmholtz_ocl_bem3d(field k, pcsurface3d gr, uint q_regular,
			    uint q_singular, basisfunctionbem3d row_basis,
			    basisfunctionbem3d col_basis, field alpha)
{
  pbem3d    bem;

  bem = new_dlp_helmholtz_bem3d(k, gr, q_regular, q_singular, row_basis,
				col_basis, alpha);

  if (row_basis == BASIS_CONSTANT_BEM3D && col_basis == BASIS_CONSTANT_BEM3D) {
    bem->nearfield = fill_dlp_cc_near_task_helmholtzbem3d;
    bem->nearfield_far = fill_dlp_cc_far_task_helmholtzbem3d;
  }

  init_helmholtzbem3d_opencl(bem);

  return bem;
}

void
del_helmholtz_ocl_bem3d(pbem3d bem)
{
  del_helmholtz_bem3d(bem);

  if (ocl_bem3d.kernels != NULL) {
    uninit_helmholtzbem3d_opencl();
  }
}

#endif
#endif
#endif

#endif
