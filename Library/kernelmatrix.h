
/* ------------------------------------------------------------
 * This is the file "kernelmatrix.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2018
 * ------------------------------------------------------------ */

/** @file kernelmatrix.h
 *  @author Steffen B&ouml;rm
 */

#ifndef KERNELMATRIX_H
#define KERNELMATRIX_H

/** @defgroup kernelmatrix kernelmatrix
 *  @brief Approximation of kernel matrices.
 *  @{ */

/** @brief Data required to approximate a kernel matrix. */
typedef struct _kernelmatrix kernelmatrix;

/** @brief Pointer to a @ref kernelmatrix object. */
typedef kernelmatrix *pkernelmatrix;

/** @brief Pointer to a constant @ref kernelmatrix object. */
typedef const kernelmatrix *pckernelmatrix;

#include "settings.h"
#include "h2matrix.h"
#include "cluster.h"
#include "clustergeometry.h"

/** @brief Representation of a kernel matrix an its approximation.
 *
 *  A kernel matrix is a matrix with entries of the form
 *  @f$g_{ij} = k(x_i,x_j)@f$, where @f$k@f$ is a kernel functions
 *  and @f$(x_i)_{i\in\Idx}@f$ are points in a suitable space.
 *
 *  If the kernel functions is locally smooth, it can be approximated
 *  by interpolation, and this gives rise to blockwise low-rank
 *  approximations. */
struct _kernelmatrix {
  /** @brief Spatial dimension. */
  uint dim;

  /** @brief Kernel function. */
  field (*kernel)(const real *xx, const real *yy, void *data);

  /** @brief Data for the kernel function. */
  void *data;

  /** @brief Number of points. */
  uint points;

  /** @brief Coordinates of points. */
  real **x;

  /** @brief Interpolation order (i.e., number of interpolation points). */
  uint m;

  /** @brief Interpolation points for the reference interval @f$[-1,1@f$. */
  real *xi_ref;
};

/** @brief Create an empty @ref kernelmatrix object.
 *
 *  @param dim Spatial dimension.
 *  @param points Number of points.
 *  @param m Interpolation order.
 *  @returns New object. */
HEADER_PREFIX pkernelmatrix
new_kernelmatrix(uint dim, uint points, uint m);

/** @brief Delete a @ref kernelmatrix object.
 *
 *  @param km Object to be deleted. */
HEADER_PREFIX void
del_kernelmatrix(pkernelmatrix km);

/** @brief Create a @ref clustergeometry object for a @ref kernelmatrix
 *  object.
 *
 *  @param km Description of the kernel matrix, particularly the points.
 *  @returns @ref clustergeometry object for the given points and dimension. */
HEADER_PREFIX pclustergeometry
creategeometry_kernelmatrix(pckernelmatrix km);

/** @brief Fill nearfield matrices.
 *
 *  @param ridx Row indices.
 *  @param cidx Column indices.
 *  @param km Description of the kernel matrix.
 *  @param N Nearfield matrix to be filled. */
HEADER_PREFIX void
fillN_kernelmatrix(const uint *ridx, const uint *cidx, pckernelmatrix km,
		   pamatrix N);

/** @brief Fill coupling matrices with interpolation coefficients.
 *
 *  @param rc Row cluster.
 *  @param cc Column cluster.
 *  @param km Description of the kernel matrix.
 *  @param S Coupling matrix to be filled. */
HEADER_PREFIX void
fillS_kernelmatrix(pccluster rc, pccluster cc,
		   pckernelmatrix km, pamatrix S);

/** @brief Fill leaf matrices for Lagrange interpolation.
 *
 *  @param tc Cluster.
 *  @param km Description of the kernel matrix.
 *  @param V Leaf matrix to be filled. */
HEADER_PREFIX void
fillV_kernelmatrix(pccluster tc,
		   pckernelmatrix km, pamatrix V);

/** @brief Fill transfer matrices for Lagrange interpolation.
 *
 *  @param sc Son cluster.
 *  @param fc Father cluster.
 *  @param km Description of the kernel matrix.
 *  @param E Transfer matrix to be filled. */
HEADER_PREFIX void
fillE_kernelmatrix(pccluster sc, pccluster fc,
		   pckernelmatrix km, pamatrix E);

/** @brief Fill a @ref clusterbasis using interpolation.
 *
 *  @param km Description of the kernel matrix.
 *  @param cb Cluster basis to be filled. */
HEADER_PREFIX void
fill_clusterbasis_kernelmatrix(pckernelmatrix km, pclusterbasis cb);

/** @brief Fill a @ref h2matrix using interpolation.
 *
 *  @param km Description of the kernel matrix.
 *  @param G Matrix to be filled. */
HEADER_PREFIX void
fill_h2matrix_kernelmatrix(pckernelmatrix km, ph2matrix G);

/** @} */

#endif
