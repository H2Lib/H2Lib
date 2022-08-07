
/* ------------------------------------------------------------
 * This is the file "ie1d.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2016
 * ------------------------------------------------------------ */

/** @file ie1d.h
 *  @author Steffen B&ouml;rm
 */

#ifndef IE1D_H
#define IE1D_H

/** @defgroup ie1d ie1d
 *  @brief Approximation of one-dimensional integral equation.
 *
 *  The @ref ie1d class is used to approximate the simple one-dimensional
 *  Fredholm integral equation
 *  @f[ -\int_0^1 \ln |x-y| u(y) \,dy = f(x) \qquad\text{ for all } x\in[0,1] .@f]
 *  @{ */

/** @brief Management object for a one-dimensional integral equation */
typedef struct _ie1d ie1d;

/** @brief Pointer to @ref ie1d object */
typedef ie1d *pie1d;

/** @brief Pointer to constant @ref ie1d object */
typedef const ie1d *pcie1d;

#include "cluster.h"
#include "hmatrix.h"
#include "h2matrix.h"

/** @brief Management object for one-dimensional integral equation */
struct _ie1d {
  /** @brief Number of basis functions */
  uint n;

  /** @brief Quadrature order */
  uint q;

  /** @brief Quadrature points in @f$[-1,1]@f$ */
  preal xq;

  /** @brief Quadrature weights */
  preal wq;

  /** @brief Maximal approximation order */
  uint mmax;

  /** @brief Interpolation points in @f$[-1,1]@f$ for all orders
   *  between <tt>0</tt> and <tt>mmax</tt> */
  preal *xi;

  /** @brief Approximation orders for all clusters */
  uint *m;

  /** @brief Callback function for nearfield entries */
  void (*nearfield)(const uint *ridx, const uint *cidx,
		    pcie1d ie, pamatrix N);

  /** @brief Callback function for admissible H-matrix blocks */
  void (*fillrk)(pccluster rc, uint rname,
		 pccluster cc, uint cname,
		 pcie1d ie, prkmatrix R);

  /** @brief Callback function for cluster basis leaf matrices */
  void (*fillV)(pccluster tc, uint tname,
		pcie1d ie, pamatrix V);

  /** @brief Callback function for cluster basis transfer matrices */
  void (*fillE)(pccluster sc, uint sname,
		pccluster fc, uint fname,
		pcie1d ie, pamatrix E);

  /** @brief Callback function for H2-matrix coupling matrices */
  void (*fillS)(pccluster rc, uint rname,
		pccluster cc, uint cname,
		pcie1d ie, pamatrix S);
};

/* ------------------------------------------------------------
 * Constructor and destructor
 * ------------------------------------------------------------ */

/** @brief Create a new @ref ie1d object.
 *
 * @remark Approximation orders and callback functions have to
 * be initialized, e.g., by @ref setup_aprx_taylor_ie1d or
 * @ref setup_aprx_interpolation_ie1d.
 *
 * @param n Number of basis functions.
 * @param q Quadrature order.
 * @returns New @ref ie1d object */
HEADER_PREFIX pie1d
new_ie1d(uint n, uint q);

/** @brief Delete a @ref ie1d object.
 *
 * @param ie Object to be deleted. */
HEADER_PREFIX void
del_ie1d(pie1d ie);

/* ------------------------------------------------------------
 * Simple cluster tree
 * ------------------------------------------------------------ */

/** @brief Create a cluster tree.
 *
 *  @param ie Problem description.
 *  @param leafsize Upper bound for the size of leaf clusters.
 *  @returns Geometrically regular cluster tree. */
HEADER_PREFIX pcluster
build_ie1d_cluster(pcie1d ie, uint leafsize);

/* ------------------------------------------------------------
 * Approximation orders
 * ------------------------------------------------------------ */

/** @brief Set up approximation orders.
 *
 *  The orders are given by @f$m_t = m_0 + m_1 (p - \level(t)@f$,
 *  where @f$p@f$ denotes the depth of the cluster tree.
 *
 *  @param ie Problem description.
 *  @param root Root cluster.
 *  @param m0 Constant term of the order function.
 *  @param m1 Linear term of the order function. */
HEADER_PREFIX void
prepare_orders_ie1d(pie1d ie, pccluster root, real m0, real m1);

/* ------------------------------------------------------------
 * Nearfield matrix
 * ------------------------------------------------------------ */

/** @brief Fill nearfield matrix.
 *
 *  @param ridx Array of row indices, at least <tt>N->rows</tt> entries.
 *  @param cidx Array of column indices, at least <tt>N->cols</tt> entries.
 *  @param ie Problem description.
 *  @param N Target matrix. */
HEADER_PREFIX void
nearfield_ie1d(const uint *ridx, const uint *cidx, pcie1d ie,
	       pamatrix N);

/* ------------------------------------------------------------
 * Approximation by Taylor expansion
 * ------------------------------------------------------------ */

/** @brief Approximate admissible submatrix of an H-matrix by Taylor expansion.
 *
 *  @param rc Row cluster.
 *  @param rname Name of row cluster.
 *  @param cc Column cluster.
 *  @param cname Name of column cluster.
 *  @param ie Problem description.
 *  @param r Target matrix. */
HEADER_PREFIX void
fillrk_taylor_ie1d(pccluster rc, uint rname,
		   pccluster cc, uint cname,
		   pcie1d ie, prkmatrix r);

/** @brief Leaf matrix for Taylor expansion.
 *
 *  @param t Cluster.
 *  @param tname Name of cluster.
 *  @param ie Problem description.
 *  @param V Target matrix. */
HEADER_PREFIX void
fillV_taylor_ie1d(pccluster t, uint tname,
		  pcie1d ie, pamatrix V);

/** @brief Transfer matrix for Taylor expansion.
 *
 *  @param s Son cluster.
 *  @param sname Name of son cluster.
 *  @param f Father cluster.
 *  @param fname Name of father cluster.
 *  @param ie Problem description.
 *  @param E Target matrix. */
HEADER_PREFIX void
fillE_taylor_ie1d(pccluster s, uint sname,
		  pccluster f, uint fname,
		  pcie1d ie, pamatrix E);

/** @brief Coupling matrix for Taylor expansion.
 *
 *  @param rc Row cluster.
 *  @param rname Name of row cluster.
 *  @param cc Column cluster.
 *  @param cname Name of column cluster.
 *  @param ie Problem description.
 *  @param S Target matrix. */
HEADER_PREFIX void
fillS_taylor_ie1d(pccluster rc, uint rname,
		  pccluster cc, uint cname,
		  pcie1d ie, pamatrix S);

/** @brief Initialize @ref ie1d object for approximation by
 *  Taylor expansion.
 *
 *  @param ie Problem description.
 *  @param root Root cluster.
 *  @param m0 Constant term of the order function.
 *  @param m1 Linear term of the order function. */
HEADER_PREFIX void
setup_aprx_taylor_ie1d(pie1d ie, pccluster root, real m0, real m1);

/* ------------------------------------------------------------
 * Approximation by interpolation
 * ------------------------------------------------------------ */

/** @brief Approximate admissible submatrix of an H-matrix by interpolation.
 *
 *  @param rc Row cluster.
 *  @param rname Name of row cluster.
 *  @param cc Column cluster.
 *  @param cname Name of column cluster.
 *  @param ie Problem description.
 *  @param r Target matrix. */
HEADER_PREFIX void
fillrk_interpolation_ie1d(pccluster rc, uint rname,
			  pccluster cc, uint cname,
			  pcie1d ie, prkmatrix r);

/** @brief Leaf matrix for interpolation.
 *
 *  @param t Cluster.
 *  @param tname Name of cluster.
 *  @param ie Problem description.
 *  @param V Target matrix. */
HEADER_PREFIX void
fillV_interpolation_ie1d(pccluster t, uint tname,
			 pcie1d ie, pamatrix V);

/** @brief Transfer matrix for interpolation.
 *
 *  @param s Son cluster.
 *  @param sname Name of son cluster.
 *  @param f Father cluster.
 *  @param fname Name of father cluster.
 *  @param ie Problem description.
 *  @param E Target matrix. */
HEADER_PREFIX void
fillE_interpolation_ie1d(pccluster s, uint sname,
			 pccluster f, uint fname,
			 pcie1d ie, pamatrix E);

/** @brief Coupling matrix for interpolation.
 *
 *  @param rc Row cluster.
 *  @param rname Name of row cluster.
 *  @param cc Column cluster.
 *  @param cname Name of column cluster.
 *  @param ie Problem description.
 *  @param S Target matrix. */
HEADER_PREFIX void
fillS_interpolation_ie1d(pccluster rc, uint rname,
			 pccluster cc, uint cname,
			 pcie1d ie, pamatrix S);

/** @brief Initialize @ref ie1d object for approximation by interpolation.
 *
 *  @param ie Problem description.
 *  @param root Root cluster.
 *  @param m0 Constant term of the order function.
 *  @param m1 Linear term of the order function. */
HEADER_PREFIX void
setup_aprx_interpolation_ie1d(pie1d ie, pccluster root, real m0, real m1);

/* ------------------------------------------------------------
 * Create matrices
 * ------------------------------------------------------------ */

/** @brief Fill an @ref hmatrix object.
 *
 *  @param ie Problem description.
 *  @param G Target matrix. */
HEADER_PREFIX void
fill_hmatrix_ie1d(pcie1d ie, phmatrix G);

/** @brief Fill a @ref clusterbasis object.
 *
 *  @param ie Problem description.
 *  @param cb Target cluster basis. */
HEADER_PREFIX void
fill_clusterbasis_ie1d(pcie1d ie, pclusterbasis cb);

/** @brief Fill a @ref h2matrix object.
 *
 *  @param ie Problem description.
 *  @param G Target matrix. */
HEADER_PREFIX void
fill_h2matrix_ie1d(pcie1d ie, ph2matrix G);

/** @} */

#endif
