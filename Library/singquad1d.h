/* ------------------------------------------------------------
 This is the file "singquad1d.h" of the H2Lib package.
 All rights reserved, Sven Christophersen 2010-2013
 ------------------------------------------------------------ */

/** @file singquad1d.h
 *  @author Sven Christophersen
 */

#ifndef SINGQUAD1D_H
#define SINGQUAD1D_H

/* C STD LIBRARY */
#include <math.h>
#include <stdio.h>
/* CORE 0 */
#include "settings.h"
/* CORE 1 */
#include "gaussquad.h"
/* CORE 2 */
/* CORE 3 */
/* SIMPLE */
/* PARTICLES */
/* BEM */

/** \defgroup singquad1d singquad1d
 *  @brief This module is responsible for standard and singular quadrature scheme
 *  that will be applied by @ref bem2d and derived modules such as @ref laplacebem2d.
 *  @{ */

/* ------------------------------------------------------------
 Definition of structs and types
 ------------------------------------------------------------ */

/**
 * @brief This struct collects all type of quadrature formulas needed by the
 * computation of matrix entries within BEM.
 */
struct _singquad1d {
	/** @brief X-component of quadrature points for singular integrals on same domain.*/
	real *x_id;
	/** @brief Y-component of quadrature points for singular integrals on same domain.*/
	real *y_id;
	/** @brief Quadrature weights for singular integrals on same domain.*/
	real *w_id;
	/** @brief Constant offset for singular integrals on same domain.*/
	real base_id;
	/** @brief Number of quadrature points for singular integrals on same domain.*/
	uint n_id;
	/**
	 * @brief X-component of quadrature points for singular integrals on domains
	 * sharing a common vertex.
	 */
	real *x_vert;
	/**
	 * @brief Y-component of quadrature points for singular integrals on domains
	 * sharing a common vertex.
	 */
	real *y_vert;
	/**
	 * @brief Quadrature weights for singular integrals on domains
	 * sharing a common vertex.
	 */
	real *w_vert;
	/**
	 * @brief Constant offset for singular integrals on domains
	 * sharing a common vertex.
	 */
	real base_vert;
	/**
	 * Number of quadrature points for singular integrals on domains
	 * sharing a common vertex.
	 */
	uint n_vert;
	/**
	 * @brief X-component of quadrature points for singular integrals on distant domains.
	 */
	real *x_dist;
	/**
	 * @brief Y-component of quadrature points for singular integrals on distant domains.
	 */
	real *y_dist;
	/** @brief Quadrature weights for singular integrals on distant domains.*/
	real *w_dist;
	/** @brief Constant offset for singular integrals on distant domains.*/
	real base_dist;
	/** Number of quadrature points for singular integrals on distant domains.*/
	uint n_dist;
	/** @brief X-component of quadrature points for integrals on a single domains.*/
	real *x_single;
	/** @brief Y-component of quadrature points for integrals on a single domains.*/
	real *y_single;
	/** @brief Quadrature weights for integrals on a single domains. */
	real *w_single;
	/** @brief Constant offset for integrals on a single domains.*/
	real base_single;
	/** Number of quadrature points for integrals on a single domains.*/
	uint n_single;
};

/**
 * @ref singquad1d is just an abbreviation for the struct _singquad1d.
 * It is necessary for the computation of singular integral arising in BEM
 * applications in 2 dimensional space.
 */
typedef struct _singquad1d singquad1d;

/**
 * Pointer to a @ref singquad1d object.
 */
typedef singquad1d* psingquad1d;

/**
 * Pointer to a constant @ref singquad1d object.
 */
typedef const singquad1d *pcsingquad1d;

/* ------------------------------------------------------------
 Constructors and destructors
 ------------------------------------------------------------ */

/**
 * @brief All quadrature rules needed for computation
 * of SLP Operator (identical egdes, edges with a common vertex, distant
 * edges)
 *
 * These quadrature rules a based on hierarchical quadrature from S. Börm and
 * W. Hackbusch.
 *
 * @param q Order of quadrature
 * @param x Gauss quadrature points in [-1,1].
 * @param w Gauss weights for above points.
 * @return Returns the newly created @ref singquad1d object.
 */
HEADER_PREFIX psingquad1d build_log_singquad1d(uint q, preal x, preal w);

/**
 * @brief All quadrature rules needed for computation
 * of DLP Operator (identical egdes, edges with a common vertex, distant
 * edges)
 *
 * These quadrature rules a based on hierarchical quadrature from S. Börm and
 * W. Hackbusch.
 *
 * @param q Order of quadrature.
 * @param x Gauss quadrature points in [-1,1].
 * @param w Gauss weights for above points.
 * @param alpha degree of homogeneity, for bem-problem value is -1.0.
 *  @return Returns the newly created @ref singquad1d object.
 */
HEADER_PREFIX psingquad1d build_pow_singquad1d(uint q, preal x, preal w,
		real alpha);

/**
 * @brief Destructor for @ref _singquad1d "singquad1d" objects.
 *
 * Delete a @ref _singquad1d "singquad1d" object.
 * @param sq @ref _singquad1d "singquad1d" object to be deleted.
 */
HEADER_PREFIX void del_singquad1d(psingquad1d sq);

/**
 * @brief This function is designed to select the correct quadrature rule for
 * a current pair of edges
 *
 * @param sq A @ref _singquad1d "singquad1d" object containing all necessary
 * quadrature rules.
 * @param tv An array defining the 2 vertices of edge @f$ t @f$.
 * @param sv An array defining the 2 vertices of edge @f$ s @f$.
 * @param tp Returning a permutation array of the vertices for @f$ t @f$.
 * @param sp Returning a permutation array of the vertices for @f$ s @f$.
 * @param x Returning the quadrature points for the edge @f$ t @f$.
 * @param y Returning the quadrature points for the edge @f$ s @f$.
 * @param w Returning the quadrature weights.
 * @param n Returning the total number of quadrature points.
 * @param base Returning a constant offset.
 * @return Returns the number of common vertices for edge @f$ t @f$ and
 * @f$ s @f$, which defines the current quadrature case.
 */
HEADER_PREFIX uint select_quadrature_singquad1d(pcsingquad1d sq, const uint *tv,
		const uint *sv, uint *tp, uint *sp, real (**x), real (**y), real **w,
		uint *n, real *base);

/** @} */

#endif
