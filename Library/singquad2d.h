/* ------------------------------------------------------------
 This is the file "singquad2d.h" of the H2Lib package.
 All rights reserved, Sven Christophersen 2011-2013
 ------------------------------------------------------------ */

/** @file singquad2d.h
 *  @author Sven Christophersen
 */

#ifndef SINGQUAD2D_H_
#define SINGQUAD2D_H_

/* C STD LIBRARY */
#include <assert.h>
/* CORE 0 */
#ifdef USE_SIMD
#include "simd.h"
#endif
#include "settings.h"
/* CORE 1 */
#include "gaussquad.h"
/* CORE 2 */
/* CORE 3 */
/* SIMPLE */
/* PARTICLES */
/* BEM */
#include "surface3d.h"

/** \defgroup singquad2d singquad2d
 *  @brief This module is responsible for standard and singular quadrature scheme
 *  that will be applied by @ref bem3d and derived modules such as @ref laplacebem3d.
 *  @{ */

/**
 * @brief This struct collects all type of quadrature formulas needed by the
 * computation of matrix entries within BEM.
 */
struct _singquad2d {
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
   * sharing a common edge.
   */
  real *x_edge;
  /**
   * @brief Y-component of quadrature points for singular integrals on domains
   * sharing a common edge.
   */
  real *y_edge;
  /**
   * @brief Quadrature weights for singular integrals on domains
   * sharing a common edge.
   */
  real *w_edge;
  /**
   * @brief Constant offset for singular integrals on domains
   * sharing a common edge.
   */
  real base_edge;
  /**
   * @brief Number of quadrature points for singular integrals on domains
   * sharing a common edge.
   */
  uint n_edge;
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
   * @brief Number of quadrature points for singular integrals on domains
   * sharing a common vertex.
   */
  uint n_vert;
  /** @brief X-component of quadrature points for singular integrals on distant domains.*/
  real *x_dist;
  /** @brief Y-component of quadrature points for singular integrals on distant domains.*/
  real *y_dist;
  /** @brief Quadrature weights for singular integrals on distant domains.*/
  real *w_dist;
  /** @brief Constant offset for singular integrals on distant domains.*/
  real base_dist;
  /** @brief Number of quadrature points for singular integrals on distant domains.*/
  uint n_dist;
  /** @brief X-component of quadrature points for a single triangle.*/
  real *x_single;
  /** @brief Y-component of quadrature points for a single triangle.*/
  real *y_single;
  /** @brief Quadrature weights for a single triangle.*/
  real *w_single;
  /** @brief Constant offset for a single triangle.*/
  real base_single;
  /** @brief Number of quadrature points for a single triangle.*/
  uint n_single;

#ifdef USE_TRIQUADPOINTS
  /**
   * @brief 3D transformed quadrature points for every triangle - x-components.
   */
  real *tri_x;

  /**
   * @brief 3D transformed quadrature points for every triangle - y-components.
   */
  real *tri_y;

  /**
   * @brief 3D transformed quadrature points for every triangle - z-components.
   */
  real *tri_z;
#endif

  /** @brief Order of basic quadrature rule for single and regular double
   integrals.*/
  uint q;
  /** @brief Order of basic quadrature rule for singular double integrals.*/
  uint q2;
  /** @brief maximal number of quadrature points.*/
  uint nmax;
};

/**
 * @ref singquad2d is just an abbreviation for the struct _singquad2d.
 * It is necessary for the computation of singular integral arising in BEM
 * applications in 3 dimensional space.
 */
typedef struct _singquad2d singquad2d;

/**
 * Pointer to a @ref singquad2d object.
 */
typedef singquad2d* psingquad2d;

/**
 * Pointer to a constant @ref singquad2d object.
 */
typedef const singquad2d *pcsingquad2d;

/* ------------------------------------------------------------
 Constructors and destructors
 ------------------------------------------------------------ */

/**
 * @brief Creates a new @ref _singquad2d "singquad2d" object containing all
 * necessary quadrature rules for singular integrals arising in BEM applications
 * in 3 dimensional space.
 *
 * @param gr Currently used geometry.
 * @param q Order of gaussian quadrature rule used for construction of single
 * triangle quadrature rule and in case of distant domains.
 * @param q2 Order of gaussian quadrature rule used for construction of singular
 * quadrature cases.
 * @return returns a new @ref _singquad2d" "singquad2d" object.
 */
HEADER_PREFIX psingquad2d
build_singquad2d(pcsurface3d gr, uint q, uint q2);

/**
 * @brief Destructor for @ref _singquad2d "singquad2d" objects.
 *
 * Delete a @ref _singquad2d "singquad2d" object.
 * @param sq @ref _singquad2d "singquad2d" object to be deleted.
 */
HEADER_PREFIX void
del_singquad2d(psingquad2d sq);

/* ------------------------------------------------------------
 Weighting quadrature rules
 ------------------------------------------------------------ */

/**
 * @brief Weighting a quadrature rule for a double integral with both linear
 * basis functions.
 *
 * @param x Points within the first triangle.
 * @param y Points within the second triangle.
 * @param w Quadrature weights. Actually this is a 10 times nq two dimensional
 * array of weights. The first 9 columns will be used for the 9 combinations of
 * linear-linear basis functions weights and the last column will be remain for
 * the standard constant-constant case.
 * @param nq Number of quadrature points and weights stored within <tt>x, y, w</tt>.
 *
 * @attention The column <tt>w[9]</tt> has to be filled with valid quadrature
 * weights before calling this function.
 */
HEADER_PREFIX void
weight_basisfunc_ll_singquad2d(real *x, real *y, real *w, uint nq);

/**
 * @brief Weighting a quadrature rule for a double integral with a combination
 * of piecewise constant and linear basis functions.
 *
 * @param x Points within the first triangle using piecewise constant basis functions.
 * @param y Points within the second triangle using linear basis functions.
 * @param w Quadrature weights. Actually this is a 4 times nq two dimensional
 * array of weights. The first 3 columns will be used for the 3 combinations of
 * constant-linear basis functions weights and the last column will be remain for
 * the standard constant-constant case.
 * @param nq Number of quadrature points and weights stored within <tt>x, y, w</tt>.
 *
 * @attention The column <tt>w[3]</tt> has to be filled with valid quadrature
 * weights before calling this function.
 */
HEADER_PREFIX void
weight_basisfunc_cl_singquad2d(real *x, real *y, real *w, uint nq);

/**
 * @brief Weighting a quadrature rule for a double integral with a combination
 * of piecewise linear and constant basis functions.
 *
 * @param x Points within the first triangle using piecewise linear basis functions.
 * @param y Points within the second triangle using piecewise constant basis functions.
 * @param w Quadrature weights. Actually this is a 4 times nq two dimensional
 * array of weights. The first 3 columns will be used for the 3 combinations of
 * linear-constant basis functions weights and the last column will be remain for
 * the standard constant-constant case.
 * @param nq Number of quadrature points and weights stored within <tt>x, y, w</tt>.
 *
 * @attention The column <tt>w[3]</tt> has to be filled with valid quadrature
 * weights before calling this function.
 */
HEADER_PREFIX void
weight_basisfunc_lc_singquad2d(real *x, real *y, real *w, uint nq);

/**
 * @brief Weighting a quadrature rule for a single integral with linear basis
 * functions.
 *
 * @param x X-component of Points within the triangle.
 * @param y Y-component of Points within the triangle.
 * @param w Quadrature weights. Actually this is a 4 times nq two dimensional
 * array of weights. The first 3 columns will be used for the 3 combinations of
 * linear basis functions weights and the last column will be remain for
 * the standard constant case.
 * @param nq Number of quadrature points and weights stored within <tt>x, y, w</tt>.
 *
 * @attention The column <tt>w[3]</tt> has to be filled with valid quadrature
 * weights before calling this function.
 */
HEADER_PREFIX void
weight_basisfunc_l_singquad2d(real *x, real *y, real*w, uint nq);

/* ------------------------------------------------------------
 Select quadrature rule
 ------------------------------------------------------------ */

/**
 * @brief Determine the number of common vertices of a pair of triangles.
 *
 * @param geo_t Array of all triangles in the geometry.
 * @param t Index for the first triangle.
 * @param s Index for the second triangle.
 * @return Number of common vertices for @p t and @p s.
 */
HEADER_PREFIX uint
fast_select_quadrature(uint (*geo_t)[3], uint t, uint s);


/**
 * @brief This function is designed to select the correct quadrature rule for
 * a current pair of triangles
 *
 * @param sq A @ref _singquad2d "singquad2d" object containing all necessary
 * quadrature rules.
 * @param tv An array defining the 3 vertices of triangle @f$ t @f$.
 * @param sv An array defining the 3 vertices of triangle @f$ s @f$.
 * @param tp Returning a permutation array of the vertices for @f$ t @f$.
 * @param sp Returning a permutation array of the vertices for @f$ s @f$.
 * @param x Returning the quadrature points for the triangle @f$ t @f$.
 * @param y Returning the quadrature points for the triangle @f$ s @f$.
 * @param w Returning the quadrature weights.
 * @param n Returning the total number of quadrature points.
 * @param base Returning a constant offset.
 * @return Returns the number of common vertices for triangle @f$ t @f$ and
 * @f$ s @f$, which defines the current quadrature case.
 */
HEADER_PREFIX uint
select_quadrature_singquad2d(pcsingquad2d sq, const uint *tv, const uint *sv,
    uint *tp, uint *sp, real **x, real **y, real **w, uint *n, real *base);

/** @} */

#endif /* SINGQUAD2D_H_ */
