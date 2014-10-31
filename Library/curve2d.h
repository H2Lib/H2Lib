/* ------------------------------------------------------------
 This is the file "curve2d.h" of the H2Lib package.
 All rights reserved, Steffen Boerm 2010
 ------------------------------------------------------------ */

/** @file curve2d.h
 *  @author Steffen B&ouml;rm */

#ifndef CURVE2D_H
#define CURVE2D_H

/**
 * @defgroup curve2d curve2d
 *
 * @brief Representation of a 2D curve by using a polygon.
 *
 * This module offers basic functions to create, load and save polygons curves.
 *
 * @{
 */

/**
 * @brief Abbreviation for struct @ref _curve2d .
 */
typedef struct _curve2d curve2d;

/**
 * @brief Abbreviation for a pointer to a @ref curve2d object.
 */
typedef curve2d *pcurve2d;

/**
 * @brief Abbreviation for a pointer to constant a @ref curve2d object.
 */
typedef const curve2d *pccurve2d;

#ifdef USE_CAIRO
#include <cairo/cairo.h>
#endif

#include "settings.h"

/**
 * @brief Representation of a polygon in 2D.
 */
struct _curve2d {
  /** @brief Number of vertices */
  uint vertices;
  /** @brief Number of edges */
  uint edges;

  /** @brief Vertex coordinates */
  real (*x)[2];
  /** @brief Edge vertices */
  uint (*e)[2];

  /** @brief Normal vectors */
  real (*n)[2];
  /** @brief Edge lengths = gram determinants*/
  preal g;
};

/* ------------------------------------------------------------
 constructors / destructors
 ------------------------------------------------------------ */

/**
 * @brief create a new @ref curve2d object with a certain number of
 * vertice and edges.
 *
 * This function will allocated storage for the members <tt>x</tt>, <tt>e</tt>,
 * <tt>n</tt> and <tt>g</tt> of the struct @ref _curve2d "curve2d".
 *
 * @param vertices Number of vertices for the new polygon object.
 * @param edges Number of edges for the new polygon object.
 * @return The newly created @ref curve2d object is returned.
 */
HEADER_PREFIX pcurve2d
new_curve2d(uint vertices, uint edges);

/**
 * @brief This function computes the normal vectors <tt>n</tt> , the gram
 * determinants <tt>g</tt> from the geometrical information.
 *
 * @param gr Geometry to be prepared.
 */
HEADER_PREFIX void
prepare_curve2d(pcurve2d gr);

/**
 * @brief Free Storage allocated for a @ref curve2d object.
 *
 * @param gr @ref curve2d "Curve2d" object to be deleted.
 */
HEADER_PREFIX void
del_curve2d(pcurve2d gr);

/* ------------------------------------------------------------
 curve creation
 ------------------------------------------------------------ */

/**
 * @brief Create a new polygon approximation of a circle with @f$ n @f$ edges
 * and radius @f$ r @f$ around the origin.
 *
 * @param edges The number of edges @f$n@f$ for the circle approximation.
 * @param r The radius @f$r@f$ of the circle approximation.
 * @return Returns a new @ref curve2d object, which is an approximation of
 * a circle of radius @f$r@f$.
 */
HEADER_PREFIX pcurve2d
new_circle_curve2d(uint edges, real r);

/**
 * @brief Create a new square with @f$ n @f$ edges
 * and edge length @f$ a @f$ around the origin.
 *
 * @param edges The number of edges @f$n@f$ for the square.
 * @param a The edge length @f$a@f$ of the square.
 * @return Returns a new @ref curve2d object, which is a square with edge length
 * @f$a@f$.
 */
HEADER_PREFIX pcurve2d
new_square_curve2d(uint edges, real a);

/**
 * @brief Create a new Hilbert-curve of order @f$ n @f$.
 * and edge length @f$ l @f$ around the origin.
 *
 * The number of edges for a Hilbert-curve of order @f$ n @f$ is @f$2^{2n} @f$.
 * The number of vertices for a Hilbert-curve of order @f$n@f$ is @f$2^{2n}-1@f$.
 * The Hilbert-curve is <b>NOT</b> closed.
 *
 * @param n The order of the Hilbert-curve.
 * @param l The edge length @f$l@f$ of the Hilbert-curve.
 * @return Returns a new @ref curve2d object, which is a Hilbert-curve
 * with edge length @f$l@f$.
 */
HEADER_PREFIX pcurve2d
new_hilbert_curve2d(uint n, real l);

/**
 * @brief Create a new polygon of a star with 16 spikes, @f$ n @f$ edges
 * and radius @f$ r @f$ around the origin.
 *
 * The spikes of the star lie on a circle with radius @f$r@f$, the inner dips
 * lie on a circle of radius @f$r/4@f$.
 *
 * @param edges The number of edges @f$n@f$ for the star. The number of edges
 * have to be a multiple of 16.
 * @param r The outer radius @f$r@f$ of the star.
 * @return Returns a new @ref curve2d object, which is a star of radius @f$r@f$.
 */
HEADER_PREFIX pcurve2d
new_star_curve2d(uint edges, real r);

/* ------------------------------------------------------------
 output functions
 ------------------------------------------------------------ */

/**
 * @brief print geometrical information to stdout.
 *
 * @param gr @ref curve2d "Curve2d" object to be printed out.
 */
HEADER_PREFIX void
print_curve2d(pccurve2d gr);

#ifdef USE_CAIRO
/**
 * @brief Draw the polygon to a cairo surface.
 *
 * @param gr @ref curve2d "Curve2d" object to be drawn.
 * @param cr Cairo surface to be drawn on.
 * @param scale Set the linewidth for the drawing.
 */
HEADER_PREFIX void
draw_curve2d(pccurve2d gr, cairo_t *cr, real scale);
#endif

/**
 * @}
 */

#endif
