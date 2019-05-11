
/* ------------------------------------------------------------
 * This is the file "gaussquad.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2010
 * ------------------------------------------------------------ */

/** @file gaussquad.h
 *  @author Steffen B&ouml;rm */

#ifndef GAUSSQUAD_H
#define GAUSSQUAD_H

#include "settings.h"

/** @defgroup gaussquad gaussquad
 *  @brief Construction of one-dimensional quadrature rules.
 *  @{ */

/** @brief Construction of quadrature points and weights for
 *  Gaussian quadrature.
 *
 *  The quadrature points are the zeros of the @f$m@f$-th order
 *  Legendre polynomial.
 *  The Legendre polynomial can be expressed as the characteristic
 *  polynomial of a tridiagonal matrix, so its zeros can be computed
 *  by solving an eigenvalue problem.
 *  The quadrature weights can be obtained from the corresponding
 *  eigenvectors.
 *
 *  @param m Required number of quadrature points.
 *    The resulting quadrature rule will be exact for polynomials
 *    of order @f$2m-1@f$.
 *  @param x Should be an array of size @f$m@f$, will be overwritten
 *    by the quadrature points in @f$[-1,1]@f$.
 *  @param w Should be an array of size @f$m@f$, will be overwritten
 *    by the quadrature weights. */
HEADER_PREFIX void
assemble_gauss(uint m, preal x, preal w);

/** @} */

#endif
