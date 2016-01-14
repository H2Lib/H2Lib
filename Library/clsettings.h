/* ------------------------------------------------------------
 This is the file "clsettings.h" of the H2Lib package.
 All rights reserved, Sven Christophersen 2015
 ------------------------------------------------------------ */

/**
 * @file clsettings.h
 * @author Sven Christophersen
 * @date 2015
 */

#ifndef CLSETTINGS_H_
#define CLSETTINGS_H_

/** @defgroup clsettings clsettings
 *  @brief Fundamental types and macros for OpenCL computations.
 *  @{ */

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#if USE_OPENCL
/** @cond DEBUG */
#define OCL_SYNTAX_HIGHLIGHT 1
/** @endcond */
#endif

#if OCL_SYNTAX_HIGHLIGHT
/** @cond DEBUG */
#define __kernel
#define __global
#define __local
#define __constant
/** @endcond */
#endif

/** @brief @ref real floating point type.
 *
 *  This type is used, e.g., for geometric coordinates, norms
 *  and diagonal elements of self-adjoint matrices. */
#ifdef USE_FLOAT
typedef float real;
#else
typedef double real;
#endif

/** @brief Field type.
 *
 *  This type is used in the linear algebra modules to represent
 *  the coefficients of matrices and vectors. */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
typedef float2 field;
#else
typedef float field;
#endif
#else
#ifdef USE_COMPLEX
typedef double2 field;
#else
typedef double field;
#endif
#endif

/** @brief Define the imaginary unit. */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
#define I ((field)(0.0f,1.0f))
#else

#endif
#else
#ifdef USE_COMPLEX
#define I ((field)(0.0,1.0))
#else
#endif
#endif

/** @brief @ref real constant zero */
__constant static real r_zero;

/** @brief @ref real constant one */
__constant static real r_one;

/** @brief @ref real constant minus one */
__constant static real r_minusone;

/** @brief @ref real constant two */
__constant static real r_two;

/** @brief @ref field constant zero */
__constant static field f_zero;

/** @brief @ref field constant one */
__constant static field f_one;

/** @brief @ref field constant minus one */
__constant static field f_minusone;

#ifdef USE_FLOAT
__constant static real r_zero = 0.0f;
__constant static real r_one = 1.0f;
__constant static real r_minusone = -1.0f;
__constant static real r_two = 2.0f;

#ifdef USE_COMPLEX
__constant static field f_zero = (field) (0.0f,0.0f);
__constant static field f_one = (field) (1.0f, 0.0f);
__constant static field f_minusone = (field) (-1.0f, 0.0f);
#else
__constant static field f_zero = 0.0f;
__constant static field f_one = 1.0f;
__constant static field f_minusone = -1.0f;
#endif
#else
__constant static real r_zero = 0.0;
__constant static real r_one = 1.0;
__constant static real r_minusone = -1.0;
__constant static real r_two = 2.0;

#ifdef USE_COMPLEX
__constant static field f_zero = (field) (0.0, 0.0);
__constant static field f_one = (field) (1.0, 0.0);
__constant static field f_minusone = (field) (-1.0, 0.0);
#else
__constant static field f_zero = 0.0;
__constant static field f_one = 1.0;
__constant static field f_minusone = -1.0;
#endif
#endif

/****************************************************
 * basis operations on complex numbers
 ****************************************************/

#ifdef USE_COMPLEX
/**
 * @brief Returns the real part of a complex number z.
 *
 * @param z Complex number given by z = a + b * I
 * @return Real part of z: Re(z) = a.
 */
inline real REAL(field z) {
  return z.x;
}

/**
 * @brief Returns the imaginary part of a complex number z.
 *
 * @param z Complex number given by z = a + b * I
 * @return Imaginary part of z: IM(z) = b.
 */
inline real IMAG(field z) {
  return z.y;
}

/**
 * @brief Performs a complex multiplication of two complex numbers x and y.
 *
 * @param x 1st factor of the complex product.
 * @param y 2nd factor of the complex product.
 * @return Returns the product @f$ x \times y@f$.
 */
inline field cmul(field x, field y) {
  return (field) (x.x * y.x - x.y * y.y, x.x * y.y + x.y * y.x);
}
#else

#endif

/**
 * @}
 */

#endif /* CLSETTINGS_H_ */
