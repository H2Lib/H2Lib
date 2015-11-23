/* ------------------------------------------------------------
 This is the file "realavector.h" of the H2Lib package.
 All rights reserved, Sven Christophersen 2015
 ------------------------------------------------------------ */

/** @file realavector.h
 @author Sven Christophersen
 */

#ifndef REALAVECTOR_H_
#define REALAVECTOR_H_

/* C STD LIBRARY */
#include <assert.h>
#include <math.h>
#include <stdlib.h>
/* CORE 0 */
#include "settings.h"
#include "basic.h"
#include "blas.h"
/* CORE 1 */
/* CORE 2 */
/* CORE 3 */
/* SIMPLE */
/* PARTICLES */
/* BEM */

/** @defgroup realavector realavector
 *  @brief Representation of a vector as an array for real coefficients.
 *
 *  The @ref realavector class is used to handle standard linear algebra
 *  operations like adding, scaling and multiplying vectors for real valued
 *  vectors.
 *  @attention If one needs to have a vector field entries from a general field,
 *             then the class @ref avector should be used instead.
 *  @{ */

/** Representation of a vector as an array. */
typedef struct _realavector realavector;

/** Pointer to a @ref avector object. */
typedef realavector *prealavector;

/** Pointer to a constant @ref avector object. */
typedef const realavector *pcrealavector;

/** Representation of a vector as an array. */
struct _realavector {
  /** @brief Vector coefficients. */
  real *v;

  /** @brief Vector dimension. */
  uint dim;

  /** @brief Points to owner of coefficient storage if this is a subvector. */
  void *owner;
};

/* ------------------------------------------------------------
 Constructors and destructors
 ------------------------------------------------------------ */

/** @brief Initialize an @ref realavector object.
 *
 *  Sets up the components of the object and allocates storage for
 *  the coefficient array.
 *
 *  @remark Should always be matched by a call to @ref uninit_realavector.
 *
 *  @param v Object to be initialized.
 *  @param dim Dimension of the new vector.
 *  @returns Initialized @ref realavector object. */
HEADER_PREFIX prealavector
init_realavector(prealavector v, uint dim);

/** @brief Initialize an @ref realavector object to represent a subvector.
 *
 *  Sets up the components of the object and uses part of the storage
 *  of another @ref realavector for the coefficient array, leading to a new
 *  vector representing a subvector of the source.
 *  Changes to the coefficients of the new vector also change
 *  coefficients of the source vector.
 *
 *  @remark Should always be matched by a call to @ref uninit_realavector that
 *  will <em>not</em> release the coefficient storage.
 *
 *  @param v Object to be initialized.
 *  @param src Source vector.
 *  @param dim Dimension of the new vector.
 *  @param off Offset in the source vector, should satisfy <tt>dim+off<=src->dim</tt>.
 *  @returns Initialized @ref realavector object. */
HEADER_PREFIX prealavector
init_sub_realavector(prealavector v, prealavector src, uint dim, uint off);
/** @brief Initialize an @ref realavector object using a given array for
 *  the coefficients.
 *
 *  Sets up the components of the object and uses the given array to
 *  represent the coefficients.
 *
 *  @remark Should always be matched by a call to @ref uninit_realavector that
 *  will <em>not</em> release the coefficient storage.
 *
 *  @param v Object to be initialized.
 *  @param src Source array, should contain at least <tt>dim</tt> elements.
 *  @param dim Dimension of the new vector.
 *  @returns Initialized @ref realavector object. */
HEADER_PREFIX prealavector
init_pointer_realavector(prealavector v, preal src, uint dim);

/** @brief Uninitialize an @ref realavector object.
 *
 *  Invalidates pointers, freeing corresponding storage if appropriate,
 *  and prepares the object for deletion.
 *
 *  @param v Object to be uninitialized. */
HEADER_PREFIX void
uninit_realavector(prealavector v);

/** @brief Create a new @ref realavector object.
 *
 *  Allocates storage for the object an sets up its components.
 *
 *  @remark Should always be matched by a call to @ref del_realavector.
 *
 *  @param dim Dimension of the new vector.
 *  @returns New @ref realavector object. */
HEADER_PREFIX prealavector
new_realavector(uint dim);

/** @brief Create a new @ref realavector object representing a subvector.
 *
 *  Allocates storage for the object, but uses part of the storage
 *  of another @ref realavector for the coefficient array, leading to a new
 *  vector representing a subvector of the source.
 *  Changes to the coefficients of the new vector also change
 *  coefficients of the source vector.
 *
 *  @remark Should always be matched by a call to @ref del_realavector that
 *  will <em>not</em> release the coefficient storage.
 *
 *  @param src Source vector.
 *  @param dim Dimension of the new vector.
 *  @param off Offset in the source vector, should satisfy <tt>dim+off<=src->dim</tt>.
 *  @returns New @ref realavector object. */
HEADER_PREFIX prealavector
new_sub_realavector(prealavector src, uint dim, uint off);

/** @brief Create a new @ref realavector object using a given array for
 *  the coefficients.
 *
 *  @param src Source array, should contain at least <tt>dim</tt> elements.
 *  @param dim Dimension of the new vector.
 *  @returns New @ref realavector object. */
HEADER_PREFIX prealavector
new_pointer_realavector(preal src, uint dim);

/** @brief Delete an @ref realavector object.
 *
 *  Releases the storage corresponding to the object.
 *
 *  @attention Make sure that there are no submatrix objects referring
 *  to this object, since they might otherwise keep using pointers
 *  to invalid storage.
 *
 *  @param v Object to be deleted. */
HEADER_PREFIX void
del_realavector(prealavector v);

/** @brief Change the dimension of an @ref realavector object without
 *  preserving its coefficients.
 *
 *  Allocates new storage for the coefficients and releases the
 *  old storage.
 *
 *  @attention Make sure that there are no submatrix objects referring
 *  to this object, since they will otherwise keep using pointers
 *  to invalid storage.
 *
 *  @param v Vector to be resized.
 *  @param dim New vector dimension. */
HEADER_PREFIX void
resize_realavector(prealavector v, uint dim);

/** @brief Reduce the dimension of an @ref realavector object without
 *  reallocating storage, preserving its coefficients.
 *
 *  @param v Vector to be resized.
 *  @param dim New vector dimension, not greater than <tt>v->dim</tt>. */
HEADER_PREFIX void
shrink_realavector(prealavector v, uint dim);

/* ------------------------------------------------------------
 Access methods
 ------------------------------------------------------------ */

#ifdef __GNUC__
INLINE_PREFIX real
getentry_realavector(pcrealavector, uint) __attribute__((const,unused));
INLINE_PREFIX void
setentry_realavector(prealavector, uint, real) __attribute__((unused));
INLINE_PREFIX real
addentry_realavector(prealavector, uint, real) __attribute__((unused));
#endif

/** @brief Read a vector entry @f$v_i@f$.
 *
 *  @param v Vector @f$v@f$.
 *  @param i Index @f$i@f$.
 *  @returns Vector entry @f$v_i@f$. */
INLINE_PREFIX real getentry_realavector(pcrealavector v, uint i) {
#ifdef FULL_DEBUG
  assert(i < v->dim);
#endif
  return v->v[i];
}

/** @brief Set a vector entry, @f$v_i \gets x@f$.
 *
 *  @param v Vector @f$v@f$.
 *  @param i Index @f$i@f$.
 *  @param x New value of @f$v_i@f$. */
INLINE_PREFIX void setentry_realavector(prealavector v, uint i, real x) {
#ifdef FULL_DEBUG
  assert(i < v->dim);
#endif

  v->v[i] = x;
}

/** @brief Add to a vector entry, @f$v_i \gets v_i + x@f$.
 *
 *  @param v Vector @f$v@f$.
 *  @param i Index @f$i@f$.
 *  @param x Summand.
 *  @returns New value of @f$v_i@f$. */
INLINE_PREFIX real addentry_realavector(prealavector v, uint i, real x) {
#ifdef FULL_DEBUG
  assert(i < v->dim);
#endif

  return (v->v[i] += x);
}

/* ------------------------------------------------------------
 Statistics
 ------------------------------------------------------------ */

/** @brief Get number of currently initialized @ref realavector objects.
 *
 *  Calls to initialization functions like @ref init_realavector and
 *  constructors like @ref new_realavector increase an internal counter,
 *  while @ref uninit_realavector and @ref del_realavector decrease it.
 *
 *  @remark Use this function to check whether a program correctly cleans
 *  up temporary variables.
 *
 *  @returns Number of currently initialized @ref realavector objects. */
HEADER_PREFIX uint
getactives_realavector();

/** @brief Get size of a given @ref realavector object.
 *
 *  Computes the size of the @ref realavector object and the storage
 *  allocated for the coefficients.
 *  If the object uses the coefficients of another object
 *  (e.g., if it was created using @ref new_sub_realavector), no coefficient
 *  storage is added.
 *
 *  @param v Vector object.
 *  @returns Size of allocated storage in bytes. */
HEADER_PREFIX size_t
getsize_realavector(pcrealavector v);

/** @brief Get heap size of a given @ref realavector object.
 *
 *  Computes the size of storage allocated for the coefficients
 *  on the heap, but not for the @ref realavector object itself.
 *  If the object uses the coefficients of another object
 *  (e.g., if it was created using @ref new_sub_realavector), no storage
 *  is required.
 *
 *  @param v Vector object.
 *  @returns Size of allocated heap storage in bytes. */
HEADER_PREFIX size_t
getsize_heap_realavector(pcrealavector v);

/* ------------------------------------------------------------
 Simple utility functions
 ------------------------------------------------------------ */

/** @brief Set a vector to zero.
 *
 *  @param v Target vector. */
HEADER_PREFIX void
clear_realavector(prealavector v);

/** @brief Set all coefficients in a vector to the same value.
 *
 *  @param v Target vector.
 *  @param x Fill value. */
HEADER_PREFIX void
fill_realavector(prealavector v, real x);

/** @brief Fill a vector with random values.
 *
 *  @param v Target vector. */
HEADER_PREFIX void
random_realavector(prealavector v);

/** @brief Copy a vector into another vector, @f$w \gets v@f$.
 *
 *  Both vectors have to be of the same dimension.
 *
 *  @remark If you want to copy between vectors of different dimension,
 *  consider using @ref copy_sub_realavector.
 *
 *  @param v Source vector.
 *  @param w Target vector. */
HEADER_PREFIX void
copy_realavector(pcrealavector v, prealavector w);

/** @brief Copy a vector into another vector, @f$w \gets v@f$.
 *
 *  If @f$w@f$ is smaller than @f$v@f$, only the first coefficients
 *  of @f$v@f$ are copied.
 *  If @f$v@f$ is smaller than @f$w@f$, only the first coefficients
 *  of @f$w@f$ are filled.
 *
 *  @param v Source vector.
 *  @param w Target vector. */
HEADER_PREFIX void
copy_sub_realavector(pcrealavector v, prealavector w);

/** @brief Print a vector.
 *
 *  @param v Vector @f$v@f$. */
HEADER_PREFIX void
print_realavector(pcrealavector v);

/* ------------------------------------------------------------
 Very basic linear algebra
 ------------------------------------------------------------ */

/** @brief Scale a vector @f$v@f$ by a factor @f$\alpha@f$,
 *  @f$v \gets \alpha v@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param v Target vector @f$v@f$. */
HEADER_PREFIX void
scale_realavector(real alpha, prealavector v);

/** @brief Compute the Euclidean norm @f$\|v\|_2@f$ of a vector @f$v@f$.
 *
 *  @param v Vector @f$v@f$.
 *  @return Returns @f$\|v\|_2@f$.
 */
HEADER_PREFIX real
norm2_realavector(pcrealavector v);

/** @brief Compute the Euclidean innner product @f$\langle x, y\rangle_2@f$
 *  of two vectors @f$x@f$ and @f$y@f$.
 *
 *  The Euclidean inner product is given by
 *  @f$\langle x, y \rangle_2 = \sum_i x_i y_i@f$.
 *
 *  @param x Vector @f$x@f$.
 *  @param y Vector @f$x@f$.
 *  @returns Euclidean inner product @f$\langle x, y\rangle_2@f$. */
HEADER_PREFIX real
dotprod_realavector(pcrealavector x, pcrealavector y);

/** @brief Add two vectors,
 *  @f$y \gets y + \alpha x@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
HEADER_PREFIX void
add_realavector(real alpha, pcrealavector x, prealavector y);

/** @} */

#endif /* REALAVECTOR_H_ */
