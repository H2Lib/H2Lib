/* ------------------------------------------------------------
 This is the file "avector.h" of the H2Lib package.
 All rights reserved, Steffen Boerm 2009
 ------------------------------------------------------------ */

/** @file avector.h
 @author Steffen B&ouml;rm
 */

#ifndef AVECTOR_H
#define AVECTOR_H

/** @defgroup avector avector
 *  @brief Representation of a vector as an array.
 *
 *  The @ref avector class is used to handle standard linear algebra
 *  operations like adding, scaling and multiplying vectors.
 *  @{ */

/** Representation of a vector as an array. */
typedef struct _avector avector;

/** Pointer to a @ref avector object. */
typedef avector *pavector;

/** Pointer to a constant @ref avector object. */
typedef const avector *pcavector;

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "basic.h"
#include "blas.h"
#include "amatrix.h"
#include "settings.h"

/** Representation of a vector as an array. */
struct _avector {
  /** @brief Vector coefficients. */
  field *v;

  /** @brief Vector dimension. */
  uint dim;

  /** @brief Points to owner of coefficient storage if this is a subvector. */
  void *owner;
};

/* ------------------------------------------------------------
 Constructors and destructors
 ------------------------------------------------------------ */

/** @brief Initialize an @ref avector object.
 *
 *  Sets up the components of the object and allocates storage for
 *  the coefficient array.
 *
 *  @remark Should always be matched by a call to @ref uninit_avector.
 *
 *  @param v Object to be initialized.
 *  @param dim Dimension of the new vector.
 *  @returns Initialized @ref avector object. */
HEADER_PREFIX pavector
init_avector(pavector v, uint dim);

/** @brief Initialize an @ref avector object to represent a subvector.
 *
 *  Sets up the components of the object and uses part of the storage
 *  of another @ref avector for the coefficient array, leading to a new
 *  vector representing a subvector of the source.
 *  Changes to the coefficients of the new vector also change
 *  coefficients of the source vector.
 *
 *  @remark Should always be matched by a call to @ref uninit_avector that
 *  will <em>not</em> release the coefficient storage.
 *
 *  @param v Object to be initialized.
 *  @param src Source vector.
 *  @param dim Dimension of the new vector.
 *  @param off Offset in the source vector, should satisfy <tt>dim+off<=src->dim</tt>.
 *  @returns Initialized @ref avector object. */
HEADER_PREFIX pavector
init_sub_avector(pavector v, pavector src, uint dim, uint off);

/** @brief Initialize an @ref avector object and set it to zero
 * 
 * Sets up the components of the object, allocates storage for the
 * coefficient array, and dets it to zero.
 * @param v Object to be initialized.
 * @param dim Dimension of the vector.
 * @returns Initialized @ref avector object.*/
HEADER_PREFIX pavector
init_zero_avector(pavector v, uint dim);

/** @brief Initialize an @ref avector object to represent a column vector
 *  of a given matrix.
 *
 *  Sets up the components of the object and uses part of the storage
 *  of an @ref amatrix corresponding to one of its columns.
 *  Changes to the coefficients of the new vector also change the column
 *  of the source matrix.
 *
 *  @remark Should always be matched by a call to @ref uninit_avector that
 *  will <em>not</em> release the coefficient storage.
 *
 *  @param v Object to be initialized.
 *  @param src Source matrix.
 *  @param col Column of the source matrix that will be used for the vector.
 *  @returns Initialized @ref avector object. */
HEADER_PREFIX pavector
init_column_avector(pavector v, pamatrix src, uint col);

/** @brief Initialize an @ref avector object using a given array for
 *  the coefficients.
 *
 *  Sets up the components of the object and uses the given array to
 *  represent the coefficients.
 *
 *  @remark Should always be matched by a call to @ref uninit_avector that
 *  will <em>not</em> release the coefficient storage.
 *
 *  @param v Object to be initialized.
 *  @param src Source array, should contain at least <tt>dim</tt> elements.
 *  @param dim Dimension of the new vector.
 *  @returns Initialized @ref avector object. */
HEADER_PREFIX pavector
init_pointer_avector(pavector v, pfield src, uint dim);

/** @brief Uninitialize an @ref avector object.
 *
 *  Invalidates pointers, freeing corresponding storage if appropriate,
 *  and prepares the object for deletion.
 *
 *  @param v Object to be uninitialized. */
HEADER_PREFIX void
uninit_avector(pavector v);

/** @brief Create a new @ref avector object.
 *
 *  Allocates storage for the object an sets up its components.
 *
 *  @remark Should always be matched by a call to @ref del_avector.
 *
 *  @param dim Dimension of the new vector.
 *  @returns New @ref avector object. */
HEADER_PREFIX pavector
new_avector(uint dim);

/** @brief Create a new @ref avector object representing a subvector.
 *
 *  Allocates storage for the object, but uses part of the storage
 *  of another @ref avector for the coefficient array, leading to a new
 *  vector representing a subvector of the source.
 *  Changes to the coefficients of the new vector also change
 *  coefficients of the source vector.
 *
 *  @remark Should always be matched by a call to @ref del_avector that
 *  will <em>not</em> release the coefficient storage.
 *
 *  @param src Source vector.
 *  @param dim Dimension of the new vector.
 *  @param off Offset in the source vector, should satisfy <tt>dim+off<=src->dim</tt>.
 *  @returns New @ref avector object. */
HEADER_PREFIX pavector
new_sub_avector(pavector src, uint dim, uint off);

/** @brief Create a new @ref avector object representing a zero vector.
 * 
 * Allocates storage for the object and sets all coefficients to 
 * zero.
 * 
 * @remark Should always be matched by a call to @ref del_avector.
 * 
 * @param dim Dimension of the new vector.
 * @returns New @ref amatrix object. */
HEADER_PREFIX pavector
new_zero_avector(uint dim);

/** @brief Create a new @ref avector object using a given array for
 *  the coefficients.
 *
 *  @param src Source array, should contain at least <tt>dim</tt> elements.
 *  @param dim Dimension of the new vector.
 *  @returns New @ref avector object. */
HEADER_PREFIX pavector
new_pointer_avector(pfield src, uint dim);

/** @brief Delete an @ref avector object.
 *
 *  Releases the storage corresponding to the object.
 *
 *  @attention Make sure that there are no submatrix objects referring
 *  to this object, since they might otherwise keep using pointers
 *  to invalid storage.
 *
 *  @param v Object to be deleted. */
HEADER_PREFIX void
del_avector(pavector v);

/** @brief Change the dimension of an @ref avector object without
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
resize_avector(pavector v, uint dim);

/** @brief Reduce the dimension of an @ref avector object without
 *  reallocating storage, preserving its coefficients.
 *
 *  @param v Vector to be resized.
 *  @param dim New vector dimension, not greater than <tt>v->dim</tt>. */
HEADER_PREFIX void
shrink_avector(pavector v, uint dim);

/* ------------------------------------------------------------
 Access methods
 ------------------------------------------------------------ */

#ifdef __GNUC__
INLINE_PREFIX field
getentry_avector(pcavector, uint) __attribute__((const,unused));
INLINE_PREFIX void
setentry_avector(pavector, uint, field) __attribute__((unused));
INLINE_PREFIX field
addentry_avector(pavector, uint, field) __attribute__((unused));
#endif

/** @brief Read a vector entry @f$v_i@f$.
 *
 *  @param v Vector @f$v@f$.
 *  @param i Index @f$i@f$.
 *  @returns Vector entry @f$v_i@f$. */
INLINE_PREFIX field getentry_avector(pcavector v, uint i) {
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
INLINE_PREFIX void setentry_avector(pavector v, uint i, field x) {
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
INLINE_PREFIX field addentry_avector(pavector v, uint i, field x) {
#ifdef FULL_DEBUG
  assert(i < v->dim);
#endif

  return (v->v[i] += x);
}

/* ------------------------------------------------------------
 Statistics
 ------------------------------------------------------------ */

/** @brief Get number of currently initialized @ref avector objects.
 *
 *  Calls to initialization functions like @ref init_avector and
 *  constructors like @ref new_avector increase an internal counter,
 *  while @ref uninit_avector and @ref del_avector decrease it.
 *
 *  @remark Use this function to check whether a program correctly cleans
 *  up temporary variables.
 *
 *  @returns Number of currently initialized @ref avector objects. */
HEADER_PREFIX uint
getactives_avector();

/** @brief Get size of a given @ref avector object.
 *
 *  Computes the size of the @ref avector object and the storage
 *  allocated for the coefficients.
 *  If the object uses the coefficients of another object
 *  (e.g., if it was created using @ref new_sub_avector), no coefficient
 *  storage is added.
 *
 *  @param v Vector object.
 *  @returns Size of allocated storage in bytes. */
HEADER_PREFIX size_t
getsize_avector(pcavector v);

/** @brief Get heap size of a given @ref avector object.
 *
 *  Computes the size of storage allocated for the coefficients
 *  on the heap, but not for the @ref avector object itself.
 *  If the object uses the coefficients of another object
 *  (e.g., if it was created using @ref new_sub_avector), no storage
 *  is required.
 *
 *  @param v Vector object.
 *  @returns Size of allocated heap storage in bytes. */
HEADER_PREFIX size_t
getsize_heap_avector(pcavector v);

/* ------------------------------------------------------------
 Simple utility functions
 ------------------------------------------------------------ */

/** @brief Set a vector to zero.
 *
 *  @param v Target vector. */
HEADER_PREFIX void
clear_avector(pavector v);

/** @brief Set all coefficients in a vector to the same value.
 *
 *  @param v Target vector.
 *  @param x Fill value. */
HEADER_PREFIX void
fill_avector(pavector v, field x);

/** @brief Fill a vector with random values.
 *
 *  @param v Target vector. */
HEADER_PREFIX void
random_avector(pavector v);

/** @brief Fill a vector with real valued random values.
 *
 *  @param v Target vector. */
HEADER_PREFIX void
random_real_avector(pavector v);

/** @brief Copy a vector into another vector, @f$w \gets v@f$.
 *
 *  Both vectors have to be of the same dimension.
 *
 *  @remark If you want to copy between vectors of different dimension,
 *  consider using @ref copy_sub_avector.
 *
 *  @param v Source vector.
 *  @param w Target vector. */
HEADER_PREFIX void
copy_avector(pcavector v, pavector w);

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
copy_sub_avector(pcavector v, pavector w);

/** @brief Print a vector.
 *
 *  @param v Vector @f$v@f$. */
HEADER_PREFIX void
print_avector(pcavector v);

/* ------------------------------------------------------------
 Very basic linear algebra
 ------------------------------------------------------------ */

/** @brief Scale a vector @f$v@f$ by a factor @f$\alpha@f$,
 *  @f$v \gets \alpha v@f$.
 *  
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param v Target vector @f$v@f$. */
HEADER_PREFIX void
scale_avector(field alpha, pavector v);

/** @brief Compute the Euclidean norm @f$\|v\|_2@f$ of a vector @f$v@f$.
 *
 *  @param v Vector @f$v@f$.
 *  @return Returns @f$\|v\|_2@f$.
 */
HEADER_PREFIX real
norm2_avector(pcavector v);

/** @brief Compute the Euclidean innner product @f$\langle x, y\rangle_2@f$
 *  of two vectors @f$x@f$ and @f$y@f$.
 *
 *  The Euclidean inner product is given by
 *  @f$\langle x, y \rangle_2 = \sum_i \bar x_i y_i@f$.
 *
 *  @param x Vector @f$x@f$.
 *  @param y Vector @f$x@f$.
 *  @returns Euclidean inner product @f$\langle x, y\rangle_2@f$. */
HEADER_PREFIX field
dotprod_avector(pcavector x, pcavector y);

/** @brief Add two vectors,
 *  @f$y \gets y + \alpha x@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
HEADER_PREFIX void
add_avector(field alpha, pcavector x, pavector y);

/** @} */

#endif
