
/* ------------------------------------------------------------
 * This is the file "duniform.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2015
 * ------------------------------------------------------------ */

/** @file duniform.h
 *  @author Steffen B&ouml;rm */

#ifndef DUNIFORM_H
#define DUNIFORM_H

/** @defgroup duniform duniform
 *  @brief Uniform matrices with directional cluster bases.
 *  @{ */

/** @brief Uniform matrices with directional cluster bases. */
typedef struct _duniform duniform;

/** @brief Pointer to @ref duniform objects. */
typedef duniform *pduniform;

/** @brief Pointer to constant @ref duniform objects. */
typedef const duniform *pcduniform;

#include "settings.h"
#include "amatrix.h"
#include "dclusterbasis.h"

/** @brief Uniform matrices with directional cluster bases.
 *
 *  Matrices are represented in the form @f$V_{t c} S_{ts} W_{s c}^*@f$,
 *  where @f$V_{t c}@f$ and @f$W_{s c}@f$ are given by directional
 *  cluster bases for the row cluster @f$t@f$, the column cluster
 *  @f$s@f$, and the direction @f$c@f$. */
struct _duniform {
  /** @brief Row cluster basis */
  pdclusterbasis rb;
  /** @brief Column cluster basis */
  pdclusterbasis cb;

  /** @brief Direction for row cluster */
  uint rd;
  /** @brief Direction for column cluster */
  uint cd;

  /** @brief Coupling matrix */
  amatrix S;
};

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

/** @brief Create a new @ref duniform object.
 *
 *  @param rb Row cluster basis.
 *  @param rd Direction for row basis.
 *  @param cb Column cluster basis.
 *  @param cd Direction for column basis.
 *  @returns New object. */
HEADER_PREFIX pduniform
new_duniform(pdclusterbasis rb, uint rd, pdclusterbasis cb, uint cd);

/** @brief Delete a @ref duniform object.
 *
 *  @param u Object to be deleted. */
HEADER_PREFIX void
del_duniform(pduniform u);

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

/** @brief Compute the storage size of a @ref duniform object.
 *
 *  @param u Object.
 *  @returns Storage size of <tt>u</tt> in bytes,
 *    not including the cluster bases. */
HEADER_PREFIX size_t
getsize_duniform(pcduniform u);

/* ------------------------------------------------------------
 * Simple utility functions
 * ------------------------------------------------------------ */

/** @brief Set a @ref duniform matrix to zero.
 *
 *  @param u Matrix. */
HEADER_PREFIX void
clear_duniform(pduniform u);

/** @brief Copy a @ref duniform matrix.
 *
 *  @param src Source matrix.
 *  @param trg Target matrix. */
HEADER_PREFIX void
copy_duniform(pcduniform src, pduniform trg);

/* ------------------------------------------------------------
 * Matrix-vector multiplication
 * ------------------------------------------------------------ */

/** @brief Fast matrix-vector multiplication.
 *
 *  Multiply the transformed coefficients in <tt>xt</tt> by
 *  the coupling matrix and add the result to <tt>yt</tt>.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param u Matrix.
 *  @param xt Transformed input vector as provided, e.g., by
 *    @ref forward_dclusterbasis. Size <tt>u->cb->ktree</tt>.
 *  @param yt Transformed output vector, functions like
 *    @ref backward_dclusterbasis can be used to obtain
 *    the final result. Size <tt>u->rb->ktree</tt>. */
HEADER_PREFIX void
fastaddeval_duniform_avector(field alpha, pcduniform u,
			     pcavector xt, pavector yt);

/** @brief Fast adjoint matrix-vector multiplication.
 *
 *  Multiply the transformed coefficients in <tt>xt</tt> by
 *  the adjoint coupling matrix and add the result to <tt>yt</tt>.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param u Matrix.
 *  @param xt Transformed input vector as provided, e.g., by
 *    @ref forward_dclusterbasis. Size <tt>u->rb->ktree</tt>
 *  @param yt Transformed output vector, functions like
 *    @ref backward_dclusterbasis can be used to obtain
 *    the final result. Size <tt>u->cb->ktree</tt>. */
HEADER_PREFIX void
fastaddevaltrans_duniform_avector(field alpha, pcduniform u,
				  pcavector xt, pavector yt);

/** @brief Slow matrix-vector multiplication @f$y \gets y + \alpha G x@f$.
 *
 *  Intended for debugging, applies forward and backward transformation
 *  for each individual block instead of sharing transformed vectors.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param u matrix.
 *  @param x Input vector @f$x@f$.
 *  @param y Output vector @f$y@f$. */
HEADER_PREFIX void
slowaddeval_duniform_avector(field alpha, pcduniform u,
			     pcavector x, pavector y);

/** @brief Slow adjoint matrix-vector multiplication
 *    @f$y \gets y + \alpha G^* x@f$.
 *
 *  Intended for debugging, applies forward and backward transformation
 *  for each individual block instead of sharing transformed vectors.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param u matrix.
 *  @param x Input vector @f$x@f$.
 *  @param y Output vector @f$y@f$. */
HEADER_PREFIX void
slowaddevaltrans_duniform_avector(field alpha, pcduniform u,
				  pcavector x, pavector y);

/* ------------------------------------------------------------
 * Conversion to a full matrix
 * ------------------------------------------------------------ */

/** @brief Convert a @ref duniform matrix to a dense matrix.
 *
 *  Adds @f$\alpha V_{t c} S_{ts} W_{s c}^*@f$ to a matrix @f$G@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param u Source matrix.
 *  @param G Target matrix. */
HEADER_PREFIX void
expand_duniform(field alpha, pcduniform u, pamatrix G);

/* ------------------------------------------------------------
 * Matrix norm 
 * ------------------------------------------------------------ */

/**
 * @brief Will compute the Euclidean norm of the coupling matrix or of the
 * product of one or two weight matrices with the coupling matrix.
 * 
 * If no directional cluster operator is given, the function will compute
 * @f$ \| S_{b} \|_{2} @f$.
 * If a directional row @f$ ro @f$ or column cluster operators @f$ co @f$ are given
 * @f[ \| C_{t, rd} S_{b} \|_{2}, \| S_{b} C_{s, cd}^{*} \|_{2} \text{ or } 
 * \| C_{t, rd} S_{b} C_{s, cd}^{*} \|_{2} ,@f]
 *  where @f$rd, cd @f$ denotes the directions corresponding to @f$ S_{b} @f$,
 *  will be computed.
 *
 * @param u @ref duniform object for the coupling matrix.
 * @param ro Optional directional row cluster operator.
 * If <tt>ro == NULL</tt>, no @f$ C_{t, rd} @f$ will be used.
 * @param co Optional directional column cluster operator.
 * If <tt>c0 == NULL</tt>, no @f$ C_{s, cd} @f$ will be used.
 *
 * @return Approximation of the Euclidean norm of the coupling matrix or of the
 * product of one or two weight matrices with the coupling matrix.
 */

HEADER_PREFIX real
norm2_fast_duniform(pcduniform u, pcdclusteroperator ro, pcdclusteroperator co);

/**
 * @brief Will compute the Frobenius norm of the coupling matrix or of the
 * product of one or two weight matrices with the coupling matrix.
 * 
 * If no directional cluster operator is given, the function will compute
 * @f$ \| S_{b} \|_{2} @f$.
 * If directional row @f$ ro @f$ or column cluster operators @f$ co @f$ are given
 * @f[ \| C_{t, rd} S_{b} \|_{2}, \| S_{b} C_{s, cd}^{*} \|_{2} \text{ or } 
 * \| C_{t, rd} S_{b} C_{s, cd}^{*} \|_{2} ,@f]
 * where @f$rd, cd @f$ denotes the directions corresponding to @f$ S_{b} @f$,
 * will be computed.
 *
 * @param u @ref duniform object for the coupling matrix.
 * @param ro Optional directional row cluster operator.
 * If <tt>ro == NULL</tt>, no @f$ C_{t, rd} @f$ will be used.
 * @param co Optional directional column cluster operator.
 * If <tt>co == NULL</tt>, no @f$ C_{s, cd} @f$ will be used.
 *
 * @return Approximation of the Frbenius norm of the coupling matrix or of the
 * product of one or two weight matrices with the coupling matrix.
 */


HEADER_PREFIX real
normfrob_fast_duniform(pcduniform u, pcdclusteroperator ro, pcdclusteroperator co);


/* ---------------------------------------------------------- 
 * Projection
 * ---------------------------------------------------------- */

/**
 * @brief Addition of two @ref duniform matrices even with different
 * directional cluster basis.
 *
 * Add to the @ref duniform matrix @f$ S_{b,\text{new}} @f$ another @ref duniform
 * matrix @f$ \tilde S @f$.
 * Where @f$ \tilde S @f$ depends on the equality of the row- and column basis.
 * The matrix can be either
 * @f[
 * \tilde S = S_{b}, \quad \tilde S = C_{t, rd} S_{b} \quad \tilde S = S_{b} C_{s, cd}^{*}
 * \text{ or } \quad \tilde S = C_{t, rd} S_{b} C_{s, cd}^{*}, @f]
 * where @f$rd, cd @f$ denotes the directions corresponding to @f$ S_{b} @f$.
 *
 * @param u Source @ref duniform matrix with corresponding @f$ S_b @f$.
 * @param ro @ref dclusteroperator describes the basis-change for the directional
 * row cluster basis with corresponding weight matrix @f$ C_{t, rd} @f$.
 * @param co @ref dclusteroperator describes the basis-change for the directional
 * column cluster basis with corresponding weight matrix @f$ C_{s, cd} @f$.
 * @param unew Source @ref duniform matrix with corresponding @f$ S_{b, \text{new}} @f$.
 */

HEADER_PREFIX void
add_projected_duniform(pcduniform u, pcdclusteroperator ro, pcdclusteroperator co, pduniform unew);

/** @} */

#endif
