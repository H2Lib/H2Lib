
/* ------------------------------------------------------------
 * This is the file "uniform.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

/** @file uniform.h
 *  @author Steffen B&ouml;rm */

#ifndef UNIFORM_H
#define UNIFORM_H

/** \defgroup uniform uniform
 *  @brief Representation of an admissible block for
 *  @f$ \mathcal H^2 @f$-matrices.
 *
 *  In case of @f$ \mathcal H^2@f$-matrices an admissible subblock is
 *  represented by
 *  @f[
 *  G|_{t \times s} = V_t S_b W_s^*,
 *  @f]
 *  where @f$ V_t, W_s @f$ are the corresponding clusterbasis.
 *
 *  A @ref _uniform "uniform" object stores the coupling-matrix @f$ S_b @f$
 *  aswell as references to the row- and column-clusterbasis.
 *  @{ */

/** @brief Representation of an admissible block for
 *  @f$ \mathcal H^2 @f$-matrices. */
typedef struct _uniform uniform;

/** @brief Pointer to @ref _uniform "uniform" object. */
typedef uniform *puniform;

/** @brief Pointer to a constant @ref _uniform "uniform" object. */
typedef const uniform *pcuniform;

/* C STD LIBRARY */
/* CORE 0 */
#include "settings.h"
/* CORE 1 */
#include "amatrix.h"
/* CORE 2 */
#include "clusterbasis.h"
#include "clusteroperator.h"
#include "rkmatrix.h"
/* CORE 3 */
/* SIMPLE */
/* PARTICLES */
/* BEM */

/** @brief Representation of an admissible block for
 *  @f$ \mathcal H^2 @f$-matrices.
 *
 *  <tt>rb</tt> and <tt>cb</tt> are references to the row- and column-
 *  @ref _clusterbasis "clusterbasis" respectively.
 *  The coupling-matrix is stored inside an @ref _amatrix "amatrix" <tt>S</tt>.
 *
 *  For convenience reasons there a two lists storing all @ref _uniform "uniform"
 *  objects belonging to the same block row defined by the row @ref _clusterbasis
 *  "clusterbasis" <tt>rb</tt> via <tt>rnext, rprev</tt> or to the same block
 *  column defined by the row @ref _clusterbasis "clusterbasis"
 *  <tt>cb</tt> via <tt>cnext, cprev</tt>.
 *
 *  If it is necessary to update these lists, please use the functions
 *  @ref ref_row_uniform, @ref ref_col_uniform,  @ref unref_row_uniform
 *  and @ref unref_col_uniform to perform this task.
 */
struct _uniform {
  /** @brief Row @ref _clusterbasis "clusterbasis" */
  pclusterbasis rb;
  /** @brief Column @ref _clusterbasis "clusterbasis" */
  pclusterbasis cb;
  /** @brief Coupling matrix */
  amatrix S;
  /** @brief Next row block in list */
  puniform rnext;
  /** @brief Previous row block in list */
  puniform rprev;
  /** @brief Next column block in list */
  puniform cnext;
  /** @brief Previous column block in list */
  puniform cprev;
};

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

/**
 * @brief Create a new @ref _uniform "uniform" object.
 *
 * Allocates storage for the object and sets the @ref _clusterbasis "clusterbasis"
 * pointers according to the input parameters.
 *
 * @remark Should always be matched by a call to @ref del_uniform.

 * @param rb Row @ref _clusterbasis "clusterbasis"
 * @param cb Column @ref _clusterbasis "clusterbasis"
 * @return Returns a pointer to a newly allocated @ref _uniform "uniform"
 * object.
 */
HEADER_PREFIX puniform
new_uniform(pclusterbasis rb, pclusterbasis cb);

/** @brief Deletes an @ref _uniform "uniform" object.
 *
 *  Releases the memory corresponding to the object.
 *  Reference inside the block row and block column lists are automatically
 *  removed by this function.
 *
 *  @param u Object to be deleted. */
HEADER_PREFIX void
del_uniform(puniform u);

/* ------------------------------------------------------------
 * Reference management
 * ------------------------------------------------------------ */

/** @brief Append this @ref _uniform "uniform" object to the block row list.
 *
 *  A reference to the @ref _uniform "uniform" object will be set into a
 *  reference list inside the row @ref _clusterbasis "clusterbasis" <tt>rb</tt>.
 *  Vice versa a reference to the current object will also be inserted into
 *  the <tt>rnext</tt> respectively <tt>rprev</tt> list.
 *
 *  @remark If one wants to remove a reference to this @ref _uniform "uniform"
 *  object, the function @ref unref_row_uniform should be called.
 *
 *  @param u @ref _uniform "Uniform" object to be referenced.
 *  @param rb Row @ref _clusterbasis "clusterbasis" to be referenced. */
HEADER_PREFIX void
ref_row_uniform(puniform u, pclusterbasis rb);

/** @brief Append this @ref _uniform "uniform" object to the block column list.
 *
 *  A reference to the @ref _uniform "uniform" object will be set into a
 *  reference list inside the column @ref _clusterbasis "clusterbasis" <tt>cb</tt>.
 *  Vice versa a reference to the current object will also be inserted into
 *  the <tt>cnext</tt> respectively <tt>cprev</tt> list.
 *
 *  @remark If one wants to remove a reference to this @ref _uniform "uniform"
 *  object, the function @ref unref_col_uniform should be called.
 *
 *  @param u @ref _uniform "Uniform" object to be referenced.
 *  @param cb Column @ref _clusterbasis "clusterbasis" to be referenced. */
HEADER_PREFIX void
ref_col_uniform(puniform u, pclusterbasis cb);

/**
 * @brief Remove block row references to an @ref _uniform "uniform" object.
 *
 * All previous made block row references to this @ref _uniform "uniform" object
 * will be removed. This is the counterpart of the function @ref ref_row_uniform.
 *
 * @param u @ref _uniform "Uniform" object for which the references should be
 * removed.
 */
HEADER_PREFIX void
unref_row_uniform(puniform u);

/**
 * @brief Remove block column references to an @ref _uniform "uniform" object.
 *
 * All previous made block column references to this @ref _uniform "uniform" object
 * will be removed. This is the counterpart of the function @ref ref_col_uniform.
 *
 * @param u @ref _uniform "Uniform" object for which the references should be
 * removed.
 */
HEADER_PREFIX void
unref_col_uniform(puniform u);

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

/** @brief Get size of a given @ref _uniform "uniform object.
 *
 *  @param u @ref _uniform "uniform" object.
 *  @returns Size of allocated storage in bytes. */
HEADER_PREFIX size_t
getsize_uniform(pcuniform u);

/* ------------------------------------------------------------
 * Simple utility functions
 * ------------------------------------------------------------ */

/**
 * @brief Sets the coupling matrix to zero.
 *
 * @param u @ref _uniform "uniform" object.
 */
HEADER_PREFIX void
clear_uniform(puniform u);

/** @brief Copy a uniform matrix @f$G@f$.
 *
 *  @param trans Set if @f$G^*@f$ should be copied instead of @f$G@f$.
 *  @param src Source matrix.
 *  @param trg target matrix. */
HEADER_PREFIX void
copy_uniform(bool trans, pcuniform src, puniform trg);

/** @brief Clone a uniform matrix.
 *
 *  @param src Source matrix.
 *  @returns Clone of the source matrix. */
HEADER_PREFIX puniform
clone_uniform(pcuniform src);

/** @brief Scale a uniform matrix.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param u Target matrix @f$S@f$, will be overwritten by @f$\alpha S@f$. */
HEADER_PREFIX void
scale_uniform(field alpha, puniform u);

/** @brief Fill a uniform matrix with random coefficients.
 *
 *  @param u Target matrix @f$S@f$. */
HEADER_PREFIX void
random_uniform(puniform u);

///**
// * @brief Computes the euclidean-norm of the coupling matrix and or the product
// * of weight matrices with the coupling matrix.
// *
// * Depending on given row- and column-@ref _clusteroperator "clusteroperators"
// * this functions computed one of the following four values via
// * vector iteration:
// * @f[
// * \lVert S_b \rVert_2, \quad \lVert C_t S_b \rVert_2, \quad
// * \lVert S_b C_s^* \rVert_2 \quad \text{or} \quad \lVert C_t S_b C_s^* \rVert_2
// * @f]
// *
// * @param u @ref _uniform "uniform" object.
// * @param rw Row @ref _clusteroperator "clusteroperator". If <tt>rw == NULL</tt>
// * then @f$ C_t @f$ is omitted in the computation of the norm.
// * @param cw Column @ref _clusteroperator "clusteroperator". If <tt>cw == NULL</tt>
// * then @f$ C_s @f$ is omitted in the computation of the norm.
// *
// * @return Approximation of the euclidean-norm of one of the above mentioned
// * matrices.
// */
//HEADER_PREFIX real
//norm2_fast_uniform(pcuniform u, pcclusteroperator rw, pcclusteroperator cw);
//
///**
// * @brief Computes the frobenius-norm of the coupling matrix and or the product
// * of weight matrices with the coupling matrix.
// *
// * Depending on given row- and column-@ref _clusteroperator "clusteroperators"
// * this functions computed one of the following four values:
// * @f[
// * \lVert S_b \rVert_F, \quad \lVert C_t S_b \rVert_F, \quad
// * \lVert S_b C_s^* \rVert_F \quad \text{or} \quad \lVert C_t S_b C_s^* \rVert_F
// * @f]
// *
// * @param u @ref _uniform "uniform" object.
// * @param rw Row @ref _clusteroperator "clusteroperator". If <tt>rw == NULL</tt>
// * then @f$ C_t @f$ is omitted in the computation of the norm.
// * @param cw Column @ref _clusteroperator "clusteroperator". If <tt>cw == NULL</tt>
// * then @f$ C_s @f$ is omitted in the computation of the norm.
// *
// * @return Approximation of the frobenius-norm of one of the above mentioned
// * matrices.
// */
//HEADER_PREFIX real
//normfrob_fast_uniform(pcuniform u, pcclusteroperator rw, pcclusteroperator cw);

/* ------------------------------------------------------------
 * Matrix-vector multiplication
 * ------------------------------------------------------------ */

/** @brief Matrix-vector multiplication
 *  @f$y \gets y + \alpha u x@f$ or @f$y \gets y + \alpha u^* x@f$.
 *
 *  The @ref _uniform "uniform" matrix or its adjoint is multiplied by the
 *  source vector @f$x@f$, the result is scaled by @f$\alpha@f$ and added
 *  to the target vector @f$y@f$. Here the @ref _uniform "uniform" matrix
 *  is represented by @f$ u = V_t S_b W_s^* @f$
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param trans Set if @f$u^*@f$ is to be used instead of @f$u@f$.
 *  @param u Matrix @f$A@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
HEADER_PREFIX void
mvm_uniform_avector(field alpha, bool trans, pcuniform u, pcavector x, pavector y);

/* ------------------------------------------------------------
 * Conversion operations
 * ------------------------------------------------------------ */

/**
 * @brief Addition of two @ref _uniform "uniform" matrices in different
 * @ref _clusterbasis "clusterbasis".
 *
 * Computes @f$  S_{b,\text{new}} \gets S_{b,\text{new}} + \hat S@f$, where
 * @f$ \hat S @f$ can be either
 * @f[
 * \hat S = S_b, \quad \hat S = C_t S_b, \quad \hat S = S_b C_s^* \quad \text{or}
 * \quad \hat S = C_t S_b C_s^*,
 * @f]
 * which depends on the equality of the row- and column @ref _clusterbasis
 * "clusterbasis".
 *
 * @param u Source @ref _uniform "uniform" matrix with the representation
 * @f$ u = V_t S_b W_s^*@f$-
 * @param ro @ref _clusteroperator "Clusteroperator" defining the basis-change
 * for the row @ref _clusterbasis "clusterbasis". Its weight matrix is
 * @f$ C_t @f$.
 * @param co @ref _clusteroperator "Clusteroperator" defining the basis-change
 * for the column @ref _clusterbasis "clusterbasis". Its weight matrix is
 * @f$ C_s @f$.
 * @param unew Source @ref _uniform "uniform" matrix with the representation
 * @f$ u_\text{new} = V_{t,\text{new}} S_{b,\text{new}} W_{s,\text{new}}^*@f$.
 */
HEADER_PREFIX void
add_projected_uniform(pcuniform u, pcclusteroperator ro, pcclusteroperator co,
    puniform unew);

/**
 * @brief Project @ref _uniform "uniform" matrix into new @ref _clusterbasis
 * "clusterbasis".
 *
 * The @ref _uniform "uniform" matrix <tt>u</tt> is represented by
 * @f$ u = V_t S_b W_S^* @f$. Depending on the equality of the old and the
 * new @ref _clusterbasis "clusterbasis" @f$ S_b @f$ is being updated to
 * @f[ S_b \gets S_b, \quad S_b \gets C_t S_b, \quad S_b \gets S_b C_s^*,
 * \quad \text{or} \quad S_b \gets C_t S_b C_s^*
 * @f]
 *
 * @param u @ref _uniform "uniform" object.
 * @param rb New row @ref _clusterbasis "clusterbasis" for <tt>u</tt>.
 * @param ro @ref _clusteroperator "Clusteroperator describing the basischange
 * between the old and new row @ref _clusterbasis "clusterbasis.
 * @param cb New column @ref _clusterbasis "clusterbasis" for <tt>u</tt>.
 * @param co @ref _clusteroperator "Clusteroperator describing the basischange
 * between the old and new column @ref _clusterbasis "clusterbasis.
 */
HEADER_PREFIX void
project_inplace_uniform(puniform u, pclusterbasis rb, pcclusteroperator ro,
    pclusterbasis cb, pcclusteroperator co);

/**
 * @brief Adds a @ref _rkmatrix "low rank matrix" to an @ref _uniform "uniform"
 * matrix.
 *
 * Updates the coupling matrix @f$ S_b @f$ of @f$ u @f$ to
 * @f$ S_b \gets S_b + V_t^* A B^* W_s @f$ for a @ref _rkmatrix "low rank matrix"
 * @f$ r = A B^* @f$ .
 *
 * @param r @ref _rkmatrix "low rank matrix" to be added to an @ref _uniform
 * "uniform" matrix.
 * @param unew @ref _uniform "uniform" matrix storing the result of the
 * computation.
 */
HEADER_PREFIX void
add_rkmatrix_uniform(pcrkmatrix r, puniform unew);

#endif

/* ------------------------------------------------------------
 * Access methods
 * ------------------------------------------------------------ */

/* Avoid problems with incomplete type definitions */
#if defined(CLUSTERBASIS_TYPE_COMPLETE) && !defined(UNIFORM_COMPLETE)
#define UNIFORM_COMPLETE

#ifdef __GNUC__
INLINE_PREFIX uint
getrows_uniform(pcuniform) __attribute__ ((const,unused));
INLINE_PREFIX uint
getcols_uniform(pcuniform) __attribute__ ((const,unused));
INLINE_PREFIX pamatrix
getS_uniform(puniform) __attribute__ ((unused));
#endif

/** @brief Get the number of rows of a @ref uniform matrix @f$G=V S W^*@f$.
 *
 *  @param u Matrix @f$G@f$.
 *  @return Number of rows of @f$G@f$, i.e., number of rows of @f$V@f$. */
INLINE_PREFIX uint
getrows_uniform(pcuniform u)
{
  return u->rb->t->size;
}

/** @brief Get the number of columns of a @ref uniform matrix @f$G=V S W^*@f$.
 *
 *  @param u Matrix @f$G@f$.
 *  @returns Number of columns of @f$G@f$, i.e., number of rows of @f$W@f$. */
INLINE_PREFIX uint
getcols_uniform(pcuniform u)
{
  return u->cb->t->size;
}

/** @brief Get the factor S of a @ref uniform matrix @f$G=V S W^*@f$.
 *
 * @param u Matrix @f$G@f$.
 * @returns Factor @f$S@f$. */
INLINE_PREFIX pamatrix
getS_uniform(puniform u)
{
  return &u->S;
}

#endif

/** @} */
