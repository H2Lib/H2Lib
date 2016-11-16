/* ------------------------------------------------------------
 This is the file "h2compression.h" of the H2Lib package.
 All rights reserved, Steffen Boerm 2011
 ------------------------------------------------------------ */

/** @file h2compression.h
 *  @author Steffen B&ouml;rm
 */

#ifndef H2COMPRESSION_H
#define H2COMPRESSION_H

/** @defgroup h2compression h2compression
 *  @brief Functions for turning an @ref amatrix, @ref hmatrix
 *  or @ref h2matrix into an @ref h2matrix.
 *
 *  The functions in this module can be used to convert array
 *  matrices and hierarchical matrices into @f$\mathcal{H}^2@f$-matrices
 *  and to recompress @f$\mathcal{H}^2@f$-matrices.
 *  The module also provides a unification algorithm that takes
 *  a block matrix consisting of independent @f$\mathcal{H}^2@f$-matrices
 *  and turns it into an @f$\mathcal{H}^2@f$-matrix with unified
 *  row and column bases.
 *  @{ */

/** @brief Cluster basis and weights for use in the unification algorithm */
typedef struct _truncblock truncblock;

/** @brief Pointer to @ref truncblock object. */
typedef truncblock *ptruncblock;

#include "h2matrix.h"
#include "hmatrix.h"
#include "truncation.h"

/* ------------------------------------------------------------
 High-level compression functions
 ------------------------------------------------------------ */

/** @brief Approximate a matrix, represented by an @ref amatrix object,
 *  by an @f$\mathcal{H}^2@f$-matrix.
 *
 *  @param G Source matrix @f$G@f$.
 *  @param b Block tree, rows and columns have to match @f$G@f$.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy.
 *  @returns @f$\mathcal{H}^2@f$-matrix approximation of @f$G@f$. */
HEADER_PREFIX ph2matrix
compress_amatrix_h2matrix(pcamatrix G, pcblock b, pctruncmode tm, real eps);

/** @brief Approximate a hierarchical matrix, represented by an
 *  @ref hmatrix object, by an @f$\mathcal{H}^2@f$-matrix.
 *
 *  @param G Source matrix @f$G@f$.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy.
 *  @returns @f$\mathcal{H}^2@f$-matrix approximation of @f$G@f$. */
HEADER_PREFIX ph2matrix
compress_hmatrix_h2matrix(pchmatrix G, pctruncmode tm, real eps);

/** @brief Approximate an @f$\mathcal{H}^2@f$-matrix, represented by an
 *  @ref h2matrix object, by a recompressed @f$\mathcal{H}^2@f$-matrix.
 *
 *  @param G Source matrix @f$G@f$.
 *  @param rbortho Set if the original row basis is orthogonal,
 *    this allows the algorithm to avoid computing row weights.
 *  @param cbortho Set if the original column basis is orthogonal,
 *    this allows the algorithm to avoid computing column weights.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy.
 *  @returns @f$\mathcal{H}^2@f$-matrix approximation of @f$G@f$. */
HEADER_PREFIX ph2matrix
compress_h2matrix_h2matrix(pch2matrix G, bool rbortho, bool cbortho,
    pctruncmode tm, real eps);

/** @brief Approximate a symmetric @f$\mathcal{H}^2@f$-matrix, represented
 *  by an @ref h2matrix object, by a recompressed @f$\mathcal{H}^2@f$-matrix.
 *
 *  Since the matrix is symmetric, we assume the row and column cluster basis
 *  are identical, i.e., <tt>G->rb==G->cb</tt>.
 *
 *  @param G Source matrix @f$G@f$.
 *  @param rbortho Set if the original row and column basis is orthogonal,
 *    this allows the algorithm to avoid computing row weights.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy.
 *  @returns @f$\mathcal{H}^2@f$-matrix approximation of @f$G@f$. */
HEADER_PREFIX ph2matrix
compress_symmetric_h2matrix_h2matrix(pch2matrix G, bool rbortho,
    pctruncmode tm, real eps);

/* ------------------------------------------------------------
 Compute local and total weights for H^2-matrices
 ------------------------------------------------------------ */

/** @brief Prepare local row weights for a given @f$\mathcal{H}^2@f$-matrix.
 *
 *  Finds matrices @f$Z^+_t@f$ with
 *  @f$\begin{pmatrix} G_{t,s_1} & \cdots & G_{t,s_\sigma} \end{pmatrix}
 *  = V_t (Z^+_t)^* P_t^*@f$, where @f$(t,s_1),\ldots,(t,s_\sigma)@f$ are the
 *  admissible blocks connected to the row cluster @f$t@f$ and @f$P_t@f$
 *  is orthogonal.
 *
 *  @remark If row and column weights have to be computed, consider
 *  using @ref localweights_h2matrix.
 *
 *  @param G Source matrix @f$G@f$.
 *  @param rbw Basis weights for the row basis, can be obtained by
 *    @ref weight_clusterbasis_clusteroperator.
 *    If this is a null pointer, the row basis is assumed to be
 *    orthogonal, so no weights are required.
 *  @param cbw Basis weights for the column basis, can be obtained by
 *    @ref weight_clusterbasis_clusteroperator.
 *    If this is a null pointer, the column basis is assumed to be
 *    orthogonal, so no weights are required.
 *  @param tm Truncation mode.
 *  @param rlw Will be filled with the local row weights, can be
 *    initialized by @ref build_from_clusterbasis_clusteroperator.
 *    If the cluster operator already contains something, the local
 *    weights will be merged with the original contents. */
HEADER_PREFIX void
rowweights_h2matrix(pch2matrix G, pcclusteroperator rbw, pcclusteroperator cbw,
    pctruncmode tm, pclusteroperator rlw);

/** @brief Prepare local column weights for a given @f$\mathcal{H}^2@f$-matrix.
 *
 *  Finds matrices @f$Z^+_s@f$ with
 *  @f$\begin{pmatrix} G_{t_1,s}\\ \vdots\\ G_{t_\tau,s} \end{pmatrix}
 *  = P_s Z^+_s W_s^*@f$, where @f$(t_1,s),\ldots,(t_\tau,s)@f$ are the
 *  admissible blocks connected to the column cluster @f$s@f$ and @f$P_s@f$
 *  is orthogonal.
 *
 *  @remark If row and column weights have to be computed, consider
 *  using @ref localweights_h2matrix.
 *
 *  @param G Source matrix @f$G@f$.
 *  @param rbw Basis weights for the row basis, can be obtained by
 *    @ref weight_clusterbasis_clusteroperator.
 *    If this is a null pointer, the row basis is assumed to be
 *    orthogonal, so no weights are required.
 *  @param cbw Basis weights for the column basis, can be obtained by
 *    @ref weight_clusterbasis_clusteroperator.
 *    If this is a null pointer, the column basis is assumed to be
 *    orthogonal, so no weights are required.
 *  @param tm Truncation mode.
 *  @param clw Will be filled with the local column weights, can be
 *    initialized by @ref build_from_clusterbasis_clusteroperator.
 *    If the cluster operator already contains something, the local
 *    weights will be merged with the original contents. */
HEADER_PREFIX void
colweights_h2matrix(pch2matrix G, pcclusteroperator rbw, pcclusteroperator cbw,
    pctruncmode tm, pclusteroperator clw);

/** @brief Prepare local weights for a given @f$\mathcal{H}^2@f$-matrix.
 *
 *  Finds matrices @f$Z^+_t@f$ with
 *  @f$\begin{pmatrix} G_{t,s_1} & \cdots & G_{t,s_\sigma} \end{pmatrix}
 *  = V_t (Z^+_t)^* P_t@f$, where @f$(t,s_1),\ldots,(t,s_\sigma)@f$ are the
 *  admissible blocks connected to the row cluster @f$t@f$ and @f$P_t@f$
 *  is orthogonal
 *  and matrices @f$Z^+_s@f$ with
 *  @f$\begin{pmatrix} G_{t_1,s}\\ \vdots\\ G_{t_\tau,s} \end{pmatrix}
 *  = Z^+_s W_s^* P_s@f$, where @f$(t_1,s),\ldots,(t_\tau,s)@f$ are the
 *  admissible blocks connected to the column cluster @f$s@f$ and @f$P_s@f$
 *  is orthogonal.
 *
 *  This algorithm is more efficient than calling @ref rowweights_h2matrix
 *  and @ref colweights_h2matrix individually, since it passes only
 *  once through the matrix.
 *
 *  @param G Source matrix @f$G@f$.
 *  @param rbw Basis weights for the row basis, can be obtained by
 *    @ref weight_clusterbasis_clusteroperator.
 *    If this is a null pointer, the row basis is assumed to be
 *    orthogonal, so no weights are required.
 *  @param cbw Basis weights for the column basis, can be obtained by
 *    @ref weight_clusterbasis_clusteroperator.
 *    If this is a null pointer, the column basis is assumed to be
 *    orthogonal, so no weights are required.
 *  @param tm Truncation mode.
 *  @param rlw Will be filled with the local row weights, can be
 *    initialized by @ref build_from_clusterbasis_clusteroperator.
 *    If the cluster operator already contains something, the local
 *    weights will be merged with the original contents.
 *  @param clw Will be filled with the local column weights, can be
 *    initialized by @ref build_from_clusterbasis_clusteroperator.
 *    If the cluster operator already contains something, the local
 *    weights will be merged with the original contents. */
HEADER_PREFIX void
localweights_h2matrix(pch2matrix G, pcclusteroperator rbw,
    pcclusteroperator cbw, pctruncmode tm, pclusteroperator rlw,
    pclusteroperator clw);

/** @brief Merge the local weights of clusters with the weights
 *  inherited from their ancestors in order to obtain total weights.
 *
 *  Finds matrices @f$Z_t@f$ with @f$P_t Z_t = \begin{pmatrix}
 *  Z^+_t\\ \zeta_{\rm age} Z_{t^+} E_t^* \end{pmatrix}@f$, where @f$t^+@f$
 *  is the father of the cluster @f$t@f$ and @f$P_t@f$ is orthogonal.
 *
 *  @param cb Cluster basis.
 *  @param tm Truncation mode, supplies @f$\zeta_{\rm age}@f$.
 *  @param lw Cluster operator, contains @f$Z^+_t@f$ and is overwritten
 *    by @f$Z_t@f$. */
HEADER_PREFIX void
accumulate_clusteroperator(pcclusterbasis cb, pctruncmode tm,
    pclusteroperator lw);

/** @brief Construct total weights for a given @f$\mathcal{H}^2@f$-matrix.
 *
 *  Finds matrices @f$Z_t@f$ with @f$G_t = V_t Z_t^* P_t^*@f$, where
 *  @f$G_t@f$ is the total cluster basis and @f$P_t@f$ is orthogonal.
 *  These are called total row weights.
 *  The total column weights are computed by applying the same procedure
 *  to @f$G^*@f$ instead of @f$G@f$.
 *
 *  @param G Source matrix @f$G@f$.
 *  @param rbortho Set if the original row basis is orthogonal,
 *    this allows the algorithm to avoid computing row weights.
 *  @param cbortho Set if the original column basis is orthogonal,
 *    this allows the algorithm to avoid computing column weights.
 *  @param tm Truncation mode.
 *  @param rw Total row weights, can be initialized by calling
 *    @ref build_from_clusterbasis_clusteroperator for the row basis.
 *  @param cw Total column weights, can be initialized by
 *    @ref build_from_clusterbasis_clusteroperator for the column basis. */
HEADER_PREFIX void
totalweights_h2matrix(pch2matrix G, bool rbortho, bool cbortho, pctruncmode tm,
    pclusteroperator rw, pclusteroperator cw);

/* ------------------------------------------------------------
 Compute truncated cluster basis
 ------------------------------------------------------------ */

/** @brief Compute a truncated cluster basis.
 *
 *  Reduces the rank of a cluster basis by computing singular value
 *  decompositions and choosing the most important left singular vectors
 *  to construct a new basis.
 *  Both total and local weights can be taken into account, and the
 *  local weights will be accumulated on the fly, making a call to
 *  @ref accumulate_clusteroperator unnecessary.
 *
 *  @param cb Original cluster basis.
 *  @param cw Total weights, ignored if null pointer.
 *  @param clw Local weights, will be accumulated, but not overwritten
 *    during the course of the algorithm, ignored if null pointer.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy.
 *  @param cbnew New cluster basis, can be initialized with
 *    @ref clonestructure_clusterbasis.
 *  @param old2new Cluster operator describing the change of basis
 *    from <tt>cb</tt> to <tt>cbnew</tt>. */
HEADER_PREFIX void
truncate_clusterbasis(pcclusterbasis cb, pcclusteroperator cw,
    pcclusteroperator clw, pctruncmode tm, real eps, pclusterbasis cbnew,
    pclusteroperator old2new);

/* ------------------------------------------------------------
 Compute adaptive cluster bases for an H^2-matrix
 ------------------------------------------------------------ */

/** @brief Construct an improved row basis for an @f$\mathcal{H}^2@f$-matrix.
 *
 *  @param G Original matrix @f$G@f$.
 *  @param rbortho Set if the original row basis is orthogonal,
 *    this allows the algorithm to avoid computing row weights.
 *  @param cbortho Set if the original column basis is orthogonal,
 *    this allows the algorithm to avoid computing column weights.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy.
 *  @param old2new Cluster operator describing the change of basis
 *    from <tt>G->rb</tt> to the new basis.
 *  @returns New row cluster basis. */
HEADER_PREFIX pclusterbasis
buildrowbasis_h2matrix(pch2matrix G, bool rbortho, bool cbortho, pctruncmode tm,
    real eps, pclusteroperator old2new);

/** @brief Construct an improved row basis for an @f$\mathcal{H}^2@f$-matrix.
 *
 *  @param G Original matrix @f$G@f$.
 *  @param rbortho Set if the original row basis is orthogonal,
 *    this allows the algorithm to avoid computing row weights.
 *  @param cbortho Set if the original column basis is orthogonal,
 *    this allows the algorithm to avoid computing column weights.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy.
 *  @param old2new Cluster operator describing the change of basis
 *    from <tt>G->cb</tt> to the new basis.
 *  @returns New column cluster basis. */
HEADER_PREFIX pclusterbasis
buildcolbasis_h2matrix(pch2matrix G, bool rbortho, bool cbortho, pctruncmode tm,
    real eps, pclusteroperator old2new);

/** @brief Construct an improved row basis for an @f$\mathcal{H}^2@f$-matrix.
 *
 *  @remark Compare to @ref buildrowbasis_h2matrix, this function
 *  constructs local weights by only one QR decomposition per cluster
 *  instead of one per block. It also handles block weights in a
 *  slightly different manner.
 *
 *  @param G Original matrix @f$G@f$.
 *  @param rbortho Set if the original row basis is orthogonal,
 *    this allows the algorithm to avoid computing row weights.
 *  @param cbortho Set if the original column basis is orthogonal,
 *    this allows the algorithm to avoid computing column weights.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy.
 *  @param rbnew New cluster basis, can be initialized with
 *    @ref clonestructure_clusterbasis.
 *  @param old2new Cluster operator describing the change of basis
 *    from <tt>G->rb</tt> to <tt>cbnew</tt>. */
HEADER_PREFIX void
truncrowbasis_h2matrix(ph2matrix G, bool rbortho, bool cbortho, pctruncmode tm,
    real eps, pclusterbasis rbnew, pclusteroperator old2new);

/** @brief Construct an improved column basis for an @f$\mathcal{H}^2@f$-matrix.
 *
 *  @remark Compare to @ref buildcolbasis_h2matrix, this function
 *  constructs local weights by only one QR decomposition per cluster
 *  instead of one per block. It also handles block weights in a
 *  slightly different manner.
 *
 *  @param G Original matrix @f$G@f$.
 *  @param rbortho Set if the original row basis is orthogonal,
 *    this allows the algorithm to avoid computing row weights.
 *  @param cbortho Set if the original column basis is orthogonal,
 *    this allows the algorithm to avoid computing column weights.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy.
 *  @param cbnew New cluster basis, can be initialized with
 *    @ref clonestructure_clusterbasis.
 *  @param old2new Cluster operator describing the change of basis
 *    from <tt>G->cb</tt> to <tt>cbnew</tt>. */
HEADER_PREFIX void
trunccolbasis_h2matrix(ph2matrix G, bool rbortho, bool cbortho, pctruncmode tm,
    real eps, pclusterbasis cbnew, pclusteroperator old2new);

/* ------------------------------------------------------------
 Approximate H^2-matrix in new cluster bases
 ------------------------------------------------------------ */

/** @brief Construct an @f$\mathcal{H}^2@f$-matrix approximation
 *  of a given @f$\mathcal{H}^2@f$-matrix in new cluster bases
 *  by blockwise projection.
 *
 *  @param G Original matrix.
 *  @param rb New row basis.
 *  @param ro Basis change from old row basis <tt>G->rb</tt> to
 *    new basis <tt>rb</tt>.
 *  @param cb New column basis.
 *  @param co Basis change from old column basis <tt>G->cb</tt> to
 *    new basis <tt>cb</tt>.
 *  @returns Approximating @f$\mathcal{H}^2@f$-matrix. */
HEADER_PREFIX ph2matrix
build_projected_h2matrix(pch2matrix G, pclusterbasis rb, pcclusteroperator ro,
    pclusterbasis cb, pcclusteroperator co);

/** @brief Switch the cluster bases of an @f$\mathcal{H}^2@f$-matrix
 *  by applying blockwise projections.
 *  @param G Original matrix, will be overwritten by projected matrix.
 *  @param rb New row basis.
 *  @param ro Basis change from old row basis <tt>G->rb</tt> to
 *    new basis <tt>rb</tt>.
 *  @param cb New column basis.
 *  @param co Basis change from old column basis <tt>G->cb</tt> to
 *    new basis <tt>cb</tt> */
HEADER_PREFIX void
project_inplace_h2matrix(ph2matrix G, pclusterbasis rb, pcclusteroperator ro,
    pclusterbasis cb, pcclusteroperator co);

/* ------------------------------------------------------------
 Specialized recompression routines for H^2-matrix arithmetic algorithms
 ------------------------------------------------------------ */

/** @brief Compute cluster weights and store them in the cluster basis.
 *
 *  The weights are the factors @f$R_t@f$ of the skinny QR decomposition
 *  @f$V_t = Q_t R_t@f$.
 *
 *  This function is similar to @ref weight_clusterbasis_clusteroperator,
 *  but stores the basis weights in <tt>cb->Z</tt> instead of in a
 *  separate @ref clusteroperator.
 *
 *  @param cb Target cluster basis, its field <tt>Z</tt> will be overwritten. */
HEADER_PREFIX void
orthoweight_clusterbasis(pclusterbasis cb);

/** @brief Compute total row weights of a matrix.
 *
 *  The row weights are the factor @f$R_t@f$ of the skinny QR decomposition
 *  @f$X_t = Q_t R_t@f$, where @f$X_t@f$ is the total row cluster basis
 *  of an @f$\mathcal{H}^2@f$-matrix @f$X@f$.
 *
 *  The matrix @f$X_t@f$ is obtained implicitly from the list
 *  starting with the pointer <tt>rlist</tt> contained in the
 *  cluster basis.
 *
 *  @param rb Row basis.
 *  @param rw Total row weight.
 *  @param tm Truncation mode. */
HEADER_PREFIX void
totalweight_row_clusteroperator(pclusterbasis rb, pclusteroperator rw,
    pctruncmode tm);

/** @brief Compute total column weights of a matrix.
 *
 *  The column weights are the factor @f$R_t@f$ of the skinny QR decomposition
 *  @f$X_t = Q_t R_t@f$, where @f$X_t@f$ is the total column cluster basis
 *  of an @f$\mathcal{H}^2@f$-matrix @f$X@f$.
 *
 *  The matrix @f$X_t@f$ is obtained implicitly from the list
 *  starting with the pointer <tt>clist</tt> contained in the
 *  cluster basis.
 *
 *  @param cb Column basis.
 *  @param cw Total column weight.
 *  @param tm Truncation mode. */
HEADER_PREFIX void
totalweight_col_clusteroperator(pclusterbasis cb, pclusteroperator cw,
    pctruncmode tm);

/** @brief Replace a cluster basis with a truncated basis.
 *
 *  This function is similar to @ref truncate_clusterbasis, but
 *  replaces the old cluster basis by the new one.
 *
 *  @param cb Old clusterbasis, will be overwritten by truncated cluster basis.
 *  @param cw Cluster weights.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy. */
HEADER_PREFIX void
truncate_inplace_clusterbasis(pclusterbasis cb, pclusteroperator cw,
    pctruncmode tm, real eps);

/** @brief Recompress an @f$\mathcal{H}^2@f$-matrix.
 *
 *  Replace an @f$\mathcal{H}^2@f$-matrix by a recompressed matrix.
 *
 *  @param G Original matrix, will be overwritten by the recompressed
 *    matrix.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy. */
HEADER_PREFIX void
recompress_inplace_h2matrix(ph2matrix G, pctruncmode tm, real eps);

/* ------------------------------------------------------------
 Unification
 ------------------------------------------------------------ */

/** @brief Description of basis and weights for one submatrix in the
 *  unification algorithm. */
struct _truncblock {
  /** @brief Cluster basis of the submatrix. */
  pcclusterbasis cb;
  /** @brief Total weights for <tt>cb</tt>. */
  pcclusteroperator cw;
  /** @brief Optional scaling factor for <tt>cw</tt>, required to
   *  apply aging factors to partial cluster bases. */
  field cw_factor;

  /** @brief Will be filled by basis change operator from <tt>cb</tt>
   *  to the new unified basis. */
  pclusteroperator old2new;

  /** @brief Temporary storage, used to handle partial cluster bases. */
  clusterbasis tmp_cb;
  /** @brief Temporary storage, used to handle partial cluster bases. */
  clusteroperator tmp_cw;

  /** @brief Next @ref truncblock in list. */
  ptruncblock next;
};

/** @brief Create a new @ref truncblock object.
 *
 *  @param cb Cluster basis of the corresponding submatrix.
 *  @param cw Total weights for <tt>cb</tt>.
 *  @param next Next @ref truncblock in list.
 *  @return Returns the newly created @ref truncblock object.
 */
HEADER_PREFIX ptruncblock
new_truncblock(pcclusterbasis cb, pcclusteroperator cw, ptruncblock next);

/** @brief Delete a list of @ref truncblock objects.
 *
 *  @param tb Head of list. */
HEADER_PREFIX void
del_truncblock(ptruncblock tb);

/** @brief Reverse the order of a list of @ref truncblock objects.
 *
 *  Since new blocks are added at the head of a list, reversing the
 *  list restores the original order.
 *
 *  @param tb Head of original list, will be changed.
 *  @returns Head of reversed list. */
HEADER_PREFIX ptruncblock
reverse_truncblock(ptruncblock tb);

/** @brief Construct a unified cluster basis.
 *
 *  Constructs a cluster basis that can approximate multiple
 *  @f$\mathcal{H}^2@f$-matrices simultaneously, e.g., multiple
 *  independent sub-@f$\mathcal{H}^2@f$-matrices that have to be
 *  merged into a large @f$\mathcal{H}^2@f$-matrix.
 *
 *  @param t Root of the corresponding cluster tree.
 *  @param tb List of @ref truncblock objects describing cluster bases
 *    and weights. After the function has completed, the fields
 *    <tt>old2new</tt> will contain the basis change operators.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy.
 *  @param cw If not a null pointer, will be set to point to a new
 *    @ref clusteroperator tree containing the total weights for the
 *    new basis. This is useful for recursive unification.
 *  @returns Unified cluster basis. */
HEADER_PREFIX pclusterbasis
unify_clusterbasis(pccluster t, ptruncblock tb, pctruncmode tm, real eps,
    pclusteroperator *cw);

/** @brief Unify @f$\mathcal{H}^2@f$-submatrices into a large
 *  @f$\mathcal{H}^2@f$-matrix.
 *
 *  Takes a block matrix containing @f$\mathcal{H}^2@f$-matrices
 *  and approximates it by a global @f$\mathcal{H}^2@f$-matrix.
 *
 *  @remark Differently from everywhere else in the library,
 *  <tt>G</tt> is <em>not</em> a proper @f$\mathcal{H}^2@f$-matrix,
 *  since the cluster bases of its immediate submatrices are allowed
 *  to differ from each other.
 *
 *  @param G Block matrix, will be overwritten by a proper
 *    @f$\mathcal{H}^2@f$-matrix approximation.
 *  @param rw1 Total row weights for all submatrices, enumerated
 *    in column-major ordering, i.e., <tt>rw1[i+j*G->rsons]</tt>
 *    corresponds to the block in row <tt>i</tt> and column <tt>j</tt>.
 *  @param cw1 Total column weights for all submatrices, enumerated
 *    in column-major ordering, i.e., <tt>rw1[i+j*G->rsons]</tt>
 *    corresponds to the block in row <tt>i</tt> and column <tt>j</tt>.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy.
 *  @param rw If not a null pointer, will be set to point to a new
 *    @ref clusteroperator tree containing the total row weights
 *    of the new @f$\mathcal{H}^2@f$-matrix.
 *  @param cw If not a null pointer, will be set to point to a new
 *    @ref clusteroperator tree containing the total column weights
 *    of the new @f$\mathcal{H}^2@f$-matrix. */
HEADER_PREFIX void
unify_h2matrix(ph2matrix G, pclusteroperator *rw1, pclusteroperator *cw1,
    pctruncmode tm, real eps, pclusteroperator *rw, pclusteroperator *cw);

/** @brief Unify @f$\mathcal{H}^2@f$-submatrices into a large
 *  @f$\mathcal{H}^2@f$-matrix, experimental parallel implementation.
 *
 *  Takes a block matrix containing @f$\mathcal{H}^2@f$-matrices
 *  and approximates it by a global @f$\mathcal{H}^2@f$-matrix.
 *
 *  @remark Differently from everywhere else in the library,
 *  <tt>G</tt> is <em>not</em> a proper @f$\mathcal{H}^2@f$-matrix,
 *  since the cluster bases of its immediate submatrices are allowed
 *  to differ from each other.
 *
 *  @param G Block matrix, will be overwritten by a proper
 *    @f$\mathcal{H}^2@f$-matrix approximation.
 *  @param pardepth Parallization depth. Parallel threads are spawned
 *    only on the next <tt>pardepth</tt> levels of the recursion.
 *  @param rw1 Total row weights for all submatrices, enumerated
 *    in column-major ordering, i.e., <tt>rw1[i+j*G->rsons]</tt>
 *    corresponds to the block in row <tt>i</tt> and column <tt>j</tt>.
 *  @param cw1 Total column weights for all submatrices, enumerated
 *    in column-major ordering, i.e., <tt>rw1[i+j*G->rsons]</tt>
 *    corresponds to the block in row <tt>i</tt> and column <tt>j</tt>.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy.
 *  @param rw If not a null pointer, will be set to point to a new
 *    @ref clusteroperator tree containing the total row weights
 *    of the new @f$\mathcal{H}^2@f$-matrix.
 *  @param cw If not a null pointer, will be set to point to a new
 *    @ref clusteroperator tree containing the total column weights
 *    of the new @f$\mathcal{H}^2@f$-matrix. */
HEADER_PREFIX void
unify_parallel_h2matrix(ph2matrix G, uint pardepth, pclusteroperator *rw1,
    pclusteroperator *cw1, pctruncmode tm, real eps, pclusteroperator *rw,
    pclusteroperator *cw);

/** @brief Converts an @ref rkmatrix into a @ref uniform matrix.
 *
 *  Returns an exact representation of an @ref rkmatrix by a new
 *  uniform matrix with new orthogonal row and column cluster bases
 *  and also computes corresponding total weight matrices.
 *  This function is intended for use with @ref unify_h2matrix in
 *  order to convert @ref rkmatrix approximations of leaves into
 *  @ref h2matrix objects.
 *
 *  @param r Original matrix.
 *  @param u Uniform matrix, will be overwritten.
 *  @param tm Truncation mode, used to set of weight matrices correctly.
 *  @param rw If not a null pointer, will be set to point to a new
 *    @ref clusteroperator tree containing the total row weights
 *    of the new uniform matrix.
 *  @param cw If not a null pointer, will be set to point to a new
 *    @ref clusteroperator tree containing the total column weights
 *    of the new uniform matrix. */
HEADER_PREFIX void
convert_rkmatrix_uniform(pcrkmatrix r, puniform u, pctruncmode tm,
    pclusteroperator *rw, pclusteroperator *cw);

/* ------------------------------------------------------------
 Compute adaptive cluster bases for H-matrices
 ------------------------------------------------------------ */

/** @brief Construct a row basis for a hierarchical matrix.
 *
 *  @param G Original matrix @f$G@f$.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy.
 *  @returns New row cluster basis. */
HEADER_PREFIX pclusterbasis
buildrowbasis_hmatrix(pchmatrix G, pctruncmode tm, real eps);

/** @brief Construct a column basis for a hierarchical matrix.
 *
 *  @param G Original matrix @f$G@f$.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy.
 *  @returns New column cluster basis. */
HEADER_PREFIX pclusterbasis
buildcolbasis_hmatrix(pchmatrix G, pctruncmode tm, real eps);

/* ------------------------------------------------------------
 Approximate H-matrix in new cluster bases
 ------------------------------------------------------------ */

/** @brief Construct an @f$\mathcal{H}^2@f$-matrix approximation
 *  of a given hierarchical matrix in given cluster bases
 *  by blockwise projection.
 *
 *  @param G Original matrix.
 *  @param rb New row basis.
 *  @param cb New column basis.
 *  @returns Approximating @f$\mathcal{H}^2@f$-matrix. */
HEADER_PREFIX ph2matrix
build_projected_hmatrix_h2matrix(pchmatrix G, pclusterbasis rb,
    pclusterbasis cb);

/* ------------------------------------------------------------
 Compute adaptive cluster bases for dense matrices
 ------------------------------------------------------------ */

/** @brief Construct a row basis for an array matrix.
 *
 *  @param G Original matrix @f$G@f$.
 *  @param b Block tree.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy.
 *  @returns New row cluster basis. */
HEADER_PREFIX pclusterbasis
buildrowbasis_amatrix(pcamatrix G, pcblock b, pctruncmode tm, real eps);

/** @brief Construct a column basis for an array matrix.
 *
 *  @param G Original matrix @f$G@f$.
 *  @param b Block tree.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy.
 *  @returns New column cluster basis. */
HEADER_PREFIX pclusterbasis
buildcolbasis_amatrix(pcamatrix G, pcblock b, pctruncmode tm, real eps);

/* ------------------------------------------------------------
 Approximate dense matrix in new cluster bases
 ------------------------------------------------------------ */

/** @brief Construct an @f$\mathcal{H}^2@f$-matrix approximation
 *  of a given array matrix in given cluster bases by blockwise projection.
 *
 *  @param G Original matrix.
 *  @param b Block tree.
 *  @param rb New row basis.
 *  @param cb New column basis.
 *  @returns Approximating @f$\mathcal{H}^2@f$-matrix. */
HEADER_PREFIX ph2matrix
build_projected_amatrix_h2matrix(pcamatrix G, pcblock b, pclusterbasis rb,
    pclusterbasis cb);

/** @} */

#endif
