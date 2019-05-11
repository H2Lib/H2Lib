
/* ------------------------------------------------------------
   This is the file "clusteroperator.h" of the H2Lib package.
   All rights reserved, Steffen Boerm 2009
   ------------------------------------------------------------ */

/** @file clusteroperator.h
 *  @author Steffen B&ouml;rm
 */

#ifndef CLUSTEROPERATOR_H
#define CLUSTEROPERATOR_H

/** @defgroup clusteroperator clusteroperator
 *  @brief Representation of cluster operators used to describe
 *  weights and transformations of cluster bases.
 *
 *  The @ref clusteroperator class represents a cluster operator
 *  @f$(Z_t)_{t\in\mathcal{T}_{\mathcal{I}}}@f$ corresponding to
 *  a cluster basis @f$(V_t)_{t\in\mathcal{T}_{\mathcal{I}}}@f$
 *  such that @f$Z_t@f$ and @f$V_t@f$ have the same number of
 *  columns.
 *
 *  Weight matrices for a cluster basis can be obtained using
 *  @ref weight_clusterbasis_clusteroperator , total weight matrices for an
 *  @ref h2matrix are available via @ref totalweights_h2matrix.
 *
 *  @{ */

/** @brief Representation of a cluster operator */
typedef struct _clusteroperator clusteroperator;

/** @brief Pointer to @ref clusteroperator object. */
typedef clusteroperator *pclusteroperator;

/** @brief Pointer to constant @ref clusterbasis object. */
typedef const clusteroperator *pcclusteroperator;

#include "cluster.h"
#include "clusterbasis.h"
#include "amatrix.h"

/** @brief Representation of a cluster operator */
struct _clusteroperator
{
  /** @brief Corresponding cluster. */
  pccluster t;
        
  /** @brief Number of rows. */
  uint krow;
  /** @brief Number of columns, usually equal to rank of corresponding cluster basis */
  uint kcol;
       
  /** @brief Coefficient matrix. */
  amatrix C;
        
  /** @brief Number of sons, either <tt>t->sons</tt> or zero. */
  uint sons;
  /** @brief Pointers to sons. */
  pclusteroperator *son;

  /** @brief References to this cluster operator. */
  uint refs;
};

#define CLUSTEROPERATOR_TYPE_COMPLETE

/* ------------------------------------------------------------
   Constructors and destructors
   ------------------------------------------------------------ */

/** @brief Initialize a @ref clusteroperator object.
 *
 *  Sets up the components of the object.
 *  If <tt>t</tt> is not a leaf, the array <tt>son</tt> is allocated,
 *  otherwise it is set to null.
 *
 *  @remark Should always be matched by a call to @ref uninit_clusteroperator.
 *
 *  @param co Object to be initialized.
 *  @param t Corresponding cluster.
 *  @returns Initialized @ref clusteroperator object. */
HEADER_PREFIX pclusteroperator
init_clusteroperator(pclusteroperator co, pccluster t);

/** @brief Initialize a @ref clusteroperator object for a leaf.
 *
 *  Sets up the components of the object.
 *  Sets <tt>son</tt> to null.
 *  If <tt>t->sons>0</tt>, this yields a partial cluster operator.
 *
 *  @remark Should always be matched by a call to @ref uninit_clusteroperator.
 *
 *  @param co Object to be initialized.
 *  @param t Corresponding cluster.
 *  @returns Initialized @ref clusteroperator object. */
HEADER_PREFIX pclusteroperator
init_leaf_clusteroperator(pclusteroperator co, pccluster t);

/** @brief Uninitializes a @ref clusteroperator object.
 *
 *  Invalidates pointers, freeing corresponding storage if appropriate,
 *  and prepares the object for deletion.
 *
 *  If this @ref clusteroperator references sons, these sons
 *  are unreferenced.
 *
 *  Only objects with <tt>cb->refs==0</tt> may be uninitialized.
 *
 *  @param co Object to be uninitialized. */
HEADER_PREFIX void
uninit_clusteroperator(pclusteroperator co);

/** @brief Create a new @ref clusteroperator object.
 *
 *  Allocates storage for the object and sets up its components.
 *
 *  @remark Should always be matched by a call to @ref del_clusteroperator.
 *  This call happens automatically if the last object referencing
 *  this object (see @ref ref_clusteroperator) is deleted.
 *
 *  @param t corresponding cluster.
 *  @return Returns the newly created @ref clusteroperator object.
 */
HEADER_PREFIX pclusteroperator
new_clusteroperator(pccluster t);

/** @brief Creates a new @ref clusteroperator object for a leaf.
 *
 *  Allocates storage for the object and sets up its components.
 *
 *  @remark Should always be matched by a call to @ref del_clusteroperator.
 *  This call happens automatically if the last object referencing
 *  this object (see @ref ref_clusteroperator) is deleted.
 *
 *  @param t Corresponding cluster.
 *  @return Returns the newly created @ref clusteroperator object.
 */
HEADER_PREFIX pclusteroperator
new_leaf_clusteroperator(pccluster t);

/** @brief Delete a @ref clusteroperator object.
 *
 *  Releases the storage corresponding to the object.
 *  If this @ref clusteroperator references sons, these sons
 *  are unreferenced.
 *
 *  Only objects with <tt>cb->refs==0</tt> may be deleted.
 *
 *  @param co Object to be deleted. */
HEADER_PREFIX void
del_clusteroperator(pclusteroperator co);

/** @brief Turn this @ref clusteroperator into a leaf.
 *
 *  Unreferences all sons, sets <tt>co->sons</tt> to zero
 *  and <tt>co->son</tt> to null.
 *
 *  @param co Target cluster operator. */
HEADER_PREFIX void
removesons_clusteroperator(pclusteroperator co);

/* ------------------------------------------------------------
   Reference counting
   ------------------------------------------------------------ */

/** @brief Set a pointer to a @ref clusterbasis object, increase its
 *  reference counter, and decrease reference counter of original
 *  pointer target.
 *
 *  @param ptr Pointer to the @ref pclusteroperator variable that will
 *         be changed.
 *  @param co @ref clusteroperator that will be referenced. */
HEADER_PREFIX void
ref_clusteroperator(pclusteroperator *ptr, pclusteroperator co);

/** @brief Reduce the reference counter of a @ref clusteroperator object.
 *
 *  If the reference counter reaches zero, the object is deleted.
 *
 *  @remark Use @ref ref_clusteroperator with <tt>cb=NULL</tt> instead,
 *  since this guarantees that the pointer is properly invalidated.
 *
 *  @param co @ref clusteroperator that will be unreferenced. */
HEADER_PREFIX void
unref_clusteroperator(pclusteroperator co);

/* ------------------------------------------------------------
   Low-level management
   ------------------------------------------------------------ */

/** @brief Updates bookkeeping information.
 *
 *  This function should be called after the sons of a
 *  @ref clusteroperator have been changed.
 *
 *  @remark This function will not update the sons recursively.
 *
 *  @param co @ref clusteroperator that will be updated. */
HEADER_PREFIX void
update_clusteroperator(pclusteroperator co);

/** @brief Change the number of rows and columns of a
 *  cluster operator and resize <tt>cb->C</tt> accordingly.
 *
 *  @param co Cluster operator that will be changed.
 *  @param krow New number of rows.
 *  @param kcol New number of columns. */
HEADER_PREFIX void
resize_clusteroperator(pclusteroperator co, uint krow, uint kcol);

/** @brief Find son @ref clusteroperator matching a given @ref cluster.
 *
 *  Find a son of <tt>cwf</tt> matching the cluster <tt>t</tt>.
 *  Returns <tt>cwf</tt> if no son was found.
 *
 *  @param cwf Father @ref clusteroperator.
 *  @param t Son @ref cluster.
 *  @returns Son of <tt>cwf</tt> with <tt>cwf->son[i]->t==t</tt>
 *         if it exists and <tt>cwf</tt> otherwise. */
HEADER_PREFIX pclusteroperator
identify_son_clusterweight_clusteroperator(pcclusteroperator cwf, pccluster t);

/* ------------------------------------------------------------
   Build cluster operator based on cluster tree of cluster basis
   ------------------------------------------------------------ */

/** @brief Construct a @ref clusteroperator matching a @ref cluster tree.
 *
 *  Constructs a @ref clusteroperator for a @ref cluster and its
 *  descendants.
 *
 *  Row and column numbers will be set to zero.
 *
 *  @param t Root @ref cluster.
 *  @returns @ref clusteroperator for <tt>t</tt> and its descendants. */
HEADER_PREFIX pclusteroperator
build_from_cluster_clusteroperator(pccluster t);

/** @brief Construct a @ref clusteroperator matching a
 *  @ref clusterbasis.
 *
 *  Constructs a @ref clusteroperator for a @ref clusterbasis and
 *  its descendants.
 *
 *  Row numbers will be set to zero, column numbers will be set
 *  to the ranks of the corresponding cluster basis.
 *
 *  @param cb Root @ref clusterbasis.
 *  @returns @ref clusteroperator for <tt>cb</tt> and its descendants. */
HEADER_PREFIX pclusteroperator
build_from_clusterbasis_clusteroperator(pcclusterbasis cb);

/* ------------------------------------------------------------
   Statistics
   ------------------------------------------------------------ */

/** @brief Get number of active @ref clusteroperator objects.
 *
 *  @returns Number of active @ref clusteroperator objects. */
HEADER_PREFIX uint
getactives_clusteroperator();

/** @brief Get the size of a given @ref clusteroperator object.
 *
 *  Yields the total size of this object and all its descendants.
 *
 *  @param co @ref clusteroperator root.
 *  @returns Size of allocated storage in bytes. */
HEADER_PREFIX size_t
getsize_clusteroperator(pcclusteroperator co);

/* ------------------------------------------------------------
   Debugging
   ------------------------------------------------------------ */

/** @brief Print the tree structure of a @ref clusteroperator.
 *
 *  @param co Root @ref clusteroperator. */
HEADER_PREFIX void
print_tree_clusteroperator(pcclusteroperator co);

/** @brief Compute and print the norms and norm differences
 *  of the coefficient matrices of two cluster operators.
 *
 *  @todo Usually a function called <tt>norm2diff</tt> should
 *  return a value, e.g., the maximal norm, instead of
 *  printing to <tt>stdout</tt>.
 *
 *  @param co1 First @ref clusteroperator.
 *  @param co2 Second @ref clusteroperator. */
HEADER_PREFIX void
norm2diff_clusteroperator(pcclusteroperator co1, pcclusteroperator co2);

/** @brief Compare two cluster weights.
 *
 *  Compare the Gramians of two cluster operators
 *  @f$(Z_{t,1})_t@f$ and @f$(Z_{2,t})_t@f$, i.e., compute
 *  @f$\frac{\|Z_{t,1}^* Z_{t,1} - Z_{t,2}^* Z_{t,2}\|_2}
 *          {\|Z_{t,1}^* Z_{t,1}\|_2\|_2}@f$ for all clusters
 *  and return the maximum.
 *  This approach makes the result invariant under orthogonal
 *  transformations applied to the rows.
 *
 *  @param co1 First @ref clusteroperator @f$Z_{1,t}@f$.
 *  @param co2 Second @ref clusteroperator @f$Z_{2,t}@f$.
 *  @returns Maximum of relative errors for all clusters in subtree. */
HEADER_PREFIX real
compareweights_clusteroperator(pcclusteroperator co1, pcclusteroperator co2);

/* ------------------------------------------------------------
   Enumeration
   ------------------------------------------------------------ */

/** @brief Enumerate @ref clusteroperator according to cluster tree.
 *
 *  The @ref clusteroperator elements are enumerated in an array
 *  of size <tt>t->desc</tt>. The enumeration starts with <tt>0</tt>
 *  assigned to the root and then proceeds with <tt>co->son[0]</tt>
 *  corresponding to the entries <tt>1</tt> to <tt>cb->son[0]->desc</tt>
 *  in the array and ending with <tt>co->son[sons-1]</tt> corresponding
 *  to the last <tt>t->son[sons-1]->desc</tt> entries.
 *
 *  @param t Cluster tree.
 *  @param co Cluster operator matching the cluster tree given by <tt>t</tt>.
 *  @returns Array of size <tt>t->desc</tt> containing pointers to
 *         the @ref clusteroperator objects corresponding to <tt>co</tt>
 *         and its descendants. */
HEADER_PREFIX pclusteroperator *
enumerate_clusteroperator(pccluster t, pclusteroperator co);

/* ------------------------------------------------------------
   Cluster basis product
   ------------------------------------------------------------ */

/** @brief Compute the cluster basis product @f$Z_t = V_{1,t}^* V_{2,t}@f$.
 *
 *  <tt>pr</tt> can only have sons if both <tt>cb1</tt> and <tt>cb2</tt>
 *  do. If <tt>pr</tt> is a leaf, currently either <tt>cb1</tt> or
 *  <tt>cb2</tt> has to be a leaf.
 *
 *  @param cb1 First cluster basis.
 *  @param cb2 Second cluster basis.
 *  @param pr Target cluster operator. */
HEADER_PREFIX void
basisproduct_clusteroperator(pcclusterbasis cb1, pcclusterbasis cb2,
			     pclusteroperator pr);

#endif

/* ------------------------------------------------------------
 * Access methods
 * ------------------------------------------------------------ */

/* Avoid problems with incomplete type definitions */
#if defined(CLUSTEROPERATOR_TYPE_COMPLETE) && !defined(CLUSTEROPERATOR_COMPLETE)
#define CLUSTEROPERATOR_COMPLETE

#ifdef __GNUC__
INLINE_PREFIX pamatrix
getC_clusteroperator(pclusteroperator) __attribute__ ((unused));
#endif

/* @brief Get the coefficient matrix @f$C@f$ of a @ref clusteroperator.
 *
 * @param co Cluster operator.
 * @returns Coefficient matrix @f$C@f$ */
INLINE_PREFIX pamatrix
getC_clusteroperator(pclusteroperator co)
{
  return &co->C;
}

#endif

/** @} */
