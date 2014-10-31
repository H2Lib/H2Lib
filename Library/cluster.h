/* ------------------------------------------------------------
 This is the file "cluster.h" of the H2Lib package.
 All rights reserved, Knut Reimer 2009
 ------------------------------------------------------------ */

/** @file cluster.h
 *  @author Knut Reimer
 */

#ifndef CLUSTER_H
#define CLUSTER_H

/** @defgroup cluster cluster
 *  @brief Representation of a cluster tree.
 * 
 * The @ref cluster class is used to represent cluster trees including 
 * coordinates of the corresponding bounding box and pointers to its sons.
 * @{*/

/** @brief Representation of a cluster tree. */
typedef struct _cluster cluster;

/** @brief Pointer to @ref cluster object. */
typedef cluster* pcluster;

/** @brief Pointer to constant @ref cluster object.*/
typedef const cluster *pccluster;

#include "settings.h"

/** @brief Representation of cluster trees.
 * 
 * Cluster trees are recursively represented labeled trees. 
 * The labels are subsets of the index set <tt>idx</tt> and the labels of the 
 * sons form a disjunct partition of the label of the father.
 * The structure provides storage for the spatial dimension and the coordinates 
 * of the axis-parallel bounding boxes. */
struct _cluster {
  /** @brief Number of indices. */
  uint size;

  /** @brief Index set. */
  uint *idx;

  /** @brief Number of sons.*/
  uint sons;

  /** @brief Pointer to son clusters.*/
  pcluster* son;

  /** @brief Spatial dimension of bounding box.*/
  uint dim;

  /** @brief Minimal coordinates of bounding box. */
  real *bmin;

  /** @brief Maximal coordinates of bounding box.*/
  real *bmax;

  /** @brief Number of descendants.*/
  uint desc;

  /** @brief Type of cluster, necessary for domain decomposition clustering.
   * 1 : domain cluster
   * 2 : interface cluster
   */
  uint type;
};

/* ------------------------------------------------------------
 Constructors and destructors
 ------------------------------------------------------------ */

/** @brief Create a new @ref cluster object.
 * 
 * Allocates storage for the object and sets the pointers to the sons to NULL.
 * 
 * @remark Should always be matched by a call to @ref del_cluster.
 * 
 * @param size Number of indices.
 * @param idx Index set.
 * @param sons Number of sons.
 * @param dim Spatial dimension of bounding box. 
 * @returns New @ref cluster object.*/
HEADER_PREFIX pcluster
new_cluster(uint size, uint *idx, uint sons, uint dim);

/** @brief Delete a @ref cluster object.
 * 
 * Releases the storage corresponding to the @ref cluster object.
 * If the cluster has sons, their storage is released too.
 * 
 * @param t Cluster object to be deleted.*/
HEADER_PREFIX void
del_cluster(pcluster t);

/** @brief Complete the initialisation of a @ref cluster object.
 *
 * Complete the initialisation of the @ref cluster object after all sons have
 * been initialised. 
 * The Number of the descendants of the cluster is counted.
 * 
 * @param t Cluster object to be completed.
 */
HEADER_PREFIX void
update_cluster(pcluster t);

/* ------------------------------------------------------------
 Statistics
 ------------------------------------------------------------ */

/** @brief Compute the depth of a @ref cluster object.
 *
 * Compute the maximal depth of a cluster tree.
 * 
 * @param t Cluster object.
 * @returns Depth of the @ref cluster tree.
 */
HEADER_PREFIX uint
getdepth_cluster(pccluster t);

/** @brief Compute the minimal level of a @ref cluster object.
 * 
 * Compute the minimal level of the leaf clusters in a cluster tree.
 * 
 * @param t Cluster object.
 * @returns Minimal level of the leaf clusters of the cluster tree.
 */
HEADER_PREFIX uint
getmindepth_cluster(pccluster t);

/* ------------------------------------------------------------
 Structure Adaptation
 ------------------------------------------------------------ */

/* Extends t until getmindepth_cluster(t) == depth */
HEADER_PREFIX void
extend_cluster(pcluster t, uint depth);

/** @brief Cut a @ref cluster object until a new depth is reached.
 * 
 * Cut a cluster tree until a given depth is reached.
 * Parts of the cluster tree, which exceed the new depth, are deleted.
 * 
 * @param t Cluster object to be cut.
 * @param depth Depth of the cluster tree after the cut.
 */
HEADER_PREFIX void
cut_cluster(pcluster t, uint depth);

/** @brief Balance a @ref cluster tree to a given depth.
 * 
 * Balances a cluster tree @f$ t @f$ in such a way that
 * @ref getdepth_cluster @f$(t) == @f$ @ref getmindepth_cluster @f$(t) 
 * == depth @f$.
 * 
 * @param t Cluster object to be balanced.
 * @param depth New depth of the cluster tree.
 */
HEADER_PREFIX void
balance_cluster(pcluster t, uint depth);

/** @brief Coarsen a @ref cluster tree to a minimal size.
 * 
 * Coarsens a cluster object @f$ t @f$ in such a way that
 * @f$ t->size \geq \text{minsize} @f$ for all clusters @f$ t @f$.
 * 
 * @param t Cluster object.
 * @param minsize minimal size for all cluster trees. 
 */
HEADER_PREFIX void
coarsen_cluster(pcluster t, uint minsize);

/** @brief Set the number of sons of a @ref cluster  tree.
 * 
 * Sets the number of sons of a given @ref cluster object to a given number.
 * All sons are initialised with NULL.
 *  
 * @param t cluster object.
 * @param sons Number of sons of the cluster object after calling this function.
 */
HEADER_PREFIX void
setsons_cluster(pcluster t, uint sons);

/* ------------------------------------------------------------
 Bounding box
 ------------------------------------------------------------ */

/** @brief Compute the euclidian diameter of the bounding box of a cluster.
 * 
 * Computes the euclidian diameter of the bounding box @f$ B_t @f$ of a 
 * @ref cluster object @f$ t @f$:
 * @f$ \text{diam}_2 (B_t) := \text{max} \left\{ \left\| x-y\right\|_2 : 
 * x,y \in B_t\right\}@f$.
 * 
 * @param t Of this cluster the euclidian diameter is computed.
 * @return Returns the euclidian diameter of the bounding box.
 */
HEADER_PREFIX real
getdiam_2_cluster(pccluster t);

/** @brief Compute the euclidian distance of two clusters.
 * 
 * Computes the euclidian distance of two bounding boxes @f$ B_t @f$ and 
 * @f$ B_s @f$ of two clusters @f$ t @f$ and @f$ s @f$:
 * @f$ \text{dist}_2(B_t, B_s) := \min \left\{ \left\| x-y \right\|_2 :
 * x \in B_t, y \in B_s \right\}.@f$
 * 
 * @param t Row cluster.
 * @param s Col cluster.
 * @return Returns the euclidian distance of two clusters.
 */
HEADER_PREFIX real
getdist_2_cluster(pccluster t, pccluster s);

/** @brief Compute the diameter of the bounding of a cluster in the maximum norm.  
 * 
 * Compute the diameter of the bounding box @f$ B_t @f$ of a @ref cluster object
 * @f$ t @f$ in the maximum norm:
 * @f$ \text{diam}_{\infty} (B_t) := \left\{  \left\| x-y \right\|_{\infty} :
 * x, y \in B_t \right\}@f$.
 * 
 * @param t Of this cluster the diameter is computed.
 * @return Returns the diameter of the bounding box in the maximum norm.
 */
HEADER_PREFIX real
getdiam_max_cluster(pccluster t);

/** @brief Compute the distance of two clusters in the maximum norm.
 * 
 * Computes the distance of the bounding boxes of two @ref cluster objects 
 * @f$ t @f$ and @f$ s @f$ in the maximum norm:
 * @f$ \text{dist}_{\infty} (B_t, B_s) := \min \left\{ \left\| x-y \right\|
 * _{\infty} : x \in B_t, y \in B_s \right\}@f$.
 * 
 * @param t Row cluster. 
 * @param s Col cluster.
 * @return Returns the distant of two clusters in the maximum norm.
 */
HEADER_PREFIX real
getdist_max_cluster(pccluster t, pccluster s);

/* ------------------------------------------------------------
 Hierarchical iterator
 ------------------------------------------------------------ */

/** @brief Hierarchical iterator for a @ref cluster tree.
 * 
 * Iterate over the cluster basis and all its descendants.
 * The <tt>pre</tt> function is called for each element before its descendants
 * are processed, the <tt>post</tt> function is called afterwards
 * The cluster tree is iterated in a depth first scheme:
 * 
 * <tt>tname</tt> will take values between <tt>tname</tt> and 
 * <tt>tname+t->sons-1</tt>, where <tt>tname</tt> is interpreted as the index of
 * the root, while the indices of the first son start at <tt>tname+1,</tt> the
 * indices of the second start at <tt>tname+1+1t->son[0]->desc,</tt> and so on.
 * 
 * @remark If the cluster tree is included in an enumerated object like 
 * an enumerated clusterbasis, the parameter tname can be used to adresss to 
 * this ordered object.
 * 
 * @param t Root of the cluster tree.
 * @param tname Identificator of the (partial) cluster tree, initially 0.
 * @param pre Callback function applied before iterating through childs.
 * @param post Callback function applied after iterating through childs.
 * @param data Auxiliary data for the callback functions <tt> pre </tt> and 
 * <tt> post </tt>.
 */
HEADER_PREFIX void
iterate_cluster(pccluster t, uint tname,
    void (*pre)(pccluster t, uint tname, void *data),
    void (*post)(pccluster t, uint tname, void *data), void *data);

/** @brief Parallel hierarchical iterator for a cluster tree.
 * 
 * Iterate over the cluster basis and all its descendants.
 * The <tt>pre</tt> function is called for each element before its descendants
 * are processed, the <tt>post</tt> function is called afterwards
 * The cluster tree is iterated in a depth first scheme:
 * 
 *  * <tt>tname</tt> will take values between <tt>tname</tt> and 
 * <tt>tname+t->sons-1</tt>, where <tt>tname</tt> is interpreted as the index of
 * the root, while the indices of the first son start at <tt>tname+1,</tt> the
 * indices of the second start at <tt>tname+1+1t->son[0]->desc,</tt> and so on.
 * 
 * @remark If the cluster tree is included in an enumerated object like 
 * an enumerated clusterbasis, the parameter tname can be used to adresss to 
 * this ordered object.
 * 
 * @param t Root of the cluster tree.
 * @param tname Identificator of the (partial) cluster tree, initially 0.
 * @param pardepth Parallelization depth.
 * @param pre  Callback function applied for iterating through childs.
 * @param post Callback function applied for iterating through childs.
 * @param data Auxiliary data for the callback funtions <tt> pre </tt> and 
 * <tt> post </tt>.
 */
HEADER_PREFIX void
iterate_parallel_cluster(pccluster t, uint tname, uint pardepth,
    void (*pre)(pccluster t, uint tname, void *data),
    void (*post)(pccluster t, uint tname, void *data), void *data);

/* ------------------------------------------------------------
 Enumeration
 ------------------------------------------------------------ */
/** @brief Enumerate a cluster tree.
 * 
 * Enumerates the @ref cluster tree @f$ t @f$ in an array of size @f$ t->desc @f$ 
 * by passing through all descendants
 * and ordering them in a pointer to a @ref cluster tree object.
 * The enumeration starts with @f$ 0 @f$ assigned to the root and then proceeds
 * with depth-first search.
 * 
 * @param t Cluster tree to be enumerated.
 * @return Returns a pointer, length @f$ t->desc @f$, to the cluster tree
 * enumerating all descendants. 
 */
HEADER_PREFIX pcluster *
enumerate_cluster(pcluster t);

/** @}*/

#endif
