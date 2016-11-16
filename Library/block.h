
/* ------------------------------------------------------------
 * This is the file "block.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

/** @file block.h
 *  @author Steffen B&ouml;rm
 */

#ifndef BLOCK_H
#define BLOCK_H

/** @defgroup block block
 *  @brief Representation of a block tree.
 * 
 * The @ref block class is used to represent block trees.
 * @{*/

/** @brief Representation of a block tree.*/
typedef struct _block block;

/** @brief Pointer to @ref block object.*/
typedef block *pblock;

/** @brief Pointer to constant @ref block object.*/
typedef const block *pcblock;

#ifdef USE_CAIRO
#include <cairo/cairo.h>
#endif

#include "cluster.h"
#include "settings.h"

/** @brief Representation of block trees.
 * 
 * A block tree @f$\mathcal{T}_{I \times J}@f$ is a labeled tree
 * for two @ref cluster trees @f$\mathcal{T}_{I}@f$ and
 * @f$\mathcal{T}_{J}@f$ and stored recursively.
 * The block tree is described by a row @ref cluster tree <tt> rc</tt> 
 * and a column @ref cluster tree <tt> cc </tt>. The admissibility flag 
 * <tt>a</tt> can be used to distinguish admissible (farfield) blocks from 
 * inadmissible (nearfield) blocks. If the block is subdivided, pointers to the 
 * son blocks can be stored, additionally the number of the row sons, column sons
 * and the number of decendents.
 * */
struct _block {
  /** @brief Row cluster.*/
  pcluster rc;

  /** @brief Column Cluster.*/
  pcluster cc;

  /** @brief Admissibility flag.*/
  bool a;

  /** @brief Son blocks, if the block is subdivided.*/
  pblock* son;

  /** @brief Number of row sons.*/
  uint rsons;

  /** @brief Number of column sons.*/
  uint csons;

  /** @brief Number of descendants.*/
  uint desc;
};

/* ------------------------------------------------------------
 * Admissibility condition
 * ------------------------------------------------------------ */

/** @brief Callback function for an admissibility condition*/
typedef bool (*admissible)(pcluster rc, pcluster cc, void *data);

/** @brief Check the euclidian (maximum) admissibility condition.
 * 
 * Checks the euclidian (maximum) admissibility of two @ref cluster trees 
 * @f$ s @f$ and @f$ t @f$.
 * The block @f$ (s,t) @f$ is admissible, if
 * @f$ \max {\text{diam}_2 (B_s)  \text{diam}_2 (B_s)} \leq \eta 
 * \text{dist}_2 (B_s, B_t) @f$. 
 * @f$ B_s @f$ and @f$ B_t @f$ are the bounding boxes for @f$ s @f$ and 
 * @f$ t @f$, @f$ \text{diam}_2 @f$ is the euclidian diameter, 
 * @f$ \text{dist}_2 @f$ the euclidian distance.
 * 
 * @param rc Row cluster @f$t@f$.
 * @param cc Col cluster @f$s@f$.
 * @param data Has to be a pointer to real eta.
 * @return TRUE, if the block @f$ (t,s) @f$ is admissible, otherwise FALSE.
 */
HEADER_PREFIX bool
admissible_2_cluster(pcluster rc, pcluster cc, void* data);

/** @brief Check the admissibility condition in the maximum norm.
 * 
 * Checks the admissibility of two @ref cluster trees 
 * @f$ s @f$ and @f$ t @f$ in the maximum norm.
 * The block @f$ (s,t) @f$ is admissible, if
 * @f$ \max {\text{diam}_{\infty} (B_s)  \text{diam}_{\infty} (B_s)} \leq \eta 
 * \text{dist}_{\infty} (B_s, B_t) @f$. 
 * @f$ B_s @f$ and @f$ B_t @f$ are the bounding boxes for @f$ s @f$ and 
 * @f$ t @f$, @f$ \text{diam}_{\infty} @f$ is the euclidian diameter, 
 * @f$ \text{dist}_{\infty} @f$ the euclidian distance.
 * 
 * @param rc Row cluster @f$t@f$.
 * @param cc Col cluster @f$s@f$.
 * @param data Has to be a pointer to real eta.
 * @returns TRUE, if the block @f$ (t,s) @f$ is admissible, otherwise FALSE.
 */
HEADER_PREFIX bool
admissible_max_cluster(pcluster rc, pcluster cc, void* data);

/** @brief Check the admissibility condition for spherical domains.
 *  
 * Checks the admissibility condition for spherical domains.
 * The block @f$ (s,t) @f$ is admissible, if
 * @f$ \max {\text{diam}_{2} (B_s)  \text{diam}_{2} (B_s)} \leq \eta 
 * \left\|c_s - c_t \right\|_2 - r_t - r_s @f$. 
 * 
 * @f$ B_s @f$ and @f$ B_t @f$ are the bounding boxes for @f$ s @f$ and 
 * @f$ t @f$, @f$ \text{diam}_{\infty} @f$ is the euclidian diameter, 
 * @f$ \text{dist}_{\infty} @f$ the euclidian distance.
 * 
 * @param rc Row cluster @f$t@f$.
 * @param cc Col cluster @f$s@f$.
 * @param data Has to be a pointer to real eta.
 * @returns TRUE, if the block @f$ (t,s) @f$ is admissible, otherwise FALSE.
 */
HEADER_PREFIX bool
admissible_sphere_cluster(pcluster rc, pcluster cc, void* data);

/** @brief Check the euclidian minimum admissibility condition.
 * 
 * Checks the euclidianminimum admissibility of two @ref cluster trees 
 * @f$ s @f$ and @f$ t @f$.
 * The block @f$ (s,t) @f$ is admissible, if
 * @f$ \min {\text{diam}_2 (B_s)  \text{diam}_2 (B_s)} \leq \eta 
 * \text{dist}_2 (B_s, B_t) @f$. 
 * @f$ B_s @f$ and @f$ B_t @f$ are the bounding boxes for @f$ s @f$ and 
 * @f$ t @f$, @f$ \text{diam}_2 @f$ is the euclidian diameter, 
 * @f$ \text{dist}_2 @f$ the euclidian distance.
 * 
 * @param rc Row cluster @f$t@f$.
 * @param cc Col cluster @f$s@f$.
 * @param data Has to be a pointer to real eta.
 * @return TRUE, if the block @f$ (t,s) @f$ is admissible, otherwise FALSE.
 */
HEADER_PREFIX bool
admissible_2_min_cluster(pcluster rc, pcluster cc, void* data);

/* ------------------------------------------------------------
 Constructors and destructors
 ------------------------------------------------------------ */

/** @brief Create a new @ref block object.
 * 
 * Allocates storage for the object and sets the pointers to the sons to NULL.
 * 
 * @remark Should always be matched by a call to @ref del_cluster.
 * 
 * @param rc Row cluster.
 * @param cc Column cluster.
 * @param a Admissibility flag.
 * @param rsons Number of sons of the row cluster.
 * @param csons Number of sons of the column cluster.
 * @returns Returns a new @ref block tree.
 */
HEADER_PREFIX pblock
new_block(pcluster rc, pcluster cc, bool a, uint rsons, uint csons);

/** @brief Delete a @ref block tree.
 * 
 * Releases the storage corresponding to the object, including the storage of
 * all its descendants.
 *
 * @param b Object to be deleted. 
 */
HEADER_PREFIX void
del_block(pblock b);

/** @brief Complete initialization of a @ref block object.
 * 
 * Completes initialization of a @ref block tree object after all sons
 * have been initialized. This function computes the number of descendants of
 * the block tree @f$ b @f$.
 * 
 * @remark Should be called after all sons of the block tree have been 
 * initialized. 
 * 
 * @param b Object to be updated.
 */
HEADER_PREFIX void
update_block(pblock b);

/* ------------------------------------------------------------
 Block clustering strategies
 ------------------------------------------------------------ */

/** @brief Build a non strict @ref block tree.
 * 
 * Builds a non strict block tree from the row @ref cluster tree 
 * @f$ rc @f$ and the column @ref cluster tree @f$ cc @f$. A block tree
 * is non-strict if all leaves @f$ (t,s) @f$ are admissible or @f$ t @f$ or
 * @f$ s @f$ have no sons.
 * 
 * @param rc Row cluster.
 * @param cc Col cluster.
 * @param data Necessary data for the admissibility condition.
 * @param admis Admissibility condition.
 * @returns Returns a non strict block tree.
 */
HEADER_PREFIX pblock
build_nonstrict_block(pcluster rc, pcluster cc, void *data, admissible admis);

/** @brief Build a non strict lower triangular @ref block tree.
 *
 * Builds a non strict lower triangular block tree from the row @ref cluster tree
 * @f$ rc @f$ and the column @ref cluster tree @f$ cc @f$. A block tree
 * is non-strict if all leaves @f$ (t,s) @f$ are admissible or @f$ t @f$ or
 * @f$ s @f$ have no sons.
 *
 * @param rc Row cluster.
 * @param cc Col cluster.
 * @param data Necessary data for the admissibility condition.
 * @param admis Admissibility condition.
 * @returns Returns a non strict lower triangular block tree.
 */
HEADER_PREFIX pblock
build_nonstrict_lower_block(pcluster rc, pcluster cc, void *data, admissible admis);

/** @brief Build a strict @ref block tree.
 * 
 * Builds a strict block tree from the row @ref cluster tree @f$ rc @f$
 * and the column @ref cluster tree @f$ cc @f$. A block tree is called 
 * strict, if all leaves @f$ (t,s) @f$ are admissible or @f$ t @f$ and @f$ s @f$
 * are leave cluster. 
 * 
 * @param rc Row cluster.
 * @param cc Col Cluster.
 * @param data Necessary data for the admissibility condition.
 * @param admis Admissibility condition.
 * @returns Returns a strict block tree.
 */
HEADER_PREFIX pblock
build_strict_block(pcluster rc, pcluster cc, void *data, admissible admis);

/** @brief Build a strict lower triangular @ref block tree.
 *
 * Builds a strict lower triangular block tree from the row
 * @ref cluster tree @f$ rc @f$ and the column @ref cluster tree @f$ cc @f$.
 * A block tree is called  strict, if all leaves @f$ (t,s) @f$ are
 * admissible or @f$ t @f$ and @f$ s @f$ are leave cluster.
 *
 * @param rc Row cluster.
 * @param cc Col Cluster.
 * @param data Necessary data for the admissibility condition.
 * @param admis Admissibility condition.
 * @returns Returns a strict lower triangular block tree.
 */
HEADER_PREFIX pblock
build_strict_lower_block(pcluster rc, pcluster cc, void *data, admissible admis);

/* ------------------------------------------------------------
 * Drawing block trees
 * ------------------------------------------------------------ */

/** @brief Draw a block tree.
 * 
 */
#ifdef USE_CAIRO
/**
 * @brief Draw a block tree to a cairo surface.
 *
 * @param cr Cairo surface to be drawn to.
 * @param b The block tree that has to be drawn.
 * @param levels Number of levels of the block tree, that should be drawn.
 *   If @p levels == 0 holds, all levels will be drawn.
 */
HEADER_PREFIX void
draw_cairo_block(cairo_t *cr, pcblock b, int levels);
#endif

/* ------------------------------------------------------------
 * Interactive visualization
 * ------------------------------------------------------------ */

/**
 * @brief Visualize a block tree in OpenGL.
 *
 * @param b The block tree that has to be drawn.
 */
HEADER_PREFIX void
view_block(pcblock b);

/* ------------------------------------------------------------
 * Hierarchical iterators
 * ------------------------------------------------------------ */

/** @brief Representation of a @ref blockentry object.*/
typedef struct _blockentry blockentry;

/**  @brief Pointer to a @ref blockentry object.*/
typedef blockentry *pblockentry;

/** @brief Pointer to a constant @ref blockentry object.*/
typedef const blockentry *pcblockentry;

/** @brief Auxiliary structure for hierarchical iterators for a @ref block
 *         cluster tree. */
struct _blockentry {
  /** @brief block tree.*/
  pcblock b;
  /** @brief Number of the block tree.*/
  uint bname;
  /** @brief Number of the row cluster of <tt>b</tt>.*/
  uint rname;
  /** @brief Number of the column cluster of <tt>b</tt>.*/
  uint cname;
  /** @brief Pointer to the blockentry object of the father of <tt>b</tt>.*/
  pcblockentry father;
  /** @brief Pointer to the blockentry object of the sucessor of <tt>b</tt>.*/
  pblockentry next;
};

/** @brief Hierarchical iterator for a @ref block tree.
 * 
 * Iterate over the blocktree and all its descendants.
 * The <tt>pre</tt> function is called for each element before its descendants
 * are processed, the <tt>post</tt> function is called afterwards
 * 
 * <tt>tname</tt> will take values between <tt>tname</tt> and 
 * <tt>tname+t->sons-1</tt>, where <tt>tname</tt> is interpreted as the index of
 * the root, while the indices of the first son start at <tt>tname+1,</tt> the
 * indices of the second start at <tt>tname+1+1t->son[0]->desc,</tt> and so on.
 * 
 * @param b Block tree.
 * @param bname Number of the block tree .
 * @param rname Number of the row cluster tree in <tt>b</tt>.
 * @param cname Number of the column cluster tree in <tt>b</tt>.  .
 * @param pre Function to be called before the descendants of <tt>b</tt> are
 *            processed.
 * @param post Function to be called after the descendants of <tt>b</tt> are 
 * 		processed.
 * @param data Auxiliary data for Callback Functions.
 */
HEADER_PREFIX void
iterate_block(pcblock b, uint bname, uint rname, uint cname,
    void (*pre)(pcblock b, uint bname, uint rname, uint cname, uint pardepth,
        void *data),
    void (*post)(pcblock b, uint bname, uint rname, uint cname, uint pardepth,
        void *data), void *data);

/** @brief Iterate through all subblocks of a @ref block tree, rowwise.
 * 
 * Iterate through all subblocks of a @ref block tree and collect all 
 * row blocks in an @ref h2matrixlist object.
 * 
 * @param b Block tree.
 * @param bname Number of the block tree.
 * @param rname Number of the row cluster tree in <tt>b</tt>.
 * @param cname Number of the column cluster tree in <tt>b</tt>.
 * @param pardepth Parallelization depth.
 * @param pre Function to be called before the sons of <tt>b</tt> are processed.
 * @param post Function to be called after the sons of <tt>b</tt> are processed.
 * @param data Auxiliary data for Callback Functions.
 */
HEADER_PREFIX void
iterate_rowlist_block(pcblock b, uint bname, uint rname, uint cname,
    uint pardepth, void (*pre)(pcblockentry pb, uint pardepth, void *data),
    void (*post)(pcblockentry pb, uint pardepth, void *data), void *data);

/** @brief Iterate through all subblocks of a @ref block tree, 
 * columnwise.
 * 
 * Iterate through all subblocks of a @ref block tree and collect all 
 * column blocks in an @ref blockentry object.
 * 
 * @param b Block tree.
 * @param bname Number of the block tree.
 * @param rname Number of the row cluster tree  in <tt>b</tt>.
 * @param cname Number of the column cluster tree in <tt>b</tt>.
 * @param pre Function to be called before the sons of <tt>b</tt> are processed.
 * @param pardepth Parallelization depth.
 * @param post Function to be called after the sons of <tt>b</tt> are processed.
 * @param data Auxiliary data for Callback Functions.
 */
HEADER_PREFIX void
iterate_collist_block(pcblock b, uint bname, uint rname, uint cname,
    uint pardepth, void (*pre)(pcblockentry pb, uint pardepth, void *data),
    void (*post)(pcblockentry pb, uint pardepth, void *data), void *data);

/** @brief Iterate through all subblocks of a @ref block tree.
 * 
 * If the iterator works with multiple threads, it guarantees that
 * threads running in parallel call the <tt>pre</tt> and <tt>post</tt>
 * functions with different row clusters.
 * 
 * @param b Block tree.
 * @param bname Number of the block tree.
 * @param rname Number of the row cluster tree in <tt>b</tt> .
 * @param cname Number of the column cluster tree in <tt>b</tt>.
 * @param pardepth Parallelization depth.
 * @param pre Function to be called before the descendants of <tt>b</tt> are
 *            processed.
 * @param post funczion to be called after the descendants of <tt>b</tt> are
 *             processed.
 * @param data Auxiliary data for Callback Functions.
 */
HEADER_PREFIX void
iterate_byrow_block(pcblock b, uint bname, uint rname, uint cname,
    uint pardepth,
    void (*pre)(pcblock b, uint bname, uint rname, uint cname, uint pardepth,
        void *data),
    void (*post)(pcblock b, uint bname, uint rname, uint cname, uint pardepth,
        void *data), void *data);

/** @brief Iterate through all subblocks of a @ref block tree.
 * 
 * If the iterator works with multiple threads, it guarantees that
 * threads running in parallel call the <tt>pre</tt> and <tt>post</tt>
 * functions with different column clusters.
 * 
 * @param b Block tree.
 * @param bname Number of the block tree.
 * @param rname Number of the row cluster tree in <tt>b</tt>.
 * @param cname Number of the column cluster tree in <tt>b</tt>.
 * @param pardepth Parallelization depth.
 * @param pre Function to be called before the descendants of <tt>b</tt> are 
 *            processed.
 * @param post Function to be called before the descendants of <tt>b</tt> are 
 *            processed.
 * @param data Auxiliary data for Callback Functions.
 */
HEADER_PREFIX void
iterate_bycol_block(pcblock b, uint bname, uint rname, uint cname,
    uint pardepth,
    void (*pre)(pcblock b, uint bname, uint rname, uint cname, uint pardepth,
        void *data),
    void (*post)(pcblock b, uint bname, uint rname, uint cname, uint pardepth,
        void *data), void *data);

/* ------------------------------------------------------------
 * Enumeration
 * ------------------------------------------------------------ */

/** @brief Enumerate a block tree.
 * 
 * Enumerates the @ref block tree @f$ t @f$ in an array  of size 
 * @f$ t->desc @f$ by passing through all descendants and ordering them in a 
 * pointer to a @ref block tree object.
 * The enumeration starts with @f$ 0 @f$ assigned to the root and then proceeds
 * with depth-first search.
 * 
 * @param t block tree to be enumerated.
 * @returns Returns a pointer, length @f$ t->desc @f$, to a block tree
 * enumerating all descendants. 
 */
HEADER_PREFIX pblock *
enumerate_block(pblock t);

/** @brief Enumerate the levels of a block tree.
 * 
 * Enumerates the @ref block tree @f$ t @f$ in an array of size 
 * @f$ t->desc @f$ by passing through all descendants, ordering them and filling
 * an array with information of the level of each descendant in the block 
 * cluster tree.
 * The array includes at position @f$ i @f$ the level of the descendant 
 * @f$ i @f$ of the block tree.
 * 
 * @param t block tree to be enumerated.
 * @returns Returns a pointer to unsigned int, length @f$ t->desc @f$, with the 
 * level of each descendant.  
 */
HEADER_PREFIX uint*
enumerate_level_block(pblock t);

/* ------------------------------------------------------------
 * Utility functions
 * ------------------------------------------------------------ */

/** @brief Compute the depth of a @ref block tree object
 *
 * Compute the maximal depth of a block tree.
 * 
 * @param b block object
 * @return Returns the maximal depth of @p b.
 */
HEADER_PREFIX uint
getdepth_block(pcblock b);

/** @brief Compute the sparsity of a @ref block object
 * 
 * Computes the @f$ C_{sp}@f$-sparsity of the block tree @f$ b @f$.
 * The @ref block tree is called @f$ C_{sp}@f$-sparse, if the number
 * of its block rows and block columns is at most @f$ C_{sp}@f$.
 * 
 * @param b Block tree.
 * @returns Returns the @f$C_{sp} @f$-spars of the blcok cluster tree.
 */
HEADER_PREFIX uint
compute_csp_block(pcblock b);

#endif

/** @}*/
