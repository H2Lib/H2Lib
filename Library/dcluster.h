
/* ------------------------------------------------------------
 * This is the file "dcluster.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2015
 * ------------------------------------------------------------ */

/** @file dcluster.h
 *  @author Steffen B&ouml;rm */

#ifndef DCLUSTER_H
#define DCLUSTER_H

/** @defgroup dcluster dcluster
 *  @brief Directional cluster trees.
 *  @{ */

/** @brief Directional cluster tree. */
typedef struct _dcluster dcluster;

/** @brief Pointer to @ref dcluster object. */
typedef dcluster *pdcluster;

/** @brief Pointer to constant @ref dcluster object. */
typedef const dcluster *pcdcluster;

/** @brief Directions for all levels of a directional cluster tree. */
typedef struct _leveldir leveldir;

/** @brief Pointer to @ref leveldir object. */
typedef leveldir *pleveldir;

/** @brief Pointer to constant @ref leveldir object. */
typedef const leveldir *pcleveldir;

#include "settings.h"
#include "cluster.h"
#include <stdlib.h>

/** @brief Directional cluster tree. */
struct _dcluster
{
  /** @brief Number of indices */
  uint size;

  /** @brief Index array,
   *    type <tt>uint idx[size]</tt> */
  uint *idx;

  /** @brief Number of sons */
  uint sons;

  /** @brief Pointers to sons,
   *    type <tt>pdcluster son[sons]</tt> */
  pdcluster *son;

  /** @brief Dimension of bounding box */
  uint dim;
  
  /** @brief Minimal coordinates of bounding boxes,
   *    type <tt>real bmin[dim]</tt> */
  real *bmin;

  /** @brief Maximal coordinates of bounding boxes,
   *    type <tt>real bmax[dim]</tt> */
  real *bmax;

  /** @brief Number of directions */
  uint directions;

  /** @brief Directions,
   *    type <tt>real dir[directions][dim]</tt>
   *
   *    This is only a pointer to an array in a @ref leveldir object,
   *    so it should not be changed or freed directly. */
  pcreal *dir;

  /** @brief Corresponding directions in sons,
   *    type <tt>uint dirson[sons][directions]</tt> */
  uint **dirson;

  /** @brief Number of descendants in tree */
  uint desc;
};

/** @brief Families of directions for all levels of the cluster tree.
 *
 *  Directional interpolation requires that row and column clusters
 *  share the same directions. This is easily achieved by using the
 *  same directions for all clusters on a given level and ensuring
 *  the the row and column clusters of blocks are always on the same
 *  level, i.e., that the block tree is level-consistent. */
struct _leveldir
{
  /** @brief Depth of the cluster tree */
  uint depth;

  /** @brief Spatial dimension */
  uint dim;

  /** @brief Diameter of largest cluster on this level */
  preal maxdiam;

  /** @brief Splitting parameters */
  uint *splits;

  /** @brief Number of directions */
  uint *directions;

  /** @brief Directions,
   *    type <tt>real dir[depth][directions][dim]</tt> */
  preal **dir;

  /** @brief Auxiliary storage */
  preal dirmem;
};

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

/** @brief Create a new directional cluster.
 *
 *  @param size Number of indices.
 *  @param idx Array of indices.
 *  @param sons Number of sons.
 *  @param dim Spatial dimension (e.g., of the bounding box).
 *  @returns New directional cluster. */
HEADER_PREFIX pdcluster
new_dcluster(uint size, uint *idx, uint sons, uint dim);

/** @brief Update a directional cluster.
 *
 *  This function is called after the sons of a cluster
 *  have been created or changed, e.g., to compute
 *  the number of descendants <tt>desc</tt>.
 *
 *  @param t Cluster to be updated. */
HEADER_PREFIX void
update_dcluster(pdcluster t);

/** @brief Delete a directional cluster tree.
 *
 *  @param t Root of the directional cluster tree to be deleted. */
HEADER_PREFIX void
del_dcluster(pdcluster t);

/** @brief Create directions for a cluster tree.
 *
 *  @param depth Depth of the cluster tree.
 *  @param dim Spatial dimension.
 *  @returns Pointer to new @ref leveldir object. */
HEADER_PREFIX pleveldir
new_leveldir(uint depth, uint dim);

/** @brief Delete directions.
 *
 *  @param ld Directions to be deleted. */
HEADER_PREFIX void
del_leveldir(pleveldir ld);

/* ------------------------------------------------------------
 * Construction of dcluster trees
 * ------------------------------------------------------------ */

/** @brief Construct a directional cluster tree.
 *
 *  @param t Root of a standard cluster tree.
 *  @returns New directional cluster tree. */
HEADER_PREFIX pdcluster
buildfromcluster_dcluster(pccluster t);

/* ------------------------------------------------------------
 * Bounding box
 * ------------------------------------------------------------ */

/** @brief Compute the diameter of the bounding box of a directional cluster.
 *
 *  @param t Directional cluster @f$t@f$.
 *  @returns Euclidean diameter of the bounding box @f$B_t@f$. */
HEADER_PREFIX real
diam_dcluster(pcdcluster t);

/** @brief Compute the distance of the bounding boxes of two
 *     directional clusters.
 *
 *  @param t Directional cluster @f$t@f$.
 *  @param s Directional cluster @f$s@f$.
 *  @returns Euclidean distance of the bounding boxes @f$B_t@f$
 *     and @f$B_s@f$. */
HEADER_PREFIX real
dist_dcluster(pcdcluster t, pcdcluster s);

/** @brief Compute the distance of the midpoints of the bounding
 *     boxes of two directional clusters.
 *
 *  @param t Directional cluster @f$t@f$.
 *  @param s Directional cluster @f$s@f$.
 *  @returns Euclidean distance of the midpoints of the
 *     bounding boxes @f$B_t@f$ and @f$B_s@f$. */
HEADER_PREFIX real
middist_dcluster(pcdcluster t, pcdcluster s);

/* ------------------------------------------------------------
 * Match directions
 * ------------------------------------------------------------ */

/** @brief Among the directions associated with a directional cluster,
 *     find the one that best matches @f$\alpha d@f$.
 *
 *  @param t Directional cluster.
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param d Direction vector @f$d@f$.
 *  @returns Index of the direction associated with @f$t@f$ that
 *     best matches @f$\alpha d@f$. */
HEADER_PREFIX uint
finddirection_dcluster(pcdcluster t, real alpha, pcreal d);

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

/** @brief Get the number of active @ref dcluster objects.
 *
 *  @returns Number of currently active @ref dcluster objects. */
HEADER_PREFIX uint
getactives_dcluster();

/** @brief Compute the storage size of a directional cluster tree.
 *
 *  @param t Root of the directional cluster tree.
 *  @returns Storage size of the directional cluster tree. */
HEADER_PREFIX size_t
getsize_dcluster(pcdcluster t);

/** @brief Compute the depth of a directional cluster tree.
 *
 *  @param t Root of the directional cluster tree.
 *  @returns Depth of the directional cluster tree. */
HEADER_PREFIX uint
getdepth_dcluster(pcdcluster t);

/** @brief Compute the number of all directions for all clusters
 *     in a directional cluster tree.
 *
 *  @param t Root of the directional cluster tree.
 *  @returns Sum of the number of directions for all clusters
 *     in the tree. */
HEADER_PREFIX uint
getalldirections_dcluster(pcdcluster t);

/* ------------------------------------------------------------
 * Construct directions for clusters
 * ------------------------------------------------------------ */

/** @brief Construct set of directions for a directional
 *     cluster tree by using a subdivision of an axis-parallel box.
 *
 *  The function computes the maximal diameters of clusters on
 *  each level of the directional cluster tree.
 *  Then the directions for each level are constructed by subdividing
 *  the surface of the cube @f$[-1,1]^3@f$ into sufficiently small
 *  boxes and projecting their centers to the unit sphere.
 *
 *  @param t Root of the directional cluster tree.
 *  @param eta1 Directional admissibility parameter @f$\eta_1@f$.
 *     The smaller @f$\eta_1@f$ is, the more directions are
 *     created.
 *  @returns Directions for all levels of the cluster tree. */
HEADER_PREFIX pleveldir
builddirections_box_dcluster(pdcluster t, real eta1);

/** @brief Find an index best matching @f$\alpha d@f$ on
 *     a given level of a directional cluster tree.
 *
 *  @param ld Directions for all levels.
 *  @param l Level.
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param d Direction vector @f$d@f$.
 *  @returns Index of direction best matching @f$\alpha d@f$ on
 *     level <tt>l</tt> of the cluster tree. */
HEADER_PREFIX uint
finddirection_leveldir(pcleveldir ld, uint l, real alpha, pcreal d);

/** @} */

#endif
