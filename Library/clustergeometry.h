
/* ------------------------------------------------------------
 * This is the file "clustergeometry.h" of the H2Lib package.
 * All rights reserved, Knut Reimer 2009
 * ------------------------------------------------------------ */

/**
 * @file clustergeometry.h
 * @author Knut Reimer
 */

#ifndef CLUSTERGEOMETRY_H
#define CLUSTERGEOMETRY_H

/**
 * @defgroup clustergeometry clustergeometry
 *  @brief Representation of a cluster geometry object used to
 *  build @ref cluster trees.
 * 
 *  The @ref clustergeometry class is used to represent geometrical
 *  structures of a domain and is an auxiliary structure to build
 *  @ref cluster trees.
 *  @{*/

/** @brief Representation of a clustergeometry.*/
typedef struct _clustergeometry clustergeometry;

/** @brief Pointer to @ref clustergeometry object.*/
typedef clustergeometry* pclustergeometry;

#include "settings.h"
#include "cluster.h"
#include "amatrix.h"
#include "eigensolvers.h"

/**
 * @brief Representation of a clustergeometry object.
 *
 *  This auxiliary structure allocates storage for different kinds of
 *  clustering strategies.
 *  Characteristic points and weights for the indices could be stored, 
 *  additionally minimal and maximal coordinates for the support
 *  bounding boxes.
 *  Also storage for internal fields is allocated. */
struct _clustergeometry {
  /** @brief Spatial dimension. */
  uint dim;

  /** @brief Number of indices.*/
  uint nidx;

  /**  @brief Characteristic points for the indices.*/
  real **x;

  /** @brief Minimal coordinates for the support bounding boxes for the indices.*/
  real **smin;

  /** @brief Maximal coordinates for the support bounding boxes for the indices.*/
  real **smax;

  /** @brief Weights for the indices. */
  real *w;

  /** @brief internal fields used to build the @ref cluster tree from the 
   clustergeometry object.*/
  real *hmin;

  /** @brief internal fields used to build the @ref cluster tree from the
   clustergeometry object.*/
  real *hmax;

  /** @brief internal fields used to build the @ref cluster tree from the
   clustergeometry object.*/
  real *buf;
};

/* ------------------------------------------------------------
 Constructors and destructors
 ------------------------------------------------------------ */

/**
 * @brief Create a new @ref clustergeometry object.
 * 
 * Allocates storage for the object and the internal fields.
 * 
 * @remark Should always be matched by a call to @ref del_clustergeometry.
 * 
 * @param dim Spatial dimension of the domain.
 * @param nidx Number of characteristic points in the domain.
 * @return Returns the newly created @ref clustergeometry object.
 */
HEADER_PREFIX pclustergeometry
new_clustergeometry(uint dim, uint nidx);

/**
 * @brief Delete a @ref clustergeometry object.
 *
 * Releases the storage corresponding to this object.
 * 
 *@param cf Object to be deleted. */
HEADER_PREFIX void
del_clustergeometry(pclustergeometry cf);

/* ------------------------------------------------------------
 Auxiliary routines
 ------------------------------------------------------------ */

/**
 * @brief Update an adaptive bounding box for an index set.
 * 
 *  Computes an adaptive bounding box for an index set and stores the result
 *  in the internal fields <tt>hmin</tt> and <tt>hmax</tt> of the 
 *  @ref clustergeometry structure.
 *  Already existing entries are overwritten.
 * 
 *  @param cf clustergeometry object, where the bounding box is stored.
 *  @param size Number of indices.
 *  @param idx Index set. */
HEADER_PREFIX void
update_point_bbox_clustergeometry(pclustergeometry cf, uint size, uint *idx);

/**
 * @brief Update a bounding box for the support of a cluster.
 * 
 *  Computes a bounding box for the support of a cluster <tt>t</tt> using only 
 *  the fields <tt>smin</tt> and <tt>smax</tt> and the index set <tt>idx</tt>
 *  of the @ref clustergeometry structure.
 *  Already existing entries are overwritten.
 * 
 *  @param cf Clustergeometry object with geometrical information.
 *  @param t In this @ref cluster tree the bounding boxes are updated. */
HEADER_PREFIX void
update_support_bbox_cluster(pclustergeometry cf, pcluster t);

/**
 * @brief Updates the bounding boxes in a @ref cluster tree object using
 *  only its sons.
 * 
 *  Computes the coordinates of the support bounding boxes of a cluster using 
 *  only the bounding boxes of its sons.
 *  Already existing values are overwritten.
 * 
 *  @param t In this cluster the bounding boxes are updated. */
HEADER_PREFIX void
update_bbox_cluster(pcluster t);

#endif

/** @}*/
