
/* ------------------------------------------------------------
 This is the file "ddcluster.h" of the H2Lib package.
 All rights reserved, Nadine Albrecht 2015
 ------------------------------------------------------------ */

/** @file ddcluster.h
 *  @author Nadine Albrecht
 */

#ifndef DDCLUSTER_H
#define DDCLUSTER_H

/*CORE 0*/
#include "settings.h"

/*CORE 1*/
#include "sparsematrix.h"

/*CORE 2*/
#include "cluster.h"
#include "clustergeometry.h"
#include "hmatrix.h"

/** @defgroup ddcluster ddcluster
 *  @brief Routines for Domain Decomposition Clustering.
 *
 * This module contains routines for regular and adaptive domain decomposition
 * clustering.
 *  @{ */

/** @brief Build a @ref cluster tree from a @ref clustergeometry object using 
 *   regular Domain Decomposition Clustering.
 * 
 * Builds a regular cluster tree basing on the geometrical informations of 
 * the clustergeometry object using a domain decomposition method.
 * The indices are splitted into three sons, two of the sons (domain clusters)
 * contain the indices correponding to two disconnected sub domains, the third 
 * one contains the indices of the separator. Therefore the sparsematrix is 
 * useful.  
 * The splitting is regular. It starts with the forwarded direction and the 
 * following directions are chosen via cycling through all possible directions.
 * If the size of the index set is greater than the leaf size, the cluster is 
 * splitted, else it is a leaf cluster.
 * During this clustering the bounding boxes are updated.
 * 
 * @param cg Clustergeometry object with geometrical information.
 * @param size Number of indices.
 * @param idx Index set.
 * @param clf Maximal leaf size.
 * @param sp Sparsematrix, used for connectivity information.
 * @param dim Dimension of Clustering.
 * @param direction Direction for the next splitting step.
 * @param flag Auxiliary array, has to be initialised with zeros.
 * @return @ref cluster tree basing on regular domain decomposition clustering.
 */
HEADER_PREFIX pcluster 
build_regular_dd_cluster(pclustergeometry cg, uint size, uint *idx, uint clf,
    psparsematrix sp, uint dim, uint direction, uint *flag);

/** @brief Build a @ref cluster tree from a @ref clustergeometry object using 
 * adaptive Domain Decomposition Clustering.
 * 
 * Builds an adaptive cluster tree basing on the geometrical informations of 
 * the clustergeometry object using a domain decomposition method.
 * The indices are splitted into three sons, two of the sons (domain cluster)
 * contain the indices correponding to two disconnected sub domains, the third 
 * one contains the indices of the separator. Therefore the sparsematrix is 
 * useful.  
 * The splitting is adaptive in the direction with the largest spatial 
 * extension and the index set is splitted accordingly, if its size is greater 
 * than the leaf size, else the cluster is a leaf cluster.
 * During this clustering the bounding boxes are updated.
 * 
 * @param cg Clustergeometry object with geometrical information.
 * @param size Number of indices.
 * @param idx Index set.
 * @param clf Maximal leaf size.
 * @param sp Sparsematrix, used for connectivity information.
 * @param dim Dimension of Clustering.
 * @param flag Auxiliary array, has to be initialised with zeros.
 * @return @ref cluster tree basing on adaptive domain decomposition clustering.
 */
HEADER_PREFIX pcluster
build_adaptive_dd_cluster(pclustergeometry cg, uint size, uint *idx, uint clf,
    psparsematrix sp, uint dim, uint *flag);

/** @brief Check the domain decomposition admissibility condition.
 * 
 * Checks the domain decomposition admissibility of two @ref cluster trees
 * @f$ s @f$ and @f$ t @f$.
 * The block @f$ (s,t) @f$ is admissible, if @f$ s @f$ and @f$ t @f$ are two 
 * disjunct domain cluster or admissible_2_min-cluster @f$ (s, t)@f$ is TRUE.
 * 
 * @param s Row cluster.
 * @param t Col cluster.
 * @param data Has to be a pointer to real.
 * @return TRUE,if the block @f$ (s,t) @f$ is admissible, otherwise FALSE.
 */
bool admissible_dd_cluster(pcluster s, pcluster t, void* data);

/** @brief Build a @ref cluster tree from a @ref clustergeometry object basing on 
 *         the clustering of the elements.
 * 
 * Builds a cluster tree basing on the geometrical informations of 
 * the clustergeometry object and the clustering of the elements stored in the 
 * @ref cluster tree object <tt>ctri</tt>.
 * The indices are splitted into three sons, two of the sons (domain cluster)
 * contain the indices correponding to two disconnected sub domains, the third 
 * one contains the indices of the separator. Therefore the sparsematrix is 
 * useful.  
 * The splitting is corresponding to the splitting of the cluster tree 
 * <tt>ctri</tt>  and the index set is splitted accordingly, if <tt>ctri</tt>
 * has sons, else the cluster is a leaf cluster.
 * During this clustering the bounding boxes are updated.
 * 
 * @remark The @ref cluster tree object <tt>ctri</tt> should be built with one 
 *         of the standard routines from @ref cluster.h. 
 * 
 * @param ctri Cluster tree object.
 * @param cg Clustergeometry object with geometrical information.
 * @param size Number of indices.
 * @param idx Index set.
 * @param clf Maximal leaf size,
 * @param dim Dimension of Clustering.
 * @return @ref cluster tree basing on the clustering of the underlying cluster
           tree object  */
HEADER_PREFIX pcluster 
build_adaptive_fromcluster_cluster(pcluster ctri, pclustergeometry cg, uint size, uint *idx, 
					    uint clf, uint dim);

/** @} */

#endif
