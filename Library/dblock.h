
/* ------------------------------------------------------------
 * This is the file "dblock.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2015
 * ------------------------------------------------------------ */

/** @file dblock.h
 *  @author Steffen B&ouml;rm */

#ifndef DBLOCK_H
#define DBLOCK_H

/** @defgroup dblock dblock
 *  @brief Directional block tree.
 *  @{ */

/** @brief Directional block tree */
typedef struct _dblock dblock;

/** @brief Pointer to @ref dblock */
typedef dblock *pdblock;

/** @brief Pointer to constant @ref dblock */
typedef const dblock *pcdblock;

/** @brief Data for directional admissibility condition */
typedef struct _diradmdata diradmdata;

/** @brief Pointer to @ref diradmdata */
typedef diradmdata *pdiradmdata;

/** @brief Pointer to constant @ref diradmdata */
typedef const diradmdata *pcdiradmdata;

#include "dcluster.h"
#include "settings.h"

#include <stdlib.h>

#ifdef USE_CAIRO
#include <cairo.h>
#endif

/** @brief Directional block tree. */
struct _dblock
{
  /** @brief Row cluster */
  pdcluster rc;
  /** @brief Column cluster */
  pdcluster cc;

  /** @brief Row direction */
  uint rd;
  /** @brief Column direction */
  uint cd;

  /** @brief Admissibility flag */
  bool adm;

  /** @brief Number of row sons */
  uint rsons;

  /** @brief Number of column sons */
  uint csons;

  /** @brief Pointers to sons.
   *
   *  <tt>son[i+j*rsons]</tt> points to the son in row <tt>i</tt>
   *  and column <tt>j</tt>. */
  pdblock *son;

  /** @brief Number of descendants, including this block */
  uint desc;
};

/** @brief Data for directional admissibility condition */
struct _diradmdata
{
  /** @brief Angular admissibility parameter
   *
   *  This parameter controls the directional resolution.
   *  Small values correspond to higher resolutions, large values
   *  to lower resolutions. */
  real eta1;

  /** @brief Parabolic admissibility parameter
   *
   *  This parameter controls how far admissible blocks
   *  have to be from the singularity.
   *  Small values correspond to a large distance, fast convergence,
   *  and high storage requirements, while large values lead to
   *  slower convergence and lower storage requirements. */
  real eta2;

  /** @brief Wave number */
  real wave_k;

  /** @brief Auxiliary vector, used internally */
  preal xy;

  /** @brief Sets of directions for all levels of the cluster tree */
  pcleveldir ld;
};

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

/** @brief Create a new directional block.
 *
 *  @param rc Row cluster.
 *  @param cc Column cluster.
 *  @param rd Row direction.
 *  @param cd Column direction.
 *  @param rsons Number of row sons.
 *  @param csons Number of column sons.
 *  @returns New directional block. */
HEADER_PREFIX pdblock
new_dblock(pdcluster rc, pdcluster cc, uint rd, uint cd,
	   uint rsons, uint csons);

/** @brief Update a directional block.
 *
 *  This function has to be called after changing the sons
 *  of a directional block, e.g., to update the <tt>desc</tt>
 *  field in the current block.
 *
 *  @param b Directional block to be updated. */
HEADER_PREFIX void
update_dblock(pdblock b);

/** @brief Delete a directional block tree.
 *
 *  @param b Directional block tree to be delete. */
HEADER_PREFIX void
del_dblock(pdblock b);

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

/** @brief Obtain the number of currently active directional blocks.
 *
 *  @returns Number of currently active directional blocks. */
HEADER_PREFIX uint
getactives_dblock();

/** @brief Compute the storage size of a directional block tree.
 *
 *  @param b Directional block tree.
 *  @returns Storage size of the directional block tree with
 *    root <tt>b</tt>. */
HEADER_PREFIX size_t
getsize_dblock(pcdblock b);

/** @brief Compute the depth of a directional block tree.
 *
 *  @param b Directional block tree.
 *  @returns Depth of the directional block tree. */
HEADER_PREFIX uint
getdepth_dblock(pcdblock b);

/* ------------------------------------------------------------
 * Drawing
 * ------------------------------------------------------------ */

#ifdef USE_CAIRO
/** @brief Draw a directional block tree.
 *
 *  @param cr Cairo context for drawing.
 *  @param b Directional block tree.
 *  @param levels Number of levels to be drawn.
 *     Zero means that all levels are drawn. */
HEADER_PREFIX void
cairodraw_dblock(cairo_t *cr, pcdblock b, int levels);
#endif

/* ------------------------------------------------------------
 * Building standard block trees
 * ------------------------------------------------------------ */

/** @brief Build a directional block tree.
 *
 *  @param rc Root row cluster.
 *  @param cc Root column cluster.
 *  @param l Root level (for accessing <tt>leveldir</tt> structures.
 *  @param admissible Admissibility condition.
 *  @param data Data for admissibility condition.
 *  @returns New directional block tree. */
HEADER_PREFIX pdblock
build_dblock(pdcluster rc, pdcluster cc, uint l,
	     bool (*admissible)(pdcluster rc, pdcluster cc, uint l,
				uint *rd, uint *cd, void *data),
	     void *data);

/* ------------------------------------------------------------
 * Parabolic admissibility condition
 * ------------------------------------------------------------ */

/** @brief Parabolic admissibility condition.
 *
 *  A pair of clusters @f$(t,s)@f$ is considered admissible if
 *  @f$\kappa \diam^2(B_t) \leq \eta_2 \dist(B_t,B_s)@f$,
 *  @f$\kappa \diam^2(B_s) \leq \eta_2 \dist(B_t,B_s)@f$,
 *  @f$\diam(B_t) \leq \eta_2 \dist(B_t,B_s)@f$, and
 *  @f$\diam(B_s) \leq \eta_2 \dist(B_t,B_s)@f$, where
 *  @f$B_t@f$ and @f$B_s@f$ are the bounding boxes of the
 *  two clusters.
 *
 *  @param rc Row cluster @f$t@f$.
 *  @param cc Column cluster @f$s@f$.
 *  @param l Level.
 *  @param rd Pointer to where the row direction will be stored.
 *  @param cd Pointer to where the column direction will be stored.
 *  @param data Additional data, has to be a pointer to an
 *      @ref diradmdata object.
 *  @returns True if @f$(t,s)@f$ is admissible. */
HEADER_PREFIX bool
parabolic_admissibility(pdcluster rc, pdcluster cc, uint l,
			uint *rd, uint *cd, void *data);

/** @brief Standard admissibility condition.
 *
 *  A pair of clusters @f$(t,s)@f$ is considered admissible if
 *  @f$\diam(B_t) \leq \eta_2 \dist(B_t,B_s)@f$, and
 *  @f$\diam(B_s) \leq \eta_2 \dist(B_t,B_s)@f$, where
 *  @f$B_t@f$ and @f$B_s@f$ are the bounding boxes of the
 *  two clusters.
 *
 *  @param rc Row cluster @f$t@f$.
 *  @param cc Column cluster @f$s@f$.
 *  @param l Level.
 *  @param rd Pointer to where the row direction will be stored.
 *  @param cd Pointer to where the column direction will be stored.
 *  @param data Additional data, has to be a pointer to an
 *      @ref diradmdata object.
 *  @returns True if @f$(t,s)@f$ is admissible. */
HEADER_PREFIX bool
standard_admissibility(pdcluster rc, pdcluster cc, uint l,
		       uint *rd, uint *cd,
		       void *data);

/** @brief Compute the maximum of @f$\diam(t)/\dist(t,s)@f$
 *  and @f$\diam(s)/\diam(t,s)@f$ for all admissible blocks.
 *
 *  @param b Directional block.
 *  @returns Effective admissibility parameter @f$\eta@f$. */
HEADER_PREFIX real
getmaxeta_dblock(pcdblock b);

/* ------------------------------------------------------------
 * Remove directions
 * ------------------------------------------------------------ */

/** 
 *  @brief Remove unused directions and set up new @ref leveldir object.
 * 
 *  Search all used directions and collect them in a new 
 *  @ref leveldir object. After that the directional cluster and 
 *  block trees will be updated.
 *  
 *  @attention During this routine the @ref leveldir <tt> lold </tt> will be destroyed.
 *  
 *  @param b Root of the directional block tree.
 *  @param t Root of the directional cluster tree which belongs to <tt> b </tt>.
 *  @param lold Leveldir object which belongs to <tt> t </tt>.
 *  @returns New @ref leveldir object which is now used for <tt> t </tt>.
 */ 

HEADER_PREFIX pleveldir
remove_unused_direction(pdblock b, pdcluster t, pleveldir lold); 

/* ------------------------------------------------------------
 * Enumeration
 * ------------------------------------------------------------ */

/** 
 * @brief Enumerates a directional block tree.
 * 
 * Enumerates a @ref dblock object <tt> b </tt> in an array of size
 * @f$ b->desc @f$. Runs through all descendants and orders them in
 * pointer to a direction block tree.
 * 
 * @param b @ref dblock object to be enumerated.
 * @returns Returns a pointer to a @ref dblock object with enumerated
 * descendants and length of @f$ b->desc @f$.
 */

HEADER_PREFIX pdblock *
enumerate_dblock(pdblock b);

/** @} */

#endif
