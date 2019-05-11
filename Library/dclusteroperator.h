/* ------------------------------------------------------------
 * This is the file "dclusteroperator.h" of the H2Lib package.
 * All rights reserved, Christina Boerst 2015/16
 * ------------------------------------------------------------ */

/** @file dclusteroperator.h
 *  @author Christina Boerst
 */
 
 #ifndef DCLUSTEROPERATOR_H
 #define DCLUSTEROPERATOR_H
 
 /** @defgroup dclusteroperator dclusteroperator
 *   @brief Representation of directional cluster operators used to describe
 *   transformations of directional cluster bases.
 *
 *  @{ */
 
 /** @brief Representation of a directional cluster operator */
 typedef struct _dclusteroperator dclusteroperator;
 
 /** @brief Pointer to @ref dclusteroperator object. */
 typedef dclusteroperator *pdclusteroperator;
 
 /** @brief Pointer to constant @ref dclusteroperator object. */
 typedef const dclusteroperator *pcdclusteroperator;
 
 #include "dcluster.h"
 #include "dclusterbasis.h"
 #include "amatrix.h"

 /** @brief Representation of a directional cluster operator */
 struct _dclusteroperator
 {
   /** @brief Corresponding directional cluster. */
   pcdcluster t;
   
   /** @brief Number of rows for each direction. */
   uint *krow;
   
   /** @brief Number of columns for each direction. */
   uint *kcol;
   
   /** @brief Number of directions. */
   uint dir;

   /** @brief Coefficient matrix for each direction. */
   pamatrix C;
   
   /** @brief Number of sons, either <tt>t->sons</tt> or zero. */
   uint sons;
   
   /** @brief Pointers to sons. */
   pdclusteroperator *son;
   
   /** @brief References to this directional cluster operator. */
   uint refs;   

 }; 
 
 
 
 /** 
  *  @brief Initialize a @ref dclusteroperator object.
  *
  *  Sets up the components of the object.
  *  If <tt> t </tt> is not a leaf, the array <tt> son </tt> is allocated,
  *  otherwise it is set to null.
  *
  *  @remark Should always be matched by a call to @ref uninit_dclusteroperator.
  *
  *  @param co @ref dclusteroperator object to be initialized.
  *  @param t Corresponding directional cluster.
  *  @returns Initialized @ref dclusteroperator object. 
  */
 
 HEADER_PREFIX pdclusteroperator
   init_dclusteroperator(pdclusteroperator co, pcdcluster t);
 
 
 /** 
 *  @brief Initialize a @ref dclusteroperator object for a leaf.
 *
 *  Sets up the components of the object.
 *  Sets <tt> son </tt> to null.
 *  If <tt> t->sons>0 </tt>, this yields a partial directional cluster operator.
 *
 *  @remark Should always be matched by a call to @ref uninit_dclusteroperator.
 *
 *  @param co @ref dclusteroperator object to be initialized.
 *  @param t Corresponding directional cluster.
 *  @returns Initialized @ref dclusteroperator object. 
 */
 
 HEADER_PREFIX pdclusteroperator
   init_leaf_dclusteroperator(pdclusteroperator co, pcdcluster t);
 
 
 /** 
 * @brief Uninitializes a @ref dclusteroperator object.
 *
 *  Invalidates pointers, freeing corresponding storage if appropriate,
 *  and prepares the object for deletion.
 *
 *  @param co @ref dclusteroperator object to be uninitialized. 
 */
 
 HEADER_PREFIX void
   uninit_dclusteroperator(pdclusteroperator co); 
 
 
 /** 
 * @brief Create a new @ref dclusteroperator object.
 *
 *  Allocates storage for the object and sets up its components.
 *
 *  @remark Should always be matched by a call to @ref del_dclusteroperator.
 *
 *  @param t corresponding directional cluster.
 *  @return Returns the new @ref dclusteroperator object.
 */
 
 HEADER_PREFIX pdclusteroperator
   new_dclusteroperator(pcdcluster t);
 
 /** 
 *  @brief Set a pointer to a @ref dclusterbasis object, increase its
 *  reference counter, and decrease reference counter of original
 *  pointer target.
 *
 *  @param ptr Pointer to the @ref pdclusteroperator variable that will
 *         be changed.
 *  @param co @ref dclusteroperator that will be referenced. 
 */
 
 HEADER_PREFIX void
   ref_dclusteroperator(pdclusteroperator *ptr, pdclusteroperator co);
 
 /** 
 *  @brief Reduce the reference counter of a @ref dclusteroperator object.
 *  
 *  If the reference counter reaches zero, the object is deleted.
 *
 *  @remark Use @ref ref_dclusteroperator with <tt>cb=NULL</tt> instead,
 *  since this guarantees that the pointer is properly invalidated.
 *
 *  @param co @ref dclusteroperator that will be unreferenced. 
 */
 
 HEADER_PREFIX void
   unref_dclusteroperator(pdclusteroperator co);  
 
 /** 
 *  @brief Creates a new @ref dclusteroperator object for a leaf.
 *  
 *  Allocates storage for the object and sets up its components.
 *
 *  @remark Should always be matched by a call to @ref del_dclusteroperator.
 *
 *  @param t Corresponding directional cluster.
 *  @return Returns the new @ref dclusteroperator object.
 */
 
 HEADER_PREFIX pdclusteroperator
   new_leaf_dclusteroperator(pcdcluster t);
 
 /** 
 *  @brief Delete a @ref dclusteroperator object.
 *
 *  Releases the storage corresponding to the object.
 *
 *  @param co @ref dclusteroperator object to be deleted. 
 */
 
 HEADER_PREFIX void
   del_dclusteroperator(pdclusteroperator co);
 
 /** 
 *  @brief Get the number of active @ref dclusteroperator objects.
 *
 *  @returns Number of active @ref dclusteroperator objects. 
 */
 
 HEADER_PREFIX uint
   getactives_dclusteroperator();
 
 
 /** 
 *  @brief Change the number of rows and columns of a
 *  directional cluster operator and resize <tt>cb->C[i]</tt> accordingly.
 *
 *  @param co Directional cluster operator that will be changed.
 *  @param krow New number of rows.
 *  @param kcol New number of columns.
 *  @param i Direction, which should be changed. 
 */
 
 HEADER_PREFIX void
   resize_dclusteroperator(pdclusteroperator co, uint krow, uint kcol, uint i);
 
 
 /** 
  * @brief Builds a @ref dclusteroperator object matching a @ref dcluster tree.
  *
  *  Constructs a directional cluster operator for a directional cluster and its
  *  descendants.
  *
  *  Row and column numbers will be set to zero.
  *
  *  @param t Root of the @ref dcluster.
  *  @returns @ref dclusteroperator for <tt> t </tt> and its descendants. 
  */
 
 HEADER_PREFIX pdclusteroperator
   build_from_dcluster_dclusteroperator(pcdcluster t);
 
 
 /** 
  *  @brief Builds a @ref dclusteroperator object matching a
  *  @ref dclusterbasis object.
  *
  *  Constructs a directional cluster operator for a directional
  *  cluster basis and its descendants.
  *
  *  Row numbers will be set to zero, column numbers will be set
  *  to the ranks of the corresponding directional cluster basis.
  *
  *  @param cb Root of the @ref dclusterbasis.
  *  @returns @ref dclusteroperator for <tt>cb</tt> and its descendants. 
  */
 
 HEADER_PREFIX pdclusteroperator
 build_from_dclusterbasis_dclusteroperator(pcdclusterbasis cb); 
 
 /**  
  *  @brief Multiplies the matrices from one directional cluster operator
  *  with the corresponding matrices from another one.
  *   
  *  For every cluster <tt> t </tt> and direction @f$ \iota @f$ 
  *  the matrix @f$ C_{t\iota} @f$ from the @ref dclusteroperator 
  *  <tt> co2 </tt> will be multiplied by the corresponding matrix from
  *  the @ref dclusteroperator <tt> co </tt> and stored in <tt> co2 </tt>.
  *  
  *  Obviously both @ref dclusteroperator objects have to belong to the
  *  same directional cluster basis and the matrix
  *  sizes have to match.
  *  
  *  @param co Directional cluster operator, only used for the multiplication.
  *  @param co2 Directional cluster operator, will be multiplied by the matrix 
  *         from <tt> co </tt> and overwritten by the product.
  */
 
 HEADER_PREFIX void
 merge_dclusteropertator(pcdclusteroperator co, pdclusteroperator co2);
 
 /** 
  *  @brief Print the tree sturcture of a given @ref dclusteroperator object.
  *
  *  @param co Root of the @ref dclusteroperator. 
  */
 
 HEADER_PREFIX void
 print_tree_dclusteroperator(pcdclusteroperator co);
 
 /* --------------------------------------- */
 /* Enumeration                             */
 /* --------------------------------------- */
 
 /** 
  *  @brief Enumerates all directional cluster operators according to a directional cluster tree.
  *
  *  @param t Directional cluster tree.
  *  @param co Directional cluster operator belongs to <tt> t </tt>. 
  *  @returns Directional cluster operator containing numeration.
  */
 
 HEADER_PREFIX pdclusteroperator *
   enumerate_dclusteroperator(pcdcluster t, pdclusteroperator co);
 
 /** @} */
 
#endif
