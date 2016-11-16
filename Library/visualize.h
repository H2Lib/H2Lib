  /*-----------------------------------------------------
  This is the file "visualize.h" of the H2Lib package.
  All rights reserved, Christina Boerst 2016.
  -------------------------------------------------------*/
	 
  /**
   * @file visualize.h
   * @author Christina Boerst
   * @date 2016
   */
 
  #ifndef VISUALIZE_H_
  #define VISUALIZE_H_
  
  #ifdef USE_FREEGLUT
  #include <GL/glut.h>
  #include <GL/freeglut.h>
  #endif

  #include "dblock.h"
  #include "block.h"
  #include "dcluster.h"
  #include "clustergeometry.h"

  #include "tri2d.h"
  #include "tet3d.h"
  #include "surface3d.h"

/** \defgroup visualize visualize
 *  @brief This module contains functions to visualize bounding boxes of
 *  a given cluster, directional cluster, block, directional block further 
 *  functions to visualize a given two or three dimensional triangulation 
 *  or a surface triangulation with freeglut.
 *  As well as functions to visualize and animate a computed solution of
 *  the Helmholtz equation.
 *  On top of this, every window content can be saved as an image ppm or with 
 *  the additional library libpng as png.
 *  @{ */


 /**
 * @brief Function for visualize one or more level of a @ref dblock object. 
 *
 * This function will display three windows. 
 * One for selecting a block level, the others will show associated column and row clusters.
 * 
 * @param b The directional block cluster tree object to visualize.
 * @param grc Clustergeometry from directional row cluster. Needed if points should be drawn too.
 *        If not, <tt> grc </tt> should be NULL.
 * @param gcc Clustergeometry from directional column cluster. Needed if points should be drawn too.
 *        If not, <tt> gcc </tt> should be NULL.
 * @param c Bool for color gradient. 
 *        If c is true, different colors for every level of the given directional block cluster tree will be used. 
 * @param argc Number of arguments given by the main function.
 * @param argv Arguments content, also given by the main function.
 */
    HEADER_PREFIX void
    visualize_dblock_bbox(pcdblock b, pclustergeometry grc, pclustergeometry gcc, bool c, int argc, char **argv);
	
  
 /**
 * @brief Function for visualize one or more level of a @ref block object. 
 *
 * This function will display three windows. 
 * One for selecting a block level, the others will show associated column and row clusters.
 *
 * @param b The block cluster tree object to visualize.
 * @param grc Clustergeometry from directional row cluster. Needed if points should be drawn too.
 * If not, <tt>grc</tt> should be NULL.
 * @param gcc Clustergeometry from directional column cluster. Needed if points should be drawn too.
 * If not, <tt>gcc</tt> should be NULL.
 * @param c Bool for coloring.
 *        If c is true, different colors for every level of the given block cluster tree will be used.
 * @param argc Number of arguments given by the main function.
 * @param argv Arguments content, also given by the main function.
 */

  HEADER_PREFIX void
  visualize_block_bbox(pcblock b, pclustergeometry grc, pclustergeometry gcc, bool c, int argc, char **argv);

 /**
 * @brief Function for visualize a certain leaf of a @ref dblock object. 
 *
 * This function will display two windows. 
 * One for selecting a certain block out of all directional leaf blocks.
 * The associated two directional clusters will be displayed in the second window.
 *
 * @param b The directional block cluster tree object to visualize.
 * @param grc Clustergeometry from directional row cluster. Needed if points should be drawn too.
 *        If not, <tt> grc </tt> should be NULL.
 * @param gcc Clustergeometry from directional column cluster. Needed if points should be drawn too.
 *        If not, <tt> gcc </tt> should be NULL.
 * @param argc Number of arguments given by the main function.
 * @param argv Arguments content, also given by the main function.
 */
  
  HEADER_PREFIX void 
  visualize_dblock_certain_bbox(pcdblock b, pclustergeometry grc, pclustergeometry gcc, int argc, char **argv);
 
 /**
 * @brief Function for visualize a certain leaf of a @ref block object. 
 *
 * This function will display two windows. 
 * One for selecting a certain block out of all leaf blocks.
 * The associated two clusters will be displayed in the second window.
 *
 * @param b The block cluster tree object to visualize.
 * @param grc Clustergeometry from directional row cluster. Needed if points should be drawn too.
 *        If not, <tt> grc </tt> should be NULL.
 * @param gcc Clustergeometry from directional column cluster. Needed if points should be drawn too.
 *        If not, <tt> gcc </tt> should be NULL. 
 * @param argc Number of arguments given by the main function.
 * @param argv Arguments content, also given by the main function.
 */


  HEADER_PREFIX void 
  visualize_block_certain_bbox(pcblock b, pclustergeometry grc, pclustergeometry gcc, int argc, char **argv);

 /**
 * @brief Function for visualize a certain leaf block of a @ref block, only leaf blocks, whose level
 *        are between l1 and l2, are selectible. 
 *
 * This function will display two windows. 
 * One for selecting a certain block out of all leaf blocks, whose level
 * are between l1 and l2.
 * The associated two clusters will be displayed in the second window.
 *
 * @param b The block cluster tree object to visualize.
 * @param grc Clustergeometry from directional row cluster. Needed if points should be drawn too.
 *        If not, <tt> grc </tt> should be NULL.
 * @param gcc Clustergeometry from directional column cluster. Needed if points should be drawn too.
 *        If not, <tt> gcc </tt> should be NULL.
 * @param l1 Uint with the block level to start drawing.
 * @param l2 Uint with the block level to stop drawing.
 * @param argc Number of arguments given by the main function.
 * @param argv Arguments content, also given by the main function.
 */

  
  HEADER_PREFIX void   
  visualize_block_level_certain_bbox(pcblock b, pclustergeometry grc, pclustergeometry gcc, uint l1, uint l2, int argc, char **argv);

 
 /**
 * @brief Function for visualize a certain leaf block of a @ref dblock, only leaf blocks, whose level
 *        are between l1 and l2, are selectible. 
 *
 * This function will display two windows. 
 * One for selecting a certain block out of all directional leaf blocks, whose level
 * are between l1 and l2.
 * The associated two directional clusters will be displayed in the second window.
 *
 * @param b The directional block cluster tree object to visualize.
 * @param grc Clustergeometry from directional row cluster. Needed if points should be drawn too.
 *        If not, <tt> grc </tt> should be NULL.
 * @param gcc Clustergeometry from directional column cluster. Needed if points should be drawn too.
 *        If not, <tt> gcc </tt> should be NULL.
 * @param l1 Uint with the directional block level to start drawing.
 * @param l2 Uint with the directional block level to stop drawing.
 * @param argc Number of arguments given by the main function.
 * @param argv Arguments content, also given by the main function.
 */

  HEADER_PREFIX void   
  visualize_dblock_level_certain_bbox(pcdblock b, pclustergeometry grc, pclustergeometry gcc, uint l1, uint l2, int argc, char **argv);

 /**
 * @brief Function for visualize @ref dcluster objects. 
 *
 * This function will display two windows. 
 * One for selecting a cluster level and the other one to display all clusters corresponding to the choosen level.
 *
 * @param t The directional cluster tree object to visualize.
 * @param gt Clustergeometry from the directional cluster. Needed if points should be drawn, else <tt>gt == NULL</tt>.
 * @param c Bool for color gradient. 
 *        If c is true, different colors for every level of the given directional cluster tree will be used.
 * @param argc Number of arguments given by the main function.
 * @param argv Arguments content, also given by the main function.
 */
  
  
  HEADER_PREFIX void 
  visualize_dcluster_bbox(pcdcluster t, pclustergeometry gt, bool c, int argc, char **argv);

 /**
 * @brief Function for visualize @ref cluster objects. 
 *
 * This function will display two windows. 
 * One for selecting a cluster level and the other one to display all clusters corresponding to the choosen level.
 *
 * @param t The cluster tree object to visualize.
 * @param gt Clustergeometry from the cluster.  Needed if points should be drawn, else <tt>gt == NULL</tt>.
 * @param c Bool for color gradient. 
 *        If c is true, different colors for every level of the given cluster tree will be used.
 * @param argc Number of arguments given by the main function.
 * @param argv Arguments content, also given by the main function.
 */

    
  HEADER_PREFIX void 
  visualize_cluster_bbox(pccluster t, pclustergeometry gt, bool c, int argc, char **argv);

 /**
 * @brief Function for visualize a complete two dimensional 
 *        triangulation @ref tri2d. 
 *
 * @param tri The two dimensional triangulation @ref tri2d that should be visualized.
 * @param argc Number of arguments given by the main function.
 * @param argv Arguments content, also given by the main function.
 */	

  HEADER_PREFIX void
  visualize_tri2d(pctri2d tri, int argc, char **argv);
	
 /**
 * @brief Function for visualize a certain triangle
 * chosen out of the complete triangulation @ref tri2d. 
 * 
 * This function will display two windows.
 * One shows a table for all triangles also including their coordinates.
 * The other one marks the choosen triangle.
 *
 * @remark This functions is made for small triangulations with not more than about one
 *        hundret triangles, otherwise the table could be confusing.
 * @param tri The two dimensional triangulation @ref tri2d that should be visualized.
 * @param argc Number of arguments given by the main function.
 * @param argv Arguments content, also given by the main function.
 */	
	
  HEADER_PREFIX void
  visualize_certain_triangle(pctri2d tri, int argc, char **argv);

 /**
 * @brief Function for visualize a complete three dimensional 
 *        triangulation @ref tet3d. 
 *
 * @param tet The three dimensional triangulation @ref tet3d that should be visualized.
 * @param argc Number of arguments given by the main function.
 * @param argv Arguments content, also given by the main function.
 */	
		
  HEADER_PREFIX void
  visualize_tet3d(pctet3d tet, int argc, char **argv);

 /**
 * @brief Function for visualize a certain tetrahedron
 *        chosen out of the complete triangulation @ref tet3d. 
 * 
 * This function will display two windows.
 * One shows a table for all tetrahedron also including their coordinates.
 * The other one marks the choosen tetrahedron.
 *
 * @remark This functions is made for small triangulations with not more than about one
 *        hundret tetrahedron, otherwise the table could be confusing.
 * @param tet The three dimensional triangulation @ref tet3d that should be visualized.
 * @param argc Number of arguments given by the main function.
 * @param argv Arguments content, also given by the main function.
 */	 	
	
  HEADER_PREFIX void
  visualize_certain_tetrahedra(pctet3d tet, int argc, char **argv);  
	
 /**
 * @brief Function for visualize a complete surface 
 *        triangulation @ref surface3d. 
 *
 * @param gr The surface3d that should be visualized.
 * @param argc Number of arguments given by the main function.
 * @param argv Arguments content, also given by the main function.
 */	

  HEADER_PREFIX void
  visualize_surface3d(pcsurface3d gr, int argc, char **argv);

 /**
 * @brief Function for visualize a certain triangle
 * chosen out of the complete surface triangulation @ref surface3d. 
 * 
 * This function will display two windows.
 * One shows a table for all triangles also including their coordinates.
 * The other one marks the choosen triangle.
 *
 * @remark This functions is made for small triangulations with not more than about one
 *        hundred triangles, otherwise the table could be confusing.
 * @param gr The surface3d object that should be visualized.
 * @param argc Number of arguments given by the main function.
 * @param argv Arguments content, also given by the main function.
 */	 
	
  HEADER_PREFIX void
  visualize_certain_surface_triangle(pcsurface3d gr, int argc,char **argv);
  
  /**
   * @brief Function for visualize a boundary values for a @ref surface3d
   * object. 
   * 
   * This function will display two windows.
   * One shows the surface grid with the boundary values by
   * using a color gradient. 
   * The other one includes a color legend.
   *
   * @attention The dimension of the vector including boundary values has to
   * match with the number of triangles, edges or vertices of the 
   * @ref surface3d object.
   * 
   * @param sur The surface3d object that should be visualized.
   * @param val @ref avector object including boundary values. 
   * @param argc Number of arguments given by the main function.
   * @param argv Arguments content, also given by the main function.
   */	 
  
  HEADER_PREFIX void
  visualize_boundaryvalue_surface_triangle(pcsurface3d sur, pcavector val, int argc, char **argv);

  /**
   * @brief Function for visualize one or more solutions of the helmholtz equation 
   *        for a @ref surface3d object. 
   * 
   * This function will display two windows.
   * One shows the surface and a plane with values of a given solution @f$ v @f$.
   * The solution has to be computed on a plane produced by a two dimensional
   * grid and a constant @f$z @f$ value!
   * During the visualization the shown solution could be changed, if more than
   * one file was given.
   * The second window includes a color legend.
   * 
   * @attention The files with the solution has to be given in a special 
   * format and the points in the right order.
   * The First line has to have 3 values, namely:
   * @f[ \#\: points\: in\: x\: direction (unsigned\: int) \quad \#\: points\: in\: y\: direction
   *  (unsigned\: int) \quad z\: coordinate\: of\: the\: plane\: (double) @f]
   * All of the other lines have to be of the form:
   * @f[ \begin{matrix} x_{1}  &  y_{1}  & \Re{v((x_{1},y_{1},z)^{T})} & \Im{v((x_{1},y_{1},z)^{T})} & |v((x_{1},y_{1},z)^{T}) | \\
   *     x_{1}  &  y_{2} & \Re{v((x_{1},y_{2},z)^{T})} & \Im{v((x_{1},y_{2},z)^{T})} & |v((x_{1},y_{2},z)^{T}) | \\
   *     \vdots &  \vdots & \vdots &  \vdots & \vdots \\    
   *     x_{2} & y_{1} &  \Re{v((x_{2},y_{1},z)^{T})} & \Im{v((x_{2},y_{1},z)^{T})} & |v((x_{2},y_{1},z)^{T}) |\\
   *     \vdots &  \vdots & \vdots &  \vdots & \vdots\\   
   *     \end{matrix} @f]
   * where @f$v((x_{i}, y_{j}, z)^{T})@f$ is the value of the solution in @f$(x_{i},y_{j}, z)^{T}@f$.      
   * There should be no empty lines and in the end @f$x\cdot y@f$ three-dimensional 
   * points to visualize. 
   * 
   * @param sur The surface3d object that should be visualized.
   * @param argc Number of arguments given by the main function.
   * @param argv Argument content, in this case files called additionally.
   */
  
  HEADER_PREFIX void
  visualize_helmholtz_solution_surface(pcsurface3d sur, int argc, char **argv);
  
  /**
   * @brief Animates solutions of the helmholtz equation in time
   *        for a @ref surface3d object. 
   * 
   * This function will display two windows.
   * One visualize a xy-plane with the solution @f$ v(x) @f$ of the Helmholtz equation
   * (switching between real, imaginary and absolute values or fade in
   * the surface is possible ).
   * A time depending solution @f$ u(x,t) @f$ of the wave equation is 
   * computed with a separation approach @f$ u(x,t) = v(x) w(t) @f$ and
   * the function in time is given by
   * @f[
   * w(t) := exp(-i \omega t) \quad \text{with } \omega = \frac{2.0*\pi}{\text{frame rate}},
   * @f]
   * where the frame rate is a visualization depending fixed value. 
   * The solution has to be computed on a xy-plane produced by a two dimensional
   * grid and a constant @f$z @f$ value!
   * During the visualization the shown solution could be changed, if more than
   * one file was given.
   * The second window includes a color legend.
   * 
   * @attention The files with the solution have to be given in a special 
   * format and the points in the right order.
   * The First line has to have 3 values, namely:
   * @f[ \#\: points\: in\: x\: direction (unsigned\: int) \quad \#\: points\: in\: y\: direction
   *  (unsigned\: int) \quad z\: coordinate\: of\: the\: plane\: (double) @f]
   * All of the other lines have to be of the form:
   * @f[ \begin{matrix} x_{1}  &  y_{1}  & \Re{v((x_{1},y_{1},z)^{T})} & \Im{v((x_{1},y_{1},z)^{T})} & |v((x_{1},y_{1},z)^{T}) | \\
   *     x_{1}  &  y_{2} & \Re{v((x_{1},y_{2},z)^{T})} & \Im{v((x_{1},y_{2},z)^{T})} & |v((x_{1},y_{2},z)^{T}) | \\
   *     \vdots &  \vdots & \vdots &  \vdots & \vdots \\    
   *     x_{2} & y_{1} &  \Re{v((x_{2},y_{1},z)^{T})} & \Im{v((x_{2},y_{1},z)^{T})} & |v((x_{2},y_{1},z)^{T}) |\\
   *     \vdots &  \vdots & \vdots &  \vdots & \vdots\\   
   *     \end{matrix} ,@f]
   * where @f$v((x_{i}, y_{j}, z)^{T})@f$ is the value of the solution in
   * @f$(x_{i},y_{j}, z)^{T}@f$.      
   * There should be no empty lines and in the end @f$x\cdot y@f$  
   * points to visualize.
   * 
   * @param sur Corresponding surface3d object.
   * @param argc Number of arguments given by the main function.
   * @param argv Argument content, in this case files with solution values.
   */
  
  HEADER_PREFIX void
  animate_helmholtz_solution(pcsurface3d sur, int argc, char **argv);


/** @}*/
  #endif
