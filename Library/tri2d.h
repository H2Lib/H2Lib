
/* ------------------------------------------------------------
 * This is the file "tri2d.h" of the H2Lib package.
 * All rights reserved, Dirk Boysen 2015
 * ------------------------------------------------------------ */

/** @file tri2d.h
 *  @author Dirk Boysen */

#ifndef TRI2D_H
#define TRI2D_H

/** @defgroup tri2d tri2d
 * @brief Two-dimensional triangular meshes.
 * @{ */

/** @brief Two-dimensional triangulat mesh*/
typedef struct _tri2d tri2d;

/** @brief Pointer to @ref tri2d object*/
typedef tri2d *ptri2d;

/** @brief Pointer to constant @ref tri2d object*/
typedef const tri2d *pctri2d;

/** Refinement data for triangular meshes*/
typedef struct _tri2dref tri2dref;

/** @brief Pointer to @ref tri2dref object*/
typedef tri2dref *ptri2dref;

/** @brief Pointer to constant @ref tri2dref object*/
typedef const tri2dref *pctri2dref;

/** @brief Tool for constructing @ref tri2d meshes*/
typedef struct _tri2dbuilder tri2dbuilder;

/** Pointer to @ref tri2dbuilder object*/
typedef tri2dbuilder *ptri2dbuilder;

#include "settings.h"

/* ------------------------------------------------------------
   Triangular mesh,
   described by vertices, edges and triangles
   ------------------------------------------------------------ */

/** @brief Representation of a two-dimensional triangular mesh.
 *
 *  Triangular meshes are represented by a hierarchy of geometric
 *  objects: vertices are at the lowest level, edges consist of
 *  vertices and triangles consist of edges.
 *  Functions like @ref getvertices_tri2d are provided to skip 
 *  levels of the hierarchy in order to obtain the vertices or
 * edges of a tetrahedron. */
struct _tri2d {
  /** @brief Number of vertices*/
  uint vertices;
  
  /** @brief Number of edges*/
  uint edges;
  
  /** @brief Number of triangles*/
  uint triangles;

  /** @brief coordinates of vertices*/
  real (*x)[2];
  
  /** @brief Start and end points of edges*/
  uint (*e)[2];			
  
  /** @brief Edge of a triangle*/
  uint (*t)[3];			
  
  /** @brief Boundary flags for vertices*/
  uint *xb;			

  /** @brief Boundary flags for edges*/
  uint *eb;		
};

/* ------------------------------------------------------------
   Constructors and destructors
   ------------------------------------------------------------ */

/** @brief Create a partially initialized @ref tri2d mesh.
 *
 *  All arrays are created, but are not filled with meaningful
 *  values.
 *
 *  @param vertices Number of vertices
 *  @param edges Number of edges
 *  @param triangles Number of tetrahedra
 *  @returns Partially initialized @ref tri2d mesh */

HEADER_PREFIX ptri2d
new_tri2d(uint vertices, uint edges, uint triangles);

/** @brief Delete a @ref tri2d object
 *
 *  @param t2 Object to be deleted */
HEADER_PREFIX void
del_tri2d(ptri2d t2);

/** @brief Create a mesh containing two triangles representing
 *  the unit square, @f$[-1,1] \times [-1,1] @f$.
 *
 *  @returns @ref tri2d mesh representing the unit square. */
HEADER_PREFIX ptri2d
new_unitsquare_tri2d();

/** @brief Create a mesh containing four triangles representing
 *  the unit square, @f$[-1,1] \times [-1,1] @f$.
 *
 *  @returns @ref tri2d mesh representing the unit square. */
HEADER_PREFIX ptri2d
new_unitcircle_tri2d();

/** @brief Create a mesh containing six triangles representing
 *  a L-shaped domain with vertices @f$(-1,-1), (0,-1), (1,-1), 
 *  (-1,0), (0,0), (1,0), (-1,1) and (0,1)@f$.
 *
 *  @returns @ref tri2d mesh representing the L-shaped domain */
HEADER_PREFIX ptri2d
new_lshape_tri2d();

/** @brief Create a mesh containing twentyfour triangles representing
 *  a U-shaped domain.
 *
 *  @returns @ref tri2d mesh representing the U-shape. */
HEADER_PREFIX ptri2d
new_ushape_tri2d();

/* ------------------------------------------------------------
   File I/O
  -------------------------------------------------------------*/

/** @brief Write a triangular mesh to a file.
 *
 *  @param t2 Source mesh
 *  @param name Filename */
HEADER_PREFIX void
write_tri2d(pctri2d t2, const char *name);

/** @brief Read a triangular mesh from a file.
 *
 *  @param name Filename
 *  @returns @ref tri2d object describing the mesh */
HEADER_PREFIX ptri2d
read_tri2d(const char *name);

/* ------------------------------------------------------------
   Get geometrical information
   ------------------------------------------------------------ */
   

/** @brief Find the vertices of a triangle.
 * 
 *  Due to the hierarchical representation of geometric objects,
 *  we cannot access the vertices of a triangle directly.
 *  This function finds the vertices and returns them in a
 *  specific order: <tt>v[i]</tt> is the vertex opposite
 *  the <tt>i</tt>-th edge.
 * @param t2 Mesh
 * @param tn Number of the triangle
 * @param v Will be filled with vertex numbers.
 * */
HEADER_PREFIX void
getvertices_tri2d(pctri2d t2, uint tn, uint v[]);

/* ------------------------------------------------------------
   Check structure for inconsistencies
   ------------------------------------------------------------ */

/** @brief Ensure that the vertices of boundary edges are in
 *  counter-clockwise order as seen from outside of the mesh.
 *
 *  @param t2 Mesh, offending boundary edges will be reordered
 *  @returns Number of edges that have been fixed. */
HEADER_PREFIX void
fixnormals_tri2d(ptri2d t2);


/** @brief Perform various consistency checks.
 *
 *  @param t2 Mesh
 *  @returns Number of inconsistencies found */
HEADER_PREFIX void
check_tri2d(pctri2d t2);


/* ------------------------------------------------------------
   Regular refinement
   ------------------------------------------------------------ */

/** @brief Representation of the refinement relationship between
 *  two meshes.
 *
 *  This structure stores for each vertex, edge and triangle
 *  how it was created during the refinement algorithm:
 *  - a vertex can be a copy of a coarse vertex or lie in a coarse edge
 *  - an edge can lie in an coarse edge or in a coarse triangle */
struct _tri2dref {
  /** @brief Father index of a vertex*/
  uint *xf;
  /** @brief Type of father for a vertex: 0 means vertex, 1 means edge*/
  uint *xt;			

  /** @brief Father index of an edge*/
  uint *ef;			
  /** @brief Type of father for an edge: 1 means edge, 2 means triangle*/
  uint *et;	

  /** @brief Father index of a triangle*/
  uint *tf;			
};

/** @brief Regular refinement of a triangular mesh.
 *
 *  The mesh is globally refined using regular refinement.
 *
 *  @param t2 Coarse mesh
 *  @param t2r If not null, the target pointer will be set to a new
 *        @ref tri2dref object describing the refinement relationship
 *  @returns Refined mesh */
HEADER_PREFIX ptri2d
refine_tri2d(pctri2d t2, ptri2dref *t2r);

/** @brief Delete a @ref tri2dref object.
 *
 *  @param t2r Object to be deleted */
HEADER_PREFIX void
del_tri2dref(ptri2dref t2r);

/* ------------------------------------------------------------
   Display grid
   ------------------------------------------------------------ */

/** @brief Draw a tri2d mesh
 * 
 * Draw a @ref tri2d mesh.
 * 
 * @param t2 This grid will be drawn.
 * @param filename In this fill the mesh will be drawn.
 * @param mark_refedges
 * @param mark_triangle 
 */

HEADER_PREFIX void
draw_cairo_tri2d(pctri2d t2, const char *filename,
	   bool mark_refedges, int mark_triangle);


/* -----------------------------------------------------------
  ------------------------------------------------------------ */
/** @brief Smooth a @ref tri2d unitcircle
 * 
 * Smooth a @ref tri2d mesh representing a unit circle.
 * The vertices of the smooth unit circle are shifted to the boundary.
 * 
 * @param t2 mesh representing the unitsquare, will be overwritten
 * with smooth mesh.
 */
HEADER_PREFIX void
smooth_unitcircle_tri2d(ptri2d t2);

/* ------------------------------------------------------------
 * Tool for constructing meshes
 * ------------------------------------------------------------ */

/** @brief Create a new @ref tri2dbuilder object.
 *
 *  @param vertices Number of vertices of the new mesh.
 *  @returns @ref tri2dbuilder object with no edges or triangles. */
HEADER_PREFIX ptri2dbuilder
new_tri2dbuilder(uint vertices);

/** @brief Delete a @ref tri2dbuilder object.
 *
 *  @param tb Object to be deleted. */
HEADER_PREFIX void
del_tri2dbuilder(ptri2dbuilder tb);

/** @brief Obtain array of vertex coordinates in @ref tri2dbuilder object.
 *
 *  Since @ref tri2dbuilder is an opaque class, this function is provided
 *  for setting the coordinates of the vertices.
 *
 *  @param tb @ref tri2dbuilder object.
 *  @returns Pointer to array of with <tt>tb->vertices</tt> entries of
 *     type <tt>real [2]</tt> representing the vertex coordinates. */
HEADER_PREFIX real
(*getx_tri2dbuilder(ptri2dbuilder tb))[2];

/** @brief Add a triangle to a @ref tri2dbuilder object.
 *
 *  @param tb Target @ref tri2dbuilder object.
 *  @param v0 Index of first vertex.
 *  @param v1 Index of second vertex.
 *  @param v2 Index of third vertex. */
HEADER_PREFIX void
addtriangle_tri2dbuilder(ptri2dbuilder tb,
			 uint v0, uint v1, uint v2);

/** @brief Create a @ref tri2d mesh from the geometrical and topological
 *  information stored in a @ref tri2dbuilder object.
 *
 *  Once the mesh has been created, the corresponding @ref tri2dbuilder
 *  object can be deleted.
 *
 *  @param tb Source @ref tri2dbuilder object.
 *  @returns @ref tri2d mesh. */
HEADER_PREFIX ptri2d
buildmesh_tri2dbuilder(ptri2dbuilder tb);

/** @}*/
#endif
