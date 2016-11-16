
/* ------------------------------------------------------------
 * This is the file "tet3d.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2015
 * ------------------------------------------------------------ */

/** @file tet3d.h
 *  @author Steffen B&ouml;rm */

#ifndef TET3D_H
#define TET3D_H

/** @defgroup tet3d tet3d
 *  @brief Three-dimensional tetrahedral meshes.
 *  @{ */

/** @brief Three-dimensional tetrahedral mesh */
typedef struct _tet3d tet3d;

/** @brief Pointer to @ref tet3d object */
typedef tet3d *ptet3d;

/** @brief Pointer to constant @ref tet3d object */
typedef const tet3d *pctet3d;

/** @brief Refinement data for tetrahedral meshes */
typedef struct _tet3dref tet3dref;

/** @brief Pointer to @ref tet3dref object */
typedef tet3dref *ptet3dref;

/** @brief Pointer to constant @ref tet3dref object */
typedef const tet3dref *pctet3dref;

/** @brief Tool for constructing @ref tet3d meshes */
typedef struct _tet3dbuilder tet3dbuilder;

/** @brief Pointer to @ref tet3dbuilder object */
typedef tet3dbuilder *ptet3dbuilder;

#include "settings.h"

/* ------------------------------------------------------------
   Tetrahedral mesh,
   described by vertices, edges, faces and tetrahedra
   ------------------------------------------------------------ */

/** @brief Representation of a three-dimensional tetrahedral mesh.
 *
 *  Tetrahedral meshes are represented by a hierarchy of geometric
 *  objects: vertices are at the lowest level, edges consist of
 *  vertices, faces consist of edges, and tetrahedra consist of faces.
 *  Functions like @ref getvertices_tet3d or @ref getedges_tet3d
 *  are provided to skip levels of the hierarchy in order to
 *  obtain the vertices or edges of a tetrahedron. */
struct _tet3d {
  /** @brief Number of vertices */
  uint vertices;

  /** @brief Number of edges */
  uint edges;

  /** @brief Number of faces */
  uint faces;

  /** @brief Number of tetrahedra */
  uint tetrahedra;

  /** @brief Coordinates of vertices */
  real (*x)[3];

  /** @brief Start and end points of edges */
  uint (*e)[2];

  /** @brief Edges on a triangular face.
   *
   *  If this is a boundary face, we assume that the faces are oriented
   *  counter-clockwise as seen from the outside of the mesh.
   *  This property is ensured by the function @ref fixnormals_tet3d
   *  and is used in the computation of outer normal vectors. */
  uint (*f)[3];			

  /** @brief Faces of a tetrahedron */
  uint (*t)[4];			

  /** @brief Boundary flags for vertices */
  uint *xb;

  /** @brief Boundary flags for edges */
  uint *eb;
  
  /** @brief Boundary flags for faces */
  uint *fb;
};

/* ------------------------------------------------------------
   Constructors and destructors
   ------------------------------------------------------------ */

/** @brief Create a partially initialized @ref tet3d mesh.
 *
 *  All arrays are created, but are not filled with meaningful
 *  values.
 *
 *  @param vertices Number of vertices
 *  @param edges Number of edges
 *  @param faces Number of faces
 *  @param tetrahedra Number of tetrahedra
 *  @returns Partially initialized @ref tet3d mesh */
HEADER_PREFIX ptet3d
new_tet3d(uint vertices, uint edges, uint faces, uint tetrahedra);

/** @brief Delete a @ref tet3d object
 *
 *  @param gr Object to be deleted */
HEADER_PREFIX void
del_tet3d(ptet3d gr);

/** @brief Create a mesh for a tetrahedron with three edges
 *  aligned with the coordinate axes.
 *
 *  @returns @ref tet3d mesh representing the tetrahedron. */
HEADER_PREFIX ptet3d
new_axis_tet3d();

/** @brief Create a mesh for a regular tetrahedron.
 *
 *  @returns @ref tet3d mesh representing a regular tetrahedron. */
HEADER_PREFIX ptet3d
new_regular_tet3d();

/** @brief Create a mesh containing six tetrahedra representing
 *  the unit cube @f$[0,1] \times [0,1] \times [0,1]@f$.
 *
 *  @returns @ref tet3d mesh representing the unit cube. */
HEADER_PREFIX ptet3d
new_unitcube_tet3d();

/* ------------------------------------------------------------
   File I/O
   ------------------------------------------------------------ */

/** @brief Write a tetrahedral mesh to a file.
 *
 *  @param gr Source mesh
 *  @param name Filename */
HEADER_PREFIX void
write_tet3d(pctet3d gr, const char *name);

/** @brief Read a tetrahedral mesh from a file.
 *
 *  @param name Filename
 *  @returns @ref tet3d object describing the mesh */
HEADER_PREFIX ptet3d
read_tet3d(const char *name);

/* ------------------------------------------------------------
   Get geometrical information
   ------------------------------------------------------------ */

/** @brief Find the vertices of a tetrahedron corresponding to
 *  a given face.
 *
 *  Due to the hierarchical representation of geometric objects,
 *  we cannot access the vertices of a tetrahedron directly.
 *  This function finds the vertices and returns them in a
 *  specific order: <tt>v[0]</tt> is the vertex opposite
 *  the local face <tt>fl</tt> of the tetrahedron, and
 *  <tt>v[1]</tt>, <tt>v[2]</tt> and <tt>v[3]</tt> are the
 *  vertices on this face such that <tt>v[i+1]</tt> is opposite
 *  the <tt>i</tt>-th edge <tt>gr->f[gr->t[tn][fl]][i]</tt>.
 *  This feature is important if the order of the edges is
 *  relevant, e.g., if we want to compute outward normal vectors
 *  for a boundary face.
 *
 *  @param gr Mesh
 *  @param tn Number of the tetrahedron
 *  @param fl Local index of the face, can only be 0, 1, 2, or 3.
 *  @param v Will be filled with the vertex numbers */
HEADER_PREFIX void
getvertices_byface_tet3d(pctet3d gr, uint tn, uint fl,
			 uint v[]);

/** @brief Find the vertices of a tetrahedron.
 *
 *  Due to the hierarchical representation of geometric objects,
 *  we cannot access the vertices of a tetrahedron directly.
 *  This function finds the vertices and returns them in a
 *  specific order: <tt>v[i]</tt> is the vertex opposite
 *  the <tt>i</tt>-th face.
 *
 *  @param gr Mesh
 *  @param tn Number of the tetrahedron
 *  @param v Will be filled with the vertex numbers */
HEADER_PREFIX void
getvertices_tet3d(pctet3d gr, uint tn, uint v[]);

/** @brief Find the edges of a tetrahedron.
 *
 *  Due to the hierarchical representation of geometric objects,
 *  we cannot access the edges of a tetrahedron directly.
 *  This function finds the edges and returns them in a
 *  specific order:
 *  - <tt>e[0]</tt> is the intersection of the faces 1 and 2.
 *  - <tt>e[1]</tt> is the intersection of the faces 2 and 3.
 *  - <tt>e[2]</tt> is the intersection of the faces 3 and 1.
 *  - <tt>e[3]</tt> is the intersection of the faces 0 and 1.
 *  - <tt>e[4]</tt> is the intersection of the faces 0 and 2.
 *  - <tt>e[5]</tt> is the intersection of the faces 0 and 3.
 *
 *  @param gr Mesh
 *  @param tn Number of the tetrahedron
 *  @param e Will be filled with the edge numbers */
HEADER_PREFIX void
getedges_tet3d(pctet3d gr, uint tn, uint e[]);

/** @brief Find the vertices of a face.
 * 
 *  Due to the hierarchical representation of geometric objects,
 *  we cannot access the vertices of a face directly.
 *  This function finds the vertices and returns them.
 * 
 *  @param t3 Mesh.
 *  @param nf Global number of the face.
 *  @param v will be filled with the vertex numbers.                    */ 
HEADER_PREFIX void
getvertices_face_tet3d(pctet3d t3, uint nf, uint v[]);

/* ------------------------------------------------------------
   Check structure for inconsistencies
   ------------------------------------------------------------ */

/** @brief Ensure that the edges of boundary faces are in
 *  counter-clockwise order as seen from outside of the mesh.
 *
 *  @param gr Mesh, offending boundary faces will be reordered
 *  @returns Number of faces that have been fixed. */
HEADER_PREFIX uint
fixnormals_tet3d(ptet3d gr);

/** @brief Perform various consistency checks.
 *
 *  @param gr Mesh*/
HEADER_PREFIX void
check_tet3d(pctet3d gr);

/** @brief Compute various statistices of a mesh.
 *
 *  @param gr Mesh
 *  @param hmin Will be overwritten by the minimal edge length
 *  @param hmax Will be overwritten by the maximal edge length
 *  @param volmin Will be overwritten by the minimal volume of a tetrahedron
 *  @param volmax Will be overwritten by the maximal volume of a tetrahedron
 *  @param relvolmin Will be overwritten by the minimal relative volume of a tetrahedron (relative with respect to the local mesh parameter)
 *  @param relvolmax Will be overwritten by the maximal relative volume of a tetrahedron (relative with respect to the local mesh parameter) */
HEADER_PREFIX void
statistics_tet3d(pctet3d gr,
		 preal hmin, preal hmax,
		 preal volmin, preal volmax,
		 preal relvolmin, preal relvolmax);

/* ------------------------------------------------------------
 * Regular refinement by Bey's algorithm
 * ------------------------------------------------------------ */

/** @brief Representation of the refinement relationship between
 *  two meshes.
 *
 *  This structure stores for each vertex, edge, face, and tetrahedron
 *  how it was created during the refinement algorithm:
 *  - a vertex can be a copy of a coarse vertex or lie in a coarse edge
 *  - an edge can lie in an coarse edge, a coarse face, or in a coarse tetrahedron
 *  - a face can lie in a coarse face or in a coarse tetrahedron */
struct _tet3dref {
  /** @brief Father index of a vertex */
  uint *xf;
  /** @brief Father type for a vertex: 0 means vertex, 1 means edges */
  uint *xt;

  /** @brief Father index of an edge */
  uint *ef;
  /** @brief Father type for an edge: 1 means edge, 2 means face, 3 means tetrahedron */
  uint *et;

  /** @brief Father index of a face */
  uint *ff;
  /** @brief Father type for a face: 2 means face, 3 means tetrahedron */
  uint *ft;

  /** @brief Father index of a tetrahedron */
  uint *tf;
};

/** @brief Regular refinement of a tetrahedral mesh.
 *
 *  The mesh is globally refined using Bey's algorithm.
 *
 *  @param gr Coarse mesh
 *  @param grr If not null, the target pointer will be set to a new
 *    @ref tet3dref object describing the refinement relationship
 *  @returns Refined mesh */
HEADER_PREFIX ptet3d
refine_tet3d(pctet3d gr, ptet3dref *grr);

/** @brief Delete a @ref tet3dref object.
 *
 *  @param grr Object to be deleted */
HEADER_PREFIX void
del_tet3dref(ptet3dref grr);

/* ------------------------------------------------------------
 * Tool for constructing meshes
 * ------------------------------------------------------------ */

/** @brief Create a new @ref tet3dbuilder object.
 *
 *  @param vertices Number of vertices of the new mesh.
 *  @returns @ref tet3dbuilder object with no edges, faces, or triangles. */
HEADER_PREFIX ptet3dbuilder
new_tet3dbuilder(uint vertices);

/** @brief Delete a @ref tet3dbuilder object.
 *
 *  @param tb Object to be deleted. */
HEADER_PREFIX void
del_tet3dbuilder(ptet3dbuilder tb);

/** @brief Obtain array of vertex coordinates in @ref tet3dbuilder object.
 *
 *  Since @ref tet3dbuilder is an opaque class, this function is provided
 *  for setting the coordinates of the vertices.
 *
 *  @param tb @ref tet3dbuilder object.
 *  @returns Pointer to array of with <tt>tb->vertices</tt> entries of
 *     type <tt>real [3]</tt> representing the vertex coordinates. */
HEADER_PREFIX real
(*getx_tet3dbuilder(ptet3dbuilder tb))[3];

/** @brief Add a tetrahedron to a @ref tet3dbuilder object.
 *
 *  @param tb Target @ref tet3dbuilder object.
 *  @param v0 Index of first vertex.
 *  @param v1 Index of second vertex.
 *  @param v2 Index of third vertex.
 *  @param v3 Index of fourth vertex. */
HEADER_PREFIX void
addtetrahedron_tet3dbuilder(ptet3dbuilder tb,
			    uint v0, uint v1, uint v2, uint v3);

/** @brief Create a @ref tet3d mesh from the geometrical and topological
 *  information stored in a @ref tet3dbuilder object.
 *
 *  Once the mesh has been created, the corresponding @ref tet3dbuilder
 *  object can be deleted.
 *
 *  @param tb Source @ref tet3dbuilder object.
 *  @returns @ref tet3d mesh. */
HEADER_PREFIX ptet3d
buildmesh_tet3dbuilder(ptet3dbuilder tb);

/** @} */

#endif
