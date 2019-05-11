/* ------------------------------------------------------------
 This is the file "macrosurface3d.h" of the H2Lib package.
 All rights reserved, Steffen Boerm 2010
 ------------------------------------------------------------ */

/** @file macrosurface3d.h
 *  @author Steffen B&ouml;rm
 */

#ifndef MACROSURFACE3D_H
#define MACROSURFACE3D_H

/** @defgroup macrosurface3d macrosurface3d
 *  @brief Representation of a surface.
 * 
 *  The @ref macrosurface3d class is used to represent piecewise
 *  parametrized surfaces in three-dimensional space.
 *  These representations can be used to construct @ref surface3d
 *  meshes for boundary element discretizations.
 *
 *  @{*/

/** @brief Representation of a parametrized surface */
typedef struct _macrosurface3d macrosurface3d;

/** @brief Pointer to a @ref macrosurface3d object */
typedef macrosurface3d *pmacrosurface3d;

/** @brief Pointer to a constant @ref macrosurface3d object */
typedef const macrosurface3d *pcmacrosurface3d;

#include "settings.h"
#include "surface3d.h"

/** @brief Representation of a parametrized surface.
 *
 *  This class is intended as a simple mesh generator for boundary
 *  element applications.
 *  A surface is described by a (hopefully small) set of triangles
 *  associated with parametrizations.
 *  
 *  The topology is described by vertices, edges, and triangles.
 *  The parametrization callback <tt>phi</tt> uses these parameters
 *  to describe the geometry.
 *
 *  The function @ref build_from_macrosurface3d_surface3d can then
 *  be used to create a standard @ref surface3d mesh. */
struct _macrosurface3d {
  /** @brief Number of vertices */
  uint vertices;

  /** @brief Number of edges */
  uint edges;

  /** @brief Number of triangles */
  uint triangles;

  /** @brief Vertex coordinates */
  real (*x)[3];

  /** @brief Edge vertices */
  uint (*e)[2];

  /** @brief Triangle vertices.
   *
   *  The vertices are ordered counter-clockwise as seen from outside
   *  of the geometry, i.e., taking the cross product of
   *  <tt>x[t[i][1]]-x[t[i][0]]</tt> and <tt>x[t[i][2]]-x[t[i][0]]</tt>
   *  is supposed to yield an outer normal vector of the (unparametrized)
   *  triangle. */
  uint (*t)[3];

  /** @brief Triangle edges.
   *
   *  Edge <tt>s[i][j]</tt> lies opposite of vertex <tt>t[i][j]</tt>. */
  uint (*s)[3];

  /** @brief Parametrization callback.
   *
   *  Given a triangle index <tt>i</tt> and coordinates <tt>xr1</tt>
   *  and <tt>xr2</tt> in the reference triangle
   *  @f$\widehat{T}=\{ (s,t)\ :\ 0\leq s,t\leq 1,\ s+t\leq 1 \}@f$,
   *  this function returns the point @f$\Phi_i(\widehat{x})@f$
   *  in the parametrized triangle in the array <tt>xt</tt>.
   *
   *  Usually, the callback function will use the geometric information
   *  in <tt>x</tt>, <tt>e</tt>, <tt>t</tt>, and <tt>s</tt> to
   *  accomplish its task.
   *  Additional information can be provided via the pointer
   *  <tt>phidata</tt>.
   *
   *  It is up to the user to ensure that the vertices and edges
   *  of a refined triangulation produced by evaluating @f$\Phi_i@f$ in
   *  the vertices and along the edges of the macro triangulation
   *  in adjacent triangles match. */
  void (*phi)(uint i, real xr1, real xr2, void *phidata, real xt[3]);

  /** @brief Pointer that will be passed to the <tt>phi</tt> callback
   *  function. */
  void *phidata;
};

/* ------------------------------------------------------------
 * Constructor and destructor
 * ------------------------------------------------------------ */

/** @brief Create a @ref macrosurface3d object.
 *
 *  @param vertices Number of vertices.
 *  @param edges Number of edges.
 *  @param triangles Number of triangles.
 *  @returns Empty @ref macrosurface3d object. */
HEADER_PREFIX pmacrosurface3d
new_macrosurface3d(uint vertices, uint edges, uint triangles);

/** @brief Delete a @ref macrosurface3d object.
 *
 *  @param mg Object to be deleted. */
HEADER_PREFIX void
del_macrosurface3d(pmacrosurface3d mg);

/* ------------------------------------------------------------
 * Examples
 * ------------------------------------------------------------ */

/** @brief Create a @ref macrosurface3d object for the unit sphere.
 *
 *  @returns A @ref macrosurface3d representation of the unit sphere. */
HEADER_PREFIX pmacrosurface3d
new_sphere_macrosurface3d();

/**
 * @brief Creates a new @ref _macrosurface3d "macrosurface3d" object for a
 * parabolic mirror.
 *
 * The geometry is derived from @ref new_sphere_macrosurface3d by deforming
 * half of the sphere parabolically towards the other have of the geometry.
 *
 * @return A @ref macrosurface3d "macrosurface3d" representation of a
 * parabolic mirror.
 */
HEADER_PREFIX pmacrosurface3d
new_parabolic_mirror_macrosurface3d();

/**
 * @brief Creates a new @ref _macrosurface3d "macrosurface3d" object for a
 * cuboid.
 *
 * The cuboid is described as
 * @f[
 * [a_x, b_x] \times [a_y, b_y] \times [a_z, b_z].
 * @f]
 *
 * @param ax Minimal extend in x-direction.
 * @param bx Maximal extend in x-direction.
 * @param ay Minimal extend in y-direction.
 * @param by Maximal extend in y-direction.
 * @param az Minimal extend in z-direction.
 * @param bz Maximal extend in z-direction.
 * @return A @ref macrosurface3d "macrosurface3d" representation of a cuboid.
 */
HEADER_PREFIX pmacrosurface3d
new_cuboid_macrosurface3d(real ax, real bx, real ay, real by, real az, real bz);

/**
 * @brief Creates a new @ref _macrosurface3d "macrosurface3d" object for a unit
 * cube.
 *
 * @return A @ref macrosurface3d "macrosurface3d" representation of a cube.
 */
HEADER_PREFIX pmacrosurface3d
new_cube_macrosurface3d();

/** @brief Create a @ref macrosurface3d object for a cylinder.
 *
 *  @attention This geometry is not fully tested yet, there might be some
 *  error in generating the mesh at some refinement levels.
 *
 *  @returns A @ref macrosurface3d representation of a cylinder. */
HEADER_PREFIX pmacrosurface3d
new_cylinder_macrosurface3d();

/* ------------------------------------------------------------
 * Polygonal approximation
 * ------------------------------------------------------------ */

/** @brief Create a triangular mesh from a @ref macrosurface3d object.
 *
 *  Each triangle in <tt>mg</tt> is refined by splitting its edges
 *  into <tt>split</tt> parts and connecting the vertices to obtain
 *  <tt>split*split</tt> plane triangles for each parametrized triangle.
 *  The resulting mesh is returned in a @ref surface3d object.
 *
 *  @param mg @ref macrosurface3d representation of the surface.
 *  @param split Number of edge refinements.
 *  @returns @ref surface3d mesh consisting of plane triangles
 *     approximating the surface. */
HEADER_PREFIX psurface3d
build_from_macrosurface3d_surface3d(pcmacrosurface3d mg, uint split);

/* ------------------------------------------------------------
 * Interactive setup of a surface3d object
 * ------------------------------------------------------------ */

/** @brief Create a @ref surface3d object interactively.
 *
 *  Current this function can only read a mesh from a file or
 *  create one for the unit sphere.
 *
 *  @returns User-defined @ref surface3d mesh. */
HEADER_PREFIX psurface3d
build_interactive_surface3d();

/** @} */

#endif
