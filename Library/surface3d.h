
/* ------------------------------------------------------------
 * This is the file "surface3d.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2010
 * ------------------------------------------------------------ */

/** @file surface3d.h
 *  @author Steffen B&ouml;rm */

#ifndef SURFACE3D_H
#define SURFACE3D_H

/**
 * @defgroup surface3d surface3d
 *
 * @brief Representation of a 3D-surface by using a triangle mesh.
 *
 * This module offers basic functions to create, load and save triangular
 * surface meshes. Also a simple red refinement is implemented via
 * @ref refine_red_surface3d .
 *
 * @{
 */

/**
 * @brief Abbreviation for the struct @ref _surface3d .
 */
typedef struct _surface3d surface3d;

/**
 * @brief Abbreviation for a pointer to a @ref surface3d object.
 */
typedef surface3d *psurface3d;

/**
 * @brief Abbreviation for a pointer to a constant @ref surface3d object.
 */
typedef const surface3d *pcsurface3d;

#include "settings.h"

/**
 * @brief Representation of a triangle surface mesh.
 */
struct _surface3d {
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
  /** @brief Triangle vertices, counter-clockwise */
  uint (*t)[3];
  /** @brief Triangle edges, s[i][j] opposite t[i][j] */
  uint (*s)[3];

  /** @brief Normal vectors */
  real (*n)[3];
  /** @brief Gram determinant, equals area doubled */
  preal g;

  /** @brief Minimal mesh size */
  real hmin;
  /** @brief Maximal mesh size */
  real hmax;
};

/* ------------------------------------------------------------
 * Constructor and destructor
 * ------------------------------------------------------------ */

/**
 * @brief create a new @ref surface3d object with a certain number of
 * vertices, edges and triangles.
 *
 * This function will allocated storage for the members <tt>x</tt>, <tt>e</tt>,
 * <tt>t</tt>, <tt>s</tt>, <tt>n</tt> and <tt>g</tt> of the struct
 * @ref _surface3d "surface3d".
 *
 * @param vertices Number of vertices for the new mesh.
 * @param edges Number of edges for the new mesh.
 * @param triangles Number of triangles for the new mesh.
 * @return A new, uninitialized @ref surface3d object.
 */
HEADER_PREFIX psurface3d
new_surface3d(uint vertices, uint edges, uint triangles);

/**
 * @brief This function computes the normal vectors <tt>n</tt> , the gram
 * determinants <tt>g</tt> and the minimal and maximal mesh size
 * <tt>hmin, hmax</tt> from the geometrical information.
 *
 * @param gr @ref surface3d "Surface3d" object to be initialized.
 */
HEADER_PREFIX void
prepare_surface3d(psurface3d gr);

/**
 * @brief Free Storage allocated for a @ref surface3d object.
 *
 * @param gr @ref surface3d "Surface3d" object to be deleted.
 */
HEADER_PREFIX void
del_surface3d(psurface3d gr);

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

/** @brief Compute geometrical properties of a surface grid.
 *
 *  @param gr Grid
 *  @param hmin Will be filled with minimal edge length.
 *  @param hmax Will be filled with maximal edge length.
 *  @param anglemin Will be filled with minimal interior angle.
 *  @param angleedge Will be will with maximal angle across edges
 *         between triangles. */
HEADER_PREFIX void
getproperties_surface3d(pcsurface3d gr,
		        preal hmin, preal hmax,
			preal anglemin, preal angleedge);

/* ------------------------------------------------------------
 * Debugging
 * ------------------------------------------------------------ */

/**
 * @brief print geometrical information to stdout.
 *
 * @param gr @ref surface3d "Surface3d" object to be printed out.
 */
HEADER_PREFIX void
print_surface3d(pcsurface3d gr);

/**
 * @brief Check if the geometrical information of a surface are consistent.
 *
 * @param gr @ref surface3d "Surface3d" object to be checked.
 * @return The total number of problems found is returned.
 */
HEADER_PREFIX uint
check_surface3d(pcsurface3d gr);

/**
 * @brief Check if the surface mesh is closed.
 *
 * @param gr @ref surface3d "Surface3d" object to be checked.
 * @return Return value is <tt>true</tt>, if surface is closed, false otherwise.
 */
HEADER_PREFIX bool
isclosed_surface3d(pcsurface3d gr);

/**
 * @brief Check if the surface is oriented correctly.
 *
 * @param gr @ref surface3d "Surface3d" object to be checked.
 @return Return value is <tt>true</tt>, if surface is oriented, false otherwise.
 */
HEADER_PREFIX bool
isoriented_surface3d(pcsurface3d gr);

/**
 * @brief Scale the geometry to a cube of given size.
 *
 * The coordinates of the geometry are scaled in a way to fit into the cube
 * defined by @f$ [a_1, b_1] \times [a_2, b_2] \times [a_3, b_3]@f$.
 *
 * @param gr @ref surface3d "Surface3d" object to be scaled.
 * @param a Minimal coordinates of scaling cube.
 * @param b Maximal coordinates of scaling cube.
 */
HEADER_PREFIX void
scale_surface3d(psurface3d gr, real *a, real *b);

/**
 * @brief Translate a geometry @p gr by a vector @f$t \in \mathbb R^3@f$.
 *
 * @param gr Geometry that should be translated.
 * @param t Translation vector @f$t@f$.
 */
HEADER_PREFIX void
translate_surface3d(psurface3d gr, real *t);

/**
 * @brief Merge to meshes into a single mesh.
 *
 * @param gr1 First mesh to be merged.
 * @param gr2 Second mesh to be merged.
 * @return New surface mesh containing all elements form @p gr1 and @p gr2.
 */
HEADER_PREFIX psurface3d
merge_surface3d(pcsurface3d gr1, pcsurface3d gr2);

/* ------------------------------------------------------------
 * File I/O
 * ------------------------------------------------------------ */

/**
 * @brief Write geometrical information of a surface mesh into a given file
 * using the H2Lib ascii representation.
 *
 * For the format description please refer to @ref read_surface3d .
 *
 * @see read_surface3d for the used file format.
 *
 * @param gr Geometry to be written to a file.
 * @param filename Filename for the geometry.
 */
HEADER_PREFIX void
write_surface3d(pcsurface3d gr, const char *filename);

/**
 * @brief Read geometrical information of a surface mesh from a given file
 * using the H2Lib ascii representation.
 *
 * The file format is the following:<br><br>
 * <code>
 * {vertices} {edges} {triangles}<br>
 * {x[0][0]} {x[0][1]} {x[0][2]}<br>
 * ...<br>
 * {x[{vertices}-1][0]} {x[{vertices}-1][1]} {x[{vertices}-1][2]}<br>
 * {e[0][0]} {e[0][1]}<br>
 * ...<br>
 * {e[{edges}-1][0]} {e[{edges}-1][1]}<br>
 * {t[0][0]} {t[0][1]} {t[0][2]}
 * {s[0][0]} {s[0][1]} {s[0][2]}<br>
 * ...<br>
 * {t[{triangles}-1][0]} {t[{triangles}-1][1]} {t[{triangles}-1][2]}
 * {s[{triangles}-1][0]} {s[{triangles}-1][1]} {s[{triangles}-1][2]}<br>
 * </code>
 *
 * @attention The normal vectors <tt>n</tt> and the gram determinants <tt>g</tt>
 * are not yet initialized. Consider calling @ref prepare_surface3d before using
 * the geometry read by this function.
 *
 * @param filename Filename for the geometry.
 * @return A new @ref surface3d object with the geometrical information read from
 * the file is returned.
 */
HEADER_PREFIX psurface3d
read_surface3d(const char *filename);

/**
 * @brief Write geometrical information of a surface mesh into a given file
 * using NetCDF.
 *
 * @param gr Geometry to be written to a file.
 * @param filename Filename for the geometry.
 */
HEADER_PREFIX void
write_nc_surface3d(pcsurface3d gr, const char *filename);

/**
 * @brief Read geometrical information of a surface mesh from a given file
 * using NetCDF.
 *
 * @attention The normal vectors <tt>n</tt> and the gram determinants <tt>g</tt>
 * are not yet initialized. Consider calling @ref prepare_surface3d before using
 * the geometry read by this function.
 *
 * @param filename Filename for the geometry.
 * @return A new @ref surface3d object with the geometrical information read from
 * the file is returned.
 */
HEADER_PREFIX psurface3d
read_nc_surface3d(const char *filename);

/**
 * @brief Read geometrical information of a surface mesh from a given file
 * using the netgen format.
 *
 * @attention The normal vectors <tt>n</tt> and the gram determinants <tt>g</tt>
 * are not yet initialized. Consider calling @ref prepare_surface3d before using
 * the geometry read by this function.
 *
 * @param filename Filename for the geometry.
 * @return A new @ref surface3d object with the geometrical information read from
 * the file is returned.
 */
HEADER_PREFIX psurface3d
read_netgen_surface3d(const char *filename);

/**
 * @brief Read geometrical information of a surface mesh from a given file
 * using the gmsh format.
 *
 * @attention The normal vectors <tt>n</tt> and the gram determinants <tt>g</tt>
 * are not yet initialized. Consider calling @ref prepare_surface3d before using
 * the geometry read by this function.
 *
 * @param filename Filename for the geometry.
 * @return A new @ref surface3d object with the geometrical information read from
 * the file is returned.
 */
HEADER_PREFIX psurface3d
read_gmsh_surface3d(const char *filename);

/**
 * @brief Read geometrical information of a surface mesh from a given file
 * using the unv format.
 *
 * @attention The normal vectors <tt>n</tt> and the gram determinants <tt>g</tt>
 * are not yet initialized. Consider calling @ref prepare_surface3d before using
 * the geometry read by this function.
 *
 * @param filename Filename for the geometry.
 * @return A new @ref surface3d object with the geometrical information read from
 * the file is returned.
 */
HEADER_PREFIX psurface3d
read_unv_surface3d(char *filename);

/* ------------------------------------------------------------
 * Mesh refinement
 * ------------------------------------------------------------ */

/**
 * @brief Apply a red refinement to a surface mesh.
 *
 * The surface mesh is globally refined with a red refinement resulting in a
 * mesh having four times as much triangles as the input mesh.
 *
 * @param in Surface to be refined.
 * @return A new @ref surface3d object is returned which is a red refinement from
 * the input surface mesh.
 */
HEADER_PREFIX psurface3d
refine_red_surface3d(psurface3d in);

/**
 * @}
 */

#endif
