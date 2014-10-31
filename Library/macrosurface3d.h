
/* ------------------------------------------------------------
   This is the file "macrosurface3d.h" of the H2Lib package.
   All rights reserved, Steffen Boerm 2010
   ------------------------------------------------------------ */

/** @file macrosurface3d.h
 *  @author Steffen B&ouml;rm
 */

#ifndef MACROSURFACE3D_H
#define MACROSURFACE3D_H

typedef struct _macrosurface3d macrosurface3d;
typedef macrosurface3d *pmacrosurface3d;
typedef const macrosurface3d *pcmacrosurface3d;

#include "settings.h"
#include "surface3d.h"

struct _macrosurface3d {
  uint vertices;		/* Number of vertices */
  uint edges;			/* Number of edges */
  uint triangles;		/* Number of triangles */

  real (*x)[3];			/* Vertex coordinates */
  uint (*e)[2];			/* Edge vertices */
  uint (*t)[3];			/* Triangle vertices, counter-clockwise */
  uint (*s)[3];			/* Triangle edges, s[i][j] opposite t[i][j] */

  void (*phi)(uint i,		/* Parametrization */
	      real xr1, real xr2, void *data,
	      real xt[3]);
  void *phidata;
};  

/* ------------------------------------------------------------
   Constructor and destructor
   ------------------------------------------------------------ */

HEADER_PREFIX pmacrosurface3d
new_macrosurface3d(uint vertices, uint edges, uint triangles);

HEADER_PREFIX void
del_macrosurface3d(pmacrosurface3d mg);

/* ------------------------------------------------------------
   Examples
   ------------------------------------------------------------ */

HEADER_PREFIX pmacrosurface3d
new_sphere_macrosurface3d();

/* ------------------------------------------------------------
   Polygonal approximation
   ------------------------------------------------------------ */

HEADER_PREFIX psurface3d
build_from_macrosurface3d_surface3d(pcmacrosurface3d mg, uint split);

/* ------------------------------------------------------------
   Interactive setup of a surface3d object
   ------------------------------------------------------------ */

HEADER_PREFIX psurface3d
build_interactive_surface3d();

#endif
