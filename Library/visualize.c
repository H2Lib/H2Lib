/* ------------------------------------------------------------
   This is the file "visualize.c" of the H2Lib package.
   All rights reserved, Christina Boerst 2015-16
   ------------------------------------------------------------ */

/**
 * @file visualize.c
 * @author Christina Boerst
 * @date 2016
 */

#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "basic.h"
#include "avector.h"
#include "realavector.h"
#include "dblock.h"
#include "block.h"
#include "dcluster.h"
#include "cluster.h"
#include "clustergeometry.h"
#include "tri2d.h"
#include "tet3d.h"
#include "surface3d.h"

#include "visualize.h"

#ifdef USE_LIBPNG
#include <png.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef USE_FREEGLUT

  /*Needed values for cluster motion and 
     projection are put together in activemovement */

typedef struct _activemovement activemovement;

typedef activemovement *pactivemovement;


struct _activemovement {

  /* mouse position values, needed for evaluating rotation angle */
  int       old_position_x;
  int       old_position_y;
  /* angles for object rotation */
  float     angle_x;
  float     angle_y;
  float     old_angle_x;
  float     old_angle_y;
  /* values for object translation */
  float     rightleft;
  float     updown;
  /* value for zooming in and out */
  float     zoom;
  /* bool, which is true, if finding block father is possible */
  bool      levelup;
  /* choosen directional block son */
  pcdblock dblock;
  /* choosen block son */
  pcblock block;
};

  /*Needed values for orthogonal projection 
     and finding selected block are put 
     together in orthomovement */

typedef struct _orthomovement orthomovement;

typedef orthomovement *porthomovement;


struct _orthomovement {

  /* scaling values from current window size and both cluster trees */
  float     orthoscaling;
  float     clustscaling_x;
  float     clustscaling_y;
  /* for rectangles */
  float     xoff;
  float     yoff;
  /* position of target block */
  float     target_x;
  float     target_y;
  /* empty space between window edge and object */
  float     offset_x;
  float     offset_y;

};

  /*Values for working with tables are put together
     in tablemovement */

typedef struct _tablemovement tablemovement;

typedef tablemovement *ptablemovement;

struct _tablemovement {

  /* scaling values from current window size */
  float     tablescaling;
  uint      pixel_x;
  uint      pixel_y;
  /* values for finding choosen triangle */
  uint      target;
  uint      page;
  uint      max_page;
  /* empty space between window edge and object */
  uint      offset_x;
  uint      offset_y;

};

  /* Struct used for a grid with solution values */
typedef struct _grid grid;

typedef grid *pgrid;

struct _grid {

  double  **point;
  uint      x;
  uint      y;

  uint      z;
};

  /* Struct collecting solutions */
typedef struct _solution solution;

typedef solution *psolution;

struct _solution {
  /* Minimum and maximum value of solution(s) */
  real    **min;
  real    **max;
  /* Bools for visualize plane with solution and model separately */
  bool      plane;
  bool      model;
  /* Vector for boundary values */
  prealavector values;
  /* Grid for solution values on a plane */
  pgrid    *gr;
  uint      current_grid_number;
  /* Real, imaginary or absolut part */
  uint      type;
  uint      grids;

  /* Time and time_step for the animation */
  bool      animate;
  real      time;
  real      time_step;
  pgrid    *timegr;
  real      scale;

};

   /*
    *   Globale variable        *       
    */

static pcdcluster t, s;		/* Cluster and block structures */
static pcdblock bl;
static pclustergeometry gt, gs;
static pccluster ct, cs;
static pcblock bbl;

static ptri2d tr;		/* Triangles and tetrahedra */
static pctri2d tr_orig;
static ptet3d th;
static pctet3d th_orig;
static psurface3d gr;
static pcsurface3d gr_orig;

static uint level, level2;	/* Color and level managment */
static bool color;
static bool version;		/* To switch between cluster end block mode */
static bool points = 1;		/* With or without points */
static bool visu = 0;		/* For showing content only after choosing a level */
static uint up;			/* Value needed to switch between drawing him or his father */
static bool direction = 0;	/* Value needed to switch between drawing directions or not */
static uint mode = 0;		/* Modes for different drawing rules */
static uint coord = 0;		/* Coordinate system */
static bool light_coord[3] = { 0, 0, 0 };	/* To light a certain coordinate */

static bool origin = 0;
static bool bbox = 0;		/* Print coordinates */
static uint animation = 0;	/* Animation value */

   /*
    *   GLUT windows    *
    *        and            *
    *      movement       *
    */

static uint ortho, clust1, clust2, one, table;
static pactivemovement move_1, move_2;
static porthomovement move_o;
static ptablemovement move_t;
static psolution sol;

   /*
    *   Values for  *
    *   snapshoots  */

#define MAX_FILE  100		/* Maximal number of snapshots */

static const GLenum format = GL_RGBA;	/* Wanted format for the pixels */
static const GLuint format_bytes = 4;	/* Needed bytes for this format */
static uint snapshots = 0;	/* Actual number of snapshots */
static bool take_a_shot[2] = { 0, 0 };	/* Boolean to decide whether it should takes the pixels or not */


   /* 
    * Constructors and destructors *
    *                              */

pactivemovement new_activemovement(bool level);

porthomovement new_orthomovement(void);

ptablemovement new_tablevement(void);

pgrid     new_grid(uint x, uint y, uint wide);

psolution new_solution(uint vectorsize, uint number, bool animated);

void
          del_activemovement(pactivemovement move);

void
          del_orthomovement(porthomovement move);

void
          del_tablemovement(ptablemovement move);

void
          del_grid(pgrid p);

void
          del_solution(psolution sol);


   /*
    *    System dependent       *
    *   transformation of       *
    coordinates           */

static void
system_transformation()
{				/* To avoid nasty behaviour of thin clients....-.-'  */

  /* glRotated(180, 1.0, 0.0, 0.0);
   */
}

   /*
    * Snapshot routines   *
    */

   /* Simple version for producing ppm data */
static void
save_current_content_ppm(char *prefix, uint paper_width, uint paper_height,
			 uint color_max, GLubyte * pixels)
{

  uint      i, j, current;
  char      filename[MAX_FILE];

  snprintf(filename, MAX_FILE, "%s%d.ppm", prefix, snapshots);	/* Create ongoing numbering */
  FILE     *file = fopen(filename, "w");
  fprintf(file, "P3\n%d %d\n%d\n", paper_width, paper_height, color_max);

  for (i = 0; i < paper_height; i++) {	/* Write pixel values in file */
    for (j = 0; j < paper_width; j++) {
      current = format_bytes * ((paper_height - i - 1) * paper_width + j);
      fprintf(file, "%3d %3d %3d ", pixels[current], pixels[current + 1],
	      pixels[current + 2]);
    }
    fprintf(file, "\n");
  }
  fclose(file);
}

  /* More complicated routine to create png data */
#ifdef USE_LIBPNG
static void
save_current_content_png(char *prefix, uint paper_width, uint paper_height,
			 GLubyte * pixels, png_byte ** png_bytes,
			 png_byte *** png_rows)
{

  uint      i, ld;
  char      filename[MAX_FILE];

  snprintf(filename, MAX_FILE, "%s%d", prefix, snapshots);	/* Create a file and open it */
  FILE     *file = fopen(filename, "wb");
  ld = format_bytes * paper_width * paper_height;
  *png_bytes = realloc(*png_bytes, sizeof(png_byte) * ld);
  *png_rows = realloc(*png_rows, sizeof(png_byte *) * paper_height);

  for (i = 0; i < ld; i++) {
    (*png_bytes)[i] = pixels[i];	/* Write our pixel in the png data */
  }

  /* And structure it row-wise */
  for (i = 0; i < paper_height; i++) {
    (*png_rows)[paper_height - i - 1] =
      &(*png_bytes)[i * paper_width * format_bytes];
  }

  /* Allocate and init png_struct to write in */
  png_structp png_snap =
    png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (!png_snap) {		/* If something goes wrong */
    printf(" No png_struct data could be allocated! \n");
    abort();
  }

  /* Allocate and init png_info object for the png_snap */
  png_infop info = png_create_info_struct(png_snap);
  if (!info) {
    printf(" No png_info data could be allocated! \n");
    png_destroy_write_struct(&png_snap, (png_infopp) NULL);
    abort();
  }
  /* Okay, now the error handling.
     If an error turns up, libpng will call setjmp(),
     so if this happen, we are out of buisness */
  if (setjmp(png_jmpbuf(png_snap))) {
    printf(" Error! File will be closed! \n");
    png_destroy_write_struct(&png_snap, &info);
    fclose(file);
    abort();
  }

  /* Now the more important part =) 
     Standard output way for opened files */
  png_init_io(png_snap, file);

  /* Set up for writing */
  png_set_IHDR(png_snap, info, paper_width, paper_height, 8,
	       PNG_COLOR_TYPE_RGBA, PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

  /* Now we can start */
  png_write_info(png_snap, info);	/* Info part */
  png_write_image(png_snap, *png_rows);	/* The image */
  png_write_end(png_snap, NULL);	/* Stops filling */

  fclose(file);

  /* Cleaning up */
  png_destroy_write_struct(&png_snap, &info);
}
#endif

static int
read_file(const char *file, double scale, uint number, psolution sol,
	  bool animate)
{

  FILE     *in;
  char      buf[80], *line;
  uint      max_x, max_y;
  double    max[3], min[3], z;
  double    lenght, maxpoint;
  int       args;
  uint      j;
  longindex l, i;

  in = fopen(file, "r");
  if (in == NULL) {
    fprintf(stderr, "No file");
  }

  /* Read size for points */
  line = fgets(buf, 80, in);

  if (line == 0) {
    fprintf(stderr, "Could not read first line of \"%s\"\n", file);
    return -1;
  }

  args = sscanf(line, "%d %d %lf", &max_x, &max_y, &z);

  if (args != 3) {
    fprintf(stderr,
	    "First line of \"%s\" doesn't contain required numbers!\n", file);
    fprintf(stderr, "It should be of the form \n");
    fprintf(stderr,
	    "# points x direction (uint)  # points y direction (uint)  z-value of the plane (double)  wave number (double)\n");
    return -1;
  }

  /* Create new point object */
  sol->gr[number] = new_grid(max_x, max_y, 8);
  l = (longindex) max_x *max_y;
  sol->gr[number]->z = z;

  /* Read point positions */

  for (i = 0; i < l; i++) {
    line = fgets(buf, 80, in);
    if (line == 0) {
      fprintf(stderr, "Could not read line %lu of \"%s\"\n", i + 1, file);
      del_grid(sol->gr[number]);
      return -1;
    }


    args =
      sscanf(line, "%lf %lf %lf %lf %lf", &(sol->gr[number]->point[i][0]),
	     &(sol->gr[number]->point[i][1]), &(sol->gr[number]->point[i][2]),
	     &(sol->gr[number]->point[i][4]),
	     &(sol->gr[number]->point[i][6]));

    if (args != 5) {
      fprintf(stderr,
	      "Could not read five coordinates in line %lu of \"%s\"\n",
	      i + 1, file);
      del_grid(sol->gr[number]);
      return -1;
    }
  }

  /* Scale and find max, min  */
  for (j = 0; j < 3; j++) {
    min[j] = 0.0;
    max[j] = 0.0;
  }

  if (animate == true) {
    sol->timegr[number] = new_grid(max_x, max_y, 8);
    for (i = 0; i < l; i++) {
      /* copy original real and imaginary parts */
      sol->timegr[number]->point[i][0] = sol->gr[number]->point[i][2];
      sol->timegr[number]->point[i][1] = sol->gr[number]->point[i][4];

      sol->gr[number]->point[i][0] /= scale;
      sol->gr[number]->point[i][1] /= scale;
      for (j = 0; j < 3; j++) {
	min[j] =
	  (sol->gr[number]->point[i][2 * (j + 1)] <
	   min[j] ? sol->gr[number]->point[i][2 * (j + 1)] : min[j]);
	max[j] =
	  (sol->gr[number]->point[i][2 * (j + 1)] >
	   max[j] ? sol->gr[number]->point[i][2 * (j + 1)] : max[j]);
	sol->gr[number]->point[i][2 * (j + 1)] += z;
	sol->gr[number]->point[i][2 * (j + 1)] /= scale;
	sol->timegr[number]->point[i][2 * (j + 1)] =
	  sol->gr[number]->point[i][2 * (j + 1)];
      }
    }
    for (j = 0; j < 3; j++) {
      sol->min[number][j] = min[j];
      sol->max[number][j] = max[j];

      /* Evaluate color values */
      lenght = scale / (sol->max[number][j] - sol->min[number][j]);
      maxpoint = (sol->max[number][j] + z) / scale;

      for (i = 0; i < l; i++) {
	/* color */
	sol->gr[number]->point[i][2 * (j + 1) + 1] =
	  1.0 +
	  (lenght * 2.0 *
	   (sol->gr[number]->point[i][2 * (j + 1)] - maxpoint));
	sol->timegr[number]->point[i][2 * (j + 1) + 1] =
	  sol->gr[number]->point[i][2 * (j + 1) + 1];
	//printf("%f %f \n", sol->gr[number]->point[i][2], sol->gr[number]->point[i][3]);  
      }
    }
  }
  else {
    for (i = 0; i < l; i++) {
      sol->gr[number]->point[i][0] /= scale;
      sol->gr[number]->point[i][1] /= scale;
      for (j = 0; j < 3; j++) {
	min[j] =
	  (sol->gr[number]->point[i][2 * (j + 1)] <
	   min[j] ? sol->gr[number]->point[i][2 * (j + 1)] : min[j]);
	max[j] =
	  (sol->gr[number]->point[i][2 * (j + 1)] >
	   max[j] ? sol->gr[number]->point[i][2 * (j + 1)] : max[j]);
	sol->gr[number]->point[i][2 * (j + 1)] += z;
	sol->gr[number]->point[i][2 * (j + 1)] /= scale;
      }
    }
    for (j = 0; j < 3; j++) {
      sol->min[number][j] = min[j];
      sol->max[number][j] = max[j];

      /* Evaluate color values */
      lenght = scale / (sol->max[number][j] - sol->min[number][j]);
      maxpoint = (sol->max[number][j] + z) / scale;

      for (i = 0; i < l; i++) {
	/* color */
	sol->gr[number]->point[i][2 * (j + 1) + 1] =
	  1.0 +
	  (lenght * 2.0 *
	   (sol->gr[number]->point[i][2 * (j + 1)] - maxpoint));
	//printf("%f %f \n", sol->gr[number]->point[i][2], sol->gr[number]->point[i][3]);  
      }
    }

  }

  fclose(in);

  return 1;
}



   /*   
    *     Useful functions        *
    *    for working with         *
    *   dblocks and dcluster  *
    */

/*    static pdblock
    buildfromblock_dblock(pcblock b, pdcluster t, pdcluster s){

    pdblock nb;
    uint i, j;
    uint rsons, csons;
    rsons = b->rsons;
    csons = b->csons;
    nb = new_dblock(t , s, 0, 0, rsons, csons);
    nb->adm = b->a;

    if(rsons > 0 && csons > 0){
  	for(i = 0; i < csons; i++){
		for(j = 0; j < rsons; j++){
			if((s->sons==0) || (t->sons == 0)){
				if(s->sons == 0){
    					nb->son[j + i * rsons] = buildfromblock_dblock(b->son[j + i * rsons], t->son[j], s);		
				}
				else{
    				nb->son[j + i * rsons] = buildfromblock_dblock(b->son[j + i * rsons], t, s->son[i]);
				}
				if((s->sons==0) && (t->sons == 0)){
					printf("Something goes wrong!\n");
				}
			}
			else{
    				nb->son[j + i * rsons] = buildfromblock_dblock(b->son[j + i * rsons], t->son[j], s->son[i]);
			}
		}
  	}
	}
    update_dblock(nb);

    return nb;
    }
*/

static    real
getdiam_2_dcluster(pcdcluster dt)
{
  real      diam2;
  uint      i;

  diam2 = 0.0;
  for (i = 0; i < dt->dim; i++) {
    diam2 += REAL_SQR(dt->bmax[i] - dt->bmin[i]);
  }
  return REAL_SQRT(diam2);
}


static    real
normmax_avactor(pcavector v)
{

  real      norm;
  uint      i;

  norm = 0.0;
  for (i = 0; i < v->dim; i++) {
    norm = (ABS(REAL(v->v[i])) > norm ? ABS(REAL(v->v[i])) : norm);
  }

  return norm;
}

   /*
    *   Functions for   *
    *       movement    *
    */

pactivemovement
new_activemovement(bool level)
{

  pactivemovement move;

  move = (pactivemovement) allocmem(sizeof(activemovement));
  move->angle_x = 0.0;
  move->angle_y = 0.0;
  move->rightleft = 0.0;
  move->updown = 0.0;
  move->zoom = 1.0;
  move->levelup = level;

  return move;
}

porthomovement
new_orthomovement()
{

  porthomovement move;
  move = (porthomovement) allocmem(sizeof(orthomovement));
  move->xoff = 0.0;
  move->yoff = 0.0;

  return move;
}

ptablemovement
new_tablemovement()
{

  ptablemovement move;
  move = (ptablemovement) allocmem(sizeof(tablemovement));

  return move;
}

pgrid
new_grid(uint x, uint y, uint wide)
{

  longindex l, i;
  pgrid     p;
  p = (pgrid) allocmem(sizeof(grid));

  p->x = x;
  p->y = y;

  l = (longindex) x *y;
  p->point = (double **) allocmem(sizeof(real *) * l);	/* Entry for evey point */

  for (i = 0; i < l; i++) {
    p->point[i] = (double *) allocmem(sizeof(real) * wide);	/* Entry (x,y,z) */
  }

  return p;
}

psolution
new_solution(uint vectorsize, uint number, bool animated)
{

  psolution sol;
  uint      i, j;

  sol = (psolution) allocmem(sizeof(solution));
  if (vectorsize > 0) {
    sol->values = new_realavector(vectorsize);
  }
  else {
    sol->values = NULL;
  }

  sol->plane = true;
  sol->model = true;
  sol->animate = animated;

  if (number > 0) {
    sol->gr = (pgrid *) allocmem(sizeof(grid) * number);
    sol->min = (real **) allocmem(sizeof(real *) * number);
    sol->max = (real **) allocmem(sizeof(real *) * number);

    sol->scale = 1.0;
    for (i = 0; i < number; i++) {
      sol->min[i] = (real *) allocmem(sizeof(real) * 3);
      sol->max[i] = (real *) allocmem(sizeof(real) * 3);
      for (j = 0; j < 3; j++) {
	sol->min[i][j] = 0.0;
	sol->max[i][j] = 0.0;
      }

      if (animated == true) {
	sol->timegr = (pgrid *) allocmem(sizeof(grid) * number);
      }
      else {
	sol->timegr = NULL;
      }
    }
  }
  else {
    sol->gr = NULL;
    sol->min = (real **) allocmem(sizeof(real *) * 1);
    sol->max = (real **) allocmem(sizeof(real *) * 1);
    sol->min[0] = (real *) allocmem(sizeof(real) * 3);
    sol->max[0] = (real *) allocmem(sizeof(real) * 3);
    for (j = 0; j < 3; j++) {
      sol->min[0][j] = 0.0;
      sol->max[0][j] = 0.0;
    }
    sol->scale = 1.0;
    sol->timegr = NULL;
  }
  sol->current_grid_number = 0;
  sol->grids = number;
  sol->type = 0;

  /* For animation */
  sol->time = 0.0;
  sol->time_step = 1.0;		/* it could be usefull to change the time step */

  return sol;
}

void
del_activemovement(pactivemovement move)
{

  free(move);
}

void
del_orthomovement(porthomovement move)
{

  free(move);
}

void
del_tablemovement(ptablemovement move)
{

  free(move);
}

void
del_grid(pgrid p)
{

  longindex l, i;

  l = (longindex) p->x * p->y;

  for (i = 0; i < l; i++) {
    freemem(p->point[i]);
  }

  freemem(p->point);
  freemem(p);

}

void
del_solution(psolution sol)
{

  uint      i;

  if (sol->values) {
    del_realavector(sol->values);
  }
  if (sol->gr) {
    for (i = 0; i < sol->grids; i++) {
      del_grid(sol->gr[i]);
    }
    freemem(sol->gr);
  }
  if (sol->timegr) {
    for (i = 0; i < sol->grids; i++) {
      del_grid(sol->timegr[i]);
    }
    freemem(sol->timegr);
  }
  for (i = 0; i < sol->grids; i++) {
    freemem(sol->min[i]);
    freemem(sol->max[i]);
  }
  freemem(sol->max);
  freemem(sol->min);

  free(sol);
}

/*						
*	Drawing directional cluster bounding boxes 		*
								                                */

    /* Case: directional cluster is one dimensional */
static void
draw_1c(pcdcluster c, uint count, uint geometry)
{

  double    lenght;
  double    stage;
  uint      i;
  /* Evaluate scaling parameter */
  if (geometry == 0) {
    lenght = t->bmax[0] - t->bmin[0];
  }
  else {
    lenght = s->bmax[0] - s->bmin[0];
  }
  lenght *= 0.9;
  stage = 1.0 * (count - level);

  glPushMatrix();
  /* Translate different cluster level on different positions */
  glTranslated(0.0, (-0.5 + 1.0 * stage / (2.0 + level2 - level)) * lenght,
	       0.0);
  /* Draw cluster */
  glPointSize(3.0);
  glBegin(GL_POINTS);
  glVertex2d((c->bmin[0]), 0.0);
  glVertex2d((c->bmax[0]), 0.0);
  glEnd();

  glBegin(GL_LINES);
  glVertex2d((c->bmin[0]), 0.0);
  glVertex2d((c->bmax[0]), 0.0);
  glEnd();

  /* Draw points */
  if (points == 1) {
    glPushMatrix();
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPointSize(5.0);
    glColor3d(0.8, 0.8, 0.8);
    glBegin(GL_POINTS);
    for (i = 0; i < c->size; i++) {
      if (geometry == 0) {
	glVertex2d((gt->x[c->idx[i]][0]), 0.0);
      }
      else {
	glVertex2d((gs->x[c->idx[i]][0]), 0.0);
      }
    }
    glEnd();
    glPopMatrix();
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
  }
  glPopMatrix();

  if (bbox == true) {
    printf("Bounding box [%.6e, %.6e]\n", c->bmin[0], c->bmax[0]);
  }
}

  /* Case: cluster is one dimensional */
static void
draw_1cc(pccluster c, uint count, uint geometry)
{

  double    lenght;
  double    stage;
  uint      i;
  /* evaluate scaling parameter */
  if (geometry == 0) {
    lenght = ct->bmax[0] - ct->bmin[0];
  }
  else {
    lenght = cs->bmax[0] - cs->bmin[0];
  }
  lenght *= 0.9;
  stage = 1.0 * (count - level);

  glPushMatrix();
  /* Translate different cluster level on different positions */
  glTranslated(0.0, (-0.5 + 1.0 * stage / (2.0 + level2 - level)) * lenght,
	       0.0);
  /* Draw cluster */
  glPointSize(3.0);
  glBegin(GL_POINTS);
  glVertex2d((c->bmin[0]), 0.0);
  glVertex2d((c->bmax[0]), 0.0);
  glEnd();

  glBegin(GL_LINES);
  glVertex2d((c->bmin[0]), 0.0);
  glVertex2d((c->bmax[0]), 0.0);
  glEnd();

  /* Draw points */
  if (points == 1) {
    glPushMatrix();
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPointSize(5.0);
    glColor3d(0.8, 0.8, 0.8);
    glBegin(GL_POINTS);
    for (i = 0; i < c->size; i++) {
      if (geometry == 0) {
	glVertex2d((gt->x[c->idx[i]][0]), 0.0);
      }
      else {
	glVertex2d((gs->x[c->idx[i]][0]), 0.0);
      }
    }
    glEnd();
    glPopMatrix();
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
  }
  glPopMatrix();

  if (bbox == true) {
    printf("Bounding box [%.6e, %.6e]\n", c->bmin[0], c->bmax[0]);
  }
}

  /*Case: directional cluster is two dimensional */
static void
draw_2c(pcdcluster c, uint geometry)
{

  uint      i;
  /* Draw cluster */
  glBegin(GL_LINE_LOOP);
  glVertex2d((c->bmin[0]), (c->bmin[1]));
  glVertex2d((c->bmax[0]), (c->bmin[1]));
  glVertex2d((c->bmax[0]), (c->bmax[1]));
  glVertex2d((c->bmin[0]), (c->bmax[1]));
  glEnd();

  /* Draw points */
  if (points == 1) {
    glPushMatrix();
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor3d(0.8, 0.8, 0.8);
    glPointSize(5.0);
    glBegin(GL_POINTS);
    for (i = 0; i < c->size; i++) {
      if (geometry == 0) {
	glVertex2d((gt->x[c->idx[i]][0]), (gt->x[c->idx[i]][1]));
      }
      else {
	glVertex2d((gs->x[c->idx[i]][0]), (gs->x[c->idx[i]][1]));
      }
    }
    glEnd();
    glPopMatrix();
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
  }

  if (bbox == true) {
    printf("Bounding box [%.6e, %.6e] x [%.6e, %.6e]\n", c->bmin[0],
	   c->bmax[0], c->bmin[1], c->bmax[1]);
  }

}

  /*Case: cluster is two dimensional */
static void
draw_2cc(pccluster c, uint geometry)
{

  uint      i;
  /* Draw cluster */
  glBegin(GL_LINE_LOOP);
  glVertex2d((c->bmin[0]), (c->bmin[1]));
  glVertex2d((c->bmax[0]), (c->bmin[1]));
  glVertex2d((c->bmax[0]), (c->bmax[1]));
  glVertex2d((c->bmin[0]), (c->bmax[1]));
  glEnd();

  /* Draw points */
  if (points == 1) {
    glPushMatrix();
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor3d(0.8, 0.8, 0.8);
    glPointSize(5.0);
    glBegin(GL_POINTS);
    for (i = 0; i < c->size; i++) {
      if (geometry == 0) {
	glVertex2d((gt->x[c->idx[i]][0]), (gt->x[c->idx[i]][1]));
      }
      else {
	glVertex2d((gs->x[c->idx[i]][0]), (gs->x[c->idx[i]][1]));
      }
    }
    glEnd();
    glPopMatrix();
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
  }

  if (bbox == true) {
    printf("Bounding box [%.6e, %.6e] x [%.6e, %.6e]\n", c->bmin[0],
	   c->bmax[0], c->bmin[1], c->bmax[1]);
  }

}

  /*Case: directional cluster is three dimensional */
static void
draw_3c(pcdcluster c, uint geometry)
{

  uint      i;
  /* Draw cluster */
  /*top */
  glBegin(GL_LINE_LOOP);
  glVertex3d((c->bmin[0]), (c->bmax[1]), (c->bmin[2]));
  glVertex3d((c->bmax[0]), (c->bmax[1]), (c->bmin[2]));
  glVertex3d((c->bmax[0]), (c->bmax[1]), (c->bmax[2]));
  glVertex3d((c->bmin[0]), (c->bmax[1]), (c->bmax[2]));
  glEnd();

  /*down */
  glBegin(GL_LINE_LOOP);
  glVertex3d((c->bmin[0]), (c->bmin[1]), (c->bmin[2]));
  glVertex3d((c->bmax[0]), (c->bmin[1]), (c->bmin[2]));
  glVertex3d((c->bmax[0]), (c->bmin[1]), (c->bmax[2]));
  glVertex3d((c->bmin[0]), (c->bmin[1]), (c->bmax[2]));
  glEnd();

  /*back left and right, front left and right */
  glBegin(GL_LINES);
  glVertex3d((c->bmin[0]), (c->bmax[1]), (c->bmin[2]));
  glVertex3d((c->bmin[0]), (c->bmin[1]), (c->bmin[2]));

  glVertex3d((c->bmax[0]), (c->bmax[1]), (c->bmin[2]));
  glVertex3d((c->bmax[0]), (c->bmin[1]), (c->bmin[2]));

  glVertex3d((c->bmin[0]), (c->bmax[1]), (c->bmax[2]));
  glVertex3d((c->bmin[0]), (c->bmin[1]), (c->bmax[2]));

  glVertex3d((c->bmax[0]), (c->bmax[1]), (c->bmax[2]));
  glVertex3d((c->bmax[0]), (c->bmin[1]), (c->bmax[2]));
  glEnd();

  /* Draw points */
  if (points == 1) {
    glPushMatrix();
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPointSize(5.0);
    glColor3d(0.8, 0.8, 0.8);
    glBegin(GL_POINTS);
    for (i = 0; i < c->size; i++) {
      if (geometry == 0) {
	glVertex3d((gt->x[c->idx[i]][0]), (gt->x[c->idx[i]][1]),
		   (gt->x[c->idx[i]][2]));
      }
      else {
	glVertex3d((gs->x[c->idx[i]][0]), (gs->x[c->idx[i]][1]),
		   (gs->x[c->idx[i]][2]));
      }
    }
    glEnd();
    glPopMatrix();
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
  }

  if (bbox == true) {
    printf("Bounding box [%.6e, %.6e] x [%.6e, %.6e] x [%.6e, %.6e]\n",
	   c->bmin[0], c->bmax[0], c->bmin[1], c->bmax[1], c->bmin[2],
	   c->bmax[2]);
  }

}

  /*Case: cluster is three dimensional */
static void
draw_3cc(pccluster c, uint geometry)
{

  uint      i;

  /* Draw cluster */
  /*top */
  glBegin(GL_LINE_LOOP);
  glVertex3d((c->bmin[0]), (c->bmax[1]), (c->bmin[2]));
  glVertex3d((c->bmax[0]), (c->bmax[1]), (c->bmin[2]));
  glVertex3d((c->bmax[0]), (c->bmax[1]), (c->bmax[2]));
  glVertex3d((c->bmin[0]), (c->bmax[1]), (c->bmax[2]));
  glEnd();

  /*down */
  glBegin(GL_LINE_LOOP);
  glVertex3d((c->bmin[0]), (c->bmin[1]), (c->bmin[2]));
  glVertex3d((c->bmax[0]), (c->bmin[1]), (c->bmin[2]));
  glVertex3d((c->bmax[0]), (c->bmin[1]), (c->bmax[2]));
  glVertex3d((c->bmin[0]), (c->bmin[1]), (c->bmax[2]));
  glEnd();

  /*back left and right, front left and right */
  glBegin(GL_LINES);
  glVertex3d((c->bmin[0]), (c->bmax[1]), (c->bmin[2]));
  glVertex3d((c->bmin[0]), (c->bmin[1]), (c->bmin[2]));

  glVertex3d((c->bmax[0]), (c->bmax[1]), (c->bmin[2]));
  glVertex3d((c->bmax[0]), (c->bmin[1]), (c->bmin[2]));

  glVertex3d((c->bmin[0]), (c->bmax[1]), (c->bmax[2]));
  glVertex3d((c->bmin[0]), (c->bmin[1]), (c->bmax[2]));

  glVertex3d((c->bmax[0]), (c->bmax[1]), (c->bmax[2]));
  glVertex3d((c->bmax[0]), (c->bmin[1]), (c->bmax[2]));
  glEnd();

  /* Draw points */
  if (points == 1) {
    glPushMatrix();
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPointSize(5.0);
    glColor3d(0.8, 0.8, 0.8);
    glBegin(GL_POINTS);
    for (i = 0; i < c->size; i++) {
      if (geometry == 0) {
	glVertex3d((gt->x[c->idx[i]][0]), (gt->x[c->idx[i]][1]),
		   (gt->x[c->idx[i]][2]));
      }
      else {
	glVertex3d((gs->x[c->idx[i]][0]), (gs->x[c->idx[i]][1]),
		   (gs->x[c->idx[i]][2]));
      }
    }
    glEnd();
    glPopMatrix();
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
  }

  if (bbox == true) {
    printf("Bounding box [%.6e, %.6e] x [%.6e, %.6e] x [%.6e, %.6e]\n",
	   c->bmin[0], c->bmax[0], c->bmin[1], c->bmax[1], c->bmin[2],
	   c->bmax[2]);
  }

}

  /* Function for drawing directions */
static void
draw_directions(pcdcluster rc, pcdcluster cc, uint dim1, uint dim2)
{

  double    scale;
  uint      j;

  assert(dim1 < 4);
  assert(dim2 < 4);

  if (dim1 == 3 && rc != NULL) {

    glPushMatrix();
    scale = rc->bmax[0] - rc->bmin[0];
    scale =
      (rc->bmax[1] - rc->bmin[1] > scale ? rc->bmax[1] - rc->bmin[1] : scale);
    scale =
      (rc->bmax[2] - rc->bmin[2] > scale ? rc->bmax[2] - rc->bmin[2] : scale);
    glTranslated(rc->bmin[0], rc->bmin[1], rc->bmin[2]);
    glTranslated((rc->bmax[0] - rc->bmin[0]) * 0.5,
		 (rc->bmax[1] - rc->bmin[1]) * 0.5,
		 (rc->bmax[2] - rc->bmin[2]) * 0.5);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glLineWidth(3);
    glColor4d(0.7, 0.7, 0.7, 0.5);
    for (j = 0; j < rc->directions; j++) {
      glBegin(GL_LINES);
      glVertex3d(0.0, 0.0, 0.0);
      glVertex3d(rc->dir[j][0] * scale, rc->dir[j][1] * scale,
		 rc->dir[j][2] * scale);
      glEnd();
    }
    glDisable(GL_BLEND);
    glPopMatrix();
  }
  if (dim2 == 3 && cc != NULL) {

    glPushMatrix();
    scale = cc->bmax[0] - cc->bmin[0];
    scale =
      (cc->bmax[1] - cc->bmin[1] > scale ? cc->bmax[1] - cc->bmin[1] : scale);
    scale =
      (cc->bmax[2] - cc->bmin[2] > scale ? cc->bmax[2] - cc->bmin[2] : scale);
    glTranslated(cc->bmin[0], cc->bmin[1], cc->bmin[2]);
    glTranslated((cc->bmax[0] - cc->bmin[0]) * 0.5,
		 (cc->bmax[1] - cc->bmin[1]) * 0.5,
		 (cc->bmax[2] - cc->bmin[2]) * 0.5);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glLineWidth(3);
    glColor4d(0.7, 0.7, 0.7, 0.5);
    for (j = 0; j < cc->directions; j++) {
      glBegin(GL_LINES);
      glVertex3d(0.0, 0.0, 0.0);
      glVertex3d(cc->dir[j][0] * scale, cc->dir[j][1] * scale,
		 cc->dir[j][2] * scale);
      glEnd();
    }
    glDisable(GL_BLEND);
    glPopMatrix();
  }
  if (dim1 == 2 && rc != NULL) {

    glPushMatrix();
    scale = rc->bmax[0] - rc->bmin[0];
    scale =
      (rc->bmax[1] - rc->bmin[1] > scale ? rc->bmax[1] - rc->bmin[1] : scale);
    glTranslated(rc->bmin[0], rc->bmin[1], 0.0);
    glTranslated((rc->bmax[0] - rc->bmin[0]) * 0.5,
		 (rc->bmax[1] - rc->bmin[1]) * 0.5, 0.0);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glLineWidth(3);
    glColor4d(0.7, 0.7, 0.7, 0.5);
    for (j = 0; j < rc->directions; j++) {
      glBegin(GL_LINES);
      glVertex2d(0.0, 0.0);
      glVertex2d(rc->dir[j][0] * scale, rc->dir[j][1] * scale);
      glEnd();
    }
    glDisable(GL_BLEND);
    glPopMatrix();
  }
  if (dim2 == 2 && cc != NULL) {

    glPushMatrix();
    scale = cc->bmax[0] - cc->bmin[0];
    scale =
      (cc->bmax[1] - cc->bmin[1] > scale ? cc->bmax[1] - cc->bmin[1] : scale);
    glTranslated(cc->bmin[0], cc->bmin[1], 0.0);
    glTranslated((cc->bmax[0] - cc->bmin[0]) * 0.5,
		 (cc->bmax[1] - cc->bmin[1]) * 0.5, 0.0);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glLineWidth(3);
    glColor4d(0.7, 0.7, 0.7, 0.5);
    for (j = 0; j < cc->directions; j++) {
      glBegin(GL_LINES);
      glVertex2d(0.0, 0.0);
      glVertex2d(cc->dir[j][0] * scale, cc->dir[j][1] * scale);
      glEnd();
    }
    glDisable(GL_BLEND);
    glPopMatrix();

  }
}

   /*
    *   Funtions for triangulations       *
    *   and surfaces, outlines and      *
    *         and midpoints                 *
    */

static double
midpoint_surface(pcsurface3d sur)
{

  double    temp[3];
  uint      i;
  double    scale, tmp;

  temp[0] = sur->x[0][0];
  temp[1] = sur->x[0][1];
  temp[2] = sur->x[0][2];

  for (i = 1; i < gr->vertices; i++) {
    temp[0] += sur->x[i][0];
    temp[1] += sur->x[i][1];
    temp[2] += sur->x[i][2];
  }

  temp[0] /= gr->vertices;
  temp[1] /= gr->vertices;
  temp[2] /= gr->vertices;

  scale = 0.0;
  for (i = 0; i < gr->vertices; i++) {
    gr->x[i][0] = sur->x[i][0] - temp[0];
    gr->x[i][1] = sur->x[i][1] - temp[1];
    gr->x[i][2] = sur->x[i][2] - temp[2];

    tmp =
      (gr->x[i][0]) * (gr->x[i][0]) + (gr->x[i][1]) * (gr->x[i][1]) +
      (gr->x[i][2]) * (gr->x[i][2]);

    scale = (tmp > scale ? tmp : scale);
  }

  scale = sqrt(scale);

  if (scale <= 0.0) {
    printf("Something is wrong with the surface!\n");
    return EXIT_FAILURE;
  }

  for (i = 0; i < gr->vertices; i++) {
    gr->x[i][0] /= scale;
    gr->x[i][1] /= scale;
    gr->x[i][2] /= scale;
  }
  return scale;
}

static void
midpoint_triangle(pctri2d tri)
{

  double    temp[2];
  uint      i;
  double    scale, tmp;

  temp[0] = tri->x[0][0];
  temp[1] = tri->x[0][1];

  for (i = 1; i < tr->vertices; i++) {
    temp[0] += tri->x[i][0];
    temp[1] += tri->x[i][1];
  }

  temp[0] /= tr->vertices;
  temp[1] /= tr->vertices;

  scale = 0.0;
  for (i = 0; i < tr->vertices; i++) {
    tr->x[i][0] = tri->x[i][0] - temp[0];
    tr->x[i][1] = tri->x[i][1] - temp[1];

    tmp = (tr->x[i][0]) * (tr->x[i][0]) + (tr->x[i][1]) * (tr->x[i][1]);

    scale = (tmp > scale ? tmp : scale);
  }

  scale = sqrt(scale);

  if (scale <= 0.0) {
    printf("Something is wrong with the triangulation!\n");
    return;
  }

  for (i = 0; i < tr->vertices; i++) {
    tr->x[i][0] /= scale;
    tr->x[i][1] /= scale;
  }

}


static void
midpoint_tetrahedra(pctet3d tet)
{

  double    temp[3];
  uint      i;
  double    scale, tmp;

  temp[0] = tet->x[0][0];
  temp[1] = tet->x[0][1];
  temp[2] = tet->x[0][2];

  for (i = 1; i < th->vertices; i++) {
    temp[0] += tet->x[i][0];
    temp[1] += tet->x[i][1];
    temp[2] += tet->x[i][2];
  }

  temp[0] /= th->vertices;
  temp[1] /= th->vertices;
  temp[2] /= th->vertices;

  scale = 0.0;
  for (i = 0; i < th->vertices; i++) {
    th->x[i][0] = tet->x[i][0] - temp[0];
    th->x[i][1] = tet->x[i][1] - temp[1];
    th->x[i][2] = tet->x[i][2] - temp[2];

    tmp =
      (th->x[i][0]) * (th->x[i][0]) + (th->x[i][1]) * (th->x[i][1]) +
      (th->x[i][2]) * (th->x[i][2]);

    scale = (tmp > scale ? tmp : scale);
  }

  scale = sqrt(scale);

  if (scale <= 0.0) {
    printf("Something is wrong with the triangulation!\n");
    return;
  }

  for (i = 0; i < th->vertices; i++) {
    th->x[i][0] /= scale;
    th->x[i][1] /= scale;
    th->x[i][2] /= scale;
  }
}

static void
coordinate_color(bool c)
{
  if (c == false) {
    glColor4d(0.8, 0.0, 0.0, 0.4);
  }
  else {
    glColor4d(0.0, 0.8, 0.0, 0.4);
  }
}

  /* Coordinate system */
static void
coordinate(void)
{

  double    radius = 0.05;
  uint      i;

  glPushMatrix();
  glTranslated(-0.75, -0.75, 1.0);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  /* x */
  coordinate_color(light_coord[0]);
  glBegin(GL_LINES);
  glVertex3d(0.0, 0.0, 0.0);
  glVertex3d(0.25, 0.0, 0.0);
  glEnd();
  glBegin(GL_TRIANGLE_FAN);
  glVertex3d(0.25, 0.0, 0.0);
  for (i = 0; i <= 5; i++) {
    glVertex3d(0.2, radius * cos(2.0 * M_PI * i / 5),
	       radius * sin(2.0 * M_PI * i / 5));
    glVertex3d(0.2, radius * cos(2.0 * M_PI * i / 5),
	       radius * sin(2.0 * M_PI * i / 5));
  }
  glEnd();

  /* y */
  coordinate_color(light_coord[1]);
  glBegin(GL_LINES);
  glVertex3d(0.0, 0.0, 0.0);
  glVertex3d(0.0, 0.25, 0.0);
  glEnd();
  glBegin(GL_TRIANGLE_FAN);
  glVertex3d(0.0, 0.25, 0.0);
  for (i = 0; i <= 5; i++) {
    glVertex3f(radius * cos(2.0 * M_PI * i / 5), 0.2,
	       radius * sin(2.0 * M_PI * i / 5));
    glVertex3f(radius * cos(2.0 * M_PI * i / 5), 0.2,
	       radius * sin(2.0 * M_PI * i / 5));
  }
  glEnd();

  /* z */
  coordinate_color(light_coord[2]);
  glBegin(GL_LINES);
  glVertex3d(0.0, 0.0, 0.0);
  glVertex3d(0.0, 0.0, 0.25);
  glEnd();
  glBegin(GL_TRIANGLE_FAN);
  glVertex3d(0.0, 0.0, 0.25);
  for (i = 0; i <= 5; i++) {
    glVertex3f(radius * cos(2.0 * M_PI * i / 5),
	       radius * sin(2.0 * M_PI * i / 5), 0.2);
    glVertex3f(radius * cos(2.0 * M_PI * i / 5),
	       radius * sin(2.0 * M_PI * i / 5), 0.2);
  }
  glEnd();

  glDisable(GL_BLEND);
  glPopMatrix();
}


  /* Outlines of a given triangulation */
static void
outline(uint version)
{

  uint      j;

  switch (version) {
  case 0:			/* Triangle */
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4d(0.6, 0.6, 0.2, 0.7);
    glBegin(GL_LINES);
    for (j = 0; j < tr->edges; j++) {
      if (tr->eb[j] == 1) {
	glVertex3d(tr->x[tr->e[j][0]][0], tr->x[tr->e[j][0]][1], 0.0);
	glVertex3d(tr->x[tr->e[j][1]][0], tr->x[tr->e[j][1]][1], 0.0);
      }
    }
    glEnd();
    glDisable(GL_BLEND);
    break;

  case 1:			/* Tetrahedra */
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4d(0.6, 0.6, 0.2, 0.7);
    glBegin(GL_LINES);
    for (j = 0; j < th->edges; j++) {
      if (th->eb[j] == 1) {
	glVertex3d(th->x[th->e[j][0]][0], th->x[th->e[j][0]][1],
		   th->x[th->e[j][0]][2]);
	glVertex3d(th->x[th->e[j][1]][0], th->x[th->e[j][1]][1],
		   th->x[th->e[j][1]][2]);
      }
    }
    glEnd();
    glDisable(GL_BLEND);
    break;

  case 2:			/* Surface */
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4d(0.6, 0.6, 0.2, 0.7);
    glBegin(GL_LINES);
    for (j = 0; j < gr->edges; j++) {
      glVertex3d(gr->x[gr->e[j][0]][0], gr->x[gr->e[j][0]][1],
		 gr->x[gr->e[j][0]][2]);
      glVertex3d(gr->x[gr->e[j][1]][0], gr->x[gr->e[j][1]][1],
		 gr->x[gr->e[j][1]][2]);
    }
    glEnd();
    glDisable(GL_BLEND);
    break;

  default:
    printf("Sorry, there is something wrong!");
    break;
  }

}


  /*    
   *      Color settings        *
   */

  /* Set color depending on current directional cluster or block level */
static void
color_level_d(uint count)
{

  uint      colorsteps;
  /* Directional block version */
  if (version == 1) {
    colorsteps = getdepth_dblock(bl);

    if (colorsteps == 0) {
      if (t->dim == s->dim) {
	glColor3d(0.0, 0.8, 0.0);
      }
      else {
	glColor3d(0.8, 0.0, 0.0);
      }
    }
    else {
      if (t->dim == s->dim) {
	glColor3d(0.0, 1.0 - (1.0 * count / (colorsteps + 1)),
		  (1.0 * count / (colorsteps + 1)));
      }
      else {
	glColor3d(1.0 - (1.0 * count / (colorsteps + 1)), 0.0,
		  (1.0 * count / (colorsteps + 1)));
      }
    }
  }
  /* Directional cluster version */
  else {
    colorsteps = getdepth_dcluster(t);

    if (colorsteps == 0) {
      if (t->dim == 1) {
	glColor3d(0.0, 0.0, 0.8);
      }
      else {
	if (t->dim == 2) {
	  glColor3d(0.0, 0.8, 0.0);
	}
	if (t->dim == 3) {
	  glColor3d(0.8, 0.0, 0.0);
	}
      }
    }
    else {
      if (t->dim == 1) {
	glColor3d(0.0, (1.0 * count / (colorsteps + 1)),
		  1.0 - (1.0 * count / (colorsteps + 1)));
      }
      if (t->dim == 2) {
	glColor3d(0.0, 1.0 - (1.0 * count / (colorsteps + 1)),
		  (1.0 * count / (colorsteps + 1)));
      }
      if (t->dim == 3) {
	glColor3d(1.0 - (1.0 * count / (colorsteps + 1)),
		  (1.0 * count / (colorsteps + 1)), 0.0);
      }
    }
  }
}

  /* Set color depending on current cluster or block level */
static void
color_level(uint count)
{

  uint      colorsteps;
  /* Block version */
  if (version == 1) {
    colorsteps = getdepth_block(bbl);

    if (colorsteps == 0) {
      if (ct->dim == cs->dim) {
	glColor3d(0.0, 0.8, 0.0);
      }
      else {
	glColor3d(0.8, 0.0, 0.0);
      }
    }
    else {
      if (ct->dim == cs->dim) {
	glColor3d(0.0, 1.0 - (1.0 * count / (colorsteps + 1)),
		  (1.0 * count / (colorsteps + 1)));
      }
      else {
	glColor3d(1.0 - (1.0 * count / (colorsteps + 1)), 0.0,
		  (1.0 * count / (colorsteps + 1)));
      }
    }
  }
  /* Cluster version */
  else {
    colorsteps = getdepth_cluster(ct);

    if (colorsteps == 0) {
      if (ct->dim == 1) {
	glColor3d(0.0, 0.0, 0.8);
      }
      else {
	if (ct->dim == 2) {
	  glColor3d(0.0, 0.8, 0.0);
	}
	if (ct->dim == 3) {
	  glColor3d(0.8, 0.0, 0.0);
	}
      }
    }
    else {
      if (ct->dim == 1) {
	glColor3d(0.0, (1.0 * count / (colorsteps + 1)),
		  1.0 - (1.0 * count / (colorsteps + 1)));
      }
      if (ct->dim == 2) {
	glColor3d(0.0, 1.0 - (1.0 * count / (colorsteps + 1)),
		  (1.0 * count / (colorsteps + 1)));
      }
      if (ct->dim == 3) {
	glColor3d(1.0 - (1.0 * count / (colorsteps + 1)),
		  (1.0 * count / (colorsteps + 1)), 0.0);
      }
    }
  }
}

  /* Set color depending on admissible flag - directional version */
static void
color_admissible_d(pcdblock b)
{

  if (b->adm == 0) {
    glColor3d(0.8, 0.0, 0.0);
  }
  else {
    glColor3d(0.0, 0.8, 0.0);
  }
}

  /* Set color depending on admissible flag */
static void
color_admissible(pcblock b)
{

  if (b->a == 0) {
    glColor3d(0.8, 0.0, 0.0);
  }
  else {
    glColor3d(0.0, 0.8, 0.0);
  }
}

  /* Find color for solution value */
static void
color_triangle(real coef)
{

  if (coef < 0.0) {
    if (coef < -1.0) {
      coef = -1.0;
    }
    glColor3d(-coef, 0.0, 1.0 + coef);
  }
  else {
    if (coef > 1.0) {
      coef = 1.0;
    }
    glColor3d(coef, coef, 1.0 - coef);
  }

}

  /*
   *    Recursive functions for drawing         *
   *     clusters and blocks                    *
   */

  /* Draw orthogonal projection of a directional block tree */
static void
draw_dblock_ortho(pcdblock b, uint count, uint roffsize, uint coffsize,
		  uint mark)
{

  uint      rsons, csons;
  uint      roff, coff;
  uint      rsize, csize;
  uint      i, j;

  rsons = b->rsons;
  csons = b->csons;

  if ((rsons + csons) != 0) {
    if (count < level2) {
      count++;
      coff = 0;
      for (j = 0; j < csons; j++) {
	roff = 0;
	for (i = 0; i < rsons; i++) {
	  roffsize += roff;
	  coffsize += coff;
	  if (visu == 1) {
	    /* Mark the target block */
	    mark = 0;
	    if ((move_o->target_x > 1.0 * roffsize)
		&& (move_o->target_x <
		    1.0 * (roffsize + b->son[i + j * rsons]->rc->size))) {
	      mark = 1;
	    }

	    if ((move_o->target_y > 1.0 * coffsize)
		&& (move_o->target_y <
		    1.0 * (coffsize + b->son[i + j * rsons]->cc->size))) {
	      mark += 1;
	    }
	  }
	  draw_dblock_ortho(b->son[i + j * rsons], count, roffsize, coffsize,
			    mark);
	  roffsize -= roff;
	  coffsize -= coff;

	  roff += b->son[i + j * rsons]->rc->size;
	}
	assert(roff == b->rc->size);
	coff += b->son[j * rsons]->cc->size;
      }
      assert(coff == b->cc->size);
      count--;
    }
  }
  else {
    /* Drawing orthogonal projection from a leaf block */
    rsize = b->rc->size;
    csize = b->cc->size;

    glPushMatrix();

    glTranslated(roffsize * 1.0, coffsize * 1.0, 0.0);

    color_admissible_d(b);
    if (visu == 1) {
      if (mark == 2) {
	glColor3d(0.0, 0.0, 1.0);
      }
    }
    glBegin(GL_TRIANGLE_STRIP);
    glVertex2d(0.0, 0.0);
    glVertex2d(rsize * 1.0, 0.0);
    glVertex2d(0.0, csize * 1.0);
    glVertex2d(rsize * 1.0, csize * 1.0);

    glEnd();

    glColor3d(0.0, 0.0, 0.0);
    glBegin(GL_LINE_STRIP);
    glVertex2d(0.0, 0.0);
    glVertex2d(rsize * 1.0, 0.0);
    glVertex2d(rsize * 1.0, csize * 1.0);
    glVertex2d(0.0, csize * 1.0);
    glVertex2d(0.0, 0.0);
    glEnd();

    glPopMatrix();
    count--;
  }
}

  /* Draw orthogonal projection of a block tree */
static void
draw_block_ortho(pcblock b, uint count, uint roffsize, uint coffsize,
		 uint mark)
{

  uint      rsons, csons;
  uint      roff, coff;
  uint      rsize, csize;
  uint      i, j;

  rsons = b->rsons;
  csons = b->csons;

  if ((rsons + csons) != 0) {
    if (count < level2) {
      count++;
      coff = 0;
      for (j = 0; j < csons; j++) {
	roff = 0;
	for (i = 0; i < rsons; i++) {
	  roffsize += roff;
	  coffsize += coff;
	  if (visu == 1) {
	    /* Mark the target block */
	    mark = 0;
	    if ((move_o->target_x > 1.0 * roffsize)
		&& (move_o->target_x <
		    1.0 * (roffsize + b->son[i + j * rsons]->rc->size))) {
	      mark = 1;
	    }

	    if ((move_o->target_y > 1.0 * coffsize)
		&& (move_o->target_y <
		    1.0 * (coffsize + b->son[i + j * rsons]->cc->size))) {
	      mark += 1;
	    }
	  }
	  draw_block_ortho(b->son[i + j * rsons], count, roffsize, coffsize,
			   mark);
	  roffsize -= roff;
	  coffsize -= coff;

	  roff += b->son[i + j * rsons]->rc->size;
	}
	assert(roff == b->rc->size);
	coff += b->son[j * rsons]->cc->size;
      }
      assert(coff == b->cc->size);
      count--;
    }
  }
  else {
    /* Drawing orthogonal projection from a leaf block */
    rsize = b->rc->size;
    csize = b->cc->size;

    glPushMatrix();

    glTranslated(roffsize * 1.0, coffsize * 1.0, 0.0);

    color_admissible(b);
    if (visu == 1) {
      if (mark == 2) {
	glColor3d(0.0, 0.0, 1.0);
      }
    }
    glBegin(GL_TRIANGLE_STRIP);
    glVertex2d(0.0, 0.0);
    glVertex2d(rsize * 1.0, 0.0);
    glVertex2d(0.0, csize * 1.0);
    glVertex2d(rsize * 1.0, csize * 1.0);

    glEnd();

    glColor3d(0.0, 0.0, 0.0);
    glBegin(GL_LINE_STRIP);
    glVertex2d(0.0, 0.0);
    glVertex2d(rsize * 1.0, 0.0);
    glVertex2d(rsize * 1.0, csize * 1.0);
    glVertex2d(0.0, csize * 1.0);
    glVertex2d(0.0, 0.0);
    glEnd();

    glPopMatrix();
    count--;
  }
}

  /* Function for drawing the father of a given directional block */
static void
draw_certain_dblock_father(pcdblock b, uint mark)
{

  /* Finding the father */
  if (up == 2) {
    if (mark == 1) {
      return;
    }
    uint      rsons, csons;
    uint      i, j;

    rsons = b->rsons;
    csons = b->csons;

    if ((rsons + csons) != 0) {
      for (j = 0; j < csons; j++) {
	for (i = 0; i < rsons; i++) {
	  if (move_1->dblock == b->son[i + j * rsons]) {
	    move_1->dblock = b;

	    color_admissible_d(b);
	    if (b->rc->dim == 1) {
	      draw_1c(b->rc, level2 - level, 0);
	    }
	    if (b->rc->dim == 2) {
	      draw_2c(b->rc, 0);
	    }
	    if (b->rc->dim == 3) {
	      draw_3c(b->rc, 0);
	    }

	    color_admissible_d(b);
	    if (b->cc->dim == 1) {
	      draw_1c(b->cc, level2 - level, 1);
	    }
	    if (b->cc->dim == 2) {
	      draw_2c(b->cc, 1);
	    }
	    if (b->cc->dim == 3) {
	      draw_3c(b->cc, 1);
	    }
	    if (direction == true) {
	      draw_directions(b->rc, b->cc, b->rc->dim, b->cc->dim);
	      glLineWidth(1);
	    }
	    return;
	  }
	}
      }
      for (j = 0; j < csons; j++) {
	for (i = 0; i < rsons; i++) {
	  draw_certain_dblock_father(b->son[i + j * rsons], mark);
	}
      }
    }
  }
  else {
    /* Just repeat drawing old father */
    pcdblock  ob = move_1->dblock;
    color_admissible_d(ob);
    if (ob->rc->dim == 1) {
      draw_1c(ob->rc, level2 - level, 0);
    }
    if (ob->rc->dim == 2) {
      draw_2c(ob->rc, 0);
    }
    if (ob->rc->dim == 3) {
      draw_3c(ob->rc, 0);
    }

    color_admissible_d(ob);
    if (ob->cc->dim == 1) {
      draw_1c(ob->cc, level2 - level, 1);
    }
    if (ob->cc->dim == 2) {
      draw_2c(ob->cc, 1);
    }
    if (ob->cc->dim == 3) {
      draw_3c(ob->cc, 1);
    }
    if (direction == true) {
      draw_directions(ob->rc, ob->cc, ob->rc->dim, ob->cc->dim);
      glLineWidth(1);
    }
  }
}

  /* Function for drawing the father of a given block */
static void
draw_certain_block_father(pcblock b, uint mark)
{

  /* Finding the father */
  if (up == 2) {
    if (mark == 1) {
      return;
    }
    uint      rsons, csons;
    uint      i, j;

    rsons = b->rsons;
    csons = b->csons;

    if ((rsons + csons) != 0) {
      for (j = 0; j < csons; j++) {
	for (i = 0; i < rsons; i++) {
	  if (move_1->block == b->son[i + j * rsons]) {
	    move_1->block = b;

	    color_admissible(b);
	    if (b->rc->dim == 1) {
	      draw_1cc(b->rc, level2 - level, 0);
	    }
	    if (b->rc->dim == 2) {
	      draw_2cc(b->rc, 0);
	    }
	    if (b->rc->dim == 3) {
	      draw_3cc(b->rc, 0);
	    }

	    color_admissible(b);
	    if (b->cc->dim == 1) {
	      draw_1cc(b->cc, level2 - level, 1);
	    }
	    if (b->cc->dim == 2) {
	      draw_2cc(b->cc, 1);
	    }
	    if (b->cc->dim == 3) {
	      draw_3cc(b->cc, 1);
	    }
	    return;
	  }
	}
      }
      for (j = 0; j < csons; j++) {
	for (i = 0; i < rsons; i++) {
	  draw_certain_block_father(b->son[i + j * rsons], mark);
	}
      }
    }
  }
  else {
    /* Just repeat drawing old father */
    pcblock   ob = move_1->block;
    color_admissible(ob);
    if (ob->rc->dim == 1) {
      draw_1cc(ob->rc, level2 - level, 0);
    }
    if (ob->rc->dim == 2) {
      draw_2cc(ob->rc, 0);
    }
    if (ob->rc->dim == 3) {
      draw_3cc(ob->rc, 0);
    }

    color_admissible(ob);
    if (ob->cc->dim == 1) {
      draw_1cc(ob->cc, level2 - level, 1);
    }
    if (ob->cc->dim == 2) {
      draw_2cc(ob->cc, 1);
    }
    if (ob->cc->dim == 3) {
      draw_3cc(ob->cc, 1);
    }
  }
}

  /* Draw a choosen leaf directional block */
static pcdblock
draw_certain_dblock(pcdblock b, uint count, uint roffsize, uint coffsize,
		    uint mark)
{

  uint      rsons, csons;
  uint      roff, coff;
  uint      i, j;
  pcdblock  found = NULL;

  rsons = b->rsons;
  csons = b->csons;

  if ((csons + rsons) != 0) {
    if (count < level2) {
      count++;
      coff = 0;
      for (j = 0; j < csons; j++) {
	roff = 0;
	for (i = 0; i < rsons; i++) {
	  roffsize += roff;
	  coffsize += coff;
	  /* Mark the target leaf block */
	  mark = 0;
	  if ((move_o->target_x > 1.0 * roffsize)
	      && (move_o->target_x <
		  1.0 * (roffsize + b->son[i + j * rsons]->rc->size))) {
	    mark = 1;
	  }

	  if ((move_o->target_y > 1.0 * coffsize)
	      && (move_o->target_y <
		  1.0 * (coffsize + b->son[i + j * rsons]->cc->size))) {
	    mark += 1;
	  }

	  found =
	    draw_certain_dblock(b->son[i + j * rsons], count, roffsize,
				coffsize, mark);

	  if (mark == 2) {
	    return found;
	  }
	  roffsize -= roff;
	  coffsize -= coff;

	  roff += b->son[i + j * rsons]->rc->size;
	}
	assert(roff == b->rc->size);
	coff += b->son[j * rsons]->cc->size;
      }
      assert(coff == b->cc->size);
      count--;
    }
  }
  else {
    /* Draw the target */
    if (mark == 2) {
      color_admissible_d(b);
      if (b->rc->dim == 1) {
	draw_1c(b->rc, level2 - level, 0);
      }
      if (b->rc->dim == 2) {
	draw_2c(b->rc, 0);
      }
      if (b->rc->dim == 3) {
	draw_3c(b->rc, 0);
      }

      color_admissible_d(b);
      if (b->cc->dim == 1) {
	draw_1c(b->cc, level2 - level, 1);
      }
      if (b->cc->dim == 2) {
	draw_2c(b->cc, 1);
      }
      if (b->cc->dim == 3) {
	draw_3c(b->cc, 1);
      }
      if (direction == true) {
	draw_directions(b->rc, b->cc, b->rc->dim, b->cc->dim);
	glLineWidth(1);
      }

      found = b;
    }
  }
  return found;
}

  /* Draw a choosen leaf block */
static pcblock
draw_certain_block(pcblock b, uint count, uint roffsize, uint coffsize,
		   uint mark)
{

  uint      rsons, csons;
  uint      roff, coff;
  uint      i, j;
  pcblock   found = NULL;

  rsons = b->rsons;
  csons = b->csons;

  if ((csons + rsons) != 0) {
    if (count < level2) {
      count++;
      coff = 0;
      for (j = 0; j < csons; j++) {
	roff = 0;
	for (i = 0; i < rsons; i++) {
	  roffsize += roff;
	  coffsize += coff;
	  /* Mark the target leaf block */
	  mark = 0;
	  if ((move_o->target_x > 1.0 * roffsize)
	      && (move_o->target_x <
		  1.0 * (roffsize + b->son[i + j * rsons]->rc->size))) {
	    mark = 1;
	  }

	  if ((move_o->target_y > 1.0 * coffsize)
	      && (move_o->target_y <
		  1.0 * (coffsize + b->son[i + j * rsons]->cc->size))) {
	    mark += 1;
	  }

	  found =
	    draw_certain_block(b->son[i + j * rsons], count, roffsize,
			       coffsize, mark);

	  if (mark == 2) {
	    return found;
	  }
	  roffsize -= roff;
	  coffsize -= coff;

	  roff += b->son[i + j * rsons]->rc->size;
	}
	assert(roff == b->rc->size);
	coff += b->son[j * rsons]->cc->size;
      }
      assert(coff == b->cc->size);
      count--;
    }
  }
  else {
    /* Draw the target */
    if (mark == 2) {
      color_admissible(b);
      if (b->rc->dim == 1) {
	draw_1cc(b->rc, level2 - level, 0);
      }
      if (b->rc->dim == 2) {
	draw_2cc(b->rc, 0);
      }
      if (b->rc->dim == 3) {
	draw_3cc(b->rc, 0);
      }

      color_admissible(b);
      if (b->cc->dim == 1) {
	draw_1cc(b->cc, level2 - level, 1);
      }
      if (b->cc->dim == 2) {
	draw_2cc(b->cc, 1);
      }
      if (b->cc->dim == 3) {
	draw_3cc(b->cc, 1);
      }

      found = b;
    }
  }
  return found;
}

  /* Drawing choosen level of a given directional cluster tree */
static void
draw_dcluster_level(pcdcluster c, uint count, uint geometry)
{

  uint      i;

  if (count < level) {
    count++;
    for (i = 0; i < c->sons; i++) {
      draw_dcluster_level(c->son[i], count, geometry);
    }
    count--;
  }
  else {
    if (count <= level2) {
      count++;
      if (t->dim == 1) {
	glColor3d(0.0, 0.8, 0.0);
	if (color == true) {
	  color_level_d(count);
	}
	draw_1c(c, count, geometry);
      }
      if (t->dim == 2) {
	glColor3d(0.0, 0.0, 0.8);
	if (color == true) {
	  color_level_d(count);
	}
	draw_2c(c, geometry);
      }
      if (t->dim == 3) {
	glColor3d(0.8, 0.0, 0.0);
	if (color == true) {
	  color_level_d(count);
	}
	draw_3c(c, geometry);
      }
      if (direction == true) {
	draw_directions(c, NULL, c->dim, 0);
	glLineWidth(1);
      }
      for (i = 0; i < c->sons; i++) {
	draw_dcluster_level(c->son[i], count, geometry);
      }
      count--;
    }
  }
}

  /* Drawing choosen level of a given cluster tree */
static void
draw_cluster_level(pccluster c, uint count, uint geometry)
{

  uint      i;

  if (count < level) {
    count++;
    for (i = 0; i < c->sons; i++) {
      draw_cluster_level(c->son[i], count, geometry);
    }
    count--;
  }
  else {
    if (count <= level2) {
      count++;
      if (ct->dim == 1) {
	glColor3d(0.0, 0.8, 0.0);
	if (color == true) {
	  color_level(count);
	}
	draw_1cc(c, count, geometry);
      }
      if (ct->dim == 2) {
	glColor3d(0.0, 0.0, 0.8);
	if (color == true) {
	  color_level(count);
	}
	draw_2cc(c, geometry);
      }
      if (ct->dim == 3) {
	glColor3d(0.8, 0.0, 0.0);
	if (color == true) {
	  color_level(count);
	}
	draw_3cc(c, geometry);
      }

      for (i = 0; i < c->sons; i++) {
	draw_cluster_level(c->son[i], count, geometry);
      }
      count--;
    }
  }
}

  /* Drawing choosen level of a directional block tree */
static void
draw_two_dcluster(pcdblock b, uint count)
{

  uint      rsons, csons;
  uint      roff, coff;
  uint      i, j;
  uint      w;
  pcdcluster q;
  uint      geometry;

  rsons = b->rsons;
  csons = b->csons;

  if (visu == 1) {
    if (count < level) {
      if ((rsons + csons) != 0) {
	count++;
	coff = 0;
	for (j = 0; j < csons; j++) {
	  roff = 0;
	  for (i = 0; i < rsons; i++) {
	    draw_two_dcluster(b->son[i + j * rsons], count);
	    roff += b->son[i + j * rsons]->rc->size;
	  }
	  assert(roff == b->rc->size);
	  coff += b->son[j * rsons]->cc->size;
	}
	assert(coff == b->cc->size);
	count--;
      }
    }
    else {
      if (count <= level2) {
	w = glutGetWindow();
	count++;
	if (w == clust2) {
	  q = b->cc;
	  geometry = 1;
	}
	else {
	  q = b->rc;
	  geometry = 0;
	}
	if (q->dim == 1) {
	  glColor3d(0.0, 0.8, 0.0);
	  if (color == true) {
	    color_level_d(count);
	  }
	  draw_1c(q, count, geometry);
	}
	if (q->dim == 2) {
	  glColor3d(0.0, 0.0, 0.8);
	  if (color == true) {
	    color_level_d(count);
	  }
	  draw_2c(q, geometry);
	}
	if (q->dim == 3) {
	  glColor3d(0.8, 0.0, 0.0);
	  if (color == true) {
	    color_level_d(count);
	  }
	  draw_3c(q, geometry);
	}
	if (direction == true) {
	  draw_directions(q, NULL, q->dim, 0);
	  glLineWidth(1);
	}

	if ((rsons + csons) != 0) {
	  coff = 0;
	  for (j = 0; j < csons; j++) {
	    roff = 0;
	    for (i = 0; i < rsons; i++) {
	      draw_two_dcluster(b->son[i + j * rsons], count);
	      roff += b->son[i + j * rsons]->rc->size;
	    }
	    assert(roff == b->rc->size);
	    coff += b->son[j * rsons]->cc->size;
	  }
	  assert(coff == b->cc->size);
	  count--;
	}
      }
    }
  }
}

  /* Drawing choosen level of a block tree */
static void
draw_two_cluster(pcblock b, uint count)
{

  uint      rsons, csons;
  uint      roff, coff;
  uint      i, j;
  uint      w;
  pccluster q;
  uint      geometry;

  rsons = b->rsons;
  csons = b->csons;

  if (visu == 1) {

    if (count < level) {
      if ((rsons + csons) != 0) {
	count++;
	coff = 0;
	for (j = 0; j < csons; j++) {
	  roff = 0;
	  for (i = 0; i < rsons; i++) {
	    draw_two_cluster(b->son[i + j * rsons], count);
	    roff += b->son[i + j * rsons]->rc->size;
	  }
	  assert(roff == b->rc->size);
	  coff += b->son[j * rsons]->cc->size;
	}
	assert(coff == b->cc->size);
	count--;
      }
    }
    else {
      if (count <= level2) {
	w = glutGetWindow();
	count++;

	if (w == clust2) {
	  q = b->cc;
	  geometry = 1;
	}
	else {
	  q = b->rc;
	  geometry = 0;
	}
	if (q->dim == 1) {
	  glColor3d(0.0, 0.8, 0.0);
	  if (color == true) {
	    color_level(count);
	  }
	  draw_1cc(q, count, geometry);
	}
	if (q->dim == 2) {
	  glColor3d(0.0, 0.0, 0.8);
	  if (color == true) {
	    color_level(count);
	  }
	  draw_2cc(q, geometry);
	}
	if (q->dim == 3) {
	  glColor3d(0.8, 0.0, 0.0);
	  if (color == true) {
	    color_level(count);
	  }
	  draw_3cc(q, geometry);
	}
	if ((rsons + csons) != 0) {
	  coff = 0;
	  for (j = 0; j < csons; j++) {
	    roff = 0;
	    for (i = 0; i < rsons; i++) {
	      draw_two_cluster(b->son[i + j * rsons], count);
	      roff += b->son[i + j * rsons]->rc->size;
	    }
	    assert(roff == b->rc->size);
	    coff += b->son[j * rsons]->cc->size;
	  }
	  assert(coff == b->cc->size);
	  count--;
	}
      }
    }
  }
}

  /*                                                            
   *    GLUT display functions  *
   */

  /* Drawing orthogonal projection of a given directional block tree */
static void
display_ortho_d()
{

  uint      count;
  uint      roffsize, coffsize;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if (bl->rc->size > bl->cc->size) {
    move_o->clustscaling_x = 1.0 / ((bl->rc->size) * 0.5);
    move_o->clustscaling_y = 1.0 / ((bl->rc->size) * 0.5);
    move_o->yoff = (bl->rc->size - bl->cc->size) * 0.5;
  }
  else {
    move_o->clustscaling_y = 1.0 / ((bl->cc->size) * 0.5);
    move_o->clustscaling_x = 1.0 / ((bl->cc->size) * 0.5);
    move_o->xoff = (bl->cc->size - bl->rc->size) * 0.5;
  }

  glPushMatrix();
  system_transformation();
  glRotated(180, 1.0, 0.0, 0.0);
  glScaled(move_o->clustscaling_x, move_o->clustscaling_y, 1.0);
  glTranslated(-(0.5 * (bl->rc->size)), -(0.5 * (bl->cc->size)), 0.0);

  count = 0;
  roffsize = 0;
  coffsize = 0;
  draw_dblock_ortho(bl, count, roffsize, coffsize, 0);

  glPopMatrix();

  glFlush();
  glutSwapBuffers();

  if (take_a_shot[0] == true) {

    static GLubyte *pixels;	/* Pointer for pixels */
    uint      paper_width = glutGet(GLUT_WINDOW_WIDTH);
    uint      paper_height = glutGet(GLUT_WINDOW_HEIGHT);

    pixels = allocmem(format_bytes * paper_height * paper_width);	/* Storage for the pixel */
    /* Reads the current buffer and save the pixels. RGBA has unsigned int type. */
    glReadPixels(0, 0, paper_width, paper_height, format, GL_UNSIGNED_BYTE,
		 pixels);

#ifdef USE_LIBPNG
    png_byte *png_bytes = NULL;	/* Memory for the png */
    png_byte **png_rows = NULL;	/* Array of memory for ever row */

    save_current_content_png("snapshot", paper_width, paper_height, pixels,
			     &png_bytes, &png_rows);
    free(png_bytes);
    free(png_rows);
#else
    save_current_content_ppm("snapshot", paper_width, paper_height, 255,
			     pixels);
#endif
    snapshots++;

    free(pixels);		/* free memory */
    take_a_shot[0] = false;	/* reset start values */
  }

}

    /* Drawing orthogonal projection of a given block tree */
static void
display_ortho()
{

  uint      count;
  uint      roffsize, coffsize;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if (bbl->rc->size > bbl->cc->size) {
    move_o->clustscaling_x = 1.0 / ((bbl->rc->size) * 0.5);
    move_o->clustscaling_y = 1.0 / ((bbl->rc->size) * 0.5);
    move_o->yoff = (bbl->rc->size - bbl->cc->size) * 0.5;
  }
  else {
    move_o->clustscaling_y = 1.0 / ((bbl->cc->size) * 0.5);
    move_o->clustscaling_x = 1.0 / ((bbl->cc->size) * 0.5);
    move_o->xoff = (bbl->cc->size - bbl->rc->size) * 0.5;
  }

  glPushMatrix();
  system_transformation();
  glRotated(180, 1.0, 0.0, 0.0);
  glScaled(move_o->clustscaling_x, move_o->clustscaling_y, 1.0);
  glTranslated(-(0.5 * (bbl->rc->size)), -(0.5 * (bbl->cc->size)), 0.0);

  count = 0;
  roffsize = 0;
  coffsize = 0;
  draw_block_ortho(bbl, count, roffsize, coffsize, 0);

  glPopMatrix();

  glFlush();
  glutSwapBuffers();

  if (take_a_shot[0] == true) {

    static GLubyte *pixels;	/* Pointer for pixels */
    uint      paper_width = glutGet(GLUT_WINDOW_WIDTH);
    uint      paper_height = glutGet(GLUT_WINDOW_HEIGHT);

    pixels = allocmem(format_bytes * paper_height * paper_width);	/* Storage for the pixel */
    /* Reads the current buffer and save the pixels. RGBA has unsigned int type. */
    glReadPixels(0, 0, paper_width, paper_height, format, GL_UNSIGNED_BYTE,
		 pixels);

#ifdef USE_LIBPNG
    png_byte *png_bytes = NULL;	/* Memory for the png */
    png_byte **png_rows = NULL;	/* Array of memory for ever row */

    save_current_content_png("snapshot", paper_width, paper_height, pixels,
			     &png_bytes, &png_rows);
    free(png_bytes);
    free(png_rows);
#else
    save_current_content_ppm("snapshot", paper_width, paper_height, 255,
			     pixels);
#endif
    snapshots++;

    free(pixels);		/* free memory */
    take_a_shot[0] = false;	/* reset start values */
  }

}

    /* Drawing level structure of the given directional block or cluster tree */
static void
display_dblock()
{

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  uint      step;

  if (version == 1) {
    step = getdepth_dblock(bl);
  }
  else {
    step = getdepth_dcluster(t);
  }
  uint      i, j;
  char      stroke[1][6];
  double    scale = 0.75 / (step + 1);
  double    eps = (2.0 / (step + 1) - scale) / 2.0;

  system_transformation();
  /* Drawing colored rectangle and write level number in it */
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslated(-1.0, 1.0, 0.0);

  glBegin(GL_TRIANGLE_STRIP);
  glColor3d(0.0, 1.0, 0.0);
  glVertex2d(0.0, 0.0);
  glVertex2d(2.0, 0.0);
  glVertex2d(0.0, -2.0 / (step + 1));
  glVertex2d(2.0, -2.0 / (step + 1));
  glEnd();

  glPushMatrix();

  glTranslatef(1.0 - 0.5 * eps, -2.0 / (step + 1.0) + eps, 0.0);
  glScalef(scale * 1.0 / 152.0, scale * 1.0 / 152.0, scale * 1.0 / 152.0);
  sprintf(stroke[0], "0");
  glColor3d(0.0, 0.0, 0.0);
  glLineWidth(2);

  for (j = 0; j < (uint) strlen(stroke[0]); j++) {
    glutStrokeCharacter(GLUT_STROKE_ROMAN, stroke[0][j]);
  }

  glPopMatrix();

  for (i = 1; i <= step; i++) {
    glBegin(GL_TRIANGLE_STRIP);
    glColor3d(0.0, 1.0 - 1.0 * i / (step + 1), 0.0 + 1.0 * i / (step + 1));
    glVertex2d(0.0, -2.0 * i / (step + 1));
    glVertex2d(1.0, -2.0 * i / (step + 1));
    glVertex2d(0.0, -2.0 * (i + 1) / (step + 1));
    glVertex2d(1.0, -2.0 * (i + 1) / (step + 1));
    glEnd();

    glPushMatrix();

    glTranslatef(0.0 + eps, -2.0 * (i + 1.0) / (step + 1.0) + eps, 0.0);
    glScalef(scale * 1.0 / 152.0, scale * 1.0 / 152.0, scale * 1.0 / 152.0);
    sprintf(stroke[0], "%d", i);
    glColor3d(0.0, 0.0, 0.0);
    glLineWidth(2);

    for (j = 0; j < (int) strlen(stroke[0]); j++) {
      glutStrokeCharacter(GLUT_STROKE_ROMAN, stroke[0][j]);
    }

    glPopMatrix();

    glBegin(GL_TRIANGLE_STRIP);
    glColor3d(0.0, 1.0 - 1.0 * i / (step + 1), 0.0 + 1.0 * i / (step + 1));
    glVertex2d(1.0, -2.0 * i / (step + 1));
    glVertex2d(2.0, -2.0 * i / (step + 1));
    glVertex2d(1.0, -2.0 * (i + 1) / (step + 1));
    glVertex2d(2.0, -2.0 * (i + 1) / (step + 1));
    glEnd();

    glPushMatrix();
    glTranslatef(1.0 + eps, -2.0 * (i + 1.0) / (step + 1.0) + eps, 0.0);
    glScalef(scale * 1.0 / 152.0, scale * 1.0 / 152.0, scale * 1.0 / 152.0);
    sprintf(stroke[0], "0-%d", i);
    glColor3d(0.0, 0.0, 0.0);

    for (j = 0; j < (int) strlen(stroke[0]); j++) {
      glutStrokeCharacter(GLUT_STROKE_ROMAN, stroke[0][j]);
    }

    glPopMatrix();
  }

  glLineWidth(1);

  glColor3d(0.0, 0.0, 0.0);
  glBegin(GL_LINES);
  glVertex2d(1.0, -2.0 / (step + 1));
  glVertex2d(1.0, -2.0);
  glEnd();

  for (i = 0; i < step; i++) {
    glBegin(GL_LINES);
    glVertex2d(0.0, -2.0 * (i + 1) / (step + 1));
    glVertex2d(2.0, -2.0 * (i + 1) / (step + 1));
    glEnd();
  }

  glFlush();
  glutSwapBuffers();

  if (take_a_shot[0] == true) {

    static GLubyte *pixels;	/* Pointer for pixels */
    uint      paper_width = glutGet(GLUT_WINDOW_WIDTH);
    uint      paper_height = glutGet(GLUT_WINDOW_HEIGHT);

    pixels = allocmem(format_bytes * paper_height * paper_width);	/* Storage for the pixel */
    /* Reads the current buffer and save the pixels. RGBA has unsigned int type. */
    glReadPixels(0, 0, paper_width, paper_height, format, GL_UNSIGNED_BYTE,
		 pixels);

#ifdef USE_LIBPNG
    png_byte *png_bytes = NULL;	/* Memory for the png */
    png_byte **png_rows = NULL;	/* Array of memory for ever row */

    save_current_content_png("snapshot", paper_width, paper_height, pixels,
			     &png_bytes, &png_rows);
    free(png_bytes);
    free(png_rows);
#else
    save_current_content_ppm("snapshot", paper_width, paper_height, 255,
			     pixels);
#endif
    snapshots++;

    free(pixels);		/* free memory */
    take_a_shot[0] = false;	/* reset start values */
  }
}

  /* Drawing level structure of the given block or cluster tree */
static void
display_block()
{

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  uint      step;

  if (version == 1) {
    step = getdepth_block(bbl);
  }
  else {
    step = getdepth_cluster(ct);
  }
  uint      i, j;
  char      stroke[1][6];
  double    scale = 0.75 / (step + 1);
  double    eps = (2.0 / (step + 1) - scale) / 2.0;

  system_transformation();
  /* Drawing colored rectangle and write level number in it */
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslated(-1.0, 1.0, 0.0);

  glBegin(GL_TRIANGLE_STRIP);
  glColor3d(0.0, 1.0, 0.0);
  glVertex2d(0.0, 0.0);
  glVertex2d(2.0, 0.0);
  glVertex2d(0.0, -2.0 / (step + 1));
  glVertex2d(2.0, -2.0 / (step + 1));
  glEnd();

  glPushMatrix();

  glTranslatef(1.0 - 0.5 * eps, -2.0 / (step + 1.0) + eps, 0.0);
  glScalef(scale * 1.0 / 152.0, scale * 1.0 / 152.0, scale * 1.0 / 152.0);
  sprintf(stroke[0], "0");
  glColor3d(0.0, 0.0, 0.0);
  glLineWidth(2);

  for (j = 0; j < (uint) strlen(stroke[0]); j++) {
    glutStrokeCharacter(GLUT_STROKE_ROMAN, stroke[0][j]);
  }

  glPopMatrix();

  for (i = 1; i <= step; i++) {
    glBegin(GL_TRIANGLE_STRIP);
    glColor3d(0.0, 1.0 - 1.0 * i / (step + 1), 0.0 + 1.0 * i / (step + 1));
    glVertex2d(0.0, -2.0 * i / (step + 1));
    glVertex2d(1.0, -2.0 * i / (step + 1));
    glVertex2d(0.0, -2.0 * (i + 1) / (step + 1));
    glVertex2d(1.0, -2.0 * (i + 1) / (step + 1));
    glEnd();

    glPushMatrix();

    glTranslatef(0.0 + eps, -2.0 * (i + 1.0) / (step + 1.0) + eps, 0.0);
    glScalef(scale * 1.0 / 152.0, scale * 1.0 / 152.0, scale * 1.0 / 152.0);
    sprintf(stroke[0], "%d", i);
    glColor3d(0.0, 0.0, 0.0);
    glLineWidth(2);

    for (j = 0; j < (int) strlen(stroke[0]); j++) {
      glutStrokeCharacter(GLUT_STROKE_ROMAN, stroke[0][j]);
    }

    glPopMatrix();

    glBegin(GL_TRIANGLE_STRIP);
    glColor3d(0.0, 1.0 - 1.0 * i / (step + 1), 0.0 + 1.0 * i / (step + 1));
    glVertex2d(1.0, -2.0 * i / (step + 1));
    glVertex2d(2.0, -2.0 * i / (step + 1));
    glVertex2d(1.0, -2.0 * (i + 1) / (step + 1));
    glVertex2d(2.0, -2.0 * (i + 1) / (step + 1));
    glEnd();

    glPushMatrix();
    glTranslatef(1.0 + eps, -2.0 * (i + 1.0) / (step + 1.0) + eps, 0.0);
    glScalef(scale * 1.0 / 152.0, scale * 1.0 / 152.0, scale * 1.0 / 152.0);
    sprintf(stroke[0], "0-%d", i);
    glColor3d(0.0, 0.0, 0.0);

    for (j = 0; j < (int) strlen(stroke[0]); j++) {
      glutStrokeCharacter(GLUT_STROKE_ROMAN, stroke[0][j]);
    }

    glPopMatrix();
  }

  glLineWidth(1);

  glColor3d(0.0, 0.0, 0.0);
  glBegin(GL_LINES);
  glVertex2d(1.0, -2.0 / (step + 1));
  glVertex2d(1.0, -2.0);
  glEnd();

  for (i = 0; i < step; i++) {
    glBegin(GL_LINES);
    glVertex2d(0.0, -2.0 * (i + 1) / (step + 1));
    glVertex2d(2.0, -2.0 * (i + 1) / (step + 1));
    glEnd();
  }

  glFlush();
  glutSwapBuffers();

  if (take_a_shot[0] == true) {

    static GLubyte *pixels;	/* Pointer for pixels */
    uint      paper_width = glutGet(GLUT_WINDOW_WIDTH);
    uint      paper_height = glutGet(GLUT_WINDOW_HEIGHT);

    pixels = allocmem(format_bytes * paper_height * paper_width);	/* Storage for the pixel */
    /* Reads the current buffer and save the pixels. RGBA has unsigned int type. */
    glReadPixels(0, 0, paper_width, paper_height, format, GL_UNSIGNED_BYTE,
		 pixels);

#ifdef USE_LIBPNG
    png_byte *png_bytes = NULL;	/* Memory for the png */
    png_byte **png_rows = NULL;	/* Array of memory for ever row */

    save_current_content_png("snapshot", paper_width, paper_height, pixels,
			     &png_bytes, &png_rows);
    free(png_bytes);
    free(png_rows);
#else
    save_current_content_ppm("snapshot", paper_width, paper_height, 255,
			     pixels);
#endif
    snapshots++;

    free(pixels);		/* free memory */
    take_a_shot[0] = false;	/* reset start values */
  }
}

  /* Drawing both cluster of a given directional block */
static void
display_two_dcluster()
{

  uint      count;
  double    scale, tmp;
  uint      w;
  w = glutGetWindow();
  pcdcluster q;
  pactivemovement move;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if (bbox == true) {
    bbox = false;
  }

  if (visu == 1) {

    if (w == clust2) {
      q = s;
      move = move_2;
    }
    else {
      q = t;
      move = move_1;
    }
    scale = (q->bmax[0] - q->bmin[0]) * 0.5;
    if (q->dim >= 2) {
      tmp = (q->bmax[1] - q->bmin[1]) * 0.5;
      if (tmp > scale) {
	scale = tmp;
      }
    }
    if (q->dim == 3) {
      tmp = (q->bmax[2] - q->bmin[2]) * 0.5;
      if (tmp > scale) {
	scale = tmp;
      }
    }
    scale = 0.95 / scale;

    glPushMatrix();
    system_transformation();
    glScaled(move->zoom * scale, move->zoom * scale, move->zoom * scale);
    glTranslated(move->rightleft, move->updown, 0.0);

    /* Put bounding box midpoint in (0,0,0) */
    if (q->dim == 1) {
      glTranslated(-(q->bmax[0] + q->bmin[0]) * 0.5, 0.0, 0.0);
    }
    else {
      if (q->dim == 2) {
	glTranslated(-(q->bmax[0] + q->bmin[0]) * 0.5,
		     -(q->bmax[1] + q->bmin[1]) * 0.5, 0.0);
      }
      else {
	glRotated(move->angle_x, 1.0, 0.0, 0.0);
	glRotated(move->angle_y, 0.0, 1.0, 0.0);

	glTranslated(-(q->bmax[0] + q->bmin[0]) * 0.5,
		     -(q->bmax[1] + q->bmin[1]) * 0.5,
		     -(q->bmax[2] + q->bmin[2]) * 0.5);
      }
    }
    count = 0;
    draw_two_dcluster(bl, count);

    glPopMatrix();

  }

  glFlush();
  glutSwapBuffers();

  if (take_a_shot[1] == true) {

    static GLubyte *pixels;	/* Pointer for pixels */
    uint      paper_width = glutGet(GLUT_WINDOW_WIDTH);
    uint      paper_height = glutGet(GLUT_WINDOW_HEIGHT);

    pixels = allocmem(format_bytes * paper_height * paper_width);	/* Storage for the pixel */
    /* Reads the current buffer and save the pixels. RGBA has unsigned int type. */
    glReadPixels(0, 0, paper_width, paper_height, format, GL_UNSIGNED_BYTE,
		 pixels);

#ifdef USE_LIBPNG
    png_byte *png_bytes = NULL;	/* Memory for the png */
    png_byte **png_rows = NULL;	/* Array of memory for ever row */

    save_current_content_png("snapshot", paper_width, paper_height, pixels,
			     &png_bytes, &png_rows);
    free(png_bytes);
    free(png_rows);
#else
    save_current_content_ppm("snapshot", paper_width, paper_height, 255,
			     pixels);
#endif
    snapshots++;

    free(pixels);		/* free memory */
    take_a_shot[1] = false;	/* reset start values */
  }

}

  /* Drawing both cluster of a given block */
static void
display_two_cluster()
{

  uint      count;
  double    scale, tmp;
  uint      w;
  w = glutGetWindow();
  pccluster q;
  pactivemovement move;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if (bbox == true) {
    bbox = false;
  }

  if (visu == 1) {
    if (w == clust2) {
      q = cs;
      move = move_2;
    }
    else {
      q = ct;
      move = move_1;
    }
    scale = (q->bmax[0] - q->bmin[0]) * 0.5;
    if (q->dim >= 2) {
      tmp = (q->bmax[1] - q->bmin[1]) * 0.5;
      if (tmp > scale) {
	scale = tmp;
      }
    }
    if (q->dim == 3) {
      tmp = (q->bmax[2] - q->bmin[2]) * 0.5;
      if (tmp > scale) {
	scale = tmp;
      }
    }
    scale = 0.95 / scale;

    glPushMatrix();
    system_transformation();
    glScaled(move->zoom * scale, move->zoom * scale, move->zoom * scale);
    glTranslated(move->rightleft, move->updown, 0.0);

    /* Put bounding box midpoint in (0,0,0) */
    if (q->dim == 1) {
      glTranslated(-(q->bmax[0] + q->bmin[0]) * 0.5, 0.0, 0.0);
    }
    else {
      if (q->dim == 2) {
	glTranslated(-(q->bmax[0] + q->bmin[0]) * 0.5,
		     -(q->bmax[1] + q->bmin[1]) * 0.5, 0.0);
      }
      else {
	glRotated(move->angle_x, 1.0, 0.0, 0.0);
	glRotated(move->angle_y, 0.0, 1.0, 0.0);

	glTranslated(-(q->bmax[0] + q->bmin[0]) * 0.5,
		     -(q->bmax[1] + q->bmin[1]) * 0.5,
		     -(q->bmax[2] + q->bmin[2]) * 0.5);
      }
    }

    count = 0;
    draw_two_cluster(bbl, count);

    glPopMatrix();

  }

  glFlush();
  glutSwapBuffers();

  if (take_a_shot[1] == true) {

    static GLubyte *pixels;	/* Pointer for pixels */
    uint      paper_width = glutGet(GLUT_WINDOW_WIDTH);
    uint      paper_height = glutGet(GLUT_WINDOW_HEIGHT);

    pixels = allocmem(format_bytes * paper_height * paper_width);	/* Storage for the pixel */
    /* Reads the current buffer and save the pixels. RGBA has unsigned int type. */
    glReadPixels(0, 0, paper_width, paper_height, format, GL_UNSIGNED_BYTE,
		 pixels);

#ifdef USE_LIBPNG
    png_byte *png_bytes = NULL;	/* Memory for the png */
    png_byte **png_rows = NULL;	/* Array of memory for ever row */

    save_current_content_png("snapshot", paper_width, paper_height, pixels,
			     &png_bytes, &png_rows);
    free(png_bytes);
    free(png_rows);
#else
    save_current_content_ppm("snapshot", paper_width, paper_height, 255,
			     pixels);
#endif
    snapshots++;

    free(pixels);		/* free memory */
    take_a_shot[1] = false;	/* reset start values */
  }
}

  /* Drawing a given directional cluster tree */
static void
display_dcluster()
{

  uint      count;
  double    scale, tmp;
  double    zoom = move_1->zoom;
  double    rightleft = move_1->rightleft;
  double    updown = move_1->updown;
  double    angle_x = move_1->angle_x;
  double    angle_y = move_1->angle_y;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if (bbox == true) {
    bbox = false;
  }

  if (visu == 1) {
    scale = (t->bmax[0] - t->bmin[0]) * 0.5;
    if (t->dim >= 2) {
      tmp = (t->bmax[1] - t->bmin[1]) * 0.5;
      if (tmp > scale) {
	scale = tmp;
      }
    }
    if (t->dim == 3) {
      tmp = (t->bmax[2] - t->bmin[2]) * 0.5;
      if (tmp > scale) {
	scale = tmp;
      }
    }
    scale = 0.95 / scale;

    glPushMatrix();
    system_transformation();
    glScaled(zoom * scale, zoom * scale, zoom * scale);
    glTranslated(rightleft, updown, 0.0);

    if (t->dim == 1) {
      glTranslated(-(t->bmax[0] + t->bmin[0]) * 0.5, 0.0, 0.0);
    }
    else {
      if (t->dim == 2) {
	glTranslated(-(t->bmax[0] + t->bmin[0]) * 0.5,
		     -(t->bmax[1] + t->bmin[1]) * 0.5, 0.0);
      }
      else {
	glRotated(angle_x, 1.0, 0.0, 0.0);
	glRotated(angle_y, 0.0, 1.0, 0.0);
	glTranslated(-(t->bmax[0] + t->bmin[0]) * 0.5,
		     -(t->bmax[1] + t->bmin[1]) * 0.5,
		     -(t->bmax[2] + t->bmin[2]) * 0.5);
      }
    }
    count = 0;
    draw_dcluster_level(t, count, 0);

    glPopMatrix();
  }

  glFlush();
  glutSwapBuffers();

  if (take_a_shot[1] == true) {

    static GLubyte *pixels;	/* Pointer for pixels */
    uint      paper_width = glutGet(GLUT_WINDOW_WIDTH);
    uint      paper_height = glutGet(GLUT_WINDOW_HEIGHT);

    pixels = allocmem(format_bytes * paper_height * paper_width);	/* Storage for the pixel */
    /* Reads the current buffer and save the pixels. RGBA has unsigned int type. */
    glReadPixels(0, 0, paper_width, paper_height, format, GL_UNSIGNED_BYTE,
		 pixels);

#ifdef USE_LIBPNG
    png_byte *png_bytes = NULL;	/* Memory for the png */
    png_byte **png_rows = NULL;	/* Array of memory for ever row */

    save_current_content_png("snapshot", paper_width, paper_height, pixels,
			     &png_bytes, &png_rows);
    free(png_bytes);
    free(png_rows);
#else
    save_current_content_ppm("snapshot", paper_width, paper_height, 255,
			     pixels);
#endif
    snapshots++;

    free(pixels);		/* free memory */
    take_a_shot[1] = false;	/* reset start values */
  }
}

  /* Drawing a given cluster tree */
static void
display_cluster()
{

  uint      count;
  double    scale, tmp;
  double    zoom = move_1->zoom;
  double    rightleft = move_1->rightleft;
  double    updown = move_1->updown;
  double    angle_x = move_1->angle_x;
  double    angle_y = move_1->angle_y;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if (bbox == true) {
    bbox = false;
  }

  if (visu == 1) {
    scale = (ct->bmax[0] - ct->bmin[0]) * 0.5;
    if (ct->dim >= 2) {
      tmp = (ct->bmax[1] - ct->bmin[1]) * 0.5;
      if (tmp > scale) {
	scale = tmp;
      }
    }
    if (ct->dim == 3) {
      tmp = (ct->bmax[2] - ct->bmin[2]) * 0.5;
      if (tmp > scale) {
	scale = tmp;
      }
    }
    scale = 0.95 / scale;

    glPushMatrix();
    system_transformation();
    glScaled(zoom * scale, zoom * scale, zoom * scale);
    glTranslated(rightleft, updown, 0.0);

    if (ct->dim == 1) {
      glTranslated(-(ct->bmax[0] + ct->bmin[0]) * 0.5, 0.0, 0.0);
    }
    else {
      if (ct->dim == 2) {
	glTranslated(-(ct->bmax[0] + ct->bmin[0]) * 0.5,
		     -(ct->bmax[1] + ct->bmin[1]) * 0.5, 0.0);
      }
      else {
	glRotated(angle_x, 1.0, 0.0, 0.0);
	glRotated(angle_y, 0.0, 1.0, 0.0);
	glTranslated(-(ct->bmax[0] + ct->bmin[0]) * 0.5,
		     -(ct->bmax[1] + ct->bmin[1]) * 0.5,
		     -(ct->bmax[2] + ct->bmin[2]) * 0.5);
      }
    }
    count = 0;
    draw_cluster_level(ct, count, 0);

    glPopMatrix();
  }

  glFlush();
  glutSwapBuffers();

  if (take_a_shot[1] == true) {

    static GLubyte *pixels;	/* Pointer for pixels */
    uint      paper_width = glutGet(GLUT_WINDOW_WIDTH);
    uint      paper_height = glutGet(GLUT_WINDOW_HEIGHT);

    pixels = allocmem(format_bytes * paper_height * paper_width);	/* Storage for the pixel */
    /* Reads the current buffer and save the pixels. RGBA has unsigned int type. */
    glReadPixels(0, 0, paper_width, paper_height, format, GL_UNSIGNED_BYTE,
		 pixels);

#ifdef USE_LIBPNG
    png_byte *png_bytes = NULL;	/* Memory for the png */
    png_byte **png_rows = NULL;	/* Array of memory for ever row */

    save_current_content_png("snapshot", paper_width, paper_height, pixels,
			     &png_bytes, &png_rows);
    free(png_bytes);
    free(png_rows);
#else
    save_current_content_ppm("snapshot", paper_width, paper_height, 255,
			     pixels);
#endif
    snapshots++;

    free(pixels);		/* free memory */
    take_a_shot[1] = false;	/* reset start values */
  }
}

  /* Drawing a single block of a given directional block tree */
static void
display_certain_dblock()
{

  uint      count;
  uint      roffsize = 0;
  uint      coffsize = 0;
  uint      mark = 0;
  pcdblock  found;
  double    scale, tmp, tmp1, tmp2;
  double    zoom = move_1->zoom;
  double    rightleft = move_1->rightleft;
  double    updown = move_1->updown;
  double    angle_x = move_1->angle_x;
  double    angle_y = move_1->angle_y;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if (visu == 1) {
    scale = (t->bmax[0] - t->bmin[0]) * 0.5;
    if (t->dim >= 2) {
      tmp = (t->bmax[1] - t->bmin[1]) * 0.5;
      if (tmp > scale) {
	scale = tmp;
      }
    }
    if (t->dim == 3) {
      tmp = (t->bmax[2] - t->bmin[2]) * 0.5;
      if (tmp > scale) {
	scale = tmp;
      }
    }
    tmp = (s->bmax[0] - s->bmin[0]) * 0.5;
    if (tmp > scale) {
      scale = tmp;
    }
    if (s->dim >= 2) {
      tmp = (s->bmax[1] - s->bmin[1]) * 0.5;
      if (tmp > scale) {
	scale = tmp;
      }
    }
    if (s->dim == 3) {
      tmp = (s->bmax[2] - s->bmin[2]) * 0.5;
      if (tmp > scale) {
	scale = tmp;
      }
    }
    scale = 0.95 / scale;

    glPushMatrix();
    system_transformation();
    glScaled(zoom * scale, zoom * scale, zoom * scale);
    glTranslated(rightleft, updown, 0.0);

    if ((t->dim == 3) || (s->dim == 3)) {
      glRotated(angle_x, 1.0, 0.0, 0.0);
      glRotated(angle_y, 0.0, 1.0, 0.0);
    }

    /* Find Cluster with highest dimension and put bounding box midpoint in (0,0,0) */

    if (t->dim == s->dim) {
      /* 1D */
      if (t->dim == 1) {
	tmp =
	  (t->bmax[0] + t->bmin[0] >
	   s->bmax[0] + s->bmin[0] ? t->bmax[0] + t->bmin[0] : s->bmax[0] +
	   s->bmin[0]);
	glTranslated(-(tmp) * 0.5, 0.0, 0.0);
      }
      else {
	/* 2D */
	if (t->dim == 2) {
	  tmp =
	    (t->bmax[0] + t->bmin[0] >
	     s->bmax[0] + s->bmin[0] ? t->bmax[0] + t->bmin[0] : s->bmax[0] +
	     s->bmin[0]);
	  tmp1 =
	    (t->bmax[1] + t->bmin[1] >
	     s->bmax[1] + s->bmin[1] ? t->bmax[1] + t->bmin[1] : s->bmax[1] +
	     s->bmin[1]);
	  glTranslated(-(tmp) * 0.5, -(tmp1) * 0.5, 0.0);
	}
	else {
	  /* 3D */
	  tmp =
	    (t->bmax[0] + t->bmin[0] >
	     s->bmax[0] + s->bmin[0] ? t->bmax[0] + t->bmin[0] : s->bmax[0] +
	     s->bmin[0]);
	  tmp1 =
	    (t->bmax[1] + t->bmin[1] >
	     s->bmax[1] + s->bmin[1] ? t->bmax[1] + t->bmin[1] : s->bmax[1] +
	     s->bmin[1]);
	  tmp2 =
	    (t->bmax[2] + t->bmin[2] >
	     s->bmax[2] + s->bmin[2] ? t->bmax[2] + t->bmin[2] : s->bmax[2] +
	     s->bmin[2]);
	  glTranslated(-(tmp) * 0.5, -(tmp1) * 0.5, -(tmp2) * 0.5);
	}
      }
    }
    else {
      if (t->dim > s->dim) {
	if (t->dim == 2) {
	  glTranslated(-(t->bmax[0] + t->bmin[0]) * 0.5,
		       -(t->bmax[1] + t->bmin[1]) * 0.5, 0.0);
	}
	else {
	  glTranslated(-(t->bmax[0] + t->bmin[0]) * 0.5,
		       -(t->bmax[1] + t->bmin[1]) * 0.5,
		       -(t->bmax[2] + t->bmin[2]) * 0.5);
	}
      }
      else {
	if (s->dim == 2) {
	  glTranslated(-(s->bmax[0] + s->bmin[0]) * 0.5,
		       -(s->bmax[1] + s->bmin[1]) * 0.5, 0.0);
	}
	else {
	  /* 3D */
	  glTranslated(-(s->bmax[0] + s->bmin[0]) * 0.5,
		       -(s->bmax[1] + s->bmin[1]) * 0.5,
		       -(s->bmax[2] + s->bmin[2]) * 0.5);
	}

      }
    }
    if (up == 0) {
      count = 0;
      found = draw_certain_dblock(bl, count, roffsize, coffsize, mark);
      if (found == NULL) {
	printf("Please choose a block to visualize.\n");
      }
      move_1->dblock = found;
    }
    else {
      if (up == 2 && move_1->dblock == bl) {
	color_admissible_d(bl);
	if (bl->rc->dim == 1) {
	  draw_1c(bl->rc, level2 - level, 0);
	}
	if (bl->rc->dim == 2) {
	  draw_2c(bl->rc, 0);
	}
	if (bl->rc->dim == 3) {
	  draw_3c(bl->rc, 0);
	}

	color_admissible_d(bl);
	if (bl->cc->dim == 1) {
	  draw_1c(bl->cc, level2 - level, 1);
	}
	if (bl->cc->dim == 2) {
	  draw_2c(bl->cc, 1);
	}
	if (bl->cc->dim == 3) {
	  draw_3c(bl->cc, 1);
	}
	up = 1;
      }
      else {
	draw_certain_dblock_father(bl, mark);
	up = 1;
      }
    }


    glPopMatrix();
  }

  glFlush();
  glutSwapBuffers();

  if (bbox == true) {
    bbox = false;
  }

  if (take_a_shot[1] == true) {

    static GLubyte *pixels;	/* Pointer for pixels */
    uint      paper_width = glutGet(GLUT_WINDOW_WIDTH);
    uint      paper_height = glutGet(GLUT_WINDOW_HEIGHT);

    pixels = allocmem(format_bytes * paper_height * paper_width);	/* Storage for the pixel */
    /* Reads the current buffer and save the pixels. RGBA has unsigned int type. */
    glReadPixels(0, 0, paper_width, paper_height, format, GL_UNSIGNED_BYTE,
		 pixels);

#ifdef USE_LIBPNG
    png_byte *png_bytes = NULL;	/* Memory for the png */
    png_byte **png_rows = NULL;	/* Array of memory for ever row */

    save_current_content_png("snapshot", paper_width, paper_height, pixels,
			     &png_bytes, &png_rows);
    free(png_bytes);
    free(png_rows);
#else
    save_current_content_ppm("snapshot", paper_width, paper_height, 255,
			     pixels);
#endif
    snapshots++;

    free(pixels);		/* free memory */
    take_a_shot[1] = false;	/* reset start values */
  }
}

  /* Drawing a single block of a given block tree */
static void
display_certain_block()
{

  uint      count;
  uint      roffsize = 0;
  uint      coffsize = 0;
  uint      mark = 0;
  pcblock   found;
  double    scale, tmp, tmp1, tmp2;
  double    zoom = move_1->zoom;
  double    rightleft = move_1->rightleft;
  double    updown = move_1->updown;
  double    angle_x = move_1->angle_x;
  double    angle_y = move_1->angle_y;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if (visu == 1) {
    scale = (ct->bmax[0] - ct->bmin[0]) * 0.5;
    if (ct->dim >= 2) {
      tmp = (ct->bmax[1] - ct->bmin[1]) * 0.5;
      if (tmp > scale) {
	scale = tmp;
      }
      if (ct->dim == 3) {
	tmp = (ct->bmax[2] - ct->bmin[2]) * 0.5;
	if (tmp > scale) {
	  scale = tmp;
	}
      }
    }
    tmp = (cs->bmax[0] - cs->bmin[0]) * 0.5;
    if (tmp > scale) {
      scale = tmp;
    }
    if (cs->dim >= 2) {
      tmp = (cs->bmax[1] - cs->bmin[1]) * 0.5;
      if (tmp > scale) {
	scale = tmp;
      }
      if (cs->dim == 3) {
	tmp = (cs->bmax[2] - cs->bmin[2]) * 0.5;
	if (tmp > scale) {
	  scale = tmp;
	}
      }
    }
    scale = 0.95 / scale;

    glPushMatrix();
    system_transformation();
    glScaled(zoom * scale, zoom * scale, zoom * scale);
    glTranslated(rightleft, updown, 0.0);


    if ((ct->dim == 3) || (cs->dim == 3)) {
      glRotated(angle_x, 1.0, 0.0, 0.0);
      glRotated(angle_y, 0.0, 1.0, 0.0);
    }

    /* Find Cluster with highest dimension and put bounding box midpoint in (0,0,0) */

    if (ct->dim == cs->dim) {
      /* 1D */
      if (ct->dim == 1) {
	tmp =
	  (ct->bmax[0] + ct->bmin[0] >
	   cs->bmax[0] + cs->bmin[0] ? ct->bmax[0] +
	   ct->bmin[0] : cs->bmax[0] + cs->bmin[0]);
	glTranslated(-(tmp) * 0.5, 0.0, 0.0);
      }
      else {
	/* 2D */
	if (ct->dim == 2) {
	  tmp =
	    (ct->bmax[0] + ct->bmin[0] >
	     cs->bmax[0] + cs->bmin[0] ? ct->bmax[0] +
	     ct->bmin[0] : cs->bmax[0] + cs->bmin[0]);
	  tmp1 =
	    (ct->bmax[1] + ct->bmin[1] >
	     cs->bmax[1] + cs->bmin[1] ? ct->bmax[1] +
	     ct->bmin[1] : cs->bmax[1] + cs->bmin[1]);
	  glTranslated(-(tmp) * 0.5, -(tmp1) * 0.5, 0.0);
	}
	else {
	  /* 3D */
	  tmp =
	    (ct->bmax[0] + ct->bmin[0] >
	     cs->bmax[0] + cs->bmin[0] ? ct->bmax[0] +
	     ct->bmin[0] : cs->bmax[0] + cs->bmin[0]);
	  tmp1 =
	    (ct->bmax[1] + ct->bmin[1] >
	     cs->bmax[1] + cs->bmin[1] ? ct->bmax[1] +
	     ct->bmin[1] : cs->bmax[1] + cs->bmin[1]);
	  tmp2 =
	    (ct->bmax[2] + ct->bmin[2] >
	     cs->bmax[2] + cs->bmin[2] ? ct->bmax[2] +
	     ct->bmin[2] : cs->bmax[2] + cs->bmin[2]);
	  glTranslated(-(tmp) * 0.5, -(tmp1) * 0.5, -(tmp2) * 0.5);
	}
      }
    }
    else {
      if (ct->dim > cs->dim) {
	if (ct->dim == 2) {
	  glTranslated(-(ct->bmax[0] + ct->bmin[0]) * 0.5,
		       -(ct->bmax[1] + ct->bmin[1]) * 0.5, 0.0);
	}
	else {
	  glTranslated(-(ct->bmax[0] + ct->bmin[0]) * 0.5,
		       -(ct->bmax[1] + ct->bmin[1]) * 0.5,
		       -(ct->bmax[2] + ct->bmin[2]) * 0.5);
	}
      }
      else {
	if (cs->dim == 2) {
	  glTranslated(-(cs->bmax[0] + cs->bmin[0]) * 0.5,
		       -(cs->bmax[1] + cs->bmin[1]) * 0.5, 0.0);
	}
	else {
	  /* 3D */
	  glTranslated(-(cs->bmax[0] + cs->bmin[0]) * 0.5,
		       -(cs->bmax[1] + cs->bmin[1]) * 0.5,
		       -(cs->bmax[2] + cs->bmin[2]) * 0.5);
	}

      }
    }

    if (up == 0) {
      count = 0;
      found = draw_certain_block(bbl, count, roffsize, coffsize, mark);
      if (found == NULL) {
	printf("Please choose a block to visualize.\n");
      }
      move_1->block = found;
    }
    else {
      if (up == 2 && move_1->block == bbl) {
	color_admissible(bbl);
	if (bbl->rc->dim == 1) {
	  draw_1cc(bbl->rc, level2 - level, 0);
	}
	if (bbl->rc->dim == 2) {
	  draw_2cc(bbl->rc, 0);
	}
	if (bbl->rc->dim == 3) {
	  draw_3cc(bbl->rc, 0);
	}

	color_admissible(bbl);
	if (bbl->cc->dim == 1) {
	  draw_1cc(bbl->cc, level2 - level, 1);
	}
	if (bbl->cc->dim == 2) {
	  draw_2cc(bbl->cc, 1);
	}
	if (bbl->cc->dim == 3) {
	  draw_3cc(bbl->cc, 1);
	}
	up = 1;
      }
      else {
	draw_certain_block_father(bbl, mark);
	up = 1;
      }
    }
    glPopMatrix();
  }

  glFlush();
  glutSwapBuffers();

  if (bbox == true) {
    bbox = false;
  }

  if (take_a_shot[1] == true) {

    static GLubyte *pixels;	/* Pointer for pixels */
    uint      paper_width = glutGet(GLUT_WINDOW_WIDTH);
    uint      paper_height = glutGet(GLUT_WINDOW_HEIGHT);

    pixels = allocmem(format_bytes * paper_height * paper_width);	/* Storage for the pixel */
    /* Reads the current buffer and save the pixels. RGBA has unsigned int type. */
    glReadPixels(0, 0, paper_width, paper_height, format, GL_UNSIGNED_BYTE,
		 pixels);

#ifdef USE_LIBPNG
    png_byte *png_bytes = NULL;	/* Memory for the png */
    png_byte **png_rows = NULL;	/* Array of memory for ever row */

    save_current_content_png("snapshot", paper_width, paper_height, pixels,
			     &png_bytes, &png_rows);
    free(png_bytes);
    free(png_rows);
#else
    save_current_content_ppm("snapshot", paper_width, paper_height, 255,
			     pixels);
#endif
    snapshots++;

    free(pixels);		/* free memory */
    take_a_shot[1] = false;	/* reset start values */
  }
}

  /* Drawing table of a two dimensional triangulation */
static void
display_table_tri()
{

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  uint      step = 5;

  uint      i, j;
  char      stroke[1][20];
  char      stroke2[3][20];
  double    scale = 0.75 / (step + 1);
  double    eps = (2.0 / (step + 1) - scale) / 2.0;
  uint      number = tr->triangles;
  uint      v[3];

  glClearColor(1.0, 1.0, 1.0, 1.0);
  system_transformation();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslated(-1.0, 1.0, 0.0);

  /* Drawing left table with number of triangle and coordinates */
  for (i = 0; i < step; i++) {
    glBegin(GL_LINE_STRIP);
    glColor3d(0.1, 0.3, 0.3);
    glVertex2d(0.0, -2.0 * i / (step + 1));
    glVertex2d(1.0, -2.0 * i / (step + 1));
    glVertex2d(1.0, -2.0 * (i + 1) / (step + 1));
    glVertex2d(0.0, -2.0 * (i + 1) / (step + 1));
    glVertex2d(0.0, -2.0 * i / (step + 1));
    glEnd();

    glPushMatrix();
    if ((i + move_t->page * 10) < number) {
      /* Number of triangle */
      glTranslatef(0.0, -2.0 * (i + 1.0) / (step + 1.0) + eps, 0.0);
      glScalef(scale * 1.0 / 152.0, scale * 1.0 / 152.0, scale * 1.0 / 152.0);
      glScalef(0.5, 0.5, 0.5);
      sprintf(stroke[0], "%d", (i + move_t->page * 10));
      if (move_t->target == (i + move_t->page * 10)) {
	glColor3d(0.0, 0.0, 0.8);
	glLineWidth(3);
      }
      else {
	glColor3d(0.8, 0.0, 0.0);
	glLineWidth(2);
      }
      for (j = 0; j < (int) strlen(stroke[0]); j++) {
	glutStrokeCharacter(GLUT_STROKE_ROMAN, stroke[0][j]);
      }
      glColor3d(0.0, 0.0, 0.0);
      glPushMatrix();
      /* Coordinate of triangle {(x, ., .) */
      glTranslatef(200.0 * eps,
		   -2.0 * (i + 1.0) / (step + 1.0) + 1500.0 * eps, 0.0);
      glScalef(0.5, 0.5, 0.5);
      getvertices_tri2d(tr, i + move_t->page, v);
      sprintf(stroke2[0], "{(%.2e, %.2e),", tr_orig->x[v[0]][0],
	      tr_orig->x[v[0]][1]);
      glLineWidth(2);

      for (j = 0; j < (int) strlen(stroke2[0]); j++) {
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, stroke2[0][j]);
      }

      glPopMatrix();

      glPushMatrix();
      /* Coordinate of triangle (., x, .) */
      glTranslatef(200.0 * eps, -2.0 * (i + 1.0) / (step + 1.0) + 500.0 * eps,
		   0.0);
      glScalef(0.5, 0.5, 0.5);
      sprintf(stroke2[1], "(%.2e, %.2e),", tr_orig->x[v[1]][0],
	      tr_orig->x[v[1]][1]);
      glLineWidth(2);

      for (j = 0; j < (int) strlen(stroke2[1]); j++) {
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, stroke2[1][j]);
      }

      glPopMatrix();

      glPushMatrix();
      /* Coordinate of triangle ( ., ., x)} */
      glTranslatef(200.0 * eps, -2.0 * (i + 1.0) / (step + 1.0) - 500.0 * eps,
		   0.0);
      glScalef(0.5, 0.5, 0.5);
      sprintf(stroke2[2], "(%.2e, %.2e)}", tr_orig->x[v[2]][0],
	      tr_orig->x[v[2]][1]);
      glLineWidth(2);

      for (j = 0; j < (int) strlen(stroke2[2]); j++) {
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, stroke2[2][j]);
      }

      glPopMatrix();
    }
    glPopMatrix();

    /* Drawing right table with number of triangle and coordinates */
    glBegin(GL_LINE_STRIP);
    glColor3d(0.1, 0.3, 0.3);
    glVertex2d(1.0, -2.0 * i / (step + 1));
    glVertex2d(2.0, -2.0 * i / (step + 1));
    glVertex2d(2.0, -2.0 * (i + 1) / (step + 1));
    glVertex2d(1.0, -2.0 * (i + 1) / (step + 1));
    glEnd();

    glPushMatrix();
    if ((i + 5 + move_t->page * 10) < number) {
      /* Number of triangle */
      glTranslatef(1.0, -2.0 * (i + 1.0) / (step + 1.0) + eps, 0.0);
      glScalef(scale * 1.0 / 152.0, scale * 1.0 / 152.0, scale * 1.0 / 152.0);
      glScalef(0.5, 0.5, 0.5);
      sprintf(stroke[0], "%d", (i + 5 + move_t->page * 10));
      if (move_t->target == (i + 5 + move_t->page * 10)) {
	glColor3d(0.0, 0.0, 0.8);
	glLineWidth(3);
      }
      else {
	glColor3d(0.8, 0.0, 0.0);
	glLineWidth(2);
      }

      for (j = 0; j < (int) strlen(stroke[0]); j++) {
	glutStrokeCharacter(GLUT_STROKE_ROMAN, stroke[0][j]);
      }
      glColor3d(0.0, 0.0, 0.0);
      glPushMatrix();
      /* Coordinate of triangle {(x, ., .) */
      glTranslatef(200.0 * eps,
		   -2.0 * (i + 1.0) / (step + 1.0) + 1500.0 * eps, 0.0);
      glScalef(0.5, 0.5, 0.5);
      getvertices_tri2d(tr, i + move_t->page, v);
      sprintf(stroke2[0], "{(%.2e, %.2e),", tr_orig->x[v[0]][0],
	      tr_orig->x[v[0]][1]);
      glLineWidth(2);

      for (j = 0; j < (int) strlen(stroke2[0]); j++) {
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, stroke2[0][j]);
      }

      glPopMatrix();

      glPushMatrix();
      /* Coordinate of triangle (., x, .) */
      glTranslatef(200.0 * eps, -2.0 * (i + 1.0) / (step + 1.0) + 500.0 * eps,
		   0.0);
      glScalef(0.5, 0.5, 0.5);
      sprintf(stroke2[1], "(%.2e, %.2e),", tr_orig->x[v[1]][0],
	      tr_orig->x[v[1]][1]);
      glLineWidth(2);

      for (j = 0; j < (int) strlen(stroke2[1]); j++) {
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, stroke2[1][j]);
      }

      glPopMatrix();

      glPushMatrix();
      /* Coordinate of triangle ( ., ., x)} */
      glTranslatef(200.0 * eps, -2.0 * (i + 1.0) / (step + 1.0) - 500.0 * eps,
		   0.0);
      glScalef(0.5, 0.5, 0.5);
      sprintf(stroke2[2], "(%.2e, %.2e)}", tr_orig->x[v[2]][0],
	      tr_orig->x[v[2]][1]);
      glLineWidth(2);

      for (j = 0; j < (int) strlen(stroke2[2]); j++) {
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, stroke2[2][j]);
      }

      glPopMatrix();
    }
    glPopMatrix();

  }

  glLineWidth(1);

  glBegin(GL_LINE_STRIP);
  glColor3d(0.1, 0.3, 0.3);
  glVertex2d(0.0, -2.0 * 5 / (step + 1));
  glVertex2d(1.0, -2.0 * 5 / (step + 1));
  glVertex2d(1.0, -2.0);
  glVertex2d(0.0, -2.0);
  glVertex2d(0.0, -2.0 * 5 / (step + 1));

  glEnd();

  glPushMatrix();
  glTranslatef(0.0 + eps, -2.0 * (5 + 1.0) / (step + 1.0) + eps, 0.0);
  glScalef(0.75 * scale * 1.0 / 152.0, 0.75 * scale * 1.0 / 152.0,
	   0.75 * scale * 1.0 / 152.0);
  sprintf(stroke[0], "triangles %d", number);
  glColor3d(0.0, 0.0, 0.0);
  glLineWidth(2);

  for (j = 0; j < (int) strlen(stroke[0]); j++) {
    glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, stroke[0][j]);
  }
  glPopMatrix();

  glLineWidth(1);

  glBegin(GL_LINE_STRIP);
  glColor3d(0.1, 0.3, 0.3);
  glVertex2d(1.0, -2.0 * 5 / (step + 1));
  glVertex2d(2.0, -2.0 * 5 / (step + 1));
  glVertex2d(2.0, -2.0);
  glVertex2d(1.0, -2.0);
  glEnd();

  glTranslatef(1.0, -2.0 * (6.0) / (step + 1.0), 0.0);

  glColor3d(0.0, 1.0, 0.0);
  glBegin(GL_TRIANGLES);
  glVertex3d(0.0, 1.0 / 6.0, 0.0);
  glVertex3d(0.5 - eps * 0.1, 0.0, 0.0);
  glVertex3d(0.5 - eps * 0.1, 1.0 / 3.0, 0.0);
  glEnd();

  glBegin(GL_TRIANGLES);
  glVertex3d(0.5 + eps * 0.1, 0.0, 0.0);
  glVertex3d(1.0, 1.0 / 6.0, 0.0);
  glVertex3d(0.5 + eps * 0.1, 1.0 / 3.0, 0.0);
  glEnd();

  glColor3d(1.0, 1.0, 1.0);
  glBegin(GL_LINE_STRIP);
  glVertex3d(0.0, 1.0 / 6.0, 0.0);
  glVertex3d(0.5 - eps * 0.1, 0.0, 0.0);
  glVertex3d(0.5 - eps * 0.1, 1.0 / 3.0, 0.0);
  glEnd();

  glBegin(GL_LINE_STRIP);
  glVertex3d(0.5 + eps * 0.1, 0.0, 0.0);
  glVertex3d(1.0, 1.0 / 6.0, 0.0);
  glVertex3d(0.5 + eps * 0.1, 1.0 / 3.0, 0.0);
  glEnd();

  glFlush();
  glutSwapBuffers();

  if (take_a_shot[0] == true) {

    static GLubyte *pixels;	/* Pointer for pixels */
    uint      paper_width = glutGet(GLUT_WINDOW_WIDTH);
    uint      paper_height = glutGet(GLUT_WINDOW_HEIGHT);

    pixels = allocmem(format_bytes * paper_height * paper_width);	/* Storage for the pixel */
    /* Reads the current buffer and save the pixels. RGBA has unsigned int type. */
    glReadPixels(0, 0, paper_width, paper_height, format, GL_UNSIGNED_BYTE,
		 pixels);

#ifdef USE_LIBPNG
    png_byte *png_bytes = NULL;	/* Memory for the png */
    png_byte **png_rows = NULL;	/* Array of memory for ever row */

    save_current_content_png("snapshot", paper_width, paper_height, pixels,
			     &png_bytes, &png_rows);
    free(png_bytes);
    free(png_rows);
#else
    save_current_content_ppm("snapshot", paper_width, paper_height, 255,
			     pixels);
#endif
    snapshots++;

    free(pixels);		/* free memory */
    take_a_shot[0] = false;	/* reset start values */
  }
}

  /* Drawing a certain triangle */
static void
display_certain_triangle()
{

  uint      v[3];

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glPushMatrix();
  system_transformation();
  glTranslated(move_1->rightleft, move_1->updown, 0.0);
  glRotated(move_1->angle_x, 1.0, 0.0, 0.0);
  glRotated(move_1->angle_y, 0.0, 1.0, 0.0);
  if (coord == true) {
    glPushMatrix();
    coordinate();
    glPopMatrix();
  }
  glScaled(move_1->zoom, move_1->zoom, move_1->zoom);
  if (origin == true) {
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPointSize(6.0);
    glColor4d(0.8, 0.8, 0.8, 0.8);
    glBegin(GL_POINTS);
    glVertex3d(0.0, 0.0, 0.0);
    glEnd();
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
  }

  if (mode == 1) {
    getvertices_tri2d(tr, move_t->target, v);
    glColor3d(0.8, 0.8, 0.8);
    glBegin(GL_LINE_LOOP);
    glVertex3d(tr->x[v[0]][0], tr->x[v[0]][1], 0.0);
    glVertex3d(tr->x[v[1]][0], tr->x[v[1]][1], 0.0);
    glVertex3d(tr->x[v[2]][0], tr->x[v[2]][1], 0.0);
    glEnd();
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4d(0.0, 0.2, 0.8, 0.3);
    glBegin(GL_TRIANGLES);
    glVertex3d(tr->x[v[0]][0], tr->x[v[0]][1], 0.0);
    glVertex3d(tr->x[v[1]][0], tr->x[v[1]][1], 0.0);
    glVertex3d(tr->x[v[2]][0], tr->x[v[2]][1], 0.0);
    glEnd();
    glDisable(GL_BLEND);
  }

  outline(0);

  glPopMatrix();

  glFlush();
  glutSwapBuffers();

  if (take_a_shot[1] == true) {

    static GLubyte *pixels;	/* Pointer for pixels */
    uint      paper_width = glutGet(GLUT_WINDOW_WIDTH);
    uint      paper_height = glutGet(GLUT_WINDOW_HEIGHT);

    pixels = allocmem(format_bytes * paper_height * paper_width);	/* Storage for the pixel */
    /* Reads the current buffer and save the pixels. RGBA has unsigned int type. */
    glReadPixels(0, 0, paper_width, paper_height, format, GL_UNSIGNED_BYTE,
		 pixels);

#ifdef USE_LIBPNG
    png_byte *png_bytes = NULL;	/* Memory for the png */
    png_byte **png_rows = NULL;	/* Array of memory for ever row */

    save_current_content_png("snapshot", paper_width, paper_height, pixels,
			     &png_bytes, &png_rows);
    free(png_bytes);
    free(png_rows);
#else
    save_current_content_ppm("snapshot", paper_width, paper_height, 255,
			     pixels);
#endif
    snapshots++;

    free(pixels);		/* free memory */
    take_a_shot[1] = false;	/* reset start values */
  }
}

  /* Drawing complete two dimensional triangulation */
static void
display_tri2d()
{

  uint      k, i;
  uint      v[3];

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glPushMatrix();
  system_transformation();
  glTranslated(move_1->rightleft, move_1->updown, 0.0);
  glRotated(move_1->angle_x, 1.0, 0.0, 0.0);
  glRotated(move_1->angle_y, 0.0, 1.0, 0.0);
  if (coord == 1) {
    glPushMatrix();
    coordinate();
    glPopMatrix();
  }
  glScaled(move_1->zoom, move_1->zoom, move_1->zoom);
  if (origin == true) {
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPointSize(6.0);
    glColor4d(0.8, 0.8, 0.8, 0.8);
    glBegin(GL_POINTS);
    glVertex3d(0.0, 0.0, 0.0);
    glEnd();
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
  }

  switch (mode) {		/* Drawing edges */
  case 0:

    glColor3d(1.0, 0.8, 0.0);
    glBegin(GL_LINES);

    for (k = 0; k < tr->edges; k++) {
      glVertex3d(tr->x[tr->e[k][0]][0], tr->x[tr->e[k][0]][1], 0.0);
      glVertex3d(tr->x[tr->e[k][1]][0], tr->x[tr->e[k][1]][1], 0.0);
    }
    glEnd();

    break;

  case 1:

    for (k = 0; k < tr->triangles; k++) {
      getvertices_tri2d(tr, k, v);
      glColor3d(0.8, 0.8, 0.8);
      glBegin(GL_LINES);
      for (i = 0; i < 3; i++) {
	glVertex3d(tr->x[v[i]][0], tr->x[v[i]][1], 0.0);
	glVertex3d(tr->x[v[(i + 1) % 3]][0], tr->x[v[(i + 1) % 3]][1], 0.0);
      }
      glEnd();
    }
    for (k = 0; k < tr->triangles; k++) {
      getvertices_tri2d(tr, k, v);
      glColor3d(0.0, 0.0, 0.8);
      glBegin(GL_TRIANGLES);
      glVertex3d(tr->x[v[0]][0], tr->x[v[0]][1], 0.0);
      glVertex3d(tr->x[v[1]][0], tr->x[v[1]][1], 0.0);
      glVertex3d(tr->x[v[2]][0], tr->x[v[2]][1], 0.0);
      glEnd();

    }

    break;
  default:
    break;

  }

  if (points == 1) {
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPointSize(6.0);
    glColor3d(0.8, 0.8, 0.8);
    glBegin(GL_POINTS);
    for (k = 0; k < tr->vertices; k++) {
      glVertex3d(tr->x[k][0], tr->x[k][1], 0.0);
    }
    glEnd();
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
  }
  glPopMatrix();

  glFlush();
  glutSwapBuffers();

  if (take_a_shot[1] == true) {

    static GLubyte *pixels;	/* Pointer for pixels */
    uint      paper_width = glutGet(GLUT_WINDOW_WIDTH);
    uint      paper_height = glutGet(GLUT_WINDOW_HEIGHT);

    pixels = allocmem(format_bytes * paper_height * paper_width);	/* Storage for the pixel */
    /* Reads the current buffer and save the pixels. RGBA has unsigned int type. */
    glReadPixels(0, 0, paper_width, paper_height, format, GL_UNSIGNED_BYTE,
		 pixels);

#ifdef USE_LIBPNG
    png_byte *png_bytes = NULL;	/* Memory for the png */
    png_byte **png_rows = NULL;	/* Array of memory for ever row */

    save_current_content_png("snapshot", paper_width, paper_height, pixels,
			     &png_bytes, &png_rows);
    free(png_bytes);
    free(png_rows);
#else
    save_current_content_ppm("snapshot", paper_width, paper_height, 255,
			     pixels);
#endif
    snapshots++;

    free(pixels);		/* free memory */
    take_a_shot[1] = false;	/* reset start values */
  }


}

  /* Drawing a certain tetrahedra */
static void
display_certain_tetrahedra()
{

  uint      v[4];
  uint      i;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glPushMatrix();
  system_transformation();
  glTranslated(move_1->rightleft, move_1->updown, 0.0);
  glRotated(move_1->angle_x, 1.0, 0.0, 0.0);
  glRotated(move_1->angle_y, 0.0, 1.0, 0.0);
  if (coord == 1) {
    glPushMatrix();
    coordinate();
    glPopMatrix();
  }
  glScaled(move_1->zoom, move_1->zoom, move_1->zoom);

  if (origin == true) {
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPointSize(6.0);
    glColor4d(0.8, 0.8, 0.8, 0.8);
    glBegin(GL_POINTS);
    glVertex3d(0.0, 0.0, 0.0);
    glEnd();
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
  }

  if (mode == 1) {

    getvertices_tet3d(th, move_t->target, v);
    glColor3d(0.8, 0.8, 0.8);
    glBegin(GL_LINES);
    for (i = 1; i < 4; i++) {
      glVertex3d(th->x[v[0]][0], th->x[v[0]][1], th->x[v[0]][2]);
      glVertex3d(th->x[v[i]][0], th->x[v[i]][1], th->x[v[i]][2]);
    }
    for (i = 2; i <= 3; i++) {
      glVertex3d(th->x[v[1]][0], th->x[v[1]][1], th->x[v[1]][2]);
      glVertex3d(th->x[v[i]][0], th->x[v[i]][1], th->x[v[i]][2]);
    }
    glVertex3d(th->x[v[3]][0], th->x[v[3]][1], th->x[v[3]][2]);
    glVertex3d(th->x[v[2]][0], th->x[v[2]][1], th->x[v[2]][2]);
    glEnd();

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4d(0.0, 0.0, 0.8, 0.5);
    for (i = 0; i < 4; i++) {
      glBegin(GL_TRIANGLES);
      glVertex3d(th->x[v[i]][0], th->x[v[i]][1], th->x[v[i]][2]);
      glVertex3d(th->x[v[(i + 1) % 4]][0], th->x[v[(i + 1) % 4]][1],
		 th->x[v[(i + 1) % 4]][2]);
      glVertex3d(th->x[v[(i + 2) % 4]][0], th->x[v[(i + 2) % 4]][1],
		 th->x[v[(i + 2) % 4]][2]);
      glEnd();
    }
    glDisable(GL_BLEND);
  }

  outline(1);

  glPopMatrix();

  glFlush();
  glutSwapBuffers();

  if (take_a_shot[1] == true) {

    static GLubyte *pixels;	/* Pointer for pixels */
    uint      paper_width = glutGet(GLUT_WINDOW_WIDTH);
    uint      paper_height = glutGet(GLUT_WINDOW_HEIGHT);

    pixels = allocmem(format_bytes * paper_height * paper_width);	/* Storage for the pixel */
    /* Reads the current buffer and save the pixels. RGBA has unsigned int type. */
    glReadPixels(0, 0, paper_width, paper_height, format, GL_UNSIGNED_BYTE,
		 pixels);

#ifdef USE_LIBPNG
    png_byte *png_bytes = NULL;	/* Memory for the png */
    png_byte **png_rows = NULL;	/* Array of memory for ever row */

    save_current_content_png("snapshot", paper_width, paper_height, pixels,
			     &png_bytes, &png_rows);
    free(png_bytes);
    free(png_rows);
#else
    save_current_content_ppm("snapshot", paper_width, paper_height, 255,
			     pixels);
#endif
    snapshots++;

    free(pixels);		/* free memory */
    take_a_shot[1] = false;	/* reset start values */
  }

}

  /* Drawing table of a three dimensional triangulation */
static void
display_table_tet()
{

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  uint      step = 4;

  uint      i, j;
  char      stroke[1][20];
  char      stroke2[4][30];
  double    scale = 0.75 / (step + 1);
  double    eps = (2.0 / (step + 1) - scale) / 2.0;
  uint      number = th->tetrahedra;
  uint      v[4];

  glClearColor(1.0, 1.0, 1.0, 1.0);
  system_transformation();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslated(-1.0, 1.0, 0.0);

  /* Drawing left table with number of tetrahedra and coordinates */
  for (i = 0; i < step; i++) {
    glBegin(GL_LINE_STRIP);
    glColor3d(0.1, 0.3, 0.3);
    glVertex2d(0.0, -2.0 * i / (step + 1));
    glVertex2d(1.0, -2.0 * i / (step + 1));
    glVertex2d(1.0, -2.0 * (i + 1) / (step + 1));
    glVertex2d(0.0, -2.0 * (i + 1) / (step + 1));
    glVertex2d(0.0, -2.0 * i / (step + 1));
    glEnd();

    glPushMatrix();
    if ((i + move_t->page * 8) < number) {
      /* Number of tetrahedra */
      glTranslatef(0.0, -2.0 * (i + 1.0) / (step + 1.0) + eps, 0.0);
      glScalef(scale * 1.0 / 152.0, scale * 1.0 / 152.0, scale * 1.0 / 152.0);
      glScalef(0.5, 0.5, 0.5);
      sprintf(stroke[0], "%d", (i + move_t->page * 8));
      if (move_t->target == (i + move_t->page * 8)) {
	glColor3d(0.0, 0.0, 0.8);
	glLineWidth(3);
      }
      else {
	glColor3d(0.8, 0.0, 0.0);
	glLineWidth(2);
      }
      for (j = 0; j < (int) strlen(stroke[0]); j++) {
	glutStrokeCharacter(GLUT_STROKE_ROMAN, stroke[0][j]);
      }
      glColor3d(0.0, 0.0, 0.0);
      glPushMatrix();
      /* Coordinate of tetrahedra {(x, ., ., .) */
      glTranslatef(200.0 * eps,
		   -2.0 * (i + 1.0) / (step + 1.0) + 1600.0 * eps, 0.0);
      glScalef(0.5, 0.5, 0.5);
      getvertices_tet3d(th, i + move_t->page, v);
      sprintf(stroke2[0], "{(%.2e, %.2e, %.2e),", th_orig->x[v[0]][0],
	      th_orig->x[v[0]][1], th_orig->x[v[0]][2]);
      glLineWidth(2);

      for (j = 0; j < (int) strlen(stroke2[0]); j++) {
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, stroke2[0][j]);
      }

      glPopMatrix();

      glPushMatrix();
      /* Coordinate of tetrahedra (., x, ., .) */
      glTranslatef(200.0 * eps, -2.0 * (i + 1.0) / (step + 1.0) + 800.0 * eps,
		   0.0);
      glScalef(0.5, 0.5, 0.5);
      sprintf(stroke2[1], "(%.2e, %.2e, %.2e),", th_orig->x[v[1]][0],
	      th_orig->x[v[1]][1], th_orig->x[v[1]][2]);
      glLineWidth(2);

      for (j = 0; j < (int) strlen(stroke2[1]); j++) {
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, stroke2[1][j]);
      }

      glPopMatrix();

      glPushMatrix();
      /* Coordinate of tetrahedra ( ., ., x, .) */
      glTranslatef(200.0 * eps, -2.0 * (i + 1.0) / (step + 1.0) + eps, 0.0);
      glScalef(0.5, 0.5, 0.5);
      sprintf(stroke2[2], "(%.2e, %.2e, %.2e),", th_orig->x[v[2]][0],
	      th_orig->x[v[2]][1], th_orig->x[v[2]][2]);
      glLineWidth(2);

      for (j = 0; j < (int) strlen(stroke2[2]); j++) {
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, stroke2[2][j]);
      }

      glPopMatrix();

      glPushMatrix();
      /* Coordinate of tetrahedra ( ., ., ., x)} */
      glTranslatef(200.0 * eps, -2.0 * (i + 1.0) / (step + 1.0) - 800.0 * eps,
		   0.0);
      glScalef(0.5, 0.5, 0.5);
      sprintf(stroke2[3], "(%.2e, %.2e, %.2e)}", th_orig->x[v[3]][0],
	      th_orig->x[v[3]][1], th_orig->x[v[3]][2]);
      glLineWidth(2);

      for (j = 0; j < (int) strlen(stroke2[3]); j++) {
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, stroke2[3][j]);
      }

      glPopMatrix();
    }
    glPopMatrix();

    /* Drawing right table with number of tetrahedra and coordinates */
    glBegin(GL_LINE_STRIP);
    glColor3d(0.1, 0.3, 0.3);
    glVertex2d(1.0, -2.0 * i / (step + 1));
    glVertex2d(2.0, -2.0 * i / (step + 1));
    glVertex2d(2.0, -2.0 * (i + 1) / (step + 1));
    glVertex2d(1.0, -2.0 * (i + 1) / (step + 1));
    glEnd();

    glPushMatrix();
    if ((i + 4 + move_t->page * 8) < number) {
      /* Number of tetrahedra */
      glTranslatef(1.0, -2.0 * (i + 1.0) / (step + 1.0) + eps, 0.0);
      glScalef(scale * 1.0 / 152.0, scale * 1.0 / 152.0, scale * 1.0 / 152.0);
      glScalef(0.5, 0.5, 0.5);
      sprintf(stroke[0], "%d", (i + 4 + move_t->page * 8));
      if (move_t->target == (i + 4 + move_t->page * 8)) {
	glColor3d(0.0, 0.0, 0.8);
	glLineWidth(3);
      }
      else {
	glColor3d(0.8, 0.0, 0.0);
	glLineWidth(2);
      }

      for (j = 0; j < (int) strlen(stroke[0]); j++) {
	glutStrokeCharacter(GLUT_STROKE_ROMAN, stroke[0][j]);
      }
      glColor3d(0.0, 0.0, 0.0);
      glPushMatrix();
      /* Coordinate of tetrahedra {(x, ., ., .) */
      glTranslatef(200.0 * eps,
		   -2.0 * (i + 1.0) / (step + 1.0) + 1600.0 * eps, 0.0);
      glScalef(0.5, 0.5, 0.5);
      getvertices_tet3d(th, i + move_t->page, v);
      sprintf(stroke2[0], "{(%.2e, %.2e, %.2e),", th_orig->x[v[0]][0],
	      th_orig->x[v[0]][1], th_orig->x[v[0]][2]);
      glLineWidth(2);

      for (j = 0; j < (int) strlen(stroke2[0]); j++) {
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, stroke2[0][j]);
      }

      glPopMatrix();

      glPushMatrix();
      /* Coordinate of tetrahedra (., x, ., .) */
      glTranslatef(200.0 * eps, -2.0 * (i + 1.0) / (step + 1.0) + 800.0 * eps,
		   0.0);
      glScalef(0.5, 0.5, 0.5);
      sprintf(stroke2[1], "(%.2e, %.2e, %.2e),", th_orig->x[v[1]][0],
	      th_orig->x[v[1]][1], th_orig->x[v[1]][2]);
      glLineWidth(2);

      for (j = 0; j < (int) strlen(stroke2[1]); j++) {
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, stroke2[1][j]);
      }

      glPopMatrix();

      glPushMatrix();
      /* Coordinate of tetrahedra ( ., ., x, .) */
      glTranslatef(200.0 * eps, -2.0 * (i + 1.0) / (step + 1.0) + eps, 0.0);
      glScalef(0.5, 0.5, 0.5);
      sprintf(stroke2[2], "(%.2e, %.2e, %.2e),", th_orig->x[v[2]][0],
	      th_orig->x[v[2]][1], th_orig->x[v[2]][2]);
      glLineWidth(2);

      for (j = 0; j < (int) strlen(stroke2[2]); j++) {
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, stroke2[2][j]);
      }

      glPopMatrix();

      glPushMatrix();
      /* Coordinate of tetrahedra ( ., ., ., x)} */
      glTranslatef(200.0 * eps, -2.0 * (i + 1.0) / (step + 1.0) - 800.0 * eps,
		   0.0);
      glScalef(0.5, 0.5, 0.5);
      sprintf(stroke2[3], "(%.2e, %.2e, %.2e)}", th_orig->x[v[3]][0],
	      th_orig->x[v[3]][1], th_orig->x[v[3]][2]);
      glLineWidth(2);

      for (j = 0; j < (int) strlen(stroke2[3]); j++) {
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, stroke2[3][j]);
      }

      glPopMatrix();
    }
    glPopMatrix();

  }

  glLineWidth(1);

  glBegin(GL_LINE_STRIP);
  glColor3d(0.1, 0.3, 0.3);
  glVertex2d(0.0, -2.0 * 4 / (step + 1));
  glVertex2d(1.0, -2.0 * 4 / (step + 1));
  glVertex2d(1.0, -2.0);
  glVertex2d(0.0, -2.0);
  glVertex2d(0.0, -2.0 * 4 / (step + 1));
  glEnd();

  glPushMatrix();
  glTranslatef(0.0 + eps, -2.0 * (4 + 1.0) / (step + 1.0) + eps, 0.0);
  glScalef(0.6 * scale * 1.0 / 152.0, 0.6 * scale * 1.0 / 152.0,
	   0.6 * scale * 1.0 / 152.0);
  sprintf(stroke[0], "tetrahedra %d", number);
  glColor3d(0.0, 0.0, 0.0);
  glLineWidth(2);

  for (j = 0; j < (int) strlen(stroke[0]); j++) {
    glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, stroke[0][j]);
  }
  glPopMatrix();

  glLineWidth(1);

  glBegin(GL_LINE_STRIP);
  glColor3d(0.1, 0.3, 0.3);
  glVertex2d(1.0, -2.0 * 4 / (step + 1));
  glVertex2d(2.0, -2.0 * 4 / (step + 1));
  glVertex2d(2.0, -2.0);
  glVertex2d(1.0, -2.0);
  glEnd();

  glTranslatef(1.0, -2.0 * (4 + 1.0) / (step + 1.0), 0.0);

  glColor3d(0.0, 1.0, 0.0);

  glBegin(GL_TRIANGLES);
  glVertex3d(0.0, 1.0 / 6.0, 0.0);
  glVertex3d(0.5 - eps * 0.1, 0.0, 0.0);
  glVertex3d(0.5 - eps * 0.1, 1.0 / 3.0, 0.0);
  glEnd();

  glBegin(GL_TRIANGLES);
  glVertex3d(0.5 + eps * 0.1, 0.0, 0.0);
  glVertex3d(1.0, 1.0 / 6.0, 0.0);
  glVertex3d(0.5 + eps * 0.1, 1.0 / 3.0, 0.0);
  glEnd();

  glColor3d(1.0, 1.0, 1.0);
  glBegin(GL_LINE_STRIP);
  glVertex3d(0.0, 1.0 / 6.0, 0.0);
  glVertex3d(0.5 - eps * 0.1, 0.0, 0.0);
  glVertex3d(0.5 - eps * 0.1, 1.0 / 3.0, 0.0);
  glEnd();

  glBegin(GL_LINE_STRIP);
  glVertex3d(0.5 + eps * 0.1, 0.0, 0.0);
  glVertex3d(1.0, 1.0 / 6.0, 0.0);
  glVertex3d(0.5 + eps * 0.1, 1.0 / 3.0, 0.0);
  glEnd();

  glFlush();
  glutSwapBuffers();

  if (take_a_shot[0] == true) {

    static GLubyte *pixels;	/* Pointer for pixels */
    uint      paper_width = glutGet(GLUT_WINDOW_WIDTH);
    uint      paper_height = glutGet(GLUT_WINDOW_HEIGHT);

    pixels = allocmem(format_bytes * paper_height * paper_width);	/* Storage for the pixel */
    /* Reads the current buffer and save the pixels. RGBA has unsigned int type. */
    glReadPixels(0, 0, paper_width, paper_height, format, GL_UNSIGNED_BYTE,
		 pixels);

#ifdef USE_LIBPNG
    png_byte *png_bytes = NULL;	/* Memory for the png */
    png_byte **png_rows = NULL;	/* Array of memory for ever row */

    save_current_content_png("snapshot", paper_width, paper_height, pixels,
			     &png_bytes, &png_rows);
    free(png_bytes);
    free(png_rows);
#else
    save_current_content_ppm("snapshot", paper_width, paper_height, 255,
			     pixels);
#endif
    snapshots++;

    free(pixels);		/* free memory */
    take_a_shot[0] = false;	/* reset start values */
  }
}


  /* Drawing complete three dimensional triangulation */
static void
display_tet3d()
{

  uint      k, i;
  uint      v[4];

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glPushMatrix();
  system_transformation();
  glTranslated(move_1->rightleft, move_1->updown, 0.0);
  glRotated(move_1->angle_x, 1.0, 0.0, 0.0);
  glRotated(move_1->angle_y, 0.0, 1.0, 0.0);
  if (coord == 1) {
    glPushMatrix();
    coordinate();
    glPopMatrix();
  }
  glScaled(move_1->zoom, move_1->zoom, move_1->zoom);
  if (origin == true) {
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPointSize(6.0);
    glColor4d(0.8, 0.8, 0.8, 0.8);
    glBegin(GL_POINTS);
    glVertex3d(0.0, 0.0, 0.0);
    glEnd();
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
  }

  switch (mode) {		/* Drawing edges */
  case 0:

    glColor3d(1.0, 0.8, 0.0);
    glBegin(GL_LINES);

    for (k = 0; k < th->edges; k++) {
      glVertex3d(th->x[th->e[k][0]][0], th->x[th->e[k][0]][1],
		 th->x[th->e[k][0]][2]);
      glVertex3d(th->x[th->e[k][1]][0], th->x[th->e[k][1]][1],
		 th->x[th->e[k][1]][2]);
    }
    glEnd();

    break;

  case 1:

    for (k = 0; k < th->tetrahedra; k++) {
      getvertices_tet3d(th, k, v);
      glColor3d(0.8, 0.8, 0.8);
      glBegin(GL_LINES);
      for (i = 1; i < 4; i++) {
	glVertex3d(th->x[v[0]][0], th->x[v[0]][1], th->x[v[0]][2]);
	glVertex3d(th->x[v[i]][0], th->x[v[i]][1], th->x[v[i]][2]);
      }
      for (i = 2; i <= 3; i++) {
	glVertex3d(th->x[v[1]][0], th->x[v[1]][1], th->x[v[1]][2]);
	glVertex3d(th->x[v[i]][0], th->x[v[i]][1], th->x[v[i]][2]);
      }
      glVertex3d(th->x[v[3]][0], th->x[v[3]][1], th->x[v[3]][2]);
      glVertex3d(th->x[v[2]][0], th->x[v[2]][1], th->x[v[2]][2]);
      glEnd();
    }

    for (k = 0; k < th->tetrahedra; k++) {
      getvertices_tet3d(th, k, v);
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glColor4d(0.0, 0.0, 0.8, 0.5);
      for (i = 0; i < 4; i++) {
	glBegin(GL_TRIANGLES);
	glVertex3d(th->x[v[i]][0], th->x[v[i]][1], th->x[v[i]][2]);
	glVertex3d(th->x[v[(i + 1) % 4]][0], th->x[v[(i + 1) % 4]][1],
		   th->x[v[(i + 1) % 4]][2]);
	glVertex3d(th->x[v[(i + 2) % 4]][0], th->x[v[(i + 2) % 4]][1],
		   th->x[v[(i + 2) % 4]][2]);
	glEnd();
      }
      glDisable(GL_BLEND);
    }

    break;

  default:
    break;

  }

  if (points == 1) {
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPointSize(6.0);
    glColor3d(0.8, 0.8, 0.8);
    glBegin(GL_POINTS);
    for (k = 0; k < th->vertices; k++) {
      glVertex3d(th->x[k][0], th->x[k][1], th->x[k][2]);
    }
    glEnd();
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
  }
  glPopMatrix();

  glFlush();
  glutSwapBuffers();

  if (take_a_shot[1] == true) {

    static GLubyte *pixels;	/* Pointer for pixels */
    uint      paper_width = glutGet(GLUT_WINDOW_WIDTH);
    uint      paper_height = glutGet(GLUT_WINDOW_HEIGHT);

    pixels = allocmem(format_bytes * paper_height * paper_width);	/* Storage for the pixel */
    /* Reads the current buffer and save the pixels. RGBA has unsigned int type. */
    glReadPixels(0, 0, paper_width, paper_height, format, GL_UNSIGNED_BYTE,
		 pixels);

#ifdef USE_LIBPNG
    png_byte *png_bytes = NULL;	/* Memory for the png */
    png_byte **png_rows = NULL;	/* Array of memory for ever row */

    save_current_content_png("snapshot", paper_width, paper_height, pixels,
			     &png_bytes, &png_rows);
    free(png_bytes);
    free(png_rows);
#else
    save_current_content_ppm("snapshot", paper_width, paper_height, 255,
			     pixels);
#endif
    snapshots++;

    free(pixels);		/* free memory */
    take_a_shot[1] = false;	/* reset start values */
  }

}

  /* Drawing complete surface triangulation */
static void
display_surface()
{

  uint      k;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glPushMatrix();
  system_transformation();
  glTranslated(move_1->rightleft, move_1->updown, 0.0);
  glRotated(move_1->angle_x, 1.0, 0.0, 0.0);
  glRotated(move_1->angle_y, 0.0, 1.0, 0.0);
  if (coord == 1) {
    glPushMatrix();
    coordinate();
    glPopMatrix();
  }
  glScaled(move_1->zoom, move_1->zoom, move_1->zoom);
  if (origin == true) {
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPointSize(6.0);
    glColor4d(0.8, 0.8, 0.8, 0.8);
    glBegin(GL_POINTS);
    glVertex3d(0.0, 0.0, 0.0);
    glEnd();
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
  }

  switch (mode) {		/* Drawing edges */
  case 0:

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4d(0.6, 0.6, 0.2, 0.7);
    glBegin(GL_LINES);
    for (k = 0; k < gr->edges; k++) {
      glVertex3d(gr->x[gr->e[k][0]][0], gr->x[gr->e[k][0]][1],
		 gr->x[gr->e[k][0]][2]);
      glVertex3d(gr->x[gr->e[k][1]][0], gr->x[gr->e[k][1]][1],
		 gr->x[gr->e[k][1]][2]);
    }
    glEnd();
    glDisable(GL_BLEND);

    break;

  case 1:

    glColor3d(0.8, 0.8, 0.8);
    glBegin(GL_LINES);
    for (k = 0; k < gr->edges; k++) {
      glVertex3d(gr->x[gr->e[k][0]][0], gr->x[gr->e[k][0]][1],
		 gr->x[gr->e[k][0]][2]);
      glVertex3d(gr->x[gr->e[k][1]][0], gr->x[gr->e[k][1]][1],
		 gr->x[gr->e[k][1]][2]);
    }
    glEnd();

    glBegin(GL_TRIANGLES);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4d(0.0, 0.0, 0.8, 0.5);
    for (k = 0; k < gr->triangles; k++) {
      glVertex3d(gr->x[gr->t[k][0]][0], gr->x[gr->t[k][0]][1],
		 gr->x[gr->t[k][0]][2]);
      glVertex3d(gr->x[gr->t[k][1]][0], gr->x[gr->t[k][1]][1],
		 gr->x[gr->t[k][1]][2]);
      glVertex3d(gr->x[gr->t[k][2]][0], gr->x[gr->t[k][2]][1],
		 gr->x[gr->t[k][2]][2]);
      glNormal3d(gr->n[k][0], gr->n[k][1], gr->n[k][2]);

    }
    glEnd();
    glDisable(GL_BLEND);

    break;

  default:
    break;

  }

  if (points == 1) {
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPointSize(6.0);
    glColor3d(0.8, 0.8, 0.8);
    glBegin(GL_POINTS);
    for (k = 0; k < gr->vertices; k++) {
      glVertex3d(gr->x[k][0], gr->x[k][1], gr->x[k][2]);
    }
    glEnd();
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
  }
  glPopMatrix();

  glDisable(GL_POLYGON_OFFSET_FILL);
  glFlush();
  glutSwapBuffers();

  if (take_a_shot[1] == true) {

    static GLubyte *pixels;	/* Pointer for pixels */
    uint      paper_width = glutGet(GLUT_WINDOW_WIDTH);
    uint      paper_height = glutGet(GLUT_WINDOW_HEIGHT);

    pixels = allocmem(format_bytes * paper_height * paper_width);	/* Storage for the pixel */
    /* Reads the current buffer and save the pixels. RGBA has unsigned int type. */
    glReadPixels(0, 0, paper_width, paper_height, format, GL_UNSIGNED_BYTE,
		 pixels);

#ifdef USE_LIBPNG
    png_byte *png_bytes = NULL;	/* Memory for the png */
    png_byte **png_rows = NULL;	/* Array of memory for ever row */

    save_current_content_png("snapshot", paper_width, paper_height, pixels,
			     &png_bytes, &png_rows);
    free(png_bytes);
    free(png_rows);
#else
    save_current_content_ppm("snapshot", paper_width, paper_height, 255,
			     pixels);
#endif
    snapshots++;

    free(pixels);		/* free memory */
    take_a_shot[1] = false;	/* reset start values */
  }

}

static void
display_certain_surface_element()
{

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glPushMatrix();
  system_transformation();
  glTranslated(move_1->rightleft, move_1->updown, 0.0);
  glRotated(move_1->angle_x, 1.0, 0.0, 0.0);
  glRotated(move_1->angle_y, 0.0, 1.0, 0.0);
  if (coord == 1) {
    glPushMatrix();
    coordinate();
    glPopMatrix();
  }
  glScaled(move_1->zoom, move_1->zoom, move_1->zoom);
  if (origin == true) {
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPointSize(6.0);
    glColor4d(0.8, 0.8, 0.8, 0.8);
    glBegin(GL_POINTS);
    glVertex3d(0.0, 0.0, 0.0);
    glEnd();
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
  }

  if (mode == 1) {

    glColor3d(0.8, 0.8, 0.8);
    glBegin(GL_LINE_LOOP);
    glVertex3d(gr->x[gr->t[move_t->target][0]][0],
	       gr->x[gr->t[move_t->target][0]][1],
	       gr->x[gr->t[move_t->target][0]][2]);
    glVertex3d(gr->x[gr->t[move_t->target][1]][0],
	       gr->x[gr->t[move_t->target][1]][1],
	       gr->x[gr->t[move_t->target][1]][2]);
    glVertex3d(gr->x[gr->t[move_t->target][2]][0],
	       gr->x[gr->t[move_t->target][2]][1],
	       gr->x[gr->t[move_t->target][2]][2]);
    glEnd();
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4d(0.0, 0.2, 0.8, 0.3);

    glBegin(GL_TRIANGLES);
    glVertex3d(gr->x[gr->t[move_t->target][0]][0],
	       gr->x[gr->t[move_t->target][0]][1],
	       gr->x[gr->t[move_t->target][0]][2]);
    glVertex3d(gr->x[gr->t[move_t->target][1]][0],
	       gr->x[gr->t[move_t->target][1]][1],
	       gr->x[gr->t[move_t->target][1]][2]);
    glVertex3d(gr->x[gr->t[move_t->target][2]][0],
	       gr->x[gr->t[move_t->target][2]][1],
	       gr->x[gr->t[move_t->target][2]][2]);
    glEnd();
    glDisable(GL_BLEND);
  }

  outline(2);

  glPopMatrix();

  glFlush();
  glutSwapBuffers();

  if (take_a_shot[1] == true) {

    static GLubyte *pixels;	/* Pointer for pixels */
    uint      paper_width = glutGet(GLUT_WINDOW_WIDTH);
    uint      paper_height = glutGet(GLUT_WINDOW_HEIGHT);

    pixels = allocmem(format_bytes * paper_height * paper_width);	/* Storage for the pixel */
    /* Reads the current buffer and save the pixels. RGBA has unsigned int type. */
    glReadPixels(0, 0, paper_width, paper_height, format, GL_UNSIGNED_BYTE,
		 pixels);

#ifdef USE_LIBPNG
    png_byte *png_bytes = NULL;	/* Memory for the png */
    png_byte **png_rows = NULL;	/* Array of memory for ever row */

    save_current_content_png("snapshot", paper_width, paper_height, pixels,
			     &png_bytes, &png_rows);
    free(png_bytes);
    free(png_rows);
#else
    save_current_content_ppm("snapshot", paper_width, paper_height, 255,
			     pixels);
#endif
    snapshots++;

    free(pixels);		/* free memory */
    take_a_shot[1] = false;	/* reset start values */
  }

}

  /* Drawing table of a surface triangulation */
static void
display_table_sur()
{

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  uint      step = 4;

  uint      i, j;
  char      stroke[1][20];
  char      stroke2[4][20];
  double    scale = 0.75 / (step + 1);
  double    eps = (2.0 / (step + 1) - scale) / 2.0;
  uint      number = gr_orig->triangles;

  glClearColor(1.0, 1.0, 1.0, 1.0);
  system_transformation();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslated(-1.0, 1.0, 0.0);

  /* Drawing left table with number of triangle and coordinates */
  for (i = 0; i < step; i++) {
    glBegin(GL_LINE_STRIP);
    glColor3d(0.1, 0.3, 0.3);
    glVertex2d(0.0, -2.0 * i / (step + 1));
    glVertex2d(1.0, -2.0 * i / (step + 1));
    glVertex2d(1.0, -2.0 * (i + 1) / (step + 1));
    glVertex2d(0.0, -2.0 * (i + 1) / (step + 1));
    glVertex2d(0.0, -2.0 * i / (step + 1));
    glEnd();

    glPushMatrix();
    if ((i + move_t->page * 8) < number) {
      /* Number of triangle */
      glTranslatef(0.0, -2.0 * (i + 1.0) / (step + 1.0) + eps, 0.0);
      glScalef(scale * 1.0 / 152.0, scale * 1.0 / 152.0, scale * 1.0 / 152.0);
      glScalef(0.5, 0.5, 0.5);
      sprintf(stroke[0], "%d", (i + move_t->page * 8));
      if (move_t->target == (i + move_t->page * 8)) {
	glColor3d(0.0, 0.0, 0.8);
	glLineWidth(3);
      }
      else {
	glColor3d(0.8, 0.0, 0.0);
	glLineWidth(2);
      }
      for (j = 0; j < (int) strlen(stroke[0]); j++) {
	glutStrokeCharacter(GLUT_STROKE_ROMAN, stroke[0][j]);
      }
      glColor3d(0.0, 0.0, 0.0);
      glPushMatrix();
      /* Coordinate of triangle {(x, ., .) */
      glTranslatef(200.0 * eps,
		   -2.0 * (i + 1.0) / (step + 1.0) + 1500.0 * eps, 0.0);
      glScalef(0.5, 0.5, 0.5);
      sprintf(stroke2[0], "{(%.2e, %.2e, %.2e),",
	      gr_orig->x[gr_orig->t[i + move_t->page * 8][0]][0],
	      gr_orig->x[gr_orig->t[i + move_t->page * 8][0]][1],
	      gr_orig->x[gr_orig->t[i + move_t->page * 8][0]][2]);
      glLineWidth(2);

      for (j = 0; j < (int) strlen(stroke2[0]); j++) {
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, stroke2[0][j]);
      }

      glPopMatrix();

      glPushMatrix();
      /* Coordinate of triangle (., x, .) */
      glTranslatef(200.0 * eps, -2.0 * (i + 1.0) / (step + 1.0) + 500.0 * eps,
		   0.0);
      glScalef(0.5, 0.5, 0.5);
      sprintf(stroke2[1], "(%.2e, %.2e, %.2e),",
	      gr_orig->x[gr_orig->t[i + move_t->page * 8][1]][0],
	      gr_orig->x[gr_orig->t[i + move_t->page * 8][1]][1],
	      gr_orig->x[gr_orig->t[i + move_t->page * 8][1]][2]);
      glLineWidth(2);

      for (j = 0; j < (int) strlen(stroke2[1]); j++) {
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, stroke2[1][j]);
      }

      glPopMatrix();

      glPushMatrix();
      /* Coordinate of triangle ( ., ., x)} */
      glTranslatef(200.0 * eps, -2.0 * (i + 1.0) / (step + 1.0) - 500.0 * eps,
		   0.0);
      glScalef(0.5, 0.5, 0.5);
      sprintf(stroke2[2], "(%.2e, %.2e, %.2e)}",
	      gr_orig->x[gr_orig->t[i + move_t->page * 8][2]][0],
	      gr_orig->x[gr_orig->t[i + move_t->page * 8][2]][1],
	      gr_orig->x[gr_orig->t[i + move_t->page * 8][2]][2]);
      glLineWidth(2);

      for (j = 0; j < (int) strlen(stroke2[2]); j++) {
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, stroke2[2][j]);
      }

      glPopMatrix();

    }
    glPopMatrix();

    /* Drawing right table with number of triangle and coordinates */
    glBegin(GL_LINE_STRIP);
    glColor3d(0.1, 0.3, 0.3);
    glVertex2d(1.0, -2.0 * i / (step + 1));
    glVertex2d(2.0, -2.0 * i / (step + 1));
    glVertex2d(2.0, -2.0 * (i + 1) / (step + 1));
    glVertex2d(1.0, -2.0 * (i + 1) / (step + 1));
    glEnd();

    glPushMatrix();
    if ((i + 4 + move_t->page * 8) < number) {
      /* Number of triangle */
      glTranslatef(1.0, -2.0 * (i + 1.0) / (step + 1.0) + eps, 0.0);
      glScalef(scale * 1.0 / 152.0, scale * 1.0 / 152.0, scale * 1.0 / 152.0);
      glScalef(0.5, 0.5, 0.5);
      sprintf(stroke[0], "%d", (i + 4 + move_t->page * 8));
      if (move_t->target == (i + 4 + move_t->page * 8)) {
	glColor3d(0.0, 0.0, 0.8);
	glLineWidth(3);
      }
      else {
	glColor3d(0.8, 0.0, 0.0);
	glLineWidth(2);
      }
      for (j = 0; j < (int) strlen(stroke[0]); j++) {
	glutStrokeCharacter(GLUT_STROKE_ROMAN, stroke[0][j]);
      }
      glColor3d(0.0, 0.0, 0.0);
      glPushMatrix();
      /* Coordinate of triangle {(x, ., .) */
      glTranslatef(200.0 * eps,
		   -2.0 * (i + 1.0) / (step + 1.0) + 1500.0 * eps, 0.0);
      glScalef(0.5, 0.5, 0.5);
      sprintf(stroke2[0], "{(%.2e, %.2e, %.2e),",
	      gr_orig->x[gr_orig->t[i + 4 + move_t->page * 8][0]][0],
	      gr_orig->x[gr_orig->t[i + 4 + move_t->page * 8][0]][1],
	      gr_orig->x[gr_orig->t[i + 4 + move_t->page * 8][0]][2]);
      glLineWidth(2);

      for (j = 0; j < (int) strlen(stroke2[0]); j++) {
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, stroke2[0][j]);
      }

      glPopMatrix();

      glPushMatrix();
      /* Coordinate of triangle (., x, .) */
      glTranslatef(200.0 * eps, -2.0 * (i + 1.0) / (step + 1.0) + 500.0 * eps,
		   0.0);
      glScalef(0.5, 0.5, 0.5);
      sprintf(stroke2[1], "(%.2e, %.2e, %.2e),",
	      gr_orig->x[gr_orig->t[i + 4 + move_t->page * 8][1]][0],
	      gr_orig->x[gr_orig->t[i + 4 + move_t->page * 8][1]][1],
	      gr_orig->x[gr_orig->t[i + 4 + move_t->page * 8][1]][2]);
      glLineWidth(2);

      for (j = 0; j < (int) strlen(stroke2[1]); j++) {
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, stroke2[1][j]);
      }

      glPopMatrix();

      glPushMatrix();
      /* Coordinate of triangle ( ., ., x)} */
      glTranslatef(200.0 * eps, -2.0 * (i + 1.0) / (step + 1.0) - 500.0 * eps,
		   0.0);
      glScalef(0.5, 0.5, 0.5);
      sprintf(stroke2[2], "(%.2e, %.2e, %.2e)}",
	      gr_orig->x[gr_orig->t[i + 4 + move_t->page * 8][2]][0],
	      gr_orig->x[gr_orig->t[i + 4 + move_t->page * 8][2]][1],
	      gr_orig->x[gr_orig->t[i + 4 + move_t->page * 8][2]][2]);
      glLineWidth(2);

      for (j = 0; j < (int) strlen(stroke2[2]); j++) {
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, stroke2[2][j]);
      }

      glPopMatrix();

    }
    glPopMatrix();
  }

  glLineWidth(1);

  glBegin(GL_LINE_STRIP);
  glColor3d(0.1, 0.3, 0.3);
  glVertex2d(0.0, -2.0 * 4 / (step + 1));
  glVertex2d(1.0, -2.0 * 4 / (step + 1));
  glVertex2d(1.0, -2.0);
  glVertex2d(0.0, -2.0);
  glVertex2d(0.0, -2.0 * 4 / (step + 1));

  glEnd();

  glPushMatrix();
  glTranslatef(0.0 + eps, -2.0 * (4 + 1.0) / (step + 1.0) + eps, 0.0);
  glScalef(0.6 * scale * 1.0 / 152.0, 0.6 * scale * 1.0 / 152.0,
	   0.6 * scale * 1.0 / 152.0);
  sprintf(stroke[0], "triangle %d", number);
  glColor3d(0.0, 0.0, 0.0);
  glLineWidth(2);

  for (j = 0; j < (int) strlen(stroke[0]); j++) {
    glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, stroke[0][j]);
  }
  glPopMatrix();

  glLineWidth(1);

  glBegin(GL_LINE_STRIP);
  glColor3d(0.1, 0.3, 0.3);
  glVertex2d(1.0, -2.0 * 4 / (step + 1));
  glVertex2d(2.0, -2.0 * 4 / (step + 1));
  glVertex2d(2.0, -2.0);
  glVertex2d(1.0, -2.0);
  glEnd();

  glTranslatef(1.0, -2.0 * (4 + 1.0) / (step + 1.0), 0.0);

  glColor3d(0.0, 1.0, 0.0);

  glBegin(GL_TRIANGLES);
  glVertex3d(0.0, 1.0 / 6.0, 0.0);
  glVertex3d(0.5 - eps * 0.1, 0.0, 0.0);
  glVertex3d(0.5 - eps * 0.1, 1.0 / 3.0, 0.0);
  glEnd();

  glBegin(GL_TRIANGLES);
  glVertex3d(0.5 + eps * 0.1, 0.0, 0.0);
  glVertex3d(1.0, 1.0 / 6.0, 0.0);
  glVertex3d(0.5 + eps * 0.1, 1.0 / 3.0, 0.0);
  glEnd();

  glColor3d(1.0, 1.0, 1.0);
  glBegin(GL_LINE_STRIP);
  glVertex3d(0.0, 1.0 / 6.0, 0.0);
  glVertex3d(0.5 - eps * 0.1, 0.0, 0.0);
  glVertex3d(0.5 - eps * 0.1, 1.0 / 3.0, 0.0);
  glEnd();

  glBegin(GL_LINE_STRIP);
  glVertex3d(0.5 + eps * 0.1, 0.0, 0.0);
  glVertex3d(1.0, 1.0 / 6.0, 0.0);
  glVertex3d(0.5 + eps * 0.1, 1.0 / 3.0, 0.0);
  glEnd();

  glFlush();
  glutSwapBuffers();

  if (take_a_shot[0] == true) {

    static GLubyte *pixels;	/* Pointer for pixels */
    uint      paper_width = glutGet(GLUT_WINDOW_WIDTH);
    uint      paper_height = glutGet(GLUT_WINDOW_HEIGHT);

    pixels = allocmem(format_bytes * paper_height * paper_width);	/* Storage for the pixel */
    /* Reads the current buffer and save the pixels. RGBA has unsigned int type. */
    glReadPixels(0, 0, paper_width, paper_height, format, GL_UNSIGNED_BYTE,
		 pixels);

#ifdef USE_LIBPNG
    png_byte *png_bytes = NULL;	/* Memory for the png */
    png_byte **png_rows = NULL;	/* Array of memory for ever row */

    save_current_content_png("snapshot", paper_width, paper_height, pixels,
			     &png_bytes, &png_rows);
    free(png_bytes);
    free(png_rows);
#else
    save_current_content_ppm("snapshot", paper_width, paper_height, 255,
			     pixels);
#endif
    snapshots++;

    free(pixels);		/* free memory */
    take_a_shot[0] = false;	/* reset start values */
  }

}



  /* Display solution on surface */
static void
display_surface_value_triangle()
{

  uint      i;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glPushMatrix();
  system_transformation();
  glTranslated(move_1->rightleft, move_1->updown, 0.0);
  glRotated(move_1->angle_x, 1.0, 0.0, 0.0);
  glRotated(move_1->angle_y, 0.0, 1.0, 0.0);
  if (coord == 1) {
    glPushMatrix();
    coordinate();
    glPopMatrix();
  }

  glScaled(move_1->zoom, move_1->zoom, move_1->zoom);
  if (origin == true) {
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPointSize(6.0);
    glColor4d(0.8, 0.8, 0.8, 0.8);
    glBegin(GL_POINTS);
    glVertex3d(0.0, 0.0, 0.0);
    glEnd();
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
  }

  for (i = 0; i < gr->triangles; i++) {

    color_triangle(sol->values->v[i]);

    glBegin(GL_TRIANGLES);
    glVertex3d(gr->x[gr->t[i][0]][0], gr->x[gr->t[i][0]][1],
	       gr->x[gr->t[i][0]][2]);
    glVertex3d(gr->x[gr->t[i][1]][0], gr->x[gr->t[i][1]][1],
	       gr->x[gr->t[i][1]][2]);
    glVertex3d(gr->x[gr->t[i][2]][0], gr->x[gr->t[i][2]][1],
	       gr->x[gr->t[i][2]][2]);
    glEnd();
  }


  glPopMatrix();
  glFlush();
  glutSwapBuffers();

  if (take_a_shot[1] == true) {

    static GLubyte *pixels;	/* Pointer for pixels */
    uint      paper_width = glutGet(GLUT_WINDOW_WIDTH);
    uint      paper_height = glutGet(GLUT_WINDOW_HEIGHT);

    pixels = allocmem(format_bytes * paper_height * paper_width);	/* Storage for the pixel */
    /* Reads the current buffer and save the pixels. RGBA has unsigned int type. */
    glReadPixels(0, 0, paper_width, paper_height, format, GL_UNSIGNED_BYTE,
		 pixels);

#ifdef USE_LIBPNG
    png_byte *png_bytes = NULL;	/* Memory for the png */
    png_byte **png_rows = NULL;	/* Array of memory for ever row */

    save_current_content_png("snapshot", paper_width, paper_height, pixels,
			     &png_bytes, &png_rows);
    free(png_bytes);
    free(png_rows);
#else
    save_current_content_ppm("snapshot", paper_width, paper_height, 255,
			     pixels);
#endif
    snapshots++;

    free(pixels);		/* free memory */
    take_a_shot[1] = false;	/* reset start values */
  }

}

  /* Display solution on surface */
static void
display_surface_value_vertices()
{

  uint      i;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glPushMatrix();
  system_transformation();
  glTranslated(move_1->rightleft, move_1->updown, 0.0);
  glRotated(move_1->angle_x, 1.0, 0.0, 0.0);
  glRotated(move_1->angle_y, 0.0, 1.0, 0.0);
  if (coord == 1) {
    glPushMatrix();
    coordinate();
    glPopMatrix();
  }

  glScaled(move_1->zoom, move_1->zoom, move_1->zoom);
  if (origin == true) {
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPointSize(6.0);
    glColor4d(0.8, 0.8, 0.8, 0.8);
    glBegin(GL_POINTS);
    glVertex3d(0.0, 0.0, 0.0);
    glEnd();
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
  }

  for (i = 0; i < gr->triangles; i++) {

    glBegin(GL_TRIANGLES);
    color_triangle(sol->values->v[gr->t[i][0]]);
    glVertex3d(gr->x[gr->t[i][0]][0], gr->x[gr->t[i][0]][1],
	       gr->x[gr->t[i][0]][2]);
    color_triangle(sol->values->v[gr->t[i][1]]);
    glVertex3d(gr->x[gr->t[i][1]][0], gr->x[gr->t[i][1]][1],
	       gr->x[gr->t[i][1]][2]);
    color_triangle(sol->values->v[gr->t[i][2]]);
    glVertex3d(gr->x[gr->t[i][2]][0], gr->x[gr->t[i][2]][1],
	       gr->x[gr->t[i][2]][2]);
    glEnd();
  }


  glPopMatrix();
  glFlush();
  glutSwapBuffers();

  if (take_a_shot[1] == true) {

    static GLubyte *pixels;	/* Pointer for pixels */
    uint      paper_width = glutGet(GLUT_WINDOW_WIDTH);
    uint      paper_height = glutGet(GLUT_WINDOW_HEIGHT);

    pixels = allocmem(format_bytes * paper_height * paper_width);	/* Storage for the pixel */
    /* Reads the current buffer and save the pixels. RGBA has unsigned int type. */
    glReadPixels(0, 0, paper_width, paper_height, format, GL_UNSIGNED_BYTE,
		 pixels);

#ifdef USE_LIBPNG
    png_byte *png_bytes = NULL;	/* Memory for the png */
    png_byte **png_rows = NULL;	/* Array of memory for ever row */

    save_current_content_png("snapshot", paper_width, paper_height, pixels,
			     &png_bytes, &png_rows);
    free(png_bytes);
    free(png_rows);
#else
    save_current_content_ppm("snapshot", paper_width, paper_height, 255,
			     pixels);
#endif
    snapshots++;

    free(pixels);		/* free memory */
    take_a_shot[1] = false;	/* reset start values */
  }

}

   /* Display solution on surface */
static void
display_surface_value_edges()
{

  uint      i;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glPushMatrix();
  system_transformation();
  glTranslated(move_1->rightleft, move_1->updown, 0.0);
  glRotated(move_1->angle_x, 1.0, 0.0, 0.0);
  glRotated(move_1->angle_y, 0.0, 1.0, 0.0);
  if (coord == 1) {
    glPushMatrix();
    coordinate();
    glPopMatrix();
  }

  glScaled(move_1->zoom, move_1->zoom, move_1->zoom);
  if (origin == true) {
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPointSize(6.0);
    glColor4d(0.8, 0.8, 0.8, 0.8);
    glBegin(GL_POINTS);
    glVertex3d(0.0, 0.0, 0.0);
    glEnd();
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
  }

  for (i = 0; i < gr->edges; i++) {

    glBegin(GL_LINES);
    color_triangle(sol->values->v[gr->e[i][0]]);
    glVertex3d(gr->x[gr->e[i][0]][0], gr->x[gr->e[i][0]][1],
	       gr->x[gr->e[i][0]][2]);
    color_triangle(sol->values->v[gr->e[i][1]]);
    glVertex3d(gr->x[gr->e[i][1]][0], gr->x[gr->e[i][1]][1],
	       gr->x[gr->e[i][1]][2]);
    glEnd();
  }


  glPopMatrix();

  glFlush();
  glutSwapBuffers();

  if (take_a_shot[1] == true) {

    static GLubyte *pixels;	/* Pointer for pixels */
    uint      paper_width = glutGet(GLUT_WINDOW_WIDTH);
    uint      paper_height = glutGet(GLUT_WINDOW_HEIGHT);

    pixels = allocmem(format_bytes * paper_height * paper_width);	/* Storage for the pixel */
    /* Reads the current buffer and save the pixels. RGBA has unsigned int type. */
    glReadPixels(0, 0, paper_width, paper_height, format, GL_UNSIGNED_BYTE,
		 pixels);

#ifdef USE_LIBPNG
    png_byte *png_bytes = NULL;	/* Memory for the png */
    png_byte **png_rows = NULL;	/* Array of memory for ever row */

    save_current_content_png("snapshot", paper_width, paper_height, pixels,
			     &png_bytes, &png_rows);
    free(png_bytes);
    free(png_rows);
#else
    save_current_content_ppm("snapshot", paper_width, paper_height, 255,
			     pixels);
#endif
    snapshots++;

    free(pixels);		/* free memory */
    take_a_shot[1] = false;	/* reset start values */
  }

}

  /* Display a legend for solution on surface */
static void
display_legend_surface()
{

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  uint      j;
  char      stroke[2][10];
  double    scale = 0.125;
  uint      current = sol->current_grid_number;
  uint      type = sol->type;

  system_transformation();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glPushMatrix();
  glTranslatef(-0.25, 0.7, 0.0);
  glScalef(scale * 1.0 / 152.0, scale * 1.0 / 152.0, scale * 1.0 / 152.0);
  sprintf(stroke[0], "max %.2e", sol->max[current][type]);
  glColor3d(1.0, 1.0, 1.0);
  glLineWidth(2);

  for (j = 0; j < (int) strlen(stroke[0]); j++) {
    glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, stroke[0][j]);
  }
  glPopMatrix();

  glPushMatrix();
  glTranslatef(-0.25, -0.8, 0.0);
  glScalef(scale * 1.0 / 152.0, scale * 1.0 / 152.0, scale * 1.0 / 152.0);
  sprintf(stroke[1], "min %.2e", sol->min[current][type]);
  glColor3d(1.0, 1.0, 1.0);
  glLineWidth(2);

  for (j = 0; j < (int) strlen(stroke[1]); j++) {
    glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, stroke[1][j]);
  }
  glPopMatrix();


  glBegin(GL_TRIANGLE_STRIP);
  glColor3d(1.0, 1.0, 0.0);
  glVertex2d(-0.5, 0.75);
  glColor3d(1.0, 1.0, 0.0);
  glVertex2d(-0.75, 0.75);
  glColor3d(0.0, 0.0, 1.0);
  glVertex2d(-0.5, 0.0);
  glColor3d(0.0, 0.0, 1.0);
  glVertex2d(-0.75, 0.0);
  glColor3d(1.0, 0.0, 0.0);
  glVertex2d(-0.5, -0.75);
  glColor3d(1.0, 0.0, 0.0);
  glVertex2d(-0.75, -0.75);
  glEnd();

  glFlush();
  glutSwapBuffers();

  if (take_a_shot[0] == true) {

    static GLubyte *pixels;	/* Pointer for pixels */
    uint      paper_width = glutGet(GLUT_WINDOW_WIDTH);
    uint      paper_height = glutGet(GLUT_WINDOW_HEIGHT);

    pixels = allocmem(format_bytes * paper_height * paper_width);	/* Storage for the pixel */
    /* Reads the current buffer and save the pixels. RGBA has unsigned int type. */
    glReadPixels(0, 0, paper_width, paper_height, format, GL_UNSIGNED_BYTE,
		 pixels);

#ifdef USE_LIBPNG
    png_byte *png_bytes = NULL;	/* Memory for the png */
    png_byte **png_rows = NULL;	/* Array of memory for ever row */

    save_current_content_png("snapshot", paper_width, paper_height, pixels,
			     &png_bytes, &png_rows);
    free(png_bytes);
    free(png_rows);
#else
    save_current_content_ppm("snapshot", paper_width, paper_height, 255,
			     pixels);
#endif
    snapshots++;

    free(pixels);		/* free memory */
    take_a_shot[0] = false;	/* reset start values */
  }

}


  /* Display surface and solution */
static void
display_solution_surface()
{

  uint      k;
  double    maxx, maxy, maxz, scale;
  longindex j;
  uint      current = sol->current_grid_number;
  uint      type = sol->type;
  pgrid     grid;
  longindex x = 0;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glPushMatrix();
  system_transformation();
  glTranslated(move_1->rightleft, move_1->updown, 0.0);
  glRotated(move_1->angle_x, 1.0, 0.0, 0.0);
  glRotated(move_1->angle_y, 0.0, 1.0, 0.0);
  if (coord == 1) {
    glPushMatrix();
    coordinate();
    glPopMatrix();
  }


  /* Find scaling for current grid */
  scale = 1.0;
  if (sol->grids > 0) {
    x = sol->gr[current]->x;
    maxy =
      ABS(sol->gr[current]->point[0][1] -
	  sol->gr[current]->point[(x - 1)][1]);
    maxx =
      ABS(sol->gr[current]->point[0][0] -
	  sol->gr[current]->point[x * (sol->gr[current]->y - 1)][0]);
    maxz = sol->max[current][type] - sol->min[current][type];

    scale = (maxx > maxy ? maxx : maxy);
    scale = (maxz > scale ? maxz : scale);

    scale *= 0.5;
    scale = (scale > 1.0 ? 1.0 / scale : 1.0);
  }

  glScaled(move_1->zoom * scale, move_1->zoom * scale, move_1->zoom * scale);
  if (origin == true) {
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPointSize(6.0);
    glColor4d(0.8, 0.8, 0.8, 0.8);
    glBegin(GL_POINTS);
    glVertex3d(0.0, 0.0, 0.0);
    glEnd();
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
  }

  /* model */
  if (sol->model == true) {
    glColor3d(0.6, 0.4, 0.1);
    glBegin(GL_TRIANGLES);
    for (k = 0; k < gr->triangles; k++) {
      glVertex3d(gr->x[gr->t[k][0]][0], gr->x[gr->t[k][0]][1],
		 gr->x[gr->t[k][0]][2]);
      glVertex3d(gr->x[gr->t[k][1]][0], gr->x[gr->t[k][1]][1],
		 gr->x[gr->t[k][1]][2]);
      glVertex3d(gr->x[gr->t[k][2]][0], gr->x[gr->t[k][2]][1],
		 gr->x[gr->t[k][2]][2]);
    }
    glEnd();
  }


  if (sol->plane == true) {
    if (sol->animate == true) {
      grid = sol->timegr[current];
    }
    else {
      grid = sol->gr[current];
    }
    j = 0;
    while (j < sol->gr[current]->y - 1) {
      glBegin(GL_TRIANGLE_STRIP);
      for (k = 0; k < x; k++) {
	color_triangle(grid->point[j * x + k][2 * (type + 1) + 1]);
	glVertex3d(sol->gr[current]->point[j * x + k][0],
		   sol->gr[current]->point[j * x + k][1],
		   grid->point[j * x + k][2 * (type + 1)]);
	color_triangle(grid->point[(j + 1) * x + k][2 * (type + 1) + 1]);
	glVertex3d(sol->gr[current]->point[(j + 1) * x + k][0],
		   sol->gr[current]->point[(j + 1) * x + k][1],
		   grid->point[(j + 1) * x + k][2 * (type + 1)]);
      }
      glEnd();
      j += 1;
    }
  }

  glPopMatrix();
  glFlush();
  glutSwapBuffers();

  if (take_a_shot[1] == true) {

    static GLubyte *pixels;	/* Pointer for pixels */
    uint      paper_width = glutGet(GLUT_WINDOW_WIDTH);
    uint      paper_height = glutGet(GLUT_WINDOW_HEIGHT);

    pixels = allocmem(format_bytes * paper_height * paper_width);	/* Storage for the pixel */
    /* Reads the current buffer and save the pixels. RGBA has unsigned int type. */
    glReadPixels(0, 0, paper_width, paper_height, format, GL_UNSIGNED_BYTE,
		 pixels);

#ifdef USE_LIBPNG
    png_byte *png_bytes = NULL;	/* Memory for the png */
    png_byte **png_rows = NULL;	/* Array of memory for ever row */

    save_current_content_png("snapshot", paper_width, paper_height, pixels,
			     &png_bytes, &png_rows);
    free(png_bytes);
    free(png_rows);
#else
    save_current_content_ppm("snapshot", paper_width, paper_height, 255,
			     pixels);
#endif
    snapshots++;

    free(pixels);		/* free memory */
    take_a_shot[1] = false;	/* reset start values */
  }

}


   /*
    *   GLUT reshape functions  *
    *                                             */

  /* Reshape of an orthogonal projection */
static void
reshape_ortho(int width, int height)
{

  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  /* Scaling */
  if (width > height) {
    glScaled(((double) height / width), 1.0, 1.0);
    move_o->orthoscaling = 2.0 / height;
    move_o->offset_x = (width - height) * 0.5;
    move_o->offset_y = 0;
  }
  else {
    glScaled(1.0, ((double) width / height), 1.0);
    move_o->orthoscaling = 2.0 / width;
    move_o->offset_x = 0;
    move_o->offset_y = (height - width) * 0.5;
  }
}

  /* Reshape of drawn directional cluster */
static void
reshape_dcluster(int width, int height)
{

  uint      w;
  w = glutGetWindow();
  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  /* Scaling */
  if (width > height) {
    glScaled(((double) height / width), 1.0, 1.0);
  }
  else {
    glScaled(1.0, ((double) width / height), 1.0);
  }

  if (clust2 != 0) {
    if (((clust1 == w) && t->dim == 3) || ((clust2 == w) && s->dim == 3)) {
      gluPerspective(35.0, 1.0, 1.0, 30.0);
      glTranslated(0.0, 0.0, -5.0);
    }
  }
  else {
    if (t->dim == 3 || s->dim == 3) {
      gluPerspective(35.0, 1.0, 1.0, 30.0);
      glTranslated(0.0, 0.0, -5.0);
    }
  }

  glEnable(GL_DEPTH_TEST);
}

  /* Reshape of drawn cluster */
static void
reshape_cluster(int width, int height)
{

  uint      w;
  w = glutGetWindow();
  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  /* Scaling */
  if (width > height) {
    glScaled(((double) height / width), 1.0, 1.0);
  }
  else {
    glScaled(1.0, ((double) width / height), 1.0);
  }

  if (clust2 != 0) {
    if (((clust1 == w) && ct->dim == 3) || ((clust2 == w) && cs->dim == 3)) {
      gluPerspective(35.0, 1.0, 1.0, 30.0);
      glTranslated(0.0, 0.0, -5.0);
    }
  }
  else {
    if (ct->dim == 3 || cs->dim == 3) {
      gluPerspective(35.0, 1.0, 1.0, 30.0);
      glTranslated(0.0, 0.0, -5.0);
    }
  }

  glEnable(GL_DEPTH_TEST);
}

  /* Reshape of a triangulation table */
static void
reshape_table(int width, int height)
{

  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  /* scaling */
  if (width > height) {
    glScaled(((double) height / width), 1.0, 1.0);
    move_t->offset_x = (uint) (width - height) * 0.5;
    move_t->offset_y = 0;
  }
  else {
    glScaled(1.0, ((double) width / height), 1.0);
    move_t->offset_x = 0;
    move_t->offset_y = (uint) (height - width) * 0.5;
  }

  /* Factors for finding choosen field */
  move_t->pixel_x = width;
  move_t->pixel_y = height;
}

  /* Reshape of a triangulation */
static void
reshape_mesh(int width, int height)
{

  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  /* Scaling */
  if (width > height) {
    glScaled(((double) height / width), 1.0, 1.0);
  }
  else {
    glScaled(1.0, ((double) width / height), 1.0);
  }

  gluPerspective(35.0, 1.0, 1.0, 30.0);
  glTranslated(0.0, 0.0, -4.0);

  glEnable(GL_DEPTH_TEST);
}

  /* Reshape for legends */
static void
reshape_legend(int width, int height)
{

  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  /* Scaling */
  if (width > height) {
    glScaled(((double) height / width), 1.0, 1.0);
  }
  else {
    glScaled(1.0, ((double) width / height), 1.0);
  }
}

 /*
  *     GLUT mouse functions    *
  */


  /* Mouse movement on a directional cluster */
static void
mouse_dcluster(int button, int state, int position_y, int position_x)
{

  uint      w;
  pactivemovement move;
  w = glutGetWindow();
  if (w == clust1) {
    move = move_1;
  }
  else {
    move = move_2;
  }
  switch (button) {
    /* Set values for roation of the given object */
  case GLUT_LEFT_BUTTON:
    if (state == GLUT_DOWN) {

      move->old_position_x = position_x;
      move->old_position_y = position_y;

      move->old_angle_x = move->angle_x;
      move->old_angle_y = move->angle_y;

    }
    break;
    /* Switch to mode for drawing the father */
  case GLUT_RIGHT_BUTTON:
    if (state == GLUT_DOWN) {
      if (move->levelup == true) {
	if (move->dblock !=NULL) {
	  up = 2;
	  glutPostWindowRedisplay(w);
	}
      }

    }
    break;

  default:
    break;
  }
}

  /* Mouse movement on a cluster */
static void
mouse_cluster(int button, int state, int position_y, int position_x)
{

  uint      w;
  pactivemovement move;
  w = glutGetWindow();
  if (w == clust1) {
    move = move_1;
  }
  else {
    move = move_2;
  }
  switch (button) {
    /* set values for roation of the given object */
  case GLUT_LEFT_BUTTON:
    if (state == GLUT_DOWN) {

      move->old_position_x = position_x;
      move->old_position_y = position_y;

      move->old_angle_x = move->angle_x;
      move->old_angle_y = move->angle_y;

    }
    break;
    /* Switch to mode for drawing the father */
  case GLUT_RIGHT_BUTTON:
    if (state == GLUT_DOWN) {
      if (move->levelup == true) {
	if (move->block !=NULL) {
	  up = 2;
	  glutPostWindowRedisplay(w);
	}
      }

    }
    break;

  default:
    break;
  }
}

  /* Mouse movement on window with level conditions   -directional version */
static void
mouse_dblock(int button, int state, int position_x, int position_y)
{

  uint      i;
  uint      step;
  if (version == 1) {
    step = getdepth_dblock(bl);
  }
  else {
    step = getdepth_dcluster(t);
  }

  switch (button) {
    /* Evaluate choosen level and change values level and level2 */
  case GLUT_LEFT_BUTTON:
    if (state == GLUT_DOWN) {

      move_o->target_x = position_x * move_o->orthoscaling;
      move_o->target_y = position_y * move_o->orthoscaling;

      i = 0;
      while (i <= step) {
	if (move_o->target_y > 2.0 * i / (step + 1)) {
	  level2 = i;
	}
	else {
	  break;
	}
	i++;
      }

      if (move_o->target_x < 1.0) {
	level = level2;
      }
      else {
	level = 0;
      }

      visu = 1;
      if (clust2 != 0) {
	glutPostWindowRedisplay(clust2);
      }
      glutPostWindowRedisplay(clust1);

    }
    break;

  default:
    break;
  }
}

  /* Mouse movement on window with level conditions */
static void
mouse_block(int button, int state, int position_x, int position_y)
{

  uint      i;
  uint      step;
  if (version == 1) {
    step = getdepth_block(bbl);
  }
  else {
    step = getdepth_cluster(ct);
  }

  switch (button) {
    /* Evaluate choosen level and change values level and level2 */
  case GLUT_LEFT_BUTTON:
    if (state == GLUT_DOWN) {

      move_o->target_x = position_x * move_o->orthoscaling;
      move_o->target_y = position_y * move_o->orthoscaling;

      i = 0;
      while (i <= step) {
	if (move_o->target_y > 2.0 * i / (step + 1)) {
	  level2 = i;
	}
	else {
	  break;
	}
	i++;
      }

      if (move_o->target_x < 1.0) {
	level = level2;
      }
      else {
	level = 0;
      }

      visu = 1;
      if (clust2 != 0) {
	glutPostWindowRedisplay(clust2);
      }
      glutPostWindowRedisplay(clust1);

    }
    break;

  default:
    break;
  }
}

  /* Mouse movement on an orthogonal projection */
static void
mouse_ortho(int button, int state, int position_x, int position_y)
{

  switch (button) {
    /* Evaluate choosen directional block */
  case GLUT_LEFT_BUTTON:
    if (state == GLUT_DOWN) {
      move_o->target_x =
	(position_x -
	 move_o->offset_x) * move_o->orthoscaling / move_o->clustscaling_x;
      move_o->target_x -= move_o->xoff;
      move_o->target_y =
	(position_y -
	 move_o->offset_y) * move_o->orthoscaling / move_o->clustscaling_y;
      move_o->target_y -= move_o->yoff;

      visu = 1;
      up = 0;
      glutPostWindowRedisplay(clust1);
      glutPostWindowRedisplay(ortho);
    }
    break;

  default:
    break;
  }
}

  /* Mouse movement on a triangulation */
static void
mouse_mesh(int button, int state, int position_y, int position_x)
{

  switch (button) {

  case GLUT_LEFT_BUTTON:
    if (state == GLUT_DOWN) {

      move_1->old_position_x = position_x;
      move_1->old_position_y = position_y;

      move_1->old_angle_x = move_1->angle_x;
      move_1->old_angle_y = move_1->angle_y;
    }
    break;

  default:
    break;
  }
}

  /* Mouse movement on a two dimensional triangulation table */
static void
mouse_table_tri(int button, int state, int position_x, int position_y)
{

  uint      h;
  uint      p;
  uint      i = 1;
  uint      last_target = move_t->target;

  switch (button) {
  case GLUT_LEFT_BUTTON:
    if (state == GLUT_DOWN) {
      if (position_x < move_t->pixel_x / 2) {
	if (position_x < move_t->offset_x) {
	  printf("Please choose a current triangle or change page.\n");
	  return;
	}
	else {
	  if ((position_y > move_t->offset_y)
	      && (position_y < move_t->pixel_y - move_t->offset_y)) {
	    h = (move_t->pixel_y - 2 * move_t->offset_y) / 6;
	    move_t->target = (move_t->page) * 10;
	    p = position_y - move_t->offset_y;
	    while ((p > h * i) && (i < 6)) {
	      move_t->target += 1;
	      i++;
	    }
	    if (i == 6) {
	      move_t->target = last_target;
	      printf("Please choose a current triangle or change page.\n");
	      return;
	    }
	  }
	  else {
	    printf("Please choose a current triangle or change page.\n");
	    return;
	  }
	}
      }
      else {
	if (position_x > (move_t->pixel_x - move_t->offset_x)) {
	  printf("Please choose a current triangle.\n");
	  return;
	}
	else {
	  if ((position_y > move_t->offset_y)
	      && (position_y < move_t->pixel_y - move_t->offset_y)) {
	    h = (move_t->pixel_y - 2 * move_t->offset_y) / 6;
	    move_t->target = (move_t->page) * 10 + 5;
	    p = position_y - move_t->offset_y;
	    while ((p > h * i) && (i < 6)) {
	      move_t->target += 1;
	      i++;
	    }
	    if (i == 6) {
	      move_t->target = last_target;
	      if (position_x >
		  (move_t->pixel_x / 2 +
		   (move_t->pixel_x - 2 * move_t->offset_x) / 4)) {
		if (move_t->max_page != 0) {
		  move_t->page =
		    ((move_t->page + 1) <=
		     (move_t->max_page) ? (move_t->page + 1) : 0);
		}
	      }
	      else {
		if (move_t->max_page != 0) {
		  move_t->page =
		    ((move_t->page) !=
		     0 ? (move_t->page - 1) : (move_t->max_page));
		}
	      }
	    }
	    glutPostWindowRedisplay(table);
	  }
	  else {
	    printf("Please choose a current triangle or change page.\n");
	    return;
	  }
	}
      }
      if (i == 6) {
	mode = 0;
      }
      else {
	if (move_t->target < tr->triangles) {
	  mode = 1;
	  glutPostWindowRedisplay(one);
	  glutPostWindowRedisplay(table);
	}
	else {
	  mode = 0;
	  printf("Please choose a current triangle or change page.\n");
	  glutPostWindowRedisplay(one);
	  glutPostWindowRedisplay(table);
	}
      }
    }
    break;

  default:
    break;
  }
}

  /* GLUT callback function for mouse movement on a three dimensional triangulation table */
static void
mouse_table_tet(int button, int state, int position_x, int position_y)
{

  uint      h;
  uint      p;
  uint      i = 1;
  uint      last_target = move_t->target;

  switch (button) {
  case GLUT_LEFT_BUTTON:
    if (state == GLUT_DOWN) {
      if (position_x < move_t->pixel_x / 2) {
	if (position_x < move_t->offset_x) {
	  printf("Please choose a current tetrahedra or change page.\n");
	  return;
	}
	else {
	  if ((position_y > move_t->offset_y)
	      && (position_y < move_t->pixel_y - move_t->offset_y)) {
	    h = (move_t->pixel_y - 2 * move_t->offset_y) / 5;
	    move_t->target = (move_t->page) * 8;
	    p = position_y - move_t->offset_y;
	    while ((p > h * i) && (i < 5)) {
	      move_t->target += 1;
	      i++;
	    }
	    if (i == 5) {
	      move_t->target = last_target;
	      printf("Please choose a current tetrahedra or change page.\n");
	      return;
	    }
	  }
	  else {
	    printf("Please choose a current tetrahedra or change page.\n");
	    return;
	  }
	}
      }
      else {
	if (position_x > (move_t->pixel_x - move_t->offset_x)) {
	  printf("Please choose a current tetrahedra.\n");
	  return;
	}
	else {
	  if ((position_y > move_t->offset_y)
	      && (position_y < move_t->pixel_y - move_t->offset_y)) {
	    h = (move_t->pixel_y - 2 * move_t->offset_y) / 5;
	    move_t->target = (move_t->page) * 8 + 4;
	    p = position_y - move_t->offset_y;
	    while ((p > h * i) && (i < 5)) {
	      move_t->target += 1;
	      i++;
	    }
	    if (i == 5) {
	      move_t->target = last_target;
	      if (position_x >
		  (move_t->pixel_x / 2 +
		   (move_t->pixel_x - 2 * move_t->offset_x) / 4)) {
		if (move_t->max_page != 0) {
		  move_t->page =
		    ((move_t->page + 1) <=
		     (move_t->max_page) ? (move_t->page + 1) : 0);
		}
	      }
	      else {
		if (move_t->max_page != 0) {
		  move_t->page =
		    ((move_t->page) !=
		     0 ? (move_t->page - 1) : (move_t->max_page));
		}
	      }
	    }
	    glutPostWindowRedisplay(table);
	  }
	  else {
	    printf("Please choose a current tetrahedra or change page.\n");
	    return;
	  }
	}
      }
      if (i == 5) {
	mode = 0;
      }
      else {
	if (move_t->target < th->tetrahedra) {
	  mode = 1;
	  glutPostWindowRedisplay(one);
	  glutPostWindowRedisplay(table);
	}
	else {
	  mode = 0;
	  printf("Please choose a current tetrahedra or change page.\n");
	  glutPostWindowRedisplay(one);
	  glutPostWindowRedisplay(table);
	}
      }
    }
    break;

  default:
    break;
  }
}

  /* Mouse movement on a surface table */
static void
mouse_table_sur(int button, int state, int position_x, int position_y)
{

  uint      h;
  uint      p;
  uint      i = 1;
  uint      last_target = move_t->target;

  switch (button) {
  case GLUT_LEFT_BUTTON:
    if (state == GLUT_DOWN) {
      if (position_x < move_t->pixel_x / 2) {
	if (position_x < move_t->offset_x) {
	  printf("Please choose a current triangle or change page.\n");
	  return;
	}
	else {
	  if ((position_y > move_t->offset_y)
	      && (position_y < move_t->pixel_y - move_t->offset_y)) {
	    h = (move_t->pixel_y - 2 * move_t->offset_y) / 5;
	    move_t->target = (move_t->page) * 8;
	    p = position_y - move_t->offset_y;
	    while ((p > h * i) && (i < 5)) {
	      move_t->target += 1;
	      i++;
	    }
	    if (i == 5) {
	      move_t->target = last_target;
	      printf("Please choose a current triangle or change page.\n");
	      return;
	    }
	  }
	  else {
	    printf("Please choose a current triangle or change page.\n");
	    return;
	  }
	}
      }
      else {
	if (position_x > (move_t->pixel_x - move_t->offset_x)) {
	  printf("Please choose a current triangle.\n");
	  return;
	}
	else {
	  if ((position_y > move_t->offset_y)
	      && (position_y < move_t->pixel_y - move_t->offset_y)) {
	    h = (move_t->pixel_y - 2 * move_t->offset_y) / 5;
	    move_t->target = (move_t->page) * 8 + 4;
	    p = position_y - move_t->offset_y;
	    while ((p > h * i) && (i < 5)) {
	      move_t->target += 1;
	      i++;
	    }
	    if (i == 5) {
	      move_t->target = last_target;
	      if (position_x >
		  (move_t->pixel_x / 2 +
		   (move_t->pixel_x - 2 * move_t->offset_x) / 4)) {
		if (move_t->max_page != 0) {
		  move_t->page =
		    ((move_t->page + 1) <=
		     (move_t->max_page) ? (move_t->page + 1) : 0);
		}
	      }
	      else {
		if (move_t->max_page != 0) {
		  move_t->page =
		    ((move_t->page) !=
		     0 ? (move_t->page - 1) : (move_t->max_page));
		}
	      }
	    }
	    glutPostWindowRedisplay(table);
	  }
	  else {
	    printf("Please choose a current triangle or change page.\n");
	    return;
	  }
	}
      }
      if (i == 5) {
	mode = 0;
      }
      else {
	if (move_t->target < gr->triangles) {
	  mode = 1;
	  glutPostWindowRedisplay(one);
	  glutPostWindowRedisplay(table);
	}
	else {
	  mode = 0;
	  printf("Please choose a current triangle or change page.\n");
	  glutPostWindowRedisplay(one);
	  glutPostWindowRedisplay(table);
	}
      }
    }
    break;

  default:
    break;
  }
}


  /*
   *    GLUT motion functions   *
   */

  /* Motion on a directional cluster object */
static void
motion_dcluster(int position_y, int position_x)
{

  uint      w;
  pactivemovement move;
  w = glutGetWindow();

  if (w == clust1) {
    move = move_1;
  }
  else {
    move = move_2;
  }

  move->angle_x = move->old_angle_x + (position_x - move->old_position_x);
  move->angle_y = move->old_angle_y + (position_y - move->old_position_y);

  glutPostRedisplay();
}

  /* Motion a cluster object */
static void
motion_cluster(int position_y, int position_x)
{

  uint      w;
  pactivemovement move;
  w = glutGetWindow();

  if (w == clust1) {
    move = move_1;
  }
  else {
    move = move_2;
  }

  move->angle_x = move->old_angle_x + (position_x - move->old_position_x);
  move->angle_y = move->old_angle_y + (position_y - move->old_position_y);

  glutPostRedisplay();
}

  /* Motion on a triangulation */
static void
motion_mesh(int position_y, int position_x)
{

  move_1->angle_x =
    move_1->old_angle_x + (position_x - move_1->old_position_x);
  move_1->angle_y =
    move_1->old_angle_y + (position_y - move_1->old_position_y);

  glutPostRedisplay();

}

  /*
   *    GLUT key functions              *
   */

  /* Key input for an orthogonal projection */
static void
key_ortho(unsigned char key, int x, int y)
{

  switch (key) {
    /*press 'esc' for close window */
  case 27:
    glutLeaveMainLoop();
    break;

    /* press 'space' for a snapshot */
  case 32:
    printf("Take a snapshot\n");
    take_a_shot[0] = true;
    glutPostRedisplay();
    break;

  default:
    break;
  }
}

  /* Key input for a directional cluster */
static void
key_dcluster(unsigned char key, int x, int y)
{

  uint      w;
  pactivemovement move;
  w = glutGetWindow();

  if (w == clust1) {
    move = move_1;
  }
  else {
    move = move_2;
  }

  switch (key) {
    /* Press 'esc' for close window */
  case 27:
    glutLeaveMainLoop();
    break;

    /* Press '+' for zooming in */
  case 43:
    move->zoom = move->zoom + 0.1;
    glutPostRedisplay();
    break;

    /* Press '-' for zooming out */
  case 45:
    move->zoom = move->zoom - 0.1;
    glutPostRedisplay();
    break;

    /* Press 'n' for remove zooming and rotation */
  case 110:
    move->zoom = 1.0;
    move->angle_x = 0.0;
    move->angle_y = 0.0;
    move->rightleft = 0.0;
    move->updown = 0.0;
    glutPostRedisplay();
    break;
    /* Press 'a' translate left */
  case 97:
    move_1->rightleft = move_1->rightleft - 0.1;
    glutPostRedisplay();
    break;
    /* Press 'd' translate right */
  case 100:
    move_1->rightleft = move_1->rightleft + 0.1;
    glutPostRedisplay();
    break;
    /* Press 'w' translate up */
  case 119:
    move_1->updown = move_1->updown + 0.1;
    glutPostRedisplay();
    break;
    /* Press 's' translate down */
  case 115:
    move_1->updown = move_1->updown - 0.1;
    glutPostRedisplay();
    break;
    /* Press 'p' to switch between visualisation with and without points */
  case 112:
    if (gt == NULL) {
      return;
    }
    points = 1 - points;
    if (clust2 != 0) {
      glutPostWindowRedisplay(clust1);
      glutPostWindowRedisplay(clust2);
    }
    else {
      glutPostRedisplay();
    }

    break;
    /* Press '>' for directions */
  case 62:
    direction = 1 - direction;
    if (clust2 != 0) {
      glutPostWindowRedisplay(clust1);
      glutPostWindowRedisplay(clust2);
    }
    else {
      glutPostRedisplay();
    }
    break;
    /* press 'space' for a snapshot */
  case 32:
    printf("Take a snapshot\n");
    take_a_shot[1] = true;
    glutPostRedisplay();
    break;
    /* press 'b' for bounding box coordinates */
  case 98:
    bbox = true;
    glutPostRedisplay();
    break;
  default:
    break;
  }
}

  /* Key input for a cluster */
static void
key_cluster(unsigned char key, int x, int y)
{

  uint      w;
  pactivemovement move;
  w = glutGetWindow();

  if (w == clust1) {
    move = move_1;
  }
  else {
    move = move_2;
  }

  switch (key) {
    /* Press 'esc' for close window */
  case 27:
    glutLeaveMainLoop();
    break;

    /* Press '+' for zooming in */
  case 43:
    move->zoom = move->zoom + 0.1;
    glutPostRedisplay();
    break;

    /* Press '-' for zooming out */
  case 45:
    move->zoom = move->zoom - 0.1;
    glutPostRedisplay();
    break;

    /* Press 'n' for remove zooming and rotation */
  case 110:
    move->zoom = 1.0;
    move->angle_x = 0.0;
    move->angle_y = 0.0;
    move->rightleft = 0.0;
    move->updown = 0.0;
    glutPostRedisplay();
    break;
    /* Press 'a' translate left */
  case 97:
    move_1->rightleft = move_1->rightleft - 0.1;
    glutPostRedisplay();
    break;
    /* Press 'd' translate right */
  case 100:
    move_1->rightleft = move_1->rightleft + 0.1;
    glutPostRedisplay();
    break;
    /* Press 'w' translate up */
  case 119:
    move_1->updown = move_1->updown + 0.1;
    glutPostRedisplay();
    break;
    /* Press 's' translate down */
  case 115:
    move_1->updown = move_1->updown - 0.1;
    glutPostRedisplay();
    break;
    /* Press 'p' to switch between visualisation with and without points */
  case 112:
    if (gt == NULL) {
      return;
    }
    points = 1 - points;
    if (clust2 != 0) {
      glutPostWindowRedisplay(clust1);
      glutPostWindowRedisplay(clust2);
    }
    else {
      glutPostRedisplay();
    }

    break;
    /* press 'space' for a snapshot */
  case 32:
    printf("Take a snapshot\n");
    take_a_shot[1] = true;
    glutPostRedisplay();
    break;
    /* press 'b' for bounding box coordinates */
  case 98:
    bbox = true;
    glutPostRedisplay();
    break;

  default:
    break;
  }
}

  /* Key input for triangulation */
static void
key_mesh(unsigned char key, int x, int y)
{

  (void) x;
  (void) y;

  switch (key) {
    /* Press 'esc' for close window */
  case 27:
    glutLeaveMainLoop();
    break;

    /* Press '+' for zooming in */
  case 43:
    move_1->zoom = move_1->zoom + 0.1;
    glutPostRedisplay();
    break;

    /* Press '-' for zooming out */
  case 45:
    move_1->zoom = move_1->zoom - 0.1;
    glutPostRedisplay();
    break;

    /* Press 'n' for remove zooming and rotation */
  case 110:
    move_1->zoom = 1.0;
    move_1->angle_x = 0.0;
    move_1->angle_y = 0.0;
    move_1->rightleft = 0.0;
    move_1->updown = 0.0;
    glutPostRedisplay();
    break;
    /* Press 'a' translate left */
  case 97:
    move_1->rightleft = move_1->rightleft - 0.1;
    glutPostRedisplay();
    break;
    /* Press 'd' translate right */
  case 100:
    move_1->rightleft = move_1->rightleft + 0.1;
    glutPostRedisplay();
    break;
    /* Press 'w' translate up */
  case 119:
    move_1->updown = move_1->updown + 0.1;
    glutPostRedisplay();
    break;
    /* Press 's' translate down */
  case 115:
    move_1->updown = move_1->updown - 0.1;
    glutPostRedisplay();
    break;
    /* Press 'p' to switch between visualisation with and without points */
  case 112:
    points = 1 - points;
    glutPostRedisplay();
    break;
    /* Press 'm' to switch between drawing triangles with edges and triangles */
  case 109:
    mode = 1 - mode;
    glutPostRedisplay();
    break;
    /* Press 'c' to switch between drawing coordinate system or not */
  case 99:
    coord = 1 - coord;
    glutPostRedisplay();
    break;
    /* Press 'x' to light x coordinate */
  case 120:
    light_coord[0] = 1;
    light_coord[1] = 0;
    light_coord[2] = 0;
    glutPostRedisplay();
    break;
    /* Press 'y' to light x coordinate */
  case 121:
    light_coord[0] = 0;
    light_coord[1] = 1;
    light_coord[2] = 0;
    glutPostRedisplay();
    break;
    /* Press 'z' to light x coordinate */
  case 122:
    light_coord[0] = 0;
    light_coord[1] = 0;
    light_coord[2] = 1;
    glutPostRedisplay();
    break;
    /* Press '0' to switch between drawing origin or not */
  case 48:
    origin = 1 - origin;
    glutPostRedisplay();
    break;
    /* press 'space' for a snapshot */
  case 32:
    printf("Take a snapshot\n");
    take_a_shot[1] = true;
    glutPostRedisplay();
    break;

  default:
    break;
  }
}

  /* Key input for solution */
static void
key_solution(unsigned char key, int x, int y)
{

  (void) x;
  (void) y;
  uint      o;

  switch (key) {
    /* Press 'esc' for close window */
  case 27:
    glutLeaveMainLoop();
    break;

    /* Press '+' for zooming in */
  case 43:
    move_1->zoom = move_1->zoom + 0.1;
    glutPostRedisplay();
    break;

    /* Press '-' for zooming out */
  case 45:
    move_1->zoom = move_1->zoom - 0.1;
    glutPostRedisplay();
    break;

    /* Press 'n' for remove zooming and rotation */
  case 110:
    move_1->zoom = 1.0;
    move_1->angle_x = 0.0;
    move_1->angle_y = 0.0;
    move_1->rightleft = 0.0;
    move_1->updown = 0.0;
    sol->time = 0.0;
    sol->time_step = 1.0;
    glutPostRedisplay();
    break;
    /* Press 'a' translate left */
  case 97:
    move_1->rightleft = move_1->rightleft - 0.1;
    glutPostRedisplay();
    break;
    /* Press 'd' translate right */
  case 100:
    move_1->rightleft = move_1->rightleft + 0.1;
    glutPostRedisplay();
    break;
    /* Press 'w' translate up */
  case 119:
    move_1->updown = move_1->updown + 0.1;
    glutPostRedisplay();
    break;
    /* Press 's' translate down */
  case 115:
    move_1->updown = move_1->updown - 0.1;
    glutPostRedisplay();
    break;
    /* Press 'p' to switch between visualisation with and without plane with solution */
  case 112:
    sol->plane = 1 - sol->plane;
    glutPostRedisplay();
    break;
    /* Press 'm' to switch between visible model or not */
  case 109:
    sol->model = 1 - sol->model;
    glutPostRedisplay();
    break;
    /* Press 'c' to switch between drawing coordinate system or not */
  case 99:
    coord = 1 - coord;
    glutPostRedisplay();
    break;
    /* Press 'x' to light x coordinate */
  case 120:
    light_coord[0] = 1;
    light_coord[1] = 0;
    light_coord[2] = 0;
    glutPostRedisplay();
    break;
    /* Press 'y' to light x coordinate */
  case 121:
    light_coord[0] = 0;
    light_coord[1] = 1;
    light_coord[2] = 0;
    glutPostRedisplay();
    break;
    /* Press 'z' to light x coordinate */
  case 122:
    light_coord[0] = 0;
    light_coord[1] = 0;
    light_coord[2] = 1;
    glutPostRedisplay();
    break;
    /* Press '0' to switch between drawing origin or not */
  case 48:
    origin = 1 - origin;
    glutPostRedisplay();
    break;
    /* press 'space' for a snapshot */
  case 32:
    printf("Take a snapshot\n");
    take_a_shot[1] = true;
    glutPostRedisplay();
    break;

    /* press '.' for next solution */
  case 46:
    if (sol->grids > 0) {
      o = sol->current_grid_number + 1;
      sol->current_grid_number = (o < sol->grids ? o : 0);
      sol->type = 0;
      sol->time = 0;
      animation = 0;
      glutPostWindowRedisplay(one);
      glutPostWindowRedisplay(table);
    }
    break;

    /* press 't' for switch types */
  case 116:
    sol->type = (sol->type + 1) % 3;
    sol->time = 0;
    animation = 0;
    glutPostWindowRedisplay(one);
    glutPostWindowRedisplay(table);
    break;

    /* press 'enter' to start animation */
  case 13:
    if (sol->animate == true) {
      animation = 1 - animation;
      glutPostWindowRedisplay(one);
    }
    break;

    /* press ' < ' for slower animation */
  case 60:
    sol->time_step -= 0.1 * sol->time_step;
    break;
    /* press ' > ' for faster animation */
  case 62:
    sol->time_step += 0.1 * sol->time_step;
    break;

  default:
    break;
  }
}

  /* Timer for evaluate time dependent solution out of Helmholtz data */
static void
timer_helmholtz(int value)
{

  longindex i, l;
  uint      current = sol->current_grid_number;
  real      scale = sol->scale;
  real      t = sol->time;
  field     e, tmp;
  real      w;
  uint      j;
  real      lenght[3], maxpoint[3];

  if (animation == true) {
    l = (sol->gr[current]->x) * (sol->gr[current]->y);
    w = 2.0 * M_PI / 25.0;	/* Angular frequence divieded by frame rate */
    // w = kappa*342.0; 
    /* Compute new solution with seperation ansatz and the solution exp(-iwt +phi) where w is the angular frequence
       and we use the simple case phi = 0 */
    e = cexp((-1.0) * I * t * w);

    /* Evaluate color values ..maybe some points aren't inside anymore */
    for (j = 0; j < 3; j++) {
      lenght[j] = scale / (sol->max[current][j] - sol->min[current][j]);
      maxpoint[j] = (sol->max[current][j] + sol->gr[current]->z) / scale;
    }


    for (i = 0; i < l; i++) {
      /* value */
      tmp =
	(sol->timegr[current]->point[i][0] +
	 I * sol->timegr[current]->point[i][1]) * e;

      sol->timegr[current]->point[i][2] = REAL(tmp);
      sol->timegr[current]->point[i][4] = IMAG(tmp);
      sol->timegr[current]->point[i][6] = ABS(tmp);
      /* color */
      sol->timegr[current]->point[i][3] =
	1.0 +
	(lenght[0] * 2.0 * (sol->timegr[current]->point[i][2] - maxpoint[0]));
      sol->timegr[current]->point[i][5] =
	1.0 +
	(lenght[1] * 2.0 * (sol->timegr[current]->point[i][4] - maxpoint[1]));
      sol->timegr[current]->point[i][7] =
	1.0 +
	(lenght[2] * 2.0 * (sol->timegr[current]->point[i][6] - maxpoint[2]));
    }
    sol->time += sol->time_step;
    glutPostWindowRedisplay(one);
  }

  glutTimerFunc(25, timer_helmholtz, 0);
}

/**********************************************************
*		        Visualize bounding boxes		                  *
***********************************************************/

void
visualize_dblock_bbox(pcdblock b, pclustergeometry grc, pclustergeometry gcc,
		      bool c, int argc, char **argv)
{

  t = b->rc;
  s = b->cc;
  gt = grc;
  gs = gcc;
  bl = b;
  color = c;
  version = 1;
  visu = 0;

  /* Checking if bounding boxes are given */
  if ((getdiam_2_dcluster(t) + getdiam_2_dcluster(s)) == 0.0) {
    printf("Drawing isn't possible, because bounding boxes are empty.");
    return;
  }

  if (gt == NULL || gs == NULL) {
    if (gt == NULL) {
      gs = gt;
    }
    if (gs == NULL) {
      gt = gs;
    }
    points = 0;
    printf("- - - - - - Caution! - - - - - -\n");
    printf("Missing clustergeometry for drawing points!\n");
  }

  if (t->dim > 3) {
    printf("Row cluster dimension is with ");
    printf("%u", t->dim);
    printf(" to hight for drawing!");

    return;
  }

  if (s->dim > 3) {
    printf("Column cluster dimension is with ");
    printf("%u", s->dim);
    printf(" to hight for drawing!");

    return;
  }

//    glutInit(&argc, argv);
  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
		GLUT_ACTION_GLUTMAINLOOP_RETURNS);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);

  glReadBuffer(GL_BACK);	/* Allowed us to read pixel out of the buffer */
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);	/* How to store the pixel */

  move_o = new_orthomovement();
  ortho = glutCreateWindow("Select level");
  glutPositionWindow(50, 300);
  glutReshapeWindow(400, 400);
  glutDisplayFunc(display_dblock);
  glutReshapeFunc(reshape_ortho);
  glutMouseFunc(mouse_dblock);
  glutKeyboardFunc(key_ortho);

  move_1 = new_activemovement(false);
  clust1 = glutCreateWindow("Visualize row cluster for certain block level");
  glutPositionWindow(550, 50);
  glutReshapeWindow(400, 400);
  glutDisplayFunc(display_two_dcluster);
  glutReshapeFunc(reshape_dcluster);
  glutMouseFunc(mouse_dcluster);
  glutMotionFunc(motion_dcluster);
  glutKeyboardFunc(key_dcluster);

  move_2 = new_activemovement(false);
  clust2 =
    glutCreateWindow("Visualize column cluster for certain block level");
  glutPositionWindow(950, 50);
  glutReshapeWindow(400, 400);
  glutDisplayFunc(display_two_dcluster);
  glutReshapeFunc(reshape_dcluster);
  glutMouseFunc(mouse_dcluster);
  glutMotionFunc(motion_dcluster);
  glutKeyboardFunc(key_dcluster);

  printf("------------------------------------\n");
  printf("Visualize Block bounding boxes!\n");
  printf("------------------------------------\n");
  printf("Use left window to select level to\n");
  printf("be drawn.\n");
  printf("Right window will show the chosen\n");
  printf("bounding boxes.\n");
  printf("If bounding boxes are 3-dimensional\n");
  printf("use mouse to chance perspective.\n");
  printf("To zoom in use '+' and '-' to zoom\n");
  printf("out.\n");
  printf("Keys 'a', 'd', 'w' and 's' translate\n");
  printf("perspective to the left, right, up\n");
  printf("and down.\n");
  printf("With 'n' perspective and zoom will\n");
  printf("be set to start values.\n");
  printf("Switch with 'p' between showing and\n");
  printf("hiding points.\n");
  printf("Fade in directions with '>'.\n");
  printf("Make snapshots of the content from one\n");
  printf("selcted window with 'space'. \n");
  printf("Leave visualization with 'esc'.\n");
  printf("------------------------------------\n");

  glutMainLoop();

  printf("Visualization has been closed.\n");
  del_activemovement(move_1);
  del_activemovement(move_2);
  del_orthomovement(move_o);
}


void
visualize_dblock_certain_bbox(pcdblock b, pclustergeometry grc,
			      pclustergeometry gcc, int argc, char **argv)
{

  visualize_dblock_level_certain_bbox(b, grc, gcc, 0, getdepth_dblock(b),
				      argc, argv);

}


void
visualize_dblock_level_certain_bbox(pcdblock b, pclustergeometry grc,
				    pclustergeometry gcc, uint l1, uint l2,
				    int argc, char **argv)
{

  t = b->rc;
  gt = grc;
  s = b->cc;
  gs = gcc;
  bl = b;
  version = 1;
  level = l1;
  level2 = l2;
  uint      depth;
  up = 0;

  assert((int) level >= 0);
  assert(level <= level2);

  /* Checking if bounding boxes are given */
  if ((getdiam_2_dcluster(t) + getdiam_2_dcluster(s)) == 0.0) {
    printf("Drawing isn't possible, because bounding boxes are empty.");
    return;
  }

  if (gt == NULL || gs == NULL) {
    if (gt == NULL) {
      gs = gt;
    }
    if (gs == NULL) {
      gt = gs;
    }
    points = 0;
    printf("- - - - - - Caution! - - - - - -\n");
    printf("Missing clustergeometry for drawing points!\n");
  }

  if (t->dim > 3) {
    printf("Row cluster dimension is with ");
    printf("%u", t->dim);
    printf(" to hight for drawing!");

    return;
  }

  if (s->dim > 3) {
    printf("Column cluster dimension is with ");
    printf("%u", s->dim);
    printf(" to hight for drawing!");

    return;
  }

  depth = getdepth_dblock(b);
  if (depth < level) {
    printf("Can't draw level ");
    printf("%u", level);
    printf(", because the block tree has only depth ");
    printf("%u.", depth);
    return;
  }

  if (depth < level2) {
    printf("Can't draw level ");
    printf("%u", level2);
    printf(", because the block tree has only depth ");
    printf("%u.", depth);
    return;
  }

//      glutInit(&argc, argv);
  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
		GLUT_ACTION_GLUTMAINLOOP_RETURNS);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);

  glReadBuffer(GL_BACK);	/* Allowed us to read pixel out of the buffer */
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);	/* How to store the pixel */

  move_o = new_orthomovement();
  ortho = glutCreateWindow("Projection blocks");
  glutPositionWindow(150, 100);
  glutReshapeWindow(600, 600);
  glutDisplayFunc(display_ortho_d);
  glutReshapeFunc(reshape_ortho);
  glutKeyboardFunc(key_ortho);
  glutMouseFunc(mouse_ortho);

  move_1 = new_activemovement(true);
  clust1 = glutCreateWindow("Visualize certain cluster");
  glutPositionWindow(800, 100);
  glutReshapeWindow(400, 400);
  glutDisplayFunc(display_certain_dblock);
  glutReshapeFunc(reshape_dcluster);
  glutMouseFunc(mouse_dcluster);
  glutMotionFunc(motion_dcluster);
  glutKeyboardFunc(key_dcluster);

  printf("------------------------------------\n");
  printf("Visualize Block bounding boxes!\n");
  printf("------------------------------------\n");
  printf("Use left window to select block to\n");
  printf("be drawn.\n");
  printf("Right window will show the chosen\n");
  printf("clusters with bounding boxes.\n");
  printf("Use rigth mouse button to switch to\n");
  printf("the father block.\n");
  printf("If bounding boxes are 3-dimensional\n");
  printf("use left mouse button to chance\n");
  printf("perspective.\n");
  printf("To zoom in use '+' and '-' to zoom\n");
  printf("out.\n");
  printf("Keys 'a', 'd', 'w' and 's' translate\n");
  printf("perspective to the left, right, up\n");
  printf("and down.\n");
  printf("With 'n' perspective and zoom will\n");
  printf("be set to start values.\n");
  printf("Switch with 'p' between showing and\n");
  printf("hiding points.\n");
  printf("Print bounding box values with 'b'.\n");
  printf("Fade in directions with '>'.\n");
  printf("Make snapshots of the content from one\n");
  printf("selcted window with 'space'. \n");
  printf("Leave visualization with 'esc'.\n");
  printf("------------------------------------\n");


  glutMainLoop();

  printf("Visualization has been closed.\n");
  del_activemovement(move_1);
  del_orthomovement(move_o);
}

void
visualize_dcluster_bbox(pcdcluster ti, pclustergeometry gti, bool c, int argc,
			char **argv)
{

  t = ti;
  s = t;
  gt = gti;
  color = c;
  version = 0;
  visu = 0;

  /* Checking if bounding boxes are given */
  if (getdiam_2_dcluster(t) == 0.0) {
    printf("Drawing isn't possible, because bounding boxes are empty.");
    return;
  }

  if (gt == NULL) {
    points = 0;
    printf("- - - - - - Caution! - - - - - -\n");
    printf("Missing clustergeometry for drawing points!\n");
  }

  if (t->dim > 3) {
    printf("Dimension is with ");
    printf("%u", t->dim);
    printf(" to hight for drawing!");

    return;
  }

//    glutInit(&argc, argv);
  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
		GLUT_ACTION_GLUTMAINLOOP_RETURNS);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);

  glReadBuffer(GL_BACK);	/* Allowed us to read pixel out of the buffer */
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);	/* How to store the pixel */

  move_o = new_orthomovement();
  ortho = glutCreateWindow("Select level");
  glutPositionWindow(50, 300);
  glutReshapeWindow(400, 400);
  glutDisplayFunc(display_dblock);
  glutReshapeFunc(reshape_ortho);
  glutMouseFunc(mouse_dblock);
  glutKeyboardFunc(key_ortho);

  move_1 = new_activemovement(false);
  clust1 = glutCreateWindow("Visualize cluster bounding boxes");
  glutPositionWindow(550, 100);
  glutReshapeWindow(600, 600);
  glutDisplayFunc(display_dcluster);
  glutReshapeFunc(reshape_dcluster);
  glutMouseFunc(mouse_dcluster);
  glutMotionFunc(motion_dcluster);
  glutKeyboardFunc(key_dcluster);

  printf("------------------------------------\n");
  printf("Visualize Cluster bounding boxes!\n");
  printf("------------------------------------\n");
  printf("Use left window to select level to\n");
  printf("be drawn.\n");
  printf("Right window will show the chosen\n");
  printf("bounding boxes.\n");
  printf("If bounding boxes are 3-dimensional\n");
  printf("use mouse to chance perspective.\n");
  printf("To zoom in use '+' and '-' to zoom\n");
  printf("out.\n");
  printf("Keys 'a', 'd', 'w' and 's' translate\n");
  printf("perspective to the left, right, up\n");
  printf("and down.\n");
  printf("With 'n' perspective and zoom will\n");
  printf("be set to start values.\n");
  printf("Switch with 'p' between showing and\n");
  printf("hiding points.\n");
  printf("Fade in directions with '>'.\n");
  printf("Make snapshots of the content from one\n");
  printf("selcted window with 'space'. \n");
  printf("Leave visualization with 'esc'.\n");
  printf("------------------------------------\n");


  glutMainLoop();

  printf("Visualization has been closed.\n");
  del_activemovement(move_1);
  del_orthomovement(move_o);
}

void
visualize_block_bbox(pcblock b, pclustergeometry grc, pclustergeometry gcc,
		     bool c, int argc, char **argv)
{

  ct = b->rc;
  cs = b->cc;
  gt = grc;
  gs = gcc;
  bbl = b;
  color = c;
  version = 1;
  visu = 0;

  /* Checking if bounding boxes are given */
  if ((getdiam_2_cluster(ct) + getdiam_2_cluster(cs)) == 0.0) {
    printf("Drawing isn't possible, because bounding boxes are empty.");
    return;
  }

  if (gt == NULL || gs == NULL) {
    if (gt == NULL) {
      gs = gt;
    }
    if (gs == NULL) {
      gt = gs;
    }
    points = 0;
    printf("- - - - - - Caution! - - - - - -\n");
    printf("Missing clustergeometry for drawing points!\n");
  }

  if (ct->dim > 3) {
    printf("Row cluster dimension is with ");
    printf("%u", ct->dim);
    printf(" to hight for drawing!");

    return;
  }

  if (cs->dim > 3) {
    printf("Column cluster dimension is with ");
    printf("%u", cs->dim);
    printf(" to hight for drawing!");

    return;
  }

//    glutInit(&argc, argv);
  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
		GLUT_ACTION_GLUTMAINLOOP_RETURNS);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);

  glReadBuffer(GL_BACK);	/* Allowed us to read pixel out of the buffer */
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);	/* How to store the pixel */

  move_o = new_orthomovement();
  ortho = glutCreateWindow("Select level");
  glutPositionWindow(50, 300);
  glutReshapeWindow(400, 400);
  glutDisplayFunc(display_block);
  glutReshapeFunc(reshape_ortho);
  glutMouseFunc(mouse_block);
  glutKeyboardFunc(key_ortho);

  move_1 = new_activemovement(false);
  clust1 = glutCreateWindow("Visualize row cluster for certain block level");
  glutPositionWindow(550, 50);
  glutReshapeWindow(400, 400);
  glutDisplayFunc(display_two_cluster);
  glutReshapeFunc(reshape_cluster);
  glutMouseFunc(mouse_cluster);
  glutMotionFunc(motion_cluster);
  glutKeyboardFunc(key_cluster);

  move_2 = new_activemovement(false);
  clust2 =
    glutCreateWindow("Visualize column cluster for certain block level");
  glutPositionWindow(950, 50);
  glutReshapeWindow(400, 400);
  glutDisplayFunc(display_two_cluster);
  glutReshapeFunc(reshape_cluster);
  glutMouseFunc(mouse_cluster);
  glutMotionFunc(motion_cluster);
  glutKeyboardFunc(key_cluster);

  printf("------------------------------------\n");
  printf("Visualize Block bounding boxes!\n");
  printf("------------------------------------\n");
  printf("Use left window to select level to\n");
  printf("be drawn.\n");
  printf("Right window will show the chosen\n");
  printf("bounding boxes.\n");
  printf("If bounding boxes are 3-dimensional\n");
  printf("use mouse to chance perspective.\n");
  printf("To zoom in use '+' and '-' to zoom\n");
  printf("out.\n");
  printf("Keys 'a', 'd', 'w' and 's' translate\n");
  printf("perspective to the left, right, up\n");
  printf("and down.\n");
  printf("With 'n' perspective and zoom will\n");
  printf("be set to start values.\n");
  printf("Switch with 'p' between showing and\n");
  printf("hiding points.\n");
  printf("Make snapshots of the content from one\n");
  printf("selcted window with 'space'. \n");
  printf("Leave visualization with 'esc'.\n");
  printf("------------------------------------\n");

  glutMainLoop();

  printf("Visualization has been closed.\n");
  del_activemovement(move_1);
  del_activemovement(move_2);
  del_orthomovement(move_o);
}


void
visualize_block_certain_bbox(pcblock b, pclustergeometry grc,
			     pclustergeometry gcc, int argc, char **argv)
{

  visualize_block_level_certain_bbox(b, grc, gcc, 0, getdepth_block(b), argc,
				     argv);

}


void
visualize_block_level_certain_bbox(pcblock b, pclustergeometry grc,
				   pclustergeometry gcc, uint l1, uint l2,
				   int argc, char **argv)
{

  ct = b->rc;
  gt = grc;
  cs = b->cc;
  gs = gcc;
  bbl = b;
  version = 1;
  level = l1;
  level2 = l2;
  uint      depth;
  up = 0;

  assert((int) level >= 0);
  assert(level <= level2);

  /* Checking if bounding boxes are given */
  if ((getdiam_2_cluster(ct) + getdiam_2_cluster(cs)) == 0.0) {
    printf("Drawing isn't possible, because bounding boxes are empty.");
    return;
  }

  if (gt == NULL || gs == NULL) {
    if (gt == NULL) {
      gs = gt;
    }
    if (gs == NULL) {
      gt = gs;
    }
    points = 0;
    printf("- - - - - - Caution! - - - - - -\n");
    printf("Missing clustergeometry for drawing points!\n");
  }

  if (ct->dim > 3) {
    printf("Row cluster dimension is with ");
    printf("%u", ct->dim);
    printf(" to hight for drawing!");

    return;
  }

  if (cs->dim > 3) {
    printf("Column cluster dimension is with ");
    printf("%u", cs->dim);
    printf(" to hight for drawing!");

    return;
  }

  depth = getdepth_block(bbl);
  if (depth < level) {
    printf("Can't draw level ");
    printf("%u", level);
    printf(", because the block tree has only depth ");
    printf("%u.", depth);
    return;
  }

  if (depth < level2) {
    printf("Can't draw level ");
    printf("%u", level2);
    printf(", because the block tree has only depth ");
    printf("%u.", depth);
    return;
  }

//      glutInit(&argc, argv);
  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
		GLUT_ACTION_GLUTMAINLOOP_RETURNS);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);

  glReadBuffer(GL_BACK);	/* Allowed us to read pixel out of the buffer */
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);	/* How to store the pixel */

  move_o = new_orthomovement();
  ortho = glutCreateWindow("Projection blocks");
  glutPositionWindow(150, 100);
  glutReshapeWindow(600, 600);
  glutDisplayFunc(display_ortho);
  glutReshapeFunc(reshape_ortho);
  glutKeyboardFunc(key_ortho);
  glutMouseFunc(mouse_ortho);

  move_1 = new_activemovement(true);
  clust1 = glutCreateWindow("Visualize certain cluster");
  glutPositionWindow(800, 100);
  glutReshapeWindow(400, 400);
  glutDisplayFunc(display_certain_block);
  glutReshapeFunc(reshape_cluster);
  glutMouseFunc(mouse_cluster);
  glutMotionFunc(motion_cluster);
  glutKeyboardFunc(key_cluster);

  printf("------------------------------------\n");
  printf("Visualize Block bounding boxes!\n");
  printf("------------------------------------\n");
  printf("Use left window to select block to\n");
  printf("be drawn.\n");
  printf("Right window will show the chosen\n");
  printf("clusters with bounding boxes.\n");
  printf("Use rigth mouse button to switch to\n");
  printf("the father block.\n");
  printf("If bounding boxes are 3-dimensional\n");
  printf("use left mouse button to chance\n");
  printf("perspective.\n");
  printf("To zoom in use '+' and '-' to zoom\n");
  printf("out.\n");
  printf("Keys 'a', 'd', 'w' and 's' translate\n");
  printf("perspective to the left, right, up\n");
  printf("and down.\n");
  printf("With 'n' perspective and zoom will\n");
  printf("be set to start values.\n");
  printf("Switch with 'p' between showing and\n");
  printf("hiding points.\n");
  printf("Print bounding box values with 'b'.\n");
  printf("Make snapshots of the content from one\n");
  printf("selcted window with 'space'. \n");
  printf("Leave visualization with 'esc'.\n");
  printf("------------------------------------\n");


  glutMainLoop();

  printf("Visualization has been closed.\n");
  del_activemovement(move_1);
  del_orthomovement(move_o);
}

void
visualize_cluster_bbox(pccluster ti, pclustergeometry gti, bool c, int argc,
		       char **argv)
{

  ct = ti;
  cs = ct;
  gt = gti;
  color = c;
  version = 0;
  visu = 0;

  /* Checking if bounding boxes are given */
  if (getdiam_2_cluster(ct) == 0.0) {
    printf("Drawing isn't possible, because bounding boxes are empty.");
    return;
  }

  if (gt == NULL) {
    points = 0;
    printf("- - - - - - Caution! - - - - - -\n");
    printf("Missing clustergeometry for drawing points!\n");
  }


  if (ct->dim > 3) {
    printf("Dimension is with ");
    printf("%u", ct->dim);
    printf(" to hight for drawing!");

    return;
  }

//    glutInit(&argc, argv);
  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
		GLUT_ACTION_GLUTMAINLOOP_RETURNS);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);

  glReadBuffer(GL_BACK);	/* Allowed us to read pixel out of the buffer */
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);	/* How to store the pixel */

  move_o = new_orthomovement();
  ortho = glutCreateWindow("Select level");
  glutPositionWindow(50, 300);
  glutReshapeWindow(400, 400);
  glutDisplayFunc(display_block);
  glutReshapeFunc(reshape_ortho);
  glutMouseFunc(mouse_block);
  glutKeyboardFunc(key_ortho);

  move_1 = new_activemovement(false);
  clust1 = glutCreateWindow("Visualize cluster bounding boxes");
  glutPositionWindow(550, 100);
  glutReshapeWindow(600, 600);
  glutDisplayFunc(display_cluster);
  glutReshapeFunc(reshape_cluster);
  glutMouseFunc(mouse_cluster);
  glutMotionFunc(motion_cluster);
  glutKeyboardFunc(key_cluster);

  printf("------------------------------------\n");
  printf("Visualize Cluster bounding boxes!\n");
  printf("------------------------------------\n");
  printf("Use left window to select level to\n");
  printf("be drawn.\n");
  printf("Right window will show the chosen\n");
  printf("bounding boxes.\n");
  printf("If bounding boxes are 3-dimensional\n");
  printf("use mouse to chance perspective.\n");
  printf("To zoom in use '+' and '-' to zoom\n");
  printf("out.\n");
  printf("Keys 'a', 'd', 'w' and 's' translate\n");
  printf("perspective to the left, right, up\n");
  printf("and down.\n");
  printf("With 'n' perspective and zoom will\n");
  printf("be set to start values.\n");
  printf("Switch with 'p' between showing and\n");
  printf("hiding points.\n");
  printf("Make snapshots of the content from one\n");
  printf("selcted window with 'space'. \n");
  printf("Leave visualization with 'esc'.\n");
  printf("------------------------------------\n");


  glutMainLoop();

  printf("Visualization has been closed.\n");
  del_activemovement(move_1);
  del_orthomovement(move_o);
}

/************************************************
*	Visualize triangulations		*
*************************************************/

void
visualize_tri2d(pctri2d tri, int argc, char **argv)
{

  uint      i;

  mode = 0;
  points = 0;
  move_1 = new_activemovement(false);
  move_t = new_tablemovement();

  /* Set up triangulation */
  tr = new_tri2d(tri->vertices, tri->edges, tri->triangles);

  if (tr->vertices < 1) {
    printf("There are no vertices to visualize.\n");
    return;
  }

  for (i = 0; i < tr->edges; i++) {
    tr->e[i][0] = tri->e[i][0];
    tr->e[i][1] = tri->e[i][1];
    tr->eb[i] = tri->eb[i];
  }

  for (i = 0; i < tr->triangles; i++) {
    tr->t[i][0] = tri->t[i][0];
    tr->t[i][1] = tri->t[i][1];
    tr->t[i][2] = tri->t[i][2];
  }

  midpoint_triangle(tri);

  move_t->target = tr->triangles + 10;
//    glutInit(&argc, argv);  
  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
		GLUT_ACTION_GLUTMAINLOOP_RETURNS);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);

  glReadBuffer(GL_BACK);	/* Allowed us to read pixel out of the buffer */
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);	/* How to store the pixel */

  one = glutCreateWindow("Mesh");
  glutSetWindow(one);
  glutPositionWindow(50, 50);
  glutReshapeWindow(300, 300);
  glutDisplayFunc(display_tri2d);
  glutReshapeFunc(reshape_mesh);
  glutKeyboardFunc(key_mesh);
  glutMouseFunc(mouse_mesh);
  glutMotionFunc(motion_mesh);

  printf("----------------------------------\n");
  printf("Visualize a complete triangulation!\n");
  printf("----------------------------------\n");
  printf("Use mouse to chance perspective.\n");
  printf("To zoom in use '+' and '-' to zoom\n");
  printf("out.\n");
  printf("Keys 'a', 'd', 'w' and 's' translate\n");
  printf("perspective to the left, right, up\n");
  printf("and down.\n");
  printf("With 'n' perspective and zoom will\n");
  printf("be set to start values.\n");
  printf("Switch with 'p' between showing and\n");
  printf("hiding points.\n");
  printf("Fade in a simple coordinate system\n");
  printf("with 'c' and mark a single direction\n");
  printf("with 'x', 'y' or 'z', to fade in the\n");
  printf("origin use '0'.\n");
  printf("Make snapshots of the content from one\n");
  printf("selcted window with 'space'. \n");
  printf("Close window with 'esc'.\n");
  printf("----------------------------------\n");

  glutMainLoop();

  del_activemovement(move_1);
  del_tablemovement(move_t);
  del_tri2d(tr);

}

void
visualize_certain_triangle(pctri2d tri, int argc, char **argv)
{

  tr_orig = tri;
  uint      i;
  uint      tmp2;

  mode = 0;
  points = 0;
  move_1 = new_activemovement(false);
  move_t = new_tablemovement();

  /* Set up triangulation */
  tr = new_tri2d(tri->vertices, tri->edges, tri->triangles);

  if (tr->vertices < 1) {
    printf("There are no vertices to visualize.\n");
    return;
  }

  for (i = 0; i < tr->edges; i++) {
    tr->e[i][0] = tri->e[i][0];
    tr->e[i][1] = tri->e[i][1];
    tr->eb[i] = tri->eb[i];
  }

  for (i = 0; i < tr->triangles; i++) {
    tr->t[i][0] = tri->t[i][0];
    tr->t[i][1] = tri->t[i][1];
    tr->t[i][2] = tri->t[i][2];
  }

  midpoint_triangle(tri);

  /* Evaluate max_page */
  move_t->max_page = 0;
  tmp2 = tr->triangles;
  while ((int) (tmp2 - 10) > 0) {
    tmp2 = tmp2 - 10;
    move_t->max_page += 1;
  }

  move_t->target = tr->triangles + 10;
  move_t->page = 0;

//    glutInit(&argc, argv);
  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
		GLUT_ACTION_GLUTMAINLOOP_RETURNS);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);

  glReadBuffer(GL_BACK);	/* Allowed us to read pixel out of the buffer */
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);	/* How to store the pixel */

  table = glutCreateWindow("Select a certain triangle");
  glutSetWindow(table);
  glutPositionWindow(50, 50);
  glutReshapeWindow(500, 500);
  glutDisplayFunc(display_table_tri);
  glutReshapeFunc(reshape_table);
  glutMouseFunc(mouse_table_tri);
  glutKeyboardFunc(key_ortho);

  one = glutCreateWindow("Certain triangle");
  glutSetWindow(one);
  glutPositionWindow(600, 50);
  glutReshapeWindow(300, 300);
  glutDisplayFunc(display_certain_triangle);
  glutReshapeFunc(reshape_mesh);
  glutKeyboardFunc(key_mesh);
  glutMouseFunc(mouse_mesh);
  glutMotionFunc(motion_mesh);

  printf("----------------------------------\n");
  printf("Visualize a certain triangle out of\n");
  printf("a complete triangulation!\n");
  printf("----------------------------------\n");
  printf("Use left window to select a certain\n");
  printf("triangle to be drawn.\n");
  printf("Right window will show the chosen\n");
  printf("one.");
  printf("Use mouse to chance perspective.\n");
  printf("To zoom in use '+' and '-' to zoom\n");
  printf("out.\n");
  printf("Keys 'a', 'd', 'w' and 's' translate\n");
  printf("perspective to the left, right, up\n");
  printf("and down.\n");
  printf("With 'n' perspective and zoom will\n");
  printf("be set to start values.\n");
  printf("Fade in a simple coordinate system\n");
  printf("with 'c' and mark a single direction\n");
  printf("with 'x', 'y' or 'z', to fade in the\n");
  printf("origin use '0'.\n");
  printf("Make snapshots of the content from one\n");
  printf("selcted window with 'space'. \n");
  printf("Close window with 'esc'.\n");
  printf("----------------------------------\n");

  glutMainLoop();

  del_activemovement(move_1);
  del_tablemovement(move_t);
  del_tri2d(tr);

}

void
visualize_tet3d(pctet3d tet, int argc, char **argv)
{

  uint      i;

  mode = 0;
  points = 0;
  move_1 = new_activemovement(false);
  move_t = new_tablemovement();

  /* Set up triangulation */
  th = new_tet3d(tet->vertices, tet->edges, tet->faces, tet->tetrahedra);

  if (th->vertices < 1) {
    printf("There are no vertices to visualize.\n");
    return;
  }

  for (i = 0; i < th->edges; i++) {
    th->e[i][0] = tet->e[i][0];
    th->e[i][1] = tet->e[i][1];
    th->eb[i] = tet->eb[i];
  }

  for (i = 0; i < th->faces; i++) {
    th->f[i][0] = tet->f[i][0];
    th->f[i][1] = tet->f[i][1];
    th->f[i][2] = tet->f[i][2];
    th->fb[i] = tet->fb[i];
  }

  for (i = 0; i < th->tetrahedra; i++) {
    th->t[i][0] = tet->t[i][0];
    th->t[i][1] = tet->t[i][1];
    th->t[i][2] = tet->t[i][2];
    th->t[i][3] = tet->t[i][3];
  }

  midpoint_tetrahedra(tet);
  move_t->target = th->tetrahedra + 8;

//    glutInit(&argc, argv);       
  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
		GLUT_ACTION_GLUTMAINLOOP_RETURNS);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);

  glReadBuffer(GL_BACK);	/* Allowed us to read pixel out of the buffer */
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);	/* How to store the pixel */

  one = glutCreateWindow("Mesh");
  glutSetWindow(one);
  glutPositionWindow(50, 50);
  glutReshapeWindow(300, 300);
  glutDisplayFunc(display_tet3d);
  glutReshapeFunc(reshape_mesh);
  glutKeyboardFunc(key_mesh);
  glutMouseFunc(mouse_mesh);
  glutMotionFunc(motion_mesh);

  printf("----------------------------------\n");
  printf("Visualize a complete triangulation!\n");
  printf("----------------------------------\n");
  printf("Use mouse to chance perspective.\n");
  printf("To zoom in use '+' and '-' to zoom\n");
  printf("out.\n");
  printf("Keys 'a', 'd', 'w' and 's' translate\n");
  printf("perspective to the left, right, up\n");
  printf("and down.");
  printf("With 'n' perspective and zoom will\n");
  printf("be set to start values.\n");
  printf("Switch with 'p' between showing and\n");
  printf("hiding points.\n");
  printf("Fade in a simple coordinate system\n");
  printf("with 'c' and mark a single direction\n");
  printf("with 'x', 'y' or 'z', to fade in the\n");
  printf("origin use '0'.\n");
  printf("Make snapshots of the content from one\n");
  printf("selcted window with 'space'. \n");
  printf("Close window with 'esc'.\n");
  printf("----------------------------------\n");

  glutMainLoop();

  del_activemovement(move_1);
  del_tablemovement(move_t);
  del_tet3d(th);

}

void
visualize_certain_tetrahedra(pctet3d tet, int argc, char **argv)
{

  th_orig = tet;
  uint      i;
  uint      tmp2;

  mode = 0;
  points = 0;
  move_1 = new_activemovement(false);
  move_t = new_tablemovement();

  /* Set up triangulation */
  th = new_tet3d(tet->vertices, tet->edges, tet->faces, tet->tetrahedra);

  if (th->vertices < 1) {
    printf("There are no vertices to visualize.\n");
    return;
  }

  for (i = 0; i < th->edges; i++) {
    th->e[i][0] = tet->e[i][0];
    th->e[i][1] = tet->e[i][1];
    th->eb[i] = tet->eb[i];
  }

  for (i = 0; i < th->faces; i++) {
    th->f[i][0] = tet->f[i][0];
    th->f[i][1] = tet->f[i][1];
    th->f[i][2] = tet->f[i][2];
    th->fb[i] = tet->fb[i];
  }

  for (i = 0; i < th->tetrahedra; i++) {
    th->t[i][0] = tet->t[i][0];
    th->t[i][1] = tet->t[i][1];
    th->t[i][2] = tet->t[i][2];
    th->t[i][3] = tet->t[i][3];
  }

  midpoint_tetrahedra(tet);
  move_t->target = th->tetrahedra + 8;

  /* Evaluate max_page */
  move_t->max_page = 0;
  tmp2 = th->tetrahedra;
  while ((int) (tmp2 - 8) > 0) {
    tmp2 = tmp2 - 8;
    move_t->max_page += 1;
  }
  move_t->page = 0;
  mode = 0;

//    glutInit(&argc, argv);
  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
		GLUT_ACTION_GLUTMAINLOOP_RETURNS);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);

  glReadBuffer(GL_BACK);	/* Allowed us to read pixel out of the buffer */
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);	/* How to store the pixel */

  table = glutCreateWindow("Select a certain tetrahedra");
  glutSetWindow(table);
  glutPositionWindow(50, 50);
  glutReshapeWindow(600, 600);
  glutDisplayFunc(display_table_tet);
  glutReshapeFunc(reshape_table);
  glutMouseFunc(mouse_table_tet);
  glutKeyboardFunc(key_ortho);

  one = glutCreateWindow("Certain triangle");
  glutSetWindow(one);
  glutPositionWindow(700, 50);
  glutReshapeWindow(300, 300);
  glutDisplayFunc(display_certain_tetrahedra);
  glutReshapeFunc(reshape_mesh);
  glutKeyboardFunc(key_mesh);
  glutMouseFunc(mouse_mesh);
  glutMotionFunc(motion_mesh);


  printf("----------------------------------\n");
  printf("Visualize a certain tetrahedra out\n");
  printf("of a complete triangulation!\n");
  printf("----------------------------------\n");
  printf("Use left window to select a certain\n");
  printf("tetrahedra to be drawn.\n");
  printf("Right window will show the chosen\n");
  printf("one.");
  printf("Use mouse to chance perspective.\n");
  printf("To zoom in use '+' and '-' to zoom\n");
  printf("out.\n");
  printf("Keys 'a', 'd', 'w' and 's' translate\n");
  printf("perspective to the left, right, up\n");
  printf("and down.\n");
  printf("With 'n' perspective and zoom will\n");
  printf("be set to start values.\n");
  printf("Fade in a simple coordinate system\n");
  printf("with 'c' and mark a single direction\n");
  printf("with 'x', 'y' or 'z', to fade in the\n");
  printf("origin use '0'.\n");
  printf("Make snapshots of the content from one\n");
  printf("selcted window with 'space'. \n");
  printf("Close window with 'esc'.\n");
  printf("----------------------------------\n");


  glutMainLoop();
  del_tet3d(th);
  del_activemovement(move_1);
  del_tablemovement(move_t);

}

void
visualize_surface3d(pcsurface3d sur, int argc, char **argv)
{

  uint      i;

  mode = 0;
  points = 0;
  move_1 = new_activemovement(false);
  move_t = new_tablemovement();


  /* Set up surface */
  gr = new_surface3d(sur->vertices, sur->edges, sur->triangles);

  if (gr->vertices < 1) {
    printf("There are no vertices to visualize.\n");
    return;
  }

  for (i = 0; i < gr->edges; i++) {
    gr->e[i][0] = sur->e[i][0];
    gr->e[i][1] = sur->e[i][1];
  }

  for (i = 0; i < gr->triangles; i++) {
    gr->t[i][0] = sur->t[i][0];
    gr->t[i][1] = sur->t[i][1];
    gr->t[i][2] = sur->t[i][2];
  }

  midpoint_surface(sur);
  move_t->target = gr->triangles + 8;

//   glutInit(&argc, argv);
  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
		GLUT_ACTION_GLUTMAINLOOP_RETURNS);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);

  glReadBuffer(GL_BACK);	/* Allowed us to read pixel out of the buffer */
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);	/* How to store the pixel */

  one = glutCreateWindow("Mesh");
  glutSetWindow(one);
  glutPositionWindow(50, 50);
  glutReshapeWindow(300, 300);
  glutDisplayFunc(display_surface);
  glutReshapeFunc(reshape_mesh);
  glutKeyboardFunc(key_mesh);
  glutMouseFunc(mouse_mesh);
  glutMotionFunc(motion_mesh);

  printf("----------------------------------\n");
  printf("Visualize a complete triangulation!\n");
  printf("----------------------------------\n");
  printf("Use mouse to chance perspective.\n");
  printf("To zoom in use '+' and '-' to zoom\n");
  printf("out.\n");
  printf("Keys 'a', 'd', 'w' and 's' translate\n");
  printf("perspective to the left, right, up\n");
  printf("and down.\n");
  printf("With 'n' perspective and zoom will\n");
  printf("be set to start values.\n");
  printf("Switch with 'p' between showing and\n");
  printf("hiding points.\n");
  printf("Fade in a simple coordinate system\n");
  printf("with 'c' and mark a single direction\n");
  printf("with 'x', 'y' or 'z', to fade in the\n");
  printf("origin use '0'.\n");
  printf("Make snapshots of the content from one\n");
  printf("selcted window with 'space'. \n");
  printf("Close window with 'esc'.\n");
  printf("----------------------------------\n");

  glutMainLoop();

  del_activemovement(move_1);
  del_tablemovement(move_t);
  del_surface3d(gr);

}

void
visualize_certain_surface_triangle(pcsurface3d sur, int argc, char **argv)
{

  gr_orig = sur;
  uint      i;
  uint      tmp2;

  mode = 0;
  points = 0;
  move_1 = new_activemovement(false);
  move_t = new_tablemovement();

  /* Set up surface */
  gr = new_surface3d(sur->vertices, sur->edges, sur->triangles);

  if (gr->vertices < 1) {
    printf("There are no vertices to visualize.\n");
    return;
  }

  for (i = 0; i < gr->edges; i++) {
    gr->e[i][0] = sur->e[i][0];
    gr->e[i][1] = sur->e[i][1];
  }

  for (i = 0; i < gr->triangles; i++) {
    gr->t[i][0] = sur->t[i][0];
    gr->t[i][1] = sur->t[i][1];
    gr->t[i][2] = sur->t[i][2];
  }

  midpoint_surface(sur);
  move_t->target = gr->triangles + 8;

  /* Evaluate max_page */
  move_t->max_page = 0;
  tmp2 = gr->triangles;
  while ((int) (tmp2 - 8) > 0) {
    tmp2 = tmp2 - 8;
    move_t->max_page += 1;
  }
  move_t->page = 0;
  mode = 0;

//  glutInit(&argc, argv);
  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
		GLUT_ACTION_GLUTMAINLOOP_RETURNS);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);

  glReadBuffer(GL_BACK);	/* Allowed us to read pixel out of the buffer */
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);	/* How to store the pixel */

  table = glutCreateWindow("Select a certain triangle");
  glutSetWindow(table);
  glutPositionWindow(50, 50);
  glutReshapeWindow(500, 500);
  glutDisplayFunc(display_table_sur);
  glutReshapeFunc(reshape_table);
  glutMouseFunc(mouse_table_sur);
  glutKeyboardFunc(key_ortho);

  one = glutCreateWindow("Certain triangle");
  glutSetWindow(one);
  glutPositionWindow(600, 50);
  glutReshapeWindow(300, 300);
  glutDisplayFunc(display_certain_surface_element);
  glutReshapeFunc(reshape_mesh);
  glutKeyboardFunc(key_mesh);
  glutMouseFunc(mouse_mesh);
  glutMotionFunc(motion_mesh);


  printf("----------------------------------\n");
  printf("Visualize a certain triangle out\n");
  printf("of a complete surface triangulation!\n");
  printf("----------------------------------\n");
  printf("Use left window to select a certain\n");
  printf("triangle to be drawn\n.");
  printf("Right window will show the chosen\n");
  printf("one.");
  printf("Use mouse to chance perspective.\n");
  printf("To zoom in use '+' and '-' to zoom\n");
  printf("out.\n");
  printf("Keys 'a', 'd', 'w' and 's' translate\n");
  printf("perspective to the left, right, up\n");
  printf("and down.\n");
  printf("With 'n' perspective and zoom will\n");
  printf("be set to start values.\n");
  printf("Fade in a simple coordinate system\n");
  printf("with 'c' and mark a single direction\n");
  printf("with 'x', 'y' or 'z', to fade in the\n");
  printf("origin use '0'.\n");
  printf("Make snapshots of the content from one\n");
  printf("selcted window with 'space'. \n");
  printf("Close window with 'esc'.\n");
  printf("----------------------------------\n");


  glutMainLoop();

  del_activemovement(move_1);
  del_tablemovement(move_t);
  del_surface3d(gr);
}


void
visualize_boundaryvalue_surface_triangle(pcsurface3d sur, pcavector val,
					 int argc, char **argv)
{

  gr_orig = sur;
  real      scale;
  real      max, min;
  uint      i;
  uint      type;

  mode = 0;
  points = 0;
  move_1 = new_activemovement(false);

  /* Set up surface */
  gr = new_surface3d(sur->vertices, sur->edges, sur->triangles);

  if (gr->vertices < 1) {
    printf("There are no vertices to visualize.\n");
    return;
  }

  for (i = 0; i < gr->edges; i++) {
    gr->e[i][0] = sur->e[i][0];
    gr->e[i][1] = sur->e[i][1];
  }

  for (i = 0; i < gr->triangles; i++) {
    gr->t[i][0] = sur->t[i][0];
    gr->t[i][1] = sur->t[i][1];
    gr->t[i][2] = sur->t[i][2];
  }

  midpoint_surface(sur);

  /* Find out if the dimension of the solution vector fits to the number of vertices, edges or triangles */

  type = 0;
  if (val->dim != gr->triangles) {
    type += 1;
    if (val->dim != gr->edges) {
      if (val->dim != gr->vertices) {
	printf("The given vector dimension doesn't match to the surface!\n");
	return;
      }
      type += 1;
    }
  }

  /* Scaling solution vector */
  scale = normmax_avactor(val);
  sol = new_solution(val->dim, 0, false);

  max = 0.0;
  min = 0.0;
  for (i = 0; i < val->dim; i++) {
    max = (REAL(val->v[i]) > max ? REAL(val->v[i]) : max);
    min = (REAL(val->v[i]) < min ? REAL(val->v[i]) : min);
    sol->values->v[i] = REAL(val->v[i]) / scale;
  }
  sol->max[0][0] = max;
  sol->min[0][0] = min;


  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
		GLUT_ACTION_GLUTMAINLOOP_RETURNS);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);

  glReadBuffer(GL_BACK);	/* Allowed us to read pixel out of the buffer */
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);	/* How to store the pixel */

  one = glutCreateWindow("Values on surface");
  glutSetWindow(one);
  glutPositionWindow(50, 50);
  glutReshapeWindow(300, 300);
  if (type == 0) {
    glutDisplayFunc(display_surface_value_triangle);
  }
  else {
    if (type == 1) {
      glutDisplayFunc(display_surface_value_edges);
    }
    else {
      glutDisplayFunc(display_surface_value_vertices);
    }
  }
  glutReshapeFunc(reshape_mesh);
  glutKeyboardFunc(key_mesh);
  glutMouseFunc(mouse_mesh);
  glutMotionFunc(motion_mesh);

  table = glutCreateWindow("Legend");
  glutSetWindow(table);
  glutPositionWindow(400, 50);
  glutReshapeWindow(300, 300);
  glutDisplayFunc(display_legend_surface);
  glutReshapeFunc(reshape_legend);
  glutKeyboardFunc(key_ortho);


  printf("----------------------------------\n");
  printf("Visualize values on surface!\n");
  printf("----------------------------------\n");
  printf("Left windows shows surface with \n");
  printf("solution, right one the legend.\n");
  printf("Use mouse to chance perspective.\n");
  printf("To zoom in use '+' and '-' to zoom\n");
  printf("out.\n");
  printf("Keys 'a', 'd', 'w' and 's' translate\n");
  printf("perspective to the left, right, up\n");
  printf("and down.\n");
  printf("With 'n' perspective and zoom will\n");
  printf("be set to start values.\n");
  printf("Fade in a simple coordinate system\n");
  printf("with 'c' and mark a single direction\n");
  printf("with 'x', 'y' or 'z', to fade in the\n");
  printf("origin use '0'.\n");
  printf("Make snapshots of the content from one\n");
  printf("selcted window with 'space'. \n");
  printf("Close window with 'esc'.\n");
  printf("----------------------------------\n");

  glutMainLoop();

  del_activemovement(move_1);
  del_surface3d(gr);
  del_solution(sol);
}

void
visualize_helmholtz_solution_surface(pcsurface3d sur, int argc, char **argv)
{

  gr_orig = sur;
  real      scale;
  uint      i, t;
  uint      type;

  mode = 0;
  points = 0;
  move_1 = new_activemovement(false);

  /* Set up surface */
  gr = new_surface3d(sur->vertices, sur->edges, sur->triangles);

  if (gr->vertices < 1) {
    printf("There are no vertices to visualize.\n");
    return;
  }

  for (i = 0; i < gr->edges; i++) {
    gr->e[i][0] = sur->e[i][0];
    gr->e[i][1] = sur->e[i][1];
  }

  for (i = 0; i < gr->triangles; i++) {
    gr->t[i][0] = sur->t[i][0];
    gr->t[i][1] = sur->t[i][1];
    gr->t[i][2] = sur->t[i][2];
  }

  scale = midpoint_surface(sur);
  sol = new_solution(0, argc - 1, false);

  printf("reading files\n");
  fflush(stdout);
  /* Read files for solutions */
  if (argc == 1) {
    sol->plane = false;
    type = 0;
  }
  else {
    printf("Reading files\n");
    type = 1;
    for (i = 1; i < argc; i++) {
      t = read_file(argv[i], scale, i - 1, sol, false);
      if (t != 1) {
	sol->plane = false;
	type = 0;
	return;
      }
    }
  }


  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
		GLUT_ACTION_GLUTMAINLOOP_RETURNS);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);

  glReadBuffer(GL_BACK);	/* Allowed us to read pixel out of the buffer */
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);	/* How to store the pixel */

  one = glutCreateWindow("Values on surface");
  glutSetWindow(one);
  glutPositionWindow(50, 50);
  glutReshapeWindow(300, 300);
  glutDisplayFunc(display_solution_surface);
  glutReshapeFunc(reshape_mesh);
  glutKeyboardFunc(key_solution);
  glutMouseFunc(mouse_mesh);
  glutMotionFunc(motion_mesh);

  if (type > 0) {
    table = glutCreateWindow("Legend");
    glutSetWindow(table);
    glutPositionWindow(400, 50);
    glutReshapeWindow(300, 300);
    glutDisplayFunc(display_legend_surface);
    glutReshapeFunc(reshape_legend);
    glutKeyboardFunc(key_ortho);

    printf("----------------------------------\n");
    printf("Visualize values on surface!\n");
    printf("----------------------------------\n");
    printf("Left windows shows surface with \n");
    printf("solution, right one the legend.\n");
    printf("Use mouse to chance perspective.\n");
    printf("To zoom in use '+' and '-' to zoom\n");
    printf("out.\n");
    printf("Keys 'a', 'd', 'w' and 's' translate\n");
    printf("perspective to the left, right, up\n");
    printf("and down.\n");
    printf("With 'n' perspective and zoom will\n");
    printf("be set to start values.\n");
    printf("Fade in a simple coordinate system\n");
    printf("with 'c' and mark a single direction\n");
    printf("with 'x', 'y' or 'z', to fade in the\n");
    printf("origin use '0'.\n");
    printf("Change visibility of surface with 'm'\n");
    printf("and solution with 'p'.\n");
    printf("If more than one file is available, \n");
    printf("change visualized file with '.'.\n");
    printf("Make snapshots of the content from one\n");
    printf("selcted window with 'space'. \n");
    printf("Close window with 'esc'.\n");
    printf("----------------------------------\n");
  }
  else {
    printf("----------------------------------\n");
    printf("Visualize values on surface!\n");
    printf("----------------------------------\n");
    printf("Missing files including data for \n");
    printf("solutions on planes!\n");
    printf("For more information take a look \n");
    printf("at the visualize.h file. \n");
  }



  glutMainLoop();

  del_activemovement(move_1);
  del_surface3d(gr);
  del_solution(sol);
}

void
animate_helmholtz_solution(pcsurface3d sur, int argc, char **argv)
{

  gr_orig = sur;
  real      scale, modify;
  uint      i, j, t;
  uint      type;

  mode = 0;
  points = 0;
  move_1 = new_activemovement(false);

  /* Set up surface */
  gr = new_surface3d(sur->vertices, sur->edges, sur->triangles);

  if (gr->vertices < 1) {
    printf("There are no vertices to visualize.\n");
    return;
  }

  for (i = 0; i < gr->edges; i++) {
    gr->e[i][0] = sur->e[i][0];
    gr->e[i][1] = sur->e[i][1];
  }

  for (i = 0; i < gr->triangles; i++) {
    gr->t[i][0] = sur->t[i][0];
    gr->t[i][1] = sur->t[i][1];
    gr->t[i][2] = sur->t[i][2];
  }

  scale = midpoint_surface(sur);
  sol = new_solution(0, argc - 1, true);

  /* Read files for solutions */
  if (argc == 1) {
    type = 0;

    printf("----------------------------------\n");
    printf("Animate Helmholtz solution in time!\n");
    printf("----------------------------------\n");
    printf("Missing files including data for \n");
    printf("solutions on planes!\n");
    printf("For more information take a look \n");
    printf("at the visualize.h file. \n");

    return;
  }
  else {
    printf("Reading files\n");
    type = 1;
    for (i = 1; i < argc; i++) {
      t = read_file(argv[i], scale, i - 1, sol, true);
      /* Modify max and min values for a animated solution */
      for (j = 0; j < 3; j++) {
	modify = (sol->max[i - 1][j] - sol->min[i - 1][j]) * 0.05;
	sol->max[i - 1][j] = sol->max[i - 1][j] + modify;
	sol->min[i - 1][j] = sol->min[i - 1][j] - modify;
      }
      if (t != 1) {
	sol->plane = false;
	type = 0;
	return;
      }
    }
  }

  /* First ignore the model */
  sol->model = false;

  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
		GLUT_ACTION_GLUTMAINLOOP_RETURNS);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);

  glReadBuffer(GL_BACK);	/* Allowed us to read pixel out of the buffer */
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);	/* How to store the pixel */

  one = glutCreateWindow("Values on surface");
  glutSetWindow(one);
  glutPositionWindow(50, 50);
  glutReshapeWindow(300, 300);
  glutDisplayFunc(display_solution_surface);
  glutReshapeFunc(reshape_mesh);
  glutKeyboardFunc(key_solution);
  glutMouseFunc(mouse_mesh);
  glutMotionFunc(motion_mesh);
  glutTimerFunc(50, timer_helmholtz, 0);

  if (type > 0) {
    table = glutCreateWindow("Legend");
    glutSetWindow(table);
    glutPositionWindow(400, 50);
    glutReshapeWindow(300, 300);
    glutDisplayFunc(display_legend_surface);
    glutReshapeFunc(reshape_legend);
    glutKeyboardFunc(key_ortho);

    printf("----------------------------------\n");
    printf("Animate Helmholtz solution in time!\n");
    printf("----------------------------------\n");
    printf("Left windows shows surface with \n");
    printf("solution, right one the legend.\n");
    printf("Use mouse to chance perspective.\n");
    printf("Start and stop animation with 'enter',\n");
    printf("faster animation '>' and slower '<' .\n");
    printf("To zoom in use '+' and '-' to zoom\n");
    printf("out.\n");
    printf("Keys 'a', 'd', 'w' and 's' translate\n");
    printf("perspective to the left, right, up\n");
    printf("and down.\n");
    printf("With 'n' perspective, zoom and\n");
    printf("animation time set to start values.\n");
    printf("Fade in a simple coordinate system\n");
    printf("with 'c' and mark a single direction\n");
    printf("with 'x', 'y' or 'z', to fade in the\n");
    printf("origin use '0'.\n");
    printf("Change visibility of surface with 'm'\n");
    printf("and solution with 'p'.\n");
    printf("If more than one file is available, \n");
    printf("change visualized file with '.'.\n");
    printf("Make snapshots of the content from one\n");
    printf("selcted window with 'space'. \n");
    printf("Close window with 'esc'.\n");
    printf("----------------------------------\n");
  }

  glutMainLoop();

  del_activemovement(move_1);
  del_surface3d(gr);
  del_solution(sol);
}


#else

void
visualize_dblock_bbox(pcdblock b, pclustergeometry grc, pclustergeometry gcc,
		      bool c, int argc, char **argv)
{

  (void) b;
  (void) c;

  fprintf(stderr, "Sorry, GLUT is not available.\n");
}

void
visualize_dblock_certain_bbox(pcdblock b, pclustergeometry grc,
			      pclustergeometry gcc, int argc, char **argv)
{

  (void) b;
  (void) grc;
  (void) gcc;

  fprintf(stderr, "Sorry, GLUT is not available.\n");
}

void
visualize_dblock_level_certain_bbox(pcdblock b, pclustergeometry grc,
				    pclustergeometry gcc, uint l1, uint l2,
				    int argc, char **argv)
{

  (void) b;
  (void) grc;
  (void) gcc;
  (void) l1;
  (void) l2;

  fprintf(stderr, "Sorry, GLUT is not available.\n");
}

void
visualize_dcluster_bbox(pcdcluster ti, pclustergeometry gti, bool c, int argc,
			char **argv)
{

  (void) gti;
  (void) ti;
  (void) c;

  fprintf(stderr, "Sorry, GLUT is not available.\n");
}


void
visualize_block_bbox(pcblock b, pclustergeometry grc, pclustergeometry gcc,
		     bool c, int argc, char **argv)
{

  (void) b;
  (void) c;

  fprintf(stderr, "Sorry, GLUT is not available.\n");
}

void
visualize_block_certain_bbox(pcblock b, pclustergeometry grc,
			     pclustergeometry gcc, int argc, char **argv)
{

  (void) b;
  (void) grc;
  (void) gcc;

  fprintf(stderr, "Sorry, GLUT is not available.\n");
}

void
visualize_block_level_certain_bbox(pcblock b, pclustergeometry grc,
				   pclustergeometry gcc, uint l1, uint l2,
				   int argc, char **argv)
{

  (void) b;
  (void) grc;
  (void) gcc;
  (void) l1;
  (void) l2;

  fprintf(stderr, "Sorry, GLUT is not available.\n");
}

void
visualize_cluster_bbox(pccluster ti, pclustergeometry gti, bool c, int argc,
		       char **argv)
{

  (void) gti;
  (void) ti;
  (void) c;

  fprintf(stderr, "Sorry, GLUT is not available.\n");
}

void
visualize_tri2d(pctri2d tri, int argc, char **argv)
{

  fprintf(stderr, "Sorry, GLUT is not available.\n");
}

void
visualize_tet3d(pctet3d tet, int argc, char **argv)
{

  fprintf(stderr, "Sorry, GLUT is not available.\n");
}

void
visualize_certain_triangle(pctri2d tri, int argc, char **argv)
{

  fprintf(stderr, "Sorry, GLUT is not available.\n");
}

void
visualize_certain_tetrahedra(pctet3d tet, int argc, char **argv)
{

  fprintf(stderr, "Sorry, GLUT is not available.\n");
}

void
visualize_surface3d(pcsurface3d sur, int argc, char **argv)
{

  fprintf(stderr, "Sorry, GLUT is not available.\n");
}

void
visualize_certain_surface_triangle(pcsurface3d sur, int argc, char **argv)
{

  fprintf(stderr, "Sorry, GLUT is not available.\n");
}

void
visualize_boundaryvalue_surface_triangle(pcsurface3d sur, pcavector val,
					 int argc, char **argv)
{

  fprintf(stderr, "Sorry, GLUT is not available.\n");
}

void
visualize_helmholtz_solution_surface(pcsurface3d sur, int argc, char **argv)
{

  fprintf(stderr, "Sorry, GLUT is not available.\n");
}

void
animate_helmholtz_solution(pcsurface3d sur, int argc, char **argv)
{

  fprintf(stderr, "Sorry, GLUT is not available.\n");
}

#endif
