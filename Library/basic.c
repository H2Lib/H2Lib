/* ------------------------------------------------------------
 This is the file "basic.c" of the H2Lib package.
 All rights reserved, Steffen Boerm 2009
 ------------------------------------------------------------ */

#include "basic.h"

#include <stdio.h>
#ifdef WIN32
#include <Windows.h>
#include <MMSystem.h>
#pragma comment(lib, "winmm")
#else
#include <sys/times.h>
#include <unistd.h>
#endif

#ifdef USE_CAIRO
#include <cairo/cairo.h>
#include <cairo/cairo-pdf.h>
#endif

#ifdef USE_OPENMP
#include <omp.h>
#endif

#ifdef USE_FREEGLUT
#include <GL/freeglut.h>
#endif

#ifdef USE_GTK3
#include <gtk/gtk.h>
#endif

int       max_pardepth = 0;

#ifdef TRACE_MEMORY
size_t    current_memory = 0;
#endif

/* ------------------------------------------------------------
 Set up the library
 ------------------------------------------------------------ */

void
init_h2lib(int *argc, char ***argv)
{
  (void) argc;
  (void) argv;

#ifdef USE_OPENMP
  char     *env;
  int       i, j;

  if (omp_in_parallel()) {
    fprintf(stderr, "Calling init_h2lib in parallel section is forbidden.\n");
    abort();
  }
  omp_set_nested(1);

  i = omp_get_num_procs();
  j = 0;
  while (i > 1) {
    i /= 2;
    j++;
  }
  max_pardepth = j + 1;

  env = getenv("H2_PARDEPTH");
  if (env)
    sscanf(env, "%d", &max_pardepth);
#else
  max_pardepth = 0;
#endif

#ifdef USE_FREEGLUT
  glutInit(argc, *argv);
#endif

#ifdef USE_GTK3
  gtk_disable_setlocale();
  gtk_init(argc, argv);
#endif
}

void
uninit_h2lib()
{

}

/* ------------------------------------------------------------
 Memory management
 ------------------------------------------------------------ */

#ifdef TRACE_MEMORY

#ifdef USE_SIMD
#define ALLOC_OFFSET (VALIGN)
#else
#define ALLOC_OFFSET (sizeof(size_t))
#endif

void     *
_h2_allocmem(size_t sz, const char *filename, int line)
{
  void     *ptr;

  sz += ALLOC_OFFSET;

  ptr = h2_malloc(sz);
  if (ptr == NULL && sz > 0) {
    (void) fprintf(stderr, "Memory allocation of %lu bytes failed in %s:%d\n",
		   (unsigned long) sz, filename, line);
    abort();
  }

  *((size_t *) ptr) = sz;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  current_memory += sz;

  return ptr + ALLOC_OFFSET;
}

uint     *
_h2_allocuint(size_t sz, const char *filename, int line)
{
  uint     *ptr;
  size_t    dsz;

  dsz = sizeof(uint) * sz + ALLOC_OFFSET;
  if ((dsz - ALLOC_OFFSET) / sizeof(uint) != sz) {
    (void) fprintf(stderr, "Integer overflow in vector allocation in %s:%d\n",
		   filename, line);
    abort();
  }

  ptr = (uint *) h2_malloc(dsz);
  if (ptr == NULL && dsz > 0) {
    (void) fprintf(stderr,
		   "Vector allocation of %lu entries failed in %s:%d\n",
		   (unsigned long) sz, filename, line);
    abort();
  }

  *((size_t *) ptr) = dsz;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  current_memory += dsz;

  return (uint *) ((void *) ptr + ALLOC_OFFSET);
}

real     *
_h2_allocreal(size_t sz, const char *filename, int line)
{
  real     *ptr;
  size_t    dsz;

  dsz = sizeof(real) * sz + ALLOC_OFFSET;
  if ((dsz - ALLOC_OFFSET) / sizeof(real) != sz) {
    (void) fprintf(stderr, "Integer overflow in vector allocation in %s:%d\n",
		   filename, line);
    abort();
  }

  ptr = (real *) h2_malloc(dsz);
  if (ptr == NULL && dsz > 0) {
    (void) fprintf(stderr,
		   "Vector allocation of %lu entries failed in %s:%d\n",
		   (unsigned long) sz, filename, line);
    abort();
  }

  *((size_t *) ptr) = dsz;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  current_memory += dsz;

  return (real *) ((void *) ptr + ALLOC_OFFSET);
}

field    *
_h2_allocfield(size_t sz, const char *filename, int line)
{
  field    *ptr;
  size_t    dsz;

  dsz = sizeof(field) * sz + ALLOC_OFFSET;
  if ((dsz - ALLOC_OFFSET) / sizeof(field) != sz) {
    (void) fprintf(stderr, "Integer overflow in vector allocation in %s:%d\n",
		   filename, line);
    abort();
  }

  ptr = (field *) h2_malloc(dsz);
  if (ptr == NULL && dsz > 0) {
    (void) fprintf(stderr,
		   "Vector allocation of %lu entries failed in %s:%d\n",
		   (unsigned long) sz, filename, line);
    abort();
  }

  *((size_t *) ptr) = dsz;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  current_memory += dsz;

  return (field *) ((void *) ptr + ALLOC_OFFSET);
}

field    *
_h2_allocmatrix(size_t rows, size_t cols, const char *filename, int line)
{
  field    *ptr;
  size_t    dsz;

  dsz = sizeof(field) * rows * cols + ALLOC_OFFSET;
  if ((dsz - ALLOC_OFFSET) / sizeof(field) != rows * cols) {
    (void) fprintf(stderr, "Integer overflow in matrix allocation in %s:%d\n",
		   filename, line);
    abort();
  }

  ptr = (field *) h2_malloc(dsz);
  if (ptr == NULL && dsz > 0) {
    (void) fprintf(stderr,
		   "Matrix allocation with %lu rows and %lu columns failed in %s:%d\n",
		   (unsigned long) rows, (unsigned long) cols, filename,
		   line);
    abort();
  }

  *((size_t *) ptr) = dsz;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  current_memory += dsz;

  return (field *) ((void *) ptr + ALLOC_OFFSET);
}

void
freemem(void *ptr)
{
  size_t    sz;

  if (ptr != NULL) {
    sz = *(size_t *) (ptr - ALLOC_OFFSET);

    h2_free(ptr - ALLOC_OFFSET);

    assert(current_memory >= sz);

#ifdef USE_OPENMP
#pragma omp atomic
#endif
    current_memory -= sz;
  }
}

size_t
get_current_memory()
{
  return current_memory;
}
#else
void     *
_h2_allocmem(size_t sz, const char *filename, int line)
{
  void     *ptr;

  ptr = h2_malloc(sz);
  if (ptr == NULL && sz > 0) {
    (void) fprintf(stderr, "Memory allocation of %lu bytes failed in %s:%d\n",
		   (unsigned long) sz, filename, line);
    abort();
  }

  return ptr;
}

uint     *
_h2_allocuint(size_t sz, const char *filename, int line)
{
  uint     *ptr;
  size_t    dsz;

  dsz = sizeof(uint) * sz;
  if (dsz / sizeof(uint) != sz) {
    (void) fprintf(stderr, "Integer overflow in vector allocation in %s:%d\n",
		   filename, line);
    abort();
  }

  ptr = (uint *) h2_malloc(dsz);
  if (ptr == NULL && dsz > 0) {
    (void) fprintf(stderr,
		   "Vector allocation of %lu entries failed in %s:%d\n",
		   (unsigned long) sz, filename, line);
    abort();
  }

  return ptr;
}

real     *
_h2_allocreal(size_t sz, const char *filename, int line)
{
  real     *ptr;
  size_t    dsz;

  dsz = sizeof(real) * sz;
  if (dsz / sizeof(real) != sz) {
    (void) fprintf(stderr, "Integer overflow in vector allocation in %s:%d\n",
		   filename, line);
    abort();
  }

  ptr = (real *) h2_malloc(dsz);
  if (ptr == NULL && dsz > 0) {
    (void) fprintf(stderr,
		   "Vector allocation of %lu entries failed in %s:%d\n",
		   (unsigned long) sz, filename, line);
    abort();
  }

  return ptr;
}

field    *
_h2_allocfield(size_t sz, const char *filename, int line)
{
  field    *ptr;
  size_t    dsz;

  dsz = sizeof(field) * sz;
  if (dsz / sizeof(field) != sz) {
    (void) fprintf(stderr, "Integer overflow in vector allocation in %s:%d\n",
		   filename, line);
    abort();
  }

  ptr = (field *) h2_malloc(dsz);
  if (ptr == NULL && dsz > 0) {
    (void) fprintf(stderr,
		   "Vector allocation of %lu entries failed in %s:%d\n",
		   (unsigned long) sz, filename, line);
    abort();
  }

  return ptr;
}

field    *
_h2_allocmatrix(size_t rows, size_t cols, const char *filename, int line)
{
  field    *ptr;
  size_t    dsz;

  dsz = sizeof(field) * rows * cols;
  if (dsz / sizeof(field) != rows * cols) {
    (void) fprintf(stderr, "Integer overflow in matrix allocation in %s:%d\n",
		   filename, line);
    abort();
  }

  ptr = (field *) h2_malloc(dsz);
  if (ptr == NULL && dsz > 0) {
    (void) fprintf(stderr,
		   "Matrix allocation with %lu rows and %lu columns failed in %s:%d\n",
		   (unsigned long) rows, (unsigned long) cols, filename,
		   line);
    abort();
  }

  return ptr;
}

void
freemem(void *ptr)
{
  if (ptr != NULL) {
    h2_free(ptr);
  }
}

size_t
get_current_memory()
{
  fprintf(stderr, "Memory tracing not available!\n"
	  "Please activate by setting TRACE_MEMORY compiler flag!\n");
  return 0;
}
#endif

/* ------------------------------------------------------------
 Sorting
 ------------------------------------------------------------ */

static void
heap_down(uint root, uint n, bool leq(uint, uint, void *),
	  void swap(uint, uint, void *), void *data)
{
  uint      child;

  child = 2 * root + 1;
  while (child < n) {
    /* Find larger child */
    if (child + 1 < n && leq(child, child + 1, data))
      child++;

    /* Check heap property */
    if (leq(root, child, data)) {
      swap(root, child, data);
      root = child;
      child = 2 * root + 1;
    }
    else
      break;
  }
}

void
_h2_heapsort(uint n, bool leq(uint, uint, void *),
	     void swap(uint, uint, void *), void *data)
{
  uint      root;

  if (n < 2)
    return;

  /* Build heap */
  root = n / 2;
  while (root > 0) {
    root--;

    heap_down(root, n, leq, swap, data);
  }

  /* Sort */
  while (n > 1) {
    swap(0, n - 1, data);
    n--;
    heap_down(0, n, leq, swap, data);
  }
}

/* ------------------------------------------------------------
 Timing
 ------------------------------------------------------------ */

struct _stopwatch {
#ifdef WIN32
  DWORD     start;
  DWORD     current;
#else
#ifdef USE_OPENMP
  double      start;
  double      current;
#else
  struct tms start;
  struct tms current;
#endif
  double      clk_tck;
#endif
};

pstopwatch
new_stopwatch()
{
  pstopwatch sw;

  sw = (pstopwatch) allocmem(sizeof(stopwatch));
#ifndef WIN32
  sw->clk_tck = sysconf(_SC_CLK_TCK);
#endif

  return sw;
}

void
del_stopwatch(pstopwatch sw)
{
  freemem(sw);
}

void
start_stopwatch(pstopwatch sw)
{
#ifdef WIN32
  sw->start = timeGetTime();
#else
#ifdef USE_OPENMP
  sw->start = omp_get_wtime();
#else
  times(&sw->start);
#endif
#endif
}

real
stop_stopwatch(pstopwatch sw)
{
#ifdef WIN32
  sw->current = timeGetTime();
  return (sw->current - sw->start) * 0.001;
#else
#ifdef USE_OPENMP
  sw->current = omp_get_wtime();
  return (sw->current - sw->start);
#else
  times(&sw->current);
  return (sw->current.tms_utime + sw->current.tms_stime - sw->start.tms_utime
	  - sw->start.tms_stime) / sw->clk_tck;
#endif
#endif
}

/* ------------------------------------------------------------
 Drawing
 ------------------------------------------------------------ */

#ifdef USE_CAIRO
cairo_t  *
new_cairopdf(const char *filename, double width, double height)
{
  cairo_t  *cr;
  cairo_surface_t *pdf;

  pdf = cairo_pdf_surface_create(filename, width, height);
  cr = cairo_create(pdf);
  cairo_scale(cr, width, height);
  cairo_set_line_width(cr, 0.5 / (width > height ? width : height));
  cairo_surface_destroy(pdf);

  return cr;
}

cairo_t  *
new_cairopng(uint width, uint height)
{
  cairo_t  *cr;
  cairo_surface_t *img;

  img = cairo_image_surface_create(CAIRO_FORMAT_RGB24, (int) width,
				   (int) height);
  cr = cairo_create(img);
  cairo_surface_destroy(img);
  cairo_set_line_width(cr, 1.0);

  cairo_save(cr);
  cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
  cairo_paint(cr);
  cairo_restore(cr);

  return cr;
}

bool
write_cairopng(cairo_t * cr, const char *filename)
{
  cairo_surface_t *img;
  cairo_status_t res;

  img = cairo_get_target(cr);
  res = cairo_surface_write_to_png(img, filename);

  return (res == CAIRO_STATUS_SUCCESS);
}
#endif
