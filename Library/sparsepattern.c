
/* ------------------------------------------------------------
   This is the file "sparsepattern.h" of the H2Lib package.
   All rights reserved, Steffen Boerm 2012
   ------------------------------------------------------------ */

#include "sparsepattern.h"

#include "basic.h"

#include <assert.h>
#include <stdio.h>

/* Constructor for non-zero entries */

static    ppatentry
create_entry(uint row, uint col, ppatentry next)
{
  ppatentry e;

  e = (ppatentry) allocmem(sizeof(patentry));

  e->row = row;
  e->col = col;
  e->next = next;

  return e;
}

/* Destructor for non-zero entries */

static void
free_entry(ppatentry e)
{
  ppatentry next;

  while (e) {
    next = e->next;
    freemem(e);
    e = next;
  }
}

/* ------------------------------------------------------------
   Constructor
   rows: Number of rows
   cols: Number of columns
   ------------------------------------------------------------ */

psparsepattern
new_sparsepattern(uint rows, uint cols)
{
  psparsepattern sp;
  uint      i;

  assert(rows > 0);
  assert(cols > 0);

  sp = (psparsepattern) allocmem(sizeof(sparsepattern));
  sp->row = (ppatentry *) allocmem(sizeof(ppatentry) * rows);

  for (i = 0; i < rows; i++)
    sp->row[i] = NULL;

  sp->rows = rows;
  sp->cols = cols;

  return sp;
}

/* ------------------------------------------------------------
   Destructor
   ------------------------------------------------------------ */

void
del_sparsepattern(psparsepattern sp)
{
  assert(sp != NULL);
  assert(sp->row != NULL);

  clear_sparsepattern(sp);
  freemem(sp->row);
  freemem(sp);
}

/* ------------------------------------------------------------
   Set a matrix to zero by removing all non-zero entries
   ------------------------------------------------------------ */

void
clear_sparsepattern(psparsepattern sp)
{
  ppatentry *row;
  uint      rows;
  uint      i;

  assert(sp != NULL);
  assert(sp->row != NULL);

  row = sp->row;
  rows = sp->rows;

  for (i = 0; i < rows; i++)
    if (row[i] != NULL) {
      free_entry(row[i]);
      row[i] = NULL;
    }
}

/* ------------------------------------------------------------
   Add a non-zero entry
   row: Row of the entry
   col: Column of the entry
   ------------------------------------------------------------ */

void
addnz_sparsepattern(psparsepattern sp, uint row, uint col)
{
  ppatentry e;

  assert(sp != NULL);
  assert(row < sp->rows);
  assert(col < sp->cols);

  e = sp->row[row];
  while (e != NULL && e->col != col)
    e = e->next;

  if (e == NULL)
    sp->row[row] = create_entry(row, col, sp->row[row]);
}

/* ------------------------------------------------------------
   Print the contents of this object
   ------------------------------------------------------------ */

void
print_sparsepattern(pcsparsepattern sp)
{
  ppatentry e;
  uint      i;

  (void) printf("sparsepattern, rows=%d, cols=%d\n", sp->rows, sp->cols);

  for (i = 0; i < sp->rows; i++) {
    (void) printf("  %u:", i);

    e = sp->row[i];
    while (e != NULL) {
      (void) printf(" %u", e->col);

      e = e->next;
    }
    (void) printf("\n");
  }
}
