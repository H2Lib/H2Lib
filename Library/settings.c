/* ------------------------------------------------------------
 This is the file "settings.h" of the H2Lib package.
 All rights reserved, Steffen Boerm 2009
 ------------------------------------------------------------ */

#include "settings.h"

#ifdef USE_FLOAT
const real r_zero = 0.0f;
const real r_one = 1.0f;
const real r_minusone = -1.0f;
#ifdef USE_COMPLEX
const field f_zero = 0.0f;
const field f_one = 1.0f;
const field f_minusone = -1.0f;
const field f_i = 0.0f + 1.0f * I;
#else
const field f_zero = 0.0f;
const field f_one = 1.0f;
const field f_minusone = -1.0f;
#endif
#else
const real r_zero = 0.0;
const real r_one = 1.0;
const real r_minusone = -1.0;
#ifdef USE_COMPLEX
const field f_zero = 0.0;
const field f_one = 1.0;
const field f_minusone = -1.0;
const field f_i = 0.0 + 1.0 * I;
#else
const field f_zero = 0.0;
const field f_one = 1.0;
const field f_minusone = -1.0;
#endif
#endif

const int i_zero = 0;
const int i_one = 1;

const uint u_zero = 0;
const uint u_one = 1;

const char *_h2_ntrans = "Not transposed";
const char *_h2_trans = "Transposed";
#ifdef USE_COMPLEX
const char *_h2_adj = "Conjugate transposed";
#else
const char *_h2_adj = "Transposed";
#endif

const char *_h2_left = "Left";
const char *_h2_right = "Right";
const char *_h2_lower = "Lower";
const char *_h2_upper = "Upper";
const char *_h2_unit = "Unit";
const char *_h2_nonunit = "Non-unit";
const char *_h2_vectors = "Vectors";
const char *_h2_skinnyvectors = "Skinny Vectors";
const char *_h2_novectors = "No vectors";
