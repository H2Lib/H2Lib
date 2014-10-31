
/* ------------------------------------------------------------
   This is the file "parameters.h" of the H2Lib package.
   All rights reserved, Steffen Boerm 2009
   ------------------------------------------------------------ */

/** @file parameters.h
 *  @author Steffen B&ouml;rm
 */

#include "settings.h"

#ifndef PARAMETERS_H
#define PARAMETERS_H

/** @defgroup parameters parameters
 *  @brief Ask for parameters.
 *
 *  This module provides functions that accept parameters either
 *  interactively from the user or from environment variables.
 *  The idea is to avoid asking too many questions by setting
 *  a couple of environment variables and only asking for those
 *  parameters that change during a series of experiments.
 *  @{ */

/** @brief Ask for an integer.
 *
 *  @param question Question to be asked.
 *  @param envname If an environment variable of this name exists, the question is not asked and the contents of the variable are used.
 *  @param deflt Default value, will be used if the user answers the question by pressing Enter.
 *  @returns Integer taken from environment, user, or default. */
HEADER_PREFIX int
askforint(const char *question, const char *envname, int deflt);

/** @brief Ask for a character.
 *
 *  @param question Question to be asked.
 *  @param envname If an environment variable of this name exists, the question is not asked and the contents of the variable are used.
 *  @param allowed String containing characters that will be accepted.
 *  @param deflt Default value, will be used if the user answers the question by pressing Enter.
 *  @returns Character taken from environment, user, or default. */
HEADER_PREFIX char
askforchar(const char *question, const char *envname, const char *allowed, char deflt);

/** @brief Ask for a real number.
 *
 *  @param question Question to be asked.
 *  @param envname If an environment variable of this name exists, the question is not asked and the contents of the variable are used.
 *  @param deflt Default value, will be used if the user answers the question by pressing Enter.
 *  @returns Real number taken from environment, user, or default. */
HEADER_PREFIX real
askforreal(const char *question, const char *envname, real deflt);

/** @brief Ask for a string.
 *
 *  @param question Question to be asked.
 *  @param envname If an environment variable of this name exists, the question is not asked and the contents of the variable are used.
 *  @param deflt Default value, will be used if the user answers the question by pressing Enter.
 *  @param buffer Input buffer for the string.
 *  @param bufsize Size of the input buffer.
 *  @returns Pointer to `buffer`. */
HEADER_PREFIX char *
askforstring(const char *question, const char *envname,
	     const char *deflt,
	     char *buffer, uint bufsize);

/** @} */

#endif
