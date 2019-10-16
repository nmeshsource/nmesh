/* nmesh.h */
/* Wolfgang Tichy, 1/2019 */

/* _XOPEN_SOURCE>=500 needed for strdup */
#define _XOPEN_SOURCE 600
/* to use _XOPEN_SOURCE we need to include nmesh.h before all other headers
   in files such as main/main/utilities.c or main/main/parameters.c */

#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "nmesh_automatic_include.h"
