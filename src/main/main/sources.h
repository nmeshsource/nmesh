/* sources.h */
/* Wolfgang Tichy, 12/2020 */

/***********************************************************************/
/* Here we specify which of the standard C and Posix sources we want.
   We can just include this file, wherever the compiler finds something
   missing. */
/***********************************************************************/

/* _XOPEN_SOURCE>=500 needed for strdup
   _XOPEN_SOURCE>=600 needed for snprintf */
#define _XOPEN_SOURCE 600
/* to use _XOPEN_SOURCE we need to include nmesh.h before all other headers
   in files such as main/main/utilities.c or main/main/parameters.c */
/* The value 600 (corresponding to the sixth revision) includes definitions
   from SUSv3. */
