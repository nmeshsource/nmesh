/* endianIO.c */
/* Wolfgang Tichy, August 2015 */

/* FIXME: According to Rob Pike the code below is stupid.
   From https://commandcenter.blogspot.com/2012/04/byte-order-fallacy.html :
   Let's say your data stream has a little-endian-encoded 32-bit integer.
   Here's how to extract it (assuming unsigned bytes):
   i = (data[0]<<0) | (data[1]<<8) | (data[2]<<16) | (data[3]<<24);
   If it's big-endian, here's how to extract it:
   i = (data[3]<<0) | (data[2]<<8) | (data[1]<<16) | (data[0]<<24);
*/


#include "nmesh.h"

/* About binary formats:
   There is /usr/include/endian.h . Is this standard?
   also see /usr/include/byteswap.h, and
   include/bits/byteswap.h for gcc/x86 optimized code */
#include <endian.h>
#if  __BYTE_ORDER == __LITTLE_ENDIAN
#define BYTE_ORDER_LITTLE 1
#else
#define BYTE_ORDER_LITTLE 0
#endif


/* use fwrite to write an array of doubles to a file in little endian format */
size_t fwrite_double_little(const double *buf, size_t nmemb, FILE *fp)
{
  int little = BYTE_ORDER_LITTLE; /* endianess */

  /* just use fwrite if we are little endian */
  if(little)
    return fwrite(buf, sizeof(double), nmemb, fp);

  /* sanity check */
  if(sizeof(char) != 1)
    errorexit("fwrite_double_little: size of char is not 1");

  /* if we get here we assume big endian */
  if(sizeof(double) == 8)
  {
    double xdouble;
    char c[8], *x;
    size_t i, count=0;

    for(i=0; i<nmemb; i++)
    {
      xdouble = buf[i];
      x = (char *) &xdouble;
      c[0] = x[7];
      c[1] = x[6];
      c[2] = x[5];
      c[3] = x[4];
      c[4] = x[3];
      c[5] = x[2];
      c[6] = x[1];
      c[7] = x[0];
      count = count + fwrite(c, sizeof(char), 8, fp);
    }
    return count/8;
  }
  errorexit("fwrite_double_little: size of double is not 8");
  return -1; /* hopefully we never get here */
}


/* use fread to read an array of doubles from a file in little endian format */
size_t fread_double_little(double *buf, size_t nmemb, FILE *fp)
{
  int little = BYTE_ORDER_LITTLE; /* endianess */

  /* just use fread if we are little endian */
  if(little)
    return fread(buf, sizeof(double), nmemb, fp);

  /* sanity check */
  if(sizeof(char) != 1)
    errorexit("fread_double_little: size of char is not 1");

  /* if we get here we assume big endian */
  if(sizeof(double) == 8)
  {
    double xdouble;
    char c[8], *x;
    size_t i, count=0;

    for(i=0; i<nmemb; i++)
    {
      count = count + fread(c, sizeof(char), 8, fp);
      x = (char *) &xdouble;
      x[0] = c[7];
      x[1] = c[6];
      x[2] = c[5];
      x[3] = c[4];
      x[4] = c[3];
      x[5] = c[2];
      x[6] = c[1];
      x[7] = c[0];
      buf[i] = xdouble;
    }
    return count/8;
  }
  errorexit("fread_double_little: size of double is not 8");
  return -1; /* hopefully we never get here */
}


/* use fwrite to write an array, but swap byte order */
size_t fwrite_swapbytes(const void *ptr, size_t size, size_t nmemb, FILE *fp)
{
  const char *buf = ptr;
  size_t i, count;
  char c[size];
  int b;

  for(count=0, i=0; i<nmemb; i++)
  {
    for(b=0; b<size; b++) c[b] = buf[i*size + size-1-b];
    count = count + fwrite(c, sizeof(char), size, fp);
  }
  return count/size;
}

/* use fread to write an array, but swap byte order */
size_t fread_swapbytes(void *ptr, size_t size, size_t nmemb, FILE *fp)
{
  char *buf = ptr;
  size_t i, count;
  char c[size];
  int b;

  for(count=0, i=0; i<nmemb; i++)
  {
    count += fread(c, sizeof(char), size, fp);
    for(b=0; b<size; b++) buf[i*size + b] = c[size-1-b];
  }
  return count/size;
}

/* use fwrite to write to a file in little endian format */
size_t fwrite_little(const void *ptr, size_t size, size_t nmemb, FILE *fp)
{
  int little = BYTE_ORDER_LITTLE; /* endianess */

  /* just use fwrite if we are little endian */
  if(little)
    return fwrite(ptr, size, nmemb, fp);
  else
    return fwrite_swapbytes(ptr, size, nmemb, fp);
}


/* use fread to read from a file in little endian format */
size_t fread_little(double *ptr, size_t size, size_t nmemb, FILE *fp)
{
  int little = BYTE_ORDER_LITTLE; /* endianess */

  /* just use fread if we are little endian */
  if(little)
    return fread(ptr, size, nmemb, fp);
  else
    return fread_swapbytes(ptr, size, nmemb, fp);
}

/* use fwrite to write to a file in big endian format */
size_t fwrite_big(const void *ptr, size_t size, size_t nmemb, FILE *fp)
{
  int little = BYTE_ORDER_LITTLE; /* endianess */

  /* just use fwrite if we are big endian */
  if(!little)
    return fwrite(ptr, size, nmemb, fp);
  else
    return fwrite_swapbytes(ptr, size, nmemb, fp);
}


/* use fread to read from a file in big endian format */
size_t fread_big(double *ptr, size_t size, size_t nmemb, FILE *fp)
{
  int little = BYTE_ORDER_LITTLE; /* endianess */

  /* just use fread if we are big endian */
  if(!little)
    return fread(ptr, size, nmemb, fp);
  else
    return fread_swapbytes(ptr, size, nmemb, fp);
}
