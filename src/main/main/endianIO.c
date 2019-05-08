/* endianIO.c */
/* Wolfgang Tichy, Feb. 2019 */

/* FIXME: According to Rob Pike the code below is stupid.
   From https://commandcenter.blogspot.com/2012/04/byte-order-fallacy.html :
   Let's say your data stream has a little-endian-encoded 32-bit integer.
   Here's how to extract it (assuming unsigned bytes):
   i = (data[0]<<0) | (data[1]<<8) | (data[2]<<16) | (data[3]<<24);
   If it's big-endian, here's how to extract it:
   i = (data[3]<<0) | (data[2]<<8) | (data[1]<<16) | (data[0]<<24);
   On the other hand what I have is general and works! */


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
size_t fread_little(void *ptr, size_t size, size_t nmemb, FILE *fp)
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
size_t fread_big(void *ptr, size_t size, size_t nmemb, FILE *fp)
{
  int little = BYTE_ORDER_LITTLE; /* endianess */

  /* just use fread if we are big endian */
  if(!little)
    return fread(ptr, size, nmemb, fp);
  else
    return fread_swapbytes(ptr, size, nmemb, fp);
}
