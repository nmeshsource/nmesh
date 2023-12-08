#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#define NPBYTES 21
#define NLOCS   ((8*NPBYTES)/3)


/* unpack ploc into loc */
void connections_loc_from_ploc(char loc[NLOCS],
                               const unsigned char ploc[NPBYTES])
{
  int i, p1,p2, b1,b2, bi1,bi2;
  unsigned char ch, c1,c2;
  //p2=0;

  /* translate ploc to loc */
  for(i=0; i<NLOCS; i++)
  {
    int i3 = i*3;
    b1 = i3;        /* bit number of 1st bit in ploc */
    b2 = i3+2;      /* bit number of 3rd bit in ploc */

    p1 = b1/8;      /* position of 1st byte in ploc */

    bi1 = b1 % 8;   /* starting bit index we need in c1 */
    bi2 = b2 % 8;   /* starting bit index we need in c2 */

    /* if all is in c1 */
    c1 = ploc[p1];
    ch = c1>>bi1;   /* shift 1st bit into correct position */
    ch = ch & 7;    /* keep only lowest 3 bits */
    //printf("ch=%o\n",ch);
    if(bi2<2) /* not all is in c1 */
    {
      p2 = b2/8;       /* position of 2nd byte in ploc */
      c2 = ploc[p2];
      //printf("c2=%o\n",c2);
      c2 = c2<<(2-bi2); /* shift 3rd bit into correct position */
      //printf("c2=%o\n",c2);
      c2 = c2 & 7;     /* keep only lowest 3 bits */
      //printf("c2=%o\n",c2);
      ch = ch | c2;
    }
    loc[i] = ch + '0'; /* set loc */

    //PRF;printf(": i=%d  b1=%d b2=%d bi1=%d bi2=%d p1=%d p2=%d %c(%d)<-%o",
    //i, b1,b2, bi1,bi2, p1,p2, loc[i],loc[i]-'0', ploc[p1]);
    //if(bi2<2) printf(",%o", ploc[p2]);
    //printf("\n");
  }
  /* do not add string end marker in loc[eloc->l], this would kill
     stuff below eloc->l which we intend to keep, because sometimes we
     just decrease l in eloc */
}

/* pack loc into ploc */
void connections_loc_to_ploc(const char loc[NLOCS],
                             unsigned char ploc[NPBYTES])
{
  int i, p1,p2, b1,b2, bi1,bi2;
  unsigned char ch, c1,c2, pc;
  //p2=0;

  ploc[0] = 0;         /* init first char in ploc */
  //ploc[NPBYTES-1] = 0; /* init last char in ploc */

  /* translate ploc to loc */
  //ploc[0] = 0; // not needed
  for(i=0; i<NLOCS; i++)
  {
    int i3 = i*3;
    b1 = i3;        /* bit number of 1st bit in ploc */
    b2 = i3+2;      /* bit number of 3rd bit in ploc */

    p1 = b1/8;      /* position of 1st byte in ploc */

    bi1 = b1 % 8;   /* starting bit index we need in c1 */
    bi2 = b2 % 8;   /* starting bit index we need in c2 */

    ch = loc[i];
    //if(ch == 0) break; //Note: do not brak like this!!!
    //ch = ch - '0';
    ch = ch & 7; // this also subtracts '0'
    c1 = ch<<bi1;   /* shift 1st bit into correct position */
    pc = ploc[p1];
    pc = pc<<(8-bi1);   /* clear all left of bi1 */
    pc = pc>>(8-bi1);   /* still needed if ploc[0]=0 */
    ploc[p1] = pc | c1; /* then add the bits form c1 */

    if(bi2<2) /* not all is in one ploc */
    {
      c2 = ch>>(2-bi2); /* shift 3rd bit into correct position */
      p2 = b2/8;        /* position of 2nd byte in ploc */
      ploc[p2] = c2;    /* this also zeros all above the bit bi2 */
    }
    //PRF;printf(": i=%d  b1=%d b2=%d bi1=%d bi2=%d p1=%d p2=%d %c(%d)->%o",
    //i, b1,b2, bi1,bi2, p1,p2, loc[i],loc[i]-'0', ploc[p1]);
    //if(bi2<2) printf(",%o", ploc[p2]);
    //printf("\n");
  }
}


/* main: get loc from l and ploc */
int main(int argc, char *argv[])
{
  char str[1000];
  const char delim[] = " ,";
  char *tok;
  int l, i;
  unsigned char ploc[NPBYTES];
  char loc[NLOCS+1];

  /* zero ploc and loc */
  for(i=0; i<NPBYTES; i++) ploc[i] = 0;
  for(i=0; i<NLOCS+1; i++) loc[i] = 0;

  printf("enter l as integer and ploc as list of integers:\n");

  printf("l = ");
  fgets(str, 999, stdin);
  l = atoi(str);

  printf("ploc = ");
  fgets(str, 999, stdin);
  for(i=0, tok = strtok(str, delim); tok!=NULL; i++)
  {
    //printf("tok=%s\n", tok);
    ploc[i] = atoi(tok);
    tok = strtok(NULL, delim);
  }
  //for(i=0; i<l; i++) printf("ploc[%d] = %d\n", i, (int) ploc[i]);

  connections_loc_from_ploc(loc, ploc);
  printf("l   = %d\n", l);
  printf("loc = ");
  for(i=0; i<l; i++) printf("%c", loc[i]);
  printf("\n");
}
