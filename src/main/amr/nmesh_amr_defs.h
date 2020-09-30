/* nmesh_amr_defs.h */
/* Wolfgang Tichy, 7/2019 */


/* Arrays on the stack instead of the heap.
   instead of:
    tArray *array = alloc_array_with_segs(n, Ne, ns);
   we can use:
    DECL_STACK_ARRAY_with_segs(array, n, Ne, ns);
 */
#define DECL_STACK_ARRAY_with_segs(array, n_, Ne_, ns_) \
  tArray array[1]; \
  array->N = (n_[0]) * (n_[1]) * (n_[2]); \
  for(int i=0; i<3; i++)  array->n[i] = n_[i]; \
  array->Ne = (Ne_); \
  array->ns = (ns_); \
  array->d_nofree = 1; \
  array->si = array->info = 0; \
  array->par = NULL; \
  double array##data_[(array->N + (Ne_)) * (ns_)]; \
  array->size = sizeof(array##data_); \
  array->d = &(array##data_[0])

/* instead of:
    tArray *array = alloc_array_with_segs(n);
   we can use:
    DECL_STACK_ARRAY(array, n);
 */
#define DECL_STACK_ARRAY(array, n_) DECL_STACK_ARRAY_with_segs(array, n_,0,1)


/* define the CONST we use for read-only multidimensional arrays */
#define CONST
/* with gcc we can use:
#define CONST const
   Because gcc has the proposed extenstion N 1923 mentioned in
   http://www.open-std.org/jtc1/sc22/wg14/www/docs/sd3.htm */
