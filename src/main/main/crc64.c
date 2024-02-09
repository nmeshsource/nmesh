/* crc64.c - compute 64 bit CRCs */

#include <stdint.h>
#include <stddef.h>

/* The hex numbers indicating the polynomials used in CRC can be given in
   serveral ways, e.g. Normal or Reversed.
   We use CRC-64-ECMA which has:
   Normal      0x42f0e1eba9ea3693
   Reversed    0xc96c5795d7870f42

   Here we use 0xc96c5795d7870f42, which is the the bit-reversed encoding of
   0x42f0e1eba9ea3693. */
#define CRC64_POLY UINT64_C(0xc96c5795d7870f42)

/* table for crc64 with 0xc96c5795d7870f42 */
static const uint64_t crc64_table[256] = {
  UINT64_C(0x0000000000000000), UINT64_C(0xb32e4cbe03a75f6f),
  UINT64_C(0xf4843657a840a05b), UINT64_C(0x47aa7ae9abe7ff34),
  UINT64_C(0x7bd0c384ff8f5e33), UINT64_C(0xc8fe8f3afc28015c),
  UINT64_C(0x8f54f5d357cffe68), UINT64_C(0x3c7ab96d5468a107),
  UINT64_C(0xf7a18709ff1ebc66), UINT64_C(0x448fcbb7fcb9e309),
  UINT64_C(0x0325b15e575e1c3d), UINT64_C(0xb00bfde054f94352),
  UINT64_C(0x8c71448d0091e255), UINT64_C(0x3f5f08330336bd3a),
  UINT64_C(0x78f572daa8d1420e), UINT64_C(0xcbdb3e64ab761d61),
  UINT64_C(0x7d9ba13851336649), UINT64_C(0xceb5ed8652943926),
  UINT64_C(0x891f976ff973c612), UINT64_C(0x3a31dbd1fad4997d),
  UINT64_C(0x064b62bcaebc387a), UINT64_C(0xb5652e02ad1b6715),
  UINT64_C(0xf2cf54eb06fc9821), UINT64_C(0x41e11855055bc74e),
  UINT64_C(0x8a3a2631ae2dda2f), UINT64_C(0x39146a8fad8a8540),
  UINT64_C(0x7ebe1066066d7a74), UINT64_C(0xcd905cd805ca251b),
  UINT64_C(0xf1eae5b551a2841c), UINT64_C(0x42c4a90b5205db73),
  UINT64_C(0x056ed3e2f9e22447), UINT64_C(0xb6409f5cfa457b28),
  UINT64_C(0xfb374270a266cc92), UINT64_C(0x48190ecea1c193fd),
  UINT64_C(0x0fb374270a266cc9), UINT64_C(0xbc9d3899098133a6),
  UINT64_C(0x80e781f45de992a1), UINT64_C(0x33c9cd4a5e4ecdce),
  UINT64_C(0x7463b7a3f5a932fa), UINT64_C(0xc74dfb1df60e6d95),
  UINT64_C(0x0c96c5795d7870f4), UINT64_C(0xbfb889c75edf2f9b),
  UINT64_C(0xf812f32ef538d0af), UINT64_C(0x4b3cbf90f69f8fc0),
  UINT64_C(0x774606fda2f72ec7), UINT64_C(0xc4684a43a15071a8),
  UINT64_C(0x83c230aa0ab78e9c), UINT64_C(0x30ec7c140910d1f3),
  UINT64_C(0x86ace348f355aadb), UINT64_C(0x3582aff6f0f2f5b4),
  UINT64_C(0x7228d51f5b150a80), UINT64_C(0xc10699a158b255ef),
  UINT64_C(0xfd7c20cc0cdaf4e8), UINT64_C(0x4e526c720f7dab87),
  UINT64_C(0x09f8169ba49a54b3), UINT64_C(0xbad65a25a73d0bdc),
  UINT64_C(0x710d64410c4b16bd), UINT64_C(0xc22328ff0fec49d2),
  UINT64_C(0x85895216a40bb6e6), UINT64_C(0x36a71ea8a7ace989),
  UINT64_C(0x0adda7c5f3c4488e), UINT64_C(0xb9f3eb7bf06317e1),
  UINT64_C(0xfe5991925b84e8d5), UINT64_C(0x4d77dd2c5823b7ba),
  UINT64_C(0x64b62bcaebc387a1), UINT64_C(0xd7986774e864d8ce),
  UINT64_C(0x90321d9d438327fa), UINT64_C(0x231c512340247895),
  UINT64_C(0x1f66e84e144cd992), UINT64_C(0xac48a4f017eb86fd),
  UINT64_C(0xebe2de19bc0c79c9), UINT64_C(0x58cc92a7bfab26a6),
  UINT64_C(0x9317acc314dd3bc7), UINT64_C(0x2039e07d177a64a8),
  UINT64_C(0x67939a94bc9d9b9c), UINT64_C(0xd4bdd62abf3ac4f3),
  UINT64_C(0xe8c76f47eb5265f4), UINT64_C(0x5be923f9e8f53a9b),
  UINT64_C(0x1c4359104312c5af), UINT64_C(0xaf6d15ae40b59ac0),
  UINT64_C(0x192d8af2baf0e1e8), UINT64_C(0xaa03c64cb957be87),
  UINT64_C(0xeda9bca512b041b3), UINT64_C(0x5e87f01b11171edc),
  UINT64_C(0x62fd4976457fbfdb), UINT64_C(0xd1d305c846d8e0b4),
  UINT64_C(0x96797f21ed3f1f80), UINT64_C(0x2557339fee9840ef),
  UINT64_C(0xee8c0dfb45ee5d8e), UINT64_C(0x5da24145464902e1),
  UINT64_C(0x1a083bacedaefdd5), UINT64_C(0xa9267712ee09a2ba),
  UINT64_C(0x955cce7fba6103bd), UINT64_C(0x267282c1b9c65cd2),
  UINT64_C(0x61d8f8281221a3e6), UINT64_C(0xd2f6b4961186fc89),
  UINT64_C(0x9f8169ba49a54b33), UINT64_C(0x2caf25044a02145c),
  UINT64_C(0x6b055fede1e5eb68), UINT64_C(0xd82b1353e242b407),
  UINT64_C(0xe451aa3eb62a1500), UINT64_C(0x577fe680b58d4a6f),
  UINT64_C(0x10d59c691e6ab55b), UINT64_C(0xa3fbd0d71dcdea34),
  UINT64_C(0x6820eeb3b6bbf755), UINT64_C(0xdb0ea20db51ca83a),
  UINT64_C(0x9ca4d8e41efb570e), UINT64_C(0x2f8a945a1d5c0861),
  UINT64_C(0x13f02d374934a966), UINT64_C(0xa0de61894a93f609),
  UINT64_C(0xe7741b60e174093d), UINT64_C(0x545a57dee2d35652),
  UINT64_C(0xe21ac88218962d7a), UINT64_C(0x5134843c1b317215),
  UINT64_C(0x169efed5b0d68d21), UINT64_C(0xa5b0b26bb371d24e),
  UINT64_C(0x99ca0b06e7197349), UINT64_C(0x2ae447b8e4be2c26),
  UINT64_C(0x6d4e3d514f59d312), UINT64_C(0xde6071ef4cfe8c7d),
  UINT64_C(0x15bb4f8be788911c), UINT64_C(0xa6950335e42fce73),
  UINT64_C(0xe13f79dc4fc83147), UINT64_C(0x521135624c6f6e28),
  UINT64_C(0x6e6b8c0f1807cf2f), UINT64_C(0xdd45c0b11ba09040),
  UINT64_C(0x9aefba58b0476f74), UINT64_C(0x29c1f6e6b3e0301b),
  UINT64_C(0xc96c5795d7870f42), UINT64_C(0x7a421b2bd420502d),
  UINT64_C(0x3de861c27fc7af19), UINT64_C(0x8ec62d7c7c60f076),
  UINT64_C(0xb2bc941128085171), UINT64_C(0x0192d8af2baf0e1e),
  UINT64_C(0x4638a2468048f12a), UINT64_C(0xf516eef883efae45),
  UINT64_C(0x3ecdd09c2899b324), UINT64_C(0x8de39c222b3eec4b),
  UINT64_C(0xca49e6cb80d9137f), UINT64_C(0x7967aa75837e4c10),
  UINT64_C(0x451d1318d716ed17), UINT64_C(0xf6335fa6d4b1b278),
  UINT64_C(0xb199254f7f564d4c), UINT64_C(0x02b769f17cf11223),
  UINT64_C(0xb4f7f6ad86b4690b), UINT64_C(0x07d9ba1385133664),
  UINT64_C(0x4073c0fa2ef4c950), UINT64_C(0xf35d8c442d53963f),
  UINT64_C(0xcf273529793b3738), UINT64_C(0x7c0979977a9c6857),
  UINT64_C(0x3ba3037ed17b9763), UINT64_C(0x888d4fc0d2dcc80c),
  UINT64_C(0x435671a479aad56d), UINT64_C(0xf0783d1a7a0d8a02),
  UINT64_C(0xb7d247f3d1ea7536), UINT64_C(0x04fc0b4dd24d2a59),
  UINT64_C(0x3886b22086258b5e), UINT64_C(0x8ba8fe9e8582d431),
  UINT64_C(0xcc0284772e652b05), UINT64_C(0x7f2cc8c92dc2746a),
  UINT64_C(0x325b15e575e1c3d0), UINT64_C(0x8175595b76469cbf),
  UINT64_C(0xc6df23b2dda1638b), UINT64_C(0x75f16f0cde063ce4),
  UINT64_C(0x498bd6618a6e9de3), UINT64_C(0xfaa59adf89c9c28c),
  UINT64_C(0xbd0fe036222e3db8), UINT64_C(0x0e21ac88218962d7),
  UINT64_C(0xc5fa92ec8aff7fb6), UINT64_C(0x76d4de52895820d9),
  UINT64_C(0x317ea4bb22bfdfed), UINT64_C(0x8250e80521188082),
  UINT64_C(0xbe2a516875702185), UINT64_C(0x0d041dd676d77eea),
  UINT64_C(0x4aae673fdd3081de), UINT64_C(0xf9802b81de97deb1),
  UINT64_C(0x4fc0b4dd24d2a599), UINT64_C(0xfceef8632775faf6),
  UINT64_C(0xbb44828a8c9205c2), UINT64_C(0x086ace348f355aad),
  UINT64_C(0x34107759db5dfbaa), UINT64_C(0x873e3be7d8faa4c5),
  UINT64_C(0xc094410e731d5bf1), UINT64_C(0x73ba0db070ba049e),
  UINT64_C(0xb86133d4dbcc19ff), UINT64_C(0x0b4f7f6ad86b4690),
  UINT64_C(0x4ce50583738cb9a4), UINT64_C(0xffcb493d702be6cb),
  UINT64_C(0xc3b1f050244347cc), UINT64_C(0x709fbcee27e418a3),
  UINT64_C(0x3735c6078c03e797), UINT64_C(0x841b8ab98fa4b8f8),
  UINT64_C(0xadda7c5f3c4488e3), UINT64_C(0x1ef430e13fe3d78c),
  UINT64_C(0x595e4a08940428b8), UINT64_C(0xea7006b697a377d7),
  UINT64_C(0xd60abfdbc3cbd6d0), UINT64_C(0x6524f365c06c89bf),
  UINT64_C(0x228e898c6b8b768b), UINT64_C(0x91a0c532682c29e4),
  UINT64_C(0x5a7bfb56c35a3485), UINT64_C(0xe955b7e8c0fd6bea),
  UINT64_C(0xaeffcd016b1a94de), UINT64_C(0x1dd181bf68bdcbb1),
  UINT64_C(0x21ab38d23cd56ab6), UINT64_C(0x9285746c3f7235d9),
  UINT64_C(0xd52f0e859495caed), UINT64_C(0x6601423b97329582),
  UINT64_C(0xd041dd676d77eeaa), UINT64_C(0x636f91d96ed0b1c5),
  UINT64_C(0x24c5eb30c5374ef1), UINT64_C(0x97eba78ec690119e),
  UINT64_C(0xab911ee392f8b099), UINT64_C(0x18bf525d915feff6),
  UINT64_C(0x5f1528b43ab810c2), UINT64_C(0xec3b640a391f4fad),
  UINT64_C(0x27e05a6e926952cc), UINT64_C(0x94ce16d091ce0da3),
  UINT64_C(0xd3646c393a29f297), UINT64_C(0x604a2087398eadf8),
  UINT64_C(0x5c3099ea6de60cff), UINT64_C(0xef1ed5546e415390),
  UINT64_C(0xa8b4afbdc5a6aca4), UINT64_C(0x1b9ae303c601f3cb),
  UINT64_C(0x56ed3e2f9e224471), UINT64_C(0xe5c372919d851b1e),
  UINT64_C(0xa26908783662e42a), UINT64_C(0x114744c635c5bb45),
  UINT64_C(0x2d3dfdab61ad1a42), UINT64_C(0x9e13b115620a452d),
  UINT64_C(0xd9b9cbfcc9edba19), UINT64_C(0x6a978742ca4ae576),
  UINT64_C(0xa14cb926613cf817), UINT64_C(0x1262f598629ba778),
  UINT64_C(0x55c88f71c97c584c), UINT64_C(0xe6e6c3cfcadb0723),
  UINT64_C(0xda9c7aa29eb3a624), UINT64_C(0x69b2361c9d14f94b),
  UINT64_C(0x2e184cf536f3067f), UINT64_C(0x9d36004b35545910),
  UINT64_C(0x2b769f17cf112238), UINT64_C(0x9858d3a9ccb67d57),
  UINT64_C(0xdff2a94067518263), UINT64_C(0x6cdce5fe64f6dd0c),
  UINT64_C(0x50a65c93309e7c0b), UINT64_C(0xe388102d33392364),
  UINT64_C(0xa4226ac498dedc50), UINT64_C(0x170c267a9b79833f),
  UINT64_C(0xdcd7181e300f9e5e), UINT64_C(0x6ff954a033a8c131),
  UINT64_C(0x28532e49984f3e05), UINT64_C(0x9b7d62f79be8616a),
  UINT64_C(0xa707db9acf80c06d), UINT64_C(0x14299724cc279f02),
  UINT64_C(0x5383edcd67c06036), UINT64_C(0xe0ada17364673f59)
};


/* return CRC if we start with crc */
uint64_t crc64_continue(uint64_t crc, const void *buf, size_t nbytes)
{
  const uint8_t *data = (const uint8_t *) buf;
  crc = crc ^ UINT64_C(0xffffffffffffffff);

  while (nbytes--)
  {
    uint32_t idx = ((uint32_t) (crc ^ *data++)) & 0xff;
    crc = crc64_table[idx] ^ (crc >> 8);
  }

  return crc ^ UINT64_C(0xffffffffffffffff);
}

/* return CRC if we start with 0 */
uint64_t crc64_0start(const void *buf, size_t nbytes)
{
  return crc64_continue(0, buf, nbytes);
}

/* write CRC into crc and also record the total number of bytes */
void crc64_continue_counters(const void *buf, size_t nbytes,
                             uint64_t *crc, size_t *ntotalbytes)
{
  *crc = crc64_continue(*crc, buf, nbytes);
  *ntotalbytes += nbytes;
}


/* table needed in crc64_x_pow_n_modP_ */
static const uint64_t crc64_x_pow_2n[64] = {
  UINT64_C(0x4000000000000000), UINT64_C(0x2000000000000000),
  UINT64_C(0x0800000000000000), UINT64_C(0x0080000000000000),
  UINT64_C(0x0000800000000000), UINT64_C(0x0000000080000000),
  UINT64_C(0xc96c5795d7870f42), UINT64_C(0x6d5f4ad7e3c3afa0),
  UINT64_C(0xd49f7e445077d8ea), UINT64_C(0x040fb02a53c216fa),
  UINT64_C(0x6bec35957b9ef3a0), UINT64_C(0xb0e3bb0658964afe),
  UINT64_C(0x218578c7a2dff638), UINT64_C(0x6dbb920f24dd5cf2),
  UINT64_C(0x7a140cfcdb4d5eb5), UINT64_C(0x41b3705ecbc4057b),
  UINT64_C(0xd46ab656accac1ea), UINT64_C(0x329beda6fc34fb73),
  UINT64_C(0x51a4fcd4350b9797), UINT64_C(0x314fa85637efae9d),
  UINT64_C(0xacf27e9a1518d512), UINT64_C(0xffe2a3388a4d8ce7),
  UINT64_C(0x48b9697e60cc2e4e), UINT64_C(0xada73cb78dd62460),
  UINT64_C(0x3ea5454d8ce5c1bb), UINT64_C(0x5e84e3a6c70feaf1),
  UINT64_C(0x90fd49b66cbd81d1), UINT64_C(0xe2943e0c1db254e8),
  UINT64_C(0xecfa6adeca8834a1), UINT64_C(0xf513e212593ee321),
  UINT64_C(0xf36ae57331040916), UINT64_C(0x63fbd333b87b6717),
  UINT64_C(0xbd60f8e152f50b8b), UINT64_C(0xa5ce4a8299c1567d),
  UINT64_C(0x0bd445f0cbdb55ee), UINT64_C(0xfdd6824e20134285),
  UINT64_C(0xcead8b6ebda2227a), UINT64_C(0xe44b17e4f5d4fb5c),
  UINT64_C(0x9b29c81ad01ca7c5), UINT64_C(0x1b4366e40fea4055),
  UINT64_C(0x27bca1551aae167b), UINT64_C(0xaa57bcd1b39a5690),
  UINT64_C(0xd7fce83fa1234db9), UINT64_C(0xcce4986efea3ff8e),
  UINT64_C(0x3602a4d9e65341f1), UINT64_C(0x722b1da2df516145),
  UINT64_C(0xecfc3ddd3a08da83), UINT64_C(0x0fb96dcca83507e6),
  UINT64_C(0x125f2fe78d70f080), UINT64_C(0x842f50b7651aa516),
  UINT64_C(0x09bc34188cd9836f), UINT64_C(0xf43666c84196d909),
  UINT64_C(0xb56feb30c0df6ccb), UINT64_C(0xaa66e04ce7f30958),
  UINT64_C(0xb7b1187e9af29547), UINT64_C(0x113255f8476495de),
  UINT64_C(0x8fb19f783095d77e), UINT64_C(0xaec4aacc7c82b133),
  UINT64_C(0xf64e6d09218428cf), UINT64_C(0x036a72ea5ac258a0),
  UINT64_C(0x5235ef12eb7aaa6a), UINT64_C(0x2fed7b1685657853),
  UINT64_C(0x8ef8951d46606fb5), UINT64_C(0x9d58c1090f034d14)
};


/* Compute (a*b) mod P
   This is done similar to   Crc Multiply(const Crc &aa, const Crc &bb)
   from https://code.google.com/p/crcutil/source/browse/code/gf_util.h */
static inline uint64_t crc64_multiply_(uint64_t a, uint64_t b)
{
  uint64_t r = 0, h = UINT64_C(1) << 63;

  if((a ^ (a-1)) < (b ^ (b-1)))
  {
    uint64_t t = a;
    a = b;
    b = t;
  }

  if(a == 0) return 0;

  for(; a != 0; a <<= 1)
  {
    if(a & h)
    {
      r ^= b;
      a ^= h;
    }

    b = (b >> 1) ^ ((b & 1) ? CRC64_POLY : 0);
  }
  return r;
}


/* Return x**n mod P */
static inline uint64_t crc64_x_pow_n_modP_(uint64_t n)
{
  uint64_t r = UINT64_C(1) << 63;
  size_t i;

  for(i = 0; n != 0; n >>= 1, i++)
    if(n & 1) r = crc64_multiply_(r, crc64_x_pow_2n[i]);

  return r;
}


/* Combine the CRCs from the 2 buffers X1 and X2 that have both been started
   from 0.
   I.e. we return the CRC of the buffer X = concatenate(X1, X2).
   In this case we have:
   CRC(X, a) = CRC(X2, CRC(X1, a))
   CRC(X, b) = CRC(X, a) + ((b-a)x^|X|) mod P. */
uint64_t crc64_0start_combine(uint64_t crc01, uint64_t crc02, size_t nbytes2)
{
  return crc02 ^ crc64_multiply_(crc01, crc64_x_pow_n_modP_(8*nbytes2));
}

/* add crc02 of len nbytes2 to crc01 and record total num of bytes */
void crc64_0start_combine_counters(uint64_t *crc01, size_t *ntotalbytes,
                                   uint64_t  crc02, size_t  nbytes2)
{
  *crc01 = crc64_0start_combine(*crc01, crc02, nbytes2);
  *ntotalbytes += nbytes2;
}


/**************************************************************************/
/* functions that collect global CRCs via MPI */
/**************************************************************************/
#include "nmesh.h"

/* Compute CRC locally with func crc_local, and then combine the results
   from all ranks.
   Note this sets *crc=0 at the start. */
void crc64_0start_global(tMesh *mesh,
                         void (*crc_local)(void *obj,
                                           uint64_t *crc, size_t *cnt),
                         void *obj, uint64_t *crc, size_t *cnt)
{
  int rk, size = nMPI_size();
  ulong sbuf[2];
  ulong *rbuf = NULL;

  if(Rank0) rbuf = malloc(sizeof(rbuf[0])*2 * size);

  /* get my local crc and put it into sbuf */
  *crc = *cnt = 0; // START WITH 0
  crc_local(obj, crc, cnt);
  sbuf[0] = *crc;
  sbuf[1] = *cnt;

  /* rank0 now gathers all CRCs */
  nMPI_Gather(sbuf,2, nMPI_UNSIGNED_LONG, rbuf,2, nMPI_UNSIGNED_LONG, 0);

  /* rank0 now combines the CRCs */
  if(Rank0)
  {
    for(rk=1; rk<size; rk++)
      crc64_0start_combine_counters(crc, cnt, rbuf[rk*2], rbuf[rk*2+1]);
    free(rbuf);
  }
}
