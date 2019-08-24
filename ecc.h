//#include "param.h"

//bch(511,264,59)
#define DATA_LEN 33
#define ECCBUF_LEN 64
#define ECC_LEN 31
#define MAX_ERROR 29
#define CODE_LEN 64
#define LOG_CODE_LEN 9
#define ECC_BITS 243



//error correction encode
int ecc_enc(const unsigned char *d, unsigned char *c);

//error corrction decode
int ecc_dec(unsigned char *d, unsigned char *c);
