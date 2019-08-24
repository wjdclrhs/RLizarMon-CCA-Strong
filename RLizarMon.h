#include <stdint.h>
#include "params.h"

#define iter 100		// iteration number for keygen & EncDec test
#define testnum 10	// repeatetion number of Enc Dec procedure in a single iteration

#define sft (sizeof(size_t) * 4 - 1)

clock_t start, finish, elapsed1, elapsed2;

int Keygen(unsigned char *pk, unsigned char *sk);

int Enc(unsigned char *c, const unsigned char *m, const unsigned char *pk);

int Dec(unsigned char *m, const unsigned char *c, const unsigned char *sk, const unsigned char *pk);


