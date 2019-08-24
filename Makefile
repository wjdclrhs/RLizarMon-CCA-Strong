CC = gcc

CFLAGS=-O3 -fomit-frame-pointer -msse2avx -mavx2 -march=native -std=c99

all :
	$(CC) $(CFLAGS) -c RLizarMon_Strong.c main.c randombytes.c fips202.c bch.c ecc.c
	$(CC) $(CFLAGS) -o RLizarMon_Strong RLizarMon_Strong.o main.o randombytes.o fips202.o bch.o ecc.o
	
run : all
	./RLizarMon_Strong

clean :
	rm -f *.o
	rm -f RLizarMon_Strong

new :
	make clean
	make all
	./RLizarMon_Strong
