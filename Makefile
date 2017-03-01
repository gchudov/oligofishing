all: fish

clean:
	rm fish	
fish: fish.c crc32.c crc32.h fasta.c fasta.h Makefile
	gcc -o fish fish.c crc32.c fasta.c -std=c99
