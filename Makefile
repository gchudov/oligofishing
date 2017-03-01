all: oligofishing

clean:
	rm oligofishing
oligofishing: fish.c crc32.c crc32.h fasta.c fasta.h Makefile
	gcc -o oligofishing fish.c crc32.c fasta.c -std=c99
