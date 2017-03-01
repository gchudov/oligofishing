all: oligofishing

clean:
	rm oligofishing
oligofishing: main.c crc32.c crc32.h fasta.c fasta.h Makefile
	gcc -o oligofishing main.c crc32.c fasta.c -std=c99
