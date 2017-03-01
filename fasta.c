#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include "fasta.h"

typedef enum
{
    fasta_state_gt,
    fasta_state_name,
    fasta_state_lf1,
    fasta_state_seq,
    fasta_state_lf2
} fasta_state_t;

int fasta_read(const char * filename, int (*callback)(fasta_item* item), void *callback_arg)
{
    int fd = open(filename, O_RDONLY);
    if (fd < 0) { perror("open failed"); return -1; }
    struct stat sb;
    if (fstat(fd, &sb) < 0) { perror("stat failed"); return -1; }
    off_t file_size = sb.st_size;
    uint8_t* fasta_data = mmap(0, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
    if (fasta_data == MAP_FAILED) { perror("mmap failed"); return -1; }

    fasta_item item;
    item.callback_arg = callback_arg;
    item.seq = NULL;
    item.seq_len = 0;
    item.name = NULL;
    item.name_len = 0;

    fasta_state_t state = fasta_state_gt;
    for (off_t i = 0; i < file_size; i++)
    {
        uint8_t c = fasta_data[i];
        switch (state)
        {
            case fasta_state_gt:
                if (c == '>') { 
                    state = fasta_state_name;
                    item.name = fasta_data + i + 1;
                    break;
                }
                fprintf(stderr, "fasta format invalid, expecting >");
                return -2;
            case fasta_state_name:
                if (c == '\r') { 
                    state = fasta_state_lf1;
                    item.name_len = fasta_data + i - item.name;
                    break;
                }
                if (c == '\n') { 
                    state = fasta_state_seq;
                    item.seq = fasta_data + i + 1;
                    break;
                }
                break;
            case fasta_state_lf1:
                if (c == '\n') {
                    state = fasta_state_seq;
                    item.seq = fasta_data + i + 1;
                    break;
                }
                fprintf(stderr, "fasta format invalid, expecting \\n");
                return -2;
            case fasta_state_seq:
                if (c == '\r') { 
                    state = fasta_state_lf2;
                    item.seq_len = fasta_data + i - item.seq;
                    break;
                }
                if (c == '\n') { 
                    state = fasta_state_gt;
                    callback(&item);
                    break;
                }
                break;
            case fasta_state_lf2:
                if (c == '\n') {
                    state = fasta_state_gt;
                    callback(&item);
                    break;
                }
                fprintf(stderr, "fasta format invalid, expecting \\n");
                return -2;
        }
    }
}
