#include <stdint.h>

typedef struct
{
    uint32_t seq_len;
    uint32_t name_len;
    uint8_t* seq;
    uint8_t* name;
    void *   callback_arg;
} fasta_item;

extern int fasta_read(const char * filename, int (*callback)(fasta_item* item), void *callback_arg);
