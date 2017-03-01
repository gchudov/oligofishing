#define _XOPEN_SOURCE
#include "crc32.h"
#include "fasta.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <ctype.h>

static const char *hooks_filename_s = "hooks.fasta";
static const char *pond_filename_s = "pond.fasta";
static long window_len = 26;
static uint32_t crc32_window_table_s[256];

typedef struct hook_t hook;
struct hook_t
{
    uint32_t crc;
    uint32_t seq_len;
    uint32_t name_len;
    uint8_t* seq;
    uint8_t* name;
    hook* next;
};
static hook* hook_table = NULL;
static int hook_alloc = 0;
static int hook_count = 0;
static hook** hook_hash = NULL;
static int hook_hash_mask = 0;

static inline uint32_t crc32_sliding(uint32_t crc, uint8_t head, uint8_t tail)
{
    return crc32_byte(crc, head) ^ crc32_window_table_s[tail];
}

int parse_hook(fasta_item* item)
{
//    printf("%.*s: %.*s\n", item->name_len, item->name, item->seq_len, item->seq);
    if (!hook_table)
    {
        hook_count = 0;
        hook_alloc = 1;
        hook_table = malloc(hook_alloc * sizeof(hook));
    }
    if (hook_count == hook_alloc)
    {
        hook_alloc <<= 1;
        hook_table = realloc(hook_table, hook_alloc * sizeof(hook));
    }
    uint32_t crc = 0;
    uint32_t len = 0;
    for (int pos = 0; pos < item->seq_len; pos++)
    {
        uint8_t c = item->seq[pos];
        if (c > 32)
        {
            crc = crc32_byte(crc, c);
            len ++;
        }
    }
    if (len != window_len) return 0;
    hook_table[hook_count].crc = crc;
    hook_table[hook_count].seq_len = len;
//    printf("CRC=%x vs %x\n", crc, crc32_bytes(0, item->seq, len));
    hook_table[hook_count].seq = item->seq;
    hook_table[hook_count].name_len = item->name_len;
    hook_table[hook_count].name = item->name;
    hook_count++;
    return 0;
}

int parse_pond(fasta_item* item)
{
//    printf("%.*s: %.*s\n", item->name_len, item->name, item->seq_len, item->seq);
    uint32_t crc = 0;
    uint32_t len = 0;
    uint8_t buf[window_len];
    uint8_t * next = item->seq;
    uint8_t * last = item->seq + item->seq_len;
    while (next < last && len < window_len)
    {
        uint8_t c = *(next++);
        if (c > 32)
        {
            crc = crc32_byte(crc, c);
            buf[len] = c;
            len ++;
        }
    }
    if (len < window_len) return 0;
    // off - position of the first byte in circular buffer buf;
    int off = 0;
    while (1)
    {
        for (hook* p = hook_hash[crc & hook_hash_mask]; p; p=p->next)
            if (p->crc == crc && !memcmp(p->seq, buf + off, window_len - off) && !memcmp(p->seq + window_len - off, buf, off))
            {
                // fprintf(stderr, "Match: hook %.*s (%.*s), fish %.*s (%.*s)\n", p->name_len, p->name, p->seq_len, p->seq, item->name_len, item->name, item->seq_len, item->seq);
                printf(">%.*s\n%.*s", item->name_len, item->name, item->seq_len, item->seq);
                return 0;
            }
        if (next >= last) break;
        uint8_t c = *(next++);
        if (c <= 32) continue;

        len++;

        crc = crc32_sliding(crc, c, buf[off]);
        buf[off] = c;
        off++;
//        if (crc32_bytes(crc32_bytes(0, buf + off, window_len - off), buf, off) != crc) abort();
        if (off == window_len) off = 0;
    };
    return 0;
}

int main(int argc, char**argv)
{
    int c;
    while ((c = getopt (argc, argv, "h:p:l:")) != -1)
        switch (c)
        {
            case 'h':
                hooks_filename_s = optarg;
                break;
            case 'p':
                pond_filename_s = optarg;
                break;
            case 'l':
                window_len = atoi(optarg);
                break;
            case '?':
                fprintf (stderr,
                    "Usage: %s -h \"%s\" -p \"%s\" -l %ld\n", argv[0], hooks_filename_s, pond_filename_s, window_len);
                return 1;
              default:
                abort ();
        }

    crc32_init();
    crc32_init_window(crc32_window_table_s, window_len);
    int rc = fasta_read(hooks_filename_s, parse_hook, NULL);
    if (rc < 0) return rc;

    while (hook_count * 4 > hook_hash_mask)
        hook_hash_mask = (hook_hash_mask << 1) | 1;
    hook_hash = malloc(sizeof(hook*) * (hook_hash_mask + 1));
    memset(hook_hash, 0, sizeof(hook*) * (hook_hash_mask + 1));
    for (int i = 0; i < hook_count; i++)
    {
        uint32_t crc = hook_table[i].crc & hook_hash_mask;
        hook_table[i].next = hook_hash[crc];
        hook_hash[crc] = &hook_table[i];
    }
    if (!hook_count)
        fprintf(stderr, "Built a hash table of size %d containing %d hooks of length %ld.\n", hook_hash_mask+1, hook_count, window_len);

    rc = fasta_read(pond_filename_s, parse_pond, NULL);
    if (rc < 0) return rc;

    return 0;
}
