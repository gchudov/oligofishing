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
static uint32_t min_window_len = 0xffffffff;
static uint32_t max_window_len = 0;
typedef uint32_t crc32_window_table_t[256];
static crc32_window_table_t* crc32_window_table_s;

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

static inline uint32_t crc32_sliding(uint32_t crc, uint32_t wnd, uint8_t head, uint8_t tail)
{
    return crc32_byte(crc, head) ^ crc32_window_table_s[wnd][tail];
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
        if (c == '>')
        {
            item->seq_len = pos;
            break;
        }
        if (c > 32)
        {
            crc = crc32_byte(crc, c);
            len ++;
        }
    }
    if (min_window_len > len) min_window_len = len;
    if (max_window_len < len) max_window_len = len;
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
    uint32_t buf_mask = 0;
    while (buf_mask <= max_window_len + 1)
        buf_mask = 1 + (buf_mask << 1);
    uint8_t buf[buf_mask + 1];
    uint32_t crc[max_window_len + 1];
    memset(crc, 0, (max_window_len + 1) * sizeof(crc[0]));
    // off_next - position of the next byte of sequence in circular buffer buf;
    uint32_t off_next = 0;
    uint8_t * next = item->seq;
    uint8_t * last = item->seq + item->seq_len;
    buf[buf_mask] = 0;
    while (next < last)
    {
        uint8_t c = *(next++);
        if (c <= 32) continue;
        if (c == '>')
        {
            last = next - 1;
            item->seq_len = last - item->seq;
            break;
        }
        buf[(off_next++) & buf_mask] = c;
        for (int wnd = min_window_len; wnd <= max_window_len; wnd++)
        {
            if (off_next < wnd)
            {
                crc[wnd] = crc32_byte(crc[wnd], c);
                continue;
            }
            uint8_t tail = buf[(off_next - 1 - wnd) & buf_mask];
            uint32_t crc_wnd = crc[wnd] = crc32_sliding(crc[wnd], wnd, c, tail);
            for (hook* p = hook_hash[crc_wnd & hook_hash_mask]; p; p=p->next)
            {
                if (p->crc != crc_wnd || p->seq_len != wnd) continue;
                int part_len = (off_next & buf_mask) > ((off_next - wnd) & buf_mask) ?
                    wnd : buf_mask + 1 - ((off_next - wnd) & buf_mask);
                if (memcmp(p->seq, buf + ((off_next - wnd) & buf_mask), part_len) || memcmp(p->seq + part_len, buf, wnd - part_len))
                    break;
                while (next < last)
                {
                    uint8_t c = *(next++);
                    if (c == '>')
                    {
                        last = next - 1;
                        item->seq_len = last - item->seq;
                        break;
                    }
                    if (c > 32) off_next++;
                }
#if 0
                fprintf(stderr, "Match: hook %.*s (%.*s), fish %.*s (%.*s)\n", 
                    p->name_len, p->name, p->seq_len, p->seq, 
                    item->name_len, item->name, item->seq_len, item->seq);
#endif
                printf(">%.*s (matched %.*s)\n%.*s", item->name_len, item->name, p->name_len, p->name, item->seq_len, item->seq);
                return 0;
            }
        }
    };
    return 0;
}

int main(int argc, char**argv)
{
    int c;
    while ((c = getopt (argc, argv, "h:p:")) != -1)
        switch (c)
        {
            case 'h':
                hooks_filename_s = optarg;
                break;
            case 'p':
                pond_filename_s = optarg;
                break;
            case '?':
                fprintf (stderr,
                    "Usage: %s -h \"%s\" -p \"%s\"\n", argv[0], hooks_filename_s, pond_filename_s);
                return 1;
              default:
                abort ();
        }

    crc32_init();
    int rc = fasta_read(hooks_filename_s, parse_hook, NULL);
    if (rc < 0) return rc;

    crc32_window_table_s = malloc(sizeof(crc32_window_table_t) * (max_window_len + 1));
    for (int wnd = min_window_len; wnd <= max_window_len; wnd++)
        crc32_init_window(crc32_window_table_s[wnd], wnd);

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
    {
        fprintf(stderr, "Built a hash table of size %d containing %d hooks of length %d..%d.\n", hook_hash_mask+1, hook_count, min_window_len, max_window_len);
        return 0;
    }

    rc = fasta_read(pond_filename_s, parse_pond, NULL);
    if (rc < 0) return rc;

    return 0;
}
