#include <stdint.h>

extern uint32_t crc32_table_s[256];

static inline uint32_t crc32_byte(uint32_t crc, uint8_t val)
{
    return (crc >> 8) ^ crc32_table_s[(crc & 0xff) ^ val];
}

extern uint32_t crc32_bytes(uint32_t crc, uint8_t* uint8_ts, long count);
extern uint32_t crc32_combine(uint32_t crc1, uint32_t crc2, long len2);
extern uint32_t crc32_substract(uint32_t crc1, uint32_t crc2, long len2);
extern void crc32_init_window(uint32_t crc32_window_table[256], long len);
extern void crc32_init(void);
