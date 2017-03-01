#include "crc32.h"
#include <stdlib.h>
#include <stdio.h>

#define GF2_DIM (32)

uint32_t crc32_table_s[256];

const uint32_t uPolynomial = 0x04c11db7;
const uint32_t uReversePolynomial = 0xedb88320;
const uint32_t uReversePolynomial2 = 0xdb710641;

static uint32_t combineTable[GF2_DIM][GF2_DIM];
static uint32_t substractTable[GF2_DIM][GF2_DIM];

uint32_t crc32_bytes(uint32_t crc, uint8_t* uint8_ts, long count)
{
    for (long i = 0; i < count; i++)
        crc = crc32_byte(crc, uint8_ts[i]);
    return crc;
}

static uint32_t Reflect(uint32_t val, int ch)
{
    uint32_t value = 0;
    // Swap bit 0 for bit 7
    // bit 1 for bit 6, etc.
    for (int i = 1; i < (ch + 1); i++)
    {
        if (0 != (val & 1))
            value |= 1U << (ch - i);
        val >>= 1;
    }
    return value;
}

static uint32_t gf2_matrix_times(uint32_t umat[GF2_DIM], uint32_t uvec)
{
    int32_t vec = (int32_t)uvec;
    int32_t* mat = (int32_t*)umat;
    return (uint32_t)(
        (*(mat++) & ((vec << 31) >> 31)) ^
        (*(mat++) & ((vec << 30) >> 31)) ^
        (*(mat++) & ((vec << 29) >> 31)) ^
        (*(mat++) & ((vec << 28) >> 31)) ^
        (*(mat++) & ((vec << 27) >> 31)) ^
        (*(mat++) & ((vec << 26) >> 31)) ^
        (*(mat++) & ((vec << 25) >> 31)) ^
        (*(mat++) & ((vec << 24) >> 31)) ^
        (*(mat++) & ((vec << 23) >> 31)) ^
        (*(mat++) & ((vec << 22) >> 31)) ^
        (*(mat++) & ((vec << 21) >> 31)) ^
        (*(mat++) & ((vec << 20) >> 31)) ^
        (*(mat++) & ((vec << 19) >> 31)) ^
        (*(mat++) & ((vec << 18) >> 31)) ^
        (*(mat++) & ((vec << 17) >> 31)) ^
        (*(mat++) & ((vec << 16) >> 31)) ^
        (*(mat++) & ((vec << 15) >> 31)) ^
        (*(mat++) & ((vec << 14) >> 31)) ^
        (*(mat++) & ((vec << 13) >> 31)) ^
        (*(mat++) & ((vec << 12) >> 31)) ^
        (*(mat++) & ((vec << 11) >> 31)) ^
        (*(mat++) & ((vec << 10) >> 31)) ^
        (*(mat++) & ((vec << 9) >> 31)) ^
        (*(mat++) & ((vec << 8) >> 31)) ^
        (*(mat++) & ((vec << 7) >> 31)) ^
        (*(mat++) & ((vec << 6) >> 31)) ^
        (*(mat++) & ((vec << 5) >> 31)) ^
        (*(mat++) & ((vec << 4) >> 31)) ^
        (*(mat++) & ((vec << 3) >> 31)) ^
        (*(mat++) & ((vec << 2) >> 31)) ^
        (*(mat++) & ((vec << 1) >> 31)) ^
        (*(mat++) & (vec >> 31)));
}

static void gf2_matrix_square(uint32_t square[GF2_DIM], uint32_t mat[GF2_DIM])
{
    for (int n = 0; n < GF2_DIM; n++)
        square[n] = gf2_matrix_times(mat, mat[n]);
}

uint32_t crc32_combine(uint32_t crc1, uint32_t crc2, long len2)
{
    /* degenerate case */
    if (len2 == 0)
        return crc1;
    if (crc1 == 0)
        return crc2;
    if (len2 < 0) 
    {
        fprintf(stderr,"crc32_combine length cannot be negative");
        abort();
    }

    int n = 3;
    do
    {
        /* apply zeros operator for this bit of len2 */
        if ((len2 & 1) != 0)
            crc1 = gf2_matrix_times(combineTable[n], crc1);
        len2 >>= 1;
        n = (n + 1) & (GF2_DIM - 1);
        /* if no more bits set, then done */
    } while (len2 != 0);

    /* return combined crc */
    crc1 ^= crc2;
    return crc1;
}

uint32_t crc32_substract(uint32_t crc1, uint32_t crc2, long len2)
{
    /* degenerate case */
    if (len2 == 0)
        return crc1;
    if (len2 < 0)
    {
        fprintf(stderr,"crc32::Substract length cannot be negative");
        abort();
    }

    crc1 ^= crc2;

    int n = 3;
    do
    {
        /* apply zeros operator for this bit of len2 */
        if ((len2 & 1) != 0)
            crc1 = gf2_matrix_times(substractTable[n], crc1);
        len2 >>= 1;
        n = (n + 1) & (GF2_DIM - 1);
        /* if no more bits set, then done */
    } while (len2 != 0);

    /* return combined crc */
    return crc1;
}

void crc32_init_window(uint32_t crc32_window_table[256], long len)
{
    for (uint32_t i = 0; i < 256; i++)
    {
        crc32_window_table[i] = crc32_combine(crc32_byte(0, i & 0xff), 0, len);
    }
}

void crc32_init(void)
{
    for (uint32_t i = 0; i < 256; i++)
    {
        crc32_table_s[i] = Reflect(i, 8) << 24;
        for (int j = 0; j < 8; j++)
            crc32_table_s[i] = (crc32_table_s[i] << 1) ^ ((crc32_table_s[i] & (1U << 31)) == 0 ? 0 : uPolynomial);
        crc32_table_s[i] = Reflect(crc32_table_s[i], 32);
    }
    combineTable[0][0] = uReversePolynomial;
    substractTable[0][31] = uReversePolynomial2;
    for (int n = 1; n < GF2_DIM; n++)
    {
        combineTable[0][n] = 1U << (n - 1);
        substractTable[0][n - 1] = 1U << n;
    }
    for (int i = 1; i < GF2_DIM; i++)
    {
        gf2_matrix_square(combineTable[i], combineTable[i - 1]);
        gf2_matrix_square(substractTable[i], substractTable[i - 1]);
    }
}
