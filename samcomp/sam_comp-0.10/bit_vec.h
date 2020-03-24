#include <stdlib.h>

struct bit_vec {
    bit_vec(int size);
    ~bit_vec();

    int  get(int idx);
    void set(int idx);
    void clear(int idx);
    void clear_range(int upto);

protected:
    uint32_t *bits;
};

bit_vec::bit_vec(int size) {
    bits = (uint32_t *)calloc(size/32+1, sizeof(*bits));
}

bit_vec::~bit_vec() {
    free(bits);
}

int bit_vec::get(int idx) {
    return bits[idx/32] & (1 << (idx%32));
}

void bit_vec::set(int idx) {
    bits[idx/32] |= (1 << (idx%32));
}

void bit_vec::clear(int idx) {
    bits[idx/32] &= ~(1 << (idx%32));
}

/* Clears *at least* upto bits */
void bit_vec::clear_range(int upto) {
    memset(bits, 0, (upto/32+1)*sizeof(*bits));
}
