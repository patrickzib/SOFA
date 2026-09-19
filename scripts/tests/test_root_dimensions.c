#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ads/calc_utils.h"

static int expect_mask(const char *name, root_mask_type actual, root_mask_type expected) {
    if (actual == expected) return 1;
    fprintf(stderr, "%s mask mismatch: expected 0x%llx, got 0x%llx\n",
            name, (unsigned long long) expected, (unsigned long long) actual);
    return 0;
}

int main(void) {
    isax_index index;
    isax_index_settings settings;
    root_mask_type bit_masks[33];
    int root_dimensions[16];
    sax_type sax[32];
    memset(&index, 0, sizeof(index));
    memset(&settings, 0, sizeof(settings));
    memset(sax, 0, sizeof(sax));
    for (int bit = 0; bit < 33; ++bit) bit_masks[bit] = (root_mask_type) 1 << bit;
    settings.n_segments = 32;
    settings.isax_index_segments = 16;
    settings.sax_bit_cardinality = 8;
    settings.bit_masks = bit_masks;
    settings.root_dimensions = root_dimensions;
    index.settings = &settings;

    for (int slot = 0; slot < 16; ++slot) root_dimensions[slot] = slot * 2;
    sax[0] = 0x80;
    sax[6] = 0x80;
    sax[30] = 0x80;
    if (!expect_mask("uniform SAX", isax_root_mask_from_sax(&index, sax, 1),
                     ((root_mask_type) 1 << 15) | ((root_mask_type) 1 << 12) | 1))
        return 1;

    double variance[32];
    for (int dimension = 0; dimension < 32; ++dimension)
        variance[dimension] = (double) dimension;
    if (configure_root_dimensions(&settings, variance, 32) != SUCCESS) return 1;
    for (int slot = 0; slot < 16; ++slot) {
        if (root_dimensions[slot] != 31 - slot) {
            fprintf(stderr, "variance rank mismatch at slot %d: %d\n",
                    slot, root_dimensions[slot]);
            return 1;
        }
    }
    memset(sax, 0, sizeof(sax));
    sax[31] = 0x80;
    sax[19] = 0x80;
    sax[16] = 0x80;
    if (!expect_mask("variance-ranked SFA", isax_root_mask_from_sax(&index, sax, 1),
                     ((root_mask_type) 1 << 15) | ((root_mask_type) 1 << 3) | 1))
        return 1;
    if (isax_root_mask_from_sax(&index, sax, 1) >= ((root_mask_type) 1 << 16)) {
        fprintf(stderr, "root mask exceeds the fixed 2^16 root table\n");
        return 1;
    }
    return 0;
}
