#include "config.h"
#include "../../globals.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "math.h"

#include "ads/calc_utils.h"


////// Utility functions ////
// Function to calculate mean of an array of floats
float calculateMean(ts_type *data, int n) {
    float sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += data[i];
    }
    return sum / n;
}

// Function to calculate standard deviation of an array of floats
float calculateStdDev(ts_type *data, int n, ts_type mean) {
    ts_type sumSquaredDiffs = 0.0;
    for (int i = 0; i < n; i++) {
        sumSquaredDiffs += (data[i] - mean) * (data[i] - mean);
    }
    return sqrt(sumSquaredDiffs / n);
}

// Function to perform zero mean normalization
void znorm(ts_type *data, int n) {
    // printf("No nor applied");
    // Calculate mean
    ts_type mean = calculateMean(data, n);

    // Calculate standard deviation
    ts_type stdDev = calculateStdDev(data, n, mean);

    // Normalize each data point
    if (stdDev < 1e-8) {
        stdDev = 1.0;
    }

    for (int i = 0; i < n; i++) {
        data[i] = (data[i] - mean) / stdDev;
    }
}

void isax_node_mbb_reset(isax_node *node, int size) {
    if (node->mbb_min == NULL) {
        node->mbb_min = malloc(sizeof(ts_type) * (size_t) size);
    }
    if (node->mbb_max == NULL) {
        node->mbb_max = malloc(sizeof(ts_type) * (size_t) size);
    }
    if (node->mbb_min == NULL || node->mbb_max == NULL) {
        fprintf(stderr, "error: could not allocate MBB arrays.\n");
        return;
    }
    for (int i = 0; i < size; ++i) {
        node->mbb_min[i] = FLT_MAX;
        node->mbb_max[i] = -FLT_MAX;
    }
    node->mbb_valid = 0;
}

void isax_node_mbb_update(isax_node *node, const ts_type *ts, int size) {
    if (node->mbb_min == NULL || node->mbb_max == NULL) {
        isax_node_mbb_reset(node, size);
    }
    node->mbb_valid = 1;
    for (int i = 0; i < size; ++i) {
        ts_type value = ts[i];
        if (value < node->mbb_min[i]) {
            node->mbb_min[i] = value;
        }
        if (value > node->mbb_max[i]) {
            node->mbb_max[i] = value;
        }
    }
}

void isax_node_mbb_update_upwards(isax_node *node, const ts_type *ts, int size) {
    for (isax_node *current = node; current != NULL; current = current->parent) {
        isax_node_mbb_update(current, ts, size);
    }
}

ts_type ts_mbb_distance_sq(const ts_type *ts, const ts_type *mbb_min, const ts_type *mbb_max,
                           int size, ts_type bound, ts_type ratio_sqrt) {
    ts_type distance = 0.0f;
    ts_type scaled_bound = bound / ratio_sqrt;

    for (int i = 0; i < size; ++i) {
        ts_type value = ts[i];
        if (value < mbb_min[i]) {
            ts_type diff = mbb_min[i] - value;
            distance += diff * diff;
        } else if (value > mbb_max[i]) {
            ts_type diff = value - mbb_max[i];
            distance += diff * diff;
        }
        if (distance >= scaled_bound) {
            return distance * ratio_sqrt;
        }
    }

    return distance * ratio_sqrt;
}
