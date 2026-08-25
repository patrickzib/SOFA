//
//  isax_node_split.c
//  isaxlib
//
//  Created by Kostas Zoumpatianos and Botao Peng, March 2020
//
#include "config.h"
#include "../../../globals.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <errno.h>

#include "ads/sax/sax.h"
#include "ads/sax/sax_breakpoints.h"
#include "ads/isax_node.h"
#include "ads/isax_index.h"
#include "ads/isax_node_split.h"
#include "ads/calc_utils.h"
#include "ads/inmemory_index_engine.h"

static int simple_split_decision(isax_node_split_data *split_data,
                                 isax_index_settings *settings);
static int informed_split_decision(isax_node_split_data *split_data,
                                   isax_index_settings *settings,
                                   isax_node_record *records_buffer,
                                   int records_buffer_size);
static int maxvar_split_decision(isax_node_split_data *split_data,
                                 isax_index_settings *settings,
                                 isax_node_record *records_buffer,
                                 int records_buffer_size);
static int maxbin_split_decision(isax_node_split_data *split_data,
                                 isax_index_settings *settings,
                                 isax_node_record *records_buffer,
                                 int records_buffer_size);

static int ensure_split_capacity(isax_node_record **split_buffer, int *capacity, int needed) {
    if (needed <= *capacity) {
        return 1;
    }
    int new_capacity = *capacity;
    while (new_capacity < needed) {
        new_capacity *= 2;
    }
    isax_node_record *resized = realloc(*split_buffer, sizeof(isax_node_record) * (size_t) new_capacity);
    if (resized == NULL) {
        fprintf(stderr, "error: could not grow split buffer.\n");
        return 0;
    }
    *split_buffer = resized;
    *capacity = new_capacity;
    return 1;
}

static int append_split_buffer(isax_node_record **split_buffer, int *split_buffer_index,
                               sax_type **sax_buffer, ts_type **ts_buffer,
                               file_position_type **pos_buffer, int *buffer_size,
                               int insertion_mode, int *capacity) {
    int size = *buffer_size;
    if (size <= 0) {
        return 1;
    }
    if (!ensure_split_capacity(split_buffer, capacity, *split_buffer_index + size)) {
        return 0;
    }
    for (int i = size - 1; i >= 0; i--) {
        (*split_buffer)[*split_buffer_index].sax = sax_buffer[i];
        (*split_buffer)[*split_buffer_index].ts = ts_buffer ? ts_buffer[i] : NULL;
        (*split_buffer)[*split_buffer_index].position = pos_buffer[i];
        (*split_buffer)[*split_buffer_index].insertion_mode = insertion_mode;
        (*split_buffer)[*split_buffer_index].destination = NULL;
        (*split_buffer_index)++;
    }
    *buffer_size = 0;
    return 1;
}

static int select_split_point(isax_node_split_data *split_data,
                              isax_index_settings *settings,
                              isax_node_record *records_buffer,
                              int records_buffer_size) {
    switch (settings->node_split_criterion) {
        case 1:
            return informed_split_decision(split_data, settings, records_buffer, records_buffer_size);
        case 2:
            return simple_split_decision(split_data, settings);
        case 3:
            return maxvar_split_decision(split_data, settings, records_buffer, records_buffer_size);
        case 4:
            return maxbin_split_decision(split_data, settings, records_buffer, records_buffer_size);
        default:
            return informed_split_decision(split_data, settings, records_buffer, records_buffer_size);
    }
}

static int append_disk_file(isax_index *index, isax_node *node,
                            isax_node_record **split_buffer, int *split_buffer_index,
                            int *capacity, const char *filename, enum insertion_mode mode) {
    FILE *file = fopen(filename, "rb");
    if (file == NULL) {
        if (errno == ENOENT) {
            return 1;
        }
        fprintf(stderr, "error: could not open split data file %s.\n", filename);
        return 0;
    }

#ifdef DEBUG
    printf("*** Splitting: %s\n\n", filename);
#endif

    const int has_ts = (mode & FULL) != 0;
    for (;;) {
        if (!ensure_split_capacity(split_buffer, capacity, *split_buffer_index + 1)) {
            fclose(file);
            return 0;
        }

        isax_node_record *record = &(*split_buffer)[*split_buffer_index];
        record->position = malloc(index->settings->position_byte_size);
        record->sax = malloc(index->settings->sax_byte_size);
        record->ts = has_ts ? malloc(index->settings->ts_byte_size) : NULL;
        record->destination = NULL;
        record->insertion_mode = mode;
        if (record->position == NULL || record->sax == NULL || (has_ts && record->ts == NULL)) {
            fprintf(stderr, "error: could not allocate split record.\n");
            free(record->position);
            free(record->sax);
            free(record->ts);
            fclose(file);
            return 0;
        }

        if (fread(record->position, sizeof(file_position_type), 1, file) != 1) {
            free(record->position);
            free(record->sax);
            free(record->ts);
            if (feof(file)) {
                break;
            }
            fprintf(stderr, "error: could not read split record position from %s.\n", filename);
            fclose(file);
            return 0;
        }
        if (fread(record->sax, sizeof(sax_type), index->settings->n_segments, file) !=
                (size_t)index->settings->n_segments ||
            (has_ts && fread(record->ts, sizeof(ts_type), index->settings->timeseries_size, file) !=
                (size_t)index->settings->timeseries_size)) {
            fprintf(stderr, "error: truncated split record in %s.\n", filename);
            free(record->position);
            free(record->sax);
            free(record->ts);
            fclose(file);
            return 0;
        }

        node->leaf_size++;
        (*split_buffer_index)++;
        index->allocated_memory += has_ts ? index->settings->full_record_size
                                          : index->settings->partial_record_size;
    }

    fclose(file);
#ifdef DEBUG
    printf("*** END OF: %s\n\n", filename);
#endif
    remove(filename);
    return 1;
}

static int append_disk_buffers(isax_index *index, isax_node *node,
                               isax_node_record **split_buffer, int *split_buffer_index,
                               int *capacity) {
    if (node->filename == NULL) {
        return 1;
    }

    size_t filename_size = strlen(node->filename) + 6;
    char *filename = malloc(filename_size);
    if (filename == NULL) {
        fprintf(stderr, "error: could not allocate split filename.\n");
        return 0;
    }

    snprintf(filename, filename_size, "%s.full", node->filename);
    int success = append_disk_file(index, node, split_buffer, split_buffer_index,
                                   capacity, filename, FULL | TMP);
    if (success) {
        snprintf(filename, filename_size, "%s.part", node->filename);
        success = append_disk_file(index, node, split_buffer, split_buffer_index,
                                   capacity, filename, PARTIAL | TMP);
    }
    free(filename);
    return success;
}

static int simple_split_decision(isax_node_split_data *split_data,
                                 isax_index_settings *settings) {
    int min_index = -1;
    for (int i = 0; i < settings->n_segments; i++) {
        if (split_data->split_mask[i] + 1 > settings->sax_bit_cardinality - 1) {
            continue;
        }
        if (min_index == -1 || split_data->split_mask[i] < split_data->split_mask[min_index]) {
            min_index = i;
        }
    }
    if (min_index == -1) {
        fprintf(stderr, "split_mask (n_segments=%d):", settings->n_segments);
        for (int i = 0; i < settings->n_segments; i++) {
            fprintf(stderr, " %u", (unsigned) split_data->split_mask[i]);
        }
        fprintf(stderr, "\n");
    }
    return min_index;
}

static int informed_split_decision(isax_node_split_data *split_data,
                                   isax_index_settings *settings,
                                   isax_node_record *records_buffer,
                                   int records_buffer_size) {
    if (records_buffer_size == 0) {
        return simple_split_decision(split_data, settings);
    }
    double *segment_mean = malloc(sizeof(double) * settings->n_segments);
    double *segment_stdev = malloc(sizeof(double) * settings->n_segments);

    int i, j;
    for (i = 0; i < settings->n_segments; i++) {
        segment_mean[i] = 0;
        segment_stdev[i] = 0;
    }
    for (i = 0; i < records_buffer_size; i++) {
        for (j = 0; j < settings->n_segments; j++) {
            segment_mean[j] += (int) records_buffer[i].sax[j];
        }
    }
    for (i = 0; i < settings->n_segments; i++) {
        segment_mean[i] /= (records_buffer_size);
    }
    for (i = 0; i < records_buffer_size; i++) {
        for (j = 0; j < settings->n_segments; j++) {
            segment_stdev[j] += pow(segment_mean[j] - (int) records_buffer[i].sax[j], 2);
        }
    }
    for (i = 0; i < settings->n_segments; i++) {
        segment_stdev[i] = sqrt(segment_stdev[i] / (records_buffer_size));
    }


    // Decide split point based on the above calculations
    int segment_to_split = -1;
    float segment_to_split_b = 0.0f;
    for (i = 0; i < settings->n_segments; i++) {
        int new_bit_cardinality = split_data->split_mask[i] + 1;
        if (new_bit_cardinality > settings->sax_bit_cardinality - 1) {
            continue;
        }

        // Calculate break point for new cardinality, a bit complex.
        int break_point_id = records_buffer[0].sax[i];
        break_point_id = (break_point_id >> ((settings->sax_bit_cardinality) - (new_bit_cardinality))) << 1;

        int new_cardinality = pow(2, new_bit_cardinality + 1);
        int right_offset = ((new_cardinality - 1) * (new_cardinality - 2)) / 2 + new_cardinality - 2;
        float b = sax_breakpoints[right_offset - break_point_id];

        if (segment_to_split == -1) {
            segment_to_split = i;
            segment_to_split_b = b;
            continue;
        }

        float left_range = segment_mean[i] - (3 * segment_stdev[i]);
        float right_range = segment_mean[i] + (3 * segment_stdev[i]);

        if (left_range <= b && b <= right_range) {
            if (fabs(segment_mean[i] - b) <= fabs(segment_mean[i] - segment_to_split_b)) {
                segment_to_split = i;
                segment_to_split_b = b;
            }
        }
    }

    free(segment_mean);
    free(segment_stdev);
    return segment_to_split;
}

static int maxvar_split_decision(isax_node_split_data *split_data,
                                 isax_index_settings *settings,
                                 isax_node_record *records_buffer,
                                 int records_buffer_size) {
    if (records_buffer_size == 0) {
        return simple_split_decision(split_data, settings);
    }
    int best_segment = -1;
    double best_variance = -1.0;

    for (int segment = 0; segment < settings->n_segments; ++segment) {
        if (split_data->split_mask[segment] + 1 > settings->sax_bit_cardinality - 1) {
            continue;
        }

        double mean = 0.0;
        for (int i = 0; i < records_buffer_size; ++i) {
            mean += (double) records_buffer[i].sax[segment];
        }
        mean /= (double) records_buffer_size;

        double variance = 0.0;
        for (int i = 0; i < records_buffer_size; ++i) {
            double diff = (double) records_buffer[i].sax[segment] - mean;
            variance += diff * diff;
        }
        variance /= (double) records_buffer_size;

        if (variance > best_variance) {
            best_variance = variance;
            best_segment = segment;
        }
    }

    return best_segment;
}

static int maxbin_split_decision(isax_node_split_data *split_data,
                                 isax_index_settings *settings,
                                 isax_node_record *records_buffer,
                                 int records_buffer_size) {
    if (records_buffer_size == 0) {
        return simple_split_decision(split_data, settings);
    }
    int best_segment = -1;
    int best_imbalance = INT_MAX;

    for (int segment = 0; segment < settings->n_segments; ++segment) {
        if (split_data->split_mask[segment] + 1 > settings->sax_bit_cardinality - 1) {
            continue;
        }

        int new_bit_cardinality = split_data->split_mask[segment] + 1;
        root_mask_type mask = settings->bit_masks[settings->sax_bit_cardinality - new_bit_cardinality - 1];

        int left_count = 0;
        int right_count = 0;
        for (int i = 0; i < records_buffer_size; ++i) {
            if (records_buffer[i].sax[segment] & mask) {
                right_count++;
            } else {
                left_count++;
            }
        }

        int imbalance = abs(left_count - right_count);
        if (imbalance < best_imbalance) {
            best_imbalance = imbalance;
            best_segment = segment;
        }
    }

    return best_segment;
}


int split_node(isax_index *index, isax_node *node, int inmemory, int kn) {
    if (node->mbb_sax_valid && node->mbb_sax_min != NULL && node->mbb_sax_max != NULL) {
        int can_split = 0;
        for (int i = 0; i < index->settings->n_segments; ++i) {
            if (node->mbb_sax_min[i] != node->mbb_sax_max[i]) {
                can_split = 1;
                break;
            }
        }
        if (!can_split) {
            return 0;
        }
    }

    // *******************************************************
    // CREATE TWO NEW NODES AND SET OLD ONE AS AN INTERMEDIATE
    // *******************************************************
    if (!node->is_leaf) {
        fprintf(stderr, "sanity error: You are trying to split something weird.\n");
        return 0;
    }

    // Create split_data for this node.
    isax_node_split_data *split_data = malloc(sizeof(isax_node_split_data));
    if (split_data == NULL) {
        fprintf(stderr, "error: could not allocate memory for node split data.\n");
        return 0;
    }

    split_data->split_mask = calloc(index->settings->n_segments, sizeof(sax_type));
    if (split_data->split_mask != NULL && node->parent != NULL) {
        memcpy(split_data->split_mask,
               node->parent->split_data->split_mask,
               sizeof(sax_type) * index->settings->n_segments);
    } else if (split_data->split_mask != NULL) {
        int root_segments = index->settings->n_segments / kn;
        for (int i = 0; i < root_segments; i++) {
            split_data->split_mask[i] = (sax_type) (kn - 1);
        }
    }
    if (split_data->split_mask == NULL) {
        fprintf(stderr, "error: could not allocate memory for split mask.\n");
        free(split_data);
        return 0;
    }

    __sync_fetch_and_add(&(index->memory_info.mem_tree_structure), 2);

    isax_node *left_child = isax_leaf_node_init(index->settings->initial_leaf_buffer_size, node);
    isax_node *right_child = isax_leaf_node_init(index->settings->initial_leaf_buffer_size, node);
    node->split_data = split_data;
    node->left_child = left_child;
    node->right_child = right_child;
    node->is_leaf = 0;
    node->leaf_size = 0;


    int split_buffer_capacity = node->buffer->full_buffer_size + node->buffer->partial_buffer_size +
                                node->buffer->tmp_full_buffer_size + node->buffer->tmp_partial_buffer_size + 1;
    if (split_buffer_capacity < 1) {
        split_buffer_capacity = 1;
    }
    isax_node_record *split_buffer = malloc(sizeof(isax_node_record) * (size_t) split_buffer_capacity);
    if (split_buffer == NULL) {
        fprintf(stderr, "fatal error: could not allocate split buffer.\n");
        exit(EXIT_FAILURE);
    }
    int split_buffer_index = 0;

    // ********************************************************
    // SPLIT SAX BUFFERS CONTAINED IN *RAM* AND PUT IN CHILDREN
    // ********************************************************
    // Split both sax and ts data and move to the new leafs

    if (!append_split_buffer(&split_buffer, &split_buffer_index,
                             node->buffer->full_sax_buffer, node->buffer->full_ts_buffer,
                             node->buffer->full_position_buffer, &node->buffer->full_buffer_size,
                             NO_TMP | FULL, &split_buffer_capacity) ||
        !append_split_buffer(&split_buffer, &split_buffer_index,
                             node->buffer->partial_sax_buffer, NULL,
                             node->buffer->partial_position_buffer, &node->buffer->partial_buffer_size,
                             NO_TMP | PARTIAL, &split_buffer_capacity) ||
        !append_split_buffer(&split_buffer, &split_buffer_index,
                             node->buffer->tmp_full_sax_buffer, node->buffer->tmp_full_ts_buffer,
                             node->buffer->tmp_full_position_buffer, &node->buffer->tmp_full_buffer_size,
                             TMP | FULL, &split_buffer_capacity) ||
        !append_split_buffer(&split_buffer, &split_buffer_index,
                             node->buffer->tmp_partial_sax_buffer, NULL,
                             node->buffer->tmp_partial_position_buffer, &node->buffer->tmp_partial_buffer_size,
                             TMP | PARTIAL, &split_buffer_capacity)) {
        free(split_buffer);
        return 0;
    }

    destroy_node_buffer(node->buffer);
    node->buffer = NULL;

    // *****************************************************
    // SPLIT BUFFERS CONTAINED ON *DISK* AND PUT IN CHILDREN
    // *****************************************************
    if (!inmemory) {
        if (!append_disk_buffers(index, node, &split_buffer, &split_buffer_index, &split_buffer_capacity)) {
            free(split_buffer);
            return 0;
        }
    }

    split_data->splitpoint = select_split_point(
        split_data,
        index->settings,
        split_buffer,
        split_buffer_index);

    if (split_data->splitpoint < 0 ||
        split_data->split_mask[split_data->splitpoint] + 1 > index->settings->sax_bit_cardinality - 1) {
        fprintf(stderr, "fallback: cannot split in depth more than %d.\n", index->settings->sax_bit_cardinality);
        int simple_split = simple_split_decision(split_data, index->settings);
        fprintf(stderr, "fallback: %d.\n", simple_split);
        if (simple_split < 0) {
            fprintf(stderr, "fallback fatal: %d.\n", simple_split);
        }
        split_data->splitpoint = simple_split;
    }

    if (split_data->splitpoint == -1 ||
        ++split_data->split_mask[split_data->splitpoint] > index->settings->sax_bit_cardinality - 1) {
        fprintf(stderr, "fatal error: cannot split in depth more than %d.\n", index->settings->sax_bit_cardinality);
        exit(-1);
    }

    root_mask_type mask = index->settings->bit_masks[
        index->settings->sax_bit_cardinality - split_data->split_mask[split_data->splitpoint] - 1];

    while (split_buffer_index > 0) {
        split_buffer_index--;
        if (mask & split_buffer[split_buffer_index].sax[split_data->splitpoint]) {
            add_record_to_node(index, right_child, &split_buffer[split_buffer_index], 1, kn);
        } else {
            add_record_to_node(index, left_child, &split_buffer[split_buffer_index], 1, kn);
        }
    }

    free(split_buffer);
    return 1;
}
