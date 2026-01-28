//
//  isax_node_split.c
//  isaxlib
//
//  Created by Kostas Zoumpatianos and Botao Peng, March 2020
//
#include "config.h"
#include "../../globals.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <pthread.h>
#include <limits.h>

#include "ads/sax/sax.h"
#include "ads/sax/sax_breakpoints.h"
#include "ads/isax_node.h"
#include "ads/isax_index.h"
#include "ads/isax_node_split.h"
#include "ads/calc_utils.h"
#include "ads/inmemory_index_engine.h"

static void append_split_buffer(isax_node_record *split_buffer, int *split_buffer_index,
                                sax_type **sax_buffer, ts_type **ts_buffer,
                                file_position_type **pos_buffer, int *buffer_size,
                                int insertion_mode) {
    int size = *buffer_size;
    for (int i = size - 1; i >= 0; i--) {
        split_buffer[*split_buffer_index].sax = sax_buffer[i];
        split_buffer[*split_buffer_index].ts = ts_buffer ? ts_buffer[i] : NULL;
        split_buffer[*split_buffer_index].position = pos_buffer[i];
        split_buffer[*split_buffer_index].insertion_mode = insertion_mode;
        (*split_buffer_index)++;
    }
    *buffer_size = 0;
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

int simple_split_decision(isax_node_split_data *split_data,
                          isax_index_settings *settings) {
    if (split_data->split_mask[0] == 0)
        return 0;

    split_data->splitpoint = 1;
    while (split_data->splitpoint < settings->n_segments) {
        if (split_data->split_mask[split_data->splitpoint] <
            split_data->split_mask[split_data->splitpoint - 1]) {
            return split_data->splitpoint;
        }
        split_data->splitpoint++;
    }

    return 0;
}

int informed_split_decision(isax_node_split_data *split_data,
                            isax_index_settings *settings,
                            isax_node_record *records_buffer,
                            int records_buffer_size) {
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
    int segment_to_split_b = -1;
    for (i = 0; i < settings->n_segments; i++) {
        int new_bit_cardinality = split_data->split_mask[i] + 1;
        if (new_bit_cardinality > settings->sax_bit_cardinality - 1) {
            continue;
        }
        // TODO: Optimize this.
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
            if (abs(segment_mean[i] - b) <= abs(segment_mean[i] - segment_to_split_b)) {
                segment_to_split = i;
                segment_to_split_b = b;
            }
        }
    }

    free(segment_mean);
    free(segment_stdev);
    return segment_to_split;
}

int maxvar_split_decision(isax_node_split_data *split_data,
                          isax_index_settings *settings,
                          isax_node_record *records_buffer,
                          int records_buffer_size) {
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

int maxbin_split_decision(isax_node_split_data *split_data,
                          isax_index_settings *settings,
                          isax_node_record *records_buffer,
                          int records_buffer_size) {
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


void split_node(isax_index *index, isax_node *node, int inmemory) {
    int sktting;
    // ******************************************************* 
    // CREATE TWO NEW NODES AND SET OLD ONE AS AN INTERMEDIATE
    // ******************************************************* 
#ifdef DEBUG
    printf("*** Splitting. ***\n\n");
#endif

#ifdef DEBUG
    if (!node->is_leaf) {
        fprintf(stderr, "sanity error: You are trying to split something weird.\n");
    }
#endif
    // Create split_data for this node.
    isax_node_split_data *split_data = malloc(sizeof(isax_node_split_data));
    if (split_data == NULL) {
        fprintf(stderr, "error: could not allocate memory for node split data.\n");
    }
    split_data->split_mask = malloc(sizeof(sax_type) * index->settings->n_segments);
    if (split_data->split_mask == NULL) {
        fprintf(stderr, "error: could not allocate memory for node split mask.\n");
    }
    if (node->parent == NULL) {
        memset(split_data->split_mask, 0, sizeof(sax_type) * index->settings->n_segments);
        split_data->splitpoint = 0;
    } else {
        memcpy(split_data->split_mask, node->parent->split_data->split_mask,
               sizeof(sax_type) * index->settings->n_segments);
    }

    int fallback_splitpoint = -1;
    for (int s = 0; s < index->settings->n_segments; s++) {
        if (split_data->split_mask[s] + 1 <= index->settings->sax_bit_cardinality - 1) {
            fallback_splitpoint = s;
            break;
        }
    }
    if (fallback_splitpoint < 0) {
        fprintf(stderr, "error 1: cannot split in depth more than %d.\n", index->settings->sax_bit_cardinality);
        return;  // no split possible anymore
    }

    node->is_leaf = 0;
    node->leaf_size = 0;

    // TODO: needed???
    isax_node_mbb_reset(node, index->settings->timeseries_size);

    __sync_fetch_and_add(&(index->memory_info.mem_tree_structure), 2);

    isax_node *left_child = isax_leaf_node_init(index->settings->initial_leaf_buffer_size);
    isax_node *right_child = isax_leaf_node_init(index->settings->initial_leaf_buffer_size);
    left_child->is_leaf = 1;
    right_child->is_leaf = 1;
    left_child->parent = node;
    right_child->parent = node;
    node->split_data = split_data;
    node->left_child = left_child;
    node->right_child = right_child;


    // ############ S P L I T   D A T A #############
    // Allocating 1 more position to cover any off-sized allocations happening due to
    // trying to load one more record from a fetched file page which does not exist.
    // e.g. line 284 ( if(!fread... )
    isax_node_record *split_buffer = malloc(sizeof(isax_node_record) *
                                            (index->settings->max_leaf_size + 1));
    int split_buffer_index = 0;

    // ********************************************************
    // SPLIT SAX BUFFERS CONTAINED IN *RAM* AND PUT IN CHILDREN
    // ******************************************************** 
    // Split both sax and ts data and move to the new leafs

    append_split_buffer(split_buffer, &split_buffer_index,
                        node->buffer->full_sax_buffer, node->buffer->full_ts_buffer,
                        node->buffer->full_position_buffer, &node->buffer->full_buffer_size,
                        NO_TMP | FULL);
    append_split_buffer(split_buffer, &split_buffer_index,
                        node->buffer->partial_sax_buffer, NULL,
                        node->buffer->partial_position_buffer, &node->buffer->partial_buffer_size,
                        NO_TMP | PARTIAL);
    append_split_buffer(split_buffer, &split_buffer_index,
                        node->buffer->tmp_full_sax_buffer, node->buffer->tmp_full_ts_buffer,
                        node->buffer->tmp_full_position_buffer, &node->buffer->tmp_full_buffer_size,
                        TMP | FULL);
    append_split_buffer(split_buffer, &split_buffer_index,
                        node->buffer->tmp_partial_sax_buffer, NULL,
                        node->buffer->tmp_partial_position_buffer, &node->buffer->tmp_partial_buffer_size,
                        TMP | PARTIAL);

    destroy_node_buffer(node->buffer);
    node->buffer = NULL;

    // *****************************************************
    // SPLIT BUFFERS CONTAINED ON *DISK* AND PUT IN CHILDREN
    // ***************************************************** 

    // File is split in two files, but it is not
    // removed from disk. It is going to be used in the end.
    if (!inmemory && node->filename != NULL) {
        char *full_fname = malloc(sizeof(char) * (strlen(node->filename) + 6));
        strcpy(full_fname, node->filename);
        strcat(full_fname, ".full");
        //COUNT_INPUT_TIME_START
        FILE *full_file = fopen(full_fname, "r");
        //COUNT_INPUT_TIME_END

        // If it can't open exit;
        if (full_file != NULL) {
#ifdef DEBUG
            printf("*** Splitting: %s\n\n", full_fname);
#endif
            //COUNT_INPUT_TIME_START
            while (!feof(full_file)) {
                split_buffer[split_buffer_index].position = malloc(index->settings->position_byte_size);
                split_buffer[split_buffer_index].sax = malloc(index->settings->sax_byte_size);
                split_buffer[split_buffer_index].ts = malloc(index->settings->ts_byte_size);
                split_buffer[split_buffer_index].insertion_mode = FULL | TMP;

                // If it can't read continue.
                if (!fread(split_buffer[split_buffer_index].position, sizeof(file_position_type),
                           1, full_file)) {
                    // Free because it is not inserted in the tree
                    free(split_buffer[split_buffer_index].position);
                } else {
                    if (!fread(split_buffer[split_buffer_index].sax, sizeof(sax_type),
                               index->settings->n_segments, full_file)) {
                        // Free because it is not inserted in the tree
                        free(split_buffer[split_buffer_index].position);
                        free(split_buffer[split_buffer_index].sax);
                        free(split_buffer[split_buffer_index].ts);
                    } else {
                        if (!fread(split_buffer[split_buffer_index].sax, sizeof(ts_type),
                                   index->settings->timeseries_size, full_file)) {
                            // Free because it is not inserted in the tree
                            free(split_buffer[split_buffer_index].position);
                            free(split_buffer[split_buffer_index].sax);
                            free(split_buffer[split_buffer_index].ts);
                        } else {
                            // Increase leaf size (from 0) so that we keep track how many raw time series we 
                            // have to move in the finalization step.
                            node->leaf_size++;
                            split_buffer_index++;
                            index->allocated_memory += index->settings->full_record_size;
                        }
                    }
                }
            }
            //COUNT_INPUT_TIME_END

#ifdef DEBUG
            printf("*** END OF: %s\n\n", full_fname);
#endif
            COUNT_OUTPUT_TIME_START
            remove(full_fname);
            COUNT_OUTPUT_TIME_END
            //COUNT_INPUT_TIME_START
            fclose(full_file);
            //COUNT_INPUT_TIME_END
        }
        free(full_fname);

        char *partial_fname = malloc(sizeof(char) * (strlen(node->filename) + 6));
        strcpy(partial_fname, node->filename);
        strcat(partial_fname, ".part");
        //COUNT_INPUT_TIME_START
        FILE *partial_file = fopen(partial_fname, "r");
        //COUNT_INPUT_TIME_END

        // If it can't open exit;
        if (partial_file != NULL) {
#ifdef DEBUG
            printf("*** Splitting: %s\n\n", partial_fname);
#endif
            //COUNT_INPUT_TIME_START

            while (!feof(partial_file)) {
                split_buffer[split_buffer_index].position = malloc(index->settings->position_byte_size);
                split_buffer[split_buffer_index].sax = malloc(index->settings->sax_byte_size);
                split_buffer[split_buffer_index].insertion_mode = PARTIAL | TMP;
                // If it can't read continue.
                if (!fread(split_buffer[split_buffer_index].position, sizeof(file_position_type),
                           1, partial_file)) {
                    // Free because it is not inserted in the tree
                    free(split_buffer[split_buffer_index].position);
                    free(split_buffer[split_buffer_index].sax);

                    continue;
                } else {
                    if (!fread(split_buffer[split_buffer_index].sax, sizeof(sax_type),
                               index->settings->n_segments, partial_file)) {
                        // Free because it is not inserted in the tree
                        free(split_buffer[split_buffer_index].position);
                        free(split_buffer[split_buffer_index].sax);
                        continue;
                    } else {
                        node->leaf_size++;
                        split_buffer_index++;
                        index->allocated_memory += index->settings->partial_record_size;
                    }
                }
            }
            //COUNT_INPUT_TIME_END
            COUNT_OUTPUT_TIME_START
            remove(partial_fname);
            COUNT_OUTPUT_TIME_END
            //COUNT_INPUT_TIME_START
            //printf("this is before sktting\n");

            sktting = fclose(partial_file);
            //printf("this is skating%d\n",sktting);
            //COUNT_INPUT_TIME_END
        }

        free(partial_fname);
    }

    //printf("sizeof split buffer: %d\n", split_buffer_index);
    int selected_splitpoint = select_split_point(
        split_data,
        index->settings,
        split_buffer,
        split_buffer_index);

    if (selected_splitpoint > -1) {
        split_data->splitpoint = selected_splitpoint;
    }
    else {
        split_data->splitpoint = fallback_splitpoint;
    }


    if (split_data->splitpoint < 0) {
        fprintf(stderr, "error 1: cannot split in depth more than %d.\n",
                index->settings->sax_bit_cardinality);
        exit(-1);
    }

    if (++split_data->split_mask[split_data->splitpoint] > index->settings->sax_bit_cardinality - 1) {
        fprintf(stderr, "error 2: cannot split in depth more than %d.\n",
                index->settings->sax_bit_cardinality);
        exit(-1);
    }

    root_mask_type mask = index->settings->bit_masks[index->settings->sax_bit_cardinality -
                                                     split_data->split_mask[split_data->splitpoint] - 1];


    while (split_buffer_index > 0) {
        split_buffer_index--;
        if (mask & split_buffer[split_buffer_index].sax[split_data->splitpoint]) {
            add_record_to_node(index, right_child, &split_buffer[split_buffer_index], 1);
        } else {
            add_record_to_node(index, left_child, &split_buffer[split_buffer_index], 1);
        }
    }

    free(split_buffer);
}

void split_node_inmemory(isax_index *index, isax_node *node) {
    split_node(index, node, 1);
}
