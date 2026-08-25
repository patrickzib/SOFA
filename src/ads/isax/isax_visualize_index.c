//
//  visualize_index.c
//  isaxlib
//
//  Created by Kostas Zoumpatianos and Botao Peng, March 2020
//
#include "config.h"
#include "../../../globals.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ads/isax_index.h"
#include "ads/isax_node.h"
#include "ads/isax_visualize_index.h"
#include "ads/sax/sax.h"

#define STRING_SIZE 256
#define BUFFER_SIZE 256


void calculate_average_depth(FILE *ifile, isax_index *index)
{
    isax_node *node = index->first_node;
    if(node == NULL){
        return;
    };

    char *buffer_number = malloc(BUFFER_SIZE * sizeof(char));
    char *buffer_depth = malloc(BUFFER_SIZE * sizeof(char));
    char *buffer_leaf_size = malloc(BUFFER_SIZE * sizeof(char));

    char logfile_out_header[STRING_SIZE] = "subtrees,average depth,average leaf size\n";
    char logfile_out_values[STRING_SIZE] = "";

    double depth = 0.0;
    unsigned long leaf_size_total = 0;
    int leaf_counter_total = 0;
    int tree_counter = 0;

    while (node != NULL) {

        double current_depth;
        double current_leaf_size;

        if(!node->is_leaf)
        {
            ++tree_counter;

            unsigned int total_depth = 0;
            unsigned int node_count = 0;

            unsigned int leaf_count = 0;
            unsigned long leaf_size = 0;

            traverse_tree(node, 0, &total_depth, &node_count, &leaf_count, &leaf_size);

            current_depth = ((double)total_depth / (double) node_count);
            current_leaf_size = (double) leaf_size / (double) leaf_count;

            depth += current_depth;
            leaf_size_total += leaf_size;
            leaf_counter_total += leaf_count;
        }

        node = node->next;
    }
    double total_depth_average = depth / (double)tree_counter;
    double total_leaf_size_average = (double)leaf_size_total / (double)leaf_counter_total;

    sprintf(buffer_number, "%d,", tree_counter);
    sprintf(buffer_depth, "%f,", total_depth_average);
    sprintf(buffer_leaf_size, "%f\n", total_leaf_size_average);

    strcat(logfile_out_values, buffer_number);
    strcat(logfile_out_values, buffer_depth);
    strcat(logfile_out_values, buffer_leaf_size);

    fprintf(ifile,"%s",logfile_out_header);
    fprintf(ifile,"%s",logfile_out_values);

    free(buffer_number);
    free(buffer_depth);
    free(buffer_leaf_size);

    return;
}

void traverse_tree(isax_node *node, unsigned int parent_depth, unsigned int *total_depth, unsigned int *node_count, unsigned int *leaf_count, unsigned long *leaf_size)
{
    if(node->is_leaf)
    {
        *leaf_count +=1;
        *leaf_size += node->leaf_size;

        *total_depth += parent_depth;
        *node_count +=1;
        return;
    }

    if(node->left_child)
    {
        traverse_tree(node->left_child, ++parent_depth, total_depth, node_count, leaf_count, leaf_size);
    }

    if(node->right_child)
    {
        traverse_tree(node->right_child, ++parent_depth, total_depth, node_count, leaf_count, leaf_size);
    }
    *total_depth += parent_depth;
    *node_count +=1;
    return;
}
