//
//  dot_exporter.h
//  isaxlib
//
//  Created by Kostas Zoumpatianos and Botao Peng, March 2020
//
#include "config.h"
#include "../../globals.h"
#include "isax_node.h"
#include "isax_index.h"

void calculate_average_depth(FILE *ifile, isax_index *index);
void traverse_tree(isax_node *node, unsigned int parent_depth, unsigned int *total_depth, unsigned int *node_count, unsigned int *leaf_count, unsigned long *leaf_size);
