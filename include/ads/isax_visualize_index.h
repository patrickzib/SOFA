//
//  dot_exporter.h
//  isaxlib
//
//  Created by Kostas Zoumpatianos and Botao Peng, March 2020
//
#include "../../config.h"
#include "../../globals.h"
#include "isax_node.h"
#include "isax_index.h"

void bst_print_dot_null(void * key, int nullcount, FILE* stream);
void bst_print_dot_aux(isax_node* node, FILE* stream);
void bst_print_dot(isax_node* tree, FILE* stream);
void isax_print_dot(isax_index* index, FILE *stream);

void calculate_average_depth(FILE *ifile, isax_index *index);
void traverse_tree(isax_node *node, unsigned int parent_depth, unsigned int *total_depth, unsigned int *node_count, unsigned int *leaf_count, unsigned long *leaf_size);