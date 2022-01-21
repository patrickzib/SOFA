//
//  visualize_index.c
//  isaxlib
//
//  Created by Kostas Zoumpatianos and Botao Peng, March 2020
//
#include "../../config.h"
#include "../../globals.h"
#include <stdio.h>
#include <string.h>

#include "ads/isax_index.h"
#include "ads/isax_node.h"
#include "ads/isax_visualize_index.h"
#include "ads/sax/sax.h"

#define STRING_SIZE 5000
#define BUFFER_SIZE 3000

void bst_print_dot_null(void * key, int nullcount, FILE* stream)
{
    fprintf(stream, "    null%d [shape=point];\n", nullcount);
    fprintf(stream, "    %lu -> null%d;\n", (size_t)key, nullcount);
}

void bst_print_dot_aux(isax_node* node, FILE* stream)
{
    static int nullcount = 0;
    
    if (node->left_child)
    {
        fprintf(stream, "    %lu -> %lu;\n", (size_t)node, (size_t)node->left_child);
        bst_print_dot_aux(node->left_child, stream);
    }
    else
        bst_print_dot_null(node, nullcount++, stream);
    
    if (node->right_child)
    {
        fprintf(stream, "    %lu -> %lu;\n", (size_t)node, (size_t)node->right_child);
        bst_print_dot_aux(node->right_child, stream);
    }
    else
        bst_print_dot_null(node, nullcount++, stream);
}

void bst_print_dot(isax_node* tree, FILE* stream)
{
    //fprintf(stream, "digraph BST {\n");
    fprintf(stream, "    node [fontname=\"Arial\"];\n");
    
    if (!tree)
        fprintf(stream, "\n");
    else if (!tree->right_child && !tree->left_child)
        fprintf(stream, "    %lu;\n", (size_t)tree);
    else
        bst_print_dot_aux(tree, stream);
    
    //fprintf(stream, "}\n");
}

void isax_print_dot(isax_index* index, FILE *stream) 
{
    isax_node *node = index->first_node;
    fprintf(stream, "digraph BST {\n");
    if(node == NULL){
        fprintf(stream, "First node is empty");
    }
    while (node != NULL) {
        if(!node->is_leaf) {
            fprintf(stream, "    %lu -> %lu;\n", (size_t)index, (size_t)node);
            bst_print_dot(node, stream);
            fprintf(stream, "\n");
            
        }
        node = node->next;
    }
    fprintf(stream, "\n}");
    
}

void calculate_average_depth(FILE *ifile, isax_index *index)
{
    isax_node *node = index->first_node;
    if(node == NULL){
        return;
    };

    char logfile_out_number[STRING_SIZE] = "subtree,";
    char logfile_out_depth[STRING_SIZE] = "average depth,";
    char logfile_out_leaf_size[STRING_SIZE] = "average leaf size,";

    char buffer_number[BUFFER_SIZE] = "";
    char buffer_depth[BUFFER_SIZE] = "";
    char buffer_leaf_size[BUFFER_SIZE] = "";

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

            sprintf(buffer_number, "%d,", tree_counter);
            sprintf(buffer_depth, "%f,", current_depth);
            sprintf(buffer_leaf_size, "%f,", current_leaf_size);

            strcat(logfile_out_number, buffer_number);
            strcat(logfile_out_depth, buffer_depth);
            strcat(logfile_out_leaf_size, buffer_leaf_size);

        }

        node = node->next;
    }
    double total_depth_average = depth / (double)tree_counter;
    double total_leaf_size_average = (double)leaf_size_total / (double)leaf_counter_total;

    sprintf(buffer_number, "total\n");
    sprintf(buffer_depth, "%f\n", total_depth_average);
    sprintf(buffer_leaf_size, "%f\n", total_leaf_size_average);

    strcat(logfile_out_number, buffer_number);
    strcat(logfile_out_depth, buffer_depth);
    strcat(logfile_out_leaf_size, buffer_leaf_size);

    fprintf(ifile,"%s",logfile_out_number);
    fprintf(ifile,"%s",logfile_out_depth);
    fprintf(ifile,"%s",logfile_out_leaf_size);

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

