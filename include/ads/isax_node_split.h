//
//  isax_node_split.h
//  isaxlib
//
//  Created by Kostas Zoumpatianos and Botao Peng, March 2020
//

#ifndef isaxlib_isax_node_split_h
#define isaxlib_isax_node_split_h
#include "config.h"
#include "../../globals.h"
#include "isax_index.h"
#include "isax_node.h"

int split_node(isax_index *index, isax_node *node, int inmemory, int kn);

#endif
