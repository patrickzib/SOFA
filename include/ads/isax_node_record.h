//
//  isax_node_record.h
//  isaxlib
//
//  Created by Kostas Zoumpatianos and Botao Peng, March 2020
//

#ifndef isaxlib_isax_node_record_h
#define isaxlib_isax_node_record_h
#include "config.h"
#include "../../globals.h"

typedef struct {
    sax_type *sax;
    ts_type *ts;
    file_position_type *position;
    enum insertion_mode insertion_mode;
    void *destination;
} isax_node_record;

#endif
