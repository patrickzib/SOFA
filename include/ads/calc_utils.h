#ifndef MESSI_SFA_CALC_UTILS_H
#define MESSI_SFA_CALC_UTILS_H

#include "../../config.h"
#include "../../globals.h"

////// Utility functions ////
// Function to calculate mean of an array of floats
float calculateMean(ts_type *data, int n);

// Function to calculate standard deviation of an array of floats
float calculateStdDev(ts_type *data, int n, ts_type mean);

// Function to perform zero mean normalization
void znorm(ts_type *data, int n);

#endif //MESSI_SFA_CALC_UTILS_H
