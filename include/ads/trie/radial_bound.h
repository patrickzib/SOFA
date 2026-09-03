#ifndef MESSI_TRIE_RADIAL_BOUND_H
#define MESSI_TRIE_RADIAL_BOUND_H

#include <float.h>
#include <math.h>

typedef struct {
    double lower;
    double upper;
} messi_distance_interval;

/* Bound the exact norm represented by a double-precision sum of squared float
 * differences.  The gamma term covers subtraction, multiplication, and
 * accumulation rounding; nextafter also covers the final square root. */
static inline messi_distance_interval messi_distance_interval_from_squared(double squared,
                                                                            int terms) {
    if (squared <= 0.0)
        return (messi_distance_interval) { 0.0, 0.0 };
    const double operations = terms > 0 ? 3.0 * (double) terms + 2.0 : 2.0;
    const double product = operations * DBL_EPSILON;
    const double gamma = product < 0.5 ? product / (1.0 - product) : 0.5;
    const double lower_squared = squared / (1.0 + gamma);
    const double upper_squared = squared / (1.0 - gamma);
    messi_distance_interval result = {
        nextafter(sqrt(lower_squared), 0.0),
        nextafter(sqrt(upper_squared), INFINITY)
    };
    return result;
}

/* The stored radius is the nearest float to the build-time double result.
 * Expand it by one float ULP before applying the reverse triangle inequality. */
static inline float messi_radial_lower_bound_squared(messi_distance_interval query_radius,
                                                      float stored_record_radius) {
    double record_lower = (double) nextafterf(stored_record_radius, -INFINITY);
    const double record_upper = (double) nextafterf(stored_record_radius, INFINITY);
    if (record_lower < 0.0) record_lower = 0.0;
    double gap = 0.0;
    if (query_radius.lower > record_upper)
        gap = query_radius.lower - record_upper;
    else if (record_lower > query_radius.upper)
        gap = record_lower - query_radius.upper;
    if (gap <= 0.0) return 0.0f;
    return nextafterf((float) (gap * gap), 0.0f);
}

#endif
