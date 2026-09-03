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
 * Expanding it by one float ULP covers both that conversion and the much
 * smaller double-precision accumulation error. */
static inline double messi_stored_radius_lower(float stored_record_radius) {
    const double lower = (double) nextafterf(stored_record_radius, -INFINITY);
    return lower > 0.0 ? lower : 0.0;
}

static inline double messi_stored_radius_upper(float stored_record_radius) {
    return (double) nextafterf(stored_record_radius, INFINITY);
}

/* Convert the squared BSF to a distance threshold, rounded outward.  Query
 * code uses this once per binary-searched window, never once per record. */
static inline double messi_bsf_radius_upper(float bsf) {
    if (!(bsf >= 0.0f) || isinf(bsf)) return INFINITY;
    return nextafter(sqrt((double) bsf), INFINITY);
}

/* Reference predicate used by numerical tests.  The query path implements
 * the same inequalities with two binary searches over sorted radii. */
static inline int messi_radial_radius_may_pass(messi_distance_interval query_radius,
                                               float stored_record_radius, float bsf) {
    const double threshold = messi_bsf_radius_upper(bsf);
    return messi_stored_radius_upper(stored_record_radius) >= query_radius.lower - threshold &&
           messi_stored_radius_lower(stored_record_radius) <= query_radius.upper + threshold;
}

#endif
