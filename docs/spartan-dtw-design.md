# SPARTAN DTW Design

## Goal

Add exact Dynamic Time Warping (DTW) search to a live in-memory SPARTAN index
without rebuilding its PCA-derived symbolic representation.  The index provides
a cheap, admissible partition-level bound; raw-domain `LB_Keogh` and banded DTW
then refine surviving records.

## Query Pipeline

For a query `Q` and Sakoe-Chiba warping radius `r`:

1. Apply the same normalization used when building the index.
2. Construct the raw DTW envelope `B_Q = [L(Q), U(Q)]` with radius `r`.
3. Convert `B_Q` to a PCA-coordinate rectangle `E_Q` using the stored SPARTAN
   PCA components and bias.
4. Traverse symbolic cells using the squared distance from `E_Q` to each cell's
   PCA interval as the partition-level lower bound.
5. For surviving records, apply raw `LB_Keogh` against `L(Q), U(Q)`.
6. Evaluate banded exact DTW only for records that survive both bounds.

All bound and DTW values are squared distances, matching the existing DTW
implementation.

## Safe PCA Envelope Projection

SPARTAN stores PCA coordinates as:

```text
z_j = b_j + sum_i(p_ji * x_i)
```

where `p_ji` is a PCA loading and `b_j` is the existing PCA bias.  Do not use
`P^T L` and `P^T U` as componentwise endpoints: PCA loadings can be negative,
so componentwise order is not preserved.

For every retained component `j`, construct the interval with sign-aware
arithmetic:

```text
e_min[j] = b[j] + sum(p_ji >= 0 ? p_ji * L[i] : p_ji * U[i])
e_max[j] = b[j] + sum(p_ji >= 0 ? p_ji * U[i] : p_ji * L[i])
```

The resulting axis-aligned box contains the projection of every raw point in
the query envelope.  For a symbolic cell interval `C_j = [c_min[j], c_max[j]]`,
use:

```text
gap(C_j, E_j) = e_min[j] - c_max[j]  if c_max[j] < e_min[j]
                 c_min[j] - e_max[j]  if e_max[j] < c_min[j]
                 0                    otherwise

LB_spartan_dtw = sum_j(gap(C_j, E_j)^2)
```

With an orthonormal-row projection `P`, this is admissible:

```text
LB_spartan_dtw
  <= dist(P^T X, P^T B_Q)
  <= dist(X, B_Q)
  =  LB_Keogh(X, Q)
  <= DTW_r(X, Q)
```

The bound is intentionally weaker than raw `LB_Keogh`; it is useful because it
can prune whole symbolic partitions before individual raw records are read.

## Repository Findings

- SPARTAN stores `pca_components` and `pca_bias`; it applies no post-PCA
  coordinate scaling.  Its existing bins are monotone intervals in the same
  PCA coordinate system.
- The PCA basis is computed in double precision but stored as `float`.  A
  future exact implementation must widen envelope endpoints outward for
  rounding and either retain a certified orthonormal basis or divide the
  projected bound by a conservative upper bound on `||P||_2`.
- A live in-memory SPARTAN trie or iSAX index can be reused: the change is to
  the query rectangle and its distance-to-cell calculation, not to indexed
  records or bins.
- `src/ads/DTWfunction.c` contains legacy iSAX-oriented envelope, `LB_Keogh`,
  and banded DTW code, but it is not compiled, linked, or routed from the
  current CLI.  Treat it as reference code rather than an active feature.

## Implementation Work

1. Add an explicit DTW query mode and validated warping-window option.
2. Add a SPARTAN helper that builds the sign-aware PCA envelope rectangle.
3. Add rectangle-to-symbolic-cell lower bounds for trie and/or iSAX traversal;
   do not use the current point-query PCA lower bound for DTW pruning.
4. Reuse or modernize the raw `LB_Keogh` and banded DTW routines behind the new
   query mode.
5. Preserve exactness when allocation or optional optimization fails by falling
   back to a complete candidate scan rather than dropping candidates.

## Validation

- Compare every DTW result with brute-force banded DTW on small float and
  integer-input datasets for several radii, including `r = 0`.
- Verify no false dismissals by checking that every node and record lower bound
  is no greater than its exact DTW distance.
- Test normalized and non-normalized inputs, nontrivial PCA bias, mixed-sign
  PCA components, and boundaries that coincide with symbolic-bin edges.
- Benchmark candidate counts and elapsed time against raw `LB_Keogh` plus a
  full scan, separately for trie and iSAX once each is implemented.

## Open Product Decisions

- Initial layout scope: trie only, iSAX only, or both.
- Initial result scope: exact 1-NN only or exact top-k as well.
- Whether the first release uses only symbolic PCA-cell bounds or additionally
  stores raw envelopes for node-level pruning.
