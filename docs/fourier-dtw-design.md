# Fourier Envelope DTW Design

## Goal

Use an SFA/Fourier symbolic index to prune exact constrained DTW search without
rebuilding a live in-memory index.  The index-level bound transforms the raw
DTW query envelope into a safe Fourier-coordinate rectangle; raw `LB_Keogh`
and banded DTW refine the remaining records.

## Query Pipeline

For a query `Q` and Sakoe-Chiba radius `r`:

1. Apply the same normalization used during index construction.
2. Build the raw envelope `B_Q = [L(Q), U(Q)]` for `r`.
3. Map that box to a sign-aware rectangle over the exact normalized Fourier
   coordinates stored by the index.
4. Prune symbolic cells by squared interval-gap distance to that rectangle.
5. Apply raw `LB_Keogh` to records in surviving cells.
6. Run banded exact DTW only for records that survive both lower bounds.

All values use squared-distance units, consistent with the current DTW code.

## Safe Fourier Envelope Projection

Let `T` be the normalized Fourier coordinate map used by the index.  For a
unitary DFT, or any subset of its coordinates, `T` is non-expansive:

```text
||T(x - y)||_2 <= ||x - y||_2
```

For every `y` in the raw envelope `B_Q`, this gives:

```text
dist(TX, T(B_Q)) <= dist(X, B_Q) = LB_Keogh(X, Q)
```

Construct an axis-aligned Fourier rectangle `H_Q` that contains `T(B_Q)`.  As
`H_Q` is a superset of the exact transformed set:

```text
dist(TX, H_Q) <= LB_Keogh(X, Q) <= DTW_r(X, Q)
```

This is not ordinary Fourier indexing of `X` and `Q`.  The DTW guarantee is
introduced first by the raw envelope; Fourier is only a non-expansive linear
representation of that envelope-derived region.

## Re-verified Contraction Argument

The Fourier claim can be stated independently of any particular symbolic
layout.  Keogh and Ratanamahatana establish that the envelope distance is a
lower bound for constrained DTW:

```text
LB_Keogh(X, Q) = dist(X, B_Q) <= DTW_r(X, Q)
```

Let `F` be the unitarily normalized DFT and `F_S` retain only the coordinate
rows used by an index.  Parseval gives `||F(x - y)||_2 = ||x - y||_2`; dropping
coordinates makes `F_S` a contraction:

```text
||F_S(x - y)||_2 <= ||x - y||_2
```

Consequently:

```text
dist(F_S X, F_S(B_Q)) <= dist(X, B_Q) = LB_Keogh(X, Q)
```

The exact transformed image `F_S(B_Q)` is generally a zonotope, not a box.
The sign-aware Fourier rectangle `H_Q` deliberately encloses that zonotope.
Distance to this superset can only decrease, which yields the final admissible
bound:

```text
dist(F_S X, H_Q) <= LB_Keogh(X, Q) <= DTW_r(X, Q)
```

This is the basis for safe symbolic-cell pruning.  Classic Fourier indexing
establishes the corresponding reduced-Fourier Euclidean lower-bound principle;
the composition with a DTW envelope follows directly from the inequalities
above.

### Per-coordinate intervals

For a normalized complex Fourier coefficient:

```text
z_k = sum_j(x_j * exp(-2 pi i k j / n)) / sqrt(n)
```

build independent intervals for its real and imaginary coordinates.  For a
linear coordinate `v = sum_j(a_j * x_j)`, use sign-aware interval arithmetic:

```text
v_min = sum(a_j >= 0 ? a_j * L[j] : a_j * U[j])
v_max = sum(a_j >= 0 ? a_j * U[j] : a_j * L[j])
```

Use `a_j = cos(2 pi k j / n) / sqrt(n)` for the real coordinate and
`a_j = -sin(2 pi k j / n) / sqrt(n)` for the imaginary coordinate.  Apply the
same formula to any selected or reordered SFA coefficient coordinate.

For a symbolic cell interval `C_d = [c_min[d], c_max[d]]` and query interval
`H_d = [h_min[d], h_max[d]]`:

```text
gap(C_d, H_d) = h_min[d] - c_max[d]  if c_max[d] < h_min[d]
                 c_min[d] - h_max[d]  if h_max[d] < c_min[d]
                 0                     otherwise

LB_fft_envelope = sum_d(gap(C_d, H_d)^2)
```

If this lower bound exceeds the current best DTW distance, the whole symbolic
cell can be safely discarded.

## Existing Fourier Representations

The repository already uses the required `1 / sqrt(n)` base normalization:

- `fft_from_ts` stores selected real and negative-imaginary coordinates, each
  scaled by `index->norm_factor = 1 / sqrt(n)`.  This is a subset of unitary
  DFT real coordinates, so it is non-expansive even though it omits conjugate
  partners and is therefore intentionally weaker.
- `fft_full_real_from_ts` uses the same edge scaling and applies `sqrt(2)` to
  non-DC/non-Nyquist real/imaginary pairs.  This preserves the full
  time-domain Euclidean norm through Parseval.

The DTW query bound must use the exact representation stored by the target
index.  In particular, SFA's selected coefficient order, DC omission for
normalized inputs, and any variance-selected coefficient list must be applied
to the envelope coordinates exactly as they are applied to records.

Do not apply an additional `1 / sqrt(n)` factor to the current SFA coordinates:
they are already normalized.  If a future index uses an unnormalized FFT,
divide its Fourier-space squared distance by `n` (or its Euclidean distance by
`sqrt(n)`) before comparing it to DTW.

For an `rfft` representation, selecting unique normalized coefficients without
compensating for omitted conjugates remains contractive and is therefore safe,
although weaker.  Apply the interior `sqrt(2)` factor only when the chosen
real-coordinate representation is explicitly designed to preserve Parseval
energy, as `fft_full_real_from_ts` is; do not add it to an otherwise unchanged
positive-frequency representation.

## Numerical Safety

- Evaluate sine/cosine and interval accumulations with enough precision, then
  round lower endpoints downward and upper endpoints upward before traversal.
- Keep the transform's actual coefficient scaling in the bound; arbitrary
  post-transform amplification invalidates the direct lower-bound comparison.
- For a non-unit operator `T`, divide the transformed squared distance by a
  conservative bound on `||T||_2^2`.
- Preserve exactness on allocation or optimization failure by falling back to
  a complete candidate scan rather than dropping candidates.

## Implementation Work

1. Add a DTW query mode and validated warping-window option.
2. Add a helper that constructs normalized Fourier-envelope intervals for the
   index's selected SFA coordinates.
3. Add rectangle-to-symbolic-cell bounds for SFA trie and/or iSAX traversal;
   do not reuse the existing point-to-cell Euclidean bound unchanged.
4. Reuse or modernize the existing raw `LB_Keogh` and banded DTW code as the
   record-level refinement stages.
5. Keep exact result semantics for 1-NN before extending the feature to top-k.

## Validation

- Compare results to brute-force banded DTW for several radii, including
  radius zero, on float and integer input files.
- Verify every Fourier-envelope node/record bound is no greater than the exact
  DTW distance for randomly generated envelopes and mixed-sign Fourier rows.
- Test DC, Nyquist, interior real/imaginary pairs, truncated coordinates,
  SFA coefficient selection, and normalized-input DC omission.
- Verify that the full real representation preserves Parseval energy and that
  the selected SFA representation remains non-expansive.
- Benchmark pruning and elapsed time against raw `LB_Keogh` plus a full scan.

## Literature Position

The individual results are established: the DTW envelope lower bound, Parseval
normalization, and reduced-Fourier lower bounds for Euclidean search.  This
document does not claim that the exact composition
`LB_Keogh envelope -> Fourier interval box -> symbolic-cell distance` has
already appeared as a published indexing construction.  Its validity here is
mathematical; publication precedent remains a separate literature-review
question.

## Open Product Decisions

- Initial layout scope: trie only, iSAX only, or both.
- Initial result scope: exact 1-NN only or exact top-k as well.
- Whether to first target SFA's selected-coordinate representation or a full
  Fourier representation with an additional symbolic index.
