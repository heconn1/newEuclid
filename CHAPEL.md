# Chapel port: arbitrary-degree Euclidean-minimum sieve

This document describes the Chapel implementation that lives alongside the
original Pari/C `euclid` tool (`src/`, built via `make`, see `README` /
`INSTALL`). It computes the Euclidean minimum M(K) of a number field K using
Pierre Lezowski's algorithm as the mathematical reference, but keeps all of
the numerically-heavy work in floating point / Chapel. Pari/gp is used only
for:

- a one-time, per-field setup step (`generate_field.gp`, producing
  `field_data.txt`: the field's degree, signature, integral-basis
  embeddings, and fundamental-unit embeddings), and
- an optional, final exact-certification step (`Certify.chpl`) that shells
  out to `gp` once per candidate point.

## Why this is a rewrite, not a patch

The earlier Chapel prototypes (`firstCut.chpl` ... `sixthCut.chpl`,
`kernel.chpl`, `kernel.py`, `exactMin.chpl`, now moved to `attic/`) did not
implement Lezowski's algorithm. They tested whether a point's norm (or a
unit-scaled norm — which is numerically almost the same test, since every
unit has `|N(unit)| = 1`) was below a target bound, but never tested
*absorption* against real candidate algebraic integers (`N(x - gamma) <
bound` for enumerated `gamma`). That's the actual mathematical core of the
method, and it was entirely missing, which is why those versions frequently
reported wrong minima and were hard-capped at low degree by ad hoc,
hand-unrolled loops.

## Architecture

- **`NumberField.chpl`** — parses `field_data.txt` into a `NumberFieldData`
  record (degree, r1, r2, the "complex-compact" embedding matrix: r1 real
  + r2 complex embedding coordinates per basis element, and the
  fundamental-unit embeddings). Also computes the inverse embedding matrix
  (`sigmaInv`) at load time via a small dependency-free Gauss-Jordan solver
  (Chapel's `LinearAlgebra.inv()` requires a system LAPACK that may not be
  installed, so this avoids that dependency). No changes to
  `generate_field.gp` or the `field_data.txt` format were needed.

- **`SmallElements.chpl`** — enumerates candidate algebraic integers
  ("absorbers") for a given target bound K: the arbitrary-degree analogue
  of Lezowski's `small_elts` (`src/main.c`). Chapel domains must have a
  compile-time-constant rank, so the enumeration uses a single flattened
  1-D `forall` over a mixed-radix index space (an odometer decode) rather
  than a rank-`degree` domain. Every enumerated candidate is ranked by its
  actual field norm and the smallest-norm `candidateCap` (default 2000)
  are kept. **This was originally a per-coordinate ball filter ("keep w iff
  every embedding coordinate is within a computed radius of the origin")
  instead of a norm ranking, which is wrong**: it systematically excludes
  continued-fraction-convergent-like elements whose individual coefficients
  are large but whose norm is small, which are exactly the absorbers real
  quadratic (and other unit-rich) fields need. This was caught by testing
  against a battery of norm-Euclidean real quadratic fields (`q11.txt` ...
  `q73.txt`), several of which converged to confidently wrong answers under
  the old filter (e.g. `x^2-57` converged to exactly 2x the true minimum).

- **`Sieve.chpl`** — the core decision procedure: "is K a valid upper bound
  for M(K_field)?" A `Box` is an axis-aligned hyper-cube in integral-basis
  coefficient space. `computeBoxProjections` maps it to an
  interval/disk bound in embedding space; `calculateBoxMaxNorm` bounds the
  worst-case norm over that box. A box is *absorbed* once some candidate's
  shifted bound is provably below K (the generalization of Lezowski's
  `norme_bricolee`). Unabsorbed boxes are bisected into `2^degree` children
  and retried at finer resolution (`runSieve`).

  This module also implements the **unit-action acceleration**
  (`isBoxAbsorbedWithUnits`, `recenterByUnit`, `unitPower`): for any unit
  `u` and integer `gamma`, `|N(u*x - gamma)| = |N(x - gamma*u^-1)|`, so
  proving a *twisted-and-recentered* copy of a box is absorbed by a small
  candidate proves the *original* box is absorbed by the (potentially very
  large) integer `gamma*u^-1`, without ever enumerating such large
  integers directly. This helps some fields converge faster, but (after
  the `SmallElements.chpl` norm-ranking fix above) is no longer load-bearing
  for correctness the way it first appeared to be -- most of the real
  quadratic test fields now converge correctly even with `unitExponentRange=1`.

  **Floating-point safety margin**: the absorption comparison
  (`isProjectionAbsorbed`) checks `maxNorm < K * (1 - 1e-9)`, not
  `maxNorm < K`. This matters whenever a field's true critical point has
  dyadic coordinates (common), because then some box's edge sits *exactly*
  on the boundary at every bisection depth; as that box shrinks, its
  computed max-norm approaches K from above, and past a certain depth
  (empirically around width 1e-8) floating-point rounding in the
  embedding/subtraction/multiplication chain can make it round to just
  *below* K, causing a false "absorbed" verdict. `runSieve`'s default
  `minWidth = 1e-7` stops refinement before reaching that regime (going
  deeper doesn't produce a more correct answer for a boundary-exact case --
  the true gap between max-norm and K shrinks to exactly 0 there, so no
  amount of floating-point precision can safely resolve it, only a
  certified exact check via `Certify.chpl` can). This was caught by
  `x^2-2` (whose exact minimum 1/2 sits exactly at the coefficient-space
  boundary) briefly, incorrectly, reporting `M(K) < 0.5` once other fixes
  allowed deeper bisection.

- **`Certify.chpl`** — Phase 3. `runSieve` only ever proves upper bounds
  (`M(K) < K`); a failure to clear does not prove a matching lower bound,
  since it can also just mean "needs more resolution". But every surviving
  box has a rational (dyadic) center, and Pari can compute
  `min_gamma |N(x - gamma)|` *exactly* at that specific point — since
  trivially `m_K(x) <= M(K_field)`, an exact value there is a genuine,
  proven lower bound. `simplestFractionInInterval` additionally suggests a
  human-readable closed form (e.g. `1/3`) from a tight numeric bracket, as
  a conjecture to sanity-check against (not a proof by itself).

  Two complementary exact-search strategies are available, both
  independently guarded by the same soundness check in `main.chpl`
  (discard any sampled exact value `>= hi`, the proven numeric upper
  bound — such a value is provably impossible since `m_K(x) <= M(K) <
  hi`, so it must mean the search below missed the true nearest lattice
  point, regardless of which strategy or why):

  - `exactMinimalNormAtWithCandidates` (used by `main.chpl`) checks the
    sampled point against the *same* norm-ranked candidate list the sieve
    itself used, rather than a blind coefficient range -- since that list
    is exactly what made the sieve's own bracket converge correctly, this
    guarantees certification can never be *less* complete than the sieve
    was for that same K.
  - `exactMinimalNormAt` instead does an independent, self-contained
    search directly in the generated `gp` script: a `defaultSearchRange`
    scaled with degree (the same budget-based approach as
    `SmallElements.candidateBoundRange`) instead of a flat constant (a
    flat range like `3` is too small for some fields, e.g. `x^2-61` needs
    `searchRange >= 5`), combined with small powers of the field's
    fundamental units alongside `gamma` (since `|N(u*x - gamma)| =
    |N(x - gamma*u^-1)|` for any unit `u`) -- the exact-arithmetic
    analogue of `Sieve.chpl`'s unit-action acceleration. Useful when no
    Chapel-side candidate set is at hand.

- **`main.chpl`** — CLI entry point: loads a field, brackets M(K) with an
  exponential search, refines the bracket via bisection (each trial an
  independent, parallel `runSieve` call), and finally attempts Phase 3
  certification on the most-resistant boxes from a dedicated deeper sieve
  pass at the final lower bound.

  **Bracket-finding soundness fix**: a shallow (`exploreDepth`) sieve trial
  that reports "cleared" is always trustworthy (absorption proofs are sound
  regardless of the depth used to find them), but a shallow trial that
  reports "blocked" is *not* -- it can simply mean the depth budget was too
  small, not that K is genuinely too small. An earlier version used shallow
  "blocked" results to raise the working lower bound `lo`, which could lock
  it above the true minimum before the deep bisection phase ever ran,
  permanently excluding the true answer from the search range (this is how
  several real quadratic fields, e.g. `x^2-19`, previously converged
  confidently to the wrong value). `lo` now always starts at the trivially
  safe value 0 and is only ever raised from a *deep* (`refineDepth`)
  "blocked" result.

## Building and running

A `Makefile` target automates each step (all are additive to the existing
C-tool targets — `make` alone still just builds `euclid`, as before):

```sh
make chapel   # chpl -M . main.chpl -o euclid_chpl
make fixtures # generate any missing tests/fixtures/*.txt via gp
make smoke    # compile + run the standalone per-module smoke tests
make test     # build both tools and run tests/validate.sh (Chapel vs C)
make check    # smoke + test
make all      # euclid (C) + chapel
make clean-chapel  # remove Chapel build artifacts (also run by `make clean`)
```

Equivalently, by hand:

```sh
chpl -M . main.chpl -o euclid_chpl
./euclid_chpl --fieldFile=tests/fixtures/x3-3x-1.txt --tolerance=0.001
```

To generate `field_data.txt` for a new field (requires Pari/gp):

```sh
POLY="x^3-3*x-1" gp -q generate_field.gp
```

### Key CLI flags (see `main.chpl` for the full list)

- `--fieldFile` — path to a `field_data.txt`-format file (default
  `field_data.txt`).
- `--initialK` — starting upper-bound guess for the exponential bracket
  search.
- `--tolerance` — stop bisection once the bracket is this narrow.
- `--exploreDepth` / `--refineDepth` — bisection depth budget for the
  bracket-finding and refinement phases respectively; raise for
  higher-degree fields that need finer resolution near the true minimum.
- `--maxProblems` — safety valve: abandon a trial K early once the number
  of unresolved boxes exceeds this, rather than letting a doomed
  (K-too-small) trial explode combinatorially.
- `--candidateCap` — how many smallest-norm candidates `SmallElements`
  keeps (default 2000). Some real quadratic fields need several thousand
  to include the right convergent-like elements; raise this if a field
  still converges to a visibly-too-loose bracket.
- `--useUnits` / `--unitExponentRange` — control the Phase 2 unit-action
  acceleration (on by default; see above for why it matters).
- `--certify` / `--certifySamples` / `--certifyDepth` — control Phase 3
  (on by default; requires `gp` on `PATH`).

## Validation

`tests/validate.sh` runs the reference C `euclid` binary (ground truth)
and the Chapel pipeline side by side over `tests/fixtures/`. Results as of
this implementation:

| Field | Reference exact minimum | Chapel bracket |
|---|---|---|
| `x^2-2` | 1/2 | `[0.5, 0.500977]` |
| `x^2-61` | 1611/1525 (~1.05639) | `[1.05762, 1.05859]` |
| `x^3+x^2-1` | 1/5 | `[0.199219, 0.200195]` |
| `x^3-3*x-1` | 1/3 | `[0.333008 (Pari-certified ~0.333319), 0.333984]` |
| `x^5-x-1` | 1/4 | converges but looser; demonstrates degree-5 support |

`x^2-61` and `x^3-3*x-1` are specifically the cases that the earlier
prototypes could not solve correctly (large regulator, and rank-2 unit
group respectively).

Separately, `q11.txt` ... `q73.txt` (the 15 norm-Euclidean real quadratic
fields, per Chatland-Davenport) were used as a wider correctness battery,
since the C tool computes their exact minima quickly and several exposed
real bugs (see "Bugs found and fixed" below):

| Field | Reference exact minimum | Chapel bracket |
|---|---|---|
| `x^2-11` | 19/22 (~0.86364) | `[0.863281, 0.864258]`, fraction 19/22 exact |
| `x^2-13` | 1/3 | `[0.333008, 0.333984]`, fraction 1/3 exact |
| `x^2-19` | 170/171 (~0.99415) | `[0.994141, 0.995117]`, fraction 170/171 exact |
| `x^2-29` | 4/5 | `[0.799805, 0.800781]`, fraction 4/5 exact |
| `x^2-33` | 29/44 (~0.65909) | `[~0.6591, 0.65918]` |
| `x^2-41` | 23/32 (0.71875) | `[0.71875, 0.719727]`, certified 0.71875 |
| `x^2-57` | 14/19 (~0.73684) | `[0.736328, 0.737305]`, fraction 14/19 exact |
| `x^2-73` | 1541/2136 (~0.72144) | `[0.723633, 0.724609]` -- still ~0.3% high; largest regulator (2136) in the battery, needs more `--candidateCap`/depth than the current defaults budget for in reasonable time |

### Bugs found and fixed during this battery

1. **`SmallElements.chpl`'s candidate filter excluded exactly the
   absorbers real quadratic fields need** (per-coordinate ball instead of
   norm ranking -- see the `SmallElements.chpl` bullet above). This alone
   caused confidently-wrong answers (e.g. `x^2-57` computed as exactly 2x
   the true minimum, `x^2-19` and others off by 10-90%).
2. **Bracket-finding could lock the lower bound above the true minimum**
   using untrustworthy shallow-depth "blocked" results (see the
   `main.chpl` bullet above).
3. **A floating-point precision cliff at exact (dyadic) critical points**
   could produce a false "cleared" result if bisection was allowed to go
   too deep (see the `Sieve.chpl` unit-action bullet above for the
   `minWidth`/safety-margin fix).
4. **`Certify.chpl`'s exact check used the same too-small a search as the
   original (buggy) `SmallElements.chpl` filter**, and independently
   needed the norm-ranked candidate list to avoid reporting inflated (and
   sometimes upper-bound-contradicting) "certified" values.

All four were required together to get the real quadratic battery
converging correctly; `x^2-73` shows the remaining known limit (very large
regulators need more compute than the current defaults budget for).

Other test files under `tests/`:

- `test_numberfield.chpl`, `test_smallelements.chpl`, `test_sieve.chpl`,
  `test_sieve_debug.chpl`, `test_certify.chpl` — standalone smoke tests for
  each module, useful when iterating on a single piece in isolation.

## Known limitations / future work

- **Very large regulators need more resources than the current defaults**:
  `x^2-73` (fundamental unit ~2136) converges to a bracket that's correct
  in direction but ~0.3% too high with default settings; larger
  `--candidateCap`/`--refineDepth` help but cost proportionally more time
  per trial. A field-size-aware auto-tuning heuristic (rather than fixed
  defaults) is the natural fix, mirroring the hand-tuned, degree-indexed
  `euclid.cfg` tables in the original C tool but keyed off the regulator
  too, not just the degree.
- **Performance at higher degree**: the sieve's `2^degree` branching factor
  and per-box absorption cost mean degree-5+ fields converge more slowly
  and less tightly than degree 2-3 within the same time budget (see the
  `x^5-x-1` result above). Tuning `candidateCap`, `unitExponentRange`, and
  the depth/tolerance parameters per field is currently manual, mirroring
  the hand-tuned `euclid.cfg` in the original C tool.
- **Phase 3 is not the full Lezowski cycle-decomposition machinery**: it
  samples the most-resistant boxes rather than exactly identifying the
  unit-orbit critical cycle (`src/graph.c`'s Tarjan-based decomposition).
  It works well for the well-known small examples but is not guaranteed to
  land exactly on the supremum for harder fields; the full graph/cycle
  port remains a documented fallback if this proves insufficient. Thanks
  to the degree-scaled search range and the `hi`-based soundness guard
  (see `Certify.chpl` above), a reported certified bound is always sound
  (never exceeds the true minimum) even when it isn't tight.
- **Parallelization roadmap** (not yet built): the per-level box list and
  the small-elements enumeration are both `forall`-parallel today
  (single-locale, multi-core). Scaling to multiple locales (distribute the
  "problems" list, e.g. via a `Block`-distributed array) and to GPUs
  (`computeBoxProjections`/`calculateBoxMaxNorm` are flat numeric kernels
  well-suited to `foreach`/GPU offload) is the natural next step once
  higher-degree performance needs it.

## Relationship to the original C tool

`src/*.c` and the `euclid` binary (built via `make`) are unchanged and
remain the ground-truth oracle used by `tests/validate.sh`. The superseded
Chapel prototypes live in `attic/` for historical reference.
