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
  than a rank-`degree` domain. Candidates are filtered by requiring every
  embedding coordinate to be within a computed radius of the origin, then
  capped to the smallest-norm `candidateCap` (default 500) so absorption
  cost per box stays bounded regardless of degree.

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
  integers directly. This is not just an optimization: fields with a large
  regulator (e.g. real quadratic fields with a large fundamental unit)
  have their useful absorbers spread out along the unit orbit rather than
  clustered near the origin, and the sieve cannot converge in practice
  without it (`x^2-61`, with fundamental unit ~39, is the empirical
  example that demonstrated this during development).

- **`Certify.chpl`** — Phase 3. `runSieve` only ever proves upper bounds
  (`M(K) < K`); a failure to clear does not prove a matching lower bound,
  since it can also just mean "needs more resolution". But every surviving
  box has a rational (dyadic) center, and Pari can compute
  `min_gamma |N(x - gamma)|` *exactly* at that specific point — since
  trivially `m_K(x) <= M(K_field)`, an exact value there is a genuine,
  proven lower bound. `simplestFractionInInterval` additionally suggests a
  human-readable closed form (e.g. `1/3`) from a tight numeric bracket, as
  a conjecture to sanity-check against (not a proof by itself).

- **`main.chpl`** — CLI entry point: loads a field, brackets M(K) with an
  exponential search, refines the bracket via bisection (each trial an
  independent, parallel `runSieve` call), and finally attempts Phase 3
  certification on the most-resistant boxes from the last "blocked" trial.

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
- `--useUnits` / `--unitExponentRange` — control the Phase 2 unit-action
  acceleration (on by default; see above for why it matters).
- `--certify` / `--certifySamples` — control Phase 3 (on by default;
  requires `gp` on `PATH`).

## Validation

`tests/validate.sh` runs the reference C `euclid` binary (ground truth)
and the Chapel pipeline side by side over `tests/fixtures/`. Results as of
this implementation:

| Field | Reference exact minimum | Chapel bracket |
|---|---|---|
| `x^2-2` | 1/2 | `[0.5, 0.500977]` |
| `x^2-61` | 1611/1525 (~1.05639) | `[1.05762, 1.05859]` |
| `x^3+x^2-1` | 1/5 | `[0.199219, 0.200195]` |
| `x^3-3*x-1` | 1/3 | `[0.332897 (Pari-certified), 0.333398]` |
| `x^5-x-1` | 1/4 | converges but looser; demonstrates degree-5 support |

`x^2-61` and `x^3-3*x-1` are specifically the cases that the earlier
prototypes could not solve correctly (large regulator, and rank-2 unit
group respectively).

Other test files under `tests/`:

- `test_numberfield.chpl`, `test_smallelements.chpl`, `test_sieve.chpl`,
  `test_sieve_debug.chpl`, `test_certify.chpl` — standalone smoke tests for
  each module, useful when iterating on a single piece in isolation.

## Known limitations / future work

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
  port remains a documented fallback if this proves insufficient.
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
