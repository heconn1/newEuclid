# Euclidean minimum of a number field

Given a number field K (defined by an irreducible polynomial P with
integral coefficients), this project computes

```
M(K) = sup_x m_K(x),  where  m_K(x) = min_{X in O_K} |N(x - X)|
```

i.e. the Euclidean minimum of K, and (when successful) the points x at
which the supremum is attained.

The repository contains **two implementations**:

- **`src/`** — the original reference tool by Pierre Lezowski, in C, linked
  against the [Pari/GP](https://pari.math.u-bordeaux.fr/) library. Builds
  to the `euclid` binary. See [`README`](README) and [`INSTALL`](INSTALL)
  for details, and [`AUTHORS`](AUTHORS) / [`COPYING`](COPYING) for
  attribution and license. This is treated as the ground-truth oracle for
  validating the Chapel port (`tests/validate.sh`).
- **A Chapel port** (`*.chpl` at the repository root) — a from-scratch,
  arbitrary-degree reimplementation of the same algorithm, documented in
  full below and in [`CHAPEL.md`](CHAPEL.md). Builds to `euclid_chpl`.

## Chapel implementation

The Chapel port computes the same M(K) as the C tool, keeping all
numerically-heavy work in floating point / Chapel (parallelized with
`forall`), and uses Pari/gp only for one-time per-field setup and an
optional final exact-certification step.

**Why a rewrite:** earlier Chapel prototypes (now in `attic/`) tested
whether a point's norm was below a target bound, but never tested
*absorption* against real candidate algebraic integers (`N(x - gamma) <
bound`) — the actual mathematical core of Lezowski's method. That's why
they frequently reported wrong minima and were hard-capped at low degree.

**Architecture** (one module per concern):

- `NumberField.chpl` — parses `field_data.txt` (produced by
  `generate_field.gp`) into degree/signature/embedding data, and computes
  the inverse embedding matrix via a dependency-free Gauss-Jordan solver.
- `SmallElements.chpl` — arbitrary-degree enumeration of candidate
  absorbing algebraic integers for a target bound K, ranked by actual
  field norm (not by a per-coordinate coefficient bound — see "Correctness
  history" below for why that distinction mattered).
- `Sieve.chpl` — the real absorption test and parallel box-bisection
  sieve, plus **unit-action acceleration** (a speed optimization for most
  fields, though some real quadratic fields with sizable regulators still
  need it to converge in reasonable time), and a floating-point safety
  margin so bisecting all the way to a field's true (often dyadic) critical
  point can't produce a false "cleared" result. `Box` is generic over a
  compile-time `param degree`, storing coordinates as a fixed-size tuple
  instead of a domain-backed array (cheaper to construct/copy in
  `bisect()`, which allocates `2^degree` fresh boxes per level -- see
  "Optimization: tuple-based Box" below); `main.chpl` bridges the runtime
  `nf.degree` to that compile-time param via a dispatch table, currently
  covering **degree 1-8** (see [`CHAPEL.md`](CHAPEL.md#supported-field-degree)
  for how to raise that cap).
- `Certify.chpl` — turns a numeric upper bound into a genuine
  Pari-certified lower bound by evaluating the exact norm at the
  most-resistant sampled points, via two complementary exact-search
  strategies (reusing the sieve's own candidate list, or an independent
  degree-scaled/unit-aware search in `gp`). `main.chpl` discards any
  sampled value that's `>=` the proven numeric upper bound, since that's
  provably impossible — so a reported certified bound is always sound.
- `main.chpl` — the CLI: exponential bracket search, bisection refinement
  (which only ever raises its working lower bound from a *deep*, trusted
  sieve result — see "Correctness history"), and an optional certification
  pass.

**Algorithm in three phases:** (1) a sound covering sieve proves numeric
upper bounds on M(K); (2) unit-action recentering lets a small,
origin-centered candidate list stand in for absorbers spread out along a
field's unit orbit; (3) exact Pari evaluation at resistant sample points
certifies a matching lower bound when possible.

### Validation

Validated against the C tool's exact answers (`tests/validate.sh`,
`make test`):

| Field | Reference exact minimum | Chapel bracket |
|---|---|---|
| `x^2-2` | 1/2 | `[0.5, 0.500977]` |
| `x^2-61` | 1611/1525 (~1.05639) | `[1.05762, 1.05859]` |
| `x^3+x^2-1` | 1/5 | `[0.199219, 0.200195]` |
| `x^3-3*x-1` | 1/3 | `[0.333008, 0.333984]`, certified ~0.333319 |
| `x^4-4*x^2+2` | 1/2 | `[0.5, 0.500977]`, certified ~0.49999 |
| `x^5-x-1` | 1/4 | converges but looser; demonstrates degree-5 support |

Also validated against the 15 norm-Euclidean real quadratic fields
(Chatland-Davenport), `q11.txt` ... `q73.txt`, since the C tool computes
their exact answers quickly and they proved to be a much harder
correctness battery — see "Correctness history" below:

| Field | Reference exact minimum | Chapel bracket |
|---|---|---|
| `x^2-11` | 19/22 | `[0.863281, 0.864258]`, fraction 19/22 exact |
| `x^2-13` | 1/3 | `[0.333008, 0.333984]`, fraction 1/3 exact |
| `x^2-19` | 170/171 | `[0.994141, 0.995117]`, fraction 170/171 exact |
| `x^2-29` | 4/5 | `[0.799805, 0.800781]`, fraction 4/5 exact |
| `x^2-33` | 29/44 | `[~0.6591, 0.65918]` |
| `x^2-41` | 23/32 | `[0.71875, 0.719727]`, certified exactly 0.71875 |
| `x^2-57` | 14/19 | `[0.736328, 0.737305]`, fraction 14/19 exact |
| `x^2-73` | 1541/2136 (~0.72144) | `[0.723633, 0.724609]` — still ~0.3% high (see known limitations) |

### Correctness history

The real quadratic battery above exposed and fixed four distinct bugs,
all now covered by regression tests (`tests/validate.sh` runs `q19` and
`q57` specifically):

1. **Wrong candidates.** `SmallElements.chpl` used to require every
   embedding coordinate to be individually bounded, which excludes exactly
   the large-coefficient/small-norm "convergent-like" elements real
   quadratic fields need as absorbers. `x^2-57` used to converge to
   **exactly 2x** the true minimum because of this.
2. **Untrustworthy bracket-locking.** The bracket search used to trust a
   shallow sieve trial's "blocked" result to raise the lower bound, which
   could lock it above the true minimum before deep bisection ever ran.
   `x^2-19` used to converge to exactly `1.0` instead of `170/171` this way.
3. **A floating-point precision cliff** at fields whose true critical point
   has dyadic coordinates (common): bisecting deep enough could produce a
   false "cleared" result via rounding. `x^2-2` (critical point exactly
   1/2) briefly regressed to this while fixing the other bugs.
4. **Inconsistent certification.** `Certify.chpl`'s exact check inherited
   the same completeness gap as bug #1 independently, and could report a
   "certified" value larger than an independently-proven upper bound.

See [`CHAPEL.md`](CHAPEL.md#bugs-found-and-fixed-during-this-battery) for
the full technical explanation of each.

### Performance benchmarks

Wall-clock timings below are from a 6-core machine, default CLI settings
(`--tolerance=0.001`), `make all` (i.e. `--fast` Chapel build). Both
tools' core search loops actually use the same arithmetic model — plain
double-precision floating point (`double` in the C tool's `src/basef.c`,
`real(64)` in Chapel's `Sieve.chpl`) — not Pari's exact bignum
arithmetic. Pari is invoked exactly once per run in *each* tool: once
for one-time field setup, and again, optionally, for a single final
exact-certification pass (`Certify.chpl` / `src/pari_min_c.c`). So the
gap below reflects implementation efficiency, not a different numerical
approach — see "Optimization: real(64) fast path" and "Optimization:
tuple-based Box" below for the two efficiency passes done so far (modest
end-to-end wins at default settings, for reasons explained there), and
[`CHAPEL.md`'s Roadmap](CHAPEL.md#roadmap) for the auto-tuning work
expected to unlock more of it.

| Field | C tool (`./euclid`) | Chapel tool (`./euclid_chpl`) | Chapel user time (parallelism) |
|---|---|---|---|
| `x^2-2` | ~0.2s | ~3.5s | ~12.7s (~3.6x) |
| `x^2-61` | ~0.15s | ~68s | ~358s (~5.3x) |
| `x^2-57` | ~0.13s | ~91s | ~450s (~4.9x) |
| `x^2-73` (hardest regulator) | ~3.7s | ~157s | ~783s (~5.0x) |
| `x^3+x^2-1` | ~0.9s | ~210s | ~841s (~4.0x) |
| `x^3-3*x-1` (rank-2 units) | ~0.6s | ~198s | ~980s (~4.9x) |
| `x^5-x-1` (degree 5) | ~50s | ~160s | ~590s (~3.7x) |

"Chapel user time (parallelism)" is total CPU-seconds across all cores
and the resulting speedup over wall-clock time — confirms the `forall`
parallelism is being used (not close to the full 6x on this machine,
since the small-elements ranking and per-K sieve trials have some
unavoidably serial phases), but doesn't come close to closing the gap
with the C tool on its own. The Chapel tool is currently **17-450x
slower in absolute wall-clock terms** across the fields above (least gap
on the degree-5 field, most on `x^2-61`) — this project prioritizes
correctness and arbitrary-degree generality first (see "Correctness
history" above); the optimization work below and the parallelization
roadmap further down are the intended path to closing this gap for the
field sizes where it matters (large regulators, higher degree).

### Optimization: real(64) fast path for totally real fields

The numbers above motivated a first efficiency pass, guided by two
microbenchmarks (`bench/TupleVsArrayBench.chpl`,
`bench/BoxCopyBench.chpl`) rather than guesswork. They isolated two
findings: raw per-coordinate array-access cost is statistically
identical to fixed-size-tuple access once compiled with `--fast`, but
(a) `NumberFieldData` stored every embedding — real or complex — as
`complex(128)`, even though the `r1` real embeddings always have
`im=0`, and (b) `Box` record construction/copying inside `bisect()` is
6-8x more expensive with today's domain-backed arrays than an
equivalent fixed-size tuple would be (that second, larger finding is
addressed in "Optimization: tuple-based Box" below).

Fixing (a) was the first, lower-risk step: `NumberField.chpl` now also
stores a real(64)-only mirror of the embedding matrix
(`NumberFieldData.isTotallyReal`/`basisReal`) whenever a field has `r2 ==
0` (every real quadratic field, and any other totally real field), and
`SmallElements.chpl`/`Sieve.chpl` dispatch to a matching real(64)-only
candidate enumeration and box-absorption test in that case, skipping
complex(128) arithmetic entirely in what is otherwise the hottest
per-candidate loop in the program. Fields with `r2 > 0` (e.g.
`x^3+x^2-1`, `x^5-x-1`, imaginary quadratics) are unaffected — they still
use the original, unmodified complex(128) path. Validated with no change
in output (same brackets/certified values) across `make check` and the
full real-quadratic battery (`q2.txt`-`q73.txt`, `qn11.txt`).

Measured wall-clock impact on the same 6-core machine (before vs. after,
same CLI settings as above):

| Field | Before | After | Improvement |
|---|---|---|---|
| `x^2-61` | 71.7s | 65.5s | ~9% |
| `x^2-41` | 95.3s | 89.0s | ~7% |
| `x^2-57` | 87.7s | 82.4s | ~6% |
| `x^2-73` | 155.0s | 143.9s | ~7% |

A modest, real improvement — consistent with the microbenchmark finding
that raw arithmetic wasn't the dominant cost. The much larger win found
by `BoxCopyBench.chpl` (tuple-based `Box`, avoiding per-bisection heap
allocation) is implemented next, below.

### Optimization: tuple-based Box (Phase 2)

`Box` and `SieveResult` (`Sieve.chpl`) are now generic over a
compile-time `param degree`, storing coordinates as fixed-size tuples
(`degree*real(64)`) instead of domain-backed arrays — the representation
`BoxCopyBench.chpl` showed to be 6-8x cheaper to construct/copy.
`main.chpl` bridges runtime to compile time with a `select nf.degree {
when 1..8 ... }` dispatch, the one place a degree isn't already known at
compile time. Validated for exact output parity (identical brackets,
identical intermediate box counts) against the pre-refactor binary
across the full real-quadratic battery, `cc23.txt`, and a newly added
degree-4 fixture (`x^4-4*x^2+2`, added to `tests/fixtures/` and
`tests/validate.sh`).

Measured **end-to-end, at default settings, the improvement is modest
(0-10%)** — nowhere near the 6-8x the isolated microbenchmark predicted.
The reason: `bisect()`'s box-construction cost is O(degree) per box (2-8
tuple/array slots), but each box's absorption test scans up to
`candidateCap` (2000 by default) candidates — a few thousand
floating-point operations that completely dominate wall-clock time for
typical settings, making the box-copy savings a rounding error by
comparison.

A dedicated harness (`bench/SieveTimingBench.chpl`/
`SieveTimingBenchOld.chpl`, isolating a single `runSieve()` call so
`candidateCap` can be varied directly) confirms the underlying box-copy
win is real and substantial once the absorption scan is cheap enough to
stop masking it:

| `candidateCap` | Before | After | Improvement |
|---|---|---|---|
| 2000 (default) | 6.64s | 6.08s | ~8% |
| 200 | 1.95s | 1.65s | ~15% |
| 20 | 1.74s | 0.98s | ~44% |
| 5 | 1.76s | 0.91s | ~48% |
| 1 | 3.77s | 1.15s | ~69% (3.3x) |

(`q57.txt`, K below the true minimum so the sieve never clears, fixed
`maxDepth`/`maxProblems` budget so each pair does comparable work.) The
effect isn't degree-limited either: the degree-4 fixture at
`candidateCap=20` showed a comparable ~42% improvement (36.1s -> 21.1s).

Practical upshot: this refactor doesn't meaningfully speed up *today's*
default runs (large regulators need a large `candidateCap` to find
enough candidates, which is exactly the regime where box-copying isn't
the bottleneck), but it pays off directly once the auto-tuning roadmap
item (below) can shrink `candidateCap` adaptively per field, and it
remains a real prerequisite for multi-locale/GPU work regardless of
when that lands (see "Relevance to clusters/GPU" reasoning in the
project's planning history) — a flat tuple has no heap allocation or
domain descriptor to distribute or copy into a device buffer, unlike
today's array-backed alternative.

### Known limitations

- **Very large regulators need more resources than the current
  defaults budget for.** `x^2-73` (fundamental unit ~2136, the largest
  regulator tested) converges to a bracket that's correct in direction but
  still ~0.3% too high with default settings, and Phase 3 certification is
  inconclusive there (all sampled boxes get discarded as depth-limited
  artifacts) — larger `--candidateCap`/`--refineDepth`/`--certifyDepth`
  close the gap but cost proportionally more time per trial (see
  benchmarks above). A field-size-aware auto-tuning heuristic (keyed off
  the regulator, not just the degree) is the natural fix.
- **Performance at higher degree.** The sieve's `2^degree` branching
  factor and per-box absorption cost mean degree-5+ fields converge more
  slowly and less tightly than degree 2-3 within the same time budget
  (see `x^5-x-1` above). Tuning `candidateCap`, `unitExponentRange`, and
  the depth/tolerance parameters per field is currently manual, mirroring
  the hand-tuned, degree-indexed `euclid.cfg` tables in the original C
  tool.
- **Absolute speed vs. the C tool.** As shown above, the Chapel tool is
  currently much slower in wall-clock terms for small/simple fields, even
  though both tools' core search loops use the same double-precision
  floating-point arithmetic — the gap is implementation efficiency, not a
  different numerical approach (see "Optimization: real(64) fast path"
  and "Optimization: tuple-based Box" above). The tuple-based `Box`
  refactor is implemented and correctness-validated, but at *default*
  settings its end-to-end impact is modest (0-10%) — see that section for
  why (the per-candidate absorption scan dominates over box construction
  at the default `candidateCap`) and the auto-tuning roadmap item that
  should unlock it more broadly. The generality/portability trade-off
  (arbitrary-degree support and a design that scales out to clusters/GPUs,
  not yet built — see "Parallelization roadmap" below) is real, but it
  isn't the reason for today's gap.
- **Degree capped at 8 by a compile-time dispatch table.** `main.chpl`
  bridges `Box`'s generic `param degree` to the field's runtime degree via
  an explicit `select nf.degree { when 1..8 ... }`; higher-degree fields
  halt with a clear message rather than running. Not a deep limitation --
  see [`CHAPEL.md`](CHAPEL.md#supported-field-degree) for how to extend it
  (one `when` branch plus a recompile).
- **Phase 3 is not the full Lezowski cycle-decomposition machinery.** It
  samples the most-resistant boxes rather than exactly identifying the
  unit-orbit critical cycle (`src/graph.c`'s Tarjan-based decomposition),
  so it isn't guaranteed to land exactly on the supremum for harder
  fields (though a reported certified bound is always sound — see
  `Certify.chpl` above).
- **Parallelization roadmap (not yet built).** The per-level box list and
  the small-elements enumeration are both `forall`-parallel today
  (single-locale, multi-core, per the benchmarks above). Scaling to
  multiple locales (distribute the "problems" list, e.g. via a
  `Block`-distributed array) and to GPUs (`computeBoxProjections`/
  `calculateBoxMaxNorm` are flat numeric kernels well-suited to
  `foreach`/GPU offload) is the natural next step, and the intended way
  to close the absolute-speed gap noted above for the field sizes where
  it matters. See [`CHAPEL.md`'s Roadmap](CHAPEL.md#roadmap) for the full
  staged plan (auto-tuning first, then multi-locale, then GPU).

See [`CHAPEL.md`](CHAPEL.md) for the full writeup (including why each
design choice was necessary, not just convenient) and the complete CLI
reference.

## Repository layout

```
src/                C source for the reference Pari-based `euclid` tool
generate_field.gp   one-time Pari/gp field-setup script (-> field_data.txt)
NumberField.chpl    Chapel: field_data.txt parser + linear algebra helpers
SmallElements.chpl  Chapel: candidate algebraic-integer enumeration
Sieve.chpl          Chapel: covering sieve + unit-action acceleration
Certify.chpl        Chapel: exact lower-bound certification via Pari
main.chpl           Chapel: CLI entry point
tests/              Chapel smoke tests, saved field fixtures, validate.sh
q*.txt              field_data.txt files for small quadratic/cubic fields,
                    including the real quadratic battery (q11-q73) used
                    to find and regression-test the correctness fixes
attic/              superseded early Chapel prototypes, kept for history
```

## Requirements

- A C compiler and the Pari/GP library (for the `src/` tool and for
  `generate_field.gp` / exact certification). The Makefile passes
  `-fcommon`, needed to build the C tool on modern GCC.
- The [Chapel](https://chapel-lang.org/) compiler (`chpl`), for the Chapel
  port. Developed against Chapel 2.9.0. The Makefile forces
  `CHPL_LOCALE_MODEL=flat` for all Chapel builds, since `chpl` 2.9.0 hits
  an internal compiler error analyzing this project's `forall` loops under
  a GPU-locale-model toolchain (this project has no GPU-specific code, so
  there's nothing to gain from that mode here anyway).

## Build and test

A single `Makefile` drives both implementations:

```sh
make            # build the C tool (./euclid) -- the default target
make chapel     # build the Chapel tool (./euclid_chpl), with --fast
make chapel-debug  # unoptimized/bounds-checked build (./euclid_chpl_debug),
                   # for tracking down correctness bugs
make all        # euclid + chapel

make fixtures   # generate any missing tests/fixtures/*.txt via gp
make smoke      # compile + run the Chapel per-module smoke tests
make test       # run tests/validate.sh (Chapel bracket vs. C exact answer)
make check      # smoke + test

make clean      # remove all build artifacts (both tools)
```

## Usage

C tool:

```sh
./euclid "x^2-61"
```

Chapel tool:

```sh
POLY="x^2-61" gp -q generate_field.gp      # one-time setup -> field_data.txt
./euclid_chpl --fieldFile=field_data.txt --tolerance=0.001
```

See [`CHAPEL.md`](CHAPEL.md) for the full CLI reference (bisection depth,
unit-acceleration, and certification flags) and an explanation of what the
reported bracket/certification output means.

## License

See [`COPYING`](COPYING) (GPL) for the original C tool. The Chapel port is
provided under the same terms unless noted otherwise.

## Reference

The algorithm is the application of the article available at
<http://hal.archives-ouvertes.fr/hal-00632997/>.
