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
  point can't produce a false "cleared" result.
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

**Known limitations:** `x^2-73` (fundamental unit ~2136, the largest
regulator tested) needs more `--candidateCap`/depth than the current
defaults budget for to close its last ~0.3% gap; degree-5+ fields converge
more slowly/loosely than degree 2-3 within the same time budget; Phase 3
samples resistant boxes rather than porting Lezowski's full unit-orbit
cycle/graph decomposition, so it isn't guaranteed to land exactly on the
supremum for harder fields (though it is always sound); multi-locale/GPU
scaling is designed for but not yet built.

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
