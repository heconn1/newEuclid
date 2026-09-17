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
  attribution and license.
- **A Chapel port** (`*.chpl` at the repository root) — a from-scratch,
  arbitrary-degree reimplementation of the same algorithm, documented in
  full below and in [`CHAPEL.md`](CHAPEL.md).

## Chapel implementation

The Chapel port computes the same M(K) as the C tool, keeping all
numerically-heavy work in floating point / Chapel (parallelized with
`forall`), and uses Pari/gp only for one-time per-field setup and an
optional final exact-certification step. It builds to the `euclid_chpl`
binary.

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
  absorbing algebraic integers for a target bound K.
- `Sieve.chpl` — the real absorption test and parallel box-bisection
  sieve, plus **unit-action acceleration**: proven necessary (not just an
  optimization) for fields with large regulators, such as `x^2-61`.
- `Certify.chpl` — turns a numeric upper bound into a genuine
  Pari-certified lower bound by evaluating the exact norm at the
  most-resistant sampled points. Its brute-force search range scales with
  degree (a flat range is provably too small for some fields, e.g.
  `x^2-61` needs `>= 5`) and also tries small unit powers, and `main.chpl`
  additionally discards any sampled value that's `>=` the proven numeric
  upper bound, since that's provably impossible — so a reported certified
  bound is always sound, even on fields the search still isn't tight for.
- `main.chpl` — the CLI: exponential bracket search, bisection
  refinement, and an optional certification pass.

**Algorithm in three phases:** (1) a sound covering sieve proves numeric
upper bounds on M(K); (2) unit-action recentering lets a small,
origin-centered candidate list stand in for absorbers spread out along a
field's unit orbit, which large-regulator fields require to converge at
all; (3) exact Pari evaluation at resistant sample points certifies a
matching lower bound when possible.

The Chapel port is validated against the C tool as its ground-truth oracle
(`tests/validate.sh`):

| Field | Reference exact minimum | Chapel bracket |
|---|---|---|
| `x^2-2` | 1/2 | `[0.5, 0.500977]` |
| `x^2-61` | 1611/1525 (~1.05639) | `[1.05762, 1.05859]` |
| `x^3+x^2-1` | 1/5 | `[0.199219, 0.200195]` |
| `x^3-3*x-1` | 1/3 | `[0.332897 (Pari-certified), 0.333398]` |
| `x^5-x-1` | 1/4 | converges but looser; demonstrates degree-5 support |

`x^2-61` and `x^3-3*x-1` are specifically the cases the earlier prototypes
could not solve correctly (large regulator, and rank-2 unit group
respectively).

**Known limitations:** degree-5+ fields converge more slowly/loosely than
degree 2-3 within the same time budget; Phase 3 samples resistant boxes
rather than porting Lezowski's full unit-orbit cycle/graph decomposition,
so it isn't guaranteed to land exactly on the supremum for harder fields
(though it is always sound — see `Certify.chpl` above); multi-locale/GPU
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
attic/              superseded early Chapel prototypes, kept for history
```

## Requirements

- A C compiler and the Pari/GP library (for the `src/` tool and for
  `generate_field.gp` / exact certification).
- The [Chapel](https://chapel-lang.org/) compiler (`chpl`), for the Chapel
  port. Developed against Chapel 2.9.

## Build and test

A single `Makefile` drives both implementations:

```sh
make            # build the C tool (./euclid) -- the default target
make chapel     # build the Chapel tool (./euclid_chpl)
make all        # build both

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
