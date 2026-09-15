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
  arbitrary-degree reimplementation of the same algorithm that keeps all
  numerically-heavy work in floating point / Chapel, parallelized with
  `forall`, and uses Pari/gp only for one-time per-field setup and an
  optional final exact-certification step. Builds to the `euclid_chpl`
  binary. See [`CHAPEL.md`](CHAPEL.md) for full architecture, algorithm,
  and usage details.

The Chapel port is validated against the C tool as its ground-truth oracle
(`tests/validate.sh`); see [`CHAPEL.md`](CHAPEL.md#validation) for results.

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
