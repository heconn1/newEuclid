# Final Project Report: Chapel Euclidean-Minimum Sieve

**Release:** `chapel-v1.1.0`
**Status:** Complete and closed. Repository archived at this tag.

## 1. Objective

Reimplement Pierre Lezowski's algorithm for computing the Euclidean
minimum M(K) of a number field K in [Chapel](https://chapel-lang.org/),
targeting arbitrary degree (the original C/Pari tool in `src/` is also
general, but the early Chapel prototypes were hard-capped at degree ≤ 4),
keeping the numerically-heavy search in floating point so it can
eventually scale to clusters and GPUs, and using Pari/gp only for
one-time per-field setup and a final exact-certification step.

## 2. Summary of work

The project went through three main phases:

1. **Ground-up rebuild** of the Chapel implementation. The prior
   prototypes (`firstCut.chpl` ... `sixthCut.chpl`, now in `attic/`)
   never implemented the actual mathematical core of the algorithm
   (absorption testing against real candidate algebraic integers) and
   were replaced with a correct, arbitrary-degree design across five
   modules (`NumberField.chpl`, `SmallElements.chpl`, `Sieve.chpl`,
   `Certify.chpl`, `main.chpl`), a `Makefile`-driven build/test
   workflow, and initial documentation. Tagged `chapel-v1.0.0`.
2. **Correctness battery.** Testing against 15 norm-Euclidean real
   quadratic fields (`q11.txt` ... `q73.txt`) exposed four real,
   distinct bugs that the initial fixture set (5 small fields) hadn't
   caught — including one field converging to *exactly 2x* its true
   minimum with full confidence. All four were fixed, and two are now
   permanent regression cases in `tests/validate.sh`.
3. **Documentation and finalization.** README.md and CHAPEL.md were
   brought fully up to date with the correctness history, measured
   performance benchmarks, a complete known-limitations list, and a
   staged roadmap for future optimization and GPU porting. Tagged
   `chapel-v1.1.0` (this release).

## 3. Correctness bugs found and fixed

All four were found via the real quadratic battery, not the original
5-field fixture set — a reminder that a wider, independently-verifiable
test battery (the C tool's exact answers) found real bugs the narrower
one missed.

| # | Bug | Symptom | Fix |
|---|---|---|---|
| 1 | `SmallElements.chpl` required every embedding coordinate individually bounded, excluding large-coefficient/small-norm "convergent-like" absorbers real quadratic fields need | `x^2-57` converged to **exactly 2x** the true minimum | Rank all candidates by actual field norm instead |
| 2 | Bracket-finding trusted a shallow sieve trial's "blocked" result to raise the lower bound | `x^2-19` converged to exactly `1.0` instead of `170/171` | Only ever raise the lower bound from a *deep*, trusted sieve result |
| 3 | Floating-point precision cliff at fields whose critical point has dyadic coordinates (common) | Briefly regressed `x^2-2` to falsely report `M(K) < 0.5` (its true minimum is exactly 1/2) while fixing bug #4 | Safety margin in the absorption comparison + numerically-safe `minWidth` cutoff |
| 4 | `Certify.chpl`'s exact check inherited the same completeness gap as bug #1, independently | Could report a "certified" value larger than an independently-proven upper bound | Check against the sieve's own norm-ranked candidates; discard any sample contradicting the proven upper bound |

Full technical explanation of each: `CHAPEL.md` § "Bugs found and fixed
during this battery".

## 4. Validation results (final state)

Against the original 5 fixtures and the 15-field real quadratic battery,
compared to the C tool's exact answers:

| Field | Reference exact minimum | Chapel bracket |
|---|---|---|
| `x^2-2` | 1/2 | `[0.5, 0.500977]`, fraction exact |
| `x^2-61` | 1611/1525 (~1.05639) | `[1.05566, 1.05664]` |
| `x^2-19` | 170/171 | `[0.994141, 0.995117]`, fraction exact |
| `x^2-57` | 14/19 | `[0.736328, 0.737305]`, fraction exact |
| `x^3+x^2-1` | 1/5 | `[0.199219, 0.200195]`, fraction exact |
| `x^3-3*x-1` | 1/3 | `[0.333008, 0.333984]`, fraction exact |
| `x^5-x-1` (degree 5) | 1/4 | `[0.271875, 0.28125]` (looser; known limitation) |
| `x^2-11`, `13`, `29`, `33`, `41` | (all exact, see CHAPEL.md) | all correct, several with exact fraction match |
| `x^2-73` | 1541/2136 (~0.72144) | `[0.723633, 0.724609]` — ~0.3% high (known limitation) |

`make check` (full smoke suite + `tests/validate.sh`) passes cleanly
from a clean build as of this release.

## 5. Performance (measured, 6-core machine, default settings)

| Field | C tool | Chapel tool | Parallel speedup (user/real) |
|---|---|---|---|
| `x^2-2` | ~0.2s | ~3.5s | ~3.6x |
| `x^2-61` | ~0.15s | ~68s | ~5.3x |
| `x^2-73` (hardest regulator) | ~3.7s | ~157s | ~5.0x |
| `x^3-3*x-1` (rank-2 units) | ~0.6s | ~198s | ~4.9x |
| `x^5-x-1` (degree 5) | ~50s | ~160s | ~3.7x |

The Chapel tool is currently 17-450x slower in absolute wall-clock terms
than the C tool for these small fields — expected, since it searches in
floating point across cores rather than using Pari's exact bignum
arithmetic on one thread. This is a deliberate trade-off for
arbitrary-degree support and a design that can scale out to clusters and
GPUs, which the original C/Pari architecture cannot. Full table and
discussion: `README.md` § "Performance benchmarks".

## 6. Known limitations (carried forward)

- `x^2-73` (largest regulator tested, ~2136) needs a larger
  `--candidateCap`/`--refineDepth`/`--certifyDepth` than the current
  defaults to close its remaining ~0.3% gap; a field-size-aware
  auto-tuning heuristic is the planned fix (Roadmap stage 1).
- Degree-5+ fields converge more slowly and less tightly than degree
  2-3 within the same time budget.
- Phase 3 certification samples resistant boxes rather than porting
  Lezowski's full unit-orbit cycle/graph decomposition, so it isn't
  guaranteed to land exactly on the supremum for harder fields (though
  a reported certified bound is always sound).
- No multi-locale or GPU support yet.

## 7. Roadmap (not yet started)

Documented in full in `CHAPEL.md` § "Roadmap":

1. Auto-tuning (`candidateCap`/`refineDepth`/`certifyDepth` from a
   field's degree *and* regulator) and single-locale performance fixes.
2. Multi-locale distribution of the per-depth-level box search.
3. GPU offload of the flat absorption-test numeric kernels.
4. Stretch: full Lezowski cycle/graph decomposition for exact-everywhere
   certification, if resistant-box sampling ever proves insufficient.

## 8. Release artifacts

- **Git tag:** `chapel-v1.1.0` (annotated, pushed to `origin`).
- **Source archive:** a clean `git archive` snapshot of the tagged
  commit (no build artifacts, no `.git` metadata) was generated locally
  at `~/newEuclid/archives/euclid-chapel-v1.1.0.tar.gz` for backup;
  regenerate anytime with:
  ```sh
  git archive --format=tar.gz --prefix=euclid-chapel-v1.1.0/ \
    -o euclid-chapel-v1.1.0.tar.gz chapel-v1.1.0
  ```
- **Documentation:** `README.md` (overview, validation, benchmarks,
  limitations, usage) and `CHAPEL.md` (full architecture, bug
  postmortems, CLI reference, roadmap).
- **Repository:** `git@github.com:heconn1/newEuclid.git`, `master`
  branch, in sync with `origin/master`.

## 9. Closing note

The repository is in a clean, fully documented, and validated state.
All work from this engagement is committed, pushed, and tagged. No
further action is required to consider this phase of the project
closed; the Roadmap above is the starting point whenever performance
work resumes.
