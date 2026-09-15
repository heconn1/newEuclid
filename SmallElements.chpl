// SmallElements.chpl
//
// Arbitrary-degree generalization of Lezowski's `small_elts` (src/main.c).
// Enumerates integer coefficient vectors (relative to the integral basis)
// that are plausible "absorbers": elements gamma in O_K whose embedding is
// small enough in at least one place that N(x - gamma) could plausibly be
// below a target bound K for some x in the fundamental domain.
//
// This is only ever used as a *sound-but-not-necessarily-complete* filter:
// omitting a genuinely useful candidate only costs efficiency (a box that
// could have been proven covered stays a "problem" and gets bisected
// further); it can never make the sieve unsound, because the final
// covering test (Sieve.chpl) independently re-verifies each candidate.
module SmallElements {
  use NumberField;
  use List;
  use Math;
  use Sort;

  record SmallElementSet {
    var n: int; // number of candidates
    var degree: int;
    var numEmbeddings: int;
    var cDom: domain(2) = {1..0, 1..0};
    var coeffs: [cDom] real(64);       // [candidate, coeff index]
    var eDom: domain(2) = {1..0, 1..0};
    var emb: [eDom] complex(128);      // [candidate, embedding index]
  }

  // Choose the largest symmetric per-coordinate integer range B such that
  // enumerating (2B+1)^degree candidates stays within `budget` -- this is
  // the formulaic analogue of the hand-tuned, degree-indexed BOUNDS/CUTTING
  // tables in euclid.cfg, generalized to arbitrary degree.
  proc candidateBoundRange(degree: int, budget: real = 200_000.0): int {
    if degree <= 0 then return 1;
    const raw = (budget ** (1.0/degree) - 1.0) / 2.0;
    var b = floor(raw): int;
    if b < 1 then b = 1;
    if b > 300 then b = 300;
    return b;
  }

  // Enumerate candidates for target bound K. `boundRange` overrides the
  // automatic budget-based choice when > 0 (useful for tuning/testing).
  // `candidateCap` bounds how many candidates are ultimately kept (the
  // smallest-norm ones survive) -- absorption cost is O(candidates) per
  // box, so this keeps per-box cost bounded regardless of degree/boundRange.
  proc smallElements(const ref nf: NumberFieldData, K: real(64), eps: real(64) = 1.0e-5,
                      boundRange: int = 0, candidateCap: int = 500): SmallElementSet {
    const degree = nf.degree;
    const numEmbeddings = nf.numEmbeddings;
    const B = if boundRange > 0 then boundRange else candidateBoundRange(degree);
    const span = 2*B + 1;

    // x = K^(1/degree) + eps; R[i] = x + sum_k |basis[i,k]| bounds how far
    // an embedding coordinate of a "plausible absorber" can be from 0.
    const x = K ** (1.0/degree) + eps;
    var R: [1..numEmbeddings] real(64);
    forall i in 1..numEmbeddings {
      var s = 0.0;
      for c in 1..degree do s += abs(nf.basis[i, c]);
      R[i] = x + s;
    }

    var total = 1: int;
    for c in 1..degree do total *= span;

    var found: list((int, [1..degree] real(64), [1..numEmbeddings] complex(128)), parSafe=true);

    forall idx in 0..#total with (ref found) {
      var w: [1..degree] real(64);
      var rem = idx;
      for c in 1..degree {
        const digit = rem % span;
        rem /= span;
        w[c] = (digit - B): real(64);
      }
      var e: [1..numEmbeddings] complex(128);
      var keep = true;
      for row in 1..numEmbeddings {
        var s: complex(128) = 0.0;
        for c in 1..degree do s += nf.basis[row, c] * w[c];
        e[row] = s;
        if abs(s) > R[row] then keep = false;
      }
      if keep then found.pushBack((0, w, e));
    }

    const pool = found.toArray();
    var keepIdx: [0..#(min(pool.size, candidateCap))] int;
    if pool.size > candidateCap {
      var keys: [pool.domain] (real(64), int);
      forall (k, idx) in zip(keys.domain, pool.domain) {
        const (_, w, e) = pool[idx];
        keys[k] = (absNorm(nf, e), idx);
      }
      sort(keys);
      for i in keepIdx.domain do keepIdx[i] = keys[i][1];
    } else {
      for i in keepIdx.domain do keepIdx[i] = pool.domain.low + i;
    }

    var result: SmallElementSet;
    result.degree = degree;
    result.numEmbeddings = numEmbeddings;
    result.n = keepIdx.size;
    result.cDom = {1..result.n, 1..degree};
    result.eDom = {1..result.n, 1..numEmbeddings};
    for i in 1..result.n {
      const (_, w, e) = pool[keepIdx[i-1]];
      for c in 1..degree do result.coeffs[i, c] = w[c];
      for row in 1..numEmbeddings do result.emb[i, row] = e[row];
    }
    return result;
  }
}
