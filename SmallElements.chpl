// SmallElements.chpl
//
// Arbitrary-degree generalization of Lezowski's `small_elts` (src/main.c).
// Enumerates integer coefficient vectors (relative to the integral basis)
// and keeps the smallest-norm ones as plausible "absorbers" gamma in O_K.
//
// IMPORTANT: the filter used to require *every* embedding coordinate to be
// individually bounded (an axis-aligned ball around the origin). That's
// wrong for fields with a nontrivial unit group: the genuinely useful small
// absorbers for e.g. a real quadratic field with a large fundamental unit
// are continued-fraction-convergent-like elements whose *coefficients* can
// be large even though their *norm* is small (that's exactly what makes an
// element "small" in the relevant sense) -- a per-coordinate ball
// systematically excludes them, which was observed to make the sieve
// converge confidently on a wrong (too-large) answer for real quadratic
// fields with sizable regulators. We instead rank all enumerated candidates
// by actual norm and keep the smallest `candidateCap` of them.
//
// This is only ever used as a *sound-but-not-necessarily-complete* filter:
// omitting a genuinely useful candidate only costs efficiency (a box that
// could have been proven covered stays a "problem" and gets bisected
// further); it can never make the sieve unsound, because the final
// covering test (Sieve.chpl) independently re-verifies each candidate.
module SmallElements {
  use NumberField;
  use Math;

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

  // Enumerate candidates in [-B,B]^degree and keep the `candidateCap`
  // smallest-norm ones. `boundRange` overrides the automatic budget-based
  // choice of B when > 0 (useful for tuning/testing); `K` is currently
  // unused for filtering (norm ranking alone decides what's kept) but is
  // kept as a parameter since callers key their candidate sets by it.
  proc smallElements(const ref nf: NumberFieldData, K: real(64), eps: real(64) = 1.0e-5,
                      boundRange: int = 0, candidateCap: int = 2000): SmallElementSet {
    const degree = nf.degree;
    const numEmbeddings = nf.numEmbeddings;
    const B = if boundRange > 0 then boundRange else candidateBoundRange(degree);
    const span = 2*B + 1;

    var total = 1: int;
    for c in 1..degree do total *= span;

    // Rank every enumerated candidate by |N(w)| and keep the smallest
    // `candidateCap`, via a single shared sorted top-cap array guarded by a
    // lock (see insertCandidate) -- avoids materializing all `total`
    // candidates before ranking, which matters once `total` is large.
    const cap = max(candidateCap, 1);
    var bestKeys: [0..#cap] real(64) = max(real(64));
    var bestCoeffs: [0..#cap, 1..degree] real(64);
    var bestEmb: [0..#cap, 1..numEmbeddings] complex(128);
    var bestCount = 0;

    forall idx in 0..#total with (ref bestKeys, ref bestCoeffs, ref bestEmb, ref bestCount) {
      var w: [1..degree] real(64);
      var rem = idx;
      for c in 1..degree {
        const digit = rem % span;
        rem /= span;
        w[c] = (digit - B): real(64);
      }
      var e: [1..numEmbeddings] complex(128);
      for row in 1..numEmbeddings {
        var s: complex(128) = 0.0;
        for c in 1..degree do s += nf.basis[row, c] * w[c];
        e[row] = s;
      }
      const n = absNorm(nf, e);
      // Insert into this task's view of the global top-cap list under a
      // lock; cap is small (hundreds), so contention is not a bottleneck
      // relative to the O(total) work above.
      insertCandidate(bestKeys, bestCoeffs, bestEmb, bestCount, n, w, e);
    }

    const keepN = min(bestCount, cap);
    var result: SmallElementSet;
    result.degree = degree;
    result.numEmbeddings = numEmbeddings;
    result.n = keepN;
    result.cDom = {1..result.n, 1..degree};
    result.eDom = {1..result.n, 1..numEmbeddings};
    for i in 1..result.n {
      for c in 1..degree do result.coeffs[i, c] = bestCoeffs[i-1, c];
      for row in 1..numEmbeddings do result.emb[i, row] = bestEmb[i-1, row];
    }
    return result;
  }

  // Maintains a small sorted (ascending by key) top-`cap` list in place,
  // guarded by a single global lock. Simpler and plenty fast enough here
  // (cap is at most a few thousand) compared to a lock-free structure.
  private var _insertLock: sync bool = true;
  private proc insertCandidate(ref keys: [] real(64), ref coeffs: [] real(64), ref emb: [] complex(128),
                                ref count: int, key: real(64), const ref w: [] real(64),
                                const ref e: [] complex(128)) {
    const cap = keys.domain.size;
    _insertLock.readFE();
    if count < cap || key < keys[cap-1] {
      var pos = min(count, cap-1);
      while pos > 0 && keys[pos-1] > key {
        keys[pos] = keys[pos-1];
        for c in coeffs.domain.dim(1) do coeffs[pos, c] = coeffs[pos-1, c];
        for r in emb.domain.dim(1) do emb[pos, r] = emb[pos-1, r];
        pos -= 1;
      }
      keys[pos] = key;
      for c in coeffs.domain.dim(1) do coeffs[pos, c] = w[c];
      for r in emb.domain.dim(1) do emb[pos, r] = e[r];
      if count < cap then count += 1;
    }
    _insertLock.writeEF(true);
  }
}
