// Sieve.chpl
//
// Arbitrary-degree, parallel port of Lezowski's core decision procedure:
// "is K an upper bound for the Euclidean minimum of this field?"
//
// A box is a coefficient-space (relative to the integral basis) hyper-cube.
// It is *absorbed* (safely discarded) once we exhibit a candidate gamma in
// O_K such that every point of the box has |N(x - gamma)| < K -- proven via
// an interval/disk bound on the box's projection into embedding space
// (computeBoxProjections/calculateBoxMaxNorm below), generalizing
// Lezowski's `norme_bricolee` (src/basef.c) to arbitrary r1,r2.
//
// Boxes that cannot be absorbed are bisected into 2^degree children and
// retried at finer resolution. If every box is eventually absorbed, K is a
// *proven* upper bound on the Euclidean minimum. Exhausting the resolution
// budget without absorbing everything means K is either too small or needs
// more depth to resolve -- see Certify.chpl (Phase 3) for how to turn a
// surviving "problem" box into a certified exact answer.
//
// Box is generic over a compile-time `param degree`, storing its
// coefficients as a fixed-size tuple rather than a domain-backed array --
// measured 6-8x cheaper to construct/copy than the array-based equivalent
// (see bench/BoxCopyBench.chpl), which matters a great deal since
// bisect() allocates 2^degree fresh boxes at every depth level for every
// unresolved box. Any proc taking `const ref box: Box(?d)` infers `d`
// automatically from the argument; only runSieve (which builds the root
// box directly from a runtime nf.degree) needs an explicit `param degree`
// formal, supplied by its caller at compile time -- see main.chpl's
// `select nf.degree { ... }` dispatch.
//
// Indexing convention: box.center/box.widths are 0-indexed tuples
// (Chapel's tuple convention), while NumberFieldData.basis/basisReal
// remain 1-indexed domain arrays (unrelated data, deliberately left
// alone). Every loop below keeps its index variable 1-based (matching
// nf.basis's convention) and only ever offsets by one at the point of
// tuple access (box.center(col-1)), never anywhere else.
module Sieve {
  use NumberField;
  use SmallElements;
  use List;
  use Math;

  record Box {
    param degree: int;
    var center: degree*real(64);
    var widths: degree*real(64);
  }

  // Converts a Box's tuple center to a 1-indexed array, for the one place
  // a box center needs to cross into Certify.chpl's still-array-based API.
  proc boxCenterToArray(const ref box: Box(?d)): [1..d] real(64) {
    var a: [1..d] real(64);
    for i in 1..d do a[i] = box.center(i-1);
    return a;
  }

  proc computeBoxProjections(const ref nf: NumberFieldData, const ref box: Box(?d),
                              ref centers: [] complex(128), ref radii: [] real(64)) {
    forall row in 1..nf.numEmbeddings {
      var c: complex(128) = 0.0;
      var r = 0.0;
      for col in 1..d {
        c += nf.basis[row, col] * box.center(col-1);
        r += abs(nf.basis[row, col]) * box.widths(col-1);
      }
      centers[row] = c;
      radii[row] = r;
    }
  }

  proc calculateBoxMaxNorm(const ref nf: NumberFieldData, const ref centers: [] complex(128),
                            const ref radii: [] real(64)): real(64) {
    var maxNorm = 1.0;
    for i in 1..nf.r1 do maxNorm *= (abs(centers[i]) + radii[i]);
    for j in 1..nf.r2 {
      const idx = nf.r1 + j;
      const m = abs(centers[idx]) + radii[idx];
      maxNorm *= m*m;
    }
    return maxNorm;
  }

  // True if some enumerated small-element candidate provably absorbs the box.
  // Dispatches to a real(64)-only fast path for totally real fields (r2 ==
  // 0, e.g. every real quadratic field), which skips complex(128)
  // arithmetic entirely in this per-candidate scan -- the hottest loop in
  // the program. Fields with r2 > 0 use the original complex(128) path.
  proc isBoxAbsorbed(const ref nf: NumberFieldData, const ref box: Box(?d),
                      const ref candidates: SmallElementSet, K: real(64)): bool {
    if nf.isTotallyReal then return isBoxAbsorbedReal(nf, box, candidates, K);
    var centers: [1..nf.numEmbeddings] complex(128);
    var radii: [1..nf.numEmbeddings] real(64);
    computeBoxProjections(nf, box, centers, radii);
    return isProjectionAbsorbed(nf, centers, radii, candidates, K);
  }

  // Small relative safety margin subtracted from K before comparing: a box
  // whose true maximum norm is *exactly* K (e.g. one edge sits precisely on
  // a genuine critical point of the field, which happens whenever that
  // critical point's coefficients are dyadic, as is common) can have its
  // computed maxNorm come out a few ULPs *below* K purely from
  // floating-point rounding in the embedding/subtraction/multiplication
  // chain -- this was observed to falsely "absorb" the exact-boundary box
  // for x^2-2 at K=0.5 once bisection went deep enough to expose it. The
  // margin is relative (scaled to K) so it stays meaningful whether K is
  // order 1 or order 1e-3.
  private const absorptionSafetyMargin = 1.0e-9;

  private proc isProjectionAbsorbed(const ref nf: NumberFieldData, const ref centers: [] complex(128),
                                     const ref radii: [] real(64), const ref candidates: SmallElementSet,
                                     K: real(64)): bool {
    const threshold = K * (1.0 - absorptionSafetyMargin);
    for i in 1..candidates.n {
      var shifted: [1..nf.numEmbeddings] complex(128);
      for row in 1..nf.numEmbeddings do shifted[row] = centers[row] - candidates.emb[i, row];
      if calculateBoxMaxNorm(nf, shifted, radii) < threshold then return true;
    }
    return false;
  }

  // Real(64)-only fast path used when nf.isTotallyReal (r2 == 0): every
  // embedding is real, so this per-candidate absorption scan never needs
  // to pay for complex(128) multiplies whose imaginary part is always
  // exactly zero. Mirrors isBoxAbsorbed/computeBoxProjections/
  // calculateBoxMaxNorm/isProjectionAbsorbed above exactly, just typed
  // real(64) throughout.
  proc isBoxAbsorbedReal(const ref nf: NumberFieldData, const ref box: Box(?d),
                          const ref candidates: SmallElementSet, K: real(64)): bool {
    var centers: [1..nf.r1] real(64);
    var radii: [1..nf.r1] real(64);
    computeBoxProjectionsReal(nf, box, centers, radii);
    return isProjectionAbsorbedReal(nf, centers, radii, candidates, K);
  }

  proc computeBoxProjectionsReal(const ref nf: NumberFieldData, const ref box: Box(?d),
                                  ref centers: [] real(64), ref radii: [] real(64)) {
    forall row in 1..nf.r1 {
      var c = 0.0;
      var r = 0.0;
      for col in 1..d {
        c += nf.basisReal[row, col] * box.center(col-1);
        r += abs(nf.basisReal[row, col]) * box.widths(col-1);
      }
      centers[row] = c;
      radii[row] = r;
    }
  }

  proc calculateBoxMaxNormReal(const ref centers: [] real(64), const ref radii: [] real(64)): real(64) {
    var maxNorm = 1.0;
    for i in centers.domain do maxNorm *= (abs(centers[i]) + radii[i]);
    return maxNorm;
  }

  private proc isProjectionAbsorbedReal(const ref nf: NumberFieldData, const ref centers: [] real(64),
                                         const ref radii: [] real(64), const ref candidates: SmallElementSet,
                                         K: real(64)): bool {
    const threshold = K * (1.0 - absorptionSafetyMargin);
    for i in 1..candidates.n {
      var shifted: [1..nf.r1] real(64);
      for row in 1..nf.r1 do shifted[row] = centers[row] - candidates.embReal[i, row];
      if calculateBoxMaxNormReal(shifted, radii) < threshold then return true;
    }
    return false;
  }

  // Unit-action acceleration (generalizes Lezowski's vecteurs_possibles2 /
  // test_des_unites, src/main.c:676-897): multiplying every point of a box
  // by a unit u and re-rounding to the fundamental domain gives a *sound*
  // way to reuse a small, origin-centered candidate list for boxes far
  // from the origin -- because for any algebraic integer gamma,
  //   |N(u*x - gamma)| = |N(x - gamma*u^-1)|,
  // so proving the twisted box is absorbed by `gamma` proves the original
  // box is absorbed by the (generally much larger/farther) integer
  // gamma*u^-1, without ever having to enumerate such far-away integers
  // directly. This is essential for fields with large regulators (e.g.
  // real quadratic fields with a large fundamental unit), where the useful
  // absorbers are spread out along the unit orbit rather than clustered
  // near the origin.
  proc unitPower(const ref nf: NumberFieldData, unitIdx: int, exponent: int): [1..nf.numEmbeddings] complex(128) {
    var p: [1..nf.numEmbeddings] complex(128) = 1.0;
    const e = abs(exponent);
    for _i in 1..e {
      forall row in 1..nf.numEmbeddings {
        if exponent > 0 then p[row] *= nf.unitEmbeddings[row, unitIdx];
        else p[row] /= nf.unitEmbeddings[row, unitIdx];
      }
    }
    return p;
  }

  proc recenterByUnit(const ref nf: NumberFieldData, const ref box: Box(?d),
                       const ref centers: [] complex(128), const ref radii: [] real(64),
                       const ref unitPow: [] complex(128)): Box(d) {
    var tCenters: [1..nf.numEmbeddings] complex(128);
    var tRadii: [1..nf.numEmbeddings] real(64);
    forall row in 1..nf.numEmbeddings {
      tCenters[row] = centers[row] * unitPow[row];
      tRadii[row] = radii[row] * abs(unitPow[row]);
    }

    const w = coefficientsFromEmbedding(nf, tCenters);
    var t: [1..nf.degree] real(64);
    for i in 1..nf.degree do t[i] = round(w[i]);

    var newCenter: d*real(64);
    for i in 1..d do newCenter(i-1) = w[i] - t[i];

    // Pull the (disk) embedding radius back into coefficient space via the
    // same triangle-inequality bound used going forward, applied to
    // sigmaInv; a disk of radius R is conservatively bounded by a
    // 2R x 2R square in the unpacked real coordinates.
    var unpackedRadius: [1..nf.degree] real(64);
    for i in 1..nf.r1 do unpackedRadius[i] = tRadii[i];
    for j in 1..nf.r2 {
      unpackedRadius[nf.r1 + j] = tRadii[nf.r1 + j];
      unpackedRadius[nf.r1 + nf.r2 + j] = tRadii[nf.r1 + j];
    }
    // Serial (not forall): writes into the local tuple newWidths, which
    // (unlike an array) has value semantics, so a parallel loop over it
    // would need an explicit `with (ref newWidths)` -- degree is at most
    // 8, so there's nothing to gain from parallelizing this anyway.
    var newWidths: d*real(64);
    for c in 1..d {
      var s = 0.0;
      for k in 1..nf.degree do s += abs(nf.sigmaInv[c, k]) * unpackedRadius[k];
      newWidths(c-1) = s;
    }

    return new Box(degree=d, center=newCenter, widths=newWidths);
  }

  // Direct absorption, then (if that fails) a sweep of fundamental unit
  // powers with recentering. `unitExponentRange` controls how many powers
  // of each fundamental unit are tried (both signs); this needs to be
  // large enough to reach whatever power of the unit is relevant for a
  // given field's regulator (see runSieve for how the default scales with
  // the field's own unit magnitudes).
  proc isBoxAbsorbedWithUnits(const ref nf: NumberFieldData, const ref box: Box(?d),
                               const ref candidates: SmallElementSet, K: real(64),
                               unitExponentRange: int = 1): bool {
    // Unit multiplication is inherently complex (unitEmbeddings are stored
    // complex(128) regardless of isTotallyReal), so recenterByUnit always
    // needs the complex projections of the *original* box. The fast
    // real(64) path only replaces the expensive per-candidate absorption
    // scan itself, both for the direct box and for each unit-twisted box.
    var centers: [1..nf.numEmbeddings] complex(128);
    var radii: [1..nf.numEmbeddings] real(64);
    computeBoxProjections(nf, box, centers, radii);
    if nf.isTotallyReal {
      if isBoxAbsorbedReal(nf, box, candidates, K) then return true;
    } else {
      if isProjectionAbsorbed(nf, centers, radii, candidates, K) then return true;
    }
    if nf.numUnits == 0 then return false;

    for u in 1..nf.numUnits {
      for e in 1..unitExponentRange {
        for sgn in (1, -1) {
          const up = unitPower(nf, u, sgn*e);
          const twisted = recenterByUnit(nf, box, centers, radii, up);
          if nf.isTotallyReal {
            if isBoxAbsorbedReal(nf, twisted, candidates, K) then return true;
          } else {
            var tc: [1..nf.numEmbeddings] complex(128);
            var tr: [1..nf.numEmbeddings] real(64);
            computeBoxProjections(nf, twisted, tc, tr);
            if isProjectionAbsorbed(nf, tc, tr, candidates, K) then return true;
          }
        }
      }
    }
    return false;
  }

  proc bisect(const ref box: Box(?d)): [0..#(1 << d)] Box(d) {
    const numChildren = 1 << d;
    var children: [0..#numChildren] Box(d);
    for bIdx in 0..#numChildren {
      var c: d*real(64);
      var w: d*real(64);
      for i in 1..d {
        w(i-1) = box.widths(i-1) * 0.5;
        const dir = if (bIdx & (1 << (i-1))) != 0 then 1.0 else -1.0;
        c(i-1) = box.center(i-1) + dir * w(i-1);
      }
      children[bIdx] = new Box(degree=d, center=c, widths=w);
    }
    return children;
  }

  record SieveResult {
    param degree: int;
    var cleared: bool;
    var depthReached: int;
    var numRemaining: int;
    var maxRemainingWidth: real(64);
    var remDom: domain(1) = {1..0};
    var remaining: [remDom] Box(degree);
  }

  private proc storeRemaining(ref res: SieveResult(?d), const ref boxes: [] Box(d)) {
    res.numRemaining = boxes.size;
    res.remDom = {1..boxes.size};
    var i = 1;
    for b in boxes {
      res.remaining[i] = b;
      i += 1;
    }
  }

  // Runs cut+absorb starting from the fundamental domain [-0.5,0.5]^degree
  // until every box is absorbed, or the resolution/depth budget is
  // exhausted. Returns the surviving ("uncleared") boxes, if any.
  // minWidth guards against a floating-point precision cliff: for a box
  // whose edge sits exactly on a genuine critical point (common, since
  // critical points are often dyadic in the integral-basis coordinates),
  // the true gap between its computed max-norm and K shrinks to zero
  // *linearly with width* as the box shrinks (not just "gets small" --
  // provably approaches exactly 0), so past some depth no floating-point
  // comparison (regardless of safety margin) can reliably tell "exactly at
  // the boundary" apart from "genuinely absorbable". Empirically this
  // starts to bite for widths below roughly 1e-7 to 1e-8; going deeper
  // than that does not yield genuinely higher-quality answers, only an
  // increasing risk of a false "cleared" result, so we stop there.
  proc runSieve(const ref nf: NumberFieldData, param degree: int, K: real(64), maxDepth: int = 30,
                minWidth: real(64) = 1.0e-7, boundRange: int = 0, verbose: bool = false,
                maxProblems: int = 200_000, useUnits: bool = false, unitExponentRange: int = 3,
                candidateCap: int = 2000): SieveResult(degree) {
    const candidates = smallElements(nf, K, boundRange=boundRange, candidateCap=candidateCap);
    if verbose then writeln("  [sieve] K=", K, " candidates=", candidates.n);
    var rootCenter: degree*real(64);
    var rootWidths: degree*real(64);
    for i in 0..degree-1 do rootWidths(i) = 0.5;
    const rootBox = new Box(degree=degree, center=rootCenter, widths=rootWidths);

    var currentDom: domain(1) = {0..#1};
    var current: [currentDom] Box(degree) = [rootBox];

    for depth in 1..maxDepth {
      if current.size == 0 {
        var res: SieveResult(degree);
        res.cleared = true; res.depthReached = depth-1; res.numRemaining = 0;
        res.maxRemainingWidth = 0.0;
        return res;
      }

      var nextList: list(Box(degree), parSafe=true);
      forall b in current with (ref nextList) {
        const absorbed = if useUnits then isBoxAbsorbedWithUnits(nf, b, candidates, K, unitExponentRange)
                                      else isBoxAbsorbed(nf, b, candidates, K);
        if !absorbed then nextList.pushBack(b);
      }

      if nextList.size == 0 {
        var res: SieveResult(degree);
        res.cleared = true; res.depthReached = depth; res.numRemaining = 0;
        res.maxRemainingWidth = 0.0;
        return res;
      }

      const next = nextList.toArray();
      if verbose then writeln("  [sieve] depth=", depth, " boxes-in=", current.size, " problems=", next.size);

      // Safety valve: if the problem count is exploding (a strong signal
      // that K is below the true minimum, so no amount of bisection will
      // ever clear this region), bail out now instead of continuing to
      // burn time/memory on an ever-growing tree (mirrors MAX_NUMBER_PB2
      // in the reference C implementation).
      if next.size > maxProblems {
        var res: SieveResult(degree);
        res.cleared = false; res.depthReached = depth;
        res.maxRemainingWidth = next[next.domain.low].widths(0);
        storeRemaining(res, next);
        return res;
      }

      // Stop refining once boxes are already far smaller than minWidth --
      // further bisection cannot change the verdict at working precision.
      const w0 = next[next.domain.low].widths(0);
      if w0 < minWidth || depth == maxDepth {
        var res: SieveResult(degree);
        res.cleared = false; res.depthReached = depth;
        res.maxRemainingWidth = w0;
        storeRemaining(res, next);
        return res;
      }

      var grownList: list(Box(degree), parSafe=true);
      forall b in next with (ref grownList) {
        for child in bisect(b) do grownList.pushBack(child);
      }
      const grown = grownList.toArray();
      currentDom = grown.domain;
      current = grown;
    }

    var res: SieveResult(degree);
    res.cleared = false; res.depthReached = maxDepth;
    if current.size > 0 then res.maxRemainingWidth = current[current.domain.low].widths(0);
    storeRemaining(res, current);
    return res;
  }
}
