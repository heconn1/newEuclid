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
module Sieve {
  use NumberField;
  use SmallElements;
  use List;
  use Math;

  record Box {
    var degree: int;
    var centerDom: domain(1);
    var center: [centerDom] real(64);
    var widths: [centerDom] real(64);
  }

  proc computeBoxProjections(const ref nf: NumberFieldData, const ref box: Box,
                              ref centers: [] complex(128), ref radii: [] real(64)) {
    forall row in 1..nf.numEmbeddings {
      var c: complex(128) = 0.0;
      var r = 0.0;
      for col in 1..nf.degree {
        c += nf.basis[row, col] * box.center[col];
        r += abs(nf.basis[row, col]) * box.widths[col];
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
  proc isBoxAbsorbed(const ref nf: NumberFieldData, const ref box: Box,
                      const ref candidates: SmallElementSet, K: real(64)): bool {
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

  proc recenterByUnit(const ref nf: NumberFieldData, const ref box: Box,
                       const ref centers: [] complex(128), const ref radii: [] real(64),
                       const ref unitPow: [] complex(128)): Box {
    var tCenters: [1..nf.numEmbeddings] complex(128);
    var tRadii: [1..nf.numEmbeddings] real(64);
    forall row in 1..nf.numEmbeddings {
      tCenters[row] = centers[row] * unitPow[row];
      tRadii[row] = radii[row] * abs(unitPow[row]);
    }

    const w = coefficientsFromEmbedding(nf, tCenters);
    var t: [1..nf.degree] real(64);
    for i in 1..nf.degree do t[i] = round(w[i]);

    var newCenter: [box.centerDom] real(64);
    for i in 1..nf.degree do newCenter[i] = w[i] - t[i];

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
    var newWidths: [box.centerDom] real(64);
    forall c in 1..nf.degree {
      var s = 0.0;
      for k in 1..nf.degree do s += abs(nf.sigmaInv[c, k]) * unpackedRadius[k];
      newWidths[c] = s;
    }

    return new Box(degree=nf.degree, centerDom=box.centerDom, center=newCenter, widths=newWidths);
  }

  // Direct absorption, then (if that fails) a sweep of fundamental unit
  // powers with recentering. `unitExponentRange` controls how many powers
  // of each fundamental unit are tried (both signs); this needs to be
  // large enough to reach whatever power of the unit is relevant for a
  // given field's regulator (see runSieve for how the default scales with
  // the field's own unit magnitudes).
  proc isBoxAbsorbedWithUnits(const ref nf: NumberFieldData, const ref box: Box,
                               const ref candidates: SmallElementSet, K: real(64),
                               unitExponentRange: int = 1): bool {
    var centers: [1..nf.numEmbeddings] complex(128);
    var radii: [1..nf.numEmbeddings] real(64);
    computeBoxProjections(nf, box, centers, radii);
    if isProjectionAbsorbed(nf, centers, radii, candidates, K) then return true;
    if nf.numUnits == 0 then return false;

    for u in 1..nf.numUnits {
      for e in 1..unitExponentRange {
        for sgn in (1, -1) {
          const up = unitPower(nf, u, sgn*e);
          const twisted = recenterByUnit(nf, box, centers, radii, up);
          var tc: [1..nf.numEmbeddings] complex(128);
          var tr: [1..nf.numEmbeddings] real(64);
          computeBoxProjections(nf, twisted, tc, tr);
          if isProjectionAbsorbed(nf, tc, tr, candidates, K) then return true;
        }
      }
    }
    return false;
  }

  proc bisect(const ref box: Box): [0..#(1 << box.degree)] Box {
    const degree = box.degree;
    const numChildren = 1 << degree;
    var children: [0..#numChildren] Box;
    for bIdx in 0..#numChildren {
      var c: [box.centerDom] real(64);
      var w: [box.centerDom] real(64);
      for d in 1..degree {
        w[d] = box.widths[d] * 0.5;
        const dir = if (bIdx & (1 << (d-1))) != 0 then 1.0 else -1.0;
        c[d] = box.center[d] + dir * w[d];
      }
      children[bIdx] = new Box(degree=degree, centerDom=box.centerDom, center=c, widths=w);
    }
    return children;
  }

  record SieveResult {
    var cleared: bool;
    var depthReached: int;
    var numRemaining: int;
    var maxRemainingWidth: real(64);
    var remDom: domain(1) = {1..0};
    var remaining: [remDom] Box;
  }

  private proc storeRemaining(ref res: SieveResult, const ref boxes: [] Box) {
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
  proc runSieve(const ref nf: NumberFieldData, K: real(64), maxDepth: int = 30,
                minWidth: real(64) = 1.0e-7, boundRange: int = 0, verbose: bool = false,
                maxProblems: int = 200_000, useUnits: bool = false, unitExponentRange: int = 3,
                candidateCap: int = 2000): SieveResult {
    const candidates = smallElements(nf, K, boundRange=boundRange, candidateCap=candidateCap);
    if verbose then writeln("  [sieve] K=", K, " candidates=", candidates.n);
    const cDom = {1..nf.degree};
    var rootCenter: [cDom] real(64) = 0.0;
    var rootWidths: [cDom] real(64) = 0.5;
    const rootBox = new Box(degree=nf.degree, centerDom=cDom, center=rootCenter, widths=rootWidths);

    var currentDom: domain(1) = {0..#1};
    var current: [currentDom] Box = [rootBox];

    for depth in 1..maxDepth {
      if current.size == 0 {
        var res: SieveResult;
        res.cleared = true; res.depthReached = depth-1; res.numRemaining = 0;
        res.maxRemainingWidth = 0.0;
        return res;
      }

      var nextList: list(Box, parSafe=true);
      forall b in current with (ref nextList) {
        const absorbed = if useUnits then isBoxAbsorbedWithUnits(nf, b, candidates, K, unitExponentRange)
                                      else isBoxAbsorbed(nf, b, candidates, K);
        if !absorbed then nextList.pushBack(b);
      }

      if nextList.size == 0 {
        var res: SieveResult;
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
        var res: SieveResult;
        res.cleared = false; res.depthReached = depth;
        res.maxRemainingWidth = next[next.domain.low].widths[1];
        storeRemaining(res, next);
        return res;
      }

      // Stop refining once boxes are already far smaller than minWidth --
      // further bisection cannot change the verdict at working precision.
      const w0 = next[next.domain.low].widths[1];
      if w0 < minWidth || depth == maxDepth {
        var res: SieveResult;
        res.cleared = false; res.depthReached = depth;
        res.maxRemainingWidth = w0;
        storeRemaining(res, next);
        return res;
      }

      var grownList: list(Box, parSafe=true);
      forall b in next with (ref grownList) {
        for child in bisect(b) do grownList.pushBack(child);
      }
      const grown = grownList.toArray();
      currentDom = grown.domain;
      current = grown;
    }

    var res: SieveResult;
    res.cleared = false; res.depthReached = maxDepth;
    if current.size > 0 then res.maxRemainingWidth = current[current.domain.low].widths[1];
    storeRemaining(res, current);
    return res;
  }
}
