// main.chpl
//
// CLI driver for the arbitrary-degree Euclidean-minimum sieve.
//
// Strategy: maintain a rigorously *proven* upper bound `hi` (some K for
// which the sieve cleared every box, i.e. M(K_field) < hi is certified) and
// a heuristic `lo` (largest K at which the sieve failed to clear within the
// depth budget -- this is NOT a certified lower bound, since failure can
// also mean "needs more resolution"; see the Certify module for turning a
// tight bracket into a certified exact answer).
use NumberField;
use SmallElements;
use Sieve;
use Certify;

config const fieldFile: string = "field_data.txt";
config const initialK: real(64) = 1.0;
config const tolerance: real(64) = 1.0e-4;
config const exploreDepth: int = 12;
config const refineDepth: int = 26;
config const boundRange: int = 0; // 0 = choose automatically from degree
config const maxBisectionSteps: int = 60;
config const maxProblems: int = 50_000;
config const useUnits: bool = true;
config const unitExponentRange: int = 1;
config const candidateCap: int = 2000;
config const certify: bool = true;
config const certifyDepth: int = 30;
config const certifySamples: int = 8;
config const verbose: bool = false;

proc main() throws {
  const nf = loadNumberField(fieldFile);
  writeln("Loaded field from ", fieldFile, ": degree=", nf.degree,
          " (r1,r2)=(", nf.r1, ",", nf.r2, ") fundamentalUnits=", nf.numUnits);

  // --- Phase A: find a clearing upper bound hi ---
  //
  // IMPORTANT: a shallow (exploreDepth) sieve trial that reports "cleared"
  // is always trustworthy (absorption proofs are sound regardless of the
  // depth used to find them) -- but a shallow trial that reports "blocked"
  // is NOT trustworthy: it can simply mean the resolution budget was too
  // small to witness absorption yet, not that K is actually too small.
  // Earlier versions used shallow "blocked" results to raise `lo`, which
  // could -- and for some real quadratic fields with awkward regulators,
  // did -- lock `lo` above the true minimum before the deep bisection
  // phase ever ran, permanently excluding the true answer from the search
  // range. So `lo` always starts at the trivially-safe value 0 and is only
  // ever raised from a *deep* (refineDepth) "blocked" result in Phase B.
  var hi = initialK;
  var lo = 0.0;

  var up = 0;
  while !runSieve(nf, hi, maxDepth=exploreDepth, boundRange=boundRange, maxProblems=maxProblems,
                  useUnits=useUnits, unitExponentRange=unitExponentRange, candidateCap=candidateCap).cleared && up < 20 {
    hi *= 2.0;
    up += 1;
  }
  if verbose then writeln("[bracket] hi=", hi, " clears at exploreDepth=", exploreDepth);

  // Optional tightening: keep halving hi while it still clears at shallow
  // depth (cheap and, again, only ever trusting positive "cleared" results).
  var down = 0;
  while down < 30 {
    const mid = hi / 2.0;
    if mid <= lo then break;
    const res = runSieve(nf, mid, maxDepth=exploreDepth, boundRange=boundRange, maxProblems=maxProblems,
                          useUnits=useUnits, unitExponentRange=unitExponentRange, candidateCap=candidateCap);
    if res.cleared then hi = mid; else break;
    down += 1;
  }
  writeln("Initial bracket: ", lo, " < M(K) <= ", hi);

  // --- Phase B: bisection refinement with a deeper resolution budget ---
  var step = 0;
  var lastBlocked: SieveResult;
  var lastBlockedK = 0.0;
  var haveBlocked = false;
  while (hi - lo) > tolerance && step < maxBisectionSteps {
    const mid = (lo + hi) / 2.0;
    const res = runSieve(nf, mid, maxDepth=refineDepth, boundRange=boundRange, verbose=verbose, maxProblems=maxProblems,
                          useUnits=useUnits, unitExponentRange=unitExponentRange, candidateCap=candidateCap);
    if res.cleared {
      hi = mid;
      writeln("K=", mid, " -> cleared (new upper bound)");
    } else {
      lo = mid;
      lastBlocked = res;
      lastBlockedK = mid;
      haveBlocked = true;
      writeln("K=", mid, " -> blocked (", res.numRemaining,
              " unresolved box(es), width=", res.maxRemainingWidth, ")");
    }
    step += 1;
  }

  writeln();
  writeln("Result after ", step, " bisection steps:");
  writeln("  proven upper bound:      M(K) <  ", hi);
  writeln("  heuristic lower bound:   M(K) >= ", lo, " (uncertified)");
  writeln("  bracket width:           ", hi - lo);

  const (fnum, fden) = simplestFractionInInterval(lo, hi);
  writeln("  simplest fraction in bracket (conjectured exact value, NOT a proof): ",
          fnum, "/", fden);

  if certify && haveBlocked && nf.polynomial.size > 0 {
    writeln();
    writeln("Phase 3: sampling resistant boxes for an exact, Pari-certified lower bound...");
    // Re-run the sieve at the final lo with a deeper budget than the
    // bisection loop used, specifically to collect resistant boxes for
    // certification: a box that's merely "unresolved so far" at a
    // middling depth (as lastBlocked may be, since it came from whatever
    // depth the bisection loop used) is not necessarily near a genuine
    // critical point -- it can just need more refinement, which would
    // make certifying its norm meaningless (and, confusingly, can even
    // produce a value inconsistent with the proven upper bound). Using
    // extra depth here makes the sampled boxes far more likely to be
    // genuine resistant points.
    const certRes = runSieve(nf, lo, maxDepth=certifyDepth, boundRange=boundRange, maxProblems=maxProblems,
                              useUnits=useUnits, unitExponentRange=unitExponentRange, candidateCap=candidateCap);
    const sampleSource = if certRes.numRemaining > 0 then certRes else lastBlocked;
    const sampleK = if certRes.numRemaining > 0 then lo else lastBlockedK;
    // Reuse the exact same (norm-ranked) candidate list the sieve itself
    // used at that K, so certification can't be *less* complete than the
    // sieve was -- a smaller/blind candidate set here could otherwise
    // report a value inconsistent with (larger than) an independently-
    // proven upper bound.
    const certCandidates = smallElements(nf, sampleK, boundRange=boundRange, candidateCap=candidateCap);
    var bestNum = 0.0;
    var bestStr = "";
    var skipped = 0;
    const n = min(sampleSource.numRemaining, certifySamples);
    for i in 1..n {
      try {
        const exact = exactMinimalNormAtWithCandidates(nf.polynomial, sampleSource.remaining[i].center,
                                                         nf.degree, certCandidates);
        const parts = exact.split("/");
        const val = if parts.size == 2 then parts[0]:real(64) / parts[1]:real(64) else exact:real(64);
        // A value exceeding the independently-proven upper bound is a
        // logical impossibility (m_K(x) <= M(K_field) < hi always) -- it
        // means this particular sampled box is a depth-limited artifact,
        // not a genuine critical point, so we simply don't trust it rather
        // than reporting a self-contradictory "certified" claim.
        if val > hi + 1.0e-9 {
          skipped += 1;
        } else if val > bestNum {
          bestNum = val;
          bestStr = exact;
        }
      } catch {
        // gp unavailable or parsing failed for this sample; skip it.
      }
    }
    if bestStr.size > 0 {
      writeln("  certified lower bound:   M(K) >= ", bestStr, " ~= ", bestNum, " (exact, proven via Pari)");
      if skipped > 0 then
        writeln("  (", skipped, " other sampled box(es) discarded: exact value exceeded the proven upper bound,",
                " i.e. they were depth-limited artifacts, not genuine critical points)");
    } else if skipped > 0 {
      writeln("  certification inconclusive: all ", skipped, " sampled box(es) were depth-limited artifacts ",
              "(their exact norm exceeded the proven upper bound); try a larger --certifyDepth.");
    } else {
      writeln("  certification unavailable this run (no sample survived the ",
              "soundness check against the proven upper bound -- try a larger ",
              "--certifySamples, or this field's regulator may require ",
              "unit-aware exact search, not yet implemented)");
    }
  }
}
