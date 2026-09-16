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
use Sieve;
use Certify;

config const fieldFile: string = "field_data.txt";
config const initialK: real(64) = 1.0;
config const tolerance: real(64) = 1.0e-4;
config const exploreDepth: int = 12;
config const refineDepth: int = 22;
config const boundRange: int = 0; // 0 = choose automatically from degree
config const maxBisectionSteps: int = 60;
config const maxProblems: int = 50_000;
config const useUnits: bool = true;
config const unitExponentRange: int = 1;
config const certify: bool = true;
config const certifySamples: int = 8;
config const verbose: bool = false;

proc main() throws {
  const nf = loadNumberField(fieldFile);
  writeln("Loaded field from ", fieldFile, ": degree=", nf.degree,
          " (r1,r2)=(", nf.r1, ",", nf.r2, ") fundamentalUnits=", nf.numUnits);

  // --- Phase A: bracket the minimum between a failing lo and a clearing hi ---
  var hi = initialK;
  var lo = 0.0;

  var up = 0;
  while !runSieve(nf, hi, maxDepth=exploreDepth, boundRange=boundRange, maxProblems=maxProblems,
                  useUnits=useUnits, unitExponentRange=unitExponentRange).cleared && up < 20 {
    lo = hi;
    hi *= 2.0;
    up += 1;
  }
  if verbose then writeln("[bracket] hi=", hi, " clears at exploreDepth=", exploreDepth);

  var down = 0;
  while down < 30 {
    const mid = hi / 2.0;
    if mid <= lo then break;
    const res = runSieve(nf, mid, maxDepth=exploreDepth, boundRange=boundRange, maxProblems=maxProblems,
                          useUnits=useUnits, unitExponentRange=unitExponentRange);
    if res.cleared then hi = mid; else { lo = mid; break; }
    down += 1;
  }
  writeln("Initial bracket: ", lo, " < M(K) <= ", hi);

  // --- Phase B: bisection refinement with a deeper resolution budget ---
  var step = 0;
  var lastBlocked: SieveResult;
  var haveBlocked = false;
  while (hi - lo) > tolerance && step < maxBisectionSteps {
    const mid = (lo + hi) / 2.0;
    const res = runSieve(nf, mid, maxDepth=refineDepth, boundRange=boundRange, verbose=verbose, maxProblems=maxProblems,
                          useUnits=useUnits, unitExponentRange=unitExponentRange);
    if res.cleared {
      hi = mid;
      writeln("K=", mid, " -> cleared (new upper bound)");
    } else {
      lo = mid;
      lastBlocked = res;
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
    var bestNum = 0.0;
    var bestStr = "";
    const n = min(lastBlocked.numRemaining, certifySamples);
    for i in 1..n {
      try {
        const exact = exactMinimalNormAt(nf.polynomial, lastBlocked.remaining[i].center, nf.degree);
        const parts = exact.split("/");
        const val = if parts.size == 2 then parts[0]:real(64) / parts[1]:real(64) else exact:real(64);
        // Soundness guard: m_K(x) <= M(K) < hi always, since hi is a proven
        // upper bound. A sample >= hi is therefore provably impossible and
        // means exactMinimalNormAt's bounded search missed the true nearest
        // lattice point (this happens for large-regulator fields, since the
        // exact search does not yet use unit-action reduction the way the
        // numeric sieve does) -- discard it rather than reporting a false
        // "certified" claim.
        if val > bestNum && val < hi { bestNum = val; bestStr = exact; }
      } catch {
        // gp unavailable or parsing failed for this sample; skip it.
      }
    }
    if bestStr.size > 0 {
      writeln("  certified lower bound:   M(K) >= ", bestStr, " ~= ", bestNum, " (exact, proven via Pari)");
    } else {
      writeln("  certification unavailable this run (no sample survived the ",
              "soundness check against the proven upper bound -- try a larger ",
              "--certifySamples, or this field's regulator may require ",
              "unit-aware exact search, not yet implemented)");
    }
  }
}
