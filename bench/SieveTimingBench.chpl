// SieveTimingBench.chpl
//
// Isolates a single runSieve() call's cost, bypassing main.chpl's outer
// bracket-search loop, so candidateCap can be varied directly to test the
// hypothesis that the per-candidate absorption scan (O(candidateCap) per
// box) dominates over box construction/copying (O(degree) per box, but
// with per-object allocation overhead for the array-based Box) for
// default settings -- and that shrinking candidateCap should let the
// tuple-Box win from bench/BoxCopyBench.chpl show up end-to-end.
//
// Usage (compile once against the OLD Sieve.chpl/main.chpl-less sources
// via git stash, and again against the NEW ones, to compare):
//   chpl --fast -M . bench/SieveTimingBench.chpl -o bench/sieve_timing_old
//   ./bench/sieve_timing_old --fieldFile=q57.txt --K=0.5 --candidateCap=2000
//   ./bench/sieve_timing_old --fieldFile=q57.txt --K=0.5 --candidateCap=20
use NumberField;
use Sieve;
use Time;

config const fieldFile: string = "q57.txt";
config const K: real(64) = 0.5; // deliberately below the true minimum for
                                 // most q*.txt fields, so the sieve never
                                 // clears and does sustained, comparable work
config const candidateCap: int = 2000;
config const maxDepth: int = 24;
config const maxProblems: int = 2_000_000;

proc main() throws {
  const nf = loadNumberField(fieldFile);
  writeln("field=", fieldFile, " degree=", nf.degree, " K=", K,
          " candidateCap=", candidateCap, " maxDepth=", maxDepth);
  var sw: stopwatch;
  sw.start();
  select nf.degree {
    when 2 {
      const res = runSieve(nf, 2, K, maxDepth=maxDepth, useUnits=true, unitExponentRange=1,
                            maxProblems=maxProblems, candidateCap=candidateCap);
      sw.stop();
      writeln("cleared=", res.cleared, " depthReached=", res.depthReached,
              " numRemaining=", res.numRemaining, " elapsed=", sw.elapsed(), "s");
    }
    when 3 {
      const res = runSieve(nf, 3, K, maxDepth=maxDepth, useUnits=true, unitExponentRange=1,
                            maxProblems=maxProblems, candidateCap=candidateCap);
      sw.stop();
      writeln("cleared=", res.cleared, " depthReached=", res.depthReached,
              " numRemaining=", res.numRemaining, " elapsed=", sw.elapsed(), "s");
    }
    when 4 {
      const res = runSieve(nf, 4, K, maxDepth=maxDepth, useUnits=true, unitExponentRange=1,
                            maxProblems=maxProblems, candidateCap=candidateCap);
      sw.stop();
      writeln("cleared=", res.cleared, " depthReached=", res.depthReached,
              " numRemaining=", res.numRemaining, " elapsed=", sw.elapsed(), "s");
    }
    when 6 {
      const res = runSieve(nf, 6, K, maxDepth=maxDepth, useUnits=true, unitExponentRange=1,
                            maxProblems=maxProblems, candidateCap=candidateCap);
      sw.stop();
      writeln("cleared=", res.cleared, " depthReached=", res.depthReached,
              " numRemaining=", res.numRemaining, " elapsed=", sw.elapsed(), "s");
    }
    otherwise halt("bench only wired for degree 2-4, 6");
  }
}
