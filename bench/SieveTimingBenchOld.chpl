// SieveTimingBenchOld.chpl -- see SieveTimingBench.chpl for the rationale.
// Matches the OLD (pre-Phase-2, array-based Box) runSieve signature, which
// takes no explicit `param degree` argument.
use NumberField;
use Sieve;
use Time;

config const fieldFile: string = "q57.txt";
config const K: real(64) = 0.5;
config const candidateCap: int = 2000;
config const maxDepth: int = 24;
config const maxProblems: int = 2_000_000;

proc main() throws {
  const nf = loadNumberField(fieldFile);
  writeln("field=", fieldFile, " degree=", nf.degree, " K=", K,
          " candidateCap=", candidateCap, " maxDepth=", maxDepth);
  var sw: stopwatch;
  sw.start();
  const res = runSieve(nf, K, maxDepth=maxDepth, useUnits=true, unitExponentRange=1,
                        maxProblems=maxProblems, candidateCap=candidateCap);
  sw.stop();
  writeln("cleared=", res.cleared, " depthReached=", res.depthReached,
          " numRemaining=", res.numRemaining, " elapsed=", sw.elapsed(), "s");
}
