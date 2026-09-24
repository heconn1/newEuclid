use NumberField;
use Sieve;
use Certify;

proc main() throws {
  const nf = loadNumberField("tests/fixtures/x3-3x-1.txt");
  // K=0.332813 was found to be "blocked" (numeric upper bound proof
  // fails) just below the true minimum 1/3 -- extract a resistant box.
  const res = runSieve(nf, 3, 0.332813, maxDepth=20, useUnits=true, unitExponentRange=1, maxProblems=50000);
  writeln("cleared=", res.cleared, " remaining=", res.numRemaining);
  if !res.cleared && res.numRemaining > 0 {
    const numToScan = min(res.numRemaining, 8);
    for i in 1..numToScan {
      const box = res.remaining[i];
      const exact = exactMinimalNormAt("x^3-3*x-1", boxCenterToArray(box), nf.degree, searchRange=3);
      writeln("box[", i, "] center=", box.center, " -> exact minimal norm = ", exact);
    }
    writeln("(expected exact minimum: 1/3)");
  }
}
