use NumberField;
use Sieve;

proc main() throws {
  const nf = loadNumberField("tests/fixtures/x2-2.txt");
  writeln("Field x^2-2 (known Euclidean minimum = 0.5)");
  for K in [0.4, 0.5, 0.6, 1.0] {
    const res = runSieve(nf, K, maxDepth=20, boundRange=8);
    writeln("K=", K, " -> cleared=", res.cleared, " depth=", res.depthReached,
            " remaining=", res.numRemaining, " maxWidth=", res.maxRemainingWidth);
  }
}
