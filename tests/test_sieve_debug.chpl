use NumberField;
use Sieve;

proc main() throws {
  const nf = loadNumberField("tests/fixtures/x2-2.txt");
  writeln("Field x^2-2 (known Euclidean minimum = 0.5), K=2.0 (should clear fast)");
  const res = runSieve(nf, 2, 2.0, maxDepth=8, boundRange=8, verbose=true);
  writeln("cleared=", res.cleared, " depth=", res.depthReached, " remaining=", res.numRemaining);
}
