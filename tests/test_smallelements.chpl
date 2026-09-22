use NumberField;
use SmallElements;

proc main() throws {
  const nf = loadNumberField("tests/fixtures/x2-2.txt");
  writeln("degree=", nf.degree, " boundRange(default)=", candidateBoundRange(nf.degree));
  const s = smallElements(nf, 1.0, boundRange = 5);
  writeln("num candidates = ", s.n, " isTotallyReal=", s.isTotallyReal);
  for i in 1..min(s.n, 10) {
    write("coeffs=[");
    for c in 1..s.degree do write(s.coeffs[i,c], " ");
    write("] emb=[");
    if s.isTotallyReal then for r in 1..s.numEmbeddings do write(s.embReal[i,r], " ");
    else for r in 1..s.numEmbeddings do write(s.emb[i,r], " ");
    writeln("]");
  }
}
