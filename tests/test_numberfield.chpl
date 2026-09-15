use NumberField;

proc main() throws {
  const nf = loadNumberField("tests/fixtures/x3-3x-1.txt");
  writeln("degree=", nf.degree, " r1=", nf.r1, " r2=", nf.r2, " numUnits=", nf.numUnits);
  writeln("basis=", nf.basis);
  writeln("units=", nf.unitEmbeddings);

  var coeffs: [1..nf.degree] real(64) = 0.5;
  const e = embed(nf, coeffs);
  writeln("embed(0.5,...)=", e);
  writeln("|N|=", absNorm(nf, e));
}
