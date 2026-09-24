// GenericBoxSyntax.chpl
//
// Throwaway syntax check for Phase 2: confirms the exact construction and
// generic-dispatch syntax for a param-degree record wrapping fixed-size
// tuples, before rewriting Sieve.chpl's real Box record around it.
record Box {
  param degree: int;
  var center: degree*real(64);
  var widths: degree*real(64);
}

// Should infer d from the argument -- no explicit param needed here.
proc bisectFirstChild(const ref box: Box(?d)): Box(d) {
  var c: d*real(64);
  var w: d*real(64);
  for i in 0..d-1 {
    w(i) = box.widths(i) * 0.5;
    c(i) = box.center(i) - w(i);
  }
  return new Box(degree=d, center=c, widths=w);
}

// The one place that DOES need an explicit param, since there's no Box
// argument to infer it from (mirrors runSieve's situation).
proc makeRootBox(param degree: int): Box(degree) {
  var c: degree*real(64);
  var w: degree*real(64);
  for i in 0..degree-1 do w(i) = 0.5;
  return new Box(degree=degree, center=c, widths=w);
}

proc main() {
  const b2 = makeRootBox(2);
  writeln("b2.degree=", b2.degree, " center=", b2.center, " widths=", b2.widths);
  const child = bisectFirstChild(b2);
  writeln("child.center=", child.center, " child.widths=", child.widths);

  const b4 = makeRootBox(4);
  writeln("b4.degree=", b4.degree, " center=", b4.center);

  // runtime-to-compile-time dispatch, mirroring main.chpl's planned `select`
  const runtimeDegree = 3;
  select runtimeDegree {
    when 1 do writeln("dispatched degree=1: ", makeRootBox(1));
    when 2 do writeln("dispatched degree=2: ", makeRootBox(2));
    when 3 do writeln("dispatched degree=3: ", makeRootBox(3));
    otherwise halt("unsupported degree");
  }
}
