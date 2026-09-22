// BoxCopyBench.chpl
//
// Targets a different hypothesis than TupleVsArrayBench.chpl: not raw
// element-access cost, but the cost of *constructing and copying* Box
// records themselves, which is what Sieve.chpl's bisect() + list-collect
// pattern actually does at every depth level:
//
//   forall b in current with (ref nextList) { ... nextList.pushBack(b); }
//   ...
//   forall b in next with (ref grownList) {
//     for child in bisect(b) do grownList.pushBack(child);   // 2^degree new Box values
//   }
//   const grown = grownList.toArray();
//
// Each Box (today) carries a domain field plus two arrays over that
// domain. Constructing/copying such a record may require heap allocation
// + array-descriptor bookkeeping per copy. A tuple-based Box has no domain
// and no heap allocation at all -- copying it is a flat memcpy of stack
// bytes. This benchmark simulates N generations of "bisect + collect" for
// a fixed degree D and measures wall-clock time for both representations.
//
// Usage:
//   chpl --fast -sD=4 BoxCopyBench.chpl -o boxcopy_d4
//   ./boxcopy_d4 --startBoxes=1 --generations=6
use Time;
use List;

config param D = 4;
config const startBoxes = 64;   // boxes present at generation 0
config const generations = 5;   // how many bisect+collect rounds to simulate
config const reps = 3;          // repeat the whole experiment this many times

record BoxArr {
  var dom: domain(1) = {1..D};
  var center: [dom] real(64);
  var widths: [dom] real(64);
}

record BoxTup {
  var center: D*real(64);
  var widths: D*real(64);
}

proc bisectArr(const ref b: BoxArr): [0..#(1<<D)] BoxArr {
  const numChildren = 1 << D;
  var children: [0..#numChildren] BoxArr;
  for bIdx in 0..#numChildren {
    var c: [b.dom] real(64);
    var w: [b.dom] real(64);
    for d in 1..D {
      w[d] = b.widths[d] * 0.5;
      const dir = if (bIdx & (1 << (d-1))) != 0 then 1.0 else -1.0;
      c[d] = b.center[d] + dir * w[d];
    }
    children[bIdx] = new BoxArr(dom=b.dom, center=c, widths=w);
  }
  return children;
}

proc bisectTup(const ref b: BoxTup): [0..#(1<<D)] BoxTup {
  const numChildren = 1 << D;
  var children: [0..#numChildren] BoxTup;
  for bIdx in 0..#numChildren {
    var c: D*real(64);
    var w: D*real(64);
    for param d in 0..D-1 {
      w(d) = b.widths(d) * 0.5;
      const dir = if (bIdx & (1 << d)) != 0 then 1.0 else -1.0;
      c(d) = b.center(d) + dir * w(d);
    }
    children[bIdx] = new BoxTup(center=c, widths=w);
  }
  return children;
}

proc runArr() {
  var currentDom: domain(1) = {0..#startBoxes};
  var current: [currentDom] BoxArr;
  for i in current.domain do current[i] = new BoxArr();

  var sw: stopwatch;
  sw.start();
  for gen in 1..generations {
    var nextList: list(BoxArr, parSafe=true);
    forall b in current with (ref nextList) {
      for child in bisectArr(b) do nextList.pushBack(child);
    }
    const grown = nextList.toArray();
    currentDom = grown.domain;
    current = grown;
  }
  sw.stop();
  writeln("array : ", sw.elapsed(), " s   (finalCount=", current.size, ")");
}

proc runTup() {
  var currentDom: domain(1) = {0..#startBoxes};
  var current: [currentDom] BoxTup;
  for i in current.domain do current[i] = new BoxTup();

  var sw: stopwatch;
  sw.start();
  for gen in 1..generations {
    var nextList: list(BoxTup, parSafe=true);
    forall b in current with (ref nextList) {
      for child in bisectTup(b) do nextList.pushBack(child);
    }
    const grown = nextList.toArray();
    currentDom = grown.domain;
    current = grown;
  }
  sw.stop();
  writeln("tuple : ", sw.elapsed(), " s   (finalCount=", current.size, ")");
}

proc main() {
  writeln("D=", D, " startBoxes=", startBoxes, " generations=", generations,
          " (final generation would hold ", startBoxes * (1<<D)**generations, " boxes)");
  for r in 1..reps {
    writeln("-- rep ", r, " --");
    runArr();
    runTup();
  }
}
