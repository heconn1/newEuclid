// TupleVsArrayBench.chpl
//
// Standalone microbenchmark isolating ONE question: for a fixed, small
// compile-time-known degree D, how much does it cost to store per-box /
// per-candidate coordinate vectors in domain-backed arrays (as Box and
// SmallElementSet do today) versus fixed-size D*real(64) tuples?
//
// The inner loop shape mirrors Sieve.chpl's isProjectionAbsorbed: for each
// of `numBoxes` box-like centers, scan `cap` candidate vectors and reduce
// over D coordinates per candidate (sum of squared differences, standing
// in for the abs()/multiply chain in calculateBoxMaxNorm).
//
// Usage:
//   chpl --fast -sD=2 TupleVsArrayBench.chpl -o bench_d2
//   ./bench_d2 --cap=2000 --numBoxes=200000
//
// Run once per degree of interest (D is a compile-time param, matching the
// "recompile if we ever need degree > 8" tradeoff discussed for tuples).
use Time;
use Random;

config param D = 4;       // compile-time degree; recompile to change
config const cap = 2000;      // candidates per box, matches candidateCap default
config const numBoxes = 200_000; // outer "boxes examined" count
config const seed = 42;

proc main() {
  writeln("D=", D, " cap=", cap, " numBoxes=", numBoxes);

  // ---- shared random candidate data, generated once ----
  var candFlat: [0..#(cap*D)] real(64);
  fillRandom(candFlat, seed);

  // ================= Array-backed version (today's approach) =================
  {
    var candArr: [1..cap, 1..D] real(64);
    for i in 1..cap do for d in 1..D do candArr[i, d] = candFlat[(i-1)*D + (d-1)];

    var total = 0.0;
    var sw: stopwatch;
    sw.start();
    forall b in 0..#numBoxes with (+ reduce total) {
      var center: [1..D] real(64);
      for d in 1..D do center[d] = ((b*7 + d*3) % 97): real(64) * 0.01;

      var best = max(real(64));
      for i in 1..cap {
        var s = 0.0;
        for d in 1..D {
          const diff = center[d] - candArr[i, d];
          s += diff * diff;
        }
        if s < best then best = s;
      }
      total += best;
    }
    sw.stop();
    writeln("array : ", sw.elapsed(), " s   (checksum=", total, ")");
  }

  // ================= Tuple-backed version (proposed change) =================
  {
    var candTup: [1..cap] D*real(64);
    for i in 1..cap {
      var t: D*real(64);
      for param d in 0..D-1 do t(d) = candFlat[(i-1)*D + d];
      candTup[i] = t;
    }

    var total = 0.0;
    var sw: stopwatch;
    sw.start();
    forall b in 0..#numBoxes with (+ reduce total) {
      var center: D*real(64);
      for param d in 0..D-1 do center(d) = ((b*7 + (d+1)*3) % 97): real(64) * 0.01;

      var best = max(real(64));
      for i in 1..cap {
        const cand = candTup[i];
        var s = 0.0;
        for param d in 0..D-1 {
          const diff = center(d) - cand(d);
          s += diff * diff;
        }
        if s < best then best = s;
      }
      total += best;
    }
    sw.stop();
    writeln("tuple : ", sw.elapsed(), " s   (checksum=", total, ")");
  }
}
