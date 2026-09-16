// Certify.chpl
//
// Phase 3: turn a numeric upper bound into a certified exact answer.
//
// The sieve (Sieve.chpl) only ever proves upper bounds on M(K): "K clears"
// means M(K_field) < K, proven in floating point. It can never certify a
// matching lower bound by itself. But every surviving "problem" box has a
// rational center (some k/2^depth in each coefficient), and at that exact
// rational point x, computing min_{gamma in O_K} |N(x - gamma)| exactly
// (via Pari, which has no floating-point error) gives a genuine, certified
// LOWER bound on M(K_field): trivially m_K(x) <= M(K_field), so if we
// exhibit an x with exact m_K(x) = v, then M(K_field) >= v.
//
// If that exact value matches the numeric upper bound from the sieve (as
// it typically will for the well-known small examples), the Euclidean
// minimum is certified exactly -- without ever porting Lezowski's
// unit-orbit cycle/graph machinery (src/graph.c).
//
// This performs exactly one Pari/gp subprocess call per certification
// attempt, keeping the "one-time Pari work" principle: Pari is used only
// to finish off a numerically-narrowed candidate, never for the bulk
// numerical search.
module Certify {
  use Subprocess;
  use IO;
  use Math;

  // Finds the fraction with the smallest denominator lying strictly
  // inside (lo, hi) (Stern-Brocot / continued-fraction mediant search).
  // Used only to *suggest* a human-readable closed form (e.g. "1/3") from
  // a tight numeric bracket; it is not itself a proof.
  proc simplestFractionInInterval(lo: real(64), hi: real(64)): (int, int) {
    if lo >= hi then return (0, 1);
    if lo <= 0.0 && hi >= 0.0 then return (0, 1);
    if lo < 0.0 then {
      const (n, d) = simplestFractionInInterval(-hi, -lo);
      return (-n, d);
    }
    var aNum = 0, aDen = 1;   // 0/1
    var bNum = 1, bDen = 0;   // 1/0 (infinity)
    for _i in 1..200 {
      const mNum = aNum + bNum;
      const mDen = aDen + bDen;
      const mVal = mNum: real(64) / mDen: real(64);
      if mVal <= lo then { aNum = mNum; aDen = mDen; }
      else if mVal >= hi then { bNum = mNum; bDen = mDen; }
      else return (mNum, mDen);
    }
    return (aNum, aDen);
  }

  // Chooses the largest symmetric per-coordinate integer range R such that
  // enumerating (2R+1)^degree exact candidates in Pari stays within
  // `budget` -- a flat small range (e.g. 3) is enough for many fields but
  // demonstrably too small for others (e.g. x^2-61 needs R>=5 to recover
  // its true exact minimum instead of an inflated, unsound value).
  proc defaultSearchRange(degree: int, budget: real = 5000.0): int {
    if degree <= 0 then return 2;
    const raw = (budget ** (1.0/degree) - 1.0) / 2.0;
    var r = floor(raw): int;
    if r < 2 then r = 2;
    if r > 60 then r = 60;
    return r;
  }

  // Returns (numerator, denominator) with denominator a power of two,
  // exact for any dyadic rational produced by our bisection (box centers
  // are always k/2^depth).
  proc toExactDyadic(x: real(64), maxBits: int = 52): (int, int) {
    for d in 0..maxBits {
      const den = 1:int << d;
      const num = round(x * den:real(64)): int;
      if abs(num:real(64) / den:real(64) - x) < 1.0e-14 then return (num, den);
    }
    const den = 1:int << maxBits;
    return (round(x * den:real(64)): int, den);
  }

  // Computes min_{gamma in O_K, u a small unit power} |N(u*x - gamma)|
  // exactly, where x has the given (rational) integral-basis coefficients,
  // by shelling out once to gp. Since |N(u)| = 1 for any unit u,
  // |N(u*x - gamma)| = |N(x - gamma*u^-1)|, so searching over small unit
  // powers alongside small gamma is exactly the exact-arithmetic analogue
  // of Sieve.chpl's unit-action acceleration: without it, this search can
  // miss the true nearest lattice point for large-regulator fields (the
  // absorbers relevant there are spread out along the unit orbit, not
  // clustered near the origin) and silently report an inflated value.
  // Returns gp's printed value (an exact rational, e.g. "1/3") on success.
  proc exactMinimalNormAt(polynomial: string, const ref coeffs: [] real(64), degree: int,
                           searchRange: int = 0, unitExponentRange: int = 1): string throws {
    var terms: string;
    for i in 1..degree {
      const (num, den) = toExactDyadic(coeffs[i]);
      if num != 0 {
        if terms.size > 0 then terms += " + ";
        terms += "(" + num:string + "/" + den:string + ")*zk[" + i:string + "]";
      }
    }
    if terms.size == 0 then terms = "0";

    const R = (if searchRange > 0 then searchRange else defaultSearchRange(degree)): string;
    const E = unitExponentRange: string;
    var script: string;
    script += "P = " + polynomial + ";\n";
    script += "nf = bnfinit(P);\n";
    script += "n = poldegree(P);\n";
    script += "zk = nf.zk;\n";
    script += "fu = nf.fu;\n";
    script += "r = length(fu);\n";
    script += "x = " + terms + ";\n";
    script += "best = -1;\n";
    script += "forvec(e=vector(r,i,[-" + E + "," + E + "]), u=1; for(i=1,r,u=u*fu[i]^e[i]); " +
              "forvec(v=vector(n,i,[-" + R + "," + R + "]), g=sum(i=1,n,v[i]*zk[i]); " +
              "m=abs(nfeltnorm(nf,u*x-g)); if(best==-1||m<best,best=m)));\n";
    script += "print(best);\n";
    script += "quit;\n";

    var sub = spawn(["gp", "-q", "-f"], stdin=pipeStyle.pipe, stdout=pipeStyle.pipe, stderr=pipeStyle.pipe);
    sub.stdin.writeln(script);
    sub.stdin.close();
    var result = "";
    var line: string;
    while sub.stdout.readLine(line) {
      const t = line.strip();
      // gp may emit "cpu time = ..." diagnostics on stdout depending on
      // the local .gprc; the actual printed value is the last real line.
      if t.size > 0 && !t.startsWith("cpu time") then result = t;
    }
    sub.wait();
    return result;
  }
}
