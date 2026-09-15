// NumberField.chpl
//
// Loads the (one-time, Pari/gp-computed) field_data.txt file produced by
// generate_field.gp into a plain data record. This is the only place that
// needs to understand the on-disk file format; everything downstream works
// purely in floating point.
//
// field_data.txt layout (see generate_field.gp):
//   degree: <n>
//   r1: <r1>
//   r2: <r2>
//   INTEGRAL_BASIS_EMBEDDINGS:
//     (degree blocks of (r1+r2) "re<TAB>im" rows, one block per integral
//      basis generator, in order)
//   FUNDAMENTAL_UNITS_COUNT: <numUnits>
//     (numUnits blocks of (r1+r2) "re<TAB>im" rows, one block per
//      fundamental unit, in order)
//
// The embedding convention (matching generate_field.gp / nfeltembed) is:
//   rows 1..r1       -> real embeddings (im part always 0)
//   rows r1+1..r1+r2 -> one representative per conjugate pair of complex
//                       embeddings (the conjugate is implied, never stored)
module NumberField {
  use IO;
  use List;

  record NumberFieldData {
    var polynomial: string;
    var degree: int;
    var r1: int;
    var r2: int;
    var numEmbeddings: int; // = r1 + r2

    // basis[row, c] = embedding #row of the c-th integral basis generator
    var basisDom: domain(2) = {1..0, 1..0};
    var basis: [basisDom] complex(128);

    var numUnits: int;
    // unitEmbeddings[row, u] = embedding #row of the u-th fundamental unit
    var unitsDom: domain(2) = {1..0, 1..0};
    var unitEmbeddings: [unitsDom] complex(128);

    // Inverse of the real (r1+2r2)x(r1+2r2) embedding matrix, i.e. maps an
    // *unpacked* real embedding vector (see unpackEmbedding) back to
    // integral-basis coefficients. Used only by the unit-action
    // acceleration (Sieve.chpl); Phase 1 absorption never needs it.
    var sigmaInvDom: domain(2) = {1..0, 1..0};
    var sigmaInv: [sigmaInvDom] real(64);
  }

  // Gauss-Jordan matrix inversion with partial pivoting -- kept dependency
  // free (Chapel's LinearAlgebra.inv() requires a system LAPACK we cannot
  // assume is installed).
  proc invertMatrix(const ref M: [] real(64), n: int): [1..n, 1..n] real(64) {
    var A: [1..n, 1..(2*n)] real(64);
    for i in 1..n {
      for j in 1..n do A[i, j] = M[i, j];
      A[i, n+i] = 1.0;
    }
    for col in 1..n {
      var piv = col;
      var best = abs(A[col, col]);
      for r in col+1..n {
        if abs(A[r, col]) > best { best = abs(A[r, col]); piv = r; }
      }
      if best < 1.0e-13 then halt("NumberField.invertMatrix: embedding matrix is singular");
      if piv != col then for j in 1..2*n {
        const t = A[col, j]; A[col, j] = A[piv, j]; A[piv, j] = t;
      }
      const d = A[col, col];
      for j in 1..2*n do A[col, j] /= d;
      for r in 1..n {
        if r != col {
          const f = A[r, col];
          if f != 0.0 then for j in 1..2*n do A[r, j] -= f * A[col, j];
        }
      }
    }
    var inv: [1..n, 1..n] real(64);
    for i in 1..n do for j in 1..n do inv[i, j] = A[i, n+j];
    return inv;
  }

  private proc isDataLine(t: string): bool {
    if t.size == 0 then return false;
    if t.startsWith("#") then return false;
    if t.startsWith("polynomial:") then return false;
    if t.startsWith("degree:") then return false;
    if t.startsWith("r1:") then return false;
    if t.startsWith("r2:") then return false;
    if t.startsWith("discriminant:") then return false;
    if t.startsWith("INTEGRAL_BASIS_EMBEDDINGS:") then return false;
    if t.startsWith("FUNDAMENTAL_UNITS_COUNT:") then return false;
    if t.startsWith("UNIT_") then return false;
    return true;
  }

  proc loadNumberField(filename: string): NumberFieldData throws {
    var file = open(filename, ioMode.r);
    var reader = file.reader(locking=false);
    var line: string;

    var degree = 0, r1 = 0, r2 = 0, numUnits = 0;
    var polynomial = "";
    var dataRe: list(real(64));
    var dataIm: list(real(64));

    while reader.readLine(line) {
      const t = line.strip();
      if t.startsWith("polynomial:") {
        const parts = t.split(":");
        polynomial = parts[1].strip();
      } else if t.startsWith("degree:") {
        const parts = t.split(":");
        const v: string = parts[1].strip();
        degree = v: int;
      } else if t.startsWith("r1:") {
        const parts = t.split(":");
        const v: string = parts[1].strip();
        r1 = v: int;
      } else if t.startsWith("r2:") {
        const parts = t.split(":");
        const v: string = parts[1].strip();
        r2 = v: int;
      } else if t.startsWith("FUNDAMENTAL_UNITS_COUNT:") {
        const parts = t.split(":");
        const v: string = parts[1].strip();
        numUnits = v: int;
      } else if isDataLine(t) {
        const parts = t.split();
        if parts.size == 2 {
          dataRe.pushBack(parts[0]: real(64));
          dataIm.pushBack(parts[1]: real(64));
        }
      }
    }
    reader.close();
    file.close();

    var nf: NumberFieldData;
    nf.polynomial = polynomial;
    nf.degree = degree;
    nf.r1 = r1;
    nf.r2 = r2;
    const numEmbeddings = r1 + r2;
    nf.numEmbeddings = numEmbeddings;
    nf.basisDom = {1..numEmbeddings, 1..degree};
    nf.numUnits = numUnits;
    nf.unitsDom = {1..numEmbeddings, 1..max(1, numUnits)};

    const expected = numEmbeddings * degree + numUnits * numEmbeddings;
    if dataRe.size != expected then
      halt("field_data.txt: expected ", expected, " embedding rows, found ", dataRe.size);

    var idx = 0;
    for c in 1..degree {
      for row in 1..numEmbeddings {
        nf.basis[row, c] = (dataRe[idx], dataIm[idx]): complex(128);
        idx += 1;
      }
    }
    for u in 1..numUnits {
      for row in 1..numEmbeddings {
        nf.unitEmbeddings[row, u] = (dataRe[idx], dataIm[idx]): complex(128);
        idx += 1;
      }
    }

    // Build the real, unpacked (degree x degree) embedding matrix and
    // invert it once. Row order matches unpackEmbedding: r1 real rows,
    // then r2 "real part of complex embedding" rows, then r2 "imaginary
    // part" rows.
    var M: [1..degree, 1..degree] real(64);
    for c in 1..degree {
      for i in 1..r1 do M[i, c] = nf.basis[i, c].re;
      for j in 1..r2 {
        M[r1+j, c] = nf.basis[r1+j, c].re;
        M[r1+r2+j, c] = nf.basis[r1+j, c].im;
      }
    }
    nf.sigmaInvDom = {1..degree, 1..degree};
    nf.sigmaInv = invertMatrix(M, degree);
    return nf;
  }

  // Unpacks a complex-compact embedding vector into a genuine real
  // (degree)-vector: r1 real embeddings, then Re() of the r2 complex
  // embeddings, then Im() of the r2 complex embeddings.
  proc unpackEmbedding(const ref nf: NumberFieldData, const ref e: [] complex(128)): [1..nf.degree] real(64) {
    var v: [1..nf.degree] real(64);
    for i in 1..nf.r1 do v[i] = e[i].re;
    for j in 1..nf.r2 {
      v[nf.r1 + j] = e[nf.r1 + j].re;
      v[nf.r1 + nf.r2 + j] = e[nf.r1 + j].im;
    }
    return v;
  }

  // Inverse of `embed` composed with unpacking: given a target embedding
  // (complex-compact), returns the (generally non-integral) coefficient
  // vector that would produce it.
  proc coefficientsFromEmbedding(const ref nf: NumberFieldData, const ref e: [] complex(128)): [1..nf.degree] real(64) {
    const v = unpackEmbedding(nf, e);
    var w: [1..nf.degree] real(64);
    forall i in 1..nf.degree {
      var s = 0.0;
      for k in 1..nf.degree do s += nf.sigmaInv[i, k] * v[k];
      w[i] = s;
    }
    return w;
  }

  // Embedding of an integer coefficient vector (length = degree) in the
  // mixed r1-real / r2-complex embedding space.
  proc embed(const ref nf: NumberFieldData, const ref coeffs: [] real(64)): [1..nf.numEmbeddings] complex(128) {
    var outVec: [1..nf.numEmbeddings] complex(128);
    forall row in 1..nf.numEmbeddings {
      var s: complex(128) = 0.0;
      for c in 1..nf.degree do s += nf.basis[row, c] * coeffs[c];
      outVec[row] = s;
    }
    return outVec;
  }

  // |N(x)| given its embedding vector.
  proc absNorm(const ref nf: NumberFieldData, const ref emb: [] complex(128)): real(64) {
    var n: real(64) = 1.0;
    for i in 1..nf.r1 do n *= abs(emb[i].re);
    for j in 1..nf.r2 {
      const z = emb[nf.r1 + j];
      n *= z.re*z.re + z.im*z.im;
    }
    return n;
  }
}
