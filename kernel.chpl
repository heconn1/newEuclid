use IO;
use LinearAlgebra;
//use Complex;

// Core configuration structure to hold our number field invariants
record LezowskiKernel {
  var degree: int;
  var r1: int;
  var r2: int;
  var num_embeddings: int; 
  
  var basisDom: domain(2);
  var basisMatrix: [basisDom] complex(128); // Standardized to 128-bit for high precision
  
  var numUnits: int;
  var unitsDom: domain(2);
  var unitsMatrix: [unitsDom] complex(128);

  // Linear transformation: projects lattice coefficients to mixed embedding space
  proc computeElementEmbeddings(ref coefficients: [] real(64)) {
    var outVector: [1..this.num_embeddings] complex(128);
    
    // Matrix-vector multiplication tailored for complex geometries
    for r in 1..this.num_embeddings {
      var sum: complex(128) = (0.0, 0.0): complex(128);
      for c in 1..this.degree {
        sum += this.basisMatrix[r, c] * coefficients[c];
      }
      outVector[r] = sum;
    }
    return outVector;
  }

  // Calculates Lezowski's generalized absolute field norm for signature (r1, r2)
  proc calculateNorm(ref embVector: [] complex(128)): real(64) {
    var norm: real(64) = 1.0;
    
    // 1. Process Real Embeddings (Product of absolute values)
    for i in 1..this.r1 {
      norm *= abs(embVector[i].re);
    }
    
    // 2. Process Complex Embeddings (Product of absolute squares |z|^2 = x^2 + y^2)
    for j in 1..this.r2 {
      var idx = this.r1 + j;
      var c_val = embVector[idx];
      norm *= (c_val.re * c_val.re) + (c_val.im * c_val.im);
    }
    
    return norm;
  }

  // Core Sieve Routine: Checks if a coordinate is bounded by target_bound
  proc testBoxCovering(ref boxCenter: [] real(64), targetBound: real(64) = 1.0): bool {
    var emb = this.computeElementEmbeddings(boxCenter);
    
    // Test base mapping against the target bound
    if this.calculateNorm(emb) < targetBound then
      return true;
      
    // Apply Fundamental Unit transformation checks
    if this.numUnits > 0 {
      for u in 1..this.numUnits {
        var shiftedEmb: [1..this.num_embeddings] complex(128);
        
        // Element-wise complex multiplication shifts the coordinates
        for r in 1..this.num_embeddings {
          shiftedEmb[r] = emb[r] * this.unitsMatrix[r, u];
        }
        
        if this.calculateNorm(shiftedEmb) < targetBound then
          return true;
      }
    }
    
    return false;
  }
}

// Standalone parser function that builds and returns the configured record cleanly
proc loadLezowskiKernel(filename: string): LezowskiKernel {
  var file = open(filename, ioMode.r);
  var reader = file.reader();
  var line: string;
  
  var local_degree: int;
  var local_r1: int;
  var local_r2: int;
  
  // Phase 1: Parse Scalars safely using explicit 0-based array indexing
  while reader.readLine(line) {
    line = line.strip();
    if line.startsWith("#") || line.size == 0 then continue;
    
    if line.startsWith("degree:") {
      var parts = line.split(":");
      local_degree = parts[1].strip(): int;
    } else if line.startsWith("r1:") {
      var parts = line.split(":");
      local_r1 = parts[1].strip(): int;
    } else if line.startsWith("r2:") {
      var parts = line.split(":");
      local_r2 = parts[1].strip(): int;
    } else if line.startsWith("INTEGRAL_BASIS_EMBEDDINGS:") {
      break; 
    }
  }
  
  var local_num_embeddings = local_r1 + local_r2;
  var local_basisDom = {1..local_num_embeddings, 1..local_degree};
  var local_basisMatrix: [local_basisDom] complex(128);
  
  // Phase 2: Read Integral Basis Embeddings Matrix
  for col in 1..local_degree {
    while reader.readLine(line) {
      line = line.strip();
      if !line.startsWith("#") && line.size > 0 then break;
    }
    
    for row in 1..local_num_embeddings {
      if row > 1 then reader.readLine(line);
      var parts = line.strip().split();
      var re = parts[0]: real(64);
      var im = parts[1]: real(64);
      // FIXED: Removed 'new' keyword and cast to matching bit-width complex primitive
      local_basisMatrix[row, col] = (re, im): complex(128);
    }
  }
  
  // Phase 3: Parse Units Invariants
  var local_numUnits: int = 0;
  while reader.readLine(line) {
    line = line.strip();
    if line.startsWith("FUNDAMENTAL_UNITS_COUNT:") {
      var parts = line.split(":");
      local_numUnits = parts[1].strip(): int;
      break;
    }
  }
  
  var local_unitsDom = {1..local_num_embeddings, 1..max(1, local_numUnits)};
  var local_unitsMatrix: [local_unitsDom] complex(128);
  
  if local_numUnits > 0 {
    for u in 1..local_numUnits {
      while reader.readLine(line) {
        if line.strip().startsWith("UNIT_") then break;
      }
      for row in 1..local_num_embeddings {
        reader.readLine(line);
        var parts = line.strip().split();
        var re = parts[0]: real(64);
        var im = parts[1]: real(64);
        // FIXED: Removed 'new' keyword and cast to matching bit-width complex primitive
        local_unitsMatrix[row, u] = (re, im): complex(128);
      }
    }
  }
  
  reader.close();
  file.close();

  // Construct and return the fully ready structural record instance directly
  return new LezowskiKernel(
    degree = local_degree,
    r1 = local_r1,
    r2 = local_r2,
    num_embeddings = local_num_embeddings,
    basisDom = local_basisDom,
    basisMatrix = local_basisMatrix,
    numUnits = local_numUnits,
    unitsDom = local_unitsDom,
    unitsMatrix = local_unitsMatrix
  );
}

// Entry point for local validation
proc main() {
  writeln("Initializing Chapel Lezowski computational kernel via standalone parser...");
  
  var kernel = loadLezowskiKernel("field_data.txt");
  writeln("Successfully loaded field data!");
  writeln("Degree: ", kernel.degree, " | Signature: (", kernel.r1, ", ", kernel.r2, ")");
  
  // Create a sample fractional test midpoint vector [0.5, ... 0.5]
  var testCoord: [1..kernel.degree] real(64);
  testCoord = 0.5;
  
  var covered = kernel.testBoxCovering(testCoord, 1.0);
  writeln("Is coordinate covered below norm bound 1.0? -> ", covered);
}
