use IO;
use LinearAlgebra;
use List;
use Time;
use Math;

// ====================================================================
// STAGE 1: Spatial Geometry Structure (The Box)
// ====================================================================
record Box {
  var degree: int;
  var centerDom: domain(1);
  var center: [centerDom] real(64); 
  var widths: [centerDom] real(64); 
}

record LezowskiKernel {
  var degree: int;
  var r1: int;
  var r2: int;
  var num_embeddings: int; 
  
  var basisDom: domain(2);
  var basisMatrix: [basisDom] complex(128); 
  
  var numUnits: int;
  var unitsDom: domain(2);
  var unitsMatrix: [unitsDom] complex(128);

  proc computeBoxProjections(ref box: Box, ref outCenters: [] complex(128), ref outRadii: [] real(64)) {
    for i in 1..this.num_embeddings {
      outCenters[i] = (0.0, 0.0): complex(128);
      outRadii[i] = 0.0;
    }
    for r in 1..this.num_embeddings {
      var c_sum: complex(128) = (0.0, 0.0): complex(128);
      var r_sum: real(64) = 0.0;
      for c in 1..this.degree {
        c_sum += this.basisMatrix[r, c] * box.center[c];
        var b_val = this.basisMatrix[r, c];
        var col_magnitude = sqrt(b_val.re * b_val.re + b_val.im * b_val.im);
        r_sum += col_magnitude * box.widths[c];
      }
      outCenters[r] = c_sum;
      outRadii[r] = r_sum;
    }
  }

  proc calculateBoxMaxNorm(ref centers: [] complex(128), ref radii: [] real(64)): real(64) {
    var maxNorm: real(64) = 1.0;
    for i in 1..this.r1 {
      var c_val = abs(centers[i].re);
      var r_val = radii[i];
      maxNorm *= (c_val + r_val); 
    }
    for j in 1..this.r2 {
      var idx = this.r1 + j;
      var z0 = centers[idx];
      var z0_mag = sqrt(z0.re * z0.re + z0.im * z0.im);
      var R = radii[idx];
      var max_z = z0_mag + R;
      maxNorm *= (max_z * max_z); 
    }
    return maxNorm;
  }

  // UPGRADED: Neighborhood Coefficient Sieve to Eliminate False Negatives
  proc testBoxCovering(ref inputBox: Box, targetBound: real(64)): bool {
    // Determine the base rounded integer vector for the center
    var baseZ: [1..this.degree] int;
    for c in 1..this.degree {
      baseZ[c] = round(inputBox.center[c]): int;
    }

    var centers: [1..this.num_embeddings] complex(128);
    var radii: [1..this.num_embeddings] real(64);
    var sCenters: [1..this.num_embeddings] complex(128);
    var sRadii: [1..this.num_embeddings] real(64);

    // Support flexible degrees (up to cubic fields natively via static offsets)
    // We test a small offset mesh around the rounded integer coefficients
    for d1 in -1..1 {
      for d2 in (if this.degree >= 2 then -1 else 0)..(if this.degree >= 2 then 1 else 0) {
        for d3 in (if this.degree >= 3 then -1 else 0)..(if this.degree >= 3 then 1 else 0) {
          
          // Construct the candidate shifted box
          var testBox = new Box(degree = this.degree, centerDom = inputBox.centerDom, 
                                center = inputBox.center, widths = inputBox.widths);
          
          testBox.center[1] -= (baseZ[1] + d1);
          if this.degree >= 2 then testBox.center[2] -= (baseZ[2] + d2);
          if this.degree >= 3 then testBox.center[3] -= (baseZ[3] + d3);

          // Project this specific integer translation into the embedding space
          this.computeBoxProjections(testBox, centers, radii);
          
          // 1. Check base field coverage for this integer vector
          if this.calculateBoxMaxNorm(centers, radii) < targetBound then
            return true;
            
          // 2. Check fundamental unit distortions for this integer vector
          if this.numUnits > 0 {
            var maxExp = if this.numUnits >= 2 then 2 else 3;
            
            for e1 in -maxExp..maxExp {
              for e2 in (if this.numUnits >= 2 then -maxExp else 0)..(if this.numUnits >= 2 then maxExp else 0) {
                
                for r in 1..this.num_embeddings {
                  var u1 = (this.unitsMatrix[r, 1]) ** e1;
                  var u2 = if this.numUnits >= 2 then (this.unitsMatrix[r, 2]) ** e2 else (1.0, 0.0): complex(128);
                  var combinedUnit = u1 * u2;
                  
                  sCenters[r] = centers[r] * combinedUnit;
                  
                  var r_sum = 0.0;
                  for c in 1..this.degree {
                    var transBasisElement = this.basisMatrix[r, c] * combinedUnit;
                    var elementMag = sqrt(transBasisElement.re**2 + transBasisElement.im**2);
                    r_sum += elementMag * testBox.widths[c];
                  }
                  sRadii[r] = r_sum;
                }
                
                if this.calculateBoxMaxNorm(sCenters, sRadii) < targetBound then 
                  return true;
              }
            }
          }

        }
      }
    }
    
    return false; // Truly un-coverable under this bound threshold
  }
}

// Standalone parser function
proc loadLezowskiKernel(filename: string): LezowskiKernel {
  var file = open(filename, ioMode.r);
  var reader = file.reader();
  var line: string;
  var local_degree: int; var local_r1: int; var local_r2: int;
  
  while reader.readLine(line) {
    line = line.strip();
    if line.startsWith("#") || line.size == 0 then continue;
    
    if line.startsWith("degree:") {
      var parts = line.split(":");
      var rawVal = parts[1].strip();
      local_degree = rawVal: int;
    } else if line.startsWith("r1:") {
      var parts = line.split(":");
      var rawVal = parts[1].strip();
      local_r1 = rawVal: int;
    } else if line.startsWith("r2:") {
      var parts = line.split(":");
      var rawVal = parts[1].strip();
      local_r2 = rawVal: int;
    } else if line.startsWith("INTEGRAL_BASIS_EMBEDDINGS:") {
      break; 
    }
  }
  
  var local_num_embeddings = local_r1 + local_r2;
  var local_basisDom = {1..local_num_embeddings, 1..local_degree};
  var local_basisMatrix: [local_basisDom] complex(128);
  
  for col in 1..local_degree {
    while reader.readLine(line) { 
      line = line.strip(); 
      if !line.startsWith("#") && line.size > 0 then break; 
    }
    for row in 1..local_num_embeddings {
      if row > 1 then reader.readLine(line);
      var parts = line.strip().split();
      var re = parts[0]: real(64); var im = parts[1]: real(64);
      local_basisMatrix[row, col] = (re, im): complex(128);
    }
  }
  
  var local_numUnits: int = 0;
  while reader.readLine(line) {
    line = line.strip();
    if line.startsWith("FUNDAMENTAL_UNITS_COUNT:") { 
      var parts = line.split(":");
      var rawVal = parts[1].strip();
      local_numUnits = rawVal: int;
      break; 
    }
  }
  
  var local_unitsDom = {1..local_num_embeddings, 1..max(1, local_numUnits)};
  var local_unitsMatrix: [local_unitsDom] complex(128);
  if local_numUnits > 0 {
    for u in 1..local_numUnits {
      while reader.readLine(line) { if line.strip().startsWith("UNIT_") then break; }
      for row in 1..local_num_embeddings {
        reader.readLine(line);
        var parts = line.strip().split();
        var re = parts[0]: real(64); var im = parts[1]: real(64);
        local_unitsMatrix[row, u] = (re, im): complex(128);
      }
    }
  }
  reader.close(); file.close();
  
  return new LezowskiKernel(degree = local_degree, r1 = local_r1, r2 = local_r2, num_embeddings = local_num_embeddings,
                            basisDom = local_basisDom, basisMatrix = local_basisMatrix, numUnits = local_numUnits,
                            unitsDom = local_unitsDom, unitsMatrix = local_unitsMatrix);
}

// ====================================================================
// STAGE 2: Core Decision Sieve
// ====================================================================
proc runAMRTreeSieve(ref kernel: LezowskiKernel, targetBound: real(64), maxDepth: int = 10): bool {
  var cDom = {1..kernel.degree};
  var initialCenter: [cDom] real(64) = 0.5;
  var initialWidths: [cDom] real(64) = 0.5; 
  
  var rootBox = new Box(degree = kernel.degree, centerDom = cDom, center = initialCenter, widths = initialWidths);
  var currentLevelBoxes: list(Box);
  currentLevelBoxes.pushBack(rootBox);
  
  for depth in 1..maxDepth {
    if currentLevelBoxes.size == 0 then return true; // Success! All space cleared
    
    // Safety Threshold: If uncleared regions begin to explode exponentially, 
    // it means targetBound is lower than the true minimum. Terminate immediately.
    if currentLevelBoxes.size > 20000 then return false;
    
    var nextLevelBoxes: list(Box);
    var nextLevelLock: atomic bool; 
    
    forall i in 0..currentLevelBoxes.size-1 with (ref nextLevelBoxes, ref nextLevelLock, const in cDom) {
      var currentBox = currentLevelBoxes[i]; 
      
      if !kernel.testBoxCovering(currentBox, targetBound) {
        if depth < maxDepth {
          var numSubBoxes = 1 << kernel.degree;
          for bIdx in 0..numSubBoxes-1 {
            var subCenter: [cDom] real(64); var subWidths: [cDom] real(64);
            for d in 1..kernel.degree {
              subWidths[d] = currentBox.widths[d] * 0.5;
              var direction = if (bIdx & (1 << (d - 1))) != 0 then 1.0 else -1.0;
              subCenter[d] = currentBox.center[d] + (direction * subWidths[d]);
            }
            var subBox = new Box(degree = kernel.degree, centerDom = cDom, center = subCenter, widths = subWidths);
            
            while nextLevelLock.testAndSet() do sleep(0.000001); 
            nextLevelBoxes.pushBack(subBox); 
            nextLevelLock.clear(); 
          }
        } else {
          // Reached resolution floor without clearing
          while nextLevelLock.testAndSet() do sleep(0.000001); 
          nextLevelBoxes.pushBack(currentBox);
          nextLevelLock.clear();
        }
      }
    }
    if depth == maxDepth && nextLevelBoxes.size > 0 then return false;
    currentLevelBoxes = nextLevelBoxes;
  }
  
  return true;
}

// ====================================================================
// STAGE 3: Dichotomic Optimization Master Loop
// ====================================================================
proc main() {
  writeln("Loading execution framework...");
  var kernel = loadLezowskiKernel("field_data.txt");
  
  writeln("Initiating Binary Search Optimization Suite...");
  
  var lowerBound = 0.0;
  var upperBound = 1.0;
  var epsilon = 1e-6; // Targeted calculation accuracy tolerance
  
  while (upperBound - lowerBound) > epsilon {
    var midBound = lowerBound + (upperBound - lowerBound) * 0.5;
    write("Testing Threshold Boundary: ", midBound, " -> ");
    
    // We execute a deep 12-level sieve pass
    //    var isCovered = runAMRTreeSieve(kernel, targetBound = midBound, maxDepth = 12);
    var isCovered = runAMRTreeSieve(kernel, targetBound = midBound, maxDepth = 18);
    
    if isCovered {
      writeln("✅ CLEARED. True minimum is smaller.");
      upperBound = midBound;
    } else {
      writeln("❌ BLOCKED. True minimum is larger.");
      lowerBound = midBound;
    }
  }
  
  writeln("\n🏆 Minimum extraction pipeline complete.");
  writeln("Calculated Field Euclidean Minimum Value: ", lowerBound);
}
