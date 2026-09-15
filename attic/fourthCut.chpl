use IO;
use LinearAlgebra;
use List;
use Time;

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

  // UPGRADED: Expanded Exponent Sieve (-2..2 Mesh) to capture rank-2 boundary regions
  proc testBoxCovering(ref inputBox: Box, targetBound): bool {
    var box = new Box(degree = this.degree, centerDom = inputBox.centerDom, 
                      center = inputBox.center, widths = inputBox.widths);
    
    for c in 1..this.degree {
      var absoluteClosestInteger = round(box.center[c]);
      box.center[c] -= absoluteClosestInteger; 
    }

    var centers: [1..this.num_embeddings] complex(128);
    var radii: [1..this.num_embeddings] real(64);
    this.computeBoxProjections(box, centers, radii);
    
    if this.calculateBoxMaxNorm(centers, radii) < targetBound then
      return true;
      
    if this.numUnits == 1 {
      for exp1 in -2..2 {
        if exp1 == 0 then continue;
        var sCenters: [1..this.num_embeddings] complex(128);
        var sRadii: [1..this.num_embeddings] real(64);
        for r in 1..this.num_embeddings {
          var u_val = (this.unitsMatrix[r, 1]) ** exp1; // Handles fractional exponents via power
          sCenters[r] = centers[r] * u_val;
          var u_mag = sqrt(u_val.re * u_val.re + u_val.im * u_val.im);
          sRadii[r] = radii[r] * u_mag;
        }
        if this.calculateBoxMaxNorm(sCenters, sRadii) < targetBound then return true;
      }
    } else if this.numUnits >= 2 {
      // Explores 5^2 = 25 transformation combinations to clear the boundaries
      for exp1 in -2..2 {
        for exp2 in -2..2 {
          if exp1 == 0 && exp2 == 0 then continue;
          var sCenters: [1..this.num_embeddings] complex(128);
          var sRadii: [1..this.num_embeddings] real(64);
          for r in 1..this.num_embeddings {
            var u1 = (this.unitsMatrix[r, 1]) ** exp1;
            var u2 = (this.unitsMatrix[r, 2]) ** exp2;
            var combinedUnit = u1 * u2;
            
            sCenters[r] = centers[r] * combinedUnit;
            var u_mag = sqrt(combinedUnit.re * combinedUnit.re + combinedUnit.im * combinedUnit.im);
            sRadii[r] = radii[r] * u_mag;
          }
          if this.calculateBoxMaxNorm(sCenters, sRadii) < targetBound then return true;
        }
      }
    }
    
    return false; 
  }
}

// Standalone parser function (Synchronized and Fixed)
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
      local_degree = parts[1].strip(): int; // FIXED: explicitly indexing element 1
    } else if line.startsWith("r1:") {
      var parts = line.split(":");
      local_r1 = parts[1].strip(): int;     // FIXED: explicitly indexing element 1
    } else if line.startsWith("r2:") {
      var parts = line.split(":");
      local_r2 = parts[1].strip(): int;     // FIXED: explicitly indexing element 1
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
      local_numUnits = parts[1].strip(): int; // FIXED: explicitly indexing element 1
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
// STAGE 3: Parallel Adaptive Mesh Refinement (AMR) Loop
// ====================================================================
proc runAMRTreeSieve(ref kernel: LezowskiKernel, targetBound, maxDepth: int = 12) {
  var cDom = {1..kernel.degree};
  var initialCenter: [cDom] real(64); var initialWidths: [cDom] real(64);
  initialCenter = 0.5; initialWidths = 0.5; 
  
  var rootBox = new Box(degree = kernel.degree, centerDom = cDom, center = initialCenter, widths = initialWidths);
  var currentLevelBoxes: list(Box);
  currentLevelBoxes.pushBack(rootBox); 
  
  writeln("Starting parallel AMR sieve strategy...");
  var totalSuspects: atomic int; // Clean thread-safe counter for max depth tracking
  
  for depth in 1..maxDepth {
    if currentLevelBoxes.size == 0 {
      writeln("🎉 Success! All branches completely covered. The Euclidean minimum is STRICTLY LESS than ", targetBound);
      return;
    }
    
    writeln("Depth level ", depth, " | Active boxes to process: ", currentLevelBoxes.size);
    var nextLevelBoxes: list(Box);
    var nextLevelLock: atomic bool; 
    
    forall i in 0..currentLevelBoxes.size-1 with (ref nextLevelBoxes, ref nextLevelLock, ref totalSuspects, const in cDom) {
      var currentBox = currentLevelBoxes[i]; 
      
      if !kernel.testBoxCovering(currentBox, targetBound) {
        if depth == maxDepth {
          totalSuspects.add(1); // Quietly increment without printing collisions
        } else {
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
        }
      }
    }
    currentLevelBoxes = nextLevelBoxes;
  }
  
  var suspectCount = totalSuspects.read();
  if suspectCount > 0 {
    writeln("\n🛑 Sieve execution completed with limits.");
    writeln("Total remaining un-cleared suspect coordinates at resolution ceiling: ", suspectCount);
    writeln("This implies the field's actual Euclidean minimum or critical point boundary is near this threshold region.");
  }
}

config const targetBound = 1.0;

proc main() {
  writeln("Loading execution framework...");
  var kernel = loadLezowskiKernel("field_data.txt");
  runAMRTreeSieve(kernel, targetBound, maxDepth = 10);
}


