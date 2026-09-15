use IO;
use LinearAlgebra;
//use Complex;

// ====================================================================
// STAGE 1: Spatial Geometry Structure (The Box)
// ====================================================================
record Box {
  var degree: int;
  
  // Coordinates are defined relative to the integral basis of O_K
  var centerDom: domain(1);
  var center: [centerDom] real(64); // Midpoint of the box along each generator axis
  var widths: [centerDom] real(64); // Half-widths (radii) along each generator axis
}

// Core configuration structure to hold our number field invariants
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

  // Helper: Computes the precise upper bound of a matrix projection for interval math
  proc computeBoxProjections(ref box: Box, ref outCenters: [] complex(128), ref outRadii: [] real(64)) {
    // Zero-out arrays first
    for i in 1..this.num_embeddings {
      outCenters[i] = (0.0, 0.0): complex(128);
      outRadii[i] = 0.0;
    }
    
    // Linearly map the center and accumulate the maximum possible geometric distortion (radius)
    for r in 1..this.num_embeddings {
      var c_sum: complex(128) = (0.0, 0.0): complex(128);
      var r_sum: real(64) = 0.0;
      
      for c in 1..this.degree {
        c_sum += this.basisMatrix[r, c] * box.center[c];
        
        // Minkowski metric bounding: calculate maximum stretch along this basis embedding column
        var b_val = this.basisMatrix[r, c];
        var col_magnitude = sqrt(b_val.re * b_val.re + b_val.im * b_val.im);
        r_sum += col_magnitude * box.widths[c];
      }
      
      outCenters[r] = c_sum;
      outRadii[r] = r_sum;
    }
  }

  // Computes the strict upper bound of the norm function over the projected mixed regions
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
  
  // ====================================================================
  // REVISED STAGE 2: Sieve Routine with Integer Translation Step
  // ====================================================================
  proc testBoxCovering(ref inputBox: Box, targetBound: real(64) = 1.0): bool {
    // 1. Create a translated copy of the box to preserve the original structure
    var box = new Box(degree = this.degree, centerDom = inputBox.centerDom, 
                      center = inputBox.center, widths = inputBox.widths);
    
    // 2. Perform O_K Integer Translation: Round the center coordinates to the nearest integers
    for c in 1..this.degree {
      var absoluteClosestInteger = round(box.center[c]);
      box.center[c] -= absoluteClosestInteger; // Translate center to the fundamental domain
    }

    var centers: [1..this.num_embeddings] complex(128);
    var radii: [1..this.num_embeddings] real(64);
    
    // 3. Project the translated algebraic box into embedding spaces
    this.computeBoxProjections(box, centers, radii);
    
    // 4. Check if the translated region is completely covered by 0 in the new domain
    if this.calculateBoxMaxNorm(centers, radii) < targetBound then
      return true;
      
    // 5. Apply Fundamental Unit transformation checks
    if this.numUnits > 0 {
      for u in 1..this.numUnits {
        var shiftedCenters: [1..this.num_embeddings] complex(128);
        var shiftedRadii: [1..this.num_embeddings] real(64);
        
        for r in 1..this.num_embeddings {
          shiftedCenters[r] = centers[r] * this.unitsMatrix[r, u];
          
          var u_val = this.unitsMatrix[r, u];
          var u_mag = sqrt(u_val.re * u_val.re + u_val.im * u_val.im);
          shiftedRadii[r] = radii[r] * u_mag;
        }
        
        if this.calculateBoxMaxNorm(shiftedCenters, shiftedRadii) < targetBound then
          return true; 
      }
    }
    
    return false; // Box must be broken down further
  }
}


// Standalone parser function (kept intact from prior fix)
proc loadLezowskiKernel(filename: string): LezowskiKernel {
  var file = open(filename, ioMode.r);
  var reader = file.reader();
  var line: string;
  
  var local_degree: int;
  var local_r1: int;
  var local_r2: int;
  
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
    } else if line.startsWith("INTEGRAL_BASIS_EMBEDDINGS:") then
      break; 
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
      var re = parts[0]: real(64);
      var im = parts[1]: real(64);
      local_basisMatrix[row, col] = (re, im): complex(128);
    }
  }
  
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
        local_unitsMatrix[row, u] = (re, im): complex(128);
      }
    }
  }
  
  reader.close();
  file.close();

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

// Validation entry point
proc main() {
  writeln("Initializing Stages 1 & 2 Sandbox...");
  var kernel = loadLezowskiKernel("field_data.txt");
  
  // Set up a structural test box: Center point at [0.5, 0.5, 0.5] with a half-width size of 0.1
  var cDom = {1..kernel.degree};
  var testCenter: [cDom] real(64);
  var testWidths: [cDom] real(64);
  
  testCenter = 0.5;
  testWidths = 0.1; // Large exploratory region footprint
  
  var testBox = new Box(degree = kernel.degree, centerDom = cDom, center = testCenter, widths = testWidths);
  
  var isCleared = kernel.testBoxCovering(testBox, targetBound = 1.0);
  writeln("Sieve check completed.");
  writeln("Can this continuous geometric box be entirely eliminated? -> ", isCleared);
}
