// ====================================================================
// UPGRADED STAGE 3: Exact Minimum Extraction Loop
// ====================================================================
proc computeExactEuclideanMinimum(ref kernel: LezowskiKernel, maxDepth: int = 12) {
  var cDom = {1..kernel.degree};
  
  var initialCenter: [cDom] real(64) = 0.5;
  var initialWidths: [cDom] real(64) = 0.5; 
  
  var rootBox = new Box(degree = kernel.degree, centerDom = cDom, center = initialCenter, widths = initialWidths);
  var currentLevelBoxes: list(Box);
  currentLevelBoxes.pushBack(rootBox);
  
  // Track the actual minimum value found across the search space
  var globalMaxLowerBound: atomic real;
  globalMaxLowerBound.write(0.0);
  
  writeln("Launching exact Euclidean minimum extraction pipeline...");
  
  for depth in 1..maxDepth {
    if currentLevelBoxes.size == 0 then break;
    
    writeln("Depth level ", depth, " | Sifting through branches: ", currentLevelBoxes.size);
    
    var nextLevelBoxes: list(Box);
    var nextLevelLock: atomic bool; 
    
    // We use a progressive threshold: any box whose maximum possible norm 
    // is less than the best local minimum found so far can be safely pruned.
    var currentThreshold = globalMaxLowerBound.read();
    // Fallback to a starting threshold if no minimum has been registered yet
    if currentThreshold == 0.0 then currentThreshold = 1.0; 
    
    forall i in 0..currentLevelBoxes.size-1 with (ref nextLevelBoxes, ref nextLevelLock, ref globalMaxLowerBound, const in cDom) {
      var currentBox = currentLevelBoxes[i]; 
      
      // If a box cannot be cleared by the current threshold, it contains a potential minimum
      if !kernel.testBoxCovering(currentBox, currentThreshold) {
        if depth == maxDepth {
          // Point Evaluation Step: Compute the exact norm at the center of the box
          var centers: [1..kernel.num_embeddings] complex(128);
          var radii: [1..kernel.num_embeddings] real(64);
          
          // Create a point-box copy with zero width to isolate the center coordinate
          var pointBox = new Box(degree = kernel.degree, centerDom = cDom, center = currentBox.center, widths = 0.0);
          kernel.computeBoxProjections(pointBox, centers, radii);
          
          var exactCenterNorm = kernel.calculateBoxMaxNorm(centers, radii);
          
          // Thread-safe update of the maximum lower bound discovered
          var existing = globalMaxLowerBound.read();
          while (exactCenterNorm > existing) {
            if globalMaxLowerBound.compareAndSwap(existing, exactCenterNorm) then break;
            existing = globalMaxLowerBound.read();
          }
        } else {
          // Bisect and branch down the tree
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
  
  writeln("\n🏆 Minimum extraction pipeline complete.");
  writeln("Calculated Field Euclidean Minimum Value: ", globalMaxLowerBound.read());
}

proc main() {
  writeln("Loading execution framework...");
  var kernel = loadLezowskiKernel("field_data.txt");
  // Run extraction tracking
  computeExactEuclideanMinimum(kernel, maxDepth = 10);
}
