// ====================================================================
// STAGE 3: Parallel Adaptive Mesh Refinement (AMR) Loop (FIXED)
// ====================================================================
proc runAMRTreeSieve(ref kernel: LezowskiKernel, targetBound: real(64) = 1.0, maxDepth: int = 12) {
  var cDom = {1..kernel.degree};
  
  // Create our initial root box representing the fundamental domain [0, 1)^degree
  var initialCenter: [cDom] real(64);
  var initialWidths: [cDom] real(64);
  initialCenter = 0.5;
  initialWidths = 0.5; // Bounds out domain cleanly
  
  var rootBox = new Box(degree = kernel.degree, centerDom = cDom, center = initialCenter, widths = initialWidths);
  
  // Maintain a worklist of boxes that need processing at the current tree depth
  var currentLevelBoxes: list(Box);
  currentLevelBoxes.push_back(rootBox); // FIXED: changed from append to push_back
  
  writeln("Starting parallel AMR sieve strategy...");
  
  for depth in 1..maxDepth {
    if currentLevelBoxes.size == 0 {
      writeln("🎉 Success! All branches completely covered. The Euclidean minimum is STRICTLY LESS than ", targetBound);
      return;
    }
    
    writeln("Depth level ", depth, " | Active boxes to process: ", currentLevelBoxes.size);
    
    // Setup a next-level thread-safe collector list for the failing sub-branches
    var nextLevelBoxes: list(Box);
    var nextLevelLock: atomic int; // Atomic primitive to safely manage parallel list appends
    
    // DATA PARALLEL LOOP: Chapel scales this natively across all CPU hardware threads
    forall i in 0..currentLevelBoxes.size-1 {
      var currentBox = currentLevelBoxes.getValue(i);
      
      // Execute our Stage 2 math sieve
      if !kernel.testBoxCovering(currentBox, targetBound) {
        // If it cannot be cleared and we reached max depth, these are likely critical points!
        if depth == maxDepth {
          // Atomically output findings without printing collisions
          writeln("⚠️ Target un-cleared at maximum depth resolution near center: ", currentBox.center);
        } else {
          // Perform 2^n coordinate bisection branching
          var numSubBoxes = 1 << kernel.degree;
          
          for bIdx in 0..numSubBoxes-1 {
            var subCenter: [cDom] real(64);
            var subWidths: [cDom] real(64);
            
            for d in 1..kernel.degree {
              subWidths[d] = currentBox.widths[d] * 0.5;
              // Binary masking shifts the center cleanly to the left or right sub-quadrant
              var direction = if (bIdx & (1 << (d - 1))) != 0 then 1.0 else -1.0;
              subCenter[d] = currentBox.center[d] + (direction * subWidths[d]);
            }
            
            var subBox = new Box(degree = kernel.degree, centerDom = cDom, center = subCenter, widths = subWidths);
            
            // Critical section protected by dynamic atomic locking spin loops
            while nextLevelLock.testAndSet() do chpl_task_yield();
            nextLevelBoxes.push_back(subBox); // FIXED: changed from append to push_back
            nextLevelLock.clear();
          }
        }
      }
    }
    
    // Advance generation parameters
    currentLevelBoxes = nextLevelBoxes;
  }
  
  if currentLevelBoxes.size > 0 {
    writeln("🛑 Algorithm reached Max Depth limits with ", currentLevelBoxes.size, " remaining suspect coordinates.");
    writeln("These points define the upper structural ceiling or critical points of your Euclidean minimum value.");
  }
}

