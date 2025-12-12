"""
Stress Test Harness for AlphaOmega Demos

Runs the demo suite 100 times to ensure reliability and catch probabilistic failures.
"""

import sys
import io
import time
import os
from contextlib import redirect_stdout

# Add project root to sys.path for core imports
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
sys.path.insert(0, project_root)

# Import demos (now in same directory)
try:
    import sedenion_sat
    import sae_demo
    import qmesh_demo
    import qmesh_routing_sim
except ImportError as e:
    print(f"Import failed: {e}")
    sys.exit(1)

def run_stress_test(module, name, iterations=100):
    print(f"Stress Testing {name} ({iterations} iterations)...")
    failures = 0
    
    start_total = time.time()
    
    for i in range(iterations):
        # Capture output to avoid spam
        f = io.StringIO()
        try:
            with redirect_stdout(f):
                # Run the demo function
                if name == "Q-Mesh":
                    success = module.run_demo()
                    if not success: raise Exception("Return False")
                elif name == "SAE":
                    module.run_demo() # SAE asserts internally
                elif name == "SAT":
                    module.run_demo() # SAT prints output but doesn't return bool, we check for crash
                elif name == "Routing":
                    module.run_simulation()
                    
        except Exception as e:
            failures += 1
            print(f"  [FAIL] Iteration {i+1}: {e}")
            # print(f.getvalue()) # Print log of failed run
            
        print(f"\r  Progress: {i+1}/{iterations} | Failures: {failures}", end="")
        
    elapsed = time.time() - start_total
    print(f"\n  Completed in {elapsed:.2f}s. Success Rate: {(iterations-failures)/iterations*100:.1f}%")
    
    return failures == 0

if __name__ == "__main__":
    print("=== ALPHAOMEGA STRESS TEST SUITE ===")
    
    all_pass = True
    
    # 1. Q-Mesh (Most flaky due to probabilistic detection)
    if not run_stress_test(qmesh_demo, "Q-Mesh", 100): all_pass = False
    
    # 2. SAE (Deterministic logic but random keys)
    if not run_stress_test(sae_demo, "SAE", 100): all_pass = False
    
    # 3. SAT (Random instance generation)
    if not run_stress_test(sedenion_sat, "SAT", 100): all_pass = False
    
    # 4. Routing (Deterministic graph algo)
    if not run_stress_test(qmesh_routing_sim, "Routing", 20): # Slower
        all_pass = False
        
    if all_pass:
        print("\n[SUCCESS] All demos passed stress testing.")
        sys.exit(0)
    else:
        print("\n[FAIL] Stability issues detected.")
        sys.exit(1)
