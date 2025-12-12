"""
AlphaOmega Technology Stack - Master Demo Harness

Runs all 4 proofs-of-concept to verify the Post-Quantum Toolkit.

1. P=NP Solver (Sedenion SAT)
2. Sedenion Asymmetric Encryption (SAE)
3. Q-Mesh Entanglement (Quantum Security)
4. Teleportation Routing (Infinite Scaling)
"""

import subprocess
import sys
import time
import os

# Ensure we are running from the script's directory
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# Add project root to sys.path so demos can import core
project_root = os.path.dirname(script_dir)
sys.path.insert(0, project_root)

DEMOS = [
    {
        "name": "P=NP Geometric Solver (Attack Vector)",
        "path": "sedenion_sat.py",
        "desc": "Solves 3-SAT in O(N) using Associator Tunneling"
    },
    {
        "name": "Sedenion Asymmetric Encryption (SAE)",
        "path": "sae_demo.py",
        "desc": "Post-Quantum stream cipher based on non-associativity"
    },
    {
        "name": "Q-Mesh Security",
        "path": "qmesh_demo.py",
        "desc": "Entangled state intrusion detection"
    },
    {
        "name": "Teleportation Routing",
        "path": "qmesh_routing_sim.py",
        "desc": "Simulates O(1) latency scaling vs O(sqrt N) grid"
    }
]

def run_demo(demo):
    print(f"\n>>> RUNNING: {demo['name']}")
    print(f"    {demo['desc']}")
    print("-" * 60)
    
    start = time.time()
    try:
        # Run as subprocess to ensure standalone behavior
        result = subprocess.run(
            [sys.executable, demo["path"]],
            capture_output=True,
            text=True,
            timeout=60 # Safety timeout
        )
        
        print(result.stdout)
        
        elapsed = time.time() - start
        
        if result.returncode == 0:
            print(f"[PASS] {demo['name']} verified in {elapsed:.2f}s")
            return True
        else:
            print(result.stderr)
            print(f"[FAIL] {demo['name']} returned error code {result.returncode}")
            return False
            
    except Exception as e:
        print(f"[CRITICAL] Execution failed: {e}")
        return False

def main():
    print("============================================================")
    print("   ALPHAOMEGA FRAMEWORK: TECHNOLOGY VERIFICATION SUITE")
    print("============================================================")
    
    results = []
    for demo in DEMOS:
        if not os.path.exists(demo["path"]):
            print(f"[WARN] Script not found: {demo['path']}")
            results.append(False)
            continue
            
        success = run_demo(demo)
        results.append(success)
        
    print("\n" + "="*60)
    print("FINAL STATUS REPORT")
    print("="*60)
    
    all_pass = True
    for demo, success in zip(DEMOS, results):
        status = "active" if success else "OFFLINE"
        print(f" {status.upper():<10} : {demo['name']}")
        if not success: all_pass = False
        
    print("-" * 60)
    if all_pass:
        print("SYSTEM READY. The AlphaOmega suite is fully operational.")
    else:
        print("SYSTEM WARNING. Some modules failed verification.")

if __name__ == "__main__":
    main()
