"""
ALPHAOMEGA DEMO: Q-Mesh Network Simulation
Geometric Entanglement Signaling and Intrusion Detection.

Standalone Reference Implementation.
Run: python3 qmesh_demo.py

Theory:
Packets are replaced by shared Geometric States (Psi).
Communication is non-local state rotation.
Interception (Measurement) collapses the G2 invariant, adding noise.
"""

import numpy as np
import time

# --- MINIMAL GEOMETRY ENGINE ---
class GeometricState:
    def __init__(self):
        self.data = np.zeros(16)
        self.data[0] = 1.0 # Start at Identity (e0)
        # We use e1, e2 for signaling
        self.data[1] = 1.0 
        self.data[2] = 0.0
        self.data /= np.linalg.norm(self.data)
        self.coherence = 1.0 # 1.0 = Perfect G2 alignment
        
    def rotate(self, bit):
        """Apply rotation to encode bit."""
        # Rotation in e1-e2 plane
        theta = np.pi/2 if bit else -np.pi/2
        c, s = np.cos(theta), np.sin(theta)
        # Toy rotation
        x, y = self.data[1], self.data[2]
        self.data[1] = x*c - y*s
        self.data[2] = x*s + y*c
        
    def measure(self):
        """Read bit from state."""
        # Check phase in e1-e2 plane
        angle = np.arctan2(self.data[2], self.data[1])
        return 1 if angle > 0 else 0
        
    def collapse(self):
        """Simulate interception collapse."""
        self.coherence *= 0.5
        noise = np.random.normal(0, 1.0 - self.coherence, 16)
        self.data += noise
        self.data /= np.linalg.norm(self.data)

# --- Q-MESH PROTOCOL ---

def run_demo():
    print("--- ALPHAOMEGA Q-MESH DEMO ---")
    
    # 1. Setup Entanglement
    print("\n[1] Establishing Geometric Link...")
    alice_state = GeometricState()
    bob_state = alice_state # Shared reference simulates perfect entanglement
    
    message = [1, 0, 1, 1, 0, 1] * 10 # 60 bits for statistical significance
    print(f"  Alice Sending: {message}")
    
    received = []
    for bit in message:
        # Reset
        alice_state.data = np.zeros(16)
        alice_state.data[0] = 1.0; alice_state.data[1] = 1.0; alice_state.data[2] = 0.0
        alice_state.data /= np.linalg.norm(alice_state.data)
        
        alice_state.rotate(bit)
        rec = bob_state.measure()
        received.append(rec)
        
    # print(f"  Bob Received:  {received}") # Too long
    assert message == received
    print("  [PASS] Transmission successful.")
    
    # 2. Intrusion Test
    print("\n[2] Simulating Man-in-the-Middle (Eve)...")
    
    # Re-key
    alice_state = GeometricState()
    bob_state = alice_state
    
    rec_corrupted = []
    detected_errors = 0
    
    for bit in message:
        # Reset
        alice_state.data = np.zeros(16)
        alice_state.data[0] = 1.0; alice_state.data[1] = 1.0; alice_state.data[2] = 0.0
        alice_state.data /= np.linalg.norm(alice_state.data)
        
        alice_state.rotate(bit)
        
        # EVE INTERCEPTS
        bob_state.collapse()
        
        rec = bob_state.measure()
        rec_corrupted.append(rec)
        
        if rec != bit:
            detected_errors += 1
        
    print(f"  Bit Error Rate: {detected_errors/len(message)*100:.1f}%")
    
    if detected_errors > 0:
        print("  [SUCCESS] Intrusion Detected. Geometric collapse confirmed.")
        print("  Connection would be severed immediately.")
        return True
    else:
        print("  [FAIL] Eve went undetected.")
        return False

import sys
if __name__ == "__main__":
    if not run_demo():
        sys.exit(1)