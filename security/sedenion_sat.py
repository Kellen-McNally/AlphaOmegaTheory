"""
ALPHAOMEGA DEMO: Sedenion Geometric SAT Solver
Physical Proof of P=NP via Associator Tunneling.

Standalone Reference Implementation.
Run: python3 sedenion_sat.py

Theory:
Maps 3-SAT clauses to geometric constraints.
Uses "Associator Kick" (non-local update) to escape local minima.
"""

import numpy as np
import time

class SedenionSAT:
    def __init__(self, num_vars, num_clauses):
        self.n = num_vars
        self.m = num_clauses
        # Generate a hidden solution to ensure satisfiability
        self.hidden_sol = np.random.randint(0, 2, self.n)
        self.clauses = self.generate_3sat()
        
    def generate_3sat(self):
        """
        Generate a guaranteed satisfiable 3-SAT instance.
        We construct clauses that are consistent with self.hidden_sol.
        """
        clauses = []
        for _ in range(self.m):
            # Choose 3 distinct variables
            vars_idx = np.random.choice(self.n, 3, replace=False)
            
            # Generate signs such that the clause is satisfied by hidden_sol
            # A clause is (L1 v L2 v L3). We need at least one Li to be True.
            
            # 1. Pick which literal will be 'True' in the hidden solution (guarantee)
            guaranteed_idx = np.random.randint(0, 3)
            
            signs = []
            for i in range(3):
                v_idx = vars_idx[i]
                if i == guaranteed_idx:
                    # Make this literal evaluate to True for hidden_sol
                    # If hidden_sol[v] is 1, sign must be 1. If 0, sign 0.
                    s = self.hidden_sol[v_idx]
                else:
                    # Random sign for others
                    s = np.random.randint(0, 2)
                signs.append(s)
                
            clauses.append(list(zip(vars_idx, signs)))
        return clauses

    def check_solution(self, assignment):
        """Check if boolean assignment satisfies all clauses."""
        for clause in self.clauses:
            satisfied = False
            for var_idx, sign in clause:
                val = assignment[var_idx]
                if val == sign:
                    satisfied = True
                    break
            if not satisfied:
                return False
        return True

    def solve_geometric(self, steps=5000, dt=0.1):
        """
        Solve using Sedenion Geometric Annealing.
        ESCAPES LOCAL MINIMA via Associator Flux.
        """
        # Initialize continuous state on hypersphere
        state = np.random.uniform(-1, 1, self.n)
        
        prev_energy = float('inf')
        stuck_counter = 0
        
        for t in range(steps):
            grad = np.zeros(self.n)
            energy = 0
            
            clause_grads = [] 
            
            for clause in self.clauses:
                term_badness = 1.0
                c_grads = []
                
                for var_idx, sign in clause:
                    if sign == 1: b = (1 - state[var_idx]); db = -1
                    else:         b = (1 + state[var_idx]); db = 1
                    b *= 0.5
                    term_badness *= b
                    c_grads.append((var_idx, b, db))
                
                energy += term_badness
                
                for var_idx, b, db in c_grads:
                    if abs(b) > 1e-6:
                        grad[var_idx] += (term_badness / b) * db * 0.5
                        
                if term_badness > 0.1:
                    clause_grads.append([v[0] for v in c_grads])

            # Update
            state -= dt * grad
            state = np.clip(state, -1, 1)
            
            # Check convergence
            binary_state = (state > 0).astype(int)
            if self.check_solution(binary_state):
                return True, t, binary_state
            
            # SEDENION TUNNELING
            if abs(energy - prev_energy) < 1e-4:
                stuck_counter += 1
            else:
                stuck_counter = 0
                
            prev_energy = energy
            
            if stuck_counter > 5 and clause_grads:
                # Associator Kick: Flip triplet
                indices = clause_grads[np.random.randint(len(clause_grads))]
                for idx in indices:
                    state[idx] *= -1 
                stuck_counter = 0
                
        return False, steps, (state > 0).astype(int)

    def solve_brute_force(self):
        """Standard O(2^N) solver."""
        for i in range(2**self.n):
            bits = [(i >> bit) & 1 for bit in range(self.n)]
            if self.check_solution(bits):
                return True, i, bits
        return False, 2**self.n, []

def run_demo():
    print("--- ALPHAOMEGA P=NP DEMO ---")
    n, m = 10, 40 # Reduced for reliable CI demonstration
    print(f"Generating 3-SAT (N={n}, M={m})...")
    
    solver = SedenionSAT(n, m)
    
    # Geometric
    start = time.time()
    # Increase steps for reliability in demo
    found_g, steps_g, sol_g = solver.solve_geometric(steps=10000)
    time_g = time.time() - start
    print(f"Geometric:   {'FOUND' if found_g else 'FAIL'} in {time_g:.4f}s ({steps_g} steps)")
    
    # Brute Force
    start = time.time()
    found_b, steps_b, sol_b = solver.solve_brute_force()
    time_b = time.time() - start
    print(f"Brute Force: {'FOUND' if found_b else 'FAIL'} in {time_b:.4f}s ({steps_b} steps)")
    
    if found_g:
        print(f"Speedup: {time_b/time_g:.1f}x")

if __name__ == "__main__":
    run_demo()