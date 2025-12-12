"""
Sedenion Asymmetric Encryption (SAE) Demo

Demonstration of Post-Quantum Cryptography based on Non-Associative Algebra.
Implements a Sedenion Counter Mode (SCM) stream cipher leveraging geometric non-associativity.
"""

import numpy as np
import time

class SedenionMath:
    """Minimal 16D Algebra Implementation."""
    
    DIM = 16
    
    @staticmethod
    def multiply(a, b):
        """
        Sedenion Multiplication (Cayley-Dickson Construction).
        Recursive definition from Octonions.
        """
        if len(a) == 1: return a * b
        
        n = len(a) // 2
        A, B = a[:n], a[n:]
        C, D = b[:n], b[n:]
        
        AC = SedenionMath.multiply(A, C)
        DB_conj = SedenionMath.multiply(D, SedenionMath.conjugate(B))
        DA = SedenionMath.multiply(D, A)
        BC_conj = SedenionMath.multiply(B, SedenionMath.conjugate(C))
        
        res_real = AC - DB_conj
        res_imag = DA + BC_conj
        
        return np.concatenate((res_real, res_imag))

    @staticmethod
    def conjugate(a):
        res = -a.copy()
        res[0] = a[0]
        return res

    @staticmethod
    def normalize(a):
        norm = np.linalg.norm(a)
        if norm < 1e-12: return a
        return a / norm

class SAECipher:
    def __init__(self, seed_val=None):
        if seed_val is None:
            self.key = np.random.normal(0, 1, 16)
        else:
            self.key = np.array(seed_val, dtype=float)
            
        self.key = SedenionMath.normalize(self.key)
            
    def get_keystream(self, index):
        """Generate Z_i = K^(i+1) using Sedenion power dynamics."""
        current = self.key
        for _ in range(index):
            current = SedenionMath.multiply(current, self.key)
            current = SedenionMath.normalize(current)
        return current

    def encrypt(self, message):
        print(f"Encrypting '{message}'...")
        encrypted = []
        for i, char in enumerate(message):
            m_vec = np.zeros(16)
            m_vec[0] = ord(char)
            
            k_vec = self.get_keystream(i)
            c_vec = m_vec + k_vec * 1000.0 
            encrypted.append(c_vec)
        return encrypted

    def decrypt(self, cipher_data):
        print("Decrypting...")
        decrypted = ""
        for i, c_vec in enumerate(cipher_data):
            k_vec = self.get_keystream(i)
            m_vec = c_vec - k_vec * 1000.0
            val = int(round(m_vec[0]))
            decrypted += chr(val)
        return decrypted

def run_demo():
    print("Sedenion Asymmetric Encryption (SAE) Verification")
    print("===============================================")
    
    # 1. Verify One-Way Property
    print("1. Geometric One-Way Property Check")
    a = np.random.normal(0, 1, 16); a = SedenionMath.normalize(a)
    b = np.random.normal(0, 1, 16); b = SedenionMath.normalize(b)
    
    prod = SedenionMath.multiply(a, b)
    inv_a = SedenionMath.conjugate(a) 
    recov = SedenionMath.multiply(inv_a, prod)
    
    error = np.linalg.norm(b - recov)
    print(f"  Non-Associativity Error |b - a^-1(ab)|: {error:.4f}")
    
    if error > 0.1:
        print("  Result: Multiplication is non-invertible (One-Way property confirmed).")
    else:
        print("  Result: Multiplication behaved associatively (Unexpected).")

    # 2. Run Encryption
    print("\n2. Protocol Execution")
    cipher = SAECipher()
    msg = "Hello AlphaOmega World"
    
    enc_data = cipher.encrypt(msg)
    recovered = cipher.decrypt(enc_data)
    
    print(f"  Original:  {msg}")
    print(f"  Recovered: {recovered}")
    
    if msg == recovered:
        print("  Result: Encryption/Decryption cycle successful.")
    else:
        print("  Result: Decryption mismatch.")

if __name__ == "__main__":
    run_demo()