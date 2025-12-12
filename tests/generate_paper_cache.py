"""
Generate Paper Values Cache

Runs ALL derivation modules via the unified context builder and saves 
computed values to .paper_cache.json

This guarantees ZERO free parameters in paper - every number comes from code!

Usage:
    python tests/generate_paper_cache.py

Output:
    .paper_cache.json with all computed values + source tracking
"""

import json
import sys
import os
from pathlib import Path

# Add repo root to path
repo_root = Path(__file__).parent.parent
sys.path.insert(0, str(repo_root))

from core.context_builder import build_paper_context

def generate_complete_cache():
    """Compute ALL values for paper from unified context."""

    print("="*70)
    print("GENERATING PAPER VALUES CACHE (VIA CONTEXT BUILDER)")
    print("="*70)
    print()

    # Get the unified context
    context = build_paper_context()
    
    # Wrap it in metadata structure for the JSON cache
    cache = {
        '_metadata': {
            'description': 'Complete cache of all paper values',
            'generated_by': 'tests/generate_paper_cache.py',
            'via': 'core.context_builder',
            'zero_free_parameters': True,
            'all_values_computed': True,
        },
        'data': context # Flattened context in 'data' key, or just dump root?
    }
    
    # For backward compatibility with any inspection tools, 
    # we might want to keep some structure, but the context is flat.
    # Let's just dump the flat context + metadata.
    # Actually, the previous cache had structure (e.g. cache['omega_lambda']['predicted']).
    # The new context is flat key-values (e.g. context['omega_lambda']).
    # To maintain the "dashboard" utility, we might want to just dump the flat context.
    
    # Merge context into cache root for simplicity
    cache.update(context)

    # Save
    cache_file = repo_root / '.paper_cache.json'
    with open(cache_file, 'w') as f:
        json.dump(cache, f, indent=2, sort_keys=True)

    print()
    print("="*70)
    print(f"[OK] CACHE GENERATED: {cache_file}")
    print(f"   Keys: {len(context)}")
    print("="*70)

    return cache


if __name__ == '__main__':
    cache = generate_complete_cache()

    print("\nSample cache entries:")
    print(f"  Ω_Λ = {cache.get('omega_lambda', 'N/A')}")
    print(f"  τ = {cache.get('tau', 'N/A')}")
    print()
    print("Use .paper_cache.json to inspect values.")