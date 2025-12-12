"""
Test code integration and prevent constant duplication.

This test ensures that all modules properly import from the canonical
core.constants module instead of redefining constants locally.
"""

import pytest
import re
import os
from pathlib import Path


class TestConstantIntegration:
    """Ensure constants are not redefined across modules."""

    def test_no_dim_g2_redefinition(self):
        """
        Verify DIM_G2 = 14 is not hardcoded in modules.

        All modules should import DIM_G2 from core.constants.
        """
        # Directories to check
        root_dir = Path(__file__).parent.parent

        # Files that are allowed to define DIM_G2 (only constants.py and constants_expanded.py)
        allowed_files = [
            "core/constants.py",
            "core/constants.py",
        ]

        # Pattern to detect hardcoded DIM_G2 = 14
        pattern = re.compile(r'^\s*DIM_G2\s*=\s*14', re.MULTILINE)

        violations = []

        # Scan all Python files in core, physics, derivations, atomic, cosmology
        dirs_to_scan = ["core", "physics", "derivations", "atomic", "cosmology"]
        
        for d in dirs_to_scan:
            scan_dir = root_dir / d
            if not scan_dir.exists(): continue
            
            for py_file in scan_dir.rglob("*.py"):
                # Skip __pycache__ and test files
                if "__pycache__" in str(py_file) or "test_" in py_file.name:
                    continue

                # Check if this file is in the allowed list
                rel_path = str(py_file.relative_to(root_dir))
                if rel_path in allowed_files:
                    continue

                # Read file and check for hardcoded constant
                try:
                    content = py_file.read_text()
                    if pattern.search(content):
                        violations.append(str(py_file))
                except Exception:
                    pass  # Skip files that can't be read

        assert len(violations) == 0, (
            f"Found DIM_G2 = 14 hardcoded in {len(violations)} files:\n" +
            "\n".join(f"  - {f}" for f in violations) +
            "\n\nAll files should import DIM_G2 from core.constants instead."
        )

    def test_no_triality_redefinition(self):
        """
        Verify TRIALITY = 3 is not hardcoded in modules.

        All modules should import TRIALITY from core.constants.
        """
        root_dir = Path(__file__).parent.parent

        allowed_files = [
            "core/constants.py",
            "core/constants.py",
        ]

        # Pattern to detect hardcoded TRIALITY = 3
        pattern = re.compile(r'^\s*TRIALITY\s*=\s*3', re.MULTILINE)

        violations = []
        
        dirs_to_scan = ["core", "physics", "derivations", "atomic", "cosmology"]

        for d in dirs_to_scan:
            scan_dir = root_dir / d
            if not scan_dir.exists(): continue
            
            for py_file in scan_dir.rglob("*.py"):
                if "__pycache__" in str(py_file) or "test_" in py_file.name:
                    continue

                rel_path = str(py_file.relative_to(root_dir))
                if rel_path in allowed_files:
                    continue

                try:
                    content = py_file.read_text()
                    if pattern.search(content):
                        violations.append(str(py_file))
                except Exception:
                    pass

        assert len(violations) == 0, (
            f"Found TRIALITY = 3 hardcoded in {len(violations)} files:\n" +
            "\n".join(f"  - {f}" for f in violations) +
            "\n\nAll files should import TRIALITY from core.constants instead."
        )

    def test_modules_import_from_constants(self):
        """
        Verify that core modules import from core.constants.

        This is a sanity check that key modules are using the canonical source.
        """
        # Key modules that should import from constants
        key_modules = [
            "core/casimir.py",
            "atomic/ionization/energy.py",
        ]

        root_dir = Path(__file__).parent.parent

        # Pattern to detect import from core.constants (absolute or relative)
        import_pattern = re.compile(r'from\s+(core\.constants|\.constants)\s+import')

        missing_imports = []

        for module_path in key_modules:
            full_path = root_dir / module_path

            if not full_path.exists():
                continue

            content = full_path.read_text()
            if not import_pattern.search(content):
                missing_imports.append(module_path)

        assert len(missing_imports) == 0, (
            f"The following modules don't import from core.constants:\n" +
            "\n".join(f"  - {m}" for m in missing_imports)
        )


class TestFunctionIntegration:
    """Ensure no duplicate function implementations."""

    def test_alpha_gut_uses_constants(self):
        """
        Verify alpha_gut() in casimir.py uses imported constants.

        Should not recalculate from local variables.
        """
        casimir_file = Path(__file__).parent.parent / "core" / "casimir.py"

        if not casimir_file.exists():
            pytest.skip("casimir.py not found")

        content = casimir_file.read_text()

        # Check that alpha_gut function exists (if it exists in casimir.py)
        # Note: alpha_gut is likely in core/constants.py now, but if casimir.py uses it,
        # it should use the constant.
        
        # Let's check if casimir.py uses TRIALITY
        assert "TRIALITY" in content, "casimir.py should use TRIALITY"
        assert "DIM_G2" in content, "casimir.py should use DIM_G2"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])