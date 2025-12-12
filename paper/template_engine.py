"""
Minimal Template Engine

A lightweight replacement for Jinja2 to support LaTeX generation
without external dependencies.

Features:
- Variable substitution: {{ variable_name }}
- Simple formatting: {{ variable_name:.2f }}
"""

import re
from typing import Dict, Any, Optional
from pathlib import Path


class TemplateEngine:
    """A simple template engine for substituting variables in text files."""

    def __init__(self, template_dir: str) -> None:
        """
        Initializes the template engine.

        Args:
            template_dir: The directory where templates are located.
        """
        self.template_dir = Path(template_dir)

    def render(self, template_name: str, context: Dict[str, Any]) -> str:
        """
        Render a template with the given context.

        Args:
            template_name: The name of the template file.
            context: A dictionary of variables to substitute into the template.

        Returns:
            The rendered template as a string.

        Raises:
            FileNotFoundError: If the template file cannot be found.
        """
        template_path = self.template_dir / template_name
        if not template_path.exists():
            raise FileNotFoundError(f"Template not found: {template_path}")

        content = template_path.read_text(encoding="utf-8")

        # 1. Handle Variable Substitution with optional formatting
        # Pattern: {{ variable_name }} or {{ variable_name:.2f }}
        def replace_var(match: re.Match) -> str:
            full_expr = match.group(1).strip()

            # Check for formatting
            if ":" in full_expr:
                var_name, format_spec = full_expr.split(":", 1)
                var_name = var_name.strip()  # Ensure no leading/trailing spaces
                format_spec = format_spec.strip()  # Ensure no leading/trailing spaces
                if var_name in context:
                    val = context[var_name]
                    try:
                        replaced_val = f"{val:{format_spec}}"
                        return replaced_val
                    except ValueError:
                        return str(val)
            else:
                var_name = full_expr
                if var_name in context:
                    replaced_val = str(context[var_name])
                    return replaced_val

            # If variable not found, return original string (or could raise error)
            return match.group(0)

        # Regex for {{ ... }}
        content = re.sub(r"\{\{(.*?)\}\}", replace_var, content)

        return content


# Global instance
_engine: Optional[TemplateEngine] = None


def render_template(template_name: str, context: Dict[str, Any]) -> str:
    """
    Renders a template using a global engine instance.

    Initializes a singleton TemplateEngine on first use.

    Args:
        template_name: The name of the template file.
        context: A dictionary of variables to substitute into the template.

    Returns:
        The rendered template as a string.
    """
    global _engine
    if _engine is None:
        # Assume templates are in paper/templates/ relative to project root
        # Adjust path logic as needed
        # If run from project root, __file__ is paper/template_engine.py
        root = Path(__file__).parent.parent
        _engine = TemplateEngine(root / "paper" / "templates")

    return _engine.render(template_name, context)
