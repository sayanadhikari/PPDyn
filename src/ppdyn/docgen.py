"""Generate the input reference (Markdown) from the schema."""
from __future__ import annotations

from .schema import SCHEMA, SECTION_ORDER, Table

DIM_TEXT = {"1": "", "SI": "SI", "L": "length", "T": "time", "V": "velocity", "A": "acceleration",
            "E": "E-field", "F": "force", "Q": "charge", "M": "mass", "TEMP": "temperature",
            "RATE": "1/time", "RATE2": "1/time²", "EGRAD": "E-field/length", "RHO": "density"}


def _fmt_default(v):
    if v is None:
        return "auto"
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, str):
        return f'"{v}"'
    return str(v)


def _table(lines, name, t: Table, level=2):
    bracket = f"[[{name}]]" if t.array else f"[{name}]"
    lines.append(f"{'#' * level} `{bracket}`\n")
    lines.append(t.doc + "\n")
    if t.keys:
        lines.append("| key | default | unit | description |")
        lines.append("|---|---|---|---|")
        for k, key in t.keys.items():
            desc = key.doc
            if key.choices:
                desc += " Choices: " + ", ".join(f"`{c}`" for c in key.choices) + "."
            if key.required:
                desc = "**Required.** " + desc
            lines.append(f"| `{k}` | `{_fmt_default(key.default)}` | {DIM_TEXT.get(key.dim, key.dim)} | {desc} |")
        lines.append("")
    for sub_name, sub in t.sub.items():
        _table(lines, f"{name}.{sub_name}", sub, level + 1)


def markdown() -> str:
    lines = ["# PPDyn input reference\n",
             "Generated from `src/ppdyn/schema.py` with `ppdyn docs`. Unknown keys are rejected with a suggestion, "
             "so typos never pass silently.\n",
             "The *unit* column gives the physical dimension. With `units.system = \"reduced\"` values are in "
             "reduced units (length a, time 1/ω₀, charge Q₀, mass m₀, temperature Q₀²/(4πε₀a) = 1/Γ). "
             "With `units.system = \"SI\"` the same keys take SI values (m, s, C, kg, K, V/m, N) and are converted "
             "internally. Keys marked *SI* are always SI.\n"]
    for name in SECTION_ORDER:
        _table(lines, name, SCHEMA[name])
    return "\n".join(lines)
