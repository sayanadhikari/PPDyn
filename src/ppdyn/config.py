"""Load, validate and complete a PPDyn TOML input file."""
from __future__ import annotations

import copy
import difflib
import os
import sys
from typing import Any

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover
    import tomli as tomllib

from .schema import SCHEMA, Key, Table


class ConfigError(ValueError):
    """Raised for invalid input files, with a readable message."""


def _suggest(name: str, options) -> str:
    close = difflib.get_close_matches(name, list(options), n=1)
    return f" Did you mean '{close[0]}'?" if close else ""


def _check_value(path: str, key: Key, value: Any) -> Any:
    k = key.kind
    if value is None:
        return None
    if k == "bool":
        if not isinstance(value, bool):
            raise ConfigError(f"{path} must be true or false, got {value!r}")
    elif k == "int":
        if isinstance(value, bool) or not isinstance(value, int):
            raise ConfigError(f"{path} must be an integer, got {value!r}")
    elif k == "float":
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ConfigError(f"{path} must be a number, got {value!r}")
        value = float(value)
    elif k == "str":
        if not isinstance(value, str):
            raise ConfigError(f"{path} must be a string, got {value!r}")
        if key.choices and value not in key.choices:
            raise ConfigError(f"{path} = {value!r} is not one of {list(key.choices)}.{_suggest(value, key.choices)}")
    elif k == "vec3":
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            value = [float(value)] * 3
        if not isinstance(value, list) or len(value) != 3:
            raise ConfigError(f"{path} must be a list of three numbers, got {value!r}")
        value = [float(v) for v in value]
    elif k == "list":
        if not isinstance(value, list):
            value = [value]
    elif k == "bc":
        vals = value if isinstance(value, list) else [value, value]
        if len(vals) != 2:
            raise ConfigError(f"{path} must be a boundary type or [low, high]")
        for v in vals:
            if v not in key.choices:
                raise ConfigError(f"{path}: {v!r} is not one of {list(key.choices)}.{_suggest(str(v), key.choices)}")
        if ("periodic" in vals) and vals[0] != vals[1]:
            raise ConfigError(f"{path}: periodic must be used on both sides")
        value = list(vals)
    return value


def _fill_table(path: str, table: Table, data: dict) -> dict:
    if not isinstance(data, dict):
        raise ConfigError(f"[{path}] must be a table")
    allowed = set(table.keys) | set(table.sub)
    for name in data:
        if name not in allowed:
            raise ConfigError(f"Unknown key '{name}' in [{path}].{_suggest(name, allowed)}")
    out = {}
    for name, key in table.keys.items():
        if name in data:
            out[name] = _check_value(f"{path}.{name}", key, data[name])
        else:
            if key.required:
                raise ConfigError(f"Missing required key '{name}' in [{path}]")
            out[name] = _check_value(f"{path}.{name}", key, copy.deepcopy(key.default))
    for name, sub in table.sub.items():
        sp = f"{path}.{name}"
        if sub.array:
            items = data.get(name, [])
            if isinstance(items, dict):
                items = [items]
            out[name] = [_fill_table(sp, sub, it) for it in items]
        elif name in data:
            out[name] = _fill_table(sp, sub, data[name])
            out[name]["enabled"] = True
        elif table is SCHEMA.get("forces"):
            out[name] = None          # force switched off
        else:
            out[name] = _fill_table(sp, sub, {})
    return out


def complete(raw: dict) -> dict:
    """Validate a raw input dictionary and fill every default."""
    for name in raw:
        if name not in SCHEMA:
            raise ConfigError(f"Unknown section [{name}].{_suggest(name, SCHEMA)}")
    cfg = {}
    for name, table in SCHEMA.items():
        if table.array:
            items = raw.get(name, [])
            if isinstance(items, dict):
                items = [items]
            if not items and not table.optional:
                raise ConfigError(f"At least one [[{name}]] table is required")
            cfg[name] = [_fill_table(name, table, it) for it in items]
        elif name == "field":
            cfg[name] = _fill_table(name, table, raw[name]) if name in raw else None
        else:
            cfg[name] = _fill_table(name, table, raw.get(name, {}))
    _cross_checks(cfg)
    return cfg


def _cross_checks(cfg: dict) -> None:
    si = cfg["units"]["system"] == "SI"
    if cfg["time"]["dt"] <= 0:
        raise ConfigError("time.dt must be positive")
    if cfg["domain"]["shape"] == "sphere" and cfg["domain"]["radius"] <= 0:
        raise ConfigError("domain.radius must be > 0 for shape = \"sphere\"")
    if cfg["neutrals"]["drag"] == "epstein":
        if not si:
            raise ConfigError("neutrals.drag = \"epstein\" needs units.system = \"SI\" (it uses pressure, gas mass and grain radius)")
        if cfg["neutrals"]["pressure"] <= 0:
            raise ConfigError("neutrals.pressure must be > 0 for Epstein drag")
    for i, sp in enumerate(cfg["species"]):
        if sp["charge"]["model"] == "potential" and not si:
            raise ConfigError(f"species[{i}].charge.model = \"potential\" needs SI units")
        if sp["count"] < 0:
            raise ConfigError(f"species[{i}].count must be >= 0")
    if cfg["output"]["every"] <= 0 or cfg["output"]["history_every"] <= 0:
        raise ConfigError("output.every and output.history_every must be positive")


def load(path: str) -> dict:
    """Read a TOML input file and return the completed configuration."""
    with open(path, "rb") as fh:
        try:
            raw = tomllib.load(fh)
        except tomllib.TOMLDecodeError as exc:
            raise ConfigError(f"{path}: {exc}") from None
    cfg = complete(raw)
    cfg["_source"] = os.path.abspath(path)
    with open(path, "r", encoding="utf-8") as fh:
        cfg["_text"] = fh.read()
    return cfg


def loads(text: str) -> dict:
    cfg = complete(tomllib.loads(text))
    cfg["_source"] = "<string>"
    cfg["_text"] = text
    return cfg
