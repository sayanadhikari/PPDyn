"""Reduced (dust-plasma) units and SI conversion.

Reduced units
-------------
length  a        (Wigner-Seitz radius by default)
charge  Q0
mass    m0
time    1/omega0,   omega0^2 = Q0^2 / (4 pi eps0 m0 a^3)

Derived: force F0 = Q0^2/(4 pi eps0 a^2), field E0 = F0/Q0, energy
F0 a, so the reduced temperature is k_B T / (F0 a) = 1/Gamma and the
Yukawa force between two reference particles is (1 + kappa r) e^(-kappa r) / r^2.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, asdict

EPS0 = 8.8541878128e-12
KB = 1.380649e-23
QE = 1.602176634e-19
AMU = 1.66053906660e-27


@dataclass
class Scales:
    length: float = 1.0    # a  [m]
    charge: float = 1.0    # Q0 [C]
    mass: float = 1.0      # m0 [kg]
    si: bool = False

    @property
    def omega(self) -> float:
        if not self.si:
            return 1.0
        return math.sqrt(self.charge ** 2 / (4 * math.pi * EPS0 * self.mass * self.length ** 3))

    @property
    def time(self) -> float:
        return 1.0 / self.omega

    @property
    def force(self) -> float:
        return self.mass * self.length * self.omega ** 2

    @property
    def efield(self) -> float:
        return self.force / self.charge

    @property
    def energy(self) -> float:
        return self.force * self.length

    @property
    def temperature(self) -> float:
        """Kelvin corresponding to one reduced temperature unit."""
        return self.energy / KB if self.si else 1.0

    def factor(self, dim: str) -> float:
        """Divide an input value of dimension ``dim`` by this to get reduced units."""
        if not self.si:
            return 1.0
        a, w = self.length, self.omega
        return {
            "1": 1.0, "SI": 1.0,
            "L": a, "T": 1.0 / w, "V": a * w, "A": a * w * w,
            "E": self.efield, "F": self.force, "Q": self.charge, "M": self.mass,
            "TEMP": self.temperature, "RATE": w, "RATE2": w * w,
            "EGRAD": self.efield / a, "RHO": self.mass / a ** 3,
        }[dim]

    def to_reduced(self, value, dim: str):
        f = self.factor(dim)
        if f == 1.0 or value is None:
            return value
        if isinstance(value, list):
            return [v / f for v in value]
        return value / f

    def as_dict(self) -> dict:
        d = asdict(self)
        d.update(omega=self.omega, time=self.time, force=self.force,
                 efield=self.efield, energy=self.energy, temperature=self.temperature)
        return d

    def describe(self) -> str:
        if not self.si:
            return "reduced units (a = Q0 = m0 = 1)"
        return (f"a = {self.length:.4g} m, Q0 = {self.charge:.4g} C ({self.charge / QE:.4g} e), "
                f"m0 = {self.mass:.4g} kg, omega0 = {self.omega:.4g} rad/s, "
                f"T unit = {self.temperature:.4g} K")


def gamma(T_reduced: float) -> float:
    """Coupling parameter for a reduced temperature."""
    return math.inf if T_reduced <= 0 else 1.0 / T_reduced


def convert_tree(node, schema_table, scales: Scales):
    """Convert a completed config table (in place) from input units to reduced units."""
    if node is None:
        return
    for name, key in schema_table.keys.items():
        if name in node and key.dim not in ("1", "SI"):
            v = node[name]
            if key.kind == "list" and isinstance(v, list):
                node[name] = [x / scales.factor(key.dim) for x in v]
            elif isinstance(v, (int, float, list)) and not isinstance(v, bool):
                node[name] = scales.to_reduced(v, key.dim)
    for name, sub in schema_table.sub.items():
        child = node.get(name)
        if child is None:
            continue
        if sub.array:
            for item in child:
                convert_tree(item, sub, scales)
        else:
            convert_tree(child, sub, scales)
