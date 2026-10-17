"""Input-file schema for PPDyn.

Every recognised key is declared here once, with its default value, type,
physical dimension and a one-line description. The same table drives
validation, unit conversion and the generated input reference
(``ppdyn docs`` / ``docs/input_reference.md``).

Dimension codes
---------------
L      length            T      time             V      velocity
A      acceleration      E      electric field   F      force
Q      charge            M      mass             TEMP   temperature
RATE   1/time            RATE2  1/time^2         EGRAD  electric field / length
RHO    mass density      1      dimensionless    SI     SI only, never scaled
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Sequence

BC_CHOICES = ("periodic", "reflecting", "absorbing", "open")
SCALE_CHOICES = ("none", "charge", "mass", "radius", "radius2")


@dataclass
class Key:
    default: Any
    kind: str                      # str|int|float|bool|vec3|list|bc|any
    doc: str
    dim: str = "1"
    choices: Sequence[str] | None = None
    required: bool = False


@dataclass
class Table:
    """A TOML table. ``array`` means it is an array of tables ([[name]])."""
    keys: dict
    doc: str
    array: bool = False
    optional: bool = True          # whole table may be absent
    sub: dict = field(default_factory=dict)


# ---------------------------------------------------------------- species
SPECIES_RADIUS = Table(doc="Particle size distribution (radius, or mass when variable = \"mass\").", keys={
    "distribution": Key("fixed", "str", "Distribution type.", choices=(
        "fixed", "normal", "lognormal", "uniform", "bimodal", "binary")),
    "variable": Key("radius", "str", "Quantity the distribution is drawn for. \"mass\" reproduces the legacy dist = True option (m ~ N(mean, std), r ~ m^(1/3)).", choices=("radius", "mass")),
    "mean": Key(0.0, "float", "Mean (fixed value for \"fixed\"). For variable = \"mass\" it is in units of species.mass.", "L"),
    "std": Key(0.0, "float", "Standard deviation (normal, lognormal: of ln r).", "L"),
    "low": Key(0.0, "float", "Lower bound (uniform).", "L"),
    "high": Key(0.0, "float", "Upper bound (uniform).", "L"),
    "values": Key([], "list", "Radii for binary/bimodal (two entries).", "L"),
    "stds": Key([], "list", "Standard deviations for bimodal (two entries).", "L"),
    "weights": Key([0.5, 0.5], "list", "Population fractions for binary/bimodal."),
    "min": Key(0.0, "float", "Samples below this are redrawn (keeps radii positive).", "L"),
})

SPECIES_CHARGE = Table(doc="Charge model. Q = value * (r / reference_radius)^exponent, or Q = 4 pi eps0 r phi.", keys={
    "model": Key("power", "str", "\"power\" (fixed when exponent = 0) or \"potential\" (SI only, OML-like Q = 4 pi eps0 r phi).", choices=("power", "potential")),
    "value": Key(1.0, "float", "Charge at the reference radius.", "Q"),
    "exponent": Key(1.0, "float", "Size exponent p (1 = constant surface potential, 0 = equal charges, 2 = legacy dist).",),
    "potential": Key(-2.0, "float", "Surface potential in volts for model = \"potential\".", "SI"),
})

SPECIES_PLACE = Table(doc="Initial placement.", keys={
    "region": Key("box", "str", "Where particles are loaded.", choices=("box", "sphere", "shell", "plane", "lattice", "file")),
    "center": Key(None, "vec3", "Region centre (defaults to domain centre).", "L"),
    "size": Key(None, "vec3", "Box size for region = box/lattice (defaults to the domain).", "L"),
    "r_min": Key(0.0, "float", "Inner radius for sphere/shell.", "L"),
    "r_max": Key(0.0, "float", "Outer radius for sphere/shell.", "L"),
    "axis": Key("z", "str", "Normal axis for region = plane.", choices=("x", "y", "z")),
    "position": Key(0.0, "float", "Plane coordinate for region = plane (e.g. top of the box for dust fall).", "L"),
    "lattice": Key("bcc", "str", "Lattice type for region = lattice.", choices=("sc", "bcc", "fcc")),
    "file": Key("", "str", "HDF5/CSV/NPY file with an (N,3) position array for region = file."),
    "jitter": Key(0.0, "float", "Random displacement added to lattice sites.", "L"),
    "avoid_objects": Key(True, "bool", "Redraw particles that land inside objects."),
})

SPECIES_VEL = Table(doc="Initial velocities.", keys={
    "distribution": Key("maxwell", "str", "Velocity loading.", choices=("maxwell", "uniform", "zero")),
    "temperature": Key(0.0, "float", "Temperature for maxwell loading.", "TEMP"),
    "vmax": Key([1.0, 1.0, 1.0], "vec3", "Full width of the uniform distribution per axis.", "V"),
    "drift": Key([0.0, 0.0, 0.0], "vec3", "Mean drift velocity added after loading.", "V"),
    "remove_net_momentum": Key(True, "bool", "Subtract the centre-of-mass velocity before adding drift."),
})

SPECIES = Table(array=True, optional=False, doc="Particle species. Repeat [[species]] for several populations.", keys={
    "name": Key("dust", "str", "Species label."),
    "count": Key(1000, "int", "Number of particles.", required=True),
    "mass": Key(1.0, "float", "Mass at the reference radius (ignored when density is set).", "M"),
    "density": Key(0.0, "float", "Material density; if > 0, m = density * 4/3 pi r^3.", "RHO"),
    "reference_radius": Key(0.0, "float", "Radius at which mass/charge values apply (0 = mean radius).", "L"),
}, sub={"radius": SPECIES_RADIUS, "charge": SPECIES_CHARGE, "placement": SPECIES_PLACE, "velocity": SPECIES_VEL})

# ---------------------------------------------------------------- forces
FORCES = Table(doc="External forces. A sub-table switches a force on.", keys={}, sub={
    "gravity": Table(doc="Uniform gravity, F = m g.", keys={
        "g": Key([0.0, 0.0, -9.81], "vec3", "Gravitational acceleration.", "A")}),
    "efield": Table(doc="Uniform electric field, F = Q E.", keys={
        "E": Key([0.0, 0.0, 0.0], "vec3", "Field vector.", "E")}),
    "linear_efield": Table(doc="Linear (sheath-like) field, E = E0 + G (x - origin) per axis, F = Q E.", keys={
        "E0": Key([0.0, 0.0, 0.0], "vec3", "Field at the origin.", "E"),
        "gradient": Key([0.0, 0.0, 0.0], "vec3", "dE_i/dx_i for each axis.", "EGRAD"),
        "origin": Key([0.0, 0.0, 0.0], "vec3", "Reference point.", "L")}),
    "harmonic": Table(doc="Harmonic confinement acceleration a = -omega2 (x - origin), independent of mass. Legacy reflecting runs used this along z.", keys={
        "omega2": Key([0.0, 0.0, 0.0], "vec3", "Squared trap frequency per axis.", "RATE2"),
        "origin": Key([0.0, 0.0, 0.0], "vec3", "Trap centre.", "L")}),
    "constant": Table(doc="Constant force, e.g. ion flow drag (dust-void f_flow).", keys={
        "F": Key([0.0, 0.0, 0.0], "vec3", "Force on a reference particle.", "F"),
        "scale_with": Key("none", "str", "Per-particle scaling of F.", choices=SCALE_CHOICES)}),
    "void": Table(doc="Radial ion-drag void force F = A exp(-r/lambda) [1 + m sin(w t)] r_hat (dust-void branch).", keys={
        "amplitude": Key(0.0, "float", "Force amplitude A on a reference particle.", "F"),
        "decay_length": Key(1.0, "float", "Decay length lambda.", "L"),
        "center": Key([0.0, 0.0, 0.0], "vec3", "Void centre.", "L"),
        "scale_with": Key("radius2", "str", "Per-particle scaling (ion drag ~ r^2).", choices=SCALE_CHOICES),
        "modulation": Key(0.0, "float", "Relative modulation depth m."),
        "omega": Key(0.0, "float", "Modulation angular frequency w.", "RATE"),
        "cylindrical": Key(False, "bool", "Use distance from the z axis instead of the centre point.")}),
    "radial_table": Table(array=True, doc="Tabulated radial force F(r) r_hat (dust-exp E-field and ion-drag profiles). Repeat for several tables.", keys={
        "file": Key("", "str", "CSV (two columns r, F) or NPZ file.", required=True),
        "r_key": Key("r", "str", "Array name for r inside an NPZ file (e.g. arr_0)."),
        "f_key": Key("F", "str", "Array name for F inside an NPZ file (e.g. arr_1)."),
        "r_unit": Key(1.0, "float", "Multiply r by this to get the input unit system."),
        "factor": Key(1.0, "float", "Multiply F by this (e.g. legacy k_fel)."),
        "kind": Key("force", "str", "Table holds a force on a reference particle or an electric field.", choices=("force", "efield")),
        "center": Key([0.0, 0.0, 0.0], "vec3", "Centre of the radial profile.", "L"),
        "scale_with": Key("none", "str", "Per-particle scaling (ignored for kind = efield, which scales with charge).", choices=SCALE_CHOICES)}),
})

SCHEMA: dict[str, Table] = {
    "run": Table(doc="Run control.", keys={
        "name": Key("ppdyn", "str", "Run name, used in output file names."),
        "threads": Key(0, "int", "Numba threads (0 = all cores)."),
        "seed": Key(12345, "int", "Random seed."),
        "progress": Key(True, "bool", "Show a progress bar."),
        "restart": Key("", "str", "Checkpoint file to restart from."),
    }),
    "units": Table(doc="Unit system. In reduced units lengths are in a, time in 1/omega0, charge in Q0, mass in m0 with omega0^2 = Q0^2 / (4 pi eps0 m0 a^3), so temperature equals 1/Gamma. With system = \"SI\" every dimensional input is SI and is converted internally.", keys={
        "system": Key("reduced", "str", "Input unit system.", choices=("reduced", "SI")),
        "length": Key(0.0, "float", "SI only: reference length a in metres (0 = Wigner-Seitz radius of the loaded density).", "SI"),
        "charge": Key(0.0, "float", "SI only: reference charge Q0 in coulomb (0 = mean |Q|).", "SI"),
        "mass": Key(0.0, "float", "SI only: reference mass m0 in kg (0 = mean mass).", "SI"),
    }),
    "domain": Table(doc="Simulation domain, centred on center.", keys={
        "shape": Key("box", "str", "Domain shape.", choices=("box", "sphere")),
        "box": Key([20.0, 20.0, 20.0], "vec3", "Box edge lengths (full lengths).", "L"),
        "radius": Key(0.0, "float", "Radius for shape = sphere.", "L"),
        "center": Key([0.0, 0.0, 0.0], "vec3", "Domain centre.", "L"),
    }),
    "boundary": Table(doc="Boundary conditions. Each axis takes one value or [low, high].", keys={
        "x": Key("periodic", "bc", "Boundary along x.", choices=BC_CHOICES),
        "y": Key("periodic", "bc", "Boundary along y.", choices=BC_CHOICES),
        "z": Key("periodic", "bc", "Boundary along z.", choices=BC_CHOICES),
        "sphere": Key("reflecting", "str", "Wall type for domain.shape = sphere.", choices=("reflecting", "absorbing", "open")),
        "on_absorb": Key("stick", "str", "Absorbed particles stick to the wall (and are frozen) or are removed.", choices=("stick", "remove")),
        "neutralize": Key(True, "bool", "Absorbed particles lose their charge (dust-fall)."),
        "restitution": Key(1.0, "float", "Normal restitution coefficient for reflecting walls."),
    }),
    "time": Table(doc="Time stepping.", keys={
        "dt": Key(0.01, "float", "Time step.", "T"),
        "t_end": Key(10.0, "float", "End time.", "T"),
        "steps": Key(0, "int", "Number of steps (overrides t_end if > 0)."),
    }),
    "interaction": Table(doc="Particle-particle interaction.", keys={
        "model": Key("yukawa", "str", "Pair potential.", choices=("yukawa", "coulomb", "none")),
        "kappa": Key(1.0, "float", "Screening parameter kappa = a / lambda_D (reduced units)."),
        "debye_length": Key(0.0, "float", "Debye length; overrides kappa when > 0 (natural input in SI).", "L"),
        "cutoff": Key(0.0, "float", "Interaction cutoff (0 = automatic from tolerance).", "L"),
        "tolerance": Key(1e-4, "float", "Automatic cutoff: pair force at the cutoff relative to its value at r = a (1e-4 gives rc = 7.4 a at kappa = 1)."),
        "method": Key("auto", "str", "auto picks cells when the cutoff fits the box, otherwise direct O(N^2).", choices=("auto", "cells", "direct")),
    }),
    "neutrals": Table(doc="Neutral gas drag and thermal noise (Langevin, integrated with BAOAB).", keys={
        "drag": Key("none", "str", "Drag model.", choices=("none", "constant", "epstein")),
        "gamma": Key(0.0, "float", "Damping rate for drag = constant.", "RATE"),
        "temperature": Key(0.0, "float", "Neutral gas temperature for thermal kicks.", "TEMP"),
        "noise": Key(True, "bool", "Add the fluctuation-dissipation thermal force."),
        "pressure": Key(0.0, "float", "Gas pressure in Pa (epstein, SI).", "SI"),
        "mass_amu": Key(39.95, "float", "Neutral atom mass in amu (epstein).", "SI"),
        "delta": Key(1.44, "float", "Epstein accommodation coefficient (1.0 to 1.44)."),
    }),
    "thermostat": Table(doc="Optional thermostat on top of neutral drag.", keys={
        "type": Key("none", "str", "Thermostat type.", choices=("none", "berendsen", "langevin")),
        "temperature": Key(0.0, "float", "Target temperature.", "TEMP"),
        "tau": Key(0.1, "float", "Berendsen coupling time.", "T"),
        "gamma": Key(1.0, "float", "Langevin damping rate.", "RATE"),
        "until": Key(-1.0, "float", "Switch the thermostat off after this time (< 0 = never).", "T"),
    }),
    "collisions": Table(doc="Hard-sphere contact collisions between finite-radius particles.", keys={
        "enabled": Key(False, "bool", "Enable contact collisions."),
        "restitution": Key(1.0, "float", "Normal restitution coefficient."),
    }),
    "forces": FORCES,
    "field": Table(doc="External electric field on a 3D grid (PINC output or a ppdyn field file), F = Q E.", keys={
        "file": Key("", "str", "Field file.", required=True),
        "format": Key("auto", "str", "File format.", choices=("auto", "ppdyn", "pinc")),
        "step": Key(-1, "int", "PINC time-step group to read (-1 = last)."),
        "scale": Key(1.0, "float", "Extra multiplier on E."),
        "outside": Key("zero", "str", "Field outside the grid.", choices=("zero", "clamp", "periodic")),
        "center_grid": Key(True, "bool", "Centre a PINC grid on the domain centre."),
        "auto_solve": Key(True, "bool", "If the file does not exist and [[objects]] are defined, solve the field with [field_solver] and save it before the run."),
    }),
    "objects": Table(array=True, doc="Solid objects (signed-distance primitives). Repeat for several; they are combined as a union.", keys={
        "type": Key("box", "str", "Primitive type.", choices=("box", "sphere", "cylinder", "cubesat")),
        "center": Key([0.0, 0.0, 0.0], "vec3", "Centre.", "L"),
        "size": Key([1.0, 1.0, 1.0], "vec3", "Edge lengths (box).", "L"),
        "radius": Key(1.0, "float", "Radius (sphere, cylinder).", "L"),
        "length": Key(1.0, "float", "Length (cylinder).", "L"),
        "axis": Key("z", "str", "Axis (cylinder).", choices=("x", "y", "z")),
        "action": Key("reflect", "str", "What happens to a particle that hits it.", choices=("reflect", "absorb")),
        "restitution": Key(1.0, "float", "Normal restitution for action = reflect."),
        "potential": Key(0.0, "float", "Conductor potential in volts used by ppdyn field.", "SI"),
        "cubesat_size": Key("1U", "str", "CubeSat form factor.", choices=("1U", "2U", "3U", "6U")),
        "wings": Key("deployed", "str", "CubeSat solar panels.", choices=("none", "body_mounted", "deployed")),
        "wing_length": Key(0.15, "float", "CubeSat deployed wing length (metres, scaled with scale).", "SI"),
        "antenna": Key(False, "bool", "Add a CubeSat antenna mast."),
        "scale": Key(1.0, "float", "Geometric scale factor applied to CubeSat dimensions (metres to input length unit)."),
    }),
    "field_solver": Table(doc="Settings for `ppdyn field`, which solves the (screened) Poisson problem around the objects and writes a field file.", keys={
        "grid": Key([64, 64, 64], "list", "Grid points per axis."),
        "debye_length": Key(0.0, "float", "Screening length (0 = pure Laplace).", "L"),
        "output": Key("field.h5", "str", "Output field file."),
        "tolerance": Key(1e-8, "float", "Linear solver relative tolerance."),
    }),
    "output": Table(doc="Output and diagnostics.", keys={
        "dir": Key("output", "str", "Output directory."),
        "every": Key(100, "int", "Snapshot interval in steps."),
        "history_every": Key(10, "int", "Energy/temperature history interval in steps."),
        "quantities": Key(["position", "velocity"], "list", "Per-particle data written each snapshot (position, velocity, charge, state, force)."),
        "precision": Key("float32", "str", "Snapshot precision.", choices=("float32", "float64")),
        "compression": Key("lzf", "str", "HDF5 compression.", choices=("none", "lzf", "gzip")),
        "xdmf": Key(True, "bool", "Write an XDMF file so ParaView opens the HDF5 directly."),
        "vtk": Key(False, "bool", "Also write legacy .vtk files per snapshot."),
        "checkpoint_every": Key(0, "int", "Checkpoint interval in steps (0 = only at the end)."),
        "overwrite": Key(True, "bool", "Replace an existing output directory."),
    }),
}

SECTION_ORDER = ["run", "units", "domain", "boundary", "time", "species", "interaction",
                 "neutrals", "thermostat", "collisions", "forces", "field", "objects",
                 "field_solver", "output"]
SCHEMA["species"] = SPECIES
