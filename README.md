<p align="center"><img src="assets/ppdyn_logo.png" alt="PPDyn" width="420"></p>

# PPDyn: Plasma Particle Dynamics

[![CI](https://github.com/sayanadhikari/PPDyn/actions/workflows/ci.yml/badge.svg)](https://github.com/sayanadhikari/PPDyn/actions/workflows/ci.yml)
[![Documentation Status](https://readthedocs.org/projects/ppdyn/badge/?version=latest)](https://ppdyn.readthedocs.io/en/latest/?badge=latest)
[![PyPI Version](https://img.shields.io/pypi/v/ppdyn.svg)](https://pypi.org/project/PPDyn/)
[![DOI](https://zenodo.org/badge/349242730.svg)](https://zenodo.org/badge/latestdoi/349242730)

Parallel molecular dynamics for dusty (complex) plasmas and charged micro-debris. Grains interact through a screened Coulomb (Yukawa) potential and feel gravity, electric fields, ion drag, neutral drag and thermal noise, solid objects and absorbing walls. Numba compiles the kernels to multithreaded machine code; a linked-cell list makes the cost scale linearly with the number of particles.

**Version 2** merges every former development branch (dust fall, dust voids, dust clouds in microgravity, debris around spacecraft) into one code. You choose the physics in the input file. The previous code is kept on the [`legacy`](https://github.com/sayanadhikari/PPDyn/tree/legacy) branch.

<p align="center"><img src="docs/media/gallery.png" alt="Example results" width="100%"></p>

## Install

```bash
pip install ppdyn                 # core
pip install "ppdyn[all]"          # plus the interactive 3D viewer and the AMG field solver
```

or from a clone with `pip install -e ".[dev]"`. Python 3.9 or newer; no compiler needed.

## Quick start

```bash
ppdyn init yukawa_crystal                 # copy a ready-to-run example
cd yukawa_crystal
ppdyn run input.toml --report             # run, then write figures and an HTML report
ppdyn view output/yukawa_crystal.h5       # interactive 3D viewer with a time slider
```

A complete input file is short and readable, and unknown or misspelt keys are rejected with a suggestion:

```toml
[domain]
box = [24.36, 24.36, 24.36]          # reduced units: a, 1/omega0, Q0, m0

[[species]]
count = 3456
[species.placement]
region = "lattice"
lattice = "bcc"
[species.velocity]
temperature = 0.0025                  # T = 1/Gamma, Gamma = 400

[interaction]
kappa = 1.0

[thermostat]
type = "langevin"
temperature = 0.01                    # heat to Gamma = 100: the crystal melts
gamma = 0.1

[time]
dt = 0.05
t_end = 150.0
```

Set `units.system = "SI"` to give every quantity in SI units (metres, seconds, kelvin, coulomb, pascal, V/m) and let PPDyn do the normalisation.

## What it can simulate

| feature | input |
|---|---|
| Yukawa or Coulomb interactions, automatic cutoff, linked cells or direct sum | `[interaction]` |
| Several species; fixed, normal, lognormal, uniform, binary or bimodal grain sizes; mass from density; charge ∝ r<sup>p</sup> or from surface potential | `[[species]]` |
| Loading in a box, sphere, shell, plane, bcc/fcc/sc lattice or from a file; Maxwellian, uniform or drifting velocities | `[species.placement]`, `[species.velocity]` |
| Periodic, reflecting, absorbing (stick or remove, optional charge neutralisation) or open walls, per face; spherical chambers | `[boundary]`, `[domain]` |
| Gravity, uniform and linear (sheath) E-fields, harmonic traps, constant ion-flow force | `[forces.*]` |
| Radial void force with optional time modulation | `[forces.void]` |
| Tabulated radial force or field profiles (measured or digitised data) | `[[forces.radial_table]]` |
| 3D electric field from PINC output or PPDyn's own field solver, trilinear interpolation | `[field]` |
| Solid objects (box, sphere, cylinder, CubeSat 1U to 6U with wings and antenna) that reflect or absorb | `[[objects]]` |
| Screened-Poisson field solver around charged objects | `ppdyn field`, `[field_solver]` |
| Epstein or constant neutral drag with exact Langevin thermal noise, Langevin or Berendsen thermostats | `[neutrals]`, `[thermostat]` |
| Hard-sphere contact collisions with restitution | `[collisions]` |
| Bit-exact checkpoint and restart | `run.restart`, `--restart` |

Every key is documented in the [input reference](docs/source/input_reference.md), which is generated from the code so it cannot drift out of date.

## Examples

`ppdyn init` lists them; each reproduces the physics of one former branch. Details and results are in [docs/source/examples.md](docs/source/examples.md).

| example | former branch | what you see |
|---|---|---|
| `yukawa_crystal` | `main` | a bcc Yukawa crystal heated through its melting point |
| `dust_fall` | `dust-fall` | charged grains fall through argon onto an absorbing, neutralising electrode; large grains land first |
| `dust_void` | `dust-void` | a monolayer dust crystal with an ion-drag void and ion flow |
| `dust_cloud` | `dust-exp` | a microgravity cloud of two grain sizes forms a void and segregates into two shells under tabulated electric and ion-drag forces |
| `debris_cubesat` | `debris`, `debris_dynamics` | charged micro-debris attracted to a charged 3U CubeSat; deposits outline the body and solar wings |

## Output and post-processing

Each run writes one HDF5 file with snapshots, a per-step history (energies, temperature per axis, particle counts, collected charge) and the full input text, plus a live `history.csv`. ParaView opens the accompanying `.xdmf` file directly.

```bash
ppdyn plot run.h5 report          # every figure and an index.html
ppdyn plot run.h5 rdf             # history, snapshot, rdf, vdf, msd, profiles, absorption
ppdyn movie run.h5 --view xz      # mp4 (or gif without ffmpeg)
ppdyn view run.h5                 # standalone interactive HTML
ppdyn export run.h5 vtk           # also csv, xdmf
ppdyn watch output/               # live plot while running
```

From Python:

```python
from ppdyn.post import Run, rdf, msd
with Run("output/yukawa_crystal.h5") as r:
    r_, g = rdf(r.position(-1), r.L, r.periodic)
```

## Performance

Time per step for a Yukawa system (κ = 1, a = 1 density) on two CPU cores; the legacy timings use the old `main` pusher on the same machine.

| particles | PPDyn 1.x (O(N²)) | PPDyn 2 (linked cells) | speed-up |
|---:|---:|---:|---:|
| 2 000 | 43 ms | 19 ms | 2.3× |
| 10 000 | 1.25 s | 70 ms | 18× |
| 30 000 | 11.1 s | 0.20 s | 56× |
| 464 000 (dust-exp size) | ≈ 44 min (extrapolated) | ≈ 3.1 s (extrapolated) | ≈ 850× |

Throughput grows with the core count; run `ppdyn bench` to measure your machine.

## Validation

The test suite (`pytest -q`, 35 tests) checks, among others:

* cell-list forces equal the direct sum to 10⁻¹⁰, including thin non-periodic slabs;
* NVE energy error 3×10⁻⁶ that falls 3.9× when dt is halved (second order);
* a Langevin gas reaches its target temperature within 2 %;
* free fall onto an absorbing floor lands at √(2h/g) to within one step, and the charge is collected;
* specular reflection from spheres, boxes and cylinders, and elastic head-on collisions;
* the field solver against the analytic screened sphere φ = V R/r e<sup>−(r−R)/λ</sup>;
* SI and reduced-unit runs of the same system give identical trajectories;
* restart from a checkpoint reproduces an uninterrupted run bit for bit;
* legacy inputs from every branch convert and run.

The numerical methods, units and the legacy bugs that were fixed on the way are described in [docs/source/numerics.md](docs/source/numerics.md).

## Coming from PPDyn 1.x

```bash
ppdyn convert input.ini     # writes input.toml with NOTE comments for every judgement call
```

See [docs/source/migration.md](docs/source/migration.md) for where each branch's features went and the behaviour changes to be aware of.

## Contributors

* [Sayan Adhikari](https://github.com/sayanadhikari), IFE, Norway
* [Rupak Mukherjee](https://github.com/RupakMukherjee), PPPL, USA
* [Gaute Holen](https://github.com/GauteHolen), UiO, Norway (dust clouds, size distributions)
* [Rinku Mishra](https://github.com/rinku-mishra), IPR, India (dust voids, debris dynamics)

Contributions are welcome; see [CONTRIBUTING.md](CONTRIBUTING.md). If you use PPDyn, please cite it using [CITATION.cff](CITATION.cff).

## License

MIT, see [LICENSE](LICENSE).
