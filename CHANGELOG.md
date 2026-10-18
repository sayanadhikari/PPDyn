# Changelog

## 2.0.0 (2026-10-18)

Complete rewrite. One code on `main` replaces the feature branches; the former `main` is preserved on `legacy`.

### Added
* New logo (`docs/logo/make_logo.py`, SVG with outlined text plus PNG, light, dark and icon variants).
* 3D gallery: PyVista renders and animations of all five examples, reproducible with `docs/media/render/render_all.py` and the `render` extra.
* TOML input with strict validation, typo suggestions and a generated input reference; SI or reduced units.
* Linked-cell Yukawa forces with cell-ordered data (O(N), 18x faster at 10^4 and 56x at 3x10^4 particles on two cores).
* BAOAB Langevin integrator: neutral drag (constant or Epstein) with exact thermal noise; Langevin and Berendsen thermostats.
* Multiple species; normal, lognormal, uniform, binary and bimodal size distributions; charge from size power law or surface potential.
* Boundaries per face (periodic, reflecting, absorbing with stick/remove and charge neutralisation, open); spherical domains.
* External forces: gravity, uniform and linear E-fields, harmonic traps, constant force, void force with time modulation, tabulated radial profiles, 3D grid fields (PINC or PPDyn format).
* Signed-distance objects (box, sphere, cylinder, CubeSat) that reflect or absorb; screened-Poisson field solver (`ppdyn field`).
* Contact collisions with restitution (race-free).
* Single HDF5 output with XDMF for ParaView, live CSV history, VTK/CSV export, bit-exact checkpoints and restart.
* Post-processing: g(r) (3D and monolayer), MSD and diffusion, velocity distributions, density profiles, absorption statistics, HTML report, movies, interactive 3D viewer.
* `ppdyn` command line tool, `ppdyn convert` for legacy input files, five examples, 35 tests, CI on Linux, macOS and Windows.

### Fixed (relative to 1.x and the branches)
* Boolean input flags were always true (`bool("False")`).
* Legacy "gravity" was a spring force; dust-void random kick was not random; void modulation used dt instead of t.
* Physical constants off by a factor of 10 in dust-exp; data race in the collision loop; unseeded mass sampling.
