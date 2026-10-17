# Migrating from PPDyn 1.x and the feature branches

PPDyn 2 merges the physics of every development branch into one code on `main`. The input file selects the features; there is no longer a branch per application. The previous `main` is preserved unchanged on the `legacy` branch, and the feature branches are kept for reference.

## Convert an old input file

```bash
ppdyn convert input.ini            # writes input.toml next to it
ppdyn check input.toml             # validate and show the run summary
```

The converter understands inputs from `main`, `dust-fall`, `dust-void`, `dust-exp` and `debris_dynamics`. Every judgement call it makes is written as a `# NOTE` comment in the generated file. Read those before running.

## Where each branch went

| branch | what it added | PPDyn 2 input |
|---|---|---|
| `main`, `dev`, `argon-mod` | Yukawa MD, periodic or reflecting box, Gaussian mass distribution with Q = m^(2/3) | `[interaction]`, `[boundary]`, `[species.radius] variable = "mass"` with `[species.charge] exponent = 2` |
| `dust-fall` | grains released from the top, absorbing and neutralising floor, fall times, collected charge | `[boundary] z = ["absorbing", "reflecting"]`, `neutralize = true`, `[species.placement] region = "plane"`; fall times in `final/absorb_time`, charge in the history |
| `dust-void` | radial void force, ion flow force, neutral drag with thermal kicks, time-modulated void | `[forces.void]` (with `modulation`, `omega`, `cylindrical`), `[forces.constant]`, `[neutrals]` |
| `dust-exp` | spherical geometry, digitised radial E-field and ion-drag profiles, binary and bimodal grain sizes, size-dependent forces, drag ∝ a² | `[domain] shape = "sphere"`, `[[forces.radial_table]]` with `scale_with`, `[species.radius] distribution = "binary"/"bimodal"/"uniform"`, `[neutrals] drag = "epstein"` |
| `debris_dynamics` | PINC E-field import with trilinear interpolation, object mask, SI normalisation, gravity | `[field] format = "pinc"`, `[[objects]]`, `units.system = "SI"`, `[forces.gravity]` |
| `debris` | Gmsh CubeSat generator, Laplace solver for the field around it | `[[objects]] type = "cubesat"`, `ppdyn field` (screened Poisson, AMG-CG) |

Ready-to-run versions of these cases are in `ppdyn init`: `yukawa_crystal`, `dust_fall`, `dust_void`, `dust_cloud` and `debris_cubesat`.

## Behaviour changes to be aware of

* **Domain.** Legacy runs used [−L, L] per axis. `domain.box` is the full edge length, centred on `domain.center` (default the origin). The converter doubles the lengths.
* **Booleans.** Legacy `bool(...)` made every flag true, including `vtkData = False`. TOML booleans are now real booleans, so old runs that relied on accidental VTK output will now honour `false`.
* **Initial temperature.** Legacy Maxwellian loading used v ~ N(0, Temp), i.e. T = Temp² for unit mass. The converter sets `temperature = Temp²`.
* **Legacy gravity.** The old reflecting/mixed runs applied −(z + L<sub>z</sub>) g, which is a spring anchored at the floor, not gravity. The converter keeps it as `[forces.harmonic]` so results stay comparable; use `[forces.gravity]` for real gravity.
* **dust-void random force.** The old "random" kick contained no random numbers. It is now a proper Langevin force; expect thermal motion that was previously missing.
* **Random numbers.** Old runs drew the particle masses from an unseeded NumPy generator, so two runs never had the same particles. Runs are now reproducible from `run.seed`, including restarts.
* **Output.** One HDF5 file per run (`<name>.h5`) instead of `particle.hdf5` plus `energy.txt`. Snapshots are in `snapshots/position` with shape (frames, N, 3), the history in `history/data` and `history.csv`. ParaView opens the `.xdmf` file directly; legacy-style `.vtk` files are available with `output.vtk = true` or `ppdyn export run.h5 vtk`.
* **Python API.** `from ppdyn import load, run` replaces importing `PPDyn.ppdyn`; `ppdyn.post.Run` reads output files.
