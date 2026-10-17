# PPDyn input reference

Generated from `src/ppdyn/schema.py` with `ppdyn docs`. Unknown keys are rejected with a suggestion, so typos never pass silently.

The *unit* column gives the physical dimension. With `units.system = "reduced"` values are in reduced units (length a, time 1/ω₀, charge Q₀, mass m₀, temperature Q₀²/(4πε₀a) = 1/Γ). With `units.system = "SI"` the same keys take SI values (m, s, C, kg, K, V/m, N) and are converted internally. Keys marked *SI* are always SI.

## `[run]`

Run control.

| key | default | unit | description |
|---|---|---|---|
| `name` | `"ppdyn"` |  | Run name, used in output file names. |
| `threads` | `0` |  | Numba threads (0 = all cores). |
| `seed` | `12345` |  | Random seed. |
| `progress` | `true` |  | Show a progress bar. |
| `restart` | `""` |  | Checkpoint file to restart from. |

## `[units]`

Unit system. In reduced units lengths are in a, time in 1/omega0, charge in Q0, mass in m0 with omega0^2 = Q0^2 / (4 pi eps0 m0 a^3), so temperature equals 1/Gamma. With system = "SI" every dimensional input is SI and is converted internally.

| key | default | unit | description |
|---|---|---|---|
| `system` | `"reduced"` |  | Input unit system. Choices: `reduced`, `SI`. |
| `length` | `0.0` | SI | SI only: reference length a in metres (0 = Wigner-Seitz radius of the loaded density). |
| `charge` | `0.0` | SI | SI only: reference charge Q0 in coulomb (0 = mean |Q|). |
| `mass` | `0.0` | SI | SI only: reference mass m0 in kg (0 = mean mass). |

## `[domain]`

Simulation domain, centred on center.

| key | default | unit | description |
|---|---|---|---|
| `shape` | `"box"` |  | Domain shape. Choices: `box`, `sphere`. |
| `box` | `[20.0, 20.0, 20.0]` | length | Box edge lengths (full lengths). |
| `radius` | `0.0` | length | Radius for shape = sphere. |
| `center` | `[0.0, 0.0, 0.0]` | length | Domain centre. |

## `[boundary]`

Boundary conditions. Each axis takes one value or [low, high].

| key | default | unit | description |
|---|---|---|---|
| `x` | `"periodic"` |  | Boundary along x. Choices: `periodic`, `reflecting`, `absorbing`, `open`. |
| `y` | `"periodic"` |  | Boundary along y. Choices: `periodic`, `reflecting`, `absorbing`, `open`. |
| `z` | `"periodic"` |  | Boundary along z. Choices: `periodic`, `reflecting`, `absorbing`, `open`. |
| `sphere` | `"reflecting"` |  | Wall type for domain.shape = sphere. Choices: `reflecting`, `absorbing`, `open`. |
| `on_absorb` | `"stick"` |  | Absorbed particles stick to the wall (and are frozen) or are removed. Choices: `stick`, `remove`. |
| `neutralize` | `true` |  | Absorbed particles lose their charge (dust-fall). |
| `restitution` | `1.0` |  | Normal restitution coefficient for reflecting walls. |

## `[time]`

Time stepping.

| key | default | unit | description |
|---|---|---|---|
| `dt` | `0.01` | time | Time step. |
| `t_end` | `10.0` | time | End time. |
| `steps` | `0` |  | Number of steps (overrides t_end if > 0). |

## `[[species]]`

Particle species. Repeat [[species]] for several populations.

| key | default | unit | description |
|---|---|---|---|
| `name` | `"dust"` |  | Species label. |
| `count` | `1000` |  | **Required.** Number of particles. |
| `mass` | `1.0` | mass | Mass at the reference radius (ignored when density is set). |
| `density` | `0.0` | density | Material density; if > 0, m = density * 4/3 pi r^3. |
| `reference_radius` | `0.0` | length | Radius at which mass/charge values apply (0 = mean radius). |

### `[species.radius]`

Particle size distribution (radius, or mass when variable = "mass").

| key | default | unit | description |
|---|---|---|---|
| `distribution` | `"fixed"` |  | Distribution type. Choices: `fixed`, `normal`, `lognormal`, `uniform`, `bimodal`, `binary`. |
| `variable` | `"radius"` |  | Quantity the distribution is drawn for. "mass" reproduces the legacy dist = True option (m ~ N(mean, std), r ~ m^(1/3)). Choices: `radius`, `mass`. |
| `mean` | `0.0` | length | Mean (fixed value for "fixed"). For variable = "mass" it is in units of species.mass. |
| `std` | `0.0` | length | Standard deviation (normal, lognormal: of ln r). |
| `low` | `0.0` | length | Lower bound (uniform). |
| `high` | `0.0` | length | Upper bound (uniform). |
| `values` | `[]` | length | Radii for binary/bimodal (two entries). |
| `stds` | `[]` | length | Standard deviations for bimodal (two entries). |
| `weights` | `[0.5, 0.5]` |  | Population fractions for binary/bimodal. |
| `min` | `0.0` | length | Samples below this are redrawn (keeps radii positive). |

### `[species.charge]`

Charge model. Q = value * (r / reference_radius)^exponent, or Q = 4 pi eps0 r phi.

| key | default | unit | description |
|---|---|---|---|
| `model` | `"power"` |  | "power" (fixed when exponent = 0) or "potential" (SI only, OML-like Q = 4 pi eps0 r phi). Choices: `power`, `potential`. |
| `value` | `1.0` | charge | Charge at the reference radius. |
| `exponent` | `1.0` |  | Size exponent p (1 = constant surface potential, 0 = equal charges, 2 = legacy dist). |
| `potential` | `-2.0` | SI | Surface potential in volts for model = "potential". |

### `[species.placement]`

Initial placement.

| key | default | unit | description |
|---|---|---|---|
| `region` | `"box"` |  | Where particles are loaded. Choices: `box`, `sphere`, `shell`, `plane`, `lattice`, `file`. |
| `center` | `auto` | length | Region centre (defaults to domain centre). |
| `size` | `auto` | length | Box size for region = box/lattice (defaults to the domain). |
| `r_min` | `0.0` | length | Inner radius for sphere/shell. |
| `r_max` | `0.0` | length | Outer radius for sphere/shell. |
| `axis` | `"z"` |  | Normal axis for region = plane. Choices: `x`, `y`, `z`. |
| `position` | `0.0` | length | Plane coordinate for region = plane (e.g. top of the box for dust fall). |
| `lattice` | `"bcc"` |  | Lattice type for region = lattice. Choices: `sc`, `bcc`, `fcc`. |
| `file` | `""` |  | HDF5/CSV/NPY file with an (N,3) position array for region = file. |
| `jitter` | `0.0` | length | Random displacement added to lattice sites. |
| `avoid_objects` | `true` |  | Redraw particles that land inside objects. |

### `[species.velocity]`

Initial velocities.

| key | default | unit | description |
|---|---|---|---|
| `distribution` | `"maxwell"` |  | Velocity loading. Choices: `maxwell`, `uniform`, `zero`. |
| `temperature` | `0.0` | temperature | Temperature for maxwell loading. |
| `vmax` | `[1.0, 1.0, 1.0]` | velocity | Full width of the uniform distribution per axis. |
| `drift` | `[0.0, 0.0, 0.0]` | velocity | Mean drift velocity added after loading. |
| `remove_net_momentum` | `true` |  | Subtract the centre-of-mass velocity before adding drift. |

## `[interaction]`

Particle-particle interaction.

| key | default | unit | description |
|---|---|---|---|
| `model` | `"yukawa"` |  | Pair potential. Choices: `yukawa`, `coulomb`, `none`. |
| `kappa` | `1.0` |  | Screening parameter kappa = a / lambda_D (reduced units). |
| `debye_length` | `0.0` | length | Debye length; overrides kappa when > 0 (natural input in SI). |
| `cutoff` | `0.0` | length | Interaction cutoff (0 = automatic from tolerance). |
| `tolerance` | `0.0001` |  | Automatic cutoff: pair force at the cutoff relative to its value at r = a (1e-4 gives rc = 7.4 a at kappa = 1). |
| `method` | `"auto"` |  | auto picks cells when the cutoff fits the box, otherwise direct O(N^2). Choices: `auto`, `cells`, `direct`. |

## `[neutrals]`

Neutral gas drag and thermal noise (Langevin, integrated with BAOAB).

| key | default | unit | description |
|---|---|---|---|
| `drag` | `"none"` |  | Drag model. Choices: `none`, `constant`, `epstein`. |
| `gamma` | `0.0` | 1/time | Damping rate for drag = constant. |
| `temperature` | `0.0` | temperature | Neutral gas temperature for thermal kicks. |
| `noise` | `true` |  | Add the fluctuation-dissipation thermal force. |
| `pressure` | `0.0` | SI | Gas pressure in Pa (epstein, SI). |
| `mass_amu` | `39.95` | SI | Neutral atom mass in amu (epstein). |
| `delta` | `1.44` |  | Epstein accommodation coefficient (1.0 to 1.44). |

## `[thermostat]`

Optional thermostat on top of neutral drag.

| key | default | unit | description |
|---|---|---|---|
| `type` | `"none"` |  | Thermostat type. Choices: `none`, `berendsen`, `langevin`. |
| `temperature` | `0.0` | temperature | Target temperature. |
| `tau` | `0.1` | time | Berendsen coupling time. |
| `gamma` | `1.0` | 1/time | Langevin damping rate. |
| `until` | `-1.0` | time | Switch the thermostat off after this time (< 0 = never). |

## `[collisions]`

Hard-sphere contact collisions between finite-radius particles.

| key | default | unit | description |
|---|---|---|---|
| `enabled` | `false` |  | Enable contact collisions. |
| `restitution` | `1.0` |  | Normal restitution coefficient. |

## `[forces]`

External forces. A sub-table switches a force on.

### `[forces.gravity]`

Uniform gravity, F = m g.

| key | default | unit | description |
|---|---|---|---|
| `g` | `[0.0, 0.0, -9.81]` | acceleration | Gravitational acceleration. |

### `[forces.efield]`

Uniform electric field, F = Q E.

| key | default | unit | description |
|---|---|---|---|
| `E` | `[0.0, 0.0, 0.0]` | E-field | Field vector. |

### `[forces.linear_efield]`

Linear (sheath-like) field, E = E0 + G (x - origin) per axis, F = Q E.

| key | default | unit | description |
|---|---|---|---|
| `E0` | `[0.0, 0.0, 0.0]` | E-field | Field at the origin. |
| `gradient` | `[0.0, 0.0, 0.0]` | E-field/length | dE_i/dx_i for each axis. |
| `origin` | `[0.0, 0.0, 0.0]` | length | Reference point. |

### `[forces.harmonic]`

Harmonic confinement acceleration a = -omega2 (x - origin), independent of mass. Legacy reflecting runs used this along z.

| key | default | unit | description |
|---|---|---|---|
| `omega2` | `[0.0, 0.0, 0.0]` | 1/time² | Squared trap frequency per axis. |
| `origin` | `[0.0, 0.0, 0.0]` | length | Trap centre. |

### `[forces.constant]`

Constant force, e.g. ion flow drag (dust-void f_flow).

| key | default | unit | description |
|---|---|---|---|
| `F` | `[0.0, 0.0, 0.0]` | force | Force on a reference particle. |
| `scale_with` | `"none"` |  | Per-particle scaling of F. Choices: `none`, `charge`, `mass`, `radius`, `radius2`. |

### `[forces.void]`

Radial ion-drag void force F = A exp(-r/lambda) [1 + m sin(w t)] r_hat (dust-void branch).

| key | default | unit | description |
|---|---|---|---|
| `amplitude` | `0.0` | force | Force amplitude A on a reference particle. |
| `decay_length` | `1.0` | length | Decay length lambda. |
| `center` | `[0.0, 0.0, 0.0]` | length | Void centre. |
| `scale_with` | `"radius2"` |  | Per-particle scaling (ion drag ~ r^2). Choices: `none`, `charge`, `mass`, `radius`, `radius2`. |
| `modulation` | `0.0` |  | Relative modulation depth m. |
| `omega` | `0.0` | 1/time | Modulation angular frequency w. |
| `cylindrical` | `false` |  | Use distance from the z axis instead of the centre point. |

### `[[forces.radial_table]]`

Tabulated radial force F(r) r_hat (dust-exp E-field and ion-drag profiles). Repeat for several tables.

| key | default | unit | description |
|---|---|---|---|
| `file` | `""` |  | **Required.** CSV (two columns r, F) or NPZ file. |
| `r_key` | `"r"` |  | Array name for r inside an NPZ file (e.g. arr_0). |
| `f_key` | `"F"` |  | Array name for F inside an NPZ file (e.g. arr_1). |
| `r_unit` | `1.0` |  | Multiply r by this to get the input unit system. |
| `factor` | `1.0` |  | Multiply F by this (e.g. legacy k_fel). |
| `kind` | `"force"` |  | Table holds a force on a reference particle or an electric field. Choices: `force`, `efield`. |
| `center` | `[0.0, 0.0, 0.0]` | length | Centre of the radial profile. |
| `scale_with` | `"none"` |  | Per-particle scaling (ignored for kind = efield, which scales with charge). Choices: `none`, `charge`, `mass`, `radius`, `radius2`. |

## `[field]`

External electric field on a 3D grid (PINC output or a ppdyn field file), F = Q E.

| key | default | unit | description |
|---|---|---|---|
| `file` | `""` |  | **Required.** Field file. |
| `format` | `"auto"` |  | File format. Choices: `auto`, `ppdyn`, `pinc`. |
| `step` | `-1` |  | PINC time-step group to read (-1 = last). |
| `scale` | `1.0` |  | Extra multiplier on E. |
| `outside` | `"zero"` |  | Field outside the grid. Choices: `zero`, `clamp`, `periodic`. |
| `center_grid` | `true` |  | Centre a PINC grid on the domain centre. |
| `auto_solve` | `true` |  | If the file does not exist and [[objects]] are defined, solve the field with [field_solver] and save it before the run. |

## `[[objects]]`

Solid objects (signed-distance primitives). Repeat for several; they are combined as a union.

| key | default | unit | description |
|---|---|---|---|
| `type` | `"box"` |  | Primitive type. Choices: `box`, `sphere`, `cylinder`, `cubesat`. |
| `center` | `[0.0, 0.0, 0.0]` | length | Centre. |
| `size` | `[1.0, 1.0, 1.0]` | length | Edge lengths (box). |
| `radius` | `1.0` | length | Radius (sphere, cylinder). |
| `length` | `1.0` | length | Length (cylinder). |
| `axis` | `"z"` |  | Axis (cylinder). Choices: `x`, `y`, `z`. |
| `action` | `"reflect"` |  | What happens to a particle that hits it. Choices: `reflect`, `absorb`. |
| `restitution` | `1.0` |  | Normal restitution for action = reflect. |
| `potential` | `0.0` | SI | Conductor potential in volts used by ppdyn field. |
| `cubesat_size` | `"1U"` |  | CubeSat form factor. Choices: `1U`, `2U`, `3U`, `6U`. |
| `wings` | `"deployed"` |  | CubeSat solar panels. Choices: `none`, `body_mounted`, `deployed`. |
| `wing_length` | `0.15` | SI | CubeSat deployed wing length (metres, scaled with scale). |
| `antenna` | `false` |  | Add a CubeSat antenna mast. |
| `scale` | `1.0` |  | Geometric scale factor applied to CubeSat dimensions (metres to input length unit). |

## `[field_solver]`

Settings for `ppdyn field`, which solves the (screened) Poisson problem around the objects and writes a field file.

| key | default | unit | description |
|---|---|---|---|
| `grid` | `[64, 64, 64]` |  | Grid points per axis. |
| `debye_length` | `0.0` | length | Screening length (0 = pure Laplace). |
| `output` | `"field.h5"` |  | Output field file. |
| `tolerance` | `1e-08` |  | Linear solver relative tolerance. |

## `[output]`

Output and diagnostics.

| key | default | unit | description |
|---|---|---|---|
| `dir` | `"output"` |  | Output directory. |
| `every` | `100` |  | Snapshot interval in steps. |
| `history_every` | `10` |  | Energy/temperature history interval in steps. |
| `quantities` | `['position', 'velocity']` |  | Per-particle data written each snapshot (position, velocity, charge, state, force). |
| `precision` | `"float32"` |  | Snapshot precision. Choices: `float32`, `float64`. |
| `compression` | `"lzf"` |  | HDF5 compression. Choices: `none`, `lzf`, `gzip`. |
| `xdmf` | `true` |  | Write an XDMF file so ParaView opens the HDF5 directly. |
| `vtk` | `false` |  | Also write legacy .vtk files per snapshot. |
| `checkpoint_every` | `0` |  | Checkpoint interval in steps (0 = only at the end). |
| `overwrite` | `true` |  | Replace an existing output directory. |
