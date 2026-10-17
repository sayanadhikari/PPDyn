# Numerics and physics model

## Equations of motion

Each grain *i* has mass m<sub>i</sub>, charge Q<sub>i</sub>, radius r<sub>i</sub> and obeys

m<sub>i</sub> dv<sub>i</sub>/dt = Σ<sub>j</sub> F<sub>ij</sub> + F<sub>ext</sub>(x<sub>i</sub>, t) − m<sub>i</sub> γ<sub>i</sub> v<sub>i</sub> + √(2 m<sub>i</sub> γ<sub>i</sub> k<sub>B</sub>T<sub>n</sub>) ξ<sub>i</sub>(t)

with the screened Coulomb (Yukawa) pair force

F<sub>ij</sub> = Q<sub>i</sub>Q<sub>j</sub> / (4πε₀) · (1 + κ r) e<sup>−κr</sup> / r² · r̂<sub>ij</sub>,  κ = a/λ<sub>D</sub>.

The last two terms are neutral drag and its fluctuation-dissipation partner (Langevin dynamics). Without drag the system is Hamiltonian (NVE).

## Reduced units

Internally everything is dimensionless:

| quantity | unit |
|---|---|
| length | a (Wigner-Seitz radius, (3/(4πn))<sup>1/3</sup>, unless set) |
| charge | Q₀ (mean \|Q\| unless set) |
| mass | m₀ (mean mass unless set) |
| time | 1/ω₀, ω₀² = Q₀² / (4πε₀ m₀ a³) |
| energy, temperature | Q₀²/(4πε₀ a), so T<sub>reduced</sub> = 1/Γ |
| electric field | Q₀/(4πε₀ a²) |

The dust plasma frequency is ω<sub>pd</sub> = √3 ω₀. With `units.system = "SI"` every dimensional input is SI and is converted on load; output stores reduced values plus the scale factors (`units` attribute), and `ppdyn plot --si` labels figures in SI.

## Time integration: BAOAB

One step of length h:

1. **B** v += (h/2) F/m
2. **A** x += (h/2) v
3. **O** v ← c₁ v + c₂ ξ, with c₁ = e<sup>−γh</sup>, c₂ = √((1 − c₁²) T/m), ξ ~ N(0, 1)
4. **A** x += (h/2) v, then walls and objects
5. force evaluation at the new positions
6. **B** v += (h/2) F/m, then contact collisions and an optional Berendsen rescale

The O step is the exact solution of the Ornstein-Uhlenbeck process, so the sampled temperature is correct for any γh, and with γ = 0 the scheme reduces to velocity Verlet (symplectic, second order). The test suite checks both: NVE energy error falls 3.9× when h is halved, and a Langevin gas reaches its target temperature within 2 %.

Neutral drag and the Langevin thermostat are two heat baths; they combine into one O step with γ = γ<sub>n</sub> + γ<sub>t</sub> and T = (γ<sub>n</sub>T<sub>n</sub> + γ<sub>t</sub>T<sub>t</sub>)/γ. Epstein drag is

γ = δ (4π/3) r² n<sub>n</sub> m<sub>n</sub> v̄<sub>th</sub> / m,  v̄<sub>th</sub> = √(8k<sub>B</sub>T<sub>n</sub>/(π m<sub>n</sub>)),

so each grain size gets its own damping rate.

## Pair forces

* **Cutoff.** With `interaction.cutoff = 0` the cutoff r<sub>c</sub> is where the pair force has fallen to `tolerance` (default 10⁻⁴) of its value at r = a: 7.4a at κ = 1, 4.3a at κ = 2.3. The potential energy is shifted so U(r<sub>c</sub>) = 0, which keeps the energy diagnostic smooth. On periodic axes r<sub>c</sub> is limited to half the box (minimum image).
* **Linked cells.** Cells of size r<sub>c</sub>/2 and a ±2-cell stencil halve the searched volume compared with r<sub>c</sub>-sized cells. Particle data are gathered in cell order before the loop so neighbour reads are contiguous in memory. Each particle sums its own force, so the parallel loop has no write conflicts and is deterministic.
* **Direct sum.** Used for unscreened Coulomb (κ = 0) or when requested. As in the legacy code there is no Ewald summation, so κ = 0 with periodic boundaries is a minimum-image approximation; prefer a small κ.

## External forces

All terms are evaluated in one parallel kernel. Each can scale per grain with charge, mass, r/r<sub>ref</sub> or (r/r<sub>ref</sub>)², which is how ion drag (∝ r²) and electric force (∝ Q ∝ r) separate grain sizes.

| input | force |
|---|---|
| `[forces.gravity]` | m g |
| `[forces.efield]` | Q E |
| `[forces.linear_efield]` | Q (E₀ + G·(x − x₀)) per axis (sheath) |
| `[forces.harmonic]` | −m ω² (x − x₀) per axis |
| `[forces.constant]` | constant F (ion flow) |
| `[forces.void]` | A e<sup>−r/λ</sup> (1 + m sin ωt) r̂ (spherical or cylindrical) |
| `[[forces.radial_table]]` | tabulated F(r) r̂ or Q E(r) r̂, linear interpolation |
| `[field]` | Q E(x) from a 3D grid, trilinear interpolation |

## Boundaries and objects

Faces can be periodic, reflecting (with restitution), absorbing or open. Absorbed grains either stick (frozen, kept in the output) or are removed; with `neutralize = true` they also lose their charge, which is the dust-fall model. The time and charge of every absorption are recorded. A spherical domain uses the same options.

Solid objects are signed-distance functions (box, sphere, cylinder, CubeSat composite). A grain found inside an object after a drift is mirrored across the surface along the SDF gradient and its normal velocity is reversed, or it is absorbed on the surface. SDFs are exact, so thin parts such as solar panels need no voxel grid.

## Contact collisions

For grains with radius > 0, overlapping and approaching pairs exchange the normal impulse of a hard-sphere collision with restitution e. Updates are computed from the pre-collision velocities of all pairs (Jacobi), which is race-free and independent of thread scheduling. A stuck grain acts as an immovable wall.

## Field solver (`ppdyn field`)

Solves (∇² − λ<sub>D</sub>⁻²) φ = 0 on a node grid with φ fixed to each object's potential inside it and φ = 0 on the outer boundary, using a 7-point stencil and conjugate gradients (smoothed-aggregation AMG preconditioner if `pyamg` is installed). E = −∇φ. The test suite compares against the analytic screened sphere φ = V R/r e<sup>−(r−R)/λ</sup>. Reading PINC `E.grid.h5` files is also supported (last or chosen time step, `Axis/Quantity denormalization factor` attributes).

## Reproducibility and restart

Thermal noise comes from a NumPy PCG64 generator seeded with `run.seed`. Checkpoints store all particle arrays and the generator state, and a restarted run reproduces an uninterrupted run bit for bit (tested).

## Legacy bugs fixed

| branch | issue | fix |
|---|---|---|
| all | `bool("False")` is `True`, so every boolean flag was always on | TOML booleans, strict validation |
| main, dust-fall | "gravity" was −(z + L<sub>z</sub>) g, a spring | uniform `[forces.gravity]`; the spring is still available as `[forces.harmonic]` |
| dust-void | random kick was a constant (no random numbers) | exact Langevin O step |
| dust-void | time dependence used sin(ω dt) instead of sin(ω t) | `[forces.void]` modulation uses t |
| dust-exp | ε₀ and k<sub>B</sub> written as 8.85·10e−12 (10× too large) | constants in `units.py` |
| debris | collision loop updated other particles' velocities inside a parallel loop (data race) | Jacobi collision update |
| all | O(N²) forces, 1.2 s per step at N = 10⁴ on 2 cores | linked cells, O(N): 70 ms per step |
