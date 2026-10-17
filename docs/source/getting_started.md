# Getting started

## Install

```bash
pip install ppdyn            # or: pip install "ppdyn[all]" for the 3D viewer and the AMG field solver
```

From a clone:

```bash
git clone https://github.com/sayanadhikari/PPDyn.git
cd PPDyn
pip install -e ".[dev]"
pytest -q
```

Only Python 3.9 or newer is needed; Numba compiles the kernels on first use and caches them.

## First run

```bash
ppdyn init                       # list the examples
ppdyn init yukawa_crystal        # copy one into ./yukawa_crystal
cd yukawa_crystal
ppdyn run input.toml --report    # simulate, then write figures and output/report/index.html
ppdyn view output/yukawa_crystal.h5   # interactive 3D viewer in the browser
```

## Everyday commands

| command | what it does |
|---|---|
| `ppdyn run input.toml` | run; `-t 16` threads, `--set time.dt=0.005` overrides, `--restart output/checkpoint.h5` continues |
| `ppdyn check input.toml` | validate the input and print the run summary without running |
| `ppdyn field input.toml` | solve the electrostatic field around `[[objects]]` |
| `ppdyn plot run.h5 [all\|report\|history\|snapshot\|rdf\|vdf\|msd\|profiles\|absorption]` | figures |
| `ppdyn movie run.h5 --view xz -c charge` | mp4 (ffmpeg) or gif |
| `ppdyn view run.h5` | self-contained interactive HTML viewer with a time slider |
| `ppdyn export run.h5 vtk\|csv\|xdmf` | ParaView/VisIt/CSV export |
| `ppdyn info run.h5` | summary of an output file |
| `ppdyn watch output/` | live energy and temperature plot while a run is going |
| `ppdyn convert input.ini` | convert a PPDyn 1.x input |
| `ppdyn bench` | performance check on this machine |
| `ppdyn docs` | print the input reference |

## Python API

```python
import ppdyn
from ppdyn.post import Run, rdf

cfg = ppdyn.load("input.toml")
cfg["interaction"]["kappa"] = 2.0          # edit before running
sim, path = ppdyn.run(cfg)

with Run(path) as r:
    r_, g = rdf(r.position(-1), r.L, r.periodic)
    T = r.history["temperature"]
```

## Parameter sweeps and clusters

`--set` overrides any input key, which makes sweeps a shell loop:

```bash
for k in 0.5 1 2 3; do ppdyn run input.toml --set interaction.kappa=$k -o out_kappa$k -q; done
```

On a cluster node set the thread count with `-t` or `NUMBA_NUM_THREADS`. Example Slurm job:

```bash
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32 --time=12:00:00
module load Python/3.11
source ~/ppdyn-env/bin/activate
ppdyn run input.toml -t $SLURM_CPUS_PER_TASK -q --report
```
