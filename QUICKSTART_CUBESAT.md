# Quick Start: Cubesat Generator for PPDyn

## Overview

This tool generates cubesat geometries and electric fields, removing the PINC dependency from PPDyn.

## Quick Start (3 Steps)

### Step 1: Install Dependencies

```bash
# Install Gmsh (choose one method)
conda install -c conda-forge gmsh
# OR
pip install gmsh

# Install other dependencies
pip install -r requirements.txt
```

### Step 2: Generate Cubesat and Electric Field

```bash
# Use default configuration
python src/generate_cubesat_fields.py

# Or specify custom config
python src/generate_cubesat_fields.py -c cubesat_config.ini -o PIC_data
```

This creates:
- `PIC_data/E.grid.h5` - Electric field data
- `PIC_data/object.grid.h5` - Object geometry mask

### Step 3: Run PPDyn

Update `input.ini`:
```ini
[directory]
picDir = PIC_data

[options]
PIC_data = True
object_data = True
```

Then run PPDyn:
```bash
python src/main.py
```

## Example: Generate 3U Cubesat

1. Edit `cubesat_config.ini`:
```ini
[cubesat]
size = 3U
wing_config = deployed
wing_length = 0.20

[domain]
Lx = 1.0
Ly = 1.0
Lz = 1.0
nx = 150
ny = 150
nz = 150

[electric_field]
floating_potential = 2.5
```

2. Generate:
```bash
python src/generate_cubesat_fields.py -c cubesat_config.ini
```

3. Use with PPDyn (update `input.ini` domain to match):
```ini
[simbox]
Lx = 1.0
Ly = 1.0
Lz = 1.0
```

## Troubleshooting

**Gmsh not found:**
```bash
pip install gmsh
# OR install system package
sudo apt-get install python3-gmsh
```

**Memory errors:**
- Reduce grid resolution (`nx`, `ny`, `nz`)
- Use smaller domain size

**Field not converging:**
- Increase domain size relative to cubesat
- Check that domain is much larger than cubesat

## Files Generated

- `PIC_data/E.grid.h5` - Electric field (compatible with PPDyn)
- `PIC_data/object.grid.h5` - Object mask (compatible with PPDyn)
- `PIC_data/cubesat.msh` - Gmsh mesh file (for visualization)

## Next Steps

See `CUBESAT_GENERATOR_README.md` for detailed documentation.
