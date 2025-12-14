# Cubesat Geometry and Electric Field Generator

This tool generates cubesat geometries and computes electric fields around them, producing HDF5 files compatible with PPDyn. It removes the dependency on PINC by providing a standalone solution.

## Features

- **Cubesat Geometry Generation**: Creates realistic cubesat geometries (1U, 2U, 3U, 6U) with solar panel wings
- **Finite Element Meshing**: Uses Gmsh for high-quality mesh generation
- **Electric Field Solver**: Solves Poisson's equation to compute electric field around charged cubesats
- **PPDyn Compatibility**: Exports to HDF5 format directly readable by PPDyn

## Installation

### Prerequisites

1. **Gmsh**: Install Gmsh library
   ```bash
   # Ubuntu/Debian
   sudo apt-get install libgmsh-dev

   # Or via conda
   conda install -c conda-forge gmsh
   ```

2. **Python Dependencies**: Install required packages
   ```bash
   pip install -r requirements.txt
   ```

   Or install manually:
   ```bash
   pip install numpy scipy gmsh h5py ini-parser
   ```

## Usage

### Basic Usage

1. **Create Configuration File** (optional - default will be created):
   ```bash
   # Edit cubesat_config.ini or use defaults
   ```

2. **Generate Cubesat and Electric Field**:
   ```bash
   python src/generate_cubesat_fields.py -c cubesat_config.ini
   ```

3. **Use with PPDyn**:
   ```bash
   # The output files will be in PIC_data/ directory
   # Update input.ini to point to this directory:
   # picDir = PIC_data

   python src/main.py
   ```

### Configuration Options

Edit `cubesat_config.ini` to customize:

#### Cubesat Parameters
- `size`: Cubesat size (`1U`, `2U`, `3U`, `6U`)
- `wing_config`: Wing configuration (`standard`, `deployed`, `none`)
- `wing_length`: Length of solar panel wings in meters
- `include_antenna`: Include antenna elements (True/False)

#### Domain Parameters
- `Lx`, `Ly`, `Lz`: Domain size in meters (should match PPDyn input.ini)
- `nx`, `ny`, `nz`: Grid resolution

#### Electric Field Parameters
- `floating_potential`: Floating potential in volts
- `solver_method`: Solver method (`sparse` or `fd`)

#### Output Parameters
- `output_dir`: Output directory for HDF5 files

### Example Configurations

#### 1U Cubesat with Standard Wings
```ini
[cubesat]
size = 1U
wing_config = standard
wing_length = 0.15
include_antenna = False

[domain]
Lx = 0.5
Ly = 0.5
Lz = 0.5
nx = 100
ny = 100
nz = 100

[electric_field]
floating_potential = 1.0
solver_method = sparse
```

#### 3U Cubesat with Deployed Wings
```ini
[cubesat]
size = 3U
wing_config = deployed
wing_length = 0.20
include_antenna = True

[domain]
Lx = 1.0
Ly = 1.0
Lz = 1.0
nx = 150
ny = 150
nz = 150

[electric_field]
floating_potential = 2.5
solver_method = sparse
```

## Output Files

The generator produces two HDF5 files compatible with PPDyn:

1. **E.grid.h5**: Electric field data
   - Attribute: `Axis denormalization factor`
   - Dataset: `timestep_XXXXXX` with shape `(nx, ny, nz, 3)` containing [Ex, Ey, Ez]

2. **object.grid.h5**: Object geometry mask
   - Attribute: `Axis denormalization factor`
   - Dataset: `Object` with shape `(nx, ny, nz)` as uint8 binary mask (1=inside, 0=outside)

## Integration with PPDyn

1. Generate cubesat geometry and field:
   ```bash
   python src/generate_cubesat_fields.py -c cubesat_config.ini
   ```

2. Update `input.ini`:
   ```ini
   [directory]
   picDir = PIC_data  ; or path to your output directory

   [options]
   PIC_data = True
   object_data = True
   ```

3. Run PPDyn:
   ```bash
   python src/main.py
   ```

## Architecture

The tool consists of several modules:

- **cubesat_generator.py**: Generates cubesat geometry using Gmsh
- **geometry_to_grid.py**: Converts unstructured mesh to regular voxel grid
- **field_solver.py**: Solves Poisson equation for electric field
- **hdf5_exporter.py**: Exports data to PPDyn-compatible HDF5 format
- **generate_cubesat_fields.py**: Main orchestration script

## Troubleshooting

### Gmsh Not Found
```bash
# Install Gmsh Python bindings
pip install gmsh

# Or install system package
sudo apt-get install python3-gmsh
```

### Memory Issues with Large Grids
- Reduce grid resolution (`nx`, `ny`, `nz`)
- Use `sparse` solver method (more memory efficient)
- Process in chunks for very large domains

### Electric Field Not Converging
- Increase domain size relative to cubesat
- Check boundary conditions
- Try different solver methods

## Future Improvements

- [ ] More accurate geometry-to-grid conversion (ray casting)
- [ ] Support for complex geometries (CAD import)
- [ ] Adaptive mesh refinement
- [ ] Time-dependent field evolution
- [ ] Multiple objects support
- [ ] GUI for geometry design

## License

Same as PPDyn (MIT License)
