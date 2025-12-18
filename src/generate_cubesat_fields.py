#!/usr/bin/env python
"""
Main script to generate cubesat geometry and electric field for PPDyn

This script:
1. Generates cubesat geometry using Gmsh
2. Converts geometry to regular grid
3. Solves Poisson equation for electric field
4. Exports to HDF5 format compatible with PPDyn
"""

import numpy as np
import argparse
import ini
import os
from os.path import join as pjoin

from cubstat_generator1 import RealisticCubesatGenerator as CubesatGenerator
# from cubstat_generator import CubesatGenerator
from geometry_to_grid import GeometryToGrid
from field_solver import FieldSolver
from hdf5_exporter import HDF5Exporter


def parse_config(config_file='cubesat_config.ini'):
    """Parse configuration file."""
    if not os.path.exists(config_file):
        # Create default config
        create_default_config(config_file)

    params = ini.parse(open(config_file).read())
    return params


def create_default_config(config_file):
    """Create default configuration file."""
    default_config = """[cubesat]
size = 1U              ; Cubesat size: 1U, 2U, 3U, 6U
wing_config = standard ; Wing configuration: standard, deployed, none
wing_length = 0.15    ; Wing length in meters
include_antenna = False ; Include antenna elements

[domain]
Lx = 0.5              ; Domain size in X (meters)
Ly = 0.5              ; Domain size in Y (meters)
Lz = 0.5              ; Domain size in Z (meters)
nx = 100              ; Grid resolution in X
ny = 100              ; Grid resolution in Y
nz = 100              ; Grid resolution in Z

[electric_field]
floating_potential = 1.0  ; Floating potential in volts
solver_method = sparse    ; Solver method: sparse, fd

[output]
output_dir = PIC_data     ; Output directory for HDF5 files
"""
    with open(config_file, 'w') as f:
        f.write(default_config)
    print(f"Created default configuration file: {config_file}")


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description='Generate cubesat geometry and electric field for PPDyn'
    )
    parser.add_argument(
        '-c', '--config',
        default='cubesat_config.ini',
        help='Configuration file (default: cubesat_config.ini)'
    )
    parser.add_argument(
        '-o', '--output',
        default=None,
        help='Output directory (overrides config)'
    )
    args = parser.parse_args()

    # Parse configuration
    params = parse_config(args.config)

    # Cubesat parameters
    size = str(params['cubesat']['size'])
    wing_config = str(params['cubesat']['wing_config'])
    wing_length = float(params['cubesat']['wing_length'])
    include_antenna = bool(params['cubesat'].get('include_antenna', False))

    # Domain parameters
    Lx = float(params['domain']['Lx'])
    Ly = float(params['domain']['Ly'])
    Lz = float(params['domain']['Lz'])
    nx = int(params['domain']['nx'])
    ny = int(params['domain']['ny'])
    nz = int(params['domain']['nz'])

    # Electric field parameters
    V_float = float(params['electric_field']['floating_potential'])
    solver_method = str(params['electric_field'].get('solver_method', 'sparse'))
    fast_mode = bool(params['electric_field'].get('fast_mode', False))
    use_iterative_str = str(params['electric_field'].get('use_iterative', 'auto')).lower()
    if use_iterative_str == 'auto':
        use_iterative = None  # Auto-detect
    elif use_iterative_str in ('true', '1', 'yes'):
        use_iterative = True
    else:
        use_iterative = False

    # Output parameters
    output_dir = args.output or str(params['output']['output_dir'])

    print("=" * 60)
    print("Cubesat Geometry and Electric Field Generator")
    print("=" * 60)
    print(f"Cubesat size: {size}")
    print(f"Wing config: {wing_config}")
    print(f"Domain: ({Lx}, {Ly}, {Lz}) m")
    print(f"Grid: ({nx}, {ny}, {nz})")
    print(f"Floating potential: {V_float} V")
    print(f"Output directory: {output_dir}")
    print("=" * 60)

    # Step 1: Generate cubesat geometry
    print("\n[1/4] Generating cubesat geometry...")
    mesh_size = min(Lx, Ly, Lz) / max(nx, ny, nz)  # Adaptive mesh size

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Generate geometry and keep Gmsh open
    import gmsh
    
    # Map wing_config: "standard" -> "body_mounted", keep "none" and "deployed" as-is
    wing_config_mapped = wing_config
    if wing_config == "standard":
        wing_config_mapped = "body_mounted"
    
    generator = CubesatGenerator(
        size=size,
        wing_config=wing_config_mapped,
        wing_length=wing_length,
        include_antenna=include_antenna
    )
    volumes = generator.generate(mesh_size=mesh_size)

    # Calculate bounding box from generator dimensions
    wx, wy, wz = generator.wx, generator.wy, generator.wz
    wing_extent = wing_length if wing_config_mapped != "none" else 0
    bbox = (-wx/2 - wing_extent, -wy/2 - wing_extent, -wz/2,
            wx/2 + wing_extent, wy/2 + wing_extent, wz/2)
    print(f"Bounding box: {bbox}")

    # Save mesh file (optional, for visualization)
    mesh_file = pjoin(output_dir, 'cubesat.msh')
    generator.save(mesh_file)
    print(f"Mesh saved to {mesh_file}")

    # Step 2: Convert geometry to regular grid
    print("\n[2/4] Converting geometry to regular grid...")

    # Create grid coordinates (matching PPDyn convention: centered at origin)
    xg = np.linspace(-Lx, Lx, nx)
    yg = np.linspace(-Ly, Ly, ny)
    zg = np.linspace(-Lz, Lz, nz)

    converter = GeometryToGrid(xg, yg, zg)

    # Use simple method for now (can be improved with ray casting)
    try:
        object_mask = converter.mesh_to_voxels_simple()
    except Exception as e:
        print(f"Warning: {e}")
        print("Using simple bounding box method")
        object_mask = converter.mesh_to_voxels_simple()

    # Finalize Gmsh after conversion
    generator.finalize()

    # Step 3: Solve for electric field
    print("\n[3/4] Solving for electric field...")
    solver = FieldSolver(
        domain_size=(Lx, Ly, Lz),
        grid_resolution=(nx, ny, nz),
        floating_potential=V_float
    )

    phi, Ex, Ey, Ez = solver.solve(
        object_mask,
        method=solver_method,
        fast_mode=fast_mode,
        use_iterative=use_iterative
    )

    # Validate results
    E_mag = np.sqrt(Ex**2 + Ey**2 + Ez**2)
    E_max = np.max(E_mag)
    E_mean = np.mean(E_mag)

    # Check for issues
    nan_count = np.sum(np.isnan(Ex) | np.isnan(Ey) | np.isnan(Ez))
    inf_count = np.sum(np.isinf(Ex) | np.isinf(Ey) | np.isinf(Ez))

    if nan_count > 0:
        print(f"  Warning: {nan_count} NaN values in electric field!")
    if inf_count > 0:
        print(f"  Warning: {inf_count} Inf values in electric field!")

    print(f"Electric field statistics:")
    print(f"  E_max = {E_max:.4e} V/m")
    print(f"  E_mean = {E_mean:.4e} V/m")
    print(f"  E_min = {np.min(E_mag):.4e} V/m")

    # Step 4: Export to HDF5
    print("\n[4/4] Exporting to HDF5...")

    # Calculate denormalization factor (grid spacing)
    denorm = xg[1] - xg[0] if len(xg) > 1 else 1.0

    exporter = HDF5Exporter(output_dir)
    exporter.export_both(
        Ex, Ey, Ez,
        object_mask,
        xg, yg, zg,
        denorm,
        timestep=0
    )

    print("\n" + "=" * 60)
    print("Generation complete!")
    print(f"Output files:")
    print(f"  - {pjoin(output_dir, 'E.grid.h5')}")
    print(f"  - {pjoin(output_dir, 'object.grid.h5')}")
    print("=" * 60)

    return 0


if __name__ == '__main__':
    exit(main())
