# Troubleshooting Guide

## NaN Energy Values in Simulation

If you see `Energy: nan` in the PPDyn output, this indicates invalid electric field values.

### Causes:
1. **Fast approximation producing invalid values** - Fixed in latest version
2. **Field normalization issues** - Division by zero or very small values
3. **NaN/Inf in exported field** - Not properly cleaned before export

### Solutions:

1. **Regenerate field with validation**:
   ```bash
   python src/generate_cubesat_fields.py -c cubesat_config.ini
   ```
   The latest version now validates and cleans NaN/Inf values.

2. **Use sparse solver instead of fast**:
   ```ini
   [electric_field]
   solver_method = sparse
   fast_mode = False
   ```
   The sparse solver is more accurate and stable.

3. **Check field file**:
   ```python
   import h5py
   import numpy as np

   with h5py.File('PIC_data/E.grid.h5', 'r') as f:
       keys = sorted(f.keys())
       E = f[keys[-1]][:]
       print(f"NaN count: {np.sum(np.isnan(E))}")
       print(f"Inf count: {np.sum(np.isinf(E))}")
   ```

4. **Verify normalization method**:
   In `input.ini`, try:
   ```ini
   E_norm_method = mean  ; Instead of 'norm'
   ```

## Gmsh Mesh Warnings

The mesh generation warning is usually harmless. It's typically about:
- Small elements
- Mesh quality
- Geometric tolerances

These warnings don't affect the simulation. To suppress:
- Already handled in code (verbosity reduced)

## Field Values Too Large

If electric field values are extremely large (>1e10 V/m):
- Check floating potential value
- Reduce domain size relative to cubesat
- Use proper normalization

## Performance Issues

See `PERFORMANCE_TIPS.md` for optimization strategies.

## Common Issues

### Issue: "No volume elements found in mesh"
**Solution**: Check that mesh_size is appropriate for domain size.

### Issue: "Field not converging"
**Solution**:
- Increase domain size
- Use iterative solver for large grids
- Check boundary conditions

### Issue: "Memory error"
**Solution**:
- Reduce grid resolution
- Use iterative solver
- Process in smaller chunks
