# Performance Tips for Cubesat Field Generator

## Fast Testing Mode

For quick testing, use the **fast approximation** mode:

```ini
[electric_field]
solver_method = fast    ; or set fast_mode = True
fast_mode = True
```

This uses an analytical approximation that is **much faster** (seconds instead of minutes) but less accurate. Perfect for:
- Testing geometry generation
- Quick validation
- Development/debugging

## Optimize Grid Resolution

Reduce grid resolution for faster solving:

```ini
[domain]
nx = 50   ; Instead of 100
ny = 50
nz = 50
```

**Trade-off**: Lower resolution = faster but less accurate field.

## Use Iterative Solver

For large grids (>50,000 points), the iterative solver is automatically used:

```ini
[electric_field]
use_iterative = auto   ; auto, True, or False
```

- **auto**: Automatically chooses best method
- **True**: Force iterative (good for large grids, uses less memory)
- **False**: Force direct solver (faster for small grids)

## Progress Monitoring

The solver now shows:
- Progress during matrix construction
- Time estimates
- Grid size information
- Convergence status

Install `tqdm` for better progress bars:
```bash
pip install tqdm
```

## Performance Comparison

| Method | Speed | Accuracy | Memory | Use Case |
|--------|-------|----------|--------|----------|
| `fast` | ⚡⚡⚡ Very Fast | ⚠️ Approximate | Low | Testing |
| `sparse` (direct) | ⚡⚡ Fast | ✅ Accurate | Medium | Small grids (<50k) |
| `sparse` (iterative) | ⚡ Moderate | ✅ Accurate | Low | Large grids (>50k) |
| `fd` | ⚡ Slow | ✅ Accurate | Low | Very small grids |

## Recommended Settings

### For Testing
```ini
[domain]
nx = 50
ny = 50
nz = 50

[electric_field]
solver_method = fast
fast_mode = True
```

### For Production
```ini
[domain]
nx = 100
ny = 100
nz = 100

[electric_field]
solver_method = sparse
fast_mode = False
use_iterative = auto
```

### For High Accuracy
```ini
[domain]
nx = 150
ny = 150
nz = 150

[electric_field]
solver_method = sparse
use_iterative = True
```

## Troubleshooting Slow Performance

1. **Check grid size**: `nx * ny * nz` should be reasonable
   - < 50,000: Fast
   - 50,000 - 500,000: Moderate (use iterative)
   - > 500,000: Slow (reduce resolution)

2. **Use fast mode for testing**: Set `fast_mode = True`

3. **Reduce domain size**: Smaller domain = fewer points to solve

4. **Check memory**: Large grids need more RAM

5. **Use iterative solver**: For grids > 50k points, iterative is faster

## Example: Quick Test Run

```bash
# Edit cubesat_config.ini
[domain]
nx = 50
ny = 50
nz = 50

[electric_field]
solver_method = fast

# Run
python src/generate_cubesat_fields.py
```

This should complete in seconds instead of minutes!
