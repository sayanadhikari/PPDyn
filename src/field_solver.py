"""
Electric Field Solver for Cubesat

Solves Poisson's equation to compute electric field around a charged cubesat.
Uses FEniCS for finite element solving, with fallback to scipy for simpler cases.
"""

import numpy as np
import h5py
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve, cg, gmres
from scipy.interpolate import griddata
import sys
import time
import threading

# Try to import tqdm for progress bars
try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    # Simple progress bar replacement
    def tqdm(iterable, desc="", total=None):
        if total is None:
            total = len(iterable)
        class SimpleProgress:
            def __init__(self, iterable, desc, total):
                self.iterable = iterable
                self.desc = desc
                self.total = total
                self.current = 0
            def __iter__(self):
                print(f"{self.desc}...", end="", flush=True)
                for item in self.iterable:
                    self.current += 1
                    if self.current % max(1, self.total // 20) == 0:
                        pct = 100 * self.current / self.total
                        print(f"\r{self.desc}... {pct:.0f}%", end="", flush=True)
                    yield item
                print(f"\r{self.desc}... 100%", flush=True)
        return SimpleProgress(iterable, desc, total)


class FieldSolver:
    """
    Solve Poisson equation: -∇²φ = ρ/ε₀
    For floating potential: φ = V_float on cubesat surface, φ = 0 at infinity
    """

    def __init__(self, domain_size, grid_resolution, floating_potential=1.0):
        """
        Initialize field solver.

        Parameters:
        -----------
        domain_size : tuple
            (Lx, Ly, Lz) - domain size in meters
        grid_resolution : tuple
            (nx, ny, nz) - grid resolution
        floating_potential : float
            Floating potential in volts
        """
        self.Lx, self.Ly, self.Lz = domain_size
        self.nx, self.ny, self.nz = grid_resolution
        self.V_float = floating_potential

        # Grid spacing (actual spacing from linspace)
        # xg goes from -Lx to +Lx with nx points, so spacing = 2*Lx/(nx-1)
        if self.nx > 1:
            self.dx = 2 * self.Lx / (self.nx - 1)
        else:
            self.dx = 2 * self.Lx
        if self.ny > 1:
            self.dy = 2 * self.Ly / (self.ny - 1)
        else:
            self.dy = 2 * self.Ly
        if self.nz > 1:
            self.dz = 2 * self.Lz / (self.nz - 1)
        else:
            self.dz = 2 * self.Lz

        # Grid coordinates (centered at origin, matching PPDyn convention)
        self.xg = np.linspace(-self.Lx, self.Lx, self.nx)
        self.yg = np.linspace(-self.Ly, self.Ly, self.ny)
        self.zg = np.linspace(-self.Lz, self.Lz, self.nz)

        # Create meshgrid
        self.X, self.Y, self.Z = np.meshgrid(self.xg, self.yg, self.zg, indexing='ij')

    def solve_poisson_fd(self, object_mask):
        """
        Solve Poisson equation using finite differences.

        Parameters:
        -----------
        object_mask : ndarray
            3D boolean array, True where object exists

        Returns:
        --------
        phi : ndarray
            Electric potential (nx, ny, nz)
        Ex, Ey, Ez : ndarray
            Electric field components (nx, ny, nz)
        """
        print("Solving Poisson equation using finite differences...")

        # Initialize potential
        phi = np.zeros((self.nx, self.ny, self.nz))

        # Set boundary conditions
        # Object surface: V = V_float
        phi[object_mask] = self.V_float

        # Domain boundaries: V = 0 (ground)
        phi[0, :, :] = 0.0
        phi[-1, :, :] = 0.0
        phi[:, 0, :] = 0.0
        phi[:, -1, :] = 0.0
        phi[:, :, 0] = 0.0
        phi[:, :, -1] = 0.0

        # Iterative solver (Jacobi method)
        max_iter = 1000
        tolerance = 1e-6

        # Create mask for points to solve (not object, not boundary)
        solve_mask = ~object_mask
        solve_mask[0, :, :] = False
        solve_mask[-1, :, :] = False
        solve_mask[:, 0, :] = False
        solve_mask[:, -1, :] = False
        solve_mask[:, :, 0] = False
        solve_mask[:, :, -1] = False

        # Coefficients for Laplacian
        cx = 1.0 / (self.dx**2)
        cy = 1.0 / (self.dy**2)
        cz = 1.0 / (self.dz**2)
        denom = 2.0 * (cx + cy + cz)

        for iteration in range(max_iter):
            phi_old = phi.copy()

            # Update interior points
            for i in range(1, self.nx-1):
                for j in range(1, self.ny-1):
                    for k in range(1, self.nz-1):
                        if solve_mask[i, j, k]:
                            phi[i, j, k] = (
                                cx * (phi_old[i+1, j, k] + phi_old[i-1, j, k]) +
                                cy * (phi_old[i, j+1, k] + phi_old[i, j-1, k]) +
                                cz * (phi_old[i, j, k+1] + phi_old[i, j, k-1])
                            ) / denom

            # Check convergence
            if np.max(np.abs(phi - phi_old)) < tolerance:
                print(f"Converged after {iteration+1} iterations")
                break

        # Compute electric field: E = -∇φ
        # Use np.gradient with proper spacing
        grad_phi = np.gradient(phi, self.dx, self.dy, self.dz)
        Ex = -grad_phi[0]  # -∂φ/∂x
        Ey = -grad_phi[1]  # -∂φ/∂y
        Ez = -grad_phi[2]  # -∂φ/∂z

        return phi, Ex, Ey, Ez

    def solve_poisson_sparse(self, object_mask, use_iterative=False, max_iter=1000, tol=1e-6):
        """
        Solve Poisson equation using sparse matrix solver (more efficient).

        Parameters:
        -----------
        object_mask : ndarray
            3D boolean array, True where object exists
        use_iterative : bool
            Use iterative solver (CG) instead of direct solver (faster for large systems)
        max_iter : int
            Maximum iterations for iterative solver
        tol : float
            Tolerance for iterative solver

        Returns:
        --------
        phi : ndarray
            Electric potential (nx, ny, nz)
        Ex, Ey, Ez : ndarray
            Electric field components (nx, ny, nz)
        """
        print("Solving Poisson equation using sparse matrix solver...")
        start_time = time.time()

        # Flatten grid
        n_total = self.nx * self.ny * self.nz
        mask_flat = object_mask.flatten()
        print(f"  Grid size: {self.nx} x {self.ny} x {self.nz} = {n_total:,} points")

        # Build sparse matrix for Laplacian
        print("  Building sparse matrix...")
        row_indices = []
        col_indices = []
        data = []
        rhs = np.zeros(n_total)

        # Coefficients
        cx = 1.0 / (self.dx**2)
        cy = 1.0 / (self.dy**2)
        cz = 1.0 / (self.dz**2)
        diag = -2.0 * (cx + cy + cz)

        idx = 0
        total_points = self.nx * self.ny * self.nz
        progress_step = max(1, total_points // 50)  # Update every 2%

        for k in tqdm(range(self.nz), desc="  Building matrix", leave=False):
            for j in range(self.ny):
                for i in range(self.nx):
                    # Boundary conditions
                    if (i == 0 or i == self.nx-1 or
                        j == 0 or j == self.ny-1 or
                        k == 0 or k == self.nz-1):
                        # Dirichlet BC: phi = 0
                        row_indices.append(idx)
                        col_indices.append(idx)
                        data.append(1.0)
                        rhs[idx] = 0.0
                    elif mask_flat[idx]:
                        # Object: phi = V_float
                        row_indices.append(idx)
                        col_indices.append(idx)
                        data.append(1.0)
                        rhs[idx] = self.V_float
                    else:
                        # Interior point: Laplacian
                        row_indices.append(idx)
                        col_indices.append(idx)
                        data.append(diag)

                        # Neighbors
                        if i > 0:
                            row_indices.append(idx)
                            col_indices.append(idx - 1)
                            data.append(cx)
                        if i < self.nx - 1:
                            row_indices.append(idx)
                            col_indices.append(idx + 1)
                            data.append(cx)
                        if j > 0:
                            row_indices.append(idx)
                            col_indices.append(idx - self.nx)
                            data.append(cy)
                        if j < self.ny - 1:
                            row_indices.append(idx)
                            col_indices.append(idx + self.nx)
                            data.append(cy)
                        if k > 0:
                            row_indices.append(idx)
                            col_indices.append(idx - self.nx * self.ny)
                            data.append(cz)
                        if k < self.nz - 1:
                            row_indices.append(idx)
                            col_indices.append(idx + self.nx * self.ny)
                            data.append(cz)

                        rhs[idx] = 0.0

                    idx += 1

        # Build sparse matrix
        print("  Assembling sparse matrix...", end="", flush=True)
        A = csr_matrix((data, (row_indices, col_indices)), shape=(n_total, n_total))
        print(f" done ({A.nnz:,} non-zero elements)")

        # Solve system
        print("  Solving linear system...", end="", flush=True)
        solve_start = time.time()

        # Progress indicator for long-running solves
        progress_active = [True]
        def show_progress(solver_type=""):
            """Show elapsed time updates during solve"""
            base_msg = f"  Solving linear system... ({solver_type})" if solver_type else "  Solving linear system..."
            while progress_active[0]:
                time.sleep(2)  # Update every 2 seconds
                if progress_active[0]:
                    elapsed = time.time() - solve_start
                    if elapsed < 60:
                        print(f"\r{base_msg} [elapsed: {elapsed:.1f}s]", end="", flush=True)
                    else:
                        mins = int(elapsed // 60)
                        secs = int(elapsed % 60)
                        print(f"\r{base_msg} [elapsed: {mins}m {secs}s]", end="", flush=True)

        if use_iterative and n_total > 50000:
            # Use iterative solver for large systems (faster, less memory)
            solver_type = "iterative CG solver"
            progress_thread = threading.Thread(target=show_progress, args=(solver_type,), daemon=True)
            progress_thread.start()
            try:
                phi_flat, info = cg(A, rhs, tol=tol, maxiter=max_iter)
                progress_active[0] = False
                if info != 0:
                    print(f"\n    Warning: CG solver did not converge (info={info}), falling back to direct solver...")
                    progress_active[0] = True
                    solve_start = time.time()  # Reset timer for direct solver
                    solver_type = "direct solver (fallback)"
                    progress_thread = threading.Thread(target=show_progress, args=(solver_type,), daemon=True)
                    progress_thread.start()
                    # Fallback to direct solver if CG fails
                    phi_flat = spsolve(A, rhs)
                    progress_active[0] = False
            except KeyboardInterrupt:
                progress_active[0] = False
                raise
        else:
            # Use direct solver for smaller systems
            solver_type = "direct solver"
            progress_thread = threading.Thread(target=show_progress, args=(solver_type,), daemon=True)
            progress_thread.start()
            try:
                phi_flat = spsolve(A, rhs)
                progress_active[0] = False
            except KeyboardInterrupt:
                progress_active[0] = False
                raise

        solve_time = time.time() - solve_start
        # Clear progress line and show final time
        print(f"\r  Solving linear system... done ({solve_time:.1f}s)                    ", flush=True)

        phi = phi_flat.reshape((self.nx, self.ny, self.nz))

        # Validate solution
        phi_max = np.max(np.abs(phi))
        phi_mean = np.mean(np.abs(phi))
        if phi_max > 1e10 * self.V_float:
            print(f"    Warning: Potential values are extremely large (max={phi_max:.2e} V)")
            print(f"      This suggests the solver may have failed. Expected range: ~{self.V_float} V")

        total_time = time.time() - start_time
        print(f"  Total solve time: {total_time:.2f} seconds")

        # Compute electric field: E = -∇φ
        print("  Computing electric field from potential...", end="", flush=True)

        # Use actual grid coordinates for gradient (more accurate)
        grad_phi = np.gradient(phi, self.xg, self.yg, self.zg)
        Ex = -grad_phi[0]  # -∂φ/∂x
        Ey = -grad_phi[1]  # -∂φ/∂y
        Ez = -grad_phi[2]  # -∂φ/∂z

        # Validate and clean field values
        Ex = np.nan_to_num(Ex, nan=0.0, posinf=0.0, neginf=0.0)
        Ey = np.nan_to_num(Ey, nan=0.0, posinf=0.0, neginf=0.0)
        Ez = np.nan_to_num(Ez, nan=0.0, posinf=0.0, neginf=0.0)

        # Clip extremely large field values to prevent numerical overflow
        # For V=1V, L=0.5m, reasonable field should be ~2-20 V/m near object
        # But allow up to 1 MV/m to account for sharp gradients near boundaries
        E_mag = np.sqrt(Ex**2 + Ey**2 + Ez**2)
        max_E = 1e6  # 1 MV/m absolute maximum
        E_max_actual = np.max(E_mag)
        E_mean_actual = np.mean(E_mag)

        if E_max_actual > max_E:
            clipped_count = np.sum(E_mag > max_E)
            print(f"\n    Warning: Clipping {clipped_count} field values ({100*clipped_count/len(E_mag.flatten()):.1f}%)")
            print(f"      Max field: {E_max_actual:.2e} V/m -> {max_E:.2e} V/m")
            print(f"      Mean field: {E_mean_actual:.2e} V/m")
            # Avoid divide by zero: only scale where E_mag > max_E and E_mag > 0
            scale = np.ones_like(E_mag)
            mask = (E_mag > max_E) & (E_mag > 0)
            scale[mask] = max_E / E_mag[mask]
            Ex = Ex * scale
            Ey = Ey * scale
            Ez = Ez * scale
        else:
            print(f" (E_max={E_max_actual:.2e} V/m, E_mean={E_mean_actual:.2e} V/m)")

        print(" done")

        return phi, Ex, Ey, Ez

    def solve_fast_approximation(self, object_mask):
        """
        Fast approximation: Use analytical solution for point charge or simple geometry.
        Good for testing but less accurate.

        Parameters:
        -----------
        object_mask : ndarray
            3D boolean array, True where object exists

        Returns:
        --------
        phi : ndarray
            Electric potential
        Ex, Ey, Ez : ndarray
            Electric field components
        """
        print("Using fast approximation (analytical solution)...")

        # Find center of object
        obj_indices = np.where(object_mask)
        if len(obj_indices[0]) == 0:
            # No object, return zero field
            phi = np.zeros((self.nx, self.ny, self.nz))
            Ex = np.zeros_like(phi)
            Ey = np.zeros_like(phi)
            Ez = np.zeros_like(phi)
            return phi, Ex, Ey, Ez

        # Center of mass of object
        cx = self.xg[obj_indices[0]].mean()
        cy = self.yg[obj_indices[1]].mean()
        cz = self.zg[obj_indices[2]].mean()

        # Approximate as point charge: V = V0 * R / r (for r > R)
        # where R is characteristic size of object
        obj_size = np.sqrt(
            (self.xg[obj_indices[0]].max() - self.xg[obj_indices[0]].min())**2 +
            (self.yg[obj_indices[1]].max() - self.yg[obj_indices[1]].min())**2 +
            (self.zg[obj_indices[2]].max() - self.zg[obj_indices[2]].min())**2
        ) / 2.0

        R = max(obj_size, self.dx)  # Minimum size

        # Compute potential using vectorized operations (much faster)
        # Create coordinate grids
        X, Y, Z = np.meshgrid(self.xg, self.yg, self.zg, indexing='ij')

        # Distance from center
        r = np.sqrt((X - cx)**2 + (Y - cy)**2 + (Z - cz)**2)

        # Avoid division by zero
        r = np.maximum(r, 1e-10)

        # Potential: V = V0 * R / r for r > R, V = V0 for r <= R
        phi = np.where(object_mask, self.V_float,
                       np.where(r > R, self.V_float * R / r, self.V_float))

        # Ensure no NaN or Inf
        phi = np.nan_to_num(phi, nan=self.V_float, posinf=self.V_float, neginf=-self.V_float)

        # Compute electric field: E = -∇φ
        # Use finite differences to avoid issues with gradient
        grad_phi = np.gradient(phi, self.dx, self.dy, self.dz)
        Ex = -grad_phi[0]
        Ey = -grad_phi[1]
        Ez = -grad_phi[2]

        # Clean up any NaN or Inf values
        Ex = np.nan_to_num(Ex, nan=0.0, posinf=0.0, neginf=0.0)
        Ey = np.nan_to_num(Ey, nan=0.0, posinf=0.0, neginf=0.0)
        Ez = np.nan_to_num(Ez, nan=0.0, posinf=0.0, neginf=0.0)

        # Validate field values
        if np.any(np.isnan(Ex)) or np.any(np.isnan(Ey)) or np.any(np.isnan(Ez)):
            print("  Warning: NaN values detected in electric field, setting to zero")
            Ex = np.nan_to_num(Ex, nan=0.0)
            Ey = np.nan_to_num(Ey, nan=0.0)
            Ez = np.nan_to_num(Ez, nan=0.0)

        if np.any(np.isinf(Ex)) or np.any(np.isinf(Ey)) or np.any(np.isinf(Ez)):
            print("  Warning: Inf values detected in electric field, setting to zero")
            Ex = np.nan_to_num(Ex, posinf=0.0, neginf=0.0)
            Ey = np.nan_to_num(Ey, posinf=0.0, neginf=0.0)
            Ez = np.nan_to_num(Ez, posinf=0.0, neginf=0.0)

        return phi, Ex, Ey, Ez

    def solve(self, object_mask, method='sparse', fast_mode=False, use_iterative=None):
        """
        Solve for electric field.

        Parameters:
        -----------
        object_mask : ndarray
            3D boolean array, True where object exists
        method : str
            'sparse', 'fd' (finite difference), or 'fast' (approximation)
        fast_mode : bool
            Use fast approximation (overrides method)
        use_iterative : bool or None
            Use iterative solver for sparse method (None = auto-detect)

        Returns:
        --------
        phi : ndarray
            Electric potential
        Ex, Ey, Ez : ndarray
            Electric field components
        """
        if fast_mode or method == 'fast':
            return self.solve_fast_approximation(object_mask)
        elif method == 'sparse':
            # Auto-detect: use iterative for large grids
            n_total = self.nx * self.ny * self.nz
            if use_iterative is None:
                use_iterative = (n_total > 50000)
            return self.solve_poisson_sparse(object_mask, use_iterative=use_iterative)
        else:
            return self.solve_poisson_fd(object_mask)
