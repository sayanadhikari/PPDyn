"""
HDF5 Exporter for PPDyn-compatible format

Exports electric field and object geometry to HDF5 files
compatible with PPDyn's expected format.
"""

import h5py
import numpy as np
from os.path import join as pjoin


class HDF5Exporter:
    """
    Export data to HDF5 format compatible with PPDyn.
    """

    def __init__(self, output_dir='PIC_data'):
        """
        Initialize exporter.

        Parameters:
        -----------
        output_dir : str
            Output directory for HDF5 files
        """
        self.output_dir = output_dir
        import os
        os.makedirs(output_dir, exist_ok=True)

    def export_electric_field(self, Ex, Ey, Ez, xg, yg, zg, denorm, timestep=0):
        """
        Export electric field to E.grid.h5

        Parameters:
        -----------
        Ex, Ey, Ez : ndarray
            Electric field components (nx, ny, nz)
        xg, yg, zg : ndarray
            Grid coordinates
        denorm : float
            Axis denormalization factor
        timestep : int
            Timestep number (for dataset key)
        """
        filename = pjoin(self.output_dir, 'E.grid.h5')

        # Validate field values before exporting
        if np.any(np.isnan(Ex)) or np.any(np.isnan(Ey)) or np.any(np.isnan(Ez)):
            print("  Warning: NaN values in electric field, cleaning...")
            Ex = np.nan_to_num(Ex, nan=0.0)
            Ey = np.nan_to_num(Ey, nan=0.0)
            Ez = np.nan_to_num(Ez, nan=0.0)

        if np.any(np.isinf(Ex)) or np.any(np.isinf(Ey)) or np.any(np.isinf(Ez)):
            print("  Warning: Inf values in electric field, cleaning...")
            Ex = np.nan_to_num(Ex, posinf=0.0, neginf=0.0)
            Ey = np.nan_to_num(Ey, posinf=0.0, neginf=0.0)
            Ez = np.nan_to_num(Ez, posinf=0.0, neginf=0.0)

        # Check for reasonable values
        E_max = np.max(np.sqrt(Ex**2 + Ey**2 + Ez**2))
        if E_max > 1e10:
            print(f"  Warning: Very large electric field values detected (max={E_max:.2e} V/m)")
            print("  This might cause numerical issues in PPDyn")

        with h5py.File(filename, 'w') as f:
            # Set attribute
            f.attrs['Axis denormalization factor'] = np.array([denorm])

            # Combine field components: shape (nx, ny, nz, 3)
            E_field = np.zeros((Ex.shape[0], Ex.shape[1], Ex.shape[2], 3), dtype=np.float64)
            E_field[:, :, :, 0] = Ex.astype(np.float64)
            E_field[:, :, :, 1] = Ey.astype(np.float64)
            E_field[:, :, :, 2] = Ez.astype(np.float64)

            # Final validation
            if np.any(np.isnan(E_field)) or np.any(np.isinf(E_field)):
                print("  Error: NaN/Inf values still present after cleaning!")
                E_field = np.nan_to_num(E_field, nan=0.0, posinf=0.0, neginf=0.0)

            # Create dataset with timestep key
            dataset_key = f'timestep_{timestep:06d}'
            f.create_dataset(dataset_key, data=E_field, compression='gzip')

            # Store field statistics as attributes
            f.attrs['E_max'] = float(np.max(np.sqrt(Ex**2 + Ey**2 + Ez**2)))
            f.attrs['E_mean'] = float(np.mean(np.sqrt(Ex**2 + Ey**2 + Ez**2)))

        print(f"Exported electric field to {filename}")
        print(f"  Field statistics: E_max = {E_max:.4e} V/m, E_mean = {np.mean(np.sqrt(Ex**2 + Ey**2 + Ez**2)):.4e} V/m")

    def export_object_grid(self, object_mask, denorm):
        """
        Export object geometry to object.grid.h5

        Parameters:
        -----------
        object_mask : ndarray
            3D boolean or uint8 array (1=inside object, 0=outside)
        denorm : float
            Axis denormalization factor
        """
        filename = pjoin(self.output_dir, 'object.grid.h5')

        # Ensure uint8 format
        if object_mask.dtype != np.uint8:
            obj_mask_uint8 = object_mask.astype(np.uint8)
        else:
            obj_mask_uint8 = object_mask

        with h5py.File(filename, 'w') as f:
            # Set attribute
            f.attrs['Axis denormalization factor'] = np.array([denorm])

            # Create dataset
            f.create_dataset('Object', data=obj_mask_uint8, compression='gzip')

        print(f"Exported object grid to {filename}")

    def export_both(self, Ex, Ey, Ez, object_mask, xg, yg, zg, denorm, timestep=0):
        """
        Export both electric field and object grid.

        Parameters:
        -----------
        Ex, Ey, Ez : ndarray
            Electric field components
        object_mask : ndarray
            Object mask
        xg, yg, zg : ndarray
            Grid coordinates
        denorm : float
            Axis denormalization factor
        timestep : int
            Timestep number
        """
        self.export_electric_field(Ex, Ey, Ez, xg, yg, zg, denorm, timestep)
        self.export_object_grid(object_mask, denorm)
