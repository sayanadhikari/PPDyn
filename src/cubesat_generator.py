"""
Cubesat Geometry Generator using Gmsh

Generates realistic cubesat geometries (1U, 2U, 3U, 6U) with solar panels/wings.
"""

import gmsh
import numpy as np
import os


class CubesatGenerator:
    """
    Generate cubesat geometries using Gmsh.

    Standard cubesat units (U):
    - 1U: 10cm x 10cm x 10cm
    - 2U: 10cm x 10cm x 20cm
    - 3U: 10cm x 10cm x 30cm
    - 6U: 10cm x 20cm x 30cm
    """

    # Standard cubesat dimensions in meters
    UNIT_SIZE = 0.1  # 10 cm

    def __init__(self, size='1U', wing_config='standard', wing_length=0.15):
        """
        Initialize cubesat generator.

        Parameters:
        -----------
        size : str
            Cubesat size: '1U', '2U', '3U', '6U'
        wing_config : str
            Wing configuration: 'standard', 'deployed', 'none'
        wing_length : float
            Length of solar panel wings in meters
        """
        self.size = size
        self.wing_config = wing_config
        self.wing_length = wing_length

        # Parse size
        if size == '1U':
            self.dimensions = (self.UNIT_SIZE, self.UNIT_SIZE, self.UNIT_SIZE)
        elif size == '2U':
            self.dimensions = (self.UNIT_SIZE, self.UNIT_SIZE, 2 * self.UNIT_SIZE)
        elif size == '3U':
            self.dimensions = (self.UNIT_SIZE, self.UNIT_SIZE, 3 * self.UNIT_SIZE)
        elif size == '6U':
            self.dimensions = (self.UNIT_SIZE, 2 * self.UNIT_SIZE, 3 * self.UNIT_SIZE)
        else:
            raise ValueError(f"Unknown cubesat size: {size}")

        self.wx, self.wy, self.wz = self.dimensions

        # Initialize Gmsh (suppress warnings for cleaner output)
        gmsh.initialize()
        # Optionally reduce verbosity
        gmsh.option.setNumber("General.Verbosity", 1)  # 0=errors only, 1=warnings, 2=info, 4=debug
        gmsh.model.add("cubesat")

    def create_body(self):
        """Create the main cubesat body."""
        # Main body box
        body = gmsh.model.occ.addBox(
            -self.wx/2, -self.wy/2, -self.wz/2,
            self.wx, self.wy, self.wz
        )
        return body

    def create_wings(self):
        """Create solar panel wings."""
        wings = []

        if self.wing_config == 'none':
            return wings

        # Wing thickness
        wing_thickness = 0.002  # 2mm

        # Standard configuration: 4 wings (one on each side)
        if self.wing_config in ['standard', 'deployed']:
            # Wing dimensions
            ww = self.wing_length  # wing width
            wh = max(self.wy, self.wz)  # wing height (matches body)

            # Wing positions (centered on each face)
            wing_positions = [
                # +X face
                (self.wx/2, 0, 0, 0, 0, 0),
                # -X face
                (-self.wx/2, 0, 0, 0, 0, np.pi),
                # +Y face
                (0, self.wy/2, 0, 0, 0, np.pi/2),
                # -Y face
                (0, -self.wy/2, 0, 0, 0, -np.pi/2),
            ]

            for x, y, z, rx, ry, rz in wing_positions:
                # Create wing box
                wing = gmsh.model.occ.addBox(
                    x - wing_thickness/2,
                    y - ww/2,
                    z - wh/2,
                    wing_thickness, ww, wh
                )

                # Rotate if needed (for deployed configuration)
                if self.wing_config == 'deployed':
                    # Deploy wings at 45 degrees
                    gmsh.model.occ.rotate(
                        [(3, wing)], x, y, z,
                        0, 0, 1, np.pi/4 if x > 0 else -np.pi/4
                    )

                wings.append(wing)

        return wings

    def create_antenna(self, length=0.05):
        """Create antenna elements (optional)."""
        antennas = []
        antenna_radius = 0.001  # 1mm radius

        # Add antenna on top (+Z)
        antenna = gmsh.model.occ.addCylinder(
            self.wx/4, self.wy/4, self.wz/2,
            0, 0, length,
            antenna_radius
        )
        antennas.append(antenna)

        return antennas

    def generate(self, include_antenna=False, mesh_size=0.01):
        """
        Generate the complete cubesat geometry.

        Parameters:
        -----------
        include_antenna : bool
            Whether to include antenna elements
        mesh_size : float
            Mesh element size in meters

        Returns:
        --------
        volumes : list
            List of volume tags
        """
        # Create body
        body = self.create_body()
        volumes = [body]

        # Create wings
        wings = self.create_wings()
        volumes.extend(wings)

        # Create antennas if requested
        if include_antenna:
            antennas = self.create_antenna()
            volumes.extend(antennas)

        # Fuse all volumes together
        if len(volumes) > 1:
            gmsh.model.occ.fuse([(3, volumes[0])], [(3, v) for v in volumes[1:]])
            gmsh.model.occ.synchronize()
        else:
            gmsh.model.occ.synchronize()

        # Create physical groups for boundary conditions
        # Get all surfaces
        surfaces = gmsh.model.getBoundary([(3, volumes[0])], oriented=False)

        # Create physical group for cubesat surface (for Dirichlet BC)
        gmsh.model.addPhysicalGroup(2, [s[1] for s in surfaces], 1)

        # Generate mesh
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)
        gmsh.model.mesh.generate(3)

        return volumes

    def get_bounding_box(self):
        """Get bounding box of the cubesat."""
        # Account for wings
        wing_extent = self.wing_length if self.wing_config != 'none' else 0

        xmin = -self.wx/2 - wing_extent
        xmax = self.wx/2 + wing_extent
        ymin = -self.wy/2 - wing_extent
        ymax = self.wy/2 + wing_extent
        zmin = -self.wz/2
        zmax = self.wz/2

        return (xmin, ymin, zmin, xmax, ymax, zmax)

    def save_geometry(self, filename):
        """Save geometry to file."""
        gmsh.write(filename)

    def finalize(self):
        """Finalize Gmsh."""
        gmsh.finalize()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.finalize()
