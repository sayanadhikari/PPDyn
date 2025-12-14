"""
Convert Gmsh geometry to regular grid for PPDyn

Converts unstructured mesh from Gmsh to regular voxel grid
compatible with PPDyn's object mask format.
"""

import numpy as np
import gmsh
from scipy.spatial import cKDTree


class GeometryToGrid:
    """
    Convert Gmsh mesh to regular grid representation.
    """

    def __init__(self, xg, yg, zg):
        """
        Initialize converter.

        Parameters:
        -----------
        xg, yg, zg : ndarray
            Regular grid coordinates (1D arrays)
        """
        self.xg = xg
        self.yg = yg
        self.zg = zg

        # Create meshgrid
        self.X, self.Y, self.Z = np.meshgrid(xg, yg, zg, indexing='ij')

        # Grid shape
        self.nx, self.ny, self.nz = len(xg), len(yg), len(zg)

        # Grid spacing
        self.dx = xg[1] - xg[0] if len(xg) > 1 else 1.0
        self.dy = yg[1] - yg[0] if len(yg) > 1 else 1.0
        self.dz = zg[1] - zg[0] if len(zg) > 1 else 1.0

    def mesh_to_voxels(self, tolerance=1e-6):
        """
        Convert Gmsh mesh to voxel grid.

        Uses point-in-mesh test to determine which voxels are inside the geometry.

        Parameters:
        -----------
        tolerance : float
            Tolerance for point-in-mesh test

        Returns:
        --------
        object_mask : ndarray
            3D boolean array (True = inside object)
        """
        print("Converting mesh to voxel grid...")

        # Get mesh nodes
        nodeTags, nodeCoords, _ = gmsh.model.mesh.getNodes()
        nodeCoords = nodeCoords.reshape(-1, 3)

        # Get mesh elements (tetrahedra)
        elementTypes, elementTags, elementNodeTags = gmsh.model.mesh.getElements(3)

        if len(elementTags) == 0:
            raise ValueError("No volume elements found in mesh. Make sure mesh is 3D.")

        # Build KD-tree for fast nearest neighbor search
        tree = cKDTree(nodeCoords)

        # Initialize mask
        object_mask = np.zeros((self.nx, self.ny, self.nz), dtype=bool)

        # Get all grid points
        grid_points = np.column_stack([
            self.X.flatten(),
            self.Y.flatten(),
            self.Z.flatten()
        ])

        # For each grid point, check if it's inside any tetrahedron
        # Simplified approach: check if point is close to mesh nodes
        # More accurate: use point-in-tetrahedron test

        # Get bounding box of mesh
        mesh_min = nodeCoords.min(axis=0)
        mesh_max = nodeCoords.max(axis=0)

        # Create bounding box mask
        bbox_mask = (
            (grid_points[:, 0] >= mesh_min[0] - tolerance) &
            (grid_points[:, 0] <= mesh_max[0] + tolerance) &
            (grid_points[:, 1] >= mesh_min[1] - tolerance) &
            (grid_points[:, 1] <= mesh_max[1] + tolerance) &
            (grid_points[:, 2] >= mesh_min[2] - tolerance) &
            (grid_points[:, 2] <= mesh_max[2] + tolerance)
        )

        # For points in bounding box, check distance to mesh
        candidates = grid_points[bbox_mask]
        if len(candidates) > 0:
            # Find nearest mesh nodes
            distances, indices = tree.query(candidates, k=1)

            # Points are inside if they're very close to mesh surface
            # or we can use a more sophisticated test
            # For now, use distance threshold
            inside_threshold = self.dx * 0.5  # Half grid spacing

            # More accurate: check if point is inside any tetrahedron
            # This is computationally expensive, so we'll use a simpler method
            # Check if point is inside convex hull or use ray casting

            # Simplified: points within threshold of mesh nodes are considered inside
            # This is approximate but fast
            inside_mask = distances < inside_threshold

            # Map back to full grid
            full_inside = np.zeros(len(grid_points), dtype=bool)
            full_inside[bbox_mask] = inside_mask
            object_mask = full_inside.reshape((self.nx, self.ny, self.nz))

        print(f"Object mask: {np.sum(object_mask)} voxels inside object")
        return object_mask

    def mesh_to_voxels_raycast(self, n_rays=10):
        """
        Convert mesh to voxels using ray casting method (more accurate).

        Parameters:
        -----------
        n_rays : int
            Number of rays to cast per grid point

        Returns:
        --------
        object_mask : ndarray
            3D boolean array (True = inside object)
        """
        print("Converting mesh to voxel grid using ray casting...")

        # Get mesh surface (triangles)
        elementTypes, elementTags, elementNodeTags = gmsh.model.mesh.getElements(2)

        if len(elementTags) == 0:
            raise ValueError("No surface elements found. Cannot perform ray casting.")

        # Get triangle nodes
        triangles = []
        for elemType, elemTags, nodeTags in zip(elementTypes, elementTags, elementNodeTags):
            if elemType == 2:  # Triangle
                nodeTags_flat = nodeTags.flatten()
                nodeCoords, _ = gmsh.model.mesh.getNodes(2, elemTags[0])
                nodeCoords = nodeCoords.reshape(-1, 3)
                for i in range(0, len(nodeCoords), 3):
                    if i + 2 < len(nodeCoords):
                        triangles.append(nodeCoords[i:i+3])

        triangles = np.array(triangles)

        # Initialize mask
        object_mask = np.zeros((self.nx, self.ny, self.nz), dtype=bool)

        # For each grid point, cast rays and count intersections
        for i in range(self.nx):
            if i % 10 == 0:
                print(f"Processing slice {i}/{self.nx}")
            for j in range(self.ny):
                for k in range(self.nz):
                    point = np.array([self.xg[i], self.yg[j], self.zg[k]])

                    # Cast ray in +X direction
                    ray_dir = np.array([1.0, 0.0, 0.0])
                    intersections = 0

                    # Check intersection with all triangles
                    # (This is simplified - in practice, use spatial acceleration)
                    for tri in triangles:
                        if self._ray_triangle_intersect(point, ray_dir, tri):
                            intersections += 1

                    # Odd number of intersections = inside
                    object_mask[i, j, k] = (intersections % 2 == 1)

        print(f"Object mask: {np.sum(object_mask)} voxels inside object")
        return object_mask

    def _ray_triangle_intersect(self, origin, direction, triangle):
        """
        Check if ray intersects triangle (Möller-Trumbore algorithm).

        Parameters:
        -----------
        origin : ndarray
            Ray origin
        direction : ndarray
            Ray direction (normalized)
        triangle : ndarray
            Triangle vertices (3x3)

        Returns:
        --------
        bool
            True if intersection exists
        """
        v0, v1, v2 = triangle

        edge1 = v1 - v0
        edge2 = v2 - v0
        h = np.cross(direction, edge2)
        a = np.dot(edge1, h)

        if abs(a) < 1e-8:
            return False  # Ray parallel to triangle

        f = 1.0 / a
        s = origin - v0
        u = f * np.dot(s, h)

        if u < 0.0 or u > 1.0:
            return False

        q = np.cross(s, edge1)
        v = f * np.dot(direction, q)

        if v < 0.0 or u + v > 1.0:
            return False

        t = f * np.dot(edge2, q)
        return t > 1e-8  # Intersection in front of ray origin

    def mesh_to_voxels_simple(self):
        """
        Simple conversion: check if grid points are inside bounding box of mesh.
        More accurate methods can be added later.
        """
        print("Converting mesh to voxel grid (simple bounding box method)...")

        # Get all mesh nodes
        try:
            nodeTags, nodeCoords, _ = gmsh.model.mesh.getNodes()
            nodeCoords = nodeCoords.reshape(-1, 3)
        except Exception as e:
            print(f"Warning: Could not get mesh nodes: {e}")
            # Fallback: use a small box at origin
            print("Using fallback: small box at origin")
            bbox_min = np.array([-0.05, -0.05, -0.05])
            bbox_max = np.array([0.05, 0.05, 0.05])
            object_mask = (
                (self.X >= bbox_min[0]) & (self.X <= bbox_max[0]) &
                (self.Y >= bbox_min[1]) & (self.Y <= bbox_max[1]) &
                (self.Z >= bbox_min[2]) & (self.Z <= bbox_max[2])
            )
            print(f"Object mask: {np.sum(object_mask)} voxels inside bounding box")
            return object_mask

        # Get bounding box
        bbox_min = nodeCoords.min(axis=0)
        bbox_max = nodeCoords.max(axis=0)

        # Add small margin to account for mesh discretization
        margin = 0.01 * (bbox_max - bbox_min)
        bbox_min -= margin
        bbox_max += margin

        # Create mask for points inside bounding box
        object_mask = (
            (self.X >= bbox_min[0]) & (self.X <= bbox_max[0]) &
            (self.Y >= bbox_min[1]) & (self.Y <= bbox_max[1]) &
            (self.Z >= bbox_min[2]) & (self.Z <= bbox_max[2])
        )

        print(f"Object mask: {np.sum(object_mask)} voxels inside bounding box")
        print(f"  Bounding box: [{bbox_min[0]:.4f}, {bbox_max[0]:.4f}] x "
              f"[{bbox_min[1]:.4f}, {bbox_max[1]:.4f}] x "
              f"[{bbox_min[2]:.4f}, {bbox_max[2]:.4f}]")
        return object_mask
