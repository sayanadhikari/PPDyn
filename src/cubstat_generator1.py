"""
Realistic CubeSat Geometry Generator using Gmsh (OCC)

Implements:
- External envelopes per Cal Poly CubeSat drawings:
  1U: 100 x 100 x 113.5 mm
  2U: 100 x 100 x 227.0 mm
  3U: 100 x 100 x 340.5 mm
  6U: 100 x 226.3 x 366.0 mm
- Rail constraint concept:
  - Minimum rail width: 8.5 mm
  - End contact pad: 6.5 mm minimum on +/-Z (modeled as "no recess" region near ends)
- Optional "tuna can" extra volume:
  - Length: 36 mm max
  - Diameter: 64 mm max

Sources: Cal Poly CubeSat Design Specification Rev 14.1 + drawings.
"""

import gmsh
import math


class RealisticCubesatGenerator:
    # Dimensions in meters derived from CDS drawings
    _SIZES = {
        "1U": (0.1000, 0.1000, 0.1135),
        "2U": (0.1000, 0.1000, 0.2270),
        "3U": (0.1000, 0.1000, 0.3405),
        "6U": (0.1000, 0.2263, 0.3660),
    }

    def __init__(
        self,
        size="1U",
        wing_config="none",        # "none", "body_mounted", "deployed"
        wing_length=0.15,
        wing_thickness=0.002,
        include_rails=True,
        rail_width=0.0085,         # 8.5 mm min
        end_contact=0.0065,        # 6.5 mm min at ends
        recess_depth=0.003,        # [Unverified] tunable recess depth
        include_tuna_can=False,
        tuna_length=0.036,         # 36 mm max
        tuna_diameter=0.064,       # 64 mm max
        include_antenna=False,
    ):
        if size not in self._SIZES:
            raise ValueError(f"Unknown cubesat size: {size}")

        self.size = size
        self.wing_config = wing_config
        self.wing_length = wing_length
        self.wing_thickness = wing_thickness

        self.include_rails = include_rails
        self.rail_width = rail_width
        self.end_contact = end_contact
        self.recess_depth = recess_depth

        self.include_tuna_can = include_tuna_can
        self.tuna_length = tuna_length
        self.tuna_diameter = tuna_diameter

        self.include_antenna = include_antenna

        self.wx, self.wy, self.wz = self._SIZES[size]

        gmsh.initialize()
        gmsh.option.setNumber("General.Verbosity", 1)
        gmsh.model.add(f"cubesat_{size}")

    # ---------- Core geometry ----------

    def _add_envelope_box(self):
        # Centered at origin
        return gmsh.model.occ.addBox(-self.wx / 2, -self.wy / 2, -self.wz / 2, self.wx, self.wy, self.wz)

    def _carve_rail_recesses(self, body_tag):
        """
        Creates the rail-based form by subtracting 4 side pockets:
        - Pockets do NOT extend into corner rail regions (rail_width)
        - Pockets do NOT extend into end contact pad regions (end_contact)
        """
        if not self.include_rails:
            return [(3, body_tag)]

        rw = self.rail_width
        ec = self.end_contact
        d = self.recess_depth

        # Guardrails to avoid inverted boxes
        if 2 * rw >= min(self.wx, self.wy):
            raise ValueError("rail_width too large for cross-section.")
        if 2 * ec >= self.wz:
            raise ValueError("end_contact too large for length.")
        if d <= 0:
            return [(3, body_tag)]

        pockets = []

        # Pocket spans in Y and Z for +/-X faces
        y_span = self.wy - 2 * rw
        z_span = self.wz - 2 * ec

        # +X pocket: subtract a slab near +X face excluding rails at corners
        pockets.append(
            gmsh.model.occ.addBox(
                self.wx / 2 - d,                 # x0
                -y_span / 2,                     # y0
                -z_span / 2,                     # z0
                d, y_span, z_span
            )
        )

        # -X pocket
        pockets.append(
            gmsh.model.occ.addBox(
                -self.wx / 2,                    # x0
                -y_span / 2,
                -z_span / 2,
                d, y_span, z_span
            )
        )

        # Pocket spans in X and Z for +/-Y faces
        x_span = self.wx - 2 * rw

        # +Y pocket
        pockets.append(
            gmsh.model.occ.addBox(
                -x_span / 2,
                self.wy / 2 - d,
                -z_span / 2,
                x_span, d, z_span
            )
        )

        # -Y pocket
        pockets.append(
            gmsh.model.occ.addBox(
                -x_span / 2,
                -self.wy / 2,
                -z_span / 2,
                x_span, d, z_span
            )
        )

        # Boolean cut
        out, _ = gmsh.model.occ.cut([(3, body_tag)], [(3, p) for p in pockets], removeObject=True, removeTool=True)
        gmsh.model.occ.synchronize()
        return out  # list of (3, tag)

    def _add_tuna_can(self):
        """
        Adds the optional extra volume (tuna can) on the -Z face, centered.
        """
        r = self.tuna_diameter / 2
        L = self.tuna_length

        # Cylinder axis along -Z
        # Base at -wz/2, extending outward to more negative z
        return gmsh.model.occ.addCylinder(0.0, 0.0, -self.wz / 2, 0.0, 0.0, -L, r)

    # ---------- Deployables ----------

    def _add_wings(self):
        """
        Solar panel wings:
        - body_mounted: thin plates flush to +/-X faces
        - deployed: two wings hinged at +/-X, rotated outward 90 deg about Y axis
        """
        wings = []
        if self.wing_config == "none":
            return wings

        t = self.wing_thickness
        L = self.wing_length

        # Panel dimensions: match body height (z) and width (y)
        panel_y = self.wy
        panel_z = self.wz

        if self.wing_config == "body_mounted":
            # Body-mounted: thin plates flush to +/-X faces
            def make_panel():
                # Create panel with thickness along x
                return gmsh.model.occ.addBox(-t / 2, -panel_y / 2, -panel_z / 2, t, panel_y, panel_z)

            # +X face panel
            w = make_panel()
            gmsh.model.occ.translate([(3, w)], self.wx / 2 + t / 2, 0, 0)
            wings.append(w)

            # -X face panel
            w = make_panel()
            gmsh.model.occ.translate([(3, w)], -self.wx / 2 - t / 2, 0, 0)
            wings.append(w)

        elif self.wing_config == "deployed":
            # Deployed: wings rotated 90 degrees outward from +/-X faces
            # Hinge line is vertical (along Z-axis) at the body face
            # Rotation about Z-axis to deploy outward in ±Y direction
            # Use offset to prevent surface overlap - wing starts slightly outside body
            hinge_offset = 0.0001  # 0.1mm offset to prevent surface overlap
            
            for sign in (+1, -1):
                # Create wing box in its initial position (extending along ±X)
                # After rotation, it will extend in ±Y direction
                # Wing dimensions: length L, width panel_y, height panel_z
                
                # Start wing slightly outside body face to avoid exact surface overlap
                wing_x0 = sign * (self.wx / 2 + hinge_offset)
                
                wing = gmsh.model.occ.addBox(
                    wing_x0,                    # x0: slightly outside body face
                    -panel_y / 2,               # y0: centered
                    -panel_z / 2,               # z0: full height
                    sign * L,                   # dx: extend outward (positive for +X, negative for -X)
                    panel_y,                    # dy: full width
                    panel_z                     # dz: full height
                )
                
                # Rotate about Z-axis through the hinge point (at body face)
                # For +X wing: rotate +90° about Z to extend in +Y direction
                # For -X wing: rotate -90° about Z to extend in -Y direction
                hinge_x = sign * (self.wx / 2)   # Hinge at body face (not at wing start)
                hinge_y = 0.0
                hinge_z = 0.0
                rotation_angle = sign * math.pi / 2  # +90° for +X, -90° for -X
                
                gmsh.model.occ.rotate(
                    [(3, wing)],
                    hinge_x, hinge_y, hinge_z,  # Rotation point (hinge at body face)
                    0, 0, 1,                    # Rotation axis (Z-axis, vertical)
                    rotation_angle
                )
                wings.append(wing)

        gmsh.model.occ.synchronize()
        return wings

    def _add_antenna(self, length=0.05, radius=0.002):
        """
        Add whip antenna extending from +Z face.
        Uses slightly larger radius (2mm) for better meshing.
        """
        # Position antenna near corner of +Z face (offset from center)
        # Start from top face, extend upward along +Z
        return gmsh.model.occ.addCylinder(
            self.wx / 4,           # x position (offset from center)
            self.wy / 4,           # y position (offset from center)
            self.wz / 2,           # z start: top of body
            0, 0, length,         # direction: along +Z axis
            radius                 # antenna radius (2mm for better meshing)
        )

    # ---------- Public API ----------

    def generate(self, mesh_size=0.01):
        """
        Generate complete CubeSat geometry with mesh.
        
        Note: For thin features (wings ~2mm, antenna ~2mm radius),
        ensure mesh_size is smaller than feature thickness for proper meshing.
        Recommended: mesh_size <= 0.001 for thin features, or use adaptive meshing.
        """
        # 1) Envelope
        body = self._add_envelope_box()

        # 2) Rails via recess carving
        final_vols = self._carve_rail_recesses(body)

        # 3) Optional extra volume ("tuna can")
        if self.include_tuna_can:
            tuna = self._add_tuna_can()
            out, _ = gmsh.model.occ.fuse(final_vols, [(3, tuna)], removeObject=True, removeTool=True)
            gmsh.model.occ.synchronize()
            final_vols = out

        # 4) Deployables (wings)
        wings = self._add_wings()
        if wings:
            # Fuse wings with body
            wing_entities = [(3, w) for w in wings]
            try:
                # Fuse wings with body
                out, _ = gmsh.model.occ.fuse(final_vols, wing_entities, removeObject=True, removeTool=True)
                gmsh.model.occ.synchronize()
                
                # Remove duplicate entities to fix any overlapping surfaces
                gmsh.model.occ.removeAllDuplicates()
                gmsh.model.occ.synchronize()
                
                final_vols = out
            except Exception as e:
                print(f"Warning: Error fusing wings: {e}")
                print("Trying alternative approach...")
                try:
                    # Alternative: try without removing tools
                    out, _ = gmsh.model.occ.fuse(final_vols, wing_entities, removeObject=False, removeTool=False)
                    gmsh.model.occ.synchronize()
                    # Remove duplicate entities
                    gmsh.model.occ.removeAllDuplicates()
                    gmsh.model.occ.synchronize()
                    final_vols = out
                except Exception as e2:
                    print(f"Warning: Alternative approach also failed: {e2}")
                    print("Continuing without wings...")
                    gmsh.model.occ.synchronize()

        # 5) Optional antenna
        if self.include_antenna:
            ant = self._add_antenna()
            # Fuse antenna with body
            try:
                out, _ = gmsh.model.occ.fuse(final_vols, [(3, ant)], removeObject=True, removeTool=True)
                gmsh.model.occ.synchronize()
                final_vols = out
            except Exception as e:
                print(f"Warning: Error fusing antenna: {e}")
                print("Continuing without antenna...")
                gmsh.model.occ.synchronize()

        # Physical group: all exterior surfaces of final solid(s)
        surfaces = gmsh.model.getBoundary(final_vols, oriented=False)
        surf_tags = [t for (dim, t) in surfaces if dim == 2]
        pg = gmsh.model.addPhysicalGroup(2, surf_tags, 1)
        gmsh.model.setPhysicalName(2, pg, "cubesat_surface")

        # Final geometry healing before meshing
        try:
            gmsh.model.occ.removeAllDuplicates()
            gmsh.model.occ.synchronize()
        except:
            pass

        # Mesh with error checking
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)
        
        # Try to generate mesh, with fallback options
        try:
            gmsh.model.mesh.generate(3)
        except Exception as e:
            error_msg = str(e)
            if "overlapping facets" in error_msg.lower() or "invalid boundary" in error_msg.lower():
                print(f"Warning: Mesh generation failed due to overlapping surfaces: {e}")
                print("Attempting to heal geometry and regenerate mesh...")
                try:
                    # Try healing with larger tolerance
                    gmsh.model.occ.removeAllDuplicates()
                    gmsh.model.occ.synchronize()
                    # Increase mesh size slightly for problematic regions
                    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size * 1.5)
                    gmsh.model.mesh.generate(3)
                    print("Mesh generation succeeded after healing.")
                except Exception as e2:
                    raise Exception(f"Mesh generation failed even after healing: {e2}")
            else:
                raise

        return final_vols

    def save(self, filename):
        gmsh.write(filename)

    def finalize(self):
        gmsh.finalize()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.finalize()


if __name__ == "__main__":
    # Example: realistic 3U with rails + tuna can + deployed wings
    with RealisticCubesatGenerator(
        size="3U",
        wing_config="deployed",
        include_rails=True,
        include_tuna_can=True,
        include_antenna=True,
        recess_depth=0.003,  # tune this if needed
    ) as gen:
        gen.generate(mesh_size=0.01)
        gen.save("cubesat_realistic.msh")
        gen.save("cubesat_realistic.step")