from init import load_object_grid
import numpy as np
from pyevtk.hl import pointsToVTK, gridToVTK
from pyevtk.vtk import VtkFile, VtkUnstructuredGrid
from os.path import join as pjoin
import h5py
import os
import config

def load_gmsh_mesh(mesh_file):
    """
    Load Gmsh mesh file and extract surface triangles.
    Returns: (points, triangles) where points is (N, 3) and triangles is (M, 3) indices
    """
    try:
        import gmsh
    except ImportError:
        return None, None

    if not os.path.exists(mesh_file):
        return None, None

    try:
        gmsh.initialize()
        gmsh.open(mesh_file)

        # Get all nodes
        node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
        # node_coords is flattened: [x1, y1, z1, x2, y2, z2, ...]
        n_nodes = len(node_tags)
        points = np.array(node_coords).reshape(n_nodes, 3)

        # Create mapping from node tag to index
        tag_to_idx = {tag: idx for idx, tag in enumerate(node_tags)}

        # Get surface triangles (2D elements)
        element_types, element_tags, element_node_tags = gmsh.model.mesh.getElements(dim=2)

        triangles = []
        for elem_type, elem_tags, node_tags in zip(element_types, element_tags, element_node_tags):
            if elem_type == 2:  # Triangle (3-node)
                # node_tags is flattened: [n1, n2, n3, n4, n5, n6, ...]
                n_elems = len(elem_tags)
                for i in range(n_elems):
                    idx = i * 3
                    n1, n2, n3 = node_tags[idx], node_tags[idx+1], node_tags[idx+2]
                    triangles.append([tag_to_idx[n1], tag_to_idx[n2], tag_to_idx[n3]])

        gmsh.finalize()

        if len(triangles) == 0:
            return None, None

        return points, np.array(triangles, dtype=np.int32)
    except Exception as e:
        print(f"Warning: Could not load Gmsh mesh: {e}")
        try:
            gmsh.finalize()
        except:
            pass
        return None, None

def vtkwrite(path):
    file_name = "particle"#"rhoNeutral" #"P"
    # Always create vtkdata directory
    os.makedirs(pjoin(path,"vtkdata"), exist_ok=True)

    # Always write object geometry (it may have updated heatmap data)
    os.makedirs(pjoin(path,"vtkdata"), exist_ok=True)

    # Try to use Gmsh mesh if available and use_mesh_directly is enabled
    mesh_file = pjoin(config.picDir, 'cubesat.msh')
    use_gmsh_mesh = getattr(config, 'use_mesh_directly', False) and os.path.exists(mesh_file)

    if use_gmsh_mesh:
        print("Using Gmsh mesh for object visualization...")
        points, triangles = load_gmsh_mesh(mesh_file)

        if points is not None and triangles is not None:
            # Write Gmsh mesh as VTK unstructured grid
            n_points = points.shape[0]
            n_cells = triangles.shape[0]

            # Flatten triangle connectivity
            connectivity = triangles.flatten()
            offsets = np.arange(3, (n_cells + 1) * 3, 3, dtype=np.int32)
            cell_types = np.ones(n_cells, dtype=np.uint8) * 5  # VTK_TRIANGLE

            # Load impact heatmap if available
            h5_file = pjoin(path, file_name+'.hdf5')
            impact_heatmap = None
            xg_obj, yg_obj, zg_obj = None, None, None
            if os.path.exists(h5_file):
                with h5py.File(h5_file, 'r') as h5_temp:
                    if 'impact_heatmap' in h5_temp.keys():
                        impact_heatmap = h5_temp['impact_heatmap'][:]
                        # Also load grid coordinates for mapping
                        from init import load_object_grid
                        xg_obj, yg_obj, zg_obj, _ = load_object_grid()

            filepath = pjoin(path, 'vtkdata', 'object_mesh')
            w = VtkFile(filepath, VtkUnstructuredGrid)
            w.openGrid()
            w.openPiece(ncells=n_cells, npoints=n_points)
            w.openElement("Points")
            w.addData("points", (points[:,0], points[:,1], points[:,2]))
            w.closeElement("Points")
            w.openElement("Cells")
            w.addData("connectivity", connectivity)
            w.addData("offsets", offsets)
            w.addData("types", cell_types)
            w.closeElement("Cells")

            # Map heatmap from voxel grid to mesh triangles
            if impact_heatmap is not None and xg_obj is not None:
                # Compute triangle centers and map to voxel grid
                triangle_centers = np.zeros((n_cells, 3))
                heatmap_values = np.zeros(n_cells)
                for i, tri in enumerate(triangles):
                    # Triangle center
                    center = np.mean(points[tri], axis=0)
                    triangle_centers[i] = center

                    # Find voxel index for this center point
                    # Use searchsorted to find grid indices
                    ix = np.searchsorted(xg_obj, center[0]) - 1
                    iy = np.searchsorted(yg_obj, center[1]) - 1
                    iz = np.searchsorted(zg_obj, center[2]) - 1

                    # Clamp to valid range
                    ix = max(0, min(ix, impact_heatmap.shape[0] - 1))
                    iy = max(0, min(iy, impact_heatmap.shape[1] - 1))
                    iz = max(0, min(iz, impact_heatmap.shape[2] - 1))

                    heatmap_values[i] = impact_heatmap[ix, iy, iz]

                w.openElement("CellData")
                w.addData("impact_heatmap", heatmap_values)
                w.closeElement("CellData")

            w.closePiece()
            w.closeGrid()

            # Append data
            pts_x = np.ascontiguousarray(points[:, 0])
            pts_y = np.ascontiguousarray(points[:, 1])
            pts_z = np.ascontiguousarray(points[:, 2])
            conn = np.ascontiguousarray(connectivity)
            offs = np.ascontiguousarray(offsets)
            ctypes = np.ascontiguousarray(cell_types)
            w.appendData((pts_x, pts_y, pts_z))
            w.appendData(conn)
            w.appendData(offs)
            w.appendData(ctypes)
            if impact_heatmap is not None and xg_obj is not None:
                heatmap_vals = np.ascontiguousarray(heatmap_values)
                w.appendData(heatmap_vals)
            w.save()
            print(f"  Exported {n_cells} triangles from Gmsh mesh")
            if impact_heatmap is not None:
                print(f"  Impact heatmap mapped to mesh (max={np.max(heatmap_values):.1f})")
        else:
            print("  Warning: Could not load Gmsh mesh, falling back to voxel grid")
            use_gmsh_mesh = False

    if not use_gmsh_mesh:
        # write object geometry as an unstructured hexahedral mesh (voxel grid)
        xg, yg, zg, mask = load_object_grid()
        # find occupied voxels
        ixs, iys, izs = np.where(mask == 1)
        ncells = ixs.size
        # allocate array for all corner points
        pts = np.zeros((ncells * 8, 3))
        for n, (i, j, k) in enumerate(zip(ixs, iys, izs)):
            # compute voxel corner coords (handle grid edges)
            x0, x1 = xg[i], xg[i+1] if i+1 < xg.size else xg[i] + (xg[i] - xg[i-1])
            y0, y1 = yg[j], yg[j+1] if j+1 < yg.size else yg[j] + (yg[j] - yg[j-1])
            z0, z1 = zg[k], zg[k+1] if k+1 < zg.size else zg[k] + (zg[k] - zg[k-1])
            corners = np.array([
                [x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
                [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]
            ])
            pts[n*8:(n+1)*8, :] = corners
        # build connectivity arrays for unstructured mesh
        connectivity = np.arange(ncells * 8, dtype=np.int32)
        offsets = np.arange(1, ncells + 1, dtype=np.int32) * 8
        cell_types = np.ones(ncells, dtype=np.uint8) * 12  # VTK_HEXAHEDRON
        # write object geometry as an unstructured hexahedral mesh
        filepath = pjoin(path, 'vtkdata', 'object_mesh')
        w = VtkFile(filepath, VtkUnstructuredGrid)
        w.openGrid()
        w.openPiece(ncells=ncells, npoints=pts.shape[0])
        w.openElement("Points")
        w.addData("points", (pts[:,0], pts[:,1], pts[:,2]))
        w.closeElement("Points")
        w.openElement("Cells")
        w.addData("connectivity", connectivity)
        w.addData("offsets", offsets)
        w.addData("types", cell_types)
        w.closeElement("Cells")
        # Add impact heatmap data if available (BEFORE closePiece!)
        # Check for impact heatmap data
        h5_file = pjoin(path, file_name+'.hdf5')
        impact_heatmap = None
        if os.path.exists(h5_file):
            with h5py.File(h5_file, 'r') as h5_temp:
                if 'impact_heatmap' in h5_temp.keys():
                    impact_heatmap = h5_temp['impact_heatmap'][:]

        if impact_heatmap is not None:
            # Extract heatmap values for occupied voxels
            heatmap_values = impact_heatmap[ixs, iys, izs]
            w.openElement("CellData")
            w.addData("impact_heatmap", heatmap_values)
            w.closeElement("CellData")
        w.closePiece()
        w.closeGrid()
        # append the binary data for points and cells
        # ensure arrays are contiguous
        pts_x = np.ascontiguousarray(pts[:, 0])
        pts_y = np.ascontiguousarray(pts[:, 1])
        pts_z = np.ascontiguousarray(pts[:, 2])
        conn = np.ascontiguousarray(connectivity)
        offs = np.ascontiguousarray(offsets)
        ctypes = np.ascontiguousarray(cell_types)
        w.appendData((pts_x, pts_y, pts_z))
        w.appendData(conn)
        w.appendData(offs)
        w.appendData(ctypes)
        if impact_heatmap is not None:
            heatmap_vals = np.ascontiguousarray(heatmap_values)
            w.appendData(heatmap_vals)
        w.save()
    # Check if HDF5 file exists
    h5_file_path = pjoin(path, file_name+'.hdf5')
    if not os.path.exists(h5_file_path):
        print(f'Warning: HDF5 file not found: {h5_file_path}')
        return
    h5 = h5py.File(pjoin(path,file_name+'.hdf5'),'r')

    # load particle radii
    radii = h5['radius'][:]  # shape (N,)

    Lx = h5.attrs["Lx"]
    Ly = h5.attrs["Ly"]
    Lz = h5.attrs["Lz"]

    dp   = h5.attrs["dp"]
    Nt   = h5.attrs["Nt"]

    data_num = np.arange(start=0, stop=Nt, step=dp, dtype=int)

    datapos = h5["/position"]

    # Check if particle_status exists in HDF5
    particle_status = None
    if 'particle_status' in h5.keys():
        particle_status = h5['particle_status'][:]

    for i in range(len(data_num)):
        datax = np.array(datapos[i,:,0])
        datay = np.array(datapos[i,:,1])
        dataz = np.array(datapos[i,:,2])

        # Filter out inactive particles (absorbed/attached) if status available
        if particle_status is not None:
            active_mask = (particle_status == 0)
            # Only include active particles
            datax = datax[active_mask]
            datay = datay[active_mask]
            dataz = dataz[active_mask]
            radii_filtered = radii[active_mask]
        else:
            # Filter out NaN positions (inactive particles)
            valid_mask = ~(np.isnan(datax) | np.isnan(datay) | np.isnan(dataz))
            datax = datax[valid_mask]
            datay = datay[valid_mask]
            dataz = dataz[valid_mask]
            radii_filtered = radii[valid_mask] if len(radii) == len(valid_mask) else radii

        # Only write if there are active particles
        if len(datax) > 0:
            pointsToVTK(
                pjoin(path,'vtkdata','points_%d'%i),
                datax, datay, dataz,
                data={'radius': radii_filtered}
            )
