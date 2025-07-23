from init import load_object_grid
import numpy as np
from pyevtk.hl import pointsToVTK, gridToVTK
from pyevtk.vtk import VtkFile, VtkUnstructuredGrid
from os.path import join as pjoin
import h5py
import os

def vtkwrite(path):
    file_name = "particle"#"rhoNeutral" #"P"
    if os.path.exists(pjoin(path,'vtkdata')) == False:
        os.mkdir(pjoin(path,'vtkdata'))
        # write object geometry as an unstructured hexahedral mesh
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
        w.save()
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

    for i in range(len(data_num)):
        datax = np.array(datapos[i,:,0])
        datay = np.array(datapos[i,:,1])
        dataz = np.array(datapos[i,:,2])
        # write points with embedded radius data
        pointsToVTK(
            pjoin(path,'vtkdata','points_%d'%i),
            datax, datay, dataz,
            data={'radius': radii}
        )
