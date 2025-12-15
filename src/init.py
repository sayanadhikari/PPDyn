from numba import jit
import numpy as np
import random
import h5py
import shutil
import os

import config


def copy_pic_data():
    # Ensure destination directory exists
    destination_dir = 'PIC_data'
    os.makedirs(destination_dir, exist_ok=True)

    # Define source files
    source_files = [
        os.path.join(config.picDir, 'E.grid.h5'),
        os.path.join(config.picDir, 'object.grid.h5')
    ]

    # Copy each file
    for src in source_files:
        if os.path.isfile(src):
            shutil.copy(src, destination_dir)
            print(f"Copied {src} to {destination_dir}")
        else:
            print(f"Warning: Source file does not exist: {src}")

def load_pic_field():
    # ------------------------------------------------------------------
    # Load PIC field from HDF5
    # ------------------------------------------------------------------
    file_path = os.path.join(config.picDir, "E.grid.h5")
    with h5py.File(file_path, "r") as f:
        denorm = f.attrs["Axis denormalization factor"][0]
        # Get the last n entry (sorted by key)
        keys = sorted(f.keys())
        last_key = keys[-1]
        # grid_shape = f[last_key].shape
        # config.Lx = 0.5 * grid_shape[0] * denorm
        # config.Ly = 0.5 * grid_shape[1] * denorm
        # config.Lz = 0.5 * grid_shape[2] * denorm
        Ex = f[last_key][:, :, :, 0]
        Ey = f[last_key][:, :, :, 1]
        Ez = f[last_key][:, :, :, 2]

        # Ey = f["/Ey"][:]
        # Ez = f["/Ez"][:]
        
        # Normalize PIC electric field by selected method: 'aOmega' or 'max'
        norm_method = getattr(config, 'E_norm_method', 'aOmega')
        print("E_norm_method:", norm_method)
    
        # E_mag = np.sqrt(Ex**2 + Ey**2 + Ez**2)
           # # grid coordinates
        dx = config.dx_a
        xg = np.arange(Ex.shape[0])*  dx - config.Lx
        yg = np.arange(Ey.shape[1])*  dx  - config.Ly
        zg = np.arange(Ez.shape[2])*  dx   - config.Lz
        # norm_method = getattr(config, 'E_norm_method', 'mean')
        # if norm_method == 'max':
        if norm_method == 'aOmega':
            Ex *= config.Escale_a
            Ey *= config.Escale_a
            Ez *= config.Escale_a
        elif norm_method == 'max':
            E_mag = np.sqrt(Ex**2 + Ey**2 + Ez**2)
            E0 = np.max(E_mag)
            if E0 != 0:
                Ex /=E0; Ey /=E0; Ez /=E0
        # elif norm_method == 'norm':
        #     # Birdsall normalization: E0 = k_B * Te / (e * lambda_D)
        #     if config.lambda_D_SI > 0 and (config.Te_J > 0):  # J/K
        #         E0 = config.Te_J / (config.qe * config.lambda_D_SI)
        # else:
        #     E_mag = np.sqrt(Ex**2 + Ey**2 + Ez**2)
        #     E0 = np.mean(E_mag)
        #     if E0 != 0:
        #         Ex /= E0; Ey /= E0; Ez /= E0

        return xg, yg, zg, Ex, Ey, Ez

def load_object_grid():
    """
    Load the object mask and its axis denormalization factor from HDF5.
    Returns:
      xg_obj, yg_obj, zg_obj: 1D arrays of grid coordinates
      obj_mask     : 3D uint8 array where 1=inside object, 0=outside
    """
    # file_path = os.path.join(config.picDir, "object.grid.h5")
    # with h5py.File(file_path, 'r') as f:
    #     dset = f['Object']              # dataset holding binary mask
    #     obj_mask    = dset[:].astype(np.uint8)
    #     axis_denorm = f.attrs['Axis denormalization factor'][0]
    #     # The object grid is assumed cubic, so all axes have the same shape
    #     nx, ny, nz = obj_mask.shape
    #     # map object grid [0, axis_denorm*nx) onto MD domain [-Lx, +Lx)
    #     xg_obj = np.arange(nx) * axis_denorm - config.Lx
    #     yg_obj = np.arange(ny) * axis_denorm - config.Ly
    #     zg_obj = np.arange(nz) * axis_denorm - config.Lz
    # return xg_obj, yg_obj, zg_obj, obj_mask

    file_path = os.path.join(config.picDir, "object.grid.h5")
    with h5py.File(file_path, 'r') as f:
        dset = f['Object']              # dataset holding binary mask
        obj_mask    = dset[:].astype(np.uint8)
        axis_denorm = f.attrs['Axis denormalization factor'][0]
            # The object grid is assumed cubic, so all axes have the same shape
        nx, ny, nz = obj_mask.shape
        dx_obj = axis_denorm/config.a
            # map object grid [0, axis_denorm*nx) onto MD domain [-Lx, +Lx)
        # xg_obj = np.arange(nx) * axis_denorm - config.Lx
        # yg_obj = np.arange(ny) * axis_denorm - config.Ly
        # zg_obj = np.arange(nz) * axis_denorm - config.Lz
        xg_obj = np.arange(nx) * dx_obj- config.Lx
        yg_obj = np.arange(ny) * dx_obj - config.Ly
        zg_obj = np.arange(nz) * dx_obj - config.Lz
    return xg_obj, yg_obj, zg_obj, obj_mask
    
@jit(nopython=True)
def initial_periodic(Q, M, xg_obj, yg_obj, zg_obj, obj_mask):
    random.seed(99999999)
    pos   = np.empty((config.N,3), dtype=np.float64)
    uvel  = np.empty((config.N,3), dtype=np.float64)
    vvel  = np.empty((config.N,3), dtype=np.float64)
    acc   = np.empty((config.N,3), dtype=np.float64)
    sv    = np.zeros((config.N,3), dtype=np.float64)

    # svx  = 0.0  # velocity sum correction term in X
    # svy  = 0.0  # velocity sum correction term in Y
    # svz  = 0.0  # velocity sum correction term in Z

    ###### Initialize time array and data dump array ######
    time  = np.linspace(0,config.tmax,config.Nt)
    data_num = np.arange(start=0, stop=config.Nt, step=config.dumpPeriod, dtype=np.int64)
    # initialize positions avoiding object geometry
    dx = xg_obj[1] - xg_obj[0]
    dy = yg_obj[1] - yg_obj[0]
    dz = zg_obj[1] - zg_obj[0]
    for i in range(config.N):
        placed = False
        while not placed:
            px = np.random.random() * 2.0 * config.Lx - config.Lx
            py = np.random.random() * 2.0 * config.Ly - config.Ly
            pz = np.random.random() * 2.0 * config.Lz - config.Lz
            # compute mask indices
            ix = int(round((px - xg_obj[0]) / dx))
            iy = int(round((py - yg_obj[0]) / dy))
            iz = int(round((pz - zg_obj[0]) / dz))
            # clamp
            if ix < 0: ix = 0
            elif ix >= obj_mask.shape[0]: ix = obj_mask.shape[0]-1
            if iy < 0: iy = 0
            elif iy >= obj_mask.shape[1]: iy = obj_mask.shape[1]-1
            if iz < 0: iz = 0
            elif iz >= obj_mask.shape[2]: iz = obj_mask.shape[2]-1
            if obj_mask[ix, iy, iz] == 0:
                pos[i,0] = px
                pos[i,1] = py
                pos[i,2] = pz
                placed = True

    # Maxwellian
    if config.maxwell_load:
        vvel[:,0] = np.random.normal(0, config.Temp, config.N)
        vvel[:,1] = np.random.normal(0, config.Temp, config.N)
        vvel[:,2] = np.random.normal(0, config.Temp, config.N)
    # Random
    else:
        vvel[:,0] = np.random.random(config.N)*config.Vxmax - config.Vxmax/2.0
        vvel[:,1] = np.random.random(config.N)*config.Vymax - config.Vymax/2.0
        vvel[:,2] = np.random.random(config.N)*config.Vzmax - config.Vzmax/2.0

        sv[:,0] = sv[:,0] + vvel[:,0]
        sv[:,1] = sv[:,1] + vvel[:,1]
        sv[:,2] = sv[:,2] + vvel[:,2]

        vvel[:,0] = vvel[:,0] - sv[:,0]/config.N
        vvel[:,1] = vvel[:,1] - sv[:,1]/config.N
        vvel[:,2] = vvel[:,2] - sv[:,2]/config.N

    # acc = 0.0*acc

    for i in range(config.N):
        acc[i,:] = 0.0
        for j in range(config.N):
            if (i != j):
                xdiff = ( pos[i,0]-pos[j,0] ) - round((pos[i,0]-pos[j,0])/(2.0*config.Lx)) * 2.0*config.Lx
                ydiff = ( pos[i,1]-pos[j,1] ) - round((pos[i,1]-pos[j,1])/(2.0*config.Ly)) * 2.0*config.Ly
                zdiff = ( pos[i,2]-pos[j,2] ) - round((pos[i,2]-pos[j,2])/(2.0*config.Lz)) * 2.0*config.Lz
                r = np.sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)
                fx = xdiff*(1+config.k*r)*np.exp(-config.k*r)*(Q[i]*Q[j])/(r*r*r)    # xdiff/(r*r*r)
                fy = ydiff*(1+config.k*r)*np.exp(-config.k*r)*(Q[i]*Q[j])/(r*r*r)    # ydiff/(r*r*r)
                fz = zdiff*(1+config.k*r)*np.exp(-config.k*r)*(Q[i]*Q[j])/(r*r*r) #+ zdiff*g + Lz*g  # zdiff/(r*r*r)
                acc[i,0] += fx/M[i]
                acc[i,1] += fy/M[i]
                acc[i,2] += fz/M[i]
    return pos,vvel,uvel,acc,time,data_num


@jit(nopython=True)
def initial_reflecting(Q,M):
    random.seed(99999999)
    pos   = np.empty((config.N,3), dtype=np.float64)
    uvel  = np.empty((config.N,3), dtype=np.float64)
    vvel  = np.empty((config.N,3), dtype=np.float64)
    acc   = np.empty((config.N,3), dtype=np.float64)
    sv    = np.zeros((config.N,3), dtype=np.float64)


    # svx  = 0.0  # velocity sum correction term in X
    # svy  = 0.0  # velocity sum correction term in Y
    # svz  = 0.0  # velocity sum correction term in Z

    ###### Initialize time array and data dump array ######
    time  = np.linspace(0,config.tmax,config.Nt)
    data_num = np.arange(start=0, stop=config.Nt, step=config.dumpPeriod, dtype=np.int64)

    pos[:,0] = np.random.random(config.N)*2.0*config.Lx - config.Lx
    pos[:,1] = np.random.random(config.N)*2.0*config.Ly - config.Ly
    pos[:,2] = np.random.random(config.N)*2.0*config.Lz - config.Lz

    # Maxwellian
    if config.maxwell_load:
        vvel[:,0] = np.random.normal(0, config.Temp, config.N)
        vvel[:,1] = np.random.normal(0, config.Temp, config.N)
        vvel[:,2] = np.random.normal(0, config.Temp, config.N)
    # Random
    else:
        vvel[:,0] = np.random.random(config.N)*config.Vxmax - config.Vxmax/2.0
        vvel[:,1] = np.random.random(config.N)*config.Vymax - config.Vymax/2.0
        vvel[:,2] = np.random.random(config.N)*config.Vzmax - config.Vzmax/2.0

        sv[:,0] = sv[:,0] + vvel[:,0]
        sv[:,1] = sv[:,1] + vvel[:,1]
        sv[:,2] = sv[:,2] + vvel[:,2]

        vvel[:,0] = vvel[:,0] - sv[:,0]/config.N
        vvel[:,1] = vvel[:,1] - sv[:,1]/config.N
        vvel[:,2] = vvel[:,2] - sv[:,2]/config.N

    for i in range(config.N):
        acc[i,0] = 0.0
        acc[i,1] = 0.0
        acc[i,2] = -(pos[i,2]+config.Lz)*config.g
        for j in range(config.N):
            if (i != j):
                xdiff = ( pos[i,0]-pos[j,0] )
                ydiff = ( pos[i,1]-pos[j,1] )
                zdiff = ( pos[i,2]-pos[j,2] )
                r = np.sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)
                fx = xdiff*(1+config.k*r)*np.exp(-config.k*r)*(Q[i]*Q[j])/(r*r*r)    # xdiff/(r*r*r)
                fy = ydiff*(1+config.k*r)*np.exp(-config.k*r)*(Q[i]*Q[j])/(r*r*r)    # ydiff/(r*r*r)
                fz = zdiff*(1+config.k*r)*np.exp(-config.k*r)*(Q[i]*Q[j])/(r*r*r) # + zdiff*g + Lz*g # zdiff/(r*r*r)
                acc[i,0] += fx/M[i]
                acc[i,1] += fy/M[i]
                acc[i,2] += fz/M[i]
    return pos,vvel,uvel,acc,time,data_num
