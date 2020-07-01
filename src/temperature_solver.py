import numpy as np
from conductionNd_serial import ConductionND
from scipy.spatial import cKDTree
import time

try: range=xrange
except: pass


def solve(self, matrix=None, rhs=None, tol=1e-15):
    """
    Construct the matrix A and vector b in Ax = b
    and solve for x

    GMRES method is default
    """
    from scipy.sparse import linalg
    if matrix is None:
        matrix = self.construct_matrix()
    if rhs is None:
        rhs = self.construct_rhs()
    res = self.temperature
    x0 = self.temperature.copy()

    T, info = linalg.cg(matrix, rhs, tol=tol, x0=x0)
    self.temperature = T
    return T

def hofmeister1999(k0, T):
    a = 2.85e-5
    b = 1.008e-8
    c = -0.384
    d = 0.0

    k_con = k0*(298.0/T)**0.45
    k_rad = 0.0175 - 1.0374e-4*T + 2.245*T**2/1e7 - 3.407*T**3/1e11
    k_pres= 1.0 + (K_0p*g_av*3350.0*np.abs(mesh.coords[:,2])*1e-9)/K_T
    k_exp = np.exp(-(a*(T - 298.0) + b*0.5*(T**2 - 88804.0) + \
                     c*(3.3557e-3 - 1.0/T) + d*(T**5 - 2.35e12)/5.0)* \
                     (gmma_T*4 + 1.0/3))
    k_new = k_con*k_exp*k_pres + k_rad
    k_new[grun_mask] = k0[grun_mask]
    k_new[kmask] = k0[kmask]
    return k_new

def nonlinear_conductivity(self, k0, tolerance):
    k = k0.copy()
    rhs = self.rhs

    error = np.array(10.0)
    i = 0

    while (error > tolerance):
        t = time.time()
        k_last = self.diffusivity.copy()
        self.diffusivity = k
        mat = self.construct_matrix()

        T = solve(self, matrix=mat, rhs=rhs)
        k = hofmeister1999(k0, T)

        error = np.absolute(k - k_last).max()
        i += 1

        print("iteration {} in {:.3f} seconds, residual = {:.2e}".format(i, time.time()-t, float(error)))


def map(lithology, lithology_index, sort_idx, *args):
    """
    Requires a tuple of vectors corresponding to an inversion variable
    these are mapped to the mesh.

    tuple(vec1, vec2, vecN) --> tuple(field1, field2, fieldN)
    """

    nf = len(args)
    nl = len(lithology_index)

    # preallocate memory
    mesh_variables = np.zeros((nf, lithology.size))

    # unpack vector to field
    for i in range(0, nl):
        idx = lithology == lithology_index[i]
        for f in range(0, nf):
            mesh_variables[f,idx] = args[f][i]

    # sort indices for local mesh
    for f in range(0, nf):
        mesh_variables[f] = mesh_variables[f][sort_idx]

    return list(mesh_variables)


systime = time.time()

## Read from file ##
# read properties from layers.info
layer_attributes = np.loadtxt('layers.info', skiprows=1, usecols=list(range(2,11)))
lithology_index  = np.loadtxt('layers.info', skiprows=1, usecols=(0,), dtype=int)
nonlinear = (layer_attributes[:,6] > 0).any() # Gruneisen parameter


# litmod_parameters = np.loadtxt('litmod3d.info')
with open('LITMOD3D.info', 'r') as litmod_parameters:
    line = litmod_parameters.readline()
    line = line.replace('D', '.')
    linesplit = line.split()
    print(line)
    
    
    T0 = float(linesplit[9])  + 273.14
    T1 = float(linesplit[10]) + 273.14
    
# T1 = 1315.0 + 273.14
# T0 = 273.14

# read thermal properties and coordinates from xyz files
directory = 'OUT_xyz/'

xyz = np.loadtxt(directory + 'Material_index.xyz')
xyz_coords = xyz[:,:3] * 1e3
lithology  = xyz[:,3].astype(int)


## Extract coordinates ##
Xcoords = np.unique(xyz_coords[:,0])
Ycoords = np.unique(xyz_coords[:,1])
Zcoords = np.unique(xyz_coords[:,2])

minX, maxX = Xcoords[0], Xcoords[-1]
minY, maxY = Ycoords[0], Ycoords[-1]
minZ, maxZ = Zcoords[0], Zcoords[-1]

nx, ny, nz = Xcoords.size, Ycoords.size, Zcoords.size



## Initialise mesh ##
mesh = ConductionND((minX, minY, minZ), (maxX, maxY, maxZ), (nx, ny, nz))

# determine how to sort fortran arrays to C arrays
tree = cKDTree(xyz_coords)
d, idx = tree.query(mesh.coords)

# map properties
rho, alpha, k0, H, beta, gmma_T, K_0p, K_T, man = map(lithology, lithology_index, idx, *layer_attributes.T)
grun_mask = gmma_T == 0
K_T[K_T == 0] = 1e99


# gravity constants
g_s = 9.81
g_400 = 9.97
depth = abs(minZ)
d_gz = (g_400 - g_s)/400e3
g_av = g_s + d_gz*depth*0.5 # Average gravity attraction for the thermal calculation



kmask = k0 == 0.0
# k2[kmask] = 10.0

# we differentiate air and asthenosphere by 20km depth
air_mask = np.logical_and(kmask, mesh.coords[:,2]>-20e3)
lab_mask = np.logical_and(kmask, mesh.coords[:,2]<-20e3)


# Cut mask
air_mask = air_mask.reshape(mesh.n)
ii0, jj0, kk0 = np.where(np.diff(air_mask.astype(int), axis=0) != 0) # find intersection
air_mask[ii0,jj0,kk0] = True
air_mask = air_mask.ravel()

lab_mask = lab_mask.reshape(mesh.n)
ii0, jj0, kk0 = np.where(np.diff(lab_mask.astype(int), axis=0) != 0) # find intersection
lab_mask[ii0+1,jj0,kk0] = True
lab_mask[ii0-1,jj0,kk0] = True
lab_mask[ii0,jj0,kk0] = True
lab_mask = lab_mask.ravel()


mesh.update_properties(k0, H)



mesh.boundary_condition('maxZ', T0, flux=False)
mesh.boundary_condition('minZ', T1, flux=False)

rhs = mesh.construct_rhs()
mesh.dirichlet_mask[air_mask] = True
mesh.dirichlet_mask[lab_mask] = True

air_idx = np.nonzero(air_mask)[0].astype(np.int32)
lab_idx = np.nonzero(lab_mask)[0].astype(np.int32)

rhs[air_idx] = np.ones_like(air_idx)*T0
rhs[lab_idx] = np.ones_like(lab_idx)*T1
mesh.rhs = rhs

# x0
mesh.temperature[:] = mesh.rhs[:]


# solve temperature
walltime = time.time()
if nonlinear:
    nonlinear_conductivity(mesh, k0, 1e-6)
else:
    solve(mesh, rhs=rhs)

print("Total solve time is {:.3f} seconds".format(time.time()-walltime))




# Adiabatic correction beneath LAB
Ta = T1 - 273.14
Tbottom = 1520.

T = mesh.temperature[:].copy().reshape(mesh.n)
zcube = mesh.coords[:,-1].reshape(mesh.n)

sublit_mask = lab_mask.reshape(mesh.n) #(T >= T1).reshape(mesh.n)
ii0, jj0, kk0 = np.where(np.diff(sublit_mask.astype(int), axis=0) != 0) # find intersection

# how many cells is 40 km
nz = mesh.n[0]
dz = np.diff(mesh.grid_coords[-1]).mean()
nc = int(40e3/dz)

z_lit = zcube[ii0,jj0,kk0]
# if mesh.grid_coords[-1].min() >= -60e3:
#     z_lit_buf = -100e3
z_lit_buf = z_lit - 40e3


# set temperature inside buffer
sublit_buffer_mask = np.zeros(mesh.n, dtype=bool)
for row in range(nc):
    # clip if buffer exceeds domain
    ii1 = np.maximum(ii0 - row, 0)
    sublit_buffer_mask[ii1,jj0,kk0] = True
    T[ii1,jj0,kk0] = (Ta + 273.14) + (zcube[ii1,jj0,kk0] - z_lit)*(1400.0-Ta)/(z_lit_buf - z_lit)
    
# set temperature below buffer
grd_ad = (1520.0 - 1400.0)/(400e3 + z_lit_buf)
grd_ad = np.clip(grd_ad, 0.35e-3, 0.6e-3)

while (ii1 > 0).any():
    row += 1
    ii1 = np.maximum(ii0 - row, 0)
    T[ii1,jj0,kk0] = (1400.0 + 273.14) - (zcube[ii1,jj0,kk0] - z_lit_buf)*grd_ad

mesh.temperature[:] = T.ravel()


# mesh.diffusivity[air_idx] = 0.0
# mesh.diffusivity[lab_idx] = 0.0


# Setup second KDTree to consistently interpolate back to original ordering
tree = cKDTree(mesh.coords)
d, idx = tree.query(xyz_coords)
global_T = mesh.temperature[idx]
global_k = mesh.diffusivity[idx]


# global_T = mesh.temperature.reshape(nz,ny,nx).ravel(order='F')
# global_k = mesh.diffusivity.reshape(nz,ny,nx).ravel(order='F')

xyz_coords *= 1e-3
xyz_temperature  = np.column_stack([xyz_coords, global_T - 273.14])
xyz_conductivity = np.column_stack([xyz_coords, global_k])

# write to file
np.savetxt(directory + 'Temperature.xyz', xyz_temperature, fmt='%13.3f %12.3f %12.3f %9.3f')
np.savetxt(directory + 'Thermal_cond.xyz', xyz_conductivity, fmt='%13.3f %12.3f %12.3f %9.2f')

print("Total system time is {:.3f} seconds".format(time.time()-systime))
