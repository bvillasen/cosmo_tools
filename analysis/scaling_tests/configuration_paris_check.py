import numpy as np



grid = np.array([ 2048, 1024, 25600  ]).astype(np.float)
mpi  = np.array([ 16, 8, 200  ])


grid = np.array([ 2048, 4096, 51200  ]).astype(np.float)
mpi  = np.array([ 8, 16, 200  ])


n_per_gpu = grid / mpi
nx, ny, nz = grid
nxp = ( nx/2 + 1)
n_slab = nxp * ny

n_mpi = mpi[0]*mpi[1]*mpi[2]
print n_mpi

print n_per_gpu
print n_slab / n_mpi 
print nz / n_mpi


