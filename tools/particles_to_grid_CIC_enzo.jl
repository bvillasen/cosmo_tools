using HDF5
current_dir = pwd()
tools_dir = current_dir
push!(LOAD_PATH, tools_dir)
using CIC_functions
using Statistics

# Input File
# dataDir = "/raid/bruno/data/"
# dataDir = "/home/bruno/Desktop/data/"
# dataDir = "/home/bruno/Desktop/data/"
dataDir = "/home/bruno/Desktop/ssd_0/data/"
inDir = dataDir * "cosmo_sims/enzo/256_hydro_50Mpc/h5_files/"
# inDir = "/home/bruno/Desktop/data/cosmo_sims/ramses/128_hydro/h5_files/"
# inDir = dataDir * "cosmo_sims/enzo/256_hydro_grackle_noUV/h5_files/"
outDir = inDir
in_base_name = "snapshot_"
out_base_name = "grid_CIC_"

Lbox = 50e3
const nPoints = 256

#Domain Parameters
const x_min = 0.0
const y_min = 0.0
const z_min = 0.0
const x_max = Lbox
const y_max = Lbox
const z_max = Lbox
const Lx = x_max - x_min
const Ly = y_max - y_min
const Lz = z_max - z_min

#Grid Properties
const nx = nPoints
const ny = nPoints
const nz = nPoints
const nCells = nx*ny*nz
const dx = Lx / nx
const dy = Ly / ny
const dz = Lz / nz


# #Cells positions in grid ( mid-point )
# c_pos_x = linspace( x_min + dx/2, x_max - dx/2, nx)
# c_pos_y = linspace( y_min + dy/2, y_max - dy/2, ny)
# c_pos_z = linspace( z_min + dz/2, z_max - dz/2, nz)


nSnap = 0
for nSnap in 0:30

  println( "\nSnapshot: $(nSnap)")
  snapKey = lpad(nSnap,3,'0')
  inFileName = inDir * in_base_name * snapKey * ".h5"
  outFileName = outDir * out_base_name * snapKey * ".h5"

  print(" Loading File: $(inFileName)\n")
  inFile = h5open( inFileName, "r")


  current_a = read( attrs(inFile), "current_a" )
  current_z = read( attrs(inFile), "current_z" )
  print(" Current redshift: $(current_z)\n")

  print(" Writing File: $(outFileName)\n")
  outFile = h5open( outFileName, "w")
  attrs(outFile)["current_a"] = current_a
  attrs(outFile)["current_z"] = current_z

  # part_type = "dm"
  for part_type in [ "dm"]
    println(" $(part_type)")
    data = inFile[part_type]
    # mass = read( data, "mass" )
    pos_x = read( data, "pos_x" )
    pos_y = read( data, "pos_y" )
    pos_z = read( data, "pos_z" )

    if part_type == "gas"
      field = read(data, "density") *dx*dy*dz
    else
      field = read( data, "mass" )
    end


    nParticles = size( field )[1]
    println( "  N parts: $(nParticles)")
    p_inside = ones( Bool, nParticles )
    get_particles_outside_CIC( p_inside, pos_x, pos_y, pos_z, x_min, x_max, y_min, y_max, z_min, z_max, dx, dy, dz  )
    dens = get_interp_CIC( p_inside, field, field, field, pos_x, pos_y, pos_z, nx, ny, nz, x_min, y_min, z_min, x_max, y_max, z_max, dx, dy, dz,  true, true )
    # dens = get_interp_CIC( p_inside, mass, pos_x, pos_y, pos_z, nx, ny, nz, x_min, y_min, z_min, x_max, y_max, z_max, dx, dy, dz, c_pos_x, c_pos_y, c_pos_z, true )
    dens_avrg = mean(dens)
    println( "  Dens mean: $(dens_avrg)")

    outFile["$(part_type)/density"] = dens
  end
  close(inFile)
  close(outFile)
end
