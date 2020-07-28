module CIC_functions

# importall Base

export get_particles_outside_CIC, get_interp_CIC


using HDF5



function get_particles_outside_CIC( p_inside, p_pos_x, p_pos_y, p_pos_z,
             x_min, x_max, y_min, y_max, z_min, z_max, dx, dy, dz )
  nParticles = size( p_pos_x )[1]
  for i in 1:nParticles
    x, y, z = p_pos_x[i], p_pos_y[i], p_pos_z[i]
    if  ( x<(x_min) || x>(x_max) )
      p_inside[i] = false
      continue
    end
    if  ( y<(y_min) || y>(y_max) )
      p_inside[i] = false
      continue
    end
    if  ( z<(z_min) || z>(z_max) )
      p_inside[i] = false
      continue
    end
  end
end

function get_indxs_CIC( dx, dy, dz, p_pos_x, p_pos_y, p_pos_z )
  idxs_x = floor( Int, (p_pos_x - 0.5*dx)/dx ) + 1
  idxs_y = floor( Int, (p_pos_y - 0.5*dy)/dy ) + 1
  idxs_z = floor( Int, (p_pos_z - 0.5*dz)/dz ) + 1
  return [ idxs_x, idxs_y, idxs_z ]
end

function get_indx_CIC(  dx, dy, dz, p_pos_x, p_pos_y, p_pos_z )
  idx_x = floor( Int, (p_pos_x - 0.5*dx)/dx ) + 1
  idx_y = floor( Int, (p_pos_y - 0.5*dy)/dy ) + 1
  idx_z = floor( Int, (p_pos_z - 0.5*dz)/dz ) + 1
  return [ idx_x, idx_y, idx_z ]
end

function get_interp_CIC( p_inside, p_field, p_dens, p_mass, p_pos_x, p_pos_y, p_pos_z, nCells_x, nCells_y, nCells_z,
  x_min, y_min, z_min, x_max, y_max, z_max, dx, dy, dz, density, periodic )
  nParticles = size( p_pos_x )[1]
  # idxs_x, idxs_y, idxs_z = get_indxs_CIC( dx, dy, dz, p_pos_x, p_pos_y, p_pos_z )
  #Interpolation array array
  interp = zeros( nCells_z+2, nCells_y+2, nCells_x+2 )
  for i in 1:nParticles
    if !p_inside[i]
      continue
    end
    m = p_field[i]
    mass = p_mass[i]
    dens = p_dens[i]
    pos_x = p_pos_x[i]
    pos_y = p_pos_y[i]
    pos_z = p_pos_z[i]
    # if !density
    #   m *= mass / dens
    # end
    idx_x, idx_y, idx_z = get_indx_CIC(  dx, dy, dz, pos_x, pos_y, pos_z )
    delta_x = 1 - (pos_x - (x_min + 0.5*dx + (idx_x-1)*dx) )/dx
    delta_y = 1 - (pos_y - (y_min + 0.5*dy + (idx_y-1)*dy) )/dy
    delta_z = 1 - (pos_z - (z_min + 0.5*dz + (idx_z-1)*dz) )/dz
    idx_x += 1
    idx_y += 1
    idx_z += 1
    # println("$(idx_x)    $(idx_y)    $(idx_z)    ")
    # println("$(delta_x)    $(delta_y)    $(delta_z)    ")
    interp[idx_z, idx_y, idx_x] += m * delta_x * delta_y * delta_z
    interp[idx_z, idx_y, idx_x+1] += m * (1-delta_x) * delta_y * delta_z
    interp[idx_z, idx_y+1, idx_x] += m * delta_x * (1-delta_y) * delta_z
    interp[idx_z+1, idx_y, idx_x] += m * delta_x * delta_y * (1-delta_z)
    interp[idx_z+1, idx_y+1, idx_x] += m * delta_x * (1-delta_y) * (1-delta_z)
    interp[idx_z, idx_y+1, idx_x+1] += m * (1-delta_x) * (1-delta_y) * delta_z
    interp[idx_z+1, idx_y, idx_x+1] += m * (1-delta_x) * delta_y * (1-delta_z)
    interp[idx_z+1, idx_y+1, idx_x+1] += m * (1-delta_x) * (1-delta_y) * (1-delta_z)
  end
  if periodic
    interp[2,:,:] += interp[end,:,:]
    interp[end-1,:,:] += interp[1,:,:]
    interp[:,2,:] += interp[:,end,:]
    interp[:,end-1,:] += interp[:,1,:]
    interp[:,:,2] += interp[:,:,end]
    interp[:,:,end-1] += interp[:,:,1]
  end
  if density
    interp /= dx*dy*dz
  end
  return interp[2:end-1,2:end-1,2:end-1]
end

end
