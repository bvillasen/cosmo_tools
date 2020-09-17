import numpy as np

pc = 3.0857e16  #m
kpc = 1e3 * pc
Mpc = 1e6 * pc
Msun = 1.98847e30  #kg
Myear = 360 * 24 * 3600 * 1e6
Gconst = 6.6740831e-11 #m3  s-2 kg-1
Gcosmo = Gconst * ( 1./kpc) * 1./1000 * 1./1000 * Msun 




planck = np.array([0.6766,0.3111,0.0497,0.8102,0.9665])

cosmo_0 = np.array([0.6835,0.3010,0.0484,0.8098,0.9722])
cosmo_1 = np.array([0.6917,0.2905,0.0477,0.8052,0.9783])
cosmo_2 = np.array([0.7001,0.2808,0.0470,0.8020,0.9846])
cosmo_3 = np.array([0.7069,0.2730,0.0465,0.7997,0.9896])

cosmos = np.array([ planck, cosmo_0, cosmo_1, cosmo_2, cosmo_3 ])

h = cosmos[:,0]
H0 = cosmos[:,0] * 100 * 1e-3
omega_m = cosmos[:,1]
omega_b = cosmos[:,2]


rho_mean = 3*H0**2/(8*np.pi* Gcosmo) * omega_m / h**2
rho_gas_mean = 3*H0**2/(8*np.pi* Gcosmo) * omega_b 

rho_plank = rho_gas_mean[0]
rho_cosmo = rho_gas_mean[1:]

diff = (rho_cosmo - rho_plank)/rho_plank