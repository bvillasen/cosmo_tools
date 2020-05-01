Gc = 4.30117902e-9 #Actually, Gc * (Msun / Mpc) in (km/s)^2
CRITICAL_DENSITY = 2.77519737e11 # 3H^2/8piG in (Msun / h) / (Mpc / h)^3
VMAX_CONST = 6.55833746e-5 #sqrt(G*(Msun/h)/(Mpc/h)) in km/s
RMAX_TO_RS = 2.1626        # Rvmax / Rs
RS_CONSTANT = 0.216216595 # ln(1+2.1626)/2.1626 - 1/(1+2.1626) */

HUBBLE_TIME_CONVERSION = 9.77813952e9 # (100 km/s/Mpc)^-1 to years */


eV_to_ergs = 1.60218e-12


#Boltazman constant
K_b = 1.38064852e-23 #m2 kg s-2 K-1

#Mass of proton
M_p = 1.6726219e-27 #kg
M_e = 9.10938356e-31 #kg
e_charge = 1.60217662e-19 # Coulombs 

c = 299792000.458 # velocity of light in m/sec
# pc = 3.086e13  #km
pc = 3.0857e16  #m
kpc = 1e3 * pc
Mpc = 1e6 * pc
Msun = 1.98847e30  #kg
Myear = 360 * 24 * 3600 * 1e6


Gconst = 6.6740831e-11 #m3  s-2 kg-1

Gcosmo = Gconst * ( 1./kpc) * 1./1000 * 1./1000 * Msun 
