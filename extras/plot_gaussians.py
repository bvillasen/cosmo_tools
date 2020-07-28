import numpy as np
import matplotlib.pyplot as plt

c = 5
b = 15

c_1 = 15
w = 5

v_min = -30
v_max = 40
nPoints = 1000
v_range = np.linspace( v_min, v_max, nPoints ) 

gaussian_func = lambda x: 1./( np.sqrt(np.pi) * b ) * np.exp( -( ( x - c)/b )**2 )

gaussian = gaussian_func( v_range )

nrows = 1
ncols = 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12*ncols,8*nrows))
fs = 13


ax.plot( v_range, gaussian )
text = r'$n_{\mathrm{HI}, i}  \,\,\,\,\, b_i \,\,\,\,\,  v_i$'
ax.text(0.5, 0.90, text, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=15)
ymax = 0.043
ax.axvline( c, ymax=gaussian_func(c)/ymax, ls='--', c='C2')
ax.set_ylim( 0, ymax)

x_c = c + c_1
x_l = c + c_1 - w
x_r = c + c_1 + w
# v_range = np.linspace( x_l, x_r, nPoints )
# ax.axvline( x_c, ymax=gaussian_func(x_c)/ymax, ls='--', c='C3')
ax.axvline( x_l, ymax=gaussian_func(x_l)/ymax, ls='--', c='C4')
ax.axvline( x_r, ymax=gaussian_func(x_r)/ymax, ls='--', c='C4')

# text = r'$v_j$'
# ax.text(0.7, 0.35, text, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=15)
# 
text = r'$v_{j-1/2}$'
ax.text(0.655, 0.59, text, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=15)

text = r'$v_{j+1/2}$'
ax.text(0.78, 0.18, text, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=15)
# ax.fill_between( v_range, gaussian_func(v_range), color='C4', alpha=0.4)


# 
# x_c = c + c_1
# x_l = c - c_1 - w
# x_r = c + c_1 + w
# v_range = np.linspace( x_l, x_r, nPoints )
# ax.fill_between( v_range, gaussian_func(v_range), color='C2', alpha=0.4)

# 
# x_c = c + c_1
# x_l = c - c_1 + w
# x_r = c + c_1 - w
# v_range = np.linspace( x_l, x_r, nPoints )
# ax.fill_between( v_range, gaussian_func(v_range), color='C2', alpha=0.4)
# 

x_c = c + c_1
x_l = c + c_1 + w
x_r = c + c_1 - w
v_range = np.linspace( x_l, x_r, nPoints )
ax.fill_between( v_range, gaussian_func(v_range), color='C2', alpha=0.4)

x_c = c + c_1
x_l = c - c_1 - w
x_r = c - c_1 + w
v_range = np.linspace( x_l, x_r, nPoints )
ax.fill_between( v_range, gaussian_func(v_range), color='C2', alpha=0.4)



ax.set_xlabel('Velocity', fontsize=fs)

fig_name = '/home/bruno/Desktop/gaussian_4.png'
fig.savefig( fig_name,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
print(( 'Saved Image: ', fig_name ))
