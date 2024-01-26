import pickle
import numpy as np
import matplotlib.pyplot as plt
import utils as ut
osx = False
time_pic = 92
time_perseus = 180
if osx:
    datadir = '/Volumes/T9/XSPL/Lasers/Outputs/Data/'
    savedir = '/Volumes/T9/XSPL/Lasers/Outputs/Plots/'
else:
    # Example: External drive is assigned drive letter 'E:'
    drive_letter = 'D:'
    data_path_on_external_drive = 'XSPL/Lasers/Outputs/Data/' 
    plot_path_on_external_drive = 'XSPL/Lasers/Outputs/Plots/' 
    # Construct the full file path
    datadir = drive_letter + '\\' + data_path_on_external_drive
    savedir = drive_letter+'\\'+plot_path_on_external_drive

fname_pic = 'Laser1D_n_ncrit_0p5_Ne_8192_vars_e3_b2_j3_time_'+str(time_pic)
fname_perseus = 'pic_comparison_nncrit_0p5_gamma_1p666_pulse_I_3p5e16_time_'+str(time_perseus)

# Load pickle from folder
with open(datadir+fname_pic+'.p', 'rb') as f:
    data_pic = pickle.load(f)

with open(datadir+fname_perseus+'.p', 'rb') as f:
    data_perseus = pickle.load(f)

x_pic = 12*np.linspace(0, 1, data_pic['x1_m'].shape[0])
x_perseus = 12*np.linspace(0, 1, data_perseus['ne'].shape[0])
xstart_pic = np.where(x_pic < 4 )[0][-1]
xstop_pic = np.where(x_pic < 8 )[0][-1] 
xstart_perseus = np.where(x_perseus < 4 )[0][-1]
xstop_perseus = np.where(x_perseus < 8 )[0][-1] 

# Define laser and plasma parameters #
m_e = 9.109383e-31
elc = 1.602177e-19
eps0 = 8.854188e-12
clight = 3e8 # m/s
wavelength = 527e-9 # m
n_over_ncrit = 0.5
omega_L = 2*np.pi*clight/wavelength # rad/s    
omega_p = omega_L*np.sqrt(n_over_ncrit) # rad/s
mu_0 = 4*np.pi*1e-7
# Compute electric field amplitude to convert from unitless to units #
amp = (elc*clight/(omega_p))/(m_e*clight**2)
l0 = clight/omega_p
e0 = 1/amp
b0 = e0/clight
j0 = b0/(mu_0*l0)

###########
# Plot E3 #
###########
fig, ax = plt.subplots(1, 1)

# Set size of plot 
fig.set_size_inches(13.385, 6.0)

ax.plot(x_pic, e0*data_pic['e3'], label='PIC', color='xkcd:light red', linewidth=4)
ax.plot(x_perseus, data_perseus['e3'], label='Perseus', color='xkcd:red', linewidth=4, linestyle='--')
ax.set_xlabel('x [\lambda_L]')
ax.set_ylabel('E')
ax.set_xlim(0, 12)
# ax.legend()
ax2 = ax.twinx()
ax2.plot(x_pic, 0.5*data_pic['x1_m'], color='k', label='PIC', linewidth=4)
ax2.plot(x_perseus, data_perseus['ne']/ut.find_critical_dens(0.527), color='xkcd:grey', label='PIC', linewidth=4, linestyle='--')

# plt.tight_layout()
fig.savefig(savedir+fname_pic+fname_perseus+'e3_comparison.png', dpi=300)
# ax.set_title('Density at t='+str(time_pic))

###########
# Plot j3 #
###########
fig, ax = plt.subplots(1, 1)

# Set size of plot 
fig.set_size_inches(13.385, 6.0)

ax.plot(x_pic[xstart_pic:xstop_pic], j0*data_pic['j3'][xstart_pic:xstop_pic], label='PIC', color='xkcd:sky blue', linewidth=4)
ax.plot(x_perseus[xstart_perseus:xstop_perseus], data_perseus['j3'][xstart_perseus:xstop_perseus], label='Perseus', color='xkcd:bright blue', linewidth=4, linestyle='--')
ax.set_xlabel('x [\lambda_L]')
ax.set_ylabel('E')
ax.set_xlim(0, 12)
# ax.legend()
ax2 = ax.twinx()
ax2.plot(x_pic, 0.5*data_pic['x1_m'], color='k', label='PIC', linewidth=4)
ax2.plot(x_perseus, data_perseus['ne']/ut.find_critical_dens(0.527), color='xkcd:grey', label='PIC', linewidth=4, linestyle='--')

# plt.tight_layout()
fig.savefig(savedir+fname_pic+fname_perseus+'j3_comparison.png', dpi=300)
# ax.set_title('Density at t='+str(time_pic))

###########
# Plot b2 #
###########
fig, ax = plt.subplots(1, 1)

# Set size of plot 
fig.set_size_inches(13.385, 6.0)

ax.plot(x_pic, b0*data_pic['b2'], label='PIC', color='xkcd:green', linewidth=4)
ax.plot(x_perseus, data_perseus['b2'], label='Perseus', color='xkcd:blue green', linewidth=4, linestyle='--')
ax.set_xlabel('x [\lambda_L]')
ax.set_ylabel('E')
ax.set_xlim(0, 12)
# ax.legend()
ax2 = ax.twinx()
ax2.plot(x_pic, 0.5*data_pic['x1_m'], color='k', label='PIC', linewidth=4)
ax2.plot(x_perseus, data_perseus['ne']/ut.find_critical_dens(0.527), color='xkcd:grey', label='PIC', linewidth=4, linestyle='--')

# plt.tight_layout()
fig.savefig(savedir+fname_pic+fname_perseus+'b2_comparison.png', dpi=300)
# ax.set_title('Density at t='+str(time_pic))

plt.show()
