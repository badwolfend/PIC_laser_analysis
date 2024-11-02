import pickle
import numpy as np
import matplotlib.pyplot as plt
import utils as ut
osx = False

# time_pic = 2471
# time_perseus = 494
time_pic = 2465
time_perseus = 494
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

fname_pic = ['Laser1D_n_ncrit_0p5_Ne_128_S_x10_long_vars_e3_b2_j3_time_'+str(time_pic), 'Laser1D_n_ncrit_0p5_Ne_8192_S_x10_long_vars_e3_b2_j3_time_'+str(time_pic)]

pic_colors = {0:{'e3':'xkcd:light red', 'j3':'xkcd:sky blue', 'b2':'xkcd:green', 'ne':'k'}, 1:{'e3':'xkcd:light red', 'j3':'xkcd:cyan', 'b2':'xkcd:green','ne':'k'}}
pic_alpha = [0.65, 0.90]
fpicstr =''

fig0, ax0 = plt.subplots(1, 1)
fig1, ax1 = plt.subplots(1, 1)
fig2, ax2 = plt.subplots(1, 1)

# Setup a list of arrays to store the data
data = []
data2 = []
for fi, fpic in enumerate(fname_pic):
    if fi == 0:
        fpicstr = fpicstr+fpic
    else:
        fpicstr = fpicstr+"_"+fpic

    # Load pickle from folder
    with open(datadir+fname_pic[fi]+'.p', 'rb') as f:
        data_pic = pickle.load(f)

    x_pic = 12*np.linspace(0, 1, data_pic['x1_m'].shape[0])
    xstart_pic = np.where(x_pic < 4 )[0][-1]
    xstop_pic = np.where(x_pic < 8 )[0][-1] 
    xstart_pic_interior = np.where(x_pic < 5 )[0][-1]
    xstop_pic_interior = np.where(x_pic < 7 )[0][-1] 
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
    maxEpic = np.max(data_pic['e3'])  
    maxJpic = np.max(data_pic['j3'])
    maxBpic = np.max(data_pic['b2'])  

    ###########
    # Plot E3 #
    ###########
    # fig0, ax0 = plt.subplots(1, 1)

    # Set size of plot 
    fig0.set_size_inches(13.385, 6.0)

    ax0.plot(x_pic, data_pic['e3'], label='PIC', color=pic_colors[fi]['e3'], linewidth=4, alpha=pic_alpha[fi])
    ax0.set_xlabel('x [\lambda_L]')
    ax0.set_ylabel('E')
    ax0.set_xlim(0, 12)
    ax0.set_ylim(-1.1, 1.1)
    # ax.legend()
    ax02 = ax0.twinx()
    ax02.plot(x_pic, 0.5*data_pic['x1_m'], color=pic_colors[fi]['ne'], label='PIC', linewidth=4, alpha=pic_alpha[fi])
    ax02.set_ylim(-0.01, 0.6)

    ###########
    # Plot j3 #
    ###########
    # fig1, ax1 = plt.subplots(1, 1)

    # Set size of plot 
    fig1.set_size_inches(13.385, 6.0)
    ax1.plot(x_pic[xstart_pic:xstop_pic], data_pic['j3'][xstart_pic:xstop_pic], label='PIC', color=pic_colors[fi]['j3'], linewidth=4, alpha=0.8)
    ax1.set_xlabel('x [\lambda_L]')
    ax1.set_ylabel('E')
    ax1.set_xlim(0, 12)
    # ax.legend()
    ax12 = ax1.twinx()
    ax12.plot(x_pic, 0.5*data_pic['x1_m'], color=pic_colors[fi]['ne'], label='PIC', linewidth=4, alpha=pic_alpha[fi])
    ax12.set_ylim(-0.01, 0.6)

    ###########
    # Plot b2 #
    ###########
    # fig2, ax2 = plt.subplots(1, 1)

    # Set size of plot 
    fig2.set_size_inches(13.385, 6.0)

    ax2.plot(x_pic, data_pic['b2']/maxBpic, label='PIC', color=pic_colors[fi]['b2'], linewidth=4, alpha=pic_alpha[fi])
    ax2.set_xlabel('x [\lambda_L]')
    ax2.set_ylabel('E')
    ax2.set_xlim(0, 12)

    # ax.legend()
    ax22 = ax2.twinx()
    ax22.plot(x_pic, 0.5*data_pic['x1_m'], color=pic_colors[fi]['ne'], label='PIC', linewidth=4, alpha=pic_alpha[fi])
    ax22.set_ylim(-0.01, 0.6)

    data.append(0.5*data_pic['x1_m'][xstart_pic_interior:xstop_pic_interior])
    data2.append(data_pic['j3'][xstart_pic:xstop_pic])
    # data2.append(data_pic['j3'][xstart_pic_interior:xstop_pic_interior])

fig0.savefig(savedir+fpicstr+'_noise_e3_comparison.png', dpi=300)
fig1.savefig(savedir+fpicstr+'_noise_j3_comparison.png', dpi=300)
fig2.savefig(savedir+fpicstr+'_noise_b2_comparison.png', dpi=300)

plt.show()

## Now compute the noise in the data ##
# Ponderomotive density fluctuations
dn = np.max(data[1]) - np.min(data[1])
dj = np.max(data2[1]) - np.min(data2[1])

# Assuming data[1] is the true density (little to no noise)
ndiff = data[0] - data[1]
relative_noise = ndiff / data[1]  # Proportional noise relative to true density
jdiff = data2[0] - data2[1]
relative_noise_j = jdiff / data2[1]  # Proportional noise relative to true density

# Compute the standard deviation of the relative noise
proportional_noise_std = np.std(relative_noise)
proportional_noise_std_j = np.std(relative_noise_j)
print("Proportional Noise Standard Deviation: ", proportional_noise_std)
print("Proportional Noise Standard Deviation j: ", proportional_noise_std_j)

# For a density of N=0.5*ncrit, the from the above value we can compute the noise in the density signal
proportional_noise_std_n0p5 = proportional_noise_std
print("Proportional Noise Standard Deviation for N=0.5*ncrit: ", proportional_noise_std_n0p5)
print("Ponderomotive Density Fluctuations: ", dn/0.5/2)
mult=1.0

# Visualization
plt.figure(figsize=(12, 6))

# Plot relative noise (jdiff / true density)
# plt.plot(mult*relative_noise, label="Relative Noise (ndiff / true density)", color="purple")
plt.plot(data2[0], label="Relative Noise (ndiff / true density)", color="purple")
plt.plot(data2[1], label="Relative Noise (ndiff / true density)", color="black")
# plt.axhline(mult*proportional_noise_std, color="red", linestyle="--", label=f"Proportional Noise Std = {mult*proportional_noise_std:.3f}")
# plt.axhline(-mult*proportional_noise_std, color="red", linestyle="--")
plt.title("Proportional Noise in Noisy Density Signal Relative to True Density")
plt.xlabel("Index")
plt.ylabel("Relative Noise")
plt.legend()

plt.show()