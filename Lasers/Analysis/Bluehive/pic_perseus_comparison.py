import pickle
import numpy as np
import matplotlib.pyplot as plt
import utils as ut
osx = False

# time_pic = 2471
# time_perseus = 494
time_pic = 2465
time_perseus = 494
time_pic = 1787
time_perseus = 492
time_pic = 1786
time_pic = 536
time_perseus = 495
time_perseus = 17

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

fname_pic = ['Laser1D_n_ncrit_0p5_Ne_8192_vars_e3_b2_j3_time_'+str(time_pic), 'Laser1D_n_ncrit_0p5_Ne_8192_S_x10_vars_e3_b2_j3_time_'+str(time_pic)]
fname_perseus = ['pic_comparison_nncrit_0p5_gamma_1p666_pulse_I_3p5e16_time_'+str(time_perseus), 'pic_comparison_nncrit_0p5_gamma_1p666_pulse_I_3p5e17_time_'+str(time_perseus)]

fname_pic = ['Laser1D_n_ncrit_0p5_Ne_8192_S_x10_long_vars_e3_b2_j3_time_'+str(time_pic)]
fname_perseus = ['pic_comparison_nncrit_0p5_gamma_1p666_pulse_I_3p5e17_0p1Ln_t_long_time_'+str(time_perseus)]

fname_pic = ['Laser1D_n_ncrit_0p5_Ne_4096_S_x10_long_vars_e3_b2_j3_time_'+str(time_pic), 'Laser1D_n_ncrit_0p5_Ne_8192_S_x10_long_vars_e3_b2_j3_time_'+str(time_pic)]
fname_perseus = ['pic_comparison_nncrit_0p5_gamma_1p666_pulse_I_3p5e17_0p1Ln_t_long_time_'+str(time_perseus),'pic_comparison_nncrit_0p5_gamma_1p666_pulse_I_3p5e17_0p1Ln_t_long_time_'+str(time_perseus)]

fname_pic = ['Laser1D_n_ncrit_0p5_Ne_128_S_x10_long_vars_e3_b2_j3_time_'+str(time_pic), 'Laser1D_n_ncrit_0p5_Ne_8192_S_x10_long_vars_e3_b2_j3_time_'+str(time_pic)]
fname_perseus = ['pic_comparison_nncrit_0p5_gamma_1p666_pulse_I_3p5e17_0p1Ln_t_long_time_'+str(time_perseus),'pic_comparison_nncrit_0p5_gamma_1p666_pulse_I_3p5e17_0p1Ln_t_long_time_'+str(time_perseus)]

fname_pic = ['Laser1D_n_ncrit_0p5_Ne_512_S_x10_long_wcoll_n0_C14_vars_e3_b2_j3_time_'+str(time_pic), 'Laser1D_n_ncrit_0p5_Ne_512_S_x10_long_wcoll_n0_C14_vars_e3_b2_j3_time_'+str(time_pic)]
fname_perseus = ['pic_comparison_nncrit_0p5_gamma_1p1_pulse_I_3p5e17_0p1Ln_t_long_time_'+str(time_perseus),'pic_comparison_nncrit_0p5_gamma_1p1_pulse_I_3p5e17_0p1Ln_t_long_time_'+str(time_perseus)]

fname_pic = ['Laser1D_n_ncrit_0p5_Ne_512_S_x10_long_wcoll_n0_He3_vars_e3_b2_j3_time_'+str(time_pic)]
fname_perseus = ['pic_comparison_nncrit_0p5_gamma_1p666_pulse_I_3p5e17_0p1Ln_t_long_He3_time_'+str(time_perseus)]

fname_pic = ['Laser1D_n_ncrit_0p5_Ne_512_S_x10_long_wcoll_n0_He3_vars_e3_b2_j3_time_'+str(time_pic)]
fname_perseus = ['pic_comparison_nncrit_0p5_gamma_1p666_pulse_I_3p5e17_0p1Ln_t_long_He3_forpic_time_'+str(time_perseus)]

fname_pic = ['Laser1D_n_ncrit_0p5_Ne_512_S_x10_long_wcoll_n0_He3_vion_vars_e3_b2_j3_time_'+str(time_pic)]
fname_perseus = ['pic_comparison_nncrit_0p5_gamma_1p666_pulse_I_3p5e17_0p1Ln_t_long_He3_cond_visc_forpic_time_'+str(time_perseus)]

fname_pic = ['Laser1D_n_ncrit_0p5_Ne_512_S_x10_long_wcoll_n0_He3_2f_Cl_1_vars_e3_b2_j3_time_'+str(time_pic)]
fname_perseus = ['pic_comparison_nncrit_0p5_gamma_1p666_pulse_I_3p5e17_1p0Ln_t_long_He3_cond_visc_forpic_time_'+str(time_perseus)]

pic_colors = {0:{'e3':'xkcd:light red', 'j3':'xkcd:sky blue', 'b2':'xkcd:green', 'ne':'xkcd:kelly green', 'ne0':'xkcd:black'}, 1:{'e3':'xkcd:light red', 'j3':'xkcd:cyan', 'b2':'xkcd:green','ne':'k'}}
pic_alpha = [0.65, 0.90]
perseus_colors = {0:{'e3':'xkcd:red', 'j3':'xkcd:bright blue', 'b2':'xkcd:blue green', 'ne':'xkcd:burnt umber'}, 1:{'e3':'xkcd:red', 'j3':'xkcd:bright blue', 'b2':'xkcd:blue green', 'ne':'xkcd:burnt umber'}}
perseus_alpha = [1, 1]
fpicstr =''

fig0, ax0 = plt.subplots(1, 1)
fig1, ax1 = plt.subplots(1, 1)
fig2, ax2 = plt.subplots(1, 1)

for fi, fpic in enumerate(fname_pic):
    if fi == 0:
        fpicstr = fpicstr+fpic+'_'+fname_perseus[fi]
    else:
        fpicstr = fpicstr+"_"+fpic

    # Load pickle from folder
    with open(datadir+fname_pic[fi]+'.p', 'rb') as f:
        data_pic = pickle.load(f)

    with open(datadir+fname_perseus[fi]+'.p', 'rb') as f:
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
    maxEperseus = np.max(data_perseus['e3'])  
    maxJperseus = np.max(data_perseus['j3'])
    maxBperseus = np.max(data_perseus['b2'])  
    ###########
    # Plot E3 #
    ###########
    # fig0, ax0 = plt.subplots(1, 1)

    # Set size of plot 
    fig0.set_size_inches(13.385, 6.0)

    ax0.plot(x_pic, e0*data_pic['e3']/maxEperseus, label='PIC', color=pic_colors[fi]['e3'], linewidth=4, alpha=pic_alpha[fi])
    ax0.plot(x_perseus, data_perseus['e3']/maxEperseus, label='Perseus', color=perseus_colors[fi]['e3'], linewidth=4, linestyle='--', alpha=perseus_alpha[fi])
    ax0.set_xlabel('x [\lambda_L]')
    ax0.set_ylabel('E')
    ax0.set_xlim(0, 12)
    ax0.set_ylim(-1.1, 1.1)
    # ax.legend()
    ax02 = ax0.twinx()
    ax02.plot(x_pic, data_pic['x1_m'], color=pic_colors[fi]['ne'], label='PIC', linewidth=4, alpha=pic_alpha[fi])
    ax02.plot(x_perseus, data_perseus['ne']/ut.find_critical_dens(0.527), color=perseus_colors[fi]['ne'], label='PIC', linewidth=4, linestyle='--', alpha=perseus_alpha[fi])
    ax02.set_ylim(-0.01, 0.6)

    # plt.tight_layout()
    # fig.savefig(savedir+fname_pic+fname_perseus+'e3_comparison.png', dpi=300)
    # ax.set_title('Density at t='+str(time_pic))

    ###########
    # Plot j3 #
    ###########
    # fig1, ax1 = plt.subplots(1, 1)

    # Set size of plot 
    fig1.set_size_inches(13.385, 6.0)
    maxJpic = np.max(j0*data_pic['j3'][xstart_pic:xstop_pic])
    ax1.plot(x_pic[xstart_pic:xstop_pic], j0*data_pic['j3'][xstart_pic:xstop_pic]/maxJperseus, label='PIC', color=pic_colors[fi]['j3'], linewidth=4, alpha=0.8)
    ax1.plot(x_perseus[xstart_perseus:xstop_perseus],data_perseus['j3'][xstart_perseus:xstop_perseus]/maxJperseus, label='Perseus', color=perseus_colors[fi]['j3'], linewidth=4, linestyle='--', alpha=perseus_alpha[fi])
    ax1.set_xlabel('x [\lambda_L]')
    ax1.set_ylabel('E')
    ax1.set_xlim(0, 12)
    # ax.legend()
    ax12 = ax1.twinx()
    ax12.plot(x_pic, data_pic['x1_m'], color=pic_colors[fi]['ne'], label='PIC', linewidth=4, alpha=pic_alpha[fi])
    ax12.plot(x_pic, data_pic['x1_m_density_t0'], color=pic_colors[fi]['ne0'], label='PIC T0', linewidth=4, alpha=1)
    ax12.plot(x_perseus, data_perseus['ne']/ut.find_critical_dens(0.527), color=perseus_colors[fi]['ne'], label='PIC', linewidth=4, linestyle='--', alpha=perseus_alpha[fi])
    ax12.set_ylim(-0.01, 0.6)

    # plt.tight_layout()
    # fig.savefig(savedir+fname_pic+fname_perseus+'j3_comparison.png', dpi=300)
    # ax12.set_title('Density at t='+str(time_pic))

    ###########
    # Plot b2 #
    ###########
    # fig2, ax2 = plt.subplots(1, 1)

    # Set size of plot 
    fig2.set_size_inches(13.385, 6.0)

    ax2.plot(x_pic, b0*data_pic['b2']/maxBperseus, label='PIC', color=pic_colors[fi]['b2'], linewidth=4, alpha=pic_alpha[fi])
    ax2.plot(x_perseus, data_perseus['b2']/maxBperseus, label='Perseus', color=perseus_colors[fi]['b2'], linewidth=4, linestyle='--', alpha=perseus_alpha[fi])
    ax2.set_xlabel('x [\lambda_L]')
    ax2.set_ylabel('E')
    ax2.set_xlim(0, 12)

    # ax.legend()
    ax22 = ax2.twinx()
    ax22.plot(x_pic, data_pic['x1_m'], color=pic_colors[fi]['ne'], label='PIC', linewidth=4, alpha=pic_alpha[fi])
    ax22.plot(x_perseus, data_perseus['ne']/ut.find_critical_dens(0.527), color=perseus_colors[fi]['ne'], label='PIC', linewidth=4, linestyle='--', alpha=perseus_alpha[fi])
    ax22.set_ylim(-0.01, 0.6)

    # plt.tight_layout()
    # fig.savefig(savedir+fname_pic+fname_perseus+'b2_comparison.png', dpi=300)
    # ax.set_title('Density at t='+str(time_pic))


fig0.savefig(savedir+fpicstr+'e3_comparison.png', dpi=300)
fig1.savefig(savedir+fpicstr+'j3_comparison.png', dpi=300)
fig2.savefig(savedir+fpicstr+'b2_comparison.png', dpi=300)
# fig0.savefig(savedir+fname_pic+fname_perseus+'e3_comparison.png', dpi=300)
# fig1.savefig(savedir+fname_pic+fname_perseus+'j3_comparison.png', dpi=300)
# fig2.savefig(savedir+fname_pic+fname_perseus+'b2_comparison.png', dpi=300)

plt.show()
