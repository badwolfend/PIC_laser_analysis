import pickle
import numpy as np
import matplotlib.pyplot as plt
import utils as ut
osx = False

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

fname_pic = ['Laser1D_n_ncrit_0p5_Ne_512_S_x10_long_wcoll_n0_C14_temperature_times']
fname_perseus = ['pic_comparison_nncrit_0p5_gamma_1p1_pulse_I_3p5e17_0p1Ln_t_long_temperature_times']

fname_pic = ['Laser1D_n_ncrit_0p5_Ne_512_S_x10_long_wcoll_n0_He3_vion_temperature_times']
fname_perseus = ['pic_comparison_nncrit_0p5_gamma_1p666_pulse_I_3p5e17_0p1Ln_t_long_He3_cond_visc_temperature_times']

fpicstr =''

fig0, ax0 = plt.subplots(1, 1)

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

    # Load the temperature values in time  
    time_pic = data_pic['time']
    te0_pic = data_pic['line_out_x_0']
    te1_pic = data_pic['line_out_x_1']

    time_perseus = data_perseus['time']
    te0_perseus = data_perseus['line_out_x_0']
    te1_perseus = data_perseus['line_out_x_1']

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
    t0 = 1/omega_p

    # Get the theory temperature
    tcode_array, Te0 = ut.get_theory_temperature(Z=1, Te0=1, ni=1.980347749E27, ne=1.980347749E27, log_A=0.1, SLas=3.5e17, tstop=time_perseus[-1])

    ###########
    # Plot E3 #
    ###########
    # fig0, ax0 = plt.subplots(1, 1)

    # Set size of plot 
    fig0.set_size_inches(13.385, 6.69)

    # ax0.plot(np.array(time_pic), te0_pic, label='PIC 0', linewidth=4, alpha=0.75, color='xkcd:sky blue')
    # ax0.plot(time_perseus, te0_perseus, label='Perseus 0', linewidth=4, alpha=0.75, color='xkcd:blue')
    ax0.plot(np.array(time_pic), te1_pic, label='PIC 1', linewidth=4, alpha=0.75, color='xkcd:rose')
    ax0.plot(time_perseus, te1_perseus, label='Perseus 1', linewidth=4, alpha=0.75, color='xkcd:blue')
    ax0.plot(tcode_array, Te0, label='Perseus 1', linewidth=4, linestyle=(0, (3, 1, 1, 1, 1, 1)), color='xkcd:black')

    # ax0.plot(np.array(time_pic), te1_pic, label='PIC 0', linewidth=4, alpha=0.75, color='xkcd:black')
    # ax0.plot(time_perseus, te1_perseus, label='Perseus 0', linewidth=4, alpha=0.75, color='xkcd:bright blue')
    ax0.set_xlabel('Time [s]')
    ax0.set_ylabel('Te')
    # ax0.set_xlim(0, 12)
    ax0.set_ylim(-10, 55)
    ax0.legend()
    fig0.tight_layout()
    fig0.savefig(savedir+fpicstr+'_temperature_comparison.png', dpi=600)
    print('Saved: '+savedir+fpicstr+'_temperature_comparison.png')
plt.show()
