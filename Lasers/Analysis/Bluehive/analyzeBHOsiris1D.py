import numpy as np
from analysis import *
import utils as ut
import utils_plots as utp

osx = False

if osx: 
    datadir = '/Volumes/T9/XSPL/Lasers/Simulations/Bluehive/OSIRIS/LasersDeck/'
    # Define save directory #
    save_dir = '/Volumes/T9/XSPL/Lasers/Outputs/'


else:

    # Example: External drive is assigned drive letter 'E:'
    drive_letter = 'D:'
    file_path_on_external_drive = 'XSPL/Lasers/Simulations/Bluehive/OSIRIS/LasersDeck/' 

    # Construct the full file path
    full_file_path = drive_letter + '\\' + file_path_on_external_drive
    datadir = full_file_path

    # Define save directory #
    save_dir = drive_letter+'\\'+'XSPL/Lasers/Outputs/'

# Define specific run now #
dirname = datadir+'Laser1D_n_ncrit_0p5_Ne_8192_S_x10'
dirname = datadir+'Laser1D_n_ncrit_0p5_Ne_8192_S_x10_long'
dirname = datadir+'Laser1D_n_ncrit_0p5_Ne_4096_S_x10_long'
# dirname = datadir+'Laser1D_n_ncrit_0p5_Ne_128_S_x10_long'
dirname = datadir+'Laser1D_n_ncrit_0p5_Ne_128_S_x10_long_wcoll_n0'
# dirname = datadir+'Laser1D_n_ncrit_0p5_Ne_128_S_x10_long_wocoll_n0'
# dirname = datadir+'Laser1D_n_ncrit_0p5_Ne_512_S_x10_long_wcoll_n0'
dirname = datadir+'Laser1D_n_ncrit_0p5_Ne_2048_S_x10_long_wcoll_n0'
dirname = datadir+'Laser1D_n_ncrit_0p5_Ne_512_S_x10_long_wcoll_n0_C14'

# dirname = datadir+'Laser1D_n_ncrit_0p5_Ne_8192'
S = 3.5e17 # W/m^2

dataset = 'e3'

# Define laser and plasma parameters #
m_e = 9.109383e-31
elc = 1.602177e-19
eps0 = 8.854188e-12
clight = 3e8 # m/s
wavelength = 527e-9 # m
wavelength_um = wavelength*1e6 # um
n_over_ncrit = 0.5
omega_L = 2*np.pi*clight/wavelength # rad/s    
omega_p = omega_L*np.sqrt(n_over_ncrit) # rad/s
# Print value in scientific notation #
print("Omega_p: "+f"{omega_p:.5e}")

# Compute electric field amplitude to convert from unitless to units #
amp = (elc*clight/(omega_p))/(m_e*clight**2)
Eamp = 1/amp
one_over_conv = (elc*clight)/(m_e*clight**2)
conv = (1/one_over_conv)*1e-9*1e-2 # Convert from 1/m to 1/cm
print("Conversion factor: "+str(conv))

# Compute the electric field given intensity and wavelength #
S_cm = S/(100**2) # W/cm^2
Emax = np.sqrt(2*S/(eps0*clight))
a0 = Emax*(elc*wavelength/(2*np.pi*m_e*clight**2))

a0_theory = 0.85*np.sqrt((wavelength_um)**2*S_cm/1e18)
print("a0: "+str(a0)+", a0 theory: "+str(a0_theory))

v_o_c = 0.0084
gamma = 1/np.sqrt(1-v_o_c**2)
print("Gamma: "+str(gamma))

# Print value in scientific notation #
print("Emax: "+f"{Emax:.2e}")

# Define geometry specific to the run deck #
n_wavelengths = 12 # How many wavelengths to span in the x1 direction
xsim = (n_wavelengths*wavelength)*omega_p/clight # Length of simulation in x1 direction unitless

ncrit_m3 = ut.find_critical_dens(0.527)
print("Critical density: "+str(ncrit_m3)+" m^-3")
ncrit_cm3 = ncrit_m3*(1e-6)
print("Half Critical density: "+str(ncrit_m3/2)+" m^-3")

# Compute units for length, time, and electric field #
mu_0 = 4*np.pi*1e-7
l0 = clight/omega_p/wavelength # Interms of wavelength
l0_m = clight/omega_p
t0 = 1/omega_p
e0 = Eamp
b0 = e0/clight
j0 = b0/l0_m/mu_0
time = 1480
time = 14790
time = 14780
time_fs = time*t0*1e15

xoutE, youtE = utp.field(rundir=dirname,dataset='e3',time=time,xlim=[0,12], tmult=t0, xmult=l0, ymult=l0, intensitymult=e0, color='RdBu', to_plot=False)
xoutB, youtB = utp.field(rundir=dirname,dataset='b2',time=time,xlim=[0,12], tmult=t0, xmult=l0, ymult=l0, intensitymult=b0, color='RdBu', to_plot=False)
# utp.make_contour2(rundir=dirname,dataset='gammax1',time=time, xlim=[0,12], tmult=10**15*t0, line_out_x = 6, xmult=clight/omega_p/wavelength, ymult=1, species='electrons', to_plot=True, to_save=True, save_dir=save_dir)
utp.make_contour2(rundir=dirname,dataset='p1x1',time=time, xlim=[0,12], tmult=10**15*t0, line_out_x = 5.75, to_fit=False, xmult=clight/omega_p/wavelength, ymult=1, species='ions', to_plot=False, to_save=True, save_dir=save_dir)
utp.make_contour2(rundir=dirname,dataset='p1x1',time=time, xlim=[0,12], tmult=10**15*t0, line_out_x = 6.1, to_fit=False, xmult=clight/omega_p/wavelength, ymult=1, species='ions', to_plot=False, to_save=True, save_dir=save_dir)
utp.make_contour2(rundir=dirname,dataset='p3x1',time=time, xlim=[0,12], tmult=10**15*t0, line_out_x = 5.75, xmult=clight/omega_p/wavelength, ymult=1, species='electrons', to_plot=False, to_save=True, save_dir=save_dir)
utp.make_contour2(rundir=dirname,dataset='p3x1',time=time, xlim=[0,12], tmult=10**15*t0, line_out_x = 6.1, xmult=clight/omega_p/wavelength, ymult=1, species='electrons', to_plot=False, to_save=True, save_dir=save_dir)
utp.fields(rundir=dirname,dataset=['e3', 'b2', 'j3'],time=time, xlim=[0,12], tmult=10**15*t0, xmult=l0, ymult=l0, intensitymult=[e0, b0, j0], colors=['r-', 'g-', 'b-'], to_normalize=True, to_plot=False, to_save=True, save_dir=save_dir)
utp.fields(rundir=dirname,dataset=['e3'],time=time, xlim=[0,12], tmult=10**15*t0, xmult=l0, ymult=l0, intensitymult=[e0], colors=['b-'], 
           to_normalize=False, to_save=True, save_dir=save_dir)

figure, ax = plt.subplots(1, 1)
ax.plot(xoutE, youtE*youtB/(4e-7*np.pi), label='S', color='r')
ax.set_xlabel('x [\lambda_L]')
ax.set_ylabel('E, B')
ax.legend()
plt.show()

# Add each plot from the field function to a frame of a movie and store this movie as a file #
if (False):
    directory_path = dirname 
    sorted_files = ut.order_files_by_number(directory=dirname, dataset=dataset)
    utp.create_movie(dataset='e3', filename=dataset+".gif", num_frames=len(sorted_files), fps=5, flist=sorted_files, xmult=clight/omega_p/wavelength, ymult=Eamp, intensitymult=1e10)
