import numpy as np
from analysis import *
import utils as ut
import utils_plots as utp

datadir = '/Volumes/T9/XSPL/Lasers/Simulations/Bluehive/OSIRIS/LasersDeck/'

# Example: External drive is assigned drive letter 'E:'
drive_letter = 'D:'
file_path_on_external_drive = 'XSPL/Lasers/Simulations/Bluehive/OSIRIS/LasersDeck/' 

# Construct the full file path
full_file_path = drive_letter + '\\' + file_path_on_external_drive
datadir = full_file_path

# List directory contents #
dirname = datadir+'Laser1D_n_ncrit_0p5_Ne_8192'
dataset = 'e3'

# Define save directory #
save_dir = drive_letter+'\\'+'XSPL/Lasers/Outputs/'
# Define laser and plasma parameters #
m_e = 9.109383e-31
elc = 1.602177e-19
eps0 = 8.854188e-12
clight = 3e8 # m/s
wavelength = 527e-9 # m
n_over_ncrit = 0.5
omega_L = 2*np.pi*clight/wavelength # rad/s    
omega_p = omega_L*np.sqrt(n_over_ncrit) # rad/s

# Compute electric field amplitude to convert from unitless to units #
amp = (elc*clight/(omega_p))/(m_e*clight**2)
Eamp = 1/amp
one_over_conv = (elc*clight)/(m_e*clight**2)
conv = (1/one_over_conv)*1e-9*1e-2 # Convert from 1/m to 1/cm
print("Conversion factor: "+str(conv))

# Compute the electric field given intensity and wavelength #
S = 3.5e16 # W/m^2
Emax = np.sqrt(2*S/(eps0*clight))
a0 = Emax*(elc*wavelength/(2*np.pi*m_e*clight**2))
print(a0)
# Print value in scientific notation #
print("Emax: "+f"{Emax:.2e}")

# Define geometry specific to the run deck #
n_wavelengths = 12 # How many wavelengths to span in the x1 direction
xsim = (n_wavelengths*wavelength)*omega_p/clight # Length of simulation in x1 direction unitless

# Example usage
directory_path = dirname 
sorted_files = ut.order_files_by_number(directory=dirname, dataset=dataset)
# print(sorted_files)

ncrit_m3 = ut.find_critical_dens(0.527)
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
time = 55
time_fs = time*t0*1e15
# make_contour(rundir=dirname,dataset='p3x1',time=time, xlim=[0,12],  xmult=clight/omega_p/wavelength, ymult=1, species='electrons', to_plot=False)
# make_contour2(rundir=dirname,dataset='p3x1',time=time, xlim=[0,12], line_out_x = 6, xmult=clight/omega_p/wavelength, ymult=1, species='electrons', to_plot=False)
# make_contour2(rundir=dirname,dataset='p3x2',time=time, xlim=[0,12], line_out_x = 6, xmult=clight/omega_p/wavelength, ymult=1, species='electrons', to_plot=False)
# utp.phasespace(rundir=dirname,dataset='x1_ene',time=time, xlim=[0,12], tmult=t0, xmult=l0, ymult=l0, species='electrons', color="Reds", to_plot=False)
# utp.phasespace(rundir=dirname,dataset='x1_m',time=time, xlim=[0,12], tmult=t0, xmult=l0, ymult=l0, species='electrons', color="copper", to_plot=False)
# utp.field(rundir=dirname,dataset='e3',time=time,xlim=[0,12], tmult=t0, xmult=l0, ymult=l0, intensitymult=Eamp, color='RdBu')
# utp.field(rundir=dirname,dataset='j3',time=time,xlim=[0,12], tmult=t0, xmult=l0, ymult=l0, intensitymult=Eamp, color='RdBu')
utp.make_contour2(rundir=dirname,dataset='p3x1',time=time, xlim=[0,12], tmult=10**15*t0, line_out_x = 6, xmult=clight/omega_p/wavelength, ymult=1, species='electrons', to_plot=True)
utp.fields(rundir=dirname,dataset=['e3', 'b2', 'j3'],time=time, xlim=[0,12], tmult=10**15*t0, xmult=l0, ymult=l0, intensitymult=[e0, b0, j0], colors=['r-', 'g-', 'b-'], to_normalize=True, to_save=True, save_dir=save_dir)
utp.fields(rundir=dirname,dataset=['j3'],time=time, xlim=[0,12], tmult=10**15*t0, xmult=l0, ymult=l0, intensitymult=[j0], colors=['b-'], 
           to_normalize=False, to_save=True, save_dir=save_dir)

# Add each plot from the field function to a frame of a movie and store this movie as a file #
if (False):
    utp.create_movie(dataset='e3', filename=dataset+".gif", num_frames=len(sorted_files), fps=5, flist=sorted_files, xmult=clight/omega_p/wavelength, ymult=Eamp, intensitymult=Emax)
