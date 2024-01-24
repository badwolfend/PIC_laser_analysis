import matplotlib.pyplot as plt
import h5py
import numpy as np
import os
from analysis import *
from matplotlib.animation import PillowWriter, FuncAnimation
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
import moviepy.editor as mpy
import re
from osh5def import H5Data, PartData, fn_rule, DataAxis, OSUnits
from scipy.optimize import curve_fit

def read_h5(filename, path=None, axis_name="AXIS/AXIS"):
    """
    HDF reader for Osiris/Visxd compatible HDF files... This will slurp in the data
    and the attributes that describe the data (e.g. title, units, scale).

    Usage:
            diag_data = read_hdf('e1-000006.h5')      # diag_data is a subclass of numpy.ndarray with extra attributes

            print(diag_data)                          # print the meta data
            print(diag_data.view(numpy.ndarray))      # print the raw data
            print(diag_data.shape)                    # prints the dimension of the raw data
            print(diag_data.run_attrs['TIME'])        # prints the simulation time associated with the hdf5 file
            diag_data.data_attrs['UNITS']             # print units of the dataset points
            list(diag_data.data_attrs)                # lists all attributes related to the data array
            list(diag_data.run_attrs)                 # lists all attributes related to the run
            print(diag_data.axes[0].attrs['UNITS'])   # prints units of X-axis
            list(diag_data.axes[0].attrs)             # lists all variables of the X-axis

            diag_data[slice(3)]
                print(rw.view(np.ndarray))

    We will convert all byte strings stored in the h5 file to strings which are easier to deal with when writing codes
    see also write_h5() function in this file

    """
    fname = filename if not path else path + '/' + filename
    data_file = h5py.File(fname, 'r')

    n_data = scan_hdf5_file_for_main_data_array(data_file)

    timestamp, name, run_attrs, data_attrs, axes, data_bundle= '', '', {}, {}, [], []
    try:
        timestamp = fn_rule.findall(os.path.basename(filename))[0]
    except IndexError:
        timestamp = '000000'

    axis_number = 1
    while True:
        try:
            # try to open up another AXIS object in the HDF's attribute directory
            #  (they are named /AXIS/AXIS1, /AXIS/AXIS2, /AXIS/AXIS3 ...)
            axis_to_look_for = axis_name + str(axis_number)
            axis = data_file[axis_to_look_for]
            # convert byte string attributes to string
            attrs = {}
            for k, v in axis.attrs.items():
                try:
                    attrs[k] = v[0].decode('utf-8') if isinstance(v[0], bytes) else v
                except IndexError:
                    attrs[k] = v.decode('utf-8') if isinstance(v, bytes) else v

            axis_min = axis[0]
            axis_max = axis[-1]
            axis_numberpoints = n_data[0].shape[-(axis_number-1)]

            data_axis = DataAxis(axis_min, axis_max, axis_numberpoints, attrs=attrs)
            axes.insert(0, data_axis)
        except KeyError:
            break
        axis_number += 1

    # we need a loop here primarily (I think) for n_ene_bin phasespace data
    the_data_hdf_object = n_data
    # for the_data_hdf_object in n_data:
        
    name = the_data_hdf_object.name[1:]  # ignore the beginning '/'

    # now read in attributes of the ROOT of the hdf5..
    #   there's lots of good info there. strip out the array if value is a string

    for key, value in data_file.attrs.items():
        try:
            run_attrs[key] = value[0].decode('utf-8') if isinstance(value[0], bytes) else value
        except IndexError:
            run_attrs[key] = value.decode('utf-8') if isinstance(value, bytes) else value
    try:
        run_attrs['TIME UNITS'] = OSUnits(run_attrs['TIME UNITS'])
    except:
        run_attrs['TIME UNITS'] = OSUnits('1 / \omega_p')
    # attach attributes assigned to the data array to
    #    the H5Data.data_attrs object, remove trivial dimension before assignment
    for key, value in the_data_hdf_object.attrs.items():
        try:
            data_attrs[key] = value[0].decode('utf-8') if isinstance(value[0], bytes) else value
        except IndexError:
            data_attrs[key] = value.decode('utf-8') if isinstance(value, bytes) else value

    # check if new data format is in use
    if not data_attrs and 'SIMULATION' in data_file:
        data_attrs['LONG_NAME'], data_attrs['UNITS'] = run_attrs.pop('LABEL', 'data'), run_attrs.pop('UNITS', OSUnits('a.u.'))
        run_attrs['SIMULATION'] = {k:v for k,v in data_file['/SIMULATION'].attrs.items()}
    # convert unit string to osunit object
    try:
        data_attrs['UNITS'] = OSUnits(data_attrs['UNITS'])
    except:
#             data_attrs['UNITS'] = OSUnits('a.u.')
        pass
    data_attrs['NAME'] = "test"

    # data_bundle.data = the_data_hdf_object[()]
    data_bundle.append(H5Data(the_data_hdf_object, timestamp=timestamp,
                                data_attrs=data_attrs, run_attrs=run_attrs, axes=axes))
    data_file.close()
    if len(data_bundle) == 1:
        return data_bundle[0]
    else:
        return data_bundle


def extract_number(filename):
    """
    Extracts the last number found in a filename using regular expression.
    If no number is found, returns None.
    """
    matches = re.findall(r'\d+', filename)
    return int(matches[-1]) if matches else None

def order_files_by_number(directory, dataset='e1'):
    """
    Orders files in the given directory based on the last number in their filenames.
    """
    odir = os.path.join(directory, 'MS', 'FLD', dataset)

    files = [file.split('.h5')[0] for file in os.listdir(odir) if file.endswith('.h5')]

    files_with_numbers = [(file, extract_number(file)) for file in files]
    # Filter out files where no number was found
    files_with_numbers = [fn for fn in files_with_numbers if fn[1] is not None]
    
    #find 
    # Sort files by the extracted number
    sorted_files = sorted(files_with_numbers, key=lambda x: x[1])
    return [odir+"/"+fn[0]+".h5" for fn in sorted_files]

def generate_plot(frame_number, flist, is1d, axes, xmult=1, ymult=1, intensitymult=1):
    """
    Function to generate each plot. You can modify this function to create different plots.
    """
    axes[0].clear()
    axes[1].clear()
    ax=axes[0]
    ax2=axes[1]
    fhere = h5py.File(flist[frame_number], 'r')
    print(frame_number)
    if(len(fhere['AXIS']) == 2):
        xaxismin = fhere['AXIS']['AXIS1'][0]
        xaxismax = fhere['AXIS']['AXIS1'][1]
        yaxismin = fhere['AXIS']['AXIS2'][0]
        yaxismax = fhere['AXIS']['AXIS2'][1]
        if is1d:
            num_pointsx = fhere[dataset].shape[1]
            num_pointsy = fhere[dataset].shape[0]

            x = xmult*np.linspace(start=xaxismin, stop=xaxismax, num=num_pointsx)
            y = ymult*np.mean(fhere[dataset], axis=0)
            y = ymult*fhere[dataset][int(num_pointsy/2),:]
            ax.plot(x, y, "r")
            ax.set_ylim(-intensitymult, intensitymult)

            try:
                fdens = flist[frame_number].split("MS")[0]+"MS/PHA/x2x1_m/electrons/x2x1_m-electrons-"+flist[frame_number].split(".h5")[0].split("-")[-1]+".h5"
                fheredens = h5py.File(fdens, "r")
                y = np.mean(fheredens["x2x1_m"], axis=0)
                ylimabs = 0.5
                ax2.plot(x,ylimabs*y)
                ax2.set_ylim(0, 0.75)

            except:
                print("No density data!")
                
        else:
            plt.imshow(fhere[dataset][:,:]+1e-12,
                aspect='auto',
                extent=[xaxismin, xaxismax, yaxismin, yaxismax], cmap="RdBu", vmin=-0.1, vmax=0.1)
            
            try:
                fdens = flist[frame_number].split("MS")[0]+"MS/PHA/x2x1_m/electrons/x2x1_m-electrons-"+flist[frame_number].split(".h5")[0].split("-")[-1]+".h5"
                fheredens = h5py.File(fdens, "r")
                plt.imshow(fheredens["x2x1_m"][:,:]+1e-12,
                extent=[xaxismin, xaxismax, yaxismin, yaxismax], alpha=0.2, cmap="grey")
            except:
                print("No density data!")
                # plt.colorbar(orientation='vertical')

def create_movie(filename, num_frames=100, fps=10, flist=[], xmult=1, ymult=1, intensitymult=1):
    """
    Function to create a movie from a sequence of plots.
    """
    # Ensure save directory exists #
    mkdir = flist[0].split("MS")[0]+"Outputs/"
    filename = mkdir+filename
    try:
        os.mkdir(mkdir)
    except:
        print("Cannot make directory to output plots and videos.")

    # Create a figure and axis for the plot
    fig, ax = plt.subplots()
    ax2 = ax.twinx()

    # Creating an animation by updating the plot for each frame
    animation = FuncAnimation(fig, generate_plot, frames=num_frames, fargs=(flist, True, [ax, ax2], xmult, ymult, intensitymult))

    # Save the animation as a GIF (you can also save it as mp4 or other formats)
    writer = PillowWriter(fps=fps)
    animation.save(filename, writer=writer)

    # Convert GIF to movie using moviepy
    clip = mpy.VideoFileClip(filename)
    clip.write_videofile(filename.replace('.gif', '.mp4'))

def find_critical_dens_ratio(wl, dens):
    # wl in micron #
    # dens in /m^3
    n_crit = (10**6)*(1.1e21)/(wl**2)
    print(n_crit)
    Rh_ratio = dens/n_crit
    return Rh_ratio

def find_critical_dens(wl):
    # wl in micron #
    # dens in /m^3
    n_crit = (10**6)*(1.1e21)/(wl**2)
    return n_crit

def field(rundir='',dataset='e1',time=0,space=-1,
    xlim=[-1,-1],ylim=[-1,-1],zlim=[-1,-1], xmult=1, ymult=1,intensitymult=1,
    plotdata=[], color=None, to_plot=True, **kwargs):

    workdir = os.getcwd()
    workdir = os.path.join(workdir, rundir)

    odir = os.path.join(workdir, 'MS', 'FLD', dataset)
    files = sorted(os.listdir(odir))

    i = 0
    for j in range(len(files)):
        fhere = h5py.File(os.path.join(odir,files[j]), 'r')
        if(fhere.attrs['TIME'] >= time):
            i = j
            break

    fhere = h5py.File(os.path.join(odir,files[i]), 'r')

    plt.figure(figsize=(12, 6))
    plt.title(dataset+' field at t = '+str(fhere.attrs['TIME']))
    plt.xlabel('$x_1 [c/\omega_p]$')
    if(len(fhere['AXIS']) == 1):
        plt.ylabel('$n [n_0]$')
    if(len(fhere['AXIS']) == 2):
        plt.ylabel('$x_2 [c/\omega_p]$')

    if(len(fhere['AXIS']) == 1):

        xaxismin = xmult*fhere['AXIS']['AXIS1'][0]
        xaxismax = xmult*fhere['AXIS']['AXIS1'][1]

        nx = len(fhere[dataset][:])
        dx = (xaxismax-xaxismin)/nx

        plt.plot(np.arange(0,xaxismax,dx),np.abs(fhere[dataset][:]))

    elif(len(fhere['AXIS']) == 2):

        xaxismin = xmult*fhere['AXIS']['AXIS1'][0]
        xaxismax = xmult*fhere['AXIS']['AXIS1'][1]
        yaxismin = ymult*fhere['AXIS']['AXIS2'][0]
        yaxismax = ymult*fhere['AXIS']['AXIS2'][1]

        if color != None:
            plt.imshow(intensitymult*fhere[dataset][:,:]+1e-12,
                       aspect='auto',
                       extent=[xaxismin, xaxismax, yaxismin, yaxismax], cmap=color)
            # plt.imshow(intensitymult*fhere[dataset][:,:]+1e-12,
            #         aspect='auto',
            #         extent=[xaxismin, xaxismax, yaxismin, yaxismax], cmap=color, vmin=-0.006, vmax=0.006)           
        else:
            plt.imshow(intensitymult*fhere[dataset][:,:]+1e-12,
                    aspect='auto',
                    extent=[xaxismin, xaxismax, yaxismin, yaxismax])
        plt.colorbar(orientation='vertical')


    if(xlim != [-1,-1]):
        plt.xlim(xlim)
    if(ylim != [-1,-1]):
        plt.ylim(ylim)
    if(zlim != [-1,-1]):
        plt.clim(zlim)

    if to_plot:
        plt.show()

def phasespace(rundir='',dataset='p1x1',species='electrons',time=0,
    xlim=[-1,-1],ylim=[-1,-1],zlim=[-1,-1], xmult=1, ymult=1,
    plotdata=[], color=None, to_plot=True):

    workdir = os.getcwd()
    workdir = os.path.join(workdir, rundir)

    odir = os.path.join(workdir, 'MS', 'PHA', dataset, species)
    files = sorted(os.listdir(odir))

    i = 0
    for j in range(len(files)):
        fhere = h5py.File(os.path.join(odir,files[j]), 'r')
        if(fhere.attrs['TIME'] >= time):
            i = j
            break

    fhere = h5py.File(os.path.join(odir,files[i]), 'r')

    plt.figure(figsize=(12, 6))
    plt.title(dataset+' phasespace at t = '+str(fhere.attrs['TIME']))
    plt.xlabel('$x_1 [c/\omega_p]$')
    if(len(fhere['AXIS']) == 1):
        plt.ylabel('$n [n_0]$')
    if(len(fhere['AXIS']) == 2):
        plt.ylabel('$x_2 [c/\omega_p]$')

    if(len(fhere['AXIS']) == 1):

        xaxismin = xmult*fhere['AXIS']['AXIS1'][0]
        xaxismax = xmult*fhere['AXIS']['AXIS1'][1]

        nx = len(fhere[dataset][:])
        dx = (xaxismax-xaxismin)/nx

        plt.plot(np.arange(0,xaxismax,dx),np.abs(fhere[dataset][:]))

    elif(len(fhere['AXIS']) == 2):

        xaxismin = xmult*fhere['AXIS']['AXIS1'][0]
        xaxismax = xmult*fhere['AXIS']['AXIS1'][1]
        yaxismin = ymult*fhere['AXIS']['AXIS2'][0]
        yaxismax = ymult*fhere['AXIS']['AXIS2'][1]

        if color != None:
            plt.imshow(fhere[dataset][:,:]+1e-12,
                       aspect='auto',
                       extent=[xaxismin, xaxismax, yaxismin, yaxismax], cmap=color)
        else:
            plt.imshow(fhere[dataset][:,:]+1e-12,
                    aspect='auto',
                    extent=[xaxismin, xaxismax, yaxismin, yaxismax])
        plt.colorbar(orientation='vertical')


    if(xlim != [-1,-1]):
        plt.xlim(xlim)
    if(ylim != [-1,-1]):
        plt.ylim(ylim)
    if(zlim != [-1,-1]):
        plt.clim(zlim)
    if to_plot:
        plt.show()

def make_contour2(rundir='',dataset='p1x1',species='electrons',time=0, line_out_x=0,
    xlim=[-1,-1],ylim=[-1,-1],zlim=[-1,-1], xmult=1, ymult=1,
    plotdata=[], color=None, to_plot=True):

    workdir = os.getcwd()
    workdir = os.path.join(workdir, rundir)

    odir = os.path.join(workdir, 'MS', 'PHA', dataset, species)
    files = sorted(os.listdir(odir))

    i = 0
    for j in range(len(files)):
        fhere = h5py.File(os.path.join(odir,files[j]), 'r')
        if(fhere.attrs['TIME'] >= time):
            i = j
            break
    print(os.path.join(odir,files[i]))
    phase_space = np.abs(read_h5(os.path.join(odir,files[i])))

    fig, (phase_plot, ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [3, 1]})

    runatts = phase_space.run_attrs
    title=runatts['NAME']
    time=phase_space.run_attrs['TIME'][0]

    # fig.suptitle('Time = '+repr(time)+'$\omega_p^{-1}$',fontsize=24)
    ext_stuff=[xmult*phase_space.axes[1].min,xmult*phase_space.axes[1].max,ymult*phase_space.axes[0].min,ymult*phase_space.axes[0].max]
    data_max=max(np.abs(np.amax(phase_space)),100)

    # phase_contour=plt.contourf(np.abs(phase_space+0.000000001),
    #             levels=[0.00001*data_max,0.0001*data_max,0.001*data_max,0.01*data_max,0.05*data_max,0.1*data_max,0.2*data_max,0.5*data_max],
    #             extent=ext_stuff,cmap='Spectral',vmin=1e-5*data_max,vmax=1.5*data_max,
    #             norm=colors.LogNorm(vmin=0.00001*data_max,vmax=1.5*data_max))
    # print(np.abs(phase_space))
    xaxis = np.linspace(ext_stuff[0], ext_stuff[1], phase_space.shape[1])
    yaxis = np.linspace(ext_stuff[2], ext_stuff[3], phase_space.shape[0])
    XX, YY = np.meshgrid(xaxis, yaxis)

    phase_contour=phase_plot.pcolormesh(XX, YY, phase_space, cmap='Spectral')
    # phase_contour=phase_plot.contourf(phase_space,extent=ext_stuff, cmap='Spectral', levels=100)
    phase_plot.set_title(species+' '+dataset+' at t = '+str(time))
    phase_plot.set_xlabel('Position [$c / \omega_{p}$]')
    phase_plot.set_ylabel('Proper Velocity $\gamma v_1$ [ c ]')
    # Mark the line-out location on the 2D plot
    phase_plot.axvline(x=line_out_x, color='red', linestyle='--')
    cbar = fig.colorbar(phase_contour, ax=phase_plot)

    # Create the side 1D plot
    line_out_index = np.argmin(np.abs(xaxis - line_out_x))
    line_out_values = phase_space[:, line_out_index]

    ## Fit to Gaussian ##
    pfit = fit_to_gaussian(yaxis, line_out_values)
    y_fit = gaussian(yaxis, *pfit)
    m_e = 9.109383e-31
    elc = 1.602177e-19
    clight = 3e8
    temp = clight**2*m_e*((pfit[2])**2)/(elc) # in eV
    print("Temperature: "+str(temp)+" eV")

    # Now also plot the 1 eV Maxwellian
    maxwellian = gaussian(yaxis, pfit[0], pfit[1], np.sqrt(elc/(m_e*1*clight**2)))
    ax2.plot(maxwellian, yaxis, color='black', linestyle='--')
    ax2.plot(line_out_values, yaxis, color='blue')
    ax2.plot(y_fit, yaxis, color='red')
    ax2.set_title("1D Line-out at x={}".format(line_out_x))
    ax2.set_xlabel("Value")
    ax2.set_ylabel("Y axis")
    ax2.yaxis.tick_right()

    plt.tight_layout()

    if to_plot:
        plt.show()  

def make_contour(rundir='',dataset='p1x1',species='electrons',time=0,
    xlim=[-1,-1],ylim=[-1,-1],zlim=[-1,-1], xmult=1, ymult=1,
    plotdata=[], color=None, to_plot=True):

    workdir = os.getcwd()
    workdir = os.path.join(workdir, rundir)

    odir = os.path.join(workdir, 'MS', 'PHA', dataset, species)
    files = sorted(os.listdir(odir))

    i = 0
    for j in range(len(files)):
        fhere = h5py.File(os.path.join(odir,files[j]), 'r')
        if(fhere.attrs['TIME'] >= time):
            i = j
            break
    print(os.path.join(odir,files[i]))
    phase_space = np.abs(read_h5(os.path.join(odir,files[i])))

    fig = plt.figure(figsize=(12, 6))
    phase_plot=plt.subplot(111 )
    runatts = phase_space.run_attrs
    title=runatts['NAME']
    time=phase_space.run_attrs['TIME'][0]

    # fig.suptitle('Time = '+repr(time)+'$\omega_p^{-1}$',fontsize=24)
    ext_stuff=[xmult*phase_space.axes[1].min,xmult*phase_space.axes[1].max,ymult*phase_space.axes[0].min,ymult*phase_space.axes[0].max]

    data_max=max(np.abs(np.amax(phase_space)),100)
    # phase_contour=plt.contourf(np.abs(phase_space+0.000000001),
    #             levels=[0.00001*data_max,0.0001*data_max,0.001*data_max,0.01*data_max,0.05*data_max,0.1*data_max,0.2*data_max,0.5*data_max],
    #             extent=ext_stuff,cmap='Spectral',vmin=1e-5*data_max,vmax=1.5*data_max,
    #             norm=colors.LogNorm(vmin=0.00001*data_max,vmax=1.5*data_max))
    # print(np.abs(phase_space))
    phase_contour=plt.contourf(phase_space,extent=ext_stuff, cmap='Spectral', levels=100)
    phase_plot.set_title(species+' '+dataset+' at t = '+str(time))
    phase_plot.set_xlabel('Position [$c / \omega_{p}$]')
    phase_plot.set_ylabel('Proper Velocity $\gamma v_1$ [ c ]')
    cbar = fig.colorbar(phase_contour)
    if to_plot:
        plt.show()  
# Define the Gaussian function
def gaussian(x, a, b, c):
    return a * np.exp(-0.5 * ((x - b) / c) ** 2)

def fit_to_gaussian(x_data, y_data):

    # Fit the Gaussian to the data
    initial_guess = [10000, 0, .01]  # Initial guess for the parameters [a, b, c]
    params, covariance = curve_fit(gaussian, x_data, y_data, p0=initial_guess)

    return params

datadir = '/Volumes/T9/XSPL/Lasers/Simulations/Bluehive/OSIRIS/LasersDeck/'
# datadir = 'D:\XSPL\Lasers\Simulations\Bluehive\OSIRIS\LasersDeck\'

# Example: External drive is assigned drive letter 'E:'
drive_letter = 'D:'
file_path_on_external_drive = 'XSPL/Lasers/Simulations/Bluehive/OSIRIS/LasersDeck/' 

# Construct the full file path
full_file_path = drive_letter + '\\' + file_path_on_external_drive
datadir = full_file_path

# List directory contents #
# print(os.listdir('../../../../../../../../../../../../../'))

# datadir = os.getcwd()+"/Lasers/Simulations/Bluehive/OSIRIS/LasersDeck/"

dirname = datadir+'Laser2D'
dirname = datadir+'Laser2D_n_ncrit_0p5'
# dirname = datadir+'Laser1D_n_ncrit_0p5'

# dirname = datadir+'Laser2D_n_ncrit_laserunits'
# dirname = datadir+'Laser2D_n_ncrit_lu_pulse'
# dirname = datadir+'langdon-fixed'

dataset = 'e3'

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
sorted_files = order_files_by_number(directory=dirname, dataset=dataset)
# print(sorted_files)

ncrit_m3 = find_critical_dens(0.532)
ncrit_cm3 = ncrit_m3*(1e-6)

time = 20
# make_contour(rundir=dirname,dataset='p3x1',time=time, xlim=[0,12],  xmult=clight/omega_p/wavelength, ymult=1, species='electrons', to_plot=False)
make_contour2(rundir=dirname,dataset='p3x1',time=time, xlim=[0,12], line_out_x = 6, xmult=clight/omega_p/wavelength, ymult=1, species='electrons', to_plot=False)
make_contour2(rundir=dirname,dataset='p3x2',time=time, xlim=[0,12], line_out_x = 6, xmult=clight/omega_p/wavelength, ymult=1, species='electrons', to_plot=False)
phasespace(rundir=dirname,dataset='x2x1_ene',time=time, xlim=[0,12], xmult=clight/omega_p/wavelength, ymult=clight/omega_p/wavelength, species='electrons', color="Reds", to_plot=False)
phasespace(rundir=dirname,dataset='x2x1_m',time=time, xlim=[0,12], xmult=clight/omega_p/wavelength, ymult=clight/omega_p/wavelength, species='electrons', color="copper", to_plot=False)
field(rundir=dirname,dataset='e3',time=time,xlim=[0,12], xmult=clight/omega_p/wavelength, ymult=clight/omega_p/wavelength, intensitymult=Eamp, color='RdBu')

# Add each plot from the field function to a frame of a movie and store this movie as a file #
if (True):
    create_movie(filename=dataset+".gif", num_frames=len(sorted_files), fps=5, flist=sorted_files, xmult=clight/omega_p/wavelength, ymult=Eamp, intensitymult=Emax)
