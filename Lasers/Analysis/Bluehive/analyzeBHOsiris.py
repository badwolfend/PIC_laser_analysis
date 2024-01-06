import matplotlib.pyplot as plt
import h5py
import numpy as np
import os
from analysis import *
from matplotlib.animation import PillowWriter, FuncAnimation
import moviepy.editor as mpya

def generate_plot(frame_number, ax):
    """
    Function to generate each plot. You can modify this function to create different plots.
    """
    ax.clear()

    ax.plot(x, y)

def create_movie(filename, num_frames=100, fps=10):
    """
    Function to create a movie from a sequence of plots.
    """
    # Create a figure and axis for the plot
    fig, ax = plt.subplots()

    # Creating an animation by updating the plot for each frame
    animation = FuncAnimation(fig, generate_plot, frames=num_frames, fargs=(ax,))

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
    xlim=[-1,-1],ylim=[-1,-1],zlim=[-1,-1],
    plotdata=[], color=None, to_plot=True, **kwargs):

    if(space != -1):
        plot_or = 1
        PATH = gen_path(rundir, plot_or)
        hdf5_data = read_hdf(PATH)
        fig, ax = plt.subplots()
        #plotme(hdf5_data.data[:,space], **kwargs)
        plotme(hdf5_data,hdf5_data.data[space,:])
        plt.title('temporal evolution of e' + str(plot_or) + ' at cell ' + str(space))
        plt.xlabel('t')
        plt.show()
        return


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

        xaxismin = fhere['AXIS']['AXIS1'][0]
        xaxismax = fhere['AXIS']['AXIS1'][1]

        nx = len(fhere[dataset][:])
        dx = (xaxismax-xaxismin)/nx

        plt.plot(np.arange(0,xaxismax,dx),np.abs(fhere[dataset][:]))

    elif(len(fhere['AXIS']) == 2):

        xaxismin = fhere['AXIS']['AXIS1'][0]
        xaxismax = fhere['AXIS']['AXIS1'][1]
        yaxismin = fhere['AXIS']['AXIS2'][0]
        yaxismax = fhere['AXIS']['AXIS2'][1]

        if color != None:
            plt.imshow(fhere[dataset][:,:]+1e-12,
                       aspect='auto',
                       extent=[xaxismin, xaxismax, yaxismin, yaxismax], cmap=color, vmin=-0.1, vmax=0.1)
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

def phasespace(rundir='',dataset='p1x1',species='electrons',time=0,
    xlim=[-1,-1],ylim=[-1,-1],zlim=[-1,-1],
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

        xaxismin = fhere['AXIS']['AXIS1'][0]
        xaxismax = fhere['AXIS']['AXIS1'][1]

        nx = len(fhere[dataset][:])
        dx = (xaxismax-xaxismin)/nx

        plt.plot(np.arange(0,xaxismax,dx),np.abs(fhere[dataset][:]))

    elif(len(fhere['AXIS']) == 2):

        xaxismin = fhere['AXIS']['AXIS1'][0]
        xaxismax = fhere['AXIS']['AXIS1'][1]
        yaxismin = fhere['AXIS']['AXIS2'][0]
        yaxismax = fhere['AXIS']['AXIS2'][1]

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

dirname = '/Users/james/Desktop/Folders/XSPL/Lasers/Simulations/Bluehive/OSIRIS/LasersDeck/Laser2D'
dirname = '/Users/james/Desktop/Folders/XSPL/Lasers/Simulations/Bluehive/OSIRIS/LasersDeck/Laser2D_n_ncrit'
dirname = '/Users/james/Desktop/Folders/XSPL/Lasers/Simulations/Bluehive/OSIRIS/LasersDeck/Laser2D_n_ncrit_laserunits'
dirname = '/Users/james/Desktop/Folders/XSPL/Lasers/Simulations/Bluehive/OSIRIS/LasersDeck/Laser2D_n_ncrit_lu_pulse'

ncrit_m3 = find_critical_dens(0.532)
ncrit_cm3 = ncrit_m3*(1e-6)
print("ncrit = ", ncrit_cm3/(4*3.14159**2), "cm^-3")
print("ncrit = ", ncrit_m3, "m^-3")

phasespace(rundir=dirname,dataset='x2x1',time=5, xlim=[0,12], species='electrons', to_plot=False)
field(rundir=dirname,dataset='e3',time=5,xlim=[0,12], color='RdBu')

# Add each plot from the field function to a frame of a movie and store this movie as a file #

