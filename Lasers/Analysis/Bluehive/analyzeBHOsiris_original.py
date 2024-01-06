import matplotlib.pyplot as plt
import h5py
import numpy as np
import os
from analysis import *

def field(rundir='',dataset='e1',time=0,space=-1,
    xlim=[-1,-1],ylim=[-1,-1],zlim=[-1,-1],
    plotdata=[], **kwargs):

    if(space != -1):
        plot_or = 1
        PATH = gen_path(rundir, plot_or)
        hdf5_data = read_hdf(PATH)
        plt.figure(figsize=(8,5))
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

    plt.figure(figsize=(6, 3.2))
    plt.title(dataset+' at t = '+str(fhere.attrs['TIME']))
    plt.xlabel('$x_1 [c/\omega_p]$')
    plt.ylabel(dataset)

    xaxismin = fhere['AXIS']['AXIS1'][0]
    xaxismax = fhere['AXIS']['AXIS1'][1]

    nx = len(fhere[dataset][:])
    dx = (xaxismax-xaxismin)/nx

    plt.plot(np.arange(0,xaxismax,dx),fhere[dataset][:])

    if(xlim != [-1,-1]):
        plt.xlim(xlim)
    if(ylim != [-1,-1]):
        plt.ylim(ylim)
    if(zlim != [-1,-1]):
        plt.clim(zlim)

    plt.show()

def phasespace(rundir='',dataset='p1x1',species='electrons',time=0,
    xlim=[-1,-1],ylim=[-1,-1],zlim=[-1,-1],
    plotdata=[]):

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

    plt.figure(figsize=(6, 3.2))
    plt.title(dataset+' phasespace at t = '+str(fhere.attrs['TIME']))
    plt.xlabel('$x_1 [c/\omega_p]$')
    if(len(fhere['AXIS']) == 1):
        plt.ylabel('$n [n_0]$')
    if(len(fhere['AXIS']) == 2):
        plt.ylabel('$p_1 [m_ec]$')

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

        plt.imshow(np.log(np.abs(fhere[dataset][:,:]+1e-12)),
                   aspect='auto',
                   extent=[xaxismin, xaxismax, yaxismin, yaxismax])
        plt.colorbar(orientation='vertical')


    if(xlim != [-1,-1]):
        plt.xlim(xlim)
    if(ylim != [-1,-1]):
        plt.ylim(ylim)
    if(zlim != [-1,-1]):
        plt.clim(zlim)

    plt.show()

dirname = '/Users/james/Desktop/Folders/XSPL/Lasers/Simulations/Bluehive/OSIRIS/LasersDeck/Laser2D'
phasespace(rundir=dirname,dataset='x2x1',time=0,xlim=[0,37.69908])