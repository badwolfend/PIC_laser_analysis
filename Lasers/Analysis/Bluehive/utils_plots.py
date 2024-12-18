from osh5def import H5Data, PartData, fn_rule, DataAxis, OSUnits
import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter, FuncAnimation
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
import moviepy.editor as mpy
import h5py
import numpy as np
from osh5def import H5Data, PartData, fn_rule, DataAxis, OSUnits
import os
import utils as ut
import osh5vis as osh5vis
import pickle
font = {'family' : 'Times',
        'size'   : 22}


def save_temperature_time_series(dataset='e3',times=[0,1],species='electrons', space=-1, xlim=[-1,-1],ylim=[-1,-1],zlim=[-1,-1], tmult=1, xmult=1, ymult=1,intensitymult=1,  plotdata=[], colors=['r'], to_plot=True, to_normalize=False, to_save=True, save_dir='./', line_out_x_0=5.75, line_out_x_1=6.0, **kwargs):
    rundir = kwargs.get('dirname', '')
    save_string = "temperature"    
    temp0_array = np.zeros(len(times))
    temp1_array = np.zeros(len(times))
    data_tosave = {}
    for i, time in enumerate(times):
        print("Time: "+str(time))
        temp_0 = make_contour2(rundir=rundir,dataset=dataset,time=time, xlim=xlim, tmult=tmult, line_out_x =line_out_x_0, xmult=xmult, ymult=ymult, species=species, to_plot=to_plot, to_save=to_save, save_dir=save_dir, to_return_temp=True, to_clear=True)
        temp_1 = make_contour2(rundir=rundir,dataset=dataset,time=time, xlim=xlim, tmult=tmult, line_out_x = line_out_x_1, xmult=xmult, ymult=ymult, species=species, to_plot=to_plot, to_save=to_save, save_dir=save_dir, to_return_temp=True, to_clear=True)
        temp0_array[i] = temp_0
        temp1_array[i] = temp_1
        data_tosave["line_out_x_0"] = temp0_array 
        data_tosave["line_out_x_1"] = temp1_array
        data_tosave["time"] = tmult*times 
        run_name = rundir.split("/")[-1]   
        save_name = save_dir+"Data/"+run_name+"_"+save_string+"_times"
        pickle.dump(data_tosave, open(save_name+".p", "wb"))  
    if to_save:
        # data_tosave["line_out_x_0"] = temp0_array 
        # data_tosave["line_out_x_1"] = temp1_array
        # data_tosave["time"] = tmult*times 
        # run_name = rundir.split("/")[-1]   
        # save_name = save_dir+"Data/"+run_name+"_"+save_string+"_times"
        # pickle.dump(data_tosave, open(save_name+".p", "wb"))  
        plt.savefig(save_dir+"Plots/"+run_name+"_"+save_string+"_time_"+str(i)+".png", dpi=600)
    if to_plot:
        plt.show()


def fields(rundir='',dataset=['e3', 'j3'],mu=1, zeff=1, time=0,space=-1,
    xlim=[-1,-1],ylim=[-1,-1],zlim=[-1,-1], tmult=1, xmult=1, ymult=1,intensitymult=[1],
    plotdata=[], colors=['r'], to_plot=True, to_normalize=False, to_save=True, save_dir='./', use_ions=False, **kwargs):
    de_0 = read_density(rundir=rundir,dataset='m',species='electrons',time=0, line_out_x=0, tmult=1, ylim=[-1,-1],zlim=[-1,-1], xmult=1, ymult=1, color=None, to_plot=False, to_save=False, save_dir=save_dir)

    if use_ions:
        delectron = read_density(rundir=rundir,dataset='m',species='ions',time=time, line_out_x=0, tmult=1, ylim=[-1,-1],zlim=[-1,-1], xmult=1, ymult=1, color=None, to_plot=False, to_save=False, save_dir=save_dir)
    else:
        delectron = read_density(rundir=rundir,dataset='m',species='electrons',time=time, line_out_x=0, tmult=1, ylim=[-1,-1],zlim=[-1,-1], xmult=1, ymult=1, color=None, to_plot=False, to_save=False, save_dir=save_dir)
    
    data_tosave = {}

    fig, ax = plt.subplots(figsize=(12, 6))
    fig.set_size_inches(13.385, 6.0)

    plt.xlabel('$x_1 [\lambda_L]$')
    
    save_string = "vars"
    for id, data in enumerate(dataset):
      save_string += "_"+data  
    #   files = ut.find_files_associated_with_dataset(rundir, dataset=data)
      files = ut.order_files_by_number(rundir, dataset=data)
      i = 0
      for j in range(len(files)):
          fhere = h5py.File(files[j], 'r')
          if(fhere.attrs['TIME'] >= time):
              i = j
              break

      fhere = h5py.File(files[i], 'r')
      plt.title(str(dataset)+' fields at t = '+str(tmult*fhere.attrs['TIME'])+", timestep = "+str(i))

      if(len(fhere['AXIS']) == 1):
          plt.ylabel('$n [n_0]$')
      if(len(fhere['AXIS']) == 2):
          plt.ylabel('$x_2 [c/\omega_p]$')

      if(len(fhere['AXIS']) == 1):

        xaxismin = xmult*fhere['AXIS']['AXIS1'][0]
        xaxismax = xmult*fhere['AXIS']['AXIS1'][1]

        nx = len(fhere[data][:])
        dx = (xaxismax-xaxismin)/nx
        xaxis = np.linspace(start=xaxismin, stop=xaxismax, num=nx) 
        xstart = np.where(xaxis < 4 )[0][-1]
        xstop = np.where(xaxis < 8 )[0][-1] 
        data_tosave[data] = fhere[data][:]

        if to_normalize:
            y = fhere[data][:]
            ynorm = y/np.max(np.abs(y))
            if data == 'j3':
                ax.plot(xaxis[xstart:xstop],ynorm[xstart:xstop], colors[id], linewidth=4, label=data)
            else:
                ax.plot(xaxis, ynorm, colors[id], linewidth=4, label=data)
            ax.set_ylim(-1.1, 1.1)
            ax.set_ylabel("Normalized Field")
        else:
            ax.plot(xaxis,intensitymult[id]*fhere[data][:], linewidth=4, label=data)
            ax.set_ylabel("Field")

        try:
            ax2 = ax.twinx()
            fdens = files[i].split("MS")[0]+"MS/PHA/x1_m/electrons/x1_m-electrons-"+files[i].split(".h5")[0].split("-")[-1]+".h5"
            ion_correction = 1  
            if use_ions:
                fdens = files[i].split("MS")[0]+"MS/PHA/x1_m/ions/x1_m-ions-"+files[i].split(".h5")[0].split("-")[-1]+".h5"
                ion_correction = zeff/(mu*1836)
            fheredens = h5py.File(fdens, "r")
            y = fheredens["x1_m"][:]
            ylimabs = 0.50*ion_correction
            # ylimabs = 1.9803477e+21
            # ylimabs = 1
            yplot = ylimabs*y
            ax2.plot(np.arange(0,xaxismax,dx),yplot, 'k-', linewidth=4, label=data)
            ax2.set_ylim(-0.01, 0.6)
            ax2.set_ylabel("N/N_0")
            data_tosave["x1_m"] = ylimabs*y
            data_tosave["x1_m_density"] = ylimabs*delectron
            data_tosave["x1_m_density_t0"] = ylimabs*de_0
           

        except:
            print("No density data!")

      elif(len(fhere['AXIS']) == 2):

          xaxismin = xmult*fhere['AXIS']['AXIS1'][0]
          xaxismax = xmult*fhere['AXIS']['AXIS1'][1]
          yaxismin = ymult*fhere['AXIS']['AXIS2'][0]
          yaxismax = ymult*fhere['AXIS']['AXIS2'][1]

          if colors != None:
              plt.imshow(intensitymult*fhere[data][:,:]+1e-12,
                        aspect='auto',
                        extent=[xaxismin, xaxismax, yaxismin, yaxismax], cmap=colors[id])
              # plt.imshow(intensitymult*fhere[dataset][:,:]+1e-12,
              #         aspect='auto',
              #         extent=[xaxismin, xaxismax, yaxismin, yaxismax], cmap=color, vmin=-0.006, vmax=0.006)           
          else:
              plt.imshow(intensitymult*fhere[data][:,:]+1e-12,
                      aspect='auto',
                      extent=[xaxismin, xaxismax, yaxismin, yaxismax])
          plt.colorbar(orientation='vertical')


      if(xlim != [-1,-1]):
          plt.xlim(xlim)
      if(ylim != [-1,-1]):
          plt.ylim(ylim)
      if(zlim != [-1,-1]):
          plt.clim(zlim)
    if to_save:
        run_name = rundir.split("/")[-1]   
        save_name = save_dir+"Data/"+run_name+"_"+save_string+"_time_"+str(i)
        pickle.dump(data_tosave, open(save_name+".p", "wb"))  
        plt.savefig(save_dir+"Plots/"+run_name+"_"+save_string+"_time_"+str(i)+".png", dpi=600)
    if to_plot:
        plt.show()

def get_all_times(rundir='',dataset='e3'):
    workdir = os.getcwd()
    workdir = os.path.join(workdir, rundir)

    odir = os.path.join(workdir, 'MS', 'FLD', dataset)
    files = sorted(os.listdir(odir))

    times = []
    for j in range(len(files)):
        fhere = h5py.File(os.path.join(odir,files[j]), 'r')
        times.append(fhere.attrs['TIME'][0])
    return times

def field(rundir='',dataset='e1',time=0,space=-1,
    xlim=[-1,-1],ylim=[-1,-1],zlim=[-1,-1], tmult=1, xmult=1, ymult=1,intensitymult=1,
    plotdata=[], color=None, to_plot=True, **kwargs):

    xout = []
    yout = []

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
    plt.title(dataset+' field at t = '+str(tmult*fhere.attrs['TIME']))
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
        xout =list(np.arange(0,xaxismax,dx))
        yout = intensitymult*fhere[dataset][:]

        plt.plot(np.arange(0,xaxismax,dx),intensitymult*fhere[dataset][:])

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
    return xout, yout

def phasespace(rundir='',dataset='p1x1',species='electrons',time=0,
    xlim=[-1,-1],ylim=[-1,-1],zlim=[-1,-1], tmult=1, xmult=1, ymult=1,
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
    plt.title(dataset+' phasespace at t = '+str(tmult*fhere.attrs['TIME']))
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


def read_density(rundir='',dataset='p1x1',species='electrons',time=0, line_out_x=0, to_fit=True, xlim=[-1,-1], tmult=1, ylim=[-1,-1],zlim=[-1,-1], xmult=1, ymult=1, plotdata=[], color=None, to_plot=False, to_save=False, save_dir='./', to_return_temp=False, to_clear=False):
    data_tosave = {}
    save_string = "vars_"+dataset+"_"+species

    workdir = os.getcwd()
    workdir = os.path.join(workdir, rundir)

    odir = os.path.join(workdir, 'MS', 'DENSITY', species, dataset)
    files = sorted(os.listdir(odir))

    i = 0
    for j in range(len(files)):
        fhere = h5py.File(os.path.join(odir,files[j]), 'r')
        if(fhere.attrs['TIME'][0] >= time):
            i = j
            break
    print(os.path.join(odir,files[i]))

    dens=ut.read_h5(os.path.join(odir,files[i]))

    dens_array = dens.view(np.ndarray)
    if to_plot:
        den_plot = plt.subplot(111)
        osh5vis.osplot(0.5*dens,title=species+' Density')
        plt.show()

    data_tosave["x1_m"] = 0.5*dens_array
    if to_save:
        run_name = rundir.split("/")[-1]   
        save_name = save_dir+"Data/"+run_name+"_"+save_string+"_time_"+str(i)
        pickle.dump(data_tosave, open(save_name+".p", "wb"))  

    return dens_array

def make_contour2(rundir='',dataset='p1x1',species='electrons',time=0, line_out_x=0, to_fit=True, xlim=[-1,-1], tmult=1, ylim=[-1,-1],zlim=[-1,-1], xmult=1, ymult=1, plotdata=[], color=None, to_plot=True, to_save=False, save_dir='./', to_return_temp=False, to_clear=False):
    
    save_string = "vars_"+dataset+"_"+species

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
    phase_space = np.abs(ut.read_h5(os.path.join(odir,files[i])))

    fig, (phase_plot, ax2) = plt.subplots(1, 2, figsize=(13.385, 6), gridspec_kw={'width_ratios': [3, 1]})

    runatts = phase_space.run_attrs
    title=runatts['NAME']
    time=tmult*phase_space.run_attrs['TIME'][0]

    ext_stuff=[xmult*phase_space.axes[1].min,xmult*phase_space.axes[1].max,ymult*phase_space.axes[0].min,ymult*phase_space.axes[0].max]
    data_max=max(np.abs(np.amax(phase_space)),100)
    print("y_min: "+str(phase_space.axes[0].min)+" y_max: "+str(phase_space.axes[0].max))
    xaxis = np.linspace(ext_stuff[0], ext_stuff[1], phase_space.shape[1])
    yaxis = np.linspace(ext_stuff[2], ext_stuff[3], phase_space.shape[0])
    XX, YY = np.meshgrid(xaxis, yaxis)

    phase_contour=phase_plot.pcolormesh(XX, YY, phase_space, cmap='twilight', vmin=0, vmax=60)
    # phase_contour=phase_plot.contourf(phase_space,extent=ext_stuff, cmap='Spectral', levels=100)
    phase_plot.set_xlabel('Position [$\lambda_{L}$]')
    phase_plot.set_ylabel('Proper Velocity $\gamma v_1$ [ c ]', labelpad=20)

    # Mark the line-out location on the 2D plot
    phase_plot.axvline(x=line_out_x, color='red', linewidth=4, linestyle='--')
    # cbar = fig.colorbar(phase_contour, ax=phase_plot)

    # Create the side 1D plot
    line_out_index = np.argmin(np.abs(xaxis - line_out_x))
    line_out_values = phase_space[:, line_out_index]

    ax2.plot(line_out_values, yaxis, color='blue', linewidth=4)
    ax2.set_title("1D Line-out at x={}".format(line_out_x))
    ax2.set_ylim(ext_stuff[2], ext_stuff[3])
    ax2.set_xlabel("Value")
    ax2.yaxis.tick_right()
    
    ## Fit to Gaussian ##
    if to_fit:
        pfit = ut.fit_to_gaussian(yaxis, line_out_values)
        y_fit = ut.gaussian(yaxis, *pfit)
        m_e = 9.109383e-31
        elc = 1.602177e-19
        clight = 3e8
        temp = clight**2*m_e*((pfit[2])**2)/(elc) # in eV
        print("x="+str(line_out_x)+", Temperature: "+str(temp)+" eV")
        phase_plot.set_title(species+' '+dataset+' at t = '+str(time)+', Temperature: '+str(temp)+' eV')

        # Now also plot the 1 eV Maxwellian
        maxwellian = ut.gaussian(yaxis, pfit[0], pfit[1], np.sqrt(elc/(m_e*1*clight**2)))
        ax2.plot(maxwellian, yaxis, color='black', linestyle='--', linewidth=4)
        ax2.plot(line_out_values, yaxis, color='blue', linewidth=4)
        ax2.plot(y_fit, yaxis, color='red', linewidth=4)
        # ax2.set_title("1D Line-out at x={}".format(line_out_x))
        # ax2.set_ylim(ext_stuff[2], ext_stuff[3])
        # ax2.set_xlabel("Value")
        # ax2.yaxis.tick_right()

    plt.tight_layout()

    if to_save:
        run_name = rundir.split("/")[-1]   
        # save_name = save_dir+"Data/"+run_name+"_"+save_string+"_time_"+str(i)
        # pickle.dump(data_tosave, open(save_name+".p", "wb"))  
        plt.savefig(save_dir+"Plots/"+run_name+"_"+save_string+"_time_"+str(i)+".png", dpi=600)
    if to_plot:
        plt.show()
    else:
        if to_clear:
            plt.close() 

    if to_return_temp:
        return temp
    

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
    phase_space = np.abs(ut.read_h5(os.path.join(odir,files[i])))

    fig = plt.figure(figsize=(12, 6))
    phase_plot=plt.subplot(111 )
    runatts = phase_space.run_attrs
    title=runatts['NAME']
    time=phase_space.run_attrs['TIME'][0]

    ext_stuff=[xmult*phase_space.axes[1].min,xmult*phase_space.axes[1].max,ymult*phase_space.axes[0].min,ymult*phase_space.axes[0].max]

    data_max=max(np.abs(np.amax(phase_space)),100)

    phase_contour=plt.contourf(phase_space,extent=ext_stuff, cmap='Spectral', levels=100)
    phase_plot.set_title(species+' '+dataset+' at t = '+str(time))
    phase_plot.set_xlabel('Position [$c / \omega_{p}$]')
    phase_plot.set_ylabel('Proper Velocity $\gamma v_1$ [ c ]')
    cbar = fig.colorbar(phase_contour)
    if to_plot:
        plt.show()  

def generate_plot(frame_number, flist, is1d, axes, dataset="e3", xmult=1, ymult=1, intensitymult=1):
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
            # y = ymult*np.mean(fhere[dataset], axis=0)
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

    if(len(fhere['AXIS']) == 1):
        xaxismin = xmult*fhere['AXIS']['AXIS1'][0]
        xaxismax = xmult*fhere['AXIS']['AXIS1'][1]

        num_pointsx = len(fhere[dataset][:])
        dx = (xaxismax-xaxismin)/num_pointsx

        # Setup the x axis #
        x = np.linspace(start=xaxismin, stop=xaxismax, num=num_pointsx)

        # Plot the data #
        y = ymult*fhere[dataset]
        ax.plot(x, y, "r")
        ax.set_ylim(-intensitymult, intensitymult)

        try:
                fdens = flist[frame_number].split("MS")[0]+"MS/PHA/x1_m/electrons/x1_m-electrons-"+flist[frame_number].split(".h5")[0].split("-")[-1]+".h5"
                fheredens = h5py.File(fdens, "r")
                y = fheredens["x1_m"][:]
                ylimabs = 0.5
                ax2.plot(x,ylimabs*y)
                ax2.set_ylim(0, 0.75)

        except:
            print("No density data!")

def create_movie(filename, num_frames=100, fps=10, flist=[], dataset='e3', xmult=1, ymult=1, intensitymult=1):
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
    animation = FuncAnimation(fig, generate_plot, frames=num_frames, fargs=(flist, True, [ax, ax2], dataset, xmult, ymult, intensitymult))

    # Save the animation as a GIF (you can also save it as mp4 or other formats)
    writer = PillowWriter(fps=fps)
    animation.save(filename, writer=writer)

    # Convert GIF to movie using moviepy
    clip = mpy.VideoFileClip(filename)
    clip.write_videofile(filename.replace('.gif', '.mp4'))