from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import h5py
import numpy as np
from osh5def import H5Data, PartData, fn_rule, DataAxis, OSUnits
import os
import re



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


# Define the Gaussian function
def gaussian(x, a, b, c):
    return a * np.exp(-0.5 * ((x - b) / c) ** 2)

def fit_to_gaussian(x_data, y_data):

    # Fit the Gaussian to the data
    initial_guess = [10000, 0, .01]  # Initial guess for the parameters [a, b, c]
    params, covariance = curve_fit(gaussian, x_data, y_data, p0=initial_guess)

    return params

def extract_number(filename):
    """
    Extracts the last number found in a filename using regular expression.
    If no number is found, returns None.
    """
    matches = re.findall(r'\d+', filename)
    return int(matches[-1]) if matches else None

def find_files_associated_with_dataset(directory, dataset='e1'):
    """
    Finds all files in the given directory that are associated with the given dataset.
    """
    workdir = os.getcwd()
    workdir = os.path.join(workdir, directory)
    odir = os.path.join(workdir, 'MS', 'FLD', dataset)

    files = [file.split('.h5')[0] for file in os.listdir(odir) if file.endswith('.h5')]
    return [odir+"/"+file+".h5" for file in files]

def scan_hdf5_file_for_main_data_array(h5_file):
    """
    Scans the hdf5 file for the main data array. This is the array that has the most dimensions.
    """
    data_arrays = []
    for key, value in h5_file.items():
        if isinstance(value, h5py.Dataset):
            data_arrays.append(value)
    return sorted(data_arrays, key=lambda x: len(x.shape), reverse=True)  


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