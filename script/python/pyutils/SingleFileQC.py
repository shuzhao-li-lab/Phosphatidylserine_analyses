# key functions

import numpy as np
import matplotlib.pyplot as plt
from pyopenms import *

def calcTIC(exp, mslevel=1):
    rt, tic = [], []
    # Iterate through all spectra of the experiment
    for spec in exp:
        # Only calculate TIC for matching (MS1) spectra
        if spec.getMSLevel() == mslevel:
            rt.append(spec.getRT())
            tic.append(sum(spec.get_peaks()[1]))

    return rt, tic

def extract_trio(exp, mslevel=1, min_intensity=1000):
    # return a list of [ (intensities, mz, rt), ...]
    trio_list = []
    for spec in exp:
        if spec.getMSLevel() == mslevel:
            this_rt = spec.getRT()
            # in OpenMS, spec.get_peaks() returns two lists: m/z and intensity in the scan.
            for m, i in zip(*spec.get_peaks()):
                if i > min_intensity:
                    # int(i) to save memory
                    trio_list.append( (int(i), m, this_rt) )
                    
    return trio_list

def get_targeted_mzrange(trio_list, mz_low, mz_high):
    # return trio_list sorted by intensity
    new = []
    for t in trio_list:
        # must in order of i, m, rt
        if mz_low < t[1] < mz_high:
            new.append(t)

    new.sort(reverse=True)
    return new

# warning: low intensity could indicate wrong peak identified
def plot_target_mzchrompeaks(full_trio_list, 
                             target_name,
                             target_mz, mz_error=0.003, 
                             target_rt=None, rt_error=10, min_intensity=10000,
                             figsize=(10,5) ):
    # rt not used right now
    range_low, range_high = target_mz-mz_error, target_mz+mz_error
    # sorted by intensity, must in order of i, m, rt
    use_data_points = get_targeted_mzrange( full_trio_list, range_low, range_high )

    # ax1 for m/z plot using log2 intensity, ax2 for chromatogram in linear intensity
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=figsize)
    
    # log2min_intensity = np.log2(min_intensity)
    RT, MZ, INS = [x[2] for x in use_data_points], [x[1] for x in use_data_points], [x[0] for x in use_data_points]
    log2INS = [np.log2(x) for x in INS]
    
    
    
    ax1.plot(MZ, log2INS, '.')
    ax1.set_xlim(range_low, range_high)
    ax1.set_title('m/z peaks')
    ax1.set_xlabel('m/z')
    ax1.set_ylabel('log2 intensity')
    #       
    ax2.plot(RT, INS, '.')
    ax2.set_title('chrom peaks')
    ax2.set_xlabel('time')
    ax2.set_ylabel('intensity')
    
    if INS and INS[0] > min_intensity:
        ppm = (MZ[0] - target_mz)/target_mz * 1e6
        fig.suptitle(target_name + '\n' +
                'top intensity at: %s, \nmass accuracy = %s ppm' %(str(MZ[0])[:9], str(ppm)[:5]))
    else:
        fig.suptitle(target_name + '\n' +
                "No signal detected above %d." %min_intensity)
                     
    plt.tight_layout()
    plt.show()
    
def plot_file(file, standard_list):
    exp = MSExperiment()
    MzMLFile().load(file, exp)
    
    # do TIC
    retention_times, intensities = calcTIC(exp)
    plt.plot(retention_times, intensities)
    plt.title(file + '\nTIC')
    plt.xlabel('time')
    plt.ylabel('intensity')
    plt.tight_layout()
    plt.show()
    
    # read trio data 
    full_trio_list = extract_trio(exp)

    # do internal standards
    for k in standard_list:
        plot_target_mzchrompeaks(full_trio_list, 
                             k[0],
                             k[1], mz_error=0.003, 
                             min_intensity=10000)
        
import matplotlib.pyplot as plt

def plot_TIC_per_file(file,savefig = False):
    exp = MSExperiment()
    MzMLFile().load(file, exp)
    
    # do TIC
    retention_times, intensities = calcTIC(exp)
    plt.plot(retention_times, intensities)
    plt.title(file + '\nTIC')
    plt.xlabel('time')
    plt.ylabel('intensity')
    plt.tight_layout()
    plt.show()
    
    
def plot_TIC_multifile(file_dict, savefig = False):
    '''
    Example of file_dict
    ====================
    {sampleID: {file_path: "path2mzML",
                color: "green",
                group: "HEU"}
                ...
                }
    '''
    for k,v in file_dict.items():
        
        exp = MSExperiment()
        MzMLFile().load(v['file_path'], exp)
        # do TIC
        retention_times, intensities = calcTIC(exp)
        plt.plot(retention_times, intensities, color= v['color'], linewidth = 0.5,alpha = 0.8)
        
    plt.title('multiple files TIC')
    plt.xlabel('time')
    plt.ylabel('intensity')
    plt.tight_layout()
    
    if savefig == True:
        plt.savefig("figure.png", dpi=600)  # save the figure with higher dpi
    plt.show()