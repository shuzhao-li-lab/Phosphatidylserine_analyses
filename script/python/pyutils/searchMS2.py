import numpy as np
from matplotlib import pyplot as plt
import pymzml
import os
import pandas as pd
import sys

# https://github.com/shuzhao-li/asari/blob/284d49db05bc95377072a865663550f59340aca3/asari/tools/plot.py
def get_potential_precursor_from_file(infile, 
                                     min_scan_number, 
                                     max_scan_number, 
                                     min_mz, 
                                     max_mz,
                                     ms_level= 2):
    '''
    input
    -----
    infile: mzML file as input
    return 
    ------
    list, [(scan_number, mz, intensity value), ...]
    Note
    ----
    The precursor charge no matter mode is always positive integer
    '''
    alldata = []
    exp = pymzml.run.Reader(infile)
    ii = 0   # scan_number starts with 0
    for spec in exp:
        if min_scan_number < ii < max_scan_number: # select the retention time region to reduce spectra of interest
            if spec.ms_level == ms_level: # select ms2 level
                precursor_dict = spec.selected_precursors[0] # maybe because it can be multiplexed so the precursor is a list
                if (min_mz < precursor_dict['mz'] < max_mz):
                    alldata.append(spec)
                
        ii += 1
    return alldata



def get_potential_precursor_from_exp_filtbyScan(exp, 
                                     min_scan_number, 
                                     max_scan_number, 
                                     min_mz, 
                                     max_mz,
                                     charge_state = 1,
                                     ms_level= 2):
    '''
    input
    -----
    infile: exp as input
    return 
    ------
    list, [(scan_number, mz, intensity value), ...]
    Note
    ----
    The precursor charge no matter mode is always positive integer
    '''
    ii = 0   # scan_number starts with 0
    matched_specs = []
    for spec in exp:
        if min_scan_number < ii < max_scan_number: # select the retention time region to reduce spectra of interest
            if spec.ms_level == ms_level: # select ms2 level
                precursor_dict = spec.selected_precursors[0] # maybe because it can be multiplexed so the precursor is a list
                if min_mz < precursor_dict['mz'] < max_mz:
                    matched_specs.append(spec)
                
        ii += 1
    return matched_specs




def get_potential_precursor_from_exp_filtbyRt(exp, 
                                     min_rt_sec, 
                                     max_rt_sec, 
                                     min_mz, 
                                     max_mz,
                                     charge_state = 1,
                                     ms_level= 2):
    '''
    input
    -----
    infile: exp as input
    return 
    ------
    list, [(scan_number, mz, intensity value), ...]
    Note
    ----
    The precursor charge no matter mode is always positive integer
    '''
    ii = 0   # scan_number starts with 0
    matched_specs = []
    for spec in exp:
        rt_sec = round(spec.scan_time[0]*60,3)
        if min_rt_sec < rt_sec < max_rt_sec: # select the retention time region to reduce spectra of interest
            if spec.ms_level == ms_level: # select ms2 level
                precursor_dict = spec.selected_precursors[0] # maybe because it can be multiplexed so the precursor is a list
                if min_mz < precursor_dict['mz'] < max_mz:
                    matched_specs.append(spec)
                
        ii += 1
    return matched_specs



import matplotlib.pyplot as plt

def plot_spectra(spectra, 
                 save_figure = False,
                 output_path = "./",
                 label = ""):
    for spectrum in spectra:
        
        mz_values = spectrum.mz
        intensity_values = spectrum.i
        
        intensity_values = intensity_values/max(intensity_values)
        
        zipped_mz_int = [(x,y) for x,y in zip(mz_values,intensity_values) if y > 0.05]
        mz_values, intensity_values = zip(*zipped_mz_int)
        
        precursor_mz = round(spectrum.selected_precursors[0]['mz'], 4)

        # Create the plot figure with more height
        fig, ax = plt.subplots(figsize=(8, 6))

        # Create the bar plot
        ax.bar(mz_values, intensity_values, width=2)

        # Add the annotation text to the top of each bar
        for x, y in zip(mz_values, intensity_values):
            ax.text(x, y, f"{round(x,2)}", ha='left', va='bottom', rotation=45, fontsize=16)

        # Set the plot title and axis labels
        ax.set_title(f"MS2 scan of {round(precursor_mz,4)}_{round(spectrum.scan_time[0]*60,1)}", fontsize = 20)
        ax.set_xlabel("m/z")
        ax.set_ylabel("Intensity")

        # Set the x-axis limits to 1-1000
        ax.set_xlim(1, 900)

        # Adjust the y-axis limits to leave some space at the top
        ax.set_ylim(top=ax.get_ylim()[1]*1.2)

        # Save the plot as PDF
        if save_figure == True:
            try:
                os.mkdir(output_path)
            except:
                None
            fig.savefig(os.path.join(output_path,f"{label}_ms2_scan_{round(precursor_mz,4)}_{round(spectrum.scan_time[0]*60,1)}.pdf"))

        # Show the plot (optional)
        plt.show()
        
        
        
def cal_ppm(mz1,
            mz2):
    return abs((mz1-mz2))*1000000/mz2

# this function only works for situation where you look at charge state = 1
def search_NL(spectra,
              NL_mz = 87.03124,
              ppm_prec = 100,
              ppm = 40):
    res_data = []

    for spec in spectra:
        if cal_ppm(spec.selected_precursors[0]['mz'],max(spec.mz)) > ppm_prec:
            sel_prec_mz = spec.selected_precursors[0]['mz'] # this m/z will not be exactly the precursor m/z    
        else:
            sel_prec_mz = max(spec.mz)
        
        sel_mz = [mz for mz in spec.mz if abs(mz - sel_prec_mz) < (np.ceil(NL_mz)+1)] # +1 is to increase a little bit the range
    
        for mz in sel_mz:
            calc_ppm = cal_ppm((sel_prec_mz - mz),NL_mz)
            if calc_ppm < ppm:
                print(calc_ppm)
                res_data.append(spec)
                break
    return res_data