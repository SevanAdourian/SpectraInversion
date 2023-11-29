import numpy as np
import matplotlib.pyplot as plt
import scipy
import sys
import re
import os
import math
import yaml
import glob
import pyasdf

from scipy.signal import periodogram, coherence
from scipy.optimize import fmin

import ProcessSeismicData as PSD

def load_data(asdf_file, tr_id):
    with pyasdf.ASDFDataSet(filename_asdf_pro, compression="gzip-3", mode='r+') as ds:
        spec_ds = ds.auxiliary_data.ProcessedSpectra[tr_id]
        sfreq = spec_ds.parameters["start_freq"]
        nfreq = spec_ds.parameters["nfreq"]
        dfreq = spec_ds.parameters["dfreq"]
        efreq = sfreq+(nfreq-1)*dfreq
        freq_array = np.arange(sfreq, efreq+dfreq, dfreq)

        return freq_array, spec_ds

def load_synthetic(conf, event, stn):

    return freq_array, synth_spec

def is_synthetic_computed(conf, event, stn):

    return synth_ok

def process_synthetic(conf, spec):

    return processed_synth_spec

def compute_misfit(processed_synth_spec, data_spec_ds):

    return adj_src

def process_adjoint_sources(conf, l2_misfit):

    return pro_adj_src


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <yaml_file>")
        sys.exit(1)

    # Get the YAML file name from the command-line argument
    yaml_file_name = sys.argv[1]
    
    # Load YAML content
    param = PSD.load_yaml_file(yaml_file_name)    

    # Generate sample data points
    f1_inp = param['minimum_frequency']
    f2_inp = param['maximum_frequency']/1.e3
    t1_inp = param['start_time_series'] * 3600.
    t2_inp = param['end_time_series'] * 3600.
    dt = 1. / (2.*f2_inp)
    
    facf = param['taper_frequency_factor']
    fact = param['taper_time_factor']

    f1 = param['low_corner']/1.e3  # low corner
    f2 = param['high_corner']/1.e3 # high corner
    
    # FCAL usage, can be unnormalized, we don't care here.
    # f_norm = np.sqrt(6.6723e-11*np.pi*5515.)
    # t_norm = 1./f_norm
    
    mex = 5
    qex = 4
    f1, f2, df, ep, nt, i1, i2 = PSD.fcal(f1_inp/f_norm, f2_inp/f_norm, dt/t_norm, t2_inp/t_norm, mex, qex)
    
    print(".... Frequency parameters:")
    print(".... ... df (mHz) =", df)
    print(".... ... nt       =", nt)
    print(".... ... f1 (mHz) =", f1)
    print(".... ... f2 (mHz) =", f2)
    
    # Computation of all parameters
    f12 = f1+facf*(f2-f1)
    f21 = f2-facf*(f2-f1)
    
    # t1 = t1_inp/t_norm
    # t2 = t2_inp/t_norm
    # t12 = t1+fact*(t2-t1)
    # t21 = t2-fact*(t2-t1)
    # dt = dt/t_norm
    # dw = df * 2. * np.pi
    
    ##################
    # ntall = nt
    # nt = int(i2 - i1)
    # x = np.linspace(f1, f2, nt)
    # x = np.concatenate((x,-np.flipud(x)[:-1]))
    
    # Load asdf file
    events    = glob.glob(f"{param['basedir']}Data/*")

    # Initialize the inventory for event
    inv = []
    
    for _i, event_path in enumerate(events):
        # Print some stuff here
        event = os.path.basename(event_path)
        print(event)
        print(f"{event_path}/asdf/{event}.h5")

        # load asdf
        filename_asdf = f"{event_path}/asdf/{event}.h5"
        with pyasdf.ASDFDataSet(filename_asdf_pro, compression="gzip-3", mode='r+') as ds:
            for stations in ds.waveforms:
                # Reconstruct inventory
                inv += stations.StationXML

                st = stations.raw_waveform
                # Select the waveform data
                for i, tr in enumerate(st):
                    NFFT = 2 ** (math.ceil(math.log(tr.stats.npts, 2)))

                    # demean/detrend
                    tr = tr.detrend('demean')
                    tr = tr.detrend('linear')
                    
                    # Remove instrument respones
                    tr = tr.remove_response(inventory=inv, output='ACC', water_level=60, plot=False,
                                            pre_filt=[0.0001, 0.00011, 0.0051, 0.0052])

                    # Removing detiding because tides out of band
                    # detiding_ds = int(1/(0.0001*2*tr.stats.delta))
                    # tr = PSD.detide(tr, param['order_detiding'], detiding_ds, plot=False, plot_obspy=False)

                    # Doing a second pass for windowing
                    tr = PSD.trim_tr(tr, param['start_time_series'], param['end_time_series'], plot=False) # Trim
                    tr = tr.taper(type=param['taper_type'], max_percentage=fact)
                    tr = tr.detrend('linear')
                                        
                    # Compute fft
                    f, acc = periodogram(tr.data, fs=tr.stats.sampling_rate, nfft=NFFT,
                                       scaling='spectrum')
                    f, acc = f[1:], acc[1:]

                    # Hann windowing in frequency
                    acc_fil = np.zeros(len(acc), dtype=float)
                    for i in range(0,len(f)):
                        hann_coeff = PSD.hann(f[i], f1, f12, f21, f2)
                        acc_fil[i] = hann_coeff*acc[i]

                    # Convert units to nm/s/s
                    acc_amp = np.sqrt(acc_fil) * 1.e9
                    f *= 1000.
                    acc_win = acc_amp[(f >= param['low_corner']) & (f <= param['high_corner'])]
                    f = f[(f >= param['low_corner']) & (f <= param['high_corner'])]
                    
                    # Check if we can do barometric corrections
                    try:
                        pr = stations.pressure[0]
                        print("This station has also pressure data, can do corrections")
                        
                        # demean/detrend
                        pr = pr.detrend('demean')
                        pr = pr.detrend('linear')
                        
                        # Windowing
                        pr = PSD.trim_tr(pr, param['start_time_series'], param['end_time_series'], plot=False) # Trim
                        pr = pr.taper(type=param['taper_type'], max_percentage=fact)
                        pr = pr.detrend('linear')
                        
                        # Compute FFT
                        f, pres = periodogram(pr.data, fs=pr.stats.sampling_rate, nfft=NFFT,
                                           scaling='spectrum')
                        f, pres = f[1:], pres[1:]

                        # Hann windowing in the frequency domain
                        pres_fil = np.zeros(len(acc), dtype=float)
                        for i in range(0,len(f)):
                            hann_coeff = PSD.hann(f[i], f1, f12, f21, f2)
                            pres_fil[i] = hann_coeff*pres[i]

                        pres_amp = np.sqrt(pres_fil)
                        f *= 1000. # We back in mHz
                        pres_win = pres_amp[(f >= param['low_corner']) & (f <= param['high_corner'])]
                        f = f[(f >= param['low_corner']) & (f <= param['high_corner'])]

