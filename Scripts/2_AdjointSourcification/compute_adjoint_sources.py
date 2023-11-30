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

def process_synthetic(conf, freq_support, spec, dt, f1, f2, f12, f21, t1, t2, t12, t21, ep):
    # From the IDSM computation, process the synthetic seismogram in the
    # same way as data

    nt = len(spec)
    # Compute inverse fourier transform to get seismogram in velocity
    seismogram = np.real(np.fft.ifft(spec))
    time_support = np.arange(0, nt) * dt

    windowed_seismogram = np.zeros(nt,dtype=complex)

    for i,t_value in np.ndenumerate(time_support):
        hann_coeff = hann(t_value, t1, t12, t21, t2)
        windowed_seismogram[i] = hann_coeff*seismogram[i]*np.exp(ep*t_value)

    processed_spec = np.fft.fft(windowed_seismogram)
    
    return processed_spec

def compute_misfit(synth_freq_array, processed_synth_spec, data_freq_array, data_spec_ds):
    # We interpolate the spectra on the same frequency support as the synthetic
    # The adjoint source should be sampled the same way in order to multiply it directly
    # for the adjoint computation

    data_spec = np.array(data_spec_ds.data)
    interpolator = interp1d(data_freq_array, data_spec, kind='linear', fill_value='extrapolate')
    interpolated_data_spec = interpolator(synth_freq_array)

    adj_src = np.zeros(
    adj_src = processed_synth_spec - interpolated_data_spec
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
    f1_inp = param['minimum_frequency']/1.e3
    f2_inp = param['maximum_frequency']/1.e3
    t1_inp = param['start_time_series'] * 3600.
    t2_inp = param['end_time_series'] * 3600.
    dt = 1. / (2.*f2_inp)
    
    facf = param['taper_frequency_factor']
    fact = param['taper_time_factor']

    f1 = param['low_corner']/1.e3  # low corner
    f2 = param['high_corner']/1.e3 # high corner

    # Computation of all parameters
    f12 = f1+facf*(f2-f1)
    f21 = f2-facf*(f2-f1)
    t12 = t1+fact*(t2-t1)
    t21 = t2-fact*(t2-t1)
        
    ##################    
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

