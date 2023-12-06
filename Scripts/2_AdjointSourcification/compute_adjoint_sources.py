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
from scipy.interpolate import interp1d

import ProcessSeismicData as PSD

class SynthSpectrum:
    def __init__(self, frequency, real_part, imaginary_part, absolute_value):
        self.frequency = frequency
        self.real = real_part
        self.imag = imaginary_part
        self.abs_val = absolute_value

class SynthSpectrumCollection:
    def __init__(self):
        self.records = {}

    def add_record(self, station, record):
        if station not in self.records:
            self.records[station] = []
        self.records[station].append(record)

def load_data(spec_ds):
    sfreq = spec_ds.parameters["start_freq"]
    nfreq = spec_ds.parameters["nfreq"]
    dfreq = spec_ds.parameters["dfreq"]
    efreq = sfreq+(nfreq-1)*dfreq
    freq_array = np.arange(sfreq, efreq+dfreq, dfreq)
    data_spec = np.array(spec_ds.data)
        
    return freq_array, data_spec

def load_synthetic(conf, event, stn):
    
    file_path = f"{conf['basedir']}Scripts/1_SyntheticComputation/results/{event}/{conf['starting_model']}/{stn}.{comp}.dat"
    try:
        frequency, real, imag, abs_val = [], [], [], []
        
        with open(file_path, 'r') as file:
            for line in file:
                if line.strip():  # Skip empty lines
                    columns = line.split()
                    frequency.append(float(columns[0]))
                    real.append(float(columns[1]))
                    imag.append(float(columns[2]))
                    abs_val.append(float(columns[3]))

                    synth_record = SynthSpectrum(np.array(frequency), np.array(real),
                                                 np.array(imag), np.array(abs_val))

    except Exception as e:
        print(f"An error occurred: {e}")
        return None
    
    return synth_record

def is_synthetic_computed(conf, event, stn, comp):

    path_to_look = f"{conf['basedir']}Scripts/1_SyntheticComputation/results/{event}/{conf['starting_model']}/{stn}.{comp}.dat"
    print(path_to_look)
    if os.path.exists(path_to_look):
        synth_ok = True
        print(f"Found synthetic for event {event} and station {stn}")
        print(f"Continuing the workflow")
    else:
        synth_ok = False
        print(f"The synthetic for event {event} and station {stn} does not exist")
        print(f"Skipping...")
    return synth_ok


def process_synthetic(conf, synth_spec, dt, dt_data,
                      t1, t2, t12, t21, ep):

    # Do ifft to get seismogram in time domain (real part)
    ns = len(synth_spec.real)
    
    synth_spec_complex = np.zeros(ns, dtype=complex)
    synth_spec_complex.real = synth_spec.real
    synth_spec_complex.imag = synth_spec.imag

    NFFT = 2 ** (math.ceil(math.log(ns, 2)))
    synth_seismo = np.fft.irfft(synth_spec_complex, n=NFFT, norm='backward')
    
    # Compute time support array
    time_support = np.arange(0, conf['end_time_series_synth']*3600, dt)
    data_time_array = np.arange(0, time_support[-1], dt_data)
    
    # Interpolate on the same time support as the data
    synth_seismo = np.real(synth_seismo[:len(time_support)])
    interpolator = interp1d(time_support, synth_seismo, kind='cubic', fill_value='extrapolate')
    synth_seismo_int = interpolator(data_time_array)

    # Correct epsilon and time domain
    synth_seismo_windowed = np.zeros(len(data_time_array))
    for i,t_value in np.ndenumerate(data_time_array):
        hann_coeff = PSD.hann(t_value, t1, t12, t21, t2)/np.pi
        synth_seismo_windowed[i] = hann_coeff*synth_seismo_int[i]*np.exp(ep*t_value)/(2*np.pi)

    # Back to spectral domain (fft)
    # Pad to make it look good
    NFFT = 2 ** (math.ceil(math.log(len(synth_seismo_windowed), 2))+3)
    processed_synth_spec = np.fft.fft(synth_seismo_windowed, n=NFFT, norm='backward')
    f = np.fft.fftfreq(NFFT, d=dt_data)

    # Window in frequency domain
    processed_synth_spec = processed_synth_spec[(f >= param['low_corner']/1e3) & (f <= (param['high_corner'])/1e3)]
    freq_support_pro = f[(f >= param['low_corner']/1e3) & (f <= (param['high_corner'])/1e3)]
    
    return freq_support_pro,processed_synth_spec


# def compute_misfit(processed_synth_spec, data_spec_ds):

#     interpolator = interp1d(data_freq_array, data_spec, kind='linear', fill_value='extrapolate')
#     interpolated_data_spec = interpolator(processed_synth_spec)

#     adj_src = np.zeros(len(processed_synth_spec))
#     adj_src = processed_synth_spec.real - data_spec_ds.real 
#     return adj_src

# def process_adjoint_sources(conf, adj_src):
    
#     return pro_adj_src


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
    t2_inp = param['end_time_series_data'] * 3600.
    dt = 1. / (2.*f2_inp)
    
    facf = param['taper_frequency_factor']
    fact = param['taper_time_factor']

    f1 = param['low_corner']/1.e3  # low corner
    f2 = param['high_corner']/1.e3 # high corner

    # Useful for rescaling synthetic - due to technical considerations of IDSM
    mex = 5
    ep = mex/(param['end_time_series_synth'] * 3600.)
        
    # Computation of all parameters
    f12 = f1+facf*(f2-f1)
    f21 = f2-facf*(f2-f1)
    t12 = t1_inp+fact*(t2_inp-t1_inp)
    t21 = t2_inp-fact*(t2_inp-t1_inp)
        
    # Load asdf file
    events    = glob.glob(f"{param['basedir']}Data/*")

    for _i, event_path in enumerate(events):
        # Print some stuff here
        event = os.path.basename(event_path)
        print(event)

        # load asdf
        filename_asdf = f"{event_path}/asdf/{event}_pro.h5"
        with pyasdf.ASDFDataSet(filename_asdf, compression="gzip-3", mode='r+') as ds:
            for spec_ds in ds.auxiliary_data.ProcessedSpectra:
                # Get station id
                stn_full = spec_ds.path
                ntw  = stn_full.split('.')[0]
                stn  = stn_full.split('.')[1]
                comp = stn_full.split('.')[3][-1]

                if stn=='KIEV':
                    # Check if we have computed synthetic for this station/event pair
                    # for now only deals with Z
                    is_synth = is_synthetic_computed(param, event, stn, comp)
                    
                    if is_synth:
                        # We have computed synthetic, let's compute adjoint sources
                        print(f"Synthetics were computed for this event/station pair")
                        print(f"Loading synthetic spectra...")
                        synth_spectrum = load_synthetic(param, event, stn)
                        
                        # Now let's load data
                        freq_data, data_spec = load_data(spec_ds)

                        # Process synthetic spectra the same way as data to be consistent
                        freq_pro, synth_spectrum_pro = process_synthetic(param, synth_spectrum, dt, 1,
                                                                           t1_inp, t2_inp, t21, t12, ep)
                        
                        plt.figure()
                        plt.plot(freq_pro*1e3, np.abs(synth_spectrum_pro))
                        plt.plot(freq_data*1e3,data_spec)
                        plt.show()
                        
                        # Compute L2-misfit between synth and data
                        # adj_src = compute_misfit(synth_spectrum_pro, freq_array, data_spec)
                        
                    else:
                        # No synthetic was computed, either data is not good, or the specific path not needed
                        print(f"No synthetic was found, skipping station {stn}")

                else:
                    pass
                # # Select the spectra data
                # for i, tr in enumerate(st):
                #     NFFT = 2 ** (math.ceil(math.log(tr.stats.npts, 2)))

                #     # demean/detrend
                #     tr = tr.detrend('demean')
                #     tr = tr.detrend('linear')
                    
                #     # Remove instrument respones
                #     tr = tr.remove_response(inventory=inv, output='ACC', water_level=60, plot=False,
                #                             pre_filt=[0.0001, 0.00011, 0.0051, 0.0052])

                #     # Removing detiding because tides out of band
                #     # detiding_ds = int(1/(0.0001*2*tr.stats.delta))
                #     # tr = PSD.detide(tr, param['order_detiding'], detiding_ds, plot=False, plot_obspy=False)

                #     # Doing a second pass for windowing
                #     tr = PSD.trim_tr(tr, param['start_time_series'], param['end_time_series'], plot=False) # Trim
                #     tr = tr.taper(type=param['taper_type'], max_percentage=fact)
                #     tr = tr.detrend('linear')
                                        
                #     # Compute fft
                #     f, acc = periodogram(tr.data, fs=tr.stats.sampling_rate, nfft=NFFT,
                #                        scaling='spectrum')
                #     f, acc = f[1:], acc[1:]

                #     # Hann windowing in frequency
                #     acc_fil = np.zeros(len(acc), dtype=float)
                #     for i in range(0,len(f)):
                #         hann_coeff = PSD.hann(f[i], f1, f12, f21, f2)
                #         acc_fil[i] = hann_coeff*acc[i]

                #     # Convert units to nm/s/s
                #     acc_amp = np.sqrt(acc_fil) * 1.e9
                #     f *= 1000.
                #     acc_win = acc_amp[(f >= param['low_corner']) & (f <= param['high_corner'])]
                #     f = f[(f >= param['low_corner']) & (f <= param['high_corner'])]
                    
                #     # Check if we can do barometric corrections
                #     try:
                #         pr = stations.pressure[0]
                #         print("This station has also pressure data, can do corrections")
                        
                #         # demean/detrend
                #         pr = pr.detrend('demean')
                #         pr = pr.detrend('linear')
                        
                #         # Windowing
                #         pr = PSD.trim_tr(pr, param['start_time_series'], param['end_time_series'], plot=False) # Trim
                #         pr = pr.taper(type=param['taper_type'], max_percentage=fact)
                #         pr = pr.detrend('linear')
                        
                #         # Compute FFT
                #         f, pres = periodogram(pr.data, fs=pr.stats.sampling_rate, nfft=NFFT,
                #                            scaling='spectrum')
                #         f, pres = f[1:], pres[1:]

                #         # Hann windowing in the frequency domain
                #         pres_fil = np.zeros(len(acc), dtype=float)
                #         for i in range(0,len(f)):
                #             hann_coeff = PSD.hann(f[i], f1, f12, f21, f2)
                #             pres_fil[i] = hann_coeff*pres[i]

                #         pres_amp = np.sqrt(pres_fil)
                #         f *= 1000. # We back in mHz
                #         pres_win = pres_amp[(f >= param['low_corner']) & (f <= param['high_corner'])]
                #         f = f[(f >= param['low_corner']) & (f <= param['high_corner'])]

