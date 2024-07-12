import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import math
import yaml
import glob
import pyasdf

from scipy.signal import periodogram, coherence
from scipy.optimize import fmin

import ProcessSeismicData as PSD

def IBPM_PressureCorrection(pressure, fs):
    # From Lepage et al., 2023. Coefficients could
    # be redefined also using the minimization technique.
    # Here we use the value defined by their paper
    # a = -3.5 # nm.s^{-2}.hPa^{-1}
    # fn = 2.25 # mHz
    # Conversion
    a = -3.5*1e-4 # m.s^{-2}.Pa^{-1}
    fn = 2.25*1e-3 # Hz
    n = len(pressure)

    # plt.plot(pressure)
    # plt.show()
    # Compute pressure fft
    pressure_f = np.fft.rfft(pressure)
    freq_array = np.fft.rfftfreq(n, d=1./fs)

    # plt.plot(freq_array, pressure_f)
    # plt.show()
    alp_f = a*(1 - (freq_array**2/fn**2))

    # Apply correction factor to pressure in frequency domain
    pressure_correction_f = alp_f * pressure_f

    # IFFT to go back to time domain
    pressure_correction_t = np.fft.irfft(pressure_correction_f)
    # plt.plot(pressure_correction_t)
    # plt.show()
    return pressure_correction_t

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <yaml_file>")
        sys.exit(1)

    # Get the YAML file name from the command-line argument
    yaml_file_name = sys.argv[1]
    
    # Load YAML content
    param = PSD.load_yaml_file(yaml_file_name)    
    
    facf = param['taper_frequency_factor']
    fact = param['taper_time_factor']

    f1 = param['low_corner']/1.e3  # low corner
    f2 = param['high_corner']/1.e3 # high corner
    f12 = f1+facf*(f2-f1)
    f21 = f2-facf*(f2-f1)

    # Load asdf files
    print(f"{param['basedir']}Data/*")
    events = glob.glob(f"{param['basedir']}Data/*")

    # Initialize the inventory for event
    inv = []

    print(f"All events are processed starting from {param['start_time_series']}h"
          f"after the origin of the earthquake and for a duration of {param['end_time_series']}h")
    # Begin loop on all events/asdf files
    print(events)
    for _i, event_path in enumerate(events):
        # Print some stuff here
        event = os.path.basename(event_path)
        print(f"Processing event {event}")

        # load asdf
        filename_asdf = f"{event_path}/asdf/{event}.h5"

        # Make a copy of asdf because we will add spectra information to it
        filename_asdf_pro = f"{event_path}/asdf/{event}_pro.h5"
        os.system(f"cp {filename_asdf} {filename_asdf_pro}")
        
        with pyasdf.ASDFDataSet(filename_asdf_pro, compression="gzip-3", mode='r+') as ds:
            for stations in ds.waveforms:
                # Reconstruct inventory
                inv2 = stations.StationXML
                inv += stations.StationXML

                st = stations.raw_waveform
                # Select the waveform data
                for i, tr in enumerate(st):
                    

                    # demean/detrend
                    tr = tr.detrend('demean')
                    tr = tr.detrend('linear')
                    
                    # Remove instrument respones
                    # Remove the response in the frequency domain as this adds noise
                    # tr = tr.remove_response(inventory=inv, output='ACC', water_level=60, plot=False,
                    #                        pre_filt=[0.0001, 0.00011, 0.0051, 0.0052])

                    # Removing detiding because tides out of band
                    # detiding_ds = int(1/(0.0001*2*tr.stats.delta))
                    # tr = PSD.detide(tr, param['order_detiding'], detiding_ds, plot=False, plot_obspy=False)

                    # Doing a second pass for windowing
                    tr = PSD.trim_tr(tr, param['start_time_series'], param['end_time_series'], plot=False) # Trim
                    tr = tr.taper(type=param['taper_type'], max_percentage=fact)
                    tr = tr.detrend('linear')
                    # tr.plot()
                    # IBPM Pressure correction here, if exists

                    # Check if we can do barometric corrections
                    try:
                        pr = stations.pressure[0]
                        print("This station has also pressure data, can do corrections")
                        
                        # demean/detrend
                        pr = pr.detrend('demean')
                        pr = pr.detrend('linear')
                        
                        # Windowing
                        pr = PSD.trim_tr(pr, param['start_time_series'], param['end_time_series'], plot=False) # Trim
                        pr_ts = IBPM_PressureCorrection(pr.data, pr.stats.sampling_rate)

                        # Applying correction to the time series
                        tr_d = tr.data[:-1] - pr_ts
                        
                        ######################
                        # length = 400
                        # tr_t = np.pad(tr.data, (
                        #     int((length * 400 * 24 / 10. - tr.stats.npts) / 2.),
                        #     int((length * 400 * 24 / 10. - tr.stats.npts) / 2.)), 'edge')
                        # NFFT = 2 ** (math.ceil(math.log(tr.stats.npts, 2)))
                        
                        # # Compute fft
                        # f, acc = periodogram(tr_t, fs=tr.stats.sampling_rate, nfft=NFFT,
                        #                      scaling='spectrum')
                        # f, acc = f[1:], acc[1:]
                        
                        # inv_resp = inv2.get_response(tr.id, tr.stats.starttime)
                        # resp, _ = inv_resp.get_evalresp_response(tr.stats.delta, NFFT, 'ACC')
                        # # resp = resp[1:-1]
                        # resp = resp[1:]
                        # acc /= np.abs(resp)
                        # acc_fil= acc
                        
                        # # Convert units to nm/s/s
                        # acc_amp = np.sqrt(acc_fil) * 1.e9
                        # f *= 1000.
                        # acc_win = acc_amp[(f >= param['low_corner']) & (f <= param['high_corner'])]
                        # f = f[(f >= param['low_corner']) & (f <= param['high_corner'])]
                        
                        # spectrum = np.abs(acc_win)
                        # frequencies = f

                        ##########################
                        # length = 400
                        # tr_t = np.pad(tr_t, (
                        #     int((length * 400 * 24 / 10. - tr.stats.npts) / 2.),
                        #     int((length * 400 * 24 / 10. - tr.stats.npts) / 2.)), 'edge')
                        # NFFT = 2 ** (math.ceil(math.log(tr.stats.npts, 2)))
                        
                        # # Compute fft
                        # f, acc = periodogram(tr_t, fs=tr.stats.sampling_rate, nfft=NFFT,
                        #                      scaling='spectrum')
                        # f, acc = f[1:], acc[1:]
                        
                        # inv_resp = inv2.get_response(tr.id, tr.stats.starttime)
                        # resp, _ = inv_resp.get_evalresp_response(tr.stats.delta, 2*NFFT, 'ACC')
                        # resp = resp[1:-1]
                        # acc /= np.abs(resp)
                         
                        # # Convert units to nm/s/s
                        # acc_amp = np.sqrt(acc) * 1.e9
                        # f *= 1000.
                        # acc_win_pc = acc_amp[(f >= param['low_corner']) & (f <= param['high_corner'])]
                        # f = f[(f >= param['low_corner']) & (f <= param['high_corner'])]
                        
                        # # Plot to check
                        # fig = plt.figure(1,figsize=(12,12))
                        # plt.plot(f,np.abs(acc_win), label='Uncorrected')
                        # plt.plot(f,np.abs(acc_win_pc), label='Corrected')
                        # plt.xlabel('Frequency (mHz)')
                        # plt.ylabel('Acceleration (nm/s/s)')
                        # plt.title((tr.id).replace('.',' '))
                        # plt.legend()
                        # # plt.savefig(f"{param['basedir']}/Figures/Spectra/{event}/{tr.id}.png",format='PNG',dpi=400)
                        # # plt.close()
                        # plt.show()
                        ##################################
                    
                    except pyasdf.WaveformNotInFileException as e:
                        print("No barometer was found for this station, no correction is applied")
                        tr_d = tr.data
                        
                    # Now back to doing usual processing of the corrected (or not) trace.
                    length = 400
                    tr_d = np.pad(tr_d, (
                        int((length * 400 * 24 / 10. - tr.stats.npts) / 2.),
                        int((length * 400 * 24 / 10. - tr.stats.npts) / 2.)), 'edge')
                    NFFT = 2 ** (math.ceil(math.log(tr.stats.npts, 2)))
                    
                    # Compute fft
                    f, acc = periodogram(tr_d, fs=tr.stats.sampling_rate, nfft=NFFT,
                                         scaling='spectrum')
                    f, acc = f[1:], acc[1:]
                    
                    inv_resp = inv2.get_response(tr.id, tr.stats.starttime)
                    resp, _ = inv_resp.get_evalresp_response(tr.stats.delta, NFFT, 'ACC')
                    resp = resp[1:]
                    acc /= np.abs(resp)
                    acc_fil= acc
                
                    # Convert units to nm/s/s
                    acc_amp = np.sqrt(acc_fil) # * 1.e9
                    # f *= 1000.
                    f1f = param['low_corner']/1.e3
                    f2f = param['high_corner']/1.e3
                    acc_win = acc_amp[(f >= f1f) & (f <= f2f)]
                    f = f[(f >= f1f) & (f <= f2f)]

                    spectrum = np.abs(acc_win)
                    frequencies = f
                    
                    # # Plot to check
                    # fig = plt.figure(1,figsize=(12,12))
                    # plt.plot(f,spectrum, label='Uncorrected')
                    # plt.xlabel('Frequency (mHz)')
                    # plt.ylabel('Acceleration (nm/s/s)')
                    # plt.title((tr.id).replace('.',' '))
                    # plt.legend()
                    # # plt.savefig(f"{param['basedir']}/Figures/Spectra/{event}/{tr.id}.png",format='PNG',dpi=400)
                    # plt.close()
                    # plt.show()
                try:
                    # Write spectra as auxiliary data
                    print(f"Writing spectrum information for station")

                    datatype = "ProcessedSpectra"
                    datapath = tr.id
                    dataparams = {
                        "start_freq": f[0],
                        "nfreq": len(f),
                        "dfreq": (f[2]-f[1])}

                    ds.add_auxiliary_data(data=spectrum, data_type=datatype,
                                      path=datapath, parameters=dataparams)
                except:
                    pass
                
            print(ds)
