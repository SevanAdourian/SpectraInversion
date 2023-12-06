import numpy as np
import os
import yaml
import sys

def load_yaml_file(file_path):
    with open(file_path, 'r') as file:
        parameter = yaml.safe_load(file)
        return parameter
    
def ensure_directory_exists(path):
    # Check if the directory exists
    if not os.path.exists(path):
        # If not, create the directory
        os.makedirs(path)
        print(f"Directory '{path}' created.")
    else:
        print(f"Directory '{path}' already exists.")

def choose_pow2_n(trace):
    orig_npts = int(trace.stats.npts)
    pow2s = [2**i for i in np.arange(35)] 
    indexes_gt_orignpts = [index for index,value in enumerate(pow2s) if value > orig_npts]
    new_npts = pow2s[indexes_gt_orignpts[3]]

    return new_npts

def detide(tr, order, detiding_ds, plot=False, plot_obspy=False):
        
    tided_trace = tr.copy()
    tr.detrend(type='spline', order=order, dspline=detiding_ds, plot=plot_obspy)

    if plot == True:
        fig, ax = plt.subplots(figsize=(12,4))
        ax = before_after_figure(tided_trace, 'Original', tr, 'Detided', 'Detiding Visualization', ax, label_axes=True)
        plt.savefig('detiding.svg')
        plt.close()

    return tr   

def trim_tr(tr, lag, duration, plot=False):
    start_trim = tr.stats.starttime + (60*60*lag) # seconds
    end_trim = start_trim + (60*60*duration)
    if end_trim > tr.stats.endtime:
                print( 'WARNING: TRIMMING INCOMPATIBLE WITH TRACE LENGTH' )
    untrimmed_trace = tr.copy()
    tr.trim(starttime=start_trim, endtime=end_trim)
    if plot == True:
        fig, ax = plt.subplots(figsize=(12,4))
        ax = before_after_figure(untrimmed_trace, 'Original', tr, 'Trimmed', 'Trimming Visualization', ax, label_axes=True)
        # plt.savefig('trimming.svg')
        plt.show()
        plt.close()
        
    return tr
        
def bandpass(tr, low_corner, high_corner, corners, zerophase, plot=False):
    
    unfiltered_trace = tr.copy()
    tr.filter('Bandpass', freqmin=low_corner, freqmax=high_corner, corners=corners, zerophase=zerophase)
    if plot == True:
        fig, ax = plt.subplots(figsize=(12,4))
        ax = before_after_figure(unfiltered_trace, 'Original', self.trace, 'Filtered', 'Filtering Visualization', ax, label_axes=True)
        plt.savefig('filtering.svg')
        plt.close()
    
    return tr


def create_station_id_string(trace, use_periods=1):
    
    sts = trace.stats
    if use_periods == 1:
        id_string = sts.network +'.'+ sts.station +'.'+ sts.location +'.'+ sts.channel
    else:
        id_string = sts.network +'_'+ sts.station +'_'+ sts.location +'_'+ sts.channel
    return id_string


def get_spectra(st):
    # Select trace
    comp_st = st.select(component=component)
    tr = comp_st[0].copy()
    detiding_ds = int(1/(0.0001*2*tr.stats.delta))

    # Detide
    tr = detide(tr, order, detiding_ds, plot=False, plot_obspy=False)
    tr = trim_tr(tr, lag, duration, plot=False) # Trim
    tr = detrend(tr, 'linear') # Detrend
    tr = taper(tr, max_perc, ttype, max_length, side) # Taper
    tr = bandpass(tr, low_corner, high_corner, corners, zerophase, plot=False) # Filter

    # Plot
    xlims = [low_corner*1000, high_corner*1000]
    trace = tr.copy()
    time_data = trace.data
    n = choose_pow2_n(trace) #131072 # int(trace.stats.npts) # length of the signal [pts] #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    df = trace.stats.sampling_rate # sampling frequency [Hz]
    k = np.arange(n)
    T = n/df
    freq_vec = k/T # two sides frequency range
    freq_vec = freq_vec[range(n//2)]

    Y = np.fft.fft(time_data, n=n)/n
    phase = np.angle(Y)
    Y = Y[range(n//2)]
    Y = Y[:int(len(Y)/10)]
    freq_vec = freq_vec[:int(len(freq_vec)/10)]*1000

    return Y**2, freq_vec


def initialize_plot():
#     fig = plt.figure(figsize=(18, 8))
    fig = plt.figure(figsize=(10, 8))
    plt.ioff()
    gs = fig.add_gridspec(3, 1)
    ax_amp = plt.subplot(gs[1:])
    ax_phase = plt.subplot(gs[0])
    return fig, ax_amp, ax_phase

def add_spectra(st, c, ls, label, return_vecs=False):
    
    # Select trace
    comp_st = st.select(component=component)
    tr = comp_st[0].copy()
    detiding_ds = int(1/(0.0001*2*tr.stats.delta))

    # Detide
    tr = detide(tr, order, detiding_ds, plot=False, plot_obspy=False)
    tr = trim_tr(tr, lag, duration, plot=False) # Trim
    tr = detrend(tr, 'linear') # Detrend
    tr = taper(tr, max_perc, ttype, max_length, side) # Taper
    tr = bandpass(tr, low_corner, high_corner, corners, zerophase, plot=False) # Filter

    # Plot
    xlims = [low_corner*1000, high_corner*1000]
    trace = tr.copy()
    time_data = trace.data
    n = choose_pow2_n(trace) #131072 # int(trace.stats.npts) # length of the signal [pts] #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    df = trace.stats.sampling_rate # sampling frequency [Hz]
    k = np.arange(n)
    T = n/df
    freq_vec = k/T # two sides frequency range
    freq_vec = freq_vec[range(n//2)]

    Y = np.fft.fft(time_data, n=n)/n
    phase = np.angle(Y)
    Y = Y[range(n//2)]
    Y = Y[:int(len(Y)/10)]
    freq_vec = freq_vec[:int(len(freq_vec)/10)]*1000
    
    ax_phase.plot(freq_vec, np.angle(Y), c=c, ls=ls, lw=4.5, alpha=0.9) # Try Y**2 and/or abs
    ax_amp.plot(freq_vec, abs(Y**2), c=c, ls=ls, lw=4.5, label=label, alpha=0.9)
    
    if return_vecs ==True:
        spec_df = pd.DataFrame(list(zip(abs(Y**2),freq_vec)), columns =['pow', 'freq'])
        return spec_df  
    
    
def fcal(f1, f2, dt, tout, mex, qex):
    fn = 0.5 / dt
    
    if fn < f2:
        raise ValueError("f2 is greater than the Nyquist frequency for the time step")
    
    ep = mex / tout
    # print(f"ep = {ep}")
    df = ep / (2.0 * np.pi * qex)
    # print(f"df = {df}")
    # df = ep / qex
    nt = int(np.ceil(1.0 / (df * dt)))
    # print(f"nt = {nt}")
        
    # df = 1.0 / (nt * dt)
    # print(f"df2 = {df}")
    i1 = np.floor(f1 / df)
    print(f"i1 = {i1}")
    f1 = (i1) * df

    i2 = np.ceil(f2 / df)
    print(f"i2 = {i2}")
    f2 = (i2) * df
    
    return f1, f2, df, ep, nt, i1, i2


def hann(t, t11, t12, t21, t22):
    if t11 == 0.0 and t12 == 0.0 and t21 == 0.0 and t22 == 0.0:
        return 0.0

    if t < t11:
        return 0.0
    elif t >= t11 and t < t12:
        hann_value = np.pi * (t - t11) / (t12 - t11)
        hann_value = 0.5 * (1.0 - np.cos(hann_value))
        return hann_value
    elif t >= t12 and t < t21:
        return 1.0
    elif t >= t21 and t < t22:
        hann_value = np.pi * (t22 - t) / (t22 - t21)
        hann_value = 0.5 * (1.0 - np.cos(hann_value))
        return hann_value
    elif t >= t22:
        return 0.0

