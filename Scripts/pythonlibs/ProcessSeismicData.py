
import numpy as np
import os
import yaml
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb

# Defining global useful constants for normalization purposes
rho_av  = 5515.0
r_earth = 6371000.
BIG_G   = 6.6723e-11 

class Parameters:
    def __init__(self, f1, f2, df, ep, nt, ntall, dt, f12, f21, t1, t2, t12, t21, fscale, tscale, anorm):
        self.f1 = f1
        self.f2 = f2
        self.df = df
        self.ep = ep
        self.nt = nt
        self.ntall = ntall
        self.dt = dt
        self.f12 = f12
        self.f21 = f21
        self.t1 = t1
        self.t2 = t2
        self.t12 = t12
        self.t21 = t21
        self.fscale = fscale
        self.tscale = tscale
        self.anorm = anorm

    def __repr__(self):
        return f"Parameters({self.__dict__})"

class Config:
    def __init__(self, file_path):
        self.file_path = file_path
        self.load_parameters()

    def load_parameters(self):
        with open(self.file_path, 'r') as file:
            for line in file:
                # Skip comments and empty lines
                line = line.strip()
                if line.startswith('#') or not line:
                    continue

                # Split line into key and value
                if '=' in line:
                    key, value = line.split('=', 1)
                    key = key.strip()
                    value = value.strip()

                    # Convert boolean-like strings to actual booleans
                    if value == '.TRUE.':
                        value = True
                    elif value == '.FALSE.':
                        value = False
                    else:
                        # Try to convert to float or int if possible
                        try:
                            if '.' in value:
                                value = float(value)
                            else:
                                value = int(value)
                        except ValueError:
                            pass  # Keep as string if conversion fails

                    # Set each key-value pair as an attribute
                    setattr(self, key, value)

        return

    def setup_run(self, conf_py):
        try:
            # conf_py = load_yaml_file(filepath_yaml)
            # Set additional parameters
            mex = 5
            qex = 4
            fscale = np.sqrt(np.pi * BIG_G * rho_av)
            tscale = 1. / fscale
            anorm = r_earth / (tscale ** 2)

            # Non-dimensionalize everything
            f1_inp = self.MINIMUM_FREQUENCY / 1000. / fscale
            f2_inp = self.MAXIMUM_FREQUENCY / 1000. / fscale
            dt_inp = 1. / (2 * f2_inp)
            t2_inp = self.LENGTH_TIME_SERIES_IN_HOURS * 3600 / tscale  # Converting hours to seconds

            # Call fcal to compute the relevant forward run parameters
            f1, f2, df, ep, nt, i1, i2 = self.fcal(f1_inp, f2_inp, dt_inp, t2_inp, mex, qex)

            # Set calculated attributes
            window_frequency_parameter = conf_py['taper_frequency_factor']
            f12 = f1 + window_frequency_parameter * (f2 - f1)
            f21 = f2 - window_frequency_parameter * (f2 - f1)

            window_time_parameter = conf_py['taper_time_factor']
            t1 = conf_py['start_time_series']*(3600. / tscale)
            t2 = conf_py['end_time_series']*(3600. / tscale)
            t12 = t1 + window_time_parameter * (t2 - t1)
            t21 = t2 - window_time_parameter * (t2 - t1)
            ntall = 2*nt - 1
            
            return Parameters(f1, f2, df, ep, nt, ntall, dt_inp, f12, f21, t1, t2, t12, t21, fscale, tscale, anorm)
        
        except AttributeError as e:
            print(f"Missing attribute: {e}")
            raise

        return

    def fcal(self, f1, f2, dt, tout, mex, qex):
        # All must be normalized here
        fn = 0.5 / dt

        if fn < f2:
            raise ValueError("f2 is greater than the Nyquist frequency for the time step")

        ep = mex / tout
        df = ep / (2.0 * np.pi * qex)
        
        i1 = np.floor(f1 / df)
        f1 = i1 * df

        i2 = np.ceil(f2 / df)
        f2 = i2 * df

        nt = int((i2-i1)+1)

        return f1, f2, df, ep, nt, i1, i2

    # def fcal(self, f1, f2, dt, tout, mex, qex):
    #     # All must be normalized here
    #     fn = 0.5 / dt

    #     if fn < f2:
    #         raise ValueError("f2 is greater than the Nyquist frequency for the time step")

    #     ep = mex / tout
    #     df = ep / (2.0 * np.pi * qex)

    #     # old fcal
    #     nt = int(np.floor(1.0 / (df * dt)))
    #     ne = int(np.log(float(nt))/np.log(2.0)+1)
    #     nt = int(2**ne)
    #     df = 1/(nt*dt)
        
    #     i1 = np.floor(f1 / df)
    #     f1 = i1 * df

    #     i2 = np.ceil(f2 / df)
    #     f2 = i2 * df

    #     nt = int((i2-i1)+1)

    #     pdb.set_trace()
    #     return f1, f2, df, ep, nt, i1, i2

    def __repr__(self):
        return f"Config({self.__dict__})"



class Spectrum:
    def __init__(self, frequency=None, real_part=None, imaginary_part=None, absolute_value=None):
        self.frequency = frequency if frequency is not None else []
        self.realpart = real_part if real_part is not None else []
        self.imagpart = imaginary_part if imaginary_part is not None else []
        self.abs_val = absolute_value if absolute_value is not None else []


    def load_synthetic(self, filepath, param):
        if not filepath:
            raise ValueError("Filepath must be provided.")
        try:
            frequency, real, imag, abs_val = [], [], [], []

            with open(filepath, 'r') as file:
                for line in file:
                    if line.strip():  # Skip empty lines
                        columns = line.split()
                        frequency.append(float(columns[0]) / 1000. / param.fscale)
                        repa = float(columns[1])
                        impa = float(columns[2])
                        # real.append(float(columns[1]))
                        # imag.append(float(columns[2]))
                        # abs_val.append(np.abs(complex(float(columns[1]), float(columns[2]))))
                        real.append(repa)
                        imag.append(impa)
                        abs_val.append(np.abs(complex(repa, impa)))
                        
            self.frequency = np.array(frequency)
            self.realpart = np.array(real)
            self.imagpart = np.array(imag)
            self.abs_val = np.array(abs_val)

        except Exception as e:
            print(f"An error occurred: {e}")

        return

    def taper_and_ifft(self, params, plotting=False, saveplot=False):

        fun = np.zeros(params.nt,dtype = complex)
        fun.real = self.realpart
        fun.imag = self.imagpart
        
        # taper in frequency
        ftapered_fun = np.zeros(params.nt,dtype=complex)
        for i in range(len(ftapered_fun)):
            hann_coeff = hann(self.frequency[i], params.f1, params.f12,
                                  params.f21, params.f2)
            ftapered_fun[i] = hann_coeff*fun[i]*params.anorm
            
        if (plotting == True):
            plt.figure()
            plt.plot(self.frequency*params.fscale*1000., np.abs(fun), label='untapered')
            plt.plot(self.frequency*params.fscale*1000., np.abs(ftapered_fun), label='tapered')
            if (saveplot == True):
                plt.savefig('tapered_spectra.png', dpi=400)
                
            plt.show()

        # constructing negative frequencies
        ftapered_fun_cplx = ftapered_fun.astype(complex)
        ftapered_fun_with_neg = np.zeros(params.ntall, dtype=complex)
        
        ftapered_fun_with_neg[:params.nt] = ftapered_fun_cplx
        ftapered_fun_with_neg[params.nt:] = np.flipud(np.conjugate(ftapered_fun_cplx))[:-1]
        
        # Calculate the inverse Fourier transform
        fun_ifft = np.fft.ifft(ftapered_fun_with_neg, norm='forward') * (1./(params.dt*params.ntall))
        
        # Undoing the exponential decay necessary during the synthetic computation
        # Doing it here because it is only specific to the synthetics
        timeseries = TimeSeries()
        
        # Re-normalization of time and data
        timeseries.times = np.arange(0, (params.ntall)*params.dt, params.dt)
        timeseries.data = np.zeros(params.ntall,dtype=complex)
        
        for i,t_value in np.ndenumerate(timeseries.times):
            timeseries.data.real[i] = fun_ifft.real[i]*np.exp(params.ep*t_value)
            
        return timeseries


class TimeSeries:
    def __init__(self):
        self.times = None
        self.data = None


    def taper_and_fft(self, params, plotting=False, saveplot=False):
        ts_tapered_fun = np.zeros(params.ntall,dtype=complex)
        
        for i,t_value in np.ndenumerate(self.times):
            hann_coeff = hann(t_value, params.t1, params.t12, params.t21, params.t2)
            ts_tapered_fun[i] = hann_coeff*self.data[i]
            
        if (plotting == True):
            plt.figure()
            plt.plot(self.times*params.tscale/3600., ts_tapered_fun, label='tapered')
            if (saveplot == True):
                plt.savefig('tapered_timeseries.png', dpi=400)
            
            plt.show()

        tapered_spectrum = np.fft.fft(ts_tapered_fun.real, norm='backward') * params.dt
        frequency_array = np.fft.fftfreq(n=params.ntall, d=params.dt*params.tscale)
        
        # tapered_spectrum = np.fft.rfft(ts_tapered_fun)
        # frequency_array = np.fft.rfftfreq(n=params.ntall, d=params.dt)
        
        processed_spectrum = Spectrum(frequency=frequency_array[0:params.nt],
                                      real_part=tapered_spectrum.real[0:params.nt],
                                      imaginary_part=tapered_spectrum.imag[0:params.nt],
                                      absolute_value=np.abs(tapered_spectrum)[0:params.nt])
        
        return processed_spectrum

    
    def detide(self, order, detiding_ds, plot=False, plot_obspy=False):
    
        detided_trace = TimeSeries()
        detided_trace.times = self.times
        detided_trace.data  = scipy.signal.detrend(self.data, type='spline', order=order, dspline=detiding_ds, plot=plot_obspy)
        
        # if plot == True:
        #     fig, ax = plt.subplots(figsize=(12,4))
        #     ax = before_after_figure(tided_trace, 'Original', tr, 'Detided', 'Detiding Visualization', ax, label_axes=True)
        #     plt.savefig('detiding.svg')
        #     plt.close()
        
        return detided_trace

    def trim(self, lag, duration, plot=False):

        times = self.times
        data = self.data
        start_trim = (60*60*lag) # seconds
        end_trim = start_trim + (60*60*duration)
        if end_trim > tr.stats.endtime:
            print( 'WARNING: TRIMMING INCOMPATIBLE WITH TRACE LENGTH' )
            
        trimmed_trace = TimeSeries()
        trimmed_trace.times = times

        bool_trimming = (times < start_trim) & (times > end_trim)
        data[bool_trimming] = 0.0
        trimmed_trace.data = data

        # if plot == True:
        #     fig, ax = plt.subplots(figsize=(12,4))
        #     ax = before_after_figure(untrimmed_trace, 'Original', tr, 'Trimmed', 'Trimming Visualization', ax, label_axes=True)
        #     # plt.savefig('trimming.svg')
        #     plt.show()
        #     plt.close()
        
        return trimmed_trace





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


def create_station_id_string(trace, use_periods=1):
    
    sts = trace.stats
    if use_periods == 1:
        id_string = sts.network +'.'+ sts.station +'.'+ sts.location +'.'+ sts.channel
    else:
        id_string = sts.network +'_'+ sts.station +'_'+ sts.location +'_'+ sts.channel
    return id_string


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


    


