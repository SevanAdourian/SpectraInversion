#!/usr/bin/env python
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client 
from scipy.signal import welch
from scipy import signal
import numpy as np
import math
from scipy.signal import periodogram, coherence
import matplotlib.pyplot as plt
from scipy.optimize import fmin

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times')
mpl.rc('font',size=22)

def trim_tr(tr, lag, duration, plot=False):
    print(lag, duration)
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

def detide(tr, order, detiding_ds, plot=False, plot_obspy=False):
        
    tided_trace = tr.copy()
    tr.detrend(type='spline', order=order, dspline=detiding_ds, plot=plot_obspy)

    if plot == True:
        fig, ax = plt.subplots(figsize=(12,4))
        ax = before_after_figure(tided_trace, 'Original', tr, 'Detided', 'Detiding Visualization', ax, label_axes=True)
        plt.savefig('detiding.svg')
        plt.close()

    return tr   

def detrend(tr, detrend_type, plot=False):
    
    trend_trace = tr.copy()
    tr.detrend(type=detrend_type)
    if plot == True:
        fig, ax = plt.subplots(figsize=(12,4))
        ax = before_after_figure(trend_trace, 'Original', tr, 'Detrended', 'Detrending Visualization', ax, label_axes=True)
        plt.savefig('detrending.svg')
        plt.close()
        
    return tr

def taper(tr, max_perc, ttype, max_length, side, plot=False, plot_envelope=False):
    
    untapered_trace = tr.copy()
    tr.taper(max_perc, type=ttype, max_length=max_length, side=side)
    if plot == True:
        fig, ax = plt.subplots(figsize=(12,4))
        ax = before_after_figure(untapered_trace, 'Original', self.trace, 'Tapered', 'Tapering Visualization', ax, label_axes=True)
        plt.savefig('tapering.svg')
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

component = 'Z'
order = 2 # 
max_perc = 0.01
ttype = 'hamming'
max_length = None
side = 'both'
corners = 4
zerophase = 'False'
low_corner = 0.2e-3
high_corner = 1.0e-3
lag = 5
duration = 96

stime = UTCDateTime('2011-070T05:46:00')
etime=stime + (lag+duration+1)*3600

mf, Mf=0.15, 1.0
sta = 'TRIS'
length = 400
client = Client('IRIS')

inv = client.get_stations(network='IU', sta=sta, 
   channel='LHZ,LDO', location='00', starttime=stime, 
               endtime=etime, level='response')


st = client.get_waveforms(network='IU', station=sta, location='00',
                          channel='LHZ,LDO', starttime=stime,
                          endtime=etime)


st.merge(fill_value=0)

# st = st.detrend('demean')
# st = st.detrend('linear')
# st = st.taper(type=ttype,max_percentage=max_perc)

# st = st.select(component=component)
detiding_ds = int(1/(0.0001*2*st[0].stats.delta))

# Detide
st = detide(st, order, detiding_ds, plot=False, plot_obspy=False)
# st = trim_tr(st, lag, duration, plot=False) # Trim
st = detrend(st, 'linear') # Detrend
st = taper(st, max_perc, ttype, max_length, side) # Taper
st = bandpass(st, low_corner, high_corner, corners, zerophase, plot=False) # Filter

# st.detrend('linear')
# st.detrend('constant')
st.sort()
print(st)

# for tr in st:
#     tr.data *= signal.get_window(('kaiser', 2. * np.pi), tr.stats.npts)

# Zero pad
# for tr in st:
#     tr.data = np.pad(tr.data, (
#         int((length * 400 * 24 / 10. - tr.stats.npts) / 2.),
#         int((length * 400 * 24 / 10. - tr.stats.npts) / 2.)), 'edge')
#     print(tr)

NFFT = 2 ** (math.ceil(math.log(st[0].stats.npts, 2)))
for idx2, tr in enumerate(st):
    if tr.stats.channel == 'LHZ':
        tr = tr.remove_response(inventory=inv, output='ACC', water_level=60, plot=False,
                                pre_filt=[0.0001, 0.00011, 0.0051, 0.0052])
        
    # fper, pper = periodogram(tr.data, fs=tr.stats.sampling_rate, nfft=NFFT,
    #                          scaling='spectrum')
    # pper, fper = pper[1:], fper[1:]

    # print(len(pper), len(fper))
    f, p = periodogram(tr.data, fs=tr.stats.sampling_rate, nfft=NFFT,
                             scaling='spectrum')
    p, f = p[1:], f[1:]

    # p = np.fft.fft(tr.data, n=NFFT)
    # f = np.fft.fftfreq(NFFT, d=1.0/tr.stats.sampling_rate)
    # p *= 1.e9
    f *= 1.e3
    # pper = np.sqrt(np.abs(pper))
    # pper *= 1.e9
    # fper *= 1.e3
    
    if tr.stats.channel == 'LHZ':
        # inv_resp = inv.get_response(tr.id, tr.stats.starttime)
        # resp, _ = inv_resp.get_evalresp_response(tr.stats.delta, NFFT, 'ACC')
        # resp = resp[1:]
        # Convert units to nm/s/s
        # p = np.sqrt(p / (np.abs(resp) ** 2)) * 10 ** 9
        p = np.sqrt(p) * 1.e9
        pLHZ = p[(f >= mf) & (f <= Mf)]
        f = f[(f >= mf) & (f <= Mf)]
        # pperLHZ = pper[(fper >= mf) & (fper <= Mf)]
        # fper = fper[(fper >= mf) & (fper <= Mf)]
        
    else:
        p = np.sqrt(p) 
        # f *= 1000.
        pLDO = p[(f >= mf) & (f <= Mf)]
        f = f[(f >= mf) & (f <= Mf)]

def presscorrt(x):
    return pLHZ - x*pLDO

def resi(x):
    val = np.sum(presscorrt(x)**2)/len(presscorrt(x))
    return val

bf = fmin(resi, [0.])

f2, Cxy  = coherence(st[0].data, st[1].data, fs=1., nperseg=NFFT/4)

f2 *= 1000.
Cxy = Cxy[(f2 >= mf) & (f2 <= Mf)]
f2 = f2[(f2 >= mf) & (f2 <= Mf)]

corr = presscorrt(bf[0])

fig = plt.figure(1,figsize=(12,12))
plt.subplot(3,1,1)
# plt.plot(f,np.sqrt(np.abs(pLHZ))-np.sqrt(np.abs(corr)), label=(st[1].id).replace('.',' '))
plt.plot(f,np.abs(pLHZ), label=(st[1].id).replace('.',' '))
plt.plot(f,np.abs(presscorrt(bf[0])), label='Corrected')
plt.ylabel('Acceleration (nm/s/s)')
plt.legend()
plt.subplot(3,1,2)
plt.plot(f,pLDO, label=(st[0].id).replace('.',' '))
plt.ylabel('Pressure')
plt.subplot(3,1,3)
plt.plot(f2,Cxy, label='Coherence')
plt.ylabel('Coherence')
plt.xlabel('Frequency (mHz)')
plt.legend()
# plt.show()
# plt.savefig('Example_' + sta + '.PNG',format='PNG',dpi=400)
# plt.close()

