import numpy as np
import matplotlib.pyplot as plt
import scipy
import sys
import re
import os
import glob
import math

# Utils
def read_complex_fortran(filename, regex):
    coord = []
    complex_array = []
    #
    with open(filename, 'r') as file:
        data = file.read()
    # 
    # Find all matches of the complex number pattern
    matches = re.findall(regex, data)
    #
    # Convert the matches into a list of complex numbers
    coord = [float(match[0]) for match in matches]
    complex_array = [complex(float(match[1]), float(match[2])) for match in matches]
    #
    return coord, complex_array

# Define the function
def fcal(f1, f2, dt, tout, mex, qex):
    fn = 0.5 / dt
    
    if fn < f2:
        raise ValueError("f2 is greater than the Nyquist frequency for the time step")
    
    ep = mex / tout
    print(f"ep = {ep}")
    df = ep / (2.0 * np.pi * qex)
    print(f"df = {df}")
    # df = ep / qex
    nt = int(np.ceil(1.0 / (df * dt)))
    print(f"nt = {nt}")
    # ne = np.floor(math.log(nt) / math.log(2)) + 1
    # print(f"ne = {ne}")
    # nt = int(2 ** ne)
    # print(f"nt2 = {nt}")
    
    # df = 1.0 / (nt * dt)
    # print(f"df2 = {df}")
    i1 = np.floor(f1 / df)
    print(f"i1 = {i1}")
    f1 = (i1) * df

    i2 = np.ceil(f2 / df)
    print(f"i2 = {i2}")
    f2 = (i2) * df
    
    return f1, f2, df, ep, nt, i1, i2


def dirac_delta(index, length, dx):
    if index < 0 or index >= length:
        raise ValueError("Index must be within the range [0, length - 1]")
    
    delta = np.zeros(length,dtype=complex)
    delta[index] = complex(0.5/dx,-0.5/dx)
    
    return delta

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


    # Generate sample data points
# num_samples = 2048
facf = 0.1
fact = 0.5
f1_inp = 0.00/1e3
f2_inp = 1.00/1e3
t1_inp = 0.
t2_inp = 128.*3600
dt = 1. / (2.*f2_inp)
index_check = 1515 - 1

# FCAL usage
f_norm = np.sqrt(6.6723e-11*np.pi*5515.)
t_norm = 1./f_norm
mex = 5
qex = 4
f1, f2, df, ep, nt, i1, i2 = fcal(f1_inp/f_norm, f2_inp/f_norm, dt/t_norm, t2_inp/t_norm, mex, qex)

print(".... Frequency parameters:")
print(".... ... df (mHz) =", df*f_norm*1000.)
print(".... ... nt       =", nt)
print(".... ... f1 (mHz) =", f1*f_norm*1000.)
print(".... ... f2 (mHz) =", f2*f_norm*1000.)
f12 = f1+facf*(f2-f1)
f21 = f2-facf*(f2-f1)
t1 = t1_inp/t_norm
t2 = t2_inp/t_norm
t12 = t1+fact*(t2-t1)
t21 = t2-fact*(t2-t1)
dt = dt/t_norm

##################
# read 3d spectra
directory_path = "/home/kalamata/sevan/DMCC/bin"
file_pattern = "spectra_3d_*.dat"
file_list = glob.glob(os.path.join(directory_path, file_pattern))
file_list.sort()
data_arrays = []
info_list = []
processed_arrays = []
time_series = []

# Dealing with 1d
spectrum_1d = np.loadtxt(f"{directory_path}/spectra_1d.dat")
data_arrays.append(spectrum_1d)
info_list.append('00')

# Now reading with 3d
for file_path in file_list:
    file_name = os.path.basename(file_path)
    info = file_name.split("_")[-1].split(".")[0]
    data = np.loadtxt(file_path)  

    data_arrays.append(data)
    info_list.append(info)

x = data_arrays[0][:,0]
x = x / f_norm / 1000.
x = np.concatenate((x,-np.flipud(x)[:-1]))
nt = len(data_arrays[0])
ntall = 2*nt-1
time_support = np.arange(0, (ntall*dt)-1e-6, dt)

for ids in range(0,len(data_arrays)):
    print(f"Processing spectrum of perturbation {ids}")
    fun = np.zeros(nt, dtype=complex)
    fun.real = data_arrays[ids][:,1]
    fun.imag = data_arrays[ids][:,2]

    # Taper in frequency
    ftapered_fun = np.zeros(nt, dtype=complex)
    
    for i in range(0,len(ftapered_fun)):
        hann_coeff = hann(x[i], f1, f12, f21, f2)
        ftapered_fun[i] = hann_coeff*fun[i]

    # Constructing the negative frequencies
    ftapered_fun_cplx = ftapered_fun.astype(complex)
    ftapered_fun_with_neg = np.zeros(ntall, dtype=complex)
    # ftapered_fun_with_neg = np.zeros(nt, dtype=complex)
    ftapered_fun_with_neg[:nt] = ftapered_fun_cplx
    ftapered_fun_with_neg[nt:] = np.flipud(np.conjugate(ftapered_fun_cplx))[:-1]

    # Calculate the inverse Fourier transform
    fun_ifft = np.fft.ifft(ftapered_fun_with_neg)

    ttapered_fun = np.zeros(ntall,dtype=complex)
    # ttapered_fun = np.zeros(nt,dtype=complex)

    for i,t_value in np.ndenumerate(time_support):
        hann_coeff = hann(t_value, t1, t12, t21, t2)
        ttapered_fun[i] = hann_coeff*fun_ifft[i]*np.exp(ep*t_value)

    time_series.append(ttapered_fun)
    fun_fft = np.fft.fft(ttapered_fun)
    processed_arrays.append(fun_fft)
    
# Plot the spectra and recover information at given index
delta_u = np.array([float(x) for x in output_fortran.split()])
delta_u = delta_u*2.
delta_adj = np.zeros(len(processed_arrays), dtype=complex)
print(f"--------------------- delta adjoint ---------------------")
for ids in range(1,len(processed_arrays)):
    delta_adj[ids] = processed_arrays[ids][index_check] - processed_arrays[0][index_check]
    print(f"{delta_adj[ids].real} / {delta_u[ids-1]}")

ratio = delta_adj[1::].real/delta_u

print(f"---------------------- ratio -----------------")
print(ratio)
# Re-dimension
x = x*f_norm*1000.
time_support = time_support*t_norm
plt.figure(figsize=(10, 6))
for ids in range(0,len(processed_arrays)):
    plt.plot(x, np.abs(processed_arrays[ids]), label=f"{info_list[ids]}")
    plt.axvline(x=x[index_check],color='black')

dper = np.arange(0.1,1.6,0.1)
plt.figure(figsize=(10,6))
plt.plot(dper,delta_adj.real[1::], label='full')
plt.plot(dper, np.array(delta_u),label='adjoint')
plt.legend()
plt.xlim(0,1.6)

# plt.figure(figsize=(10, 6))
# for ids in range(0,len(processed_arrays)):
#     plt.plot(time_support, time_series[ids], label=f"{info_list[ids]}")

plt.legend()
plt.tight_layout()
plt.show()
