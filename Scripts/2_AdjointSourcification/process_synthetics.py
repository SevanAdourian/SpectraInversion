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
import pickle

from scipy.signal import periodogram, coherence, get_window, detrend
from scipy.optimize import fmin
from scipy.interpolate import interp1d

import ProcessSeismicData as PSD

import pdb


# That's temporary, there will be a proper tree later on
file_path1 = '/data/SA/IDSM_forward/build/AFI.Z.dat'


def process_synthetic(conf, spectrum):
    # time_series = taper_and_ifft(conf, spectrum, plotting=True)
    time_series = spectrum.taper_and_ifft(conf, plotting=True)
    processed_spectrum = time_series.taper_and_fft(conf, plotting=True)

    return processed_spectrum

############################################################


# Main
if __name__ == "__main__":

    config = PSD.Config("/data/SA/IDSM_forward/setup/PARAMETERS")

    params = config.setup_run("./conf.yaml")

    # Read in event and stations to process -- to do later

    # Load synthetics
    synthetic_spectrum = PSD.Spectrum()
    synthetic_spectrum.load_synthetic(file_path1, params) # In the future will take event and station as argument

    # time_series = taper_and_ifft(params, synthetic_spectrum, plotting=True)
    pro_spec = process_synthetic(params, synthetic_spectrum)
    
    sys.exit()
    # Read data -- Also will be done later
    
    # Read in processing parameters for each event/station pair -- for now just doing one at a time, hardcoded
    # Later will be read from the metadata of the file.
    
    
    # Probably here do equivalent of picking to accept/reject spectrum

    # Compute misfit

    # Backprocess adjoint source

    # Write in format that is accepted by IDSM-adjoint
    


# ##############################################3


# Load associate data
# file = open("C060994A_IU_AFI_60.pkl", 'rb')
# obj_file = pickle.load(file)

