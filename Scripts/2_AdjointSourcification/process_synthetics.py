import numpy as np
import matplotlib.pyplot as plt
import sys
import yaml

import ProcessSeismicData as PSD

# That's temporary, there will be a proper tree later on
file_path1     = '/data/SA/IDSM_forward/build/ANMO.Z.dat'
filepath_yaml = './process_spectra.yaml'

def process_synthetic(conf, spectrum):
    # time_series = taper_and_ifft(conf, spectrum, plotting=True)
    time_series = spectrum.taper_and_ifft(conf, plotting=False)
    processed_spectrum = time_series.taper_and_fft(conf, plotting=False)

    return processed_spectrum

############################################################


# Main
if __name__ == "__main__":

    conf_py = PSD.load_yaml_file(filepath_yaml)
    config = PSD.Config("../1_SyntheticComputation/PARAMETERS")
        
    params = config.setup_run(conf_py)
    
    # Read in event and stations to process -- to do later

    # Load synthetics
    synthetic_spectrum = PSD.Spectrum()
    synthetic_spectrum.load_synthetic(file_path1, params) # In the future will take event and station as argument

    # Process spectrum
    pro_spec = process_synthetic(params, synthetic_spectrum)

    # Plot processed spectrum
    plt.figure(figsize=(12,6))
    plt.plot(pro_spec.frequency*1000.0, pro_spec.abs_val, label='IDSM')
    plt.xlabel(r'$f$ (mHz)')
    plt.ylabel('amplitude (m/s/s)')
    plt.legend()
    plt.xlim(0.3,2.7)
    plt.show()
