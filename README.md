# SpectraInversion
workflow tentative for spectra inversion using IDSM and adjoint method - work in progress

# Install required packages
conda env create -f environment.yml

conda activate specinv

cd ./Scripts/pythonlibs

conda develop . # Adds the local libraries to the PYTHONPATH

cd ../../

# Run data download
cd Scripts/0_DataCollection/

python3 0_dl_time_domain_seismic_data.py download_events.yaml

It should download all mseeds and xmls for accelerations (currently for LHZ) and barometer data, if available (channels LDO and LDI). Data are stored in <basedir>/Data/<event_ID>/

1 asdf file will be created per event containing the traces and metadata <event_ID.h5>.

# Run data processing and spectra computations
python3 1_process_spectra.py process_spectra.yaml

This will process the data and compute amplitude spectra. These spectra will be appended to another asdf file <event_ID>_pro.h5 in an auxiliary_data field named ProcessedSpectra organized in trace.id

Example to plot a spectra from the asdf file once it was created:

```
with pyasdf.ASDFDataSet(<event_ID>_pro.h5, compression="gzip-3", mode='r+') as ds:
    spec_ds = ds.auxiliary_data.ProcessedSpectra[tr_id]
    sfreq = spec_ds.parameters["start_freq"]
    nfreq = spec_ds.parameters["nfreq"]
    dfreq = spec_ds.parameters["dfreq"]
    efreq = sfreq+(nfreq-1)*dfreq
    freq_array = np.arange(sfreq, efreq+dfreq, dfreq)

    plt.figure()
    plt.plot(freq_array, np.array(spec_ds.data))
```