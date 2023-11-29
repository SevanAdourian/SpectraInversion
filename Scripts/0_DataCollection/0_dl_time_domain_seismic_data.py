import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import glob
import yaml
import pyasdf

from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.mass_downloader import GlobalDomain, Restrictions, MassDownloader

def load_yaml_file(file_path):
    with open(file_path, 'r') as file:
        return yaml.safe_load(file)

def ensure_directory_exists(path):
    # Check if the directory exists
    if not os.path.exists(path):
        # If not, create the directory
        os.makedirs(path)
        print(f"Directory '{path}' created.")
    else:
        print(f"Directory '{path}' already exists.")
        

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <yaml_file>")
        sys.exit(1)

    # Get the YAML file name from the command-line argument
    yaml_file_name = sys.argv[1]
    
    # Load YAML content
    param = load_yaml_file(yaml_file_name)    
    
    waveform_dir = f"{param['basedir']}/Data/{param['event_ID']}/mseeds/"
    stations_dir = f"{param['basedir']}/Data/{param['event_ID']}/xmls/"
    asdf_dir     = f"{param['basedir']}/Data/{param['event_ID']}/asdf/"
    
    ensure_directory_exists(waveform_dir)
    ensure_directory_exists(stations_dir)
    ensure_directory_exists(asdf_dir)
    
    # Moved to yaml file
    origin_time = UTCDateTime(param['time_onset'])
    
    restrictions = Restrictions(starttime=origin_time,
                                endtime=origin_time + (param['lag']+param['duration']+1)*3600,
                                reject_channels_with_gaps=True,
                                network="G,IU,II,",
                                channel="LHZ,LDO,LDI",
                                sanitize=True,
                                minimum_length=0.95,
                                channel_priorities=["LDO","LDI"],
                                location_priorities=["00"]
                                )

    domain = GlobalDomain() # CircularDomain, GlobalDomain, RectangularDomain
    
    c = Client("IRIS")
    inv = c.get_stations(network="G,IU,II", channel="LHZ,LDO,LDI", 
                         level="response", 
                         starttime=origin_time, endtime=origin_time + \
                         (param['lag']+param['duration']+1)*3600)
    
    mdl = MassDownloader(providers=["IRIS"])
    mdl.download(domain, restrictions, mseed_storage=waveform_dir, stationxml_storage=stations_dir)

    # USE INV TO PRINT RELEVANT INFO FOR DOWNLOADED DATA
    
    # Writing to an asdf file
    filename = f"{asdf_dir}/{param['event_ID']}.h5"
    with pyasdf.ASDFDataSet(filename, compression="gzip-3") as ds:
        waveform_files = glob.glob(f"{waveform_dir}*LHZ*.mseed")
        for _i, filename in enumerate(waveform_files):
            print("Adding file %i of %i ..." % (_i + 1, len(waveform_files)))
            ds.add_waveforms(filename, tag="raw_waveform", event_id=param['event_ID'])

        pressure_files = glob.glob(f"{waveform_dir}*LD*.mseed")
        for _i, filename in enumerate(pressure_files):
            print("Adding file %i of %i ..." % (_i + 1, len(pressure_files)))
            ds.add_waveforms(filename, tag="pressure", event_id=param['event_ID'])

        station_files = glob.glob(f"{stations_dir}/*.xml")
        for _i, filename in enumerate(station_files):
            print("Adding station %i of %i ..." % (_i + 1, len(station_files)))
            ds.add_stationxml(filename)

