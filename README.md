# SpectraInversion
workflow tentative for spectra inversion using IDSM and adjoint method - work in progress

# Install required packages
conda env create -f environment.yml
conda activate specinv
cd ./Scripts/pythonlibs
conda develop . # Adds the local libraries to the PYTHONPATH
cd ../../

# Run data selection
cd Scripts/0_DataCollection/