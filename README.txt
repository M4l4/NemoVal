This package provides quick and hopefully painless plotting of data from the NAA models.

There are three scripts that can be run. The first is single_plots.py, which is used to make plots of single variables
and anomalies. The second is three_model.py, which is used to make a three graph comparison between models with
optional anomalies. The final one is transect.py, which is used to make a transect.

All of three files are configured near the top of the files, and have instructions for use there.

Each file outputs in to the ./plots directory, and destroys it every run. If you want to keep the files, make a copy.

Requirements can be found in environment.yml, which can be easily imported in to a conda environment by running the
command 'conda env create -f environment.yml'. CDO and NetCDF must be installed separately from the python bindings.

Known issues:
    This file could use some work.
    The transect line is not ideal, this will be fixed soon.
    Scale factors need to be updated in both transect.py and three_model.py to differentiate between CMOC and PISCES,
this will be fixed soon.
    single_plots.py has no scaling applied to any data.
    Taylor plots in single_plots.py do not work.