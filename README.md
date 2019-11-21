# PMCtrack
Polar Mesoscale Cyclones (PMC) tracking algorithm

## Dependencies
* cmake
* netCDF with Fortran enabled
* Fortran-90 compiler

### Using conda environment
The easiest way to install all the dependencies is to use [conda](https://docs.conda.io/en/latest/).
Once you have downloaded conda ([Miniconda](https://docs.conda.io/en/latest/miniconda.html) is recommended), create an isolated environment by running
```bash
conda create -n pmctrack -c conda-forge cmake libgfortran fortran-compiler netcdf-fortran
```
This should fetch all the necessary libraries which then can be used for compilation.

## Installation
### Unix
Run the script
```bash
./INSTALL.sh
```
All options are stored in `CMakeLists.txt` file. You can add shortcut for your platform.

## Running the program
```bash
./track.x [settings.conf]
```

## Papers
### New version
* Sergeev, D., I. A. Renfrew, T. Spengler, A. Terpstra, and S. I. Watanabe, 2019: North Atlantic polar mesoscale cyclones in ERA5 and ERA-Interim reanalyses. _in prep._.

### Based on the original version (prior to changes committed here)
* Watanabe, S. I., H. Niino, and W. Yanase, 2018: Composite Analysis of Polar Mesocyclones over the Western Part of the Sea of Japan. _Monthly Weather Review_, 146, 985–1004, https://doi.org/10.1175/MWR-D-17-0107.1.

* Watanabe, S. I., H. Niino, and W. Yanase, 2017: Structure and Environment of Polar Mesocyclones over the Northeastern Part of the Sea of Japan. _Monthly Weather Review_, 145 (6), 2217–2233, DOI:10.1175/MWR-D-16-0342.1.

* Watanabe, S. I., H. Niino, and W. Yanase, 2016: Climatology of Polar Mesocyclones over the Sea of Japan Using a New Objective Tracking Method. _Monthly Weather Review_, 144 (7), 2503-2515, DOI:10.1175/MWR-D-15-0349.1.

## Analysing output
To accompany the tracking algorithm itself, a software package for postprocessing and analysis has been written in Python: [octant](https://github.com/dennissergeev/octant).
