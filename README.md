# PMCtrack
Polar Mesoscale Cyclones (PMC) tracking algorithm

## Dependencies
* cmake
* netCDF
* Fortran-90 compiler

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
Coming soon

### Based on the original version (prior to changes shown here)
* Watanabe, S. I., H. Niino, and W. Yanase, 2018: Composite Analysis of Polar Mesocyclones over the Western Part of the Sea of Japan. Monthly Weather Review, 146, 985–1004, https://doi.org/10.1175/MWR-D-17-0107.1.

* Watanabe, S. I., H. Niino, and W. Yanase, 2017: Structure and Environment of Polar Mesocyclones over the Northeastern Part of the Sea of Japan. Monthly Weather Review, 145 (6), 2217–2233, DOI:10.1175/MWR-D-16-0342.1.

* Watanabe, S. I., H. Niino, and W. Yanase, 2016: Climatology of Polar Mesocyclones over the Sea of Japan Using a New Objective Tracking Method. Monthly Weather Review, 144 (7), 2503-2515, DOI:10.1175/MWR-D-15-0349.1.

## Analysing output
To accompany the tracking algorithm itself, a software package for postprocessing and analysis has been written in Python: [octant](https://github.com/dennissergeev/octant).
