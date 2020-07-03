
# PMT Analysis repository #
pmt_analysis provides the basic tools used for the study of the response of the PMTs of the ICARUS experiment and their calibration. This project is currently under develelopment.

## Getting Started ##
### Get the source ###
Clone the reposotory on your personal space using `git` or fork the project to your personal GitHub page.
### External software dependencies ###
The core of the projects depends on:
* ROOT 6.15 or higher: https://root.cern.ch/downloading-root  
* Eigen: http://eigen.tuxfamily.org/index.php?title=Main_Page
* CMake: https://cmake.org/
Some task-specific scripts may requires some common python packages such as numpy, pandas, and uproot.

On an experiment GPVM machine at FNAL, all required packages are found on CVMFS and configurable as ups products thanks to  `./configure.sh`. Alternatively, the user should edit the `configure.sh` scripts pointing to their local installation of the required packages.

### Set up your personal space ###
On a GPVM machine at FNAL, the environment is configured by `configure.sh`. On a personal machine, the user has to edit `configure.sh` in the most appropriate fashion. To prepare the environment:
``` bash
chmod u+x configure.sh
./configure.sh
```
The second command has to be repeated at any new login.

## Use macros ##
During the data exploratory phase or for quick analysis tasks which doesn't involve using many files at each times, it is possible to use ROOT macros linked to the main project libraries. To do so:
``` bash
cd macros/
root -l loadLib.cc my_macro.cc
```
The script `loadLib.cc` must always be called before a root macro including any of the project libraries inside `inc/` and `src/`. Its job is to create a ROOT dictionary for the project custom libraries.

Examples macros are provided to show to the user how to structure a simple analysis project or quick data visualization.

## Compile pmt_analysis ##
For the processing of large pmt datasets  it is necessary to compile pmt_analysis. This part of the project is not yet completed, however a small example is provided in `waveform/` . More detailed instruction will come in a near future.

Compiling requires a local installation of CMake or the correct UPS product sourced by `configure.sh`. Then in the main project directory do:
``` bash
mkdir builddir/
cd builddir
cmake ../
cmake --build .
cd ../
```
If every step is successful a command called like the compiled script ( expect for the .cc part ) is created. In case the CMakeList.txt file is not changed, it is possible to run only the command `cmake --build /path/to/bulddir` targeting the build directory. 
