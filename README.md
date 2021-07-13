# Example codes
The repository contains a few examples of the codes written by me.  Some of the codes require calibration files, some packages such as CFTITSIO, PIL, so you may not be able to run/compile those codes. Below is the list of available codes,

1. SSMSI:- The tool written for AstroSat Scanning Sky Monitor image reconstruction. The code performs forward fitting using Singular Value Decomposition and gives the available sources in the FOV. This code is implemented using MPI, CUDA and pthreads to run the reconstruction in the parallel environment. The code has the functionality to use multiple GPU cards. It is extensively used on a cluster having 4 nodes and 12 GPU cards (3 cards/node).

2. SimEvents: This code simulates a list of events for AstroSat CZTI using Ray Tracing and writes Detector Plane Histogram and event file. Later simulated files are used for further analysis. This code is implemented using MPI and CUDA (4 GPU cards). 

3. GetCorrelationGPU:- This code generates a cross-correlation image given an event file, finds available sources and compute their intensities. The code is written in CUDA

4. CUDA_DBNN:- The code is a GPU version of the code written for galaxies, star classification written by Prof. Sajeeth. This code is developed in CUDA.

5. mpi_autodbnn:- The code is MPI version of the code written for galaxies, star classification written by Prof. Sajeeth. This code is developed in MPI.  

6. CRTS:- This tool is used to compute different statistical parameters of 500 million CRTS light curves. This code is developed in MPI.

7. czti_pipeline:- CZTI level-1 to level-2 data reduction pipeline. The pipeline is written in C/C++. The pipeline is publically available at http://astrosat-ssc.iucaa.in/uploads/czti/czti_pipeline_20180308_V2.1.tar 

8. Bayesian_spectrum:- Bayesian-based spectrum reconstruction tool used for AstroSat CZTI. The code is written in C.


Below are the links for a few websites, web-based tools developed/managed by me.

1. A tool to convert time in various format to AstroSat time and vice versa: http://astrosat-ssc.iucaa.in:8080/astrosattime/

2. A tool to generate orbit files for AstroSat CZTI: http://astrosat-ssc.iucaa.in:8080/orbitgen/

3. A tool to check visibility of a source in given time period: http://astrosat-ssc.iucaa.in:8080/AstroVisCal/

4. AstroSat Science Support Cell Website: http://astrosat-ssc.iucaa.in

5. AstroSat Website: http://astrosat.iucaa.in

6. IUCAA HPC website: http://hpc.iucaa.in



