ARSER is a Python package for identifying
periodic expression profiles in analyzing circadian microarray data
and has been released under the GPL

Table of Contents
=================
- Pre-installation
- Usage
- Examples
- Input/Output Files
- Additional Tools
- Citations
- Additional Information


Pre-installation
================
ARSER package is implemented by Python calling R program. Before using the package, 
please install the following software and packages first:

- Programing environments:    
    1. Python version 2.5 or later
    2. R version 2.7 or later
- Packages:
    3. scipy version 0.7 or later
    4. numpy version 1.1 or later
    5. matplotlib version 0.99 or later
    6. Rpy
- Tips:
    To avoid wading through all the details (and potential complications) on Installation, 
    the easiest thing for you to do is use one of the pre-packaged python distributions 
    that already provide scipy/numpy/matplotlib built in. The Enthought Python Distribution (EPD) 
    for Windows, OS X or Redhat is an excellent choice. Another alternative for Windows users is Python (x, y).
    
Usage
=====
Command-line running:  

    python arser.py input_file_name output_file_name start(optional) end(optional) default_period(optional)
Options:

    ARSER searches period in the range [start, end]
    start: period searching range start, default 20h
    end: period searching range end, default 28h
    default_period: default period used by ARSER for searching, default 24h
example: 

    $ python arser.py data.txt output.txt >& log.txt

Input/Output Files
==================
Note: Sample input and output files can be found in the *examples/* subdirectory.

* Input: 
	* Microarray data file with a header line which records the time-points. The 1st column is probesets, other columns are expression values over time. It is assumed that the samples are linearly spaced (e.g., one point every 4 hrs, etc). *The current version of ARSER does NOT allow for non-linear sampling.*
         
* Output:
	* The 1st column is probesets, other columns are values of parameters as followed:

            mean        -> mean value for raw y values
            period      -> period identified by ARSER
            amplitude   -> amplitude for single cosine model
            phase       -> phase for single cosine model
            R2          -> R square of regression curve
            R2adj       -> adjusted R square of regression curve
            coef_var    -> (standard deviation) / mean
            pvalue 			-> F test for testing significant regression model
            qvalue      -> FDR by stroy method
            FDR_BH      -> FDR by BH method
            filter_type -> filtering for noise by ARSER
                            0 -- no filtering
                            1 --  filtering
            ar_method   -> methods for autogressive model fitting
                            'mle' -- maximum likelihood estimate
                            'burg' -- burg algorithm
                            'yule-walker' -- yule-walker equations
                            'default' -- harmonic regression with 24h
            period_number -> number of cycles in time series
                 
Additional Tools
================
See the README file in the tools subdirectory.

Citations
=======================
Rendong Yang and Zhen Su, Analyzing circadian expression data by harmonic regression based on autoregressive spectral estimation *Bioinformatics*. 2012 Jun 15;26(12):i168-74.

Additional Information
=======================
ARSER website at 
http://bioinformatics.cau.edu.cn/ARSER

A webserver of ARSER algorithm can be found at
http://bioinfo.cau.edu.cn/BioClock

Please contact me if you have suggestions for improvement or 
if any problem arises in the use of the program, 
or the interpretation of the results.

Rendong Yang

cauyrd@gmail.com
