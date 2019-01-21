# pmu_plots



The code in this package contains a set of functions that plot statistical graphics acording to 15 minutes of PMU measurements (frequency, voltange phase angle, voltage magnitude).

Input:
 - PMU measurements, one csv file per each PMU sensor, as the examples in the folder 'data_example'
 - 2-d locations of sensors, as in 'coords.csv'
 - Nominal voltage of each sensor's bus, as in 'nominalvolts.csv'
 - Modification of 'conf' file with the correct inputs
 
 Output:
 - Folder with set of plots created by the functions
 
 To Run:
 > python plot_sigmaPCA.py conf
 or
 > python plot_autocorr.py conf



%   Copyright (c) 2019 by Mauro Escobar, IEOR Columbia University   www.columbia.edu/~me2533
