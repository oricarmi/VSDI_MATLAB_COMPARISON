# VSDI_MATLAB_COMPARISON

To run all of the methods on a specific experiment, please contact us and we will provide a link to download the raw data requested.
Once the raw data is available, run the first section of "Real_data_comparison.m" file.
The path to the raw .mat file needs to be provided in fname. 
** To get the same results as in the article, make sure the "paramsori.csv" file, which is a csv with all the user defined parameters for each algorithm, is the same as Table2 in the supp. material.

To plot the figures from the provided results (outputs) after running all experiments: 
load an output .mat file and then run the file plot_figures (need to load each output .mat file separately to matlab and then run).
Take note that there are results for the simulation and experimental data, and each figure is either for simulation or experimental so the appropriate results should be selected.
All of the output files are in the DATA folder of this repo (uploaded via git lfs). 

For any questions, feel free to contact oricarmi92@gmail.com
