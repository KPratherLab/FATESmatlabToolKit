function PeakDataTMP = read_peak(PeakFile)
% READ_PEAK reads in the files containing the peak/spectra (currently pkl)
% data into matrices. These matrices are used
% by parse_part to populate the PEAKMat data tables.  This file may be altered to  
% accomodate other file formats (other than .pkl).  Efforts should
% be made when altering the code to use the fastest data import MATLAB
% function appropiate, as this step of study creation can often be a
% bottleneck. 

fid = fopen(PeakFile);
PeakDataTMP = textscan(fid,'%f %f %f %f %f %f %f %f\n','Delimiter',',');
fclose(fid);