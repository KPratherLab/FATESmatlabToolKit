function PeakDataTMP = read_peak(PeakFile)
% READ_PEAK reads in the files containing the peak/spectra (currently pkl)
% data into matrices. These matrices are used
% by parse_part to populate the PEAKMat data tables.  This file may be altered to  
% accomodate other file formats (other than .pkl).  Efforts should
% be made when altering the code to use the fastest data import MATLAB
% function appropiate, as this step of study creation can often be a
% bottleneck.  PEAKFile is the full path to the spectra file.  PEAKDataTMP is a
% matrix where each row contains data for a single peak.  Each column
% contains a single data type (m/z, area, etc).  PEAKDataTMP must be output
% in this format to be compatible with parse_part.

fid = fopen(PeakFile); %open file
PeakDataTMP = textscan(fid,'%f %f %f %f %f %f %f %f\n','Delimiter',','); %read file
PeakDataTMP = cell2mat(PeakDataTMP); %convert from cell array to matrix, each row contains data for one peak
fclose(fid); %close file