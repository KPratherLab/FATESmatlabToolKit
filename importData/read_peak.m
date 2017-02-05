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

%% file format for Prather ATOFMS
%comment this out or change if not using Prather ATOFMS pkl file format
fid = fopen(PeakFile); %open file
PeakDataTMP = textscan(fid,'%f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32\n','Delimiter',','); %read file
PeakDataTMP = cell2mat(PeakDataTMP); %convert from cell array to matrix, each row contains data for one peak
fclose(fid); %close file

%% file format for COMMERCIAL ATOFMS
%comment this out or change if not using commerical ATOFMS pkl file format
% 
% %read in spectra data
% fid = fopen(PeakFile);
% PeakDataTMP = textscan(fid,'%*s%s%*s','Delimiter',{'{','}'}); %reads in spectra as single string
% PeakDataTMP = PeakDataTMP{1}; 
% PeakDataTMP = replace(PeakDataTMP,{'(' ')' ','},''); %take out any special characters
% fclose(fid);
% 
% %read in number of peaks
% fid = fopen(PeakFile);
% numPeaks = textscan(fid,'%*s %*s %*s %*f %*f %d  %*[^\n]','Delimiter',',');
% numPeaks = numPeaks{1}; %need numPeaks for each spectra to be able to reshape spectra matrix later
% numSpectra = length(numPeaks);
% numPeaks = num2cell(numPeaks);
% fclose(fid);
% 
% %convert spectra string to numbers and reshape
% %each row will contain info for a single peak
% PeakDataTMP = cellfun(@(x,y) (sscanf(x,'%f',[5 y]))',PeakDataTMP,numPeaks,'UniformOutput',0);
% for i = 1:numSpectra
%     PeakDataTMP{i}(:,6) = i; %add particle ID to each peak row
% end
% PeakDataTMP = cell2mat(PeakDataTMP);
% 
% %add spectra polarity (pos or neg) to each peak row
% PeakDataTMP(:,7) = 0; 
% PeakDataTMP(PeakDataTMP(:,1) > 0,7) = 1;

%% file format for ALABAMA
%comment this out or change if not using ALABAMA pkl file format
% 
% %format string to read in spectra from HitParticle file
% formatHit = '%*s %*f';
% for i = 1:500
%     formatHit = strcat(formatHit,' %f32');
% end
% 
% %read in spectra
% fid = fopen(PeakFile);
% PeakDataTMP = textscan(fid,formatHit,'Delimiter','\t');
% PeakDataTMP = cell2mat(PeakDataTMP); %make matrix of spectra data
% fclose(fid);
% sizeHit = size(PeakDataTMP,1);
% PeakDataTMP = PeakDataTMP';
% 
% %set up some variables to add additional peak info
% %mz values for a spectra
% MZValue = single([-1*(1:250) 1:250]);
% MZValue = MZValue';
% numMZ = length(MZValue);
% 
% %spectra polarity values
% Polarity = single([zeros(numMZ/2,1); ones(numMZ/2,1)]);
% 
% %particle index values
% tmp = single(1:sizeHit);
% tmp = repmat(tmp,numMZ,1);
% tmp = tmp(:);
% 
% %Combine all Data
% PeakDataTMP = PeakDataTMP(:); %area
% PeakDataTMP(:,2) = repmat(MZValue,sizeHit,1); %m/z
% PeakDataTMP(:,3) = repmat(Polarity,sizeHit,1); %polarity
% PeakDataTMP(:,4) = tmp; %particle identifier
% % 
