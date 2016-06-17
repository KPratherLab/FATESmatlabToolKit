function [Spectrum] = get_spectrum(PID,Polarity,ColList)
% GET_SPECTRUM returns spectra for particles 
% Call as [Spectrum] = get_spectrum(PID,Polarity,ColList)
% 
% Where PID is set of paritcle identifiers stored as a nx2 matrix
%   PID(:,1) = InstID
%   PID(:,2) = PartID
%
% Polarity specifies the spectrum polarity as
%   Polarity = 0 - negative spectra
%   Polarity = 1 - positive spectra
%   Polarity = 2 - negative and positive spectra
% 
% Spectrum is a cell matrix with one row for each particle
% that has these columns
%   Spectrum{i,1} = MZ 
% 
% Additional columns contain data from columns in ColList.  ColList is
% optional and defaults to
% {'Area','RelArea','Height','BlowScale'}.  So the default
% cell matrix has the form
%   Spectrum{i,1} = MZ 
%   Spectrum{i,2} = AREA 
%   Spectrum{i,3} = RELAREA 
%   Spectrum{i,4} = HEIGHT 
%   Spectrum{i,5} = BLOWSCALE 
%
% The contents of each element are vectors with one element for each 
% peak in the spectrum.  The vectors are empty for particles which do
% not have a spectrum of Polarity.  Note that the lengths of these 
% vectors differ among particles.
%
% See also GET_COLUMN, GET_INT_SPECTRUM

% YAADA - Software Toolkit to Analyze Single-Particle Mass Spectral Data
%
% Copyright (C) 1999-2000 California Institute of Technology
% Copyright (C) 2001-2008 Arizona State University
% Copyright (C) 2008 Jonathan O. Allen

% Jonathan O. Allen  03 Jul 00

% Added ColList to allow changes in data structure.  This also
% improves performance.
% JOA  06 Oct 01

% Removed Width from default columns returned
% Fixed bugs
% - id must match columns
% - data returned in correct order for unsorted IDs
% JOA  2003-01-07

% Added Verbose as parameter
% JOA  2003-08-28

% Fix for @column/subsref returns column vector
% JOA  2004-04-07

% Removed flag from intersect(A,B,'sorted')
% JOA  2005-09-15

% Commented out writing the peakid list to the output Spectrum to improve
% runtime.  Made various other coding changes to improve runtime
% Camille Sultana May 2014

global PEAKFlds YAADA PEAKMat STUDY

%verify inputs
if nargin < 2 || nargin > 4
  error('Call as [Spectrum] = get_spectrum(PID,Polarity,ColList,Verbose)');
end
if (~(size(PID,2)==2))  %PID is now Nx2 matrix of instid, partids
  error('Invalid PID');
end
if ~isnumeric(Polarity) || length(Polarity) > 1
  error('Expecting scalar for Polarity');
end
if Polarity < 0 || Polarity > 2
  error('Expecting 0 (negative), 1 (positive), or 2 (both) for Polarity');
end
if ~exist('ColList','var')
  ColList = {'AREA','RELAREA','HEIGHT','BLOWSCALE'};
else
  if ischar(ColList)
    ColList = {ColList};
  end
  if ~iscell(ColList)
    error('Expecting string or cell of strings for ColList');
  end
end

% verify that requested columns exist and find index of requested column in
% PEAKFlds
NumCol = length(ColList);
indexCol = [];
for k = 1:NumCol
  n = strcmpi(ColList{k},fieldnames(PEAKFlds)); %find if column exists
  if ~any(n) %column does not exist
    error('Column %s not found',ColList{k});
   else
   indexCol = [indexCol find(n)]; %find index of requested column
  end
end

%sort and find unique PIDs
[PID,PartIdx]     = sortrows(PID); 
[~,PartIdxIdx] = sort(PartIdx); 
uPID              = unique(PID,'rows'); %,'legacy'); %find unique PIDs in list

%check for duplicates in PID list
if length(uPID) < length(PID)
  Duplicates = 1;
else
  Duplicates = 0;
end  

%set up commmand and run it with helper script
curr_nr =size(PEAKMat,1);
if (curr_nr<STUDY.NumPkRows)
    %then current in memory PEAK Matrix is everything, no need to load
    %PEAKMat
    rows_set2read = peak_findrows_formany_instpart(PID(:,1),PID(:,2));
    command_2do  ='get_spectrum_peak_script'; %Fill up Spectrum
    peak_commhelper_script;  %exec the command 2do over PeakMat matrix file parts
else
    get_spectrum_peak_script;
end
                 
% copy data for duplicate PIDs
% this will be slow, so only run as needed
if Duplicates
  % data is in first occurance in sorted list  
  FirstPID =  find([1; any(diff(PID),2)]); %find location of 1st unique PID
  LastPID = [FirstKID(2:end)-1; length(PID)]; %find location of last PID in set of unique PIDs
  for i = 1:length(FirstPID) 
      Spectrum(FirstPID(i):LastPID(i),:) = Spectrum(FirstPID,:);
  end
end

% reorder spectral data to match original PID
Spectrum = Spectrum(PartIdxIdx,:);

return
