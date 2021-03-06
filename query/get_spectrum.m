function [Spectrum] = get_spectrum(PID,Polarity,ColList)
% GET_SPECTRUM returns spectra for particles. All peaks stored in the
% external binary file are written into Spectrum without any binning.
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
% Spectrum is a nx1 cell matrix with each cell containing spectra data for
% a single particle
% Within each cell spectral data is organized in a MxK matrix with M being
% the number of peaks.  Each row contains data for a single peak. Each
% column pertains to a specific data type, with the first column always being m/z. 
%   Spectrum{i}(:,1) = MZ
%
% Additional columns contain data from columns in ColList.  ColList is
% optional and defaults to
% {'Area','RelArea','Height','BlowScale'}.  So the default matrix has the form
%   Spectrum{i}(:,1) = MZ
%   Spectrum{i}(:,2) = AREA
%   Spectrum{i}(:,3) = RELAREA
%   Spectrum{i}(:,4) = HEIGHT
%
% The matrix is empty for particles which do
% not have a spectrum of Polarity.  Note that the number of rows of each matrix
% differs among particles.

global PEAKFlds FATES PEAKMat STUDY LastPIDcompare tmpPID tmpTracker

%% verify inputs
if nargin < 2 || nargin > 3
    error('Call as [Spectrum] = get_spectrum(PID,Polarity,ColList)');
end
if (~(size(PID,2)==2))  %PID is a Nx2 matrix of instid, partids
    error('Invalid PID');
end
if ~isnumeric(Polarity) || length(Polarity) > 1
    error('Expecting scalar for Polarity');
end
if Polarity < 0 || Polarity > 2
    error('Expecting 0 (negative), 1 (positive), or 2 (both) for Polarity');
end
if ~exist('ColList','var')
    ColList = {'AREA','RELAREA','HEIGHT'};
else
    if ischar(ColList)
        ColList = {ColList};
    end
    if ~iscell(ColList)
        error('Expecting string or cell of strings for ColList');
    end
end

%% verify that requested columns exist and find index of requested column in
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

%sort and get unique PIDs 
%must have a unique PID list to be able to use
%ismemberR2012a_Get_spectrum_peak_script later.  Performing the unique
%check external to ismember saves a bunch of time
[PID,~,idxOrig] = unique(PID,'rows');

%% set up commmand and run it with helper script
curr_nr =size(PEAKMat,1);
if (curr_nr<STUDY.NumPkRows)
    %then current in memory PEAK Matrix is everything, no need to load
    %PEAKMat
    rows_set2read = peak_findrows_formany_instpart(PID(:,1),PID(:,2));
    command_2do  ='get_spectrum_peak_script'; %Fill up Spectrum
    LastPIDcompare = [0 0];
    tmpPID = PID;
    tmpTracker = 1:(size(PID,1));
    peak_commhelper_script;  %exec the command 2do over PeakMat matrix file parts
else
    get_spectrum_peak_script;
end

% reorder spectral data to match original PID
Spectrum = Spectrum(idxOrig,:);

clearvars -global LastPIDcompare tmpPID tmpTracker
return

