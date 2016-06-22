%set up commmand and run it with helper script
% ---------------------------------------------------
% this script is just the code to be run inside the loop that reads PEAKMat
% e.g. this script is run in the eval statement as the peak mat is read in 
%rows_set2read=[];   this is performed outside of script
%comm_results=[];
% ---------------------------------------------------
% this script is called by peak_commhelper_script within get_spectrum

global LastPIDcompare tmpPID tmpTracker

%retrieve desired data from PEAKMat
if Polarity == 2 %for all rows
    mzData = PEAKMat(:,[PEAKFlds.MZ indexCol]); %get mz and requested column data
    tempUKID = PEAKMat(:,[PEAKFlds.INSTID PEAKFlds.PARTID]); %get INSTID and PARTID   
else %for only selected polarity
   PolIdx = PEAKMat(:,PEAKFlds.SPECID) == Polarity;
   mzData = PEAKMat(PolIdx,[PEAKFlds.MZ indexCol]); %get mz and requested column data
   tempUKID = PEAKMat(PolIdx,[PEAKFlds.INSTID PEAKFlds.PARTID]); %for specific polarity
end

%find unique PIDs from PEAKMat data
getdiff = diff(tempUKID); 
Uidx = [true;any(getdiff,2)]; 
UKID     = tempUKID(Uidx,:); %list of unique PIDs
FirstPIDcompare = UKID(1,:);
LastKID = [FirstKID(2:end)-1; length(tempUKID)]; %index of last PID in set of unique PIDS

%the next three lines are a way of getting around the built in intersect
%function which becomes slower in a non linear function as the PID list
%gets longer due to some internal checks (unique) that we don't need to run
%as we already know that UKID has already been filtered for
%uniqueness. tmpPID does not have to be unique.  However if ismemberR2012a_Get_spectrum_peak_script becomes
%too combersome to maintain can just replace the next three lines with the
%intersect line commented out below them.  
%find index of PIDs that have requested spectra
[MatchSIDIdx,MatchUKIDIdx] = ismemberR2012a_Get_spectrum_peak_script(tmpPID,UKID,'rows'); %determine what and where UKID are in tmpPID
MatchUKIDIdx = MatchUKIDIdx(MatchSIDIdx); %get index of matches in UKID
MatchSIDIdx = find(MatchSIDIdx); %get index of matches in tmpPID
% [~,MatchSIDIdx,MatchUKIDIdx] = intersect(tmpPID,UKID,'rows'); %SAVE this, find index of PIDs that have requested spectra

%set up output variables
if exist('firstChunk','var')
else
    SpecCols = size(mzData,2); %number of data columns
    Spectrum = cell(length(PID),1); %create empty cell array to hold spectra
    firstChunk = 1;
end

if ~isempty(MatchSIDIdx)
%for the first pid have to check to make sure that there aren't already
%peaks written from previous call to binary file
PeakIdx = FirstKID(MatchUKIDIdx(1)):LastKID(MatchUKIDIdx(1)); %find index for all peaks for a certain PID
if ~isequal(LastPIDcompare, FirstPIDcompare)%make sure that this particle doesn't have peak info from previous section of binary file read in    
    Spectrum{tmpTracker(MatchSIDIdx(1)),1} = mzData(PeakIdx,:); %no previous peak info found
else
    Spectrum{tmpTracker(MatchSIDIdx(1)),1} = [Spectrum{tmpTracker(MatchSIDIdx(1)),1}; mzData(PeakIdx,:)]; %if partial spectra already there append data
end

%copy data into the cells
for j = 2:length(MatchSIDIdx)
    PeakIdx = FirstKID(MatchUKIDIdx(j)):LastKID(MatchUKIDIdx(j)); %find index for all peaks for a certain PID
    Spectrum{tmpTracker(MatchSIDIdx(j)),1} = mzData(PeakIdx,:); % copy peak data for a certain PID into Spectrum array
end

%alter tracker metrics 
LastPIDcompare = UKID(end,:); %last PID added into Spectrum 
tmpPID(MatchSIDIdx(1:(end-1)),:) = [];  %remove PIDs from comparison list where data has already been found
tmpTracker(MatchSIDIdx(1:(end-1))) = []; %remove indexes
end

clear mzData MatchSIDIdx MatchUKIDIdx Uidx UKID FirstKID LastKID getdiff tempUKID PolIdx








