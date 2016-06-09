%set up commmand and run it with helper script
% ---------------------------------------------------
% this script is just the code to be run inside the loop that reads PEAKMat
% e.g. this script is run in the eval statement as the peak mat is read in 
%rows_set2read=[];   this is performed outside of script
%comm_results=[];
% ---------------------------------------------------


if Polarity == 2 %for all rows
    mzData = PEAKMat(:,[PEAKFlds.MZ indexCol]); %get mz and requested column data
    tempUKID = PEAKMat(:,[PEAKFlds.INSTID PEAKFlds.PARTID]); %get INSTID and PARTID   
else %for only selected polarity
   PolIdx = PEAKMat(:,PEAKFlds.SPECID) == Polarity;
   mzData = PEAKMat(PolIdx,[PEAKFlds.MZ indexCol]); %get mz and requested column data
   tempUKID = PEAKMat(PolIdx,[PEAKFlds.INSTID PEAKFlds.PARTID]); %for specific polarity
end

%find unique PIDs
getdiff = diff(tempUKID); 
Uidx = [true;any(getdiff,2)]; 
UKID     = tempUKID(Uidx,:); %list of unique PIDs
FirstKID = find(Uidx); %index of first PID in set of unique PIDS
LastKID = [FirstKID(2:end)-1; length(tempUKID)]; %index of last PID in set of unique PIDS
[~,MatchSIDIdx,MatchUKIDIdx] = intersect(PID,UKID,'rows'); %find index of PIDs that have requested spectra
  
%set up variables
if exist('firstChunk','var')
else
    SpecCols = size(mzData,2); %number of data columns
    Spectrum = cell(length(PID),SpecCols); %create empty cell array to hold spectra
    firstChunk = 1;
end

%copy data into the cells
for j = 1:length(MatchSIDIdx)
    PeakIdx = FirstKID(MatchUKIDIdx(j)):LastKID(MatchUKIDIdx(j)); %find index for all peaks for a certain PID
    if isempty(Spectrum{MatchSIDIdx(j),1}) %have to make sure incomplete spectra from previous Spectrum hasn't been already written
        for k = 1:SpecCols
            Spectrum{MatchSIDIdx(j),k} = double(mzData(PeakIdx,k)); % copy peak data for a certain PID into Spectrum array
        end
    else
        for k = 1:SpecCols
            Spectrum{MatchSIDIdx(j),k} = [Spectrum{MatchSIDIdx(j),k}; double(mzData(PeakIdx,k))]; %if partial spectra already there append data
        end
    end  
end

clear mzData MatchSIDIdx MatchUKIDIdx Uidx UKID FirstKID LastKID getdiff tempUKID PolIdx








