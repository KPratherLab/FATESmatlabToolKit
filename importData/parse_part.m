function parse_part(HitPartFile,MissedPartFile,PeakFile,InstID)
% PARSE_PART reads files with particle and spectral data 
% and populates FATES study data structures and external data files
% Call as parse_part(HitPartFile,MissedPartFile,PeakFile,InstID)
% Particle data are contained in the HitPartFile and MissedPartFile
% Spectra data is in the PeakFile
% The InstID provided will be applied to all data imported into FATES supplied in the raw data files 

global INST PEAK STUDY PEAKFlds PARTidFlds PARTidMat PARTdataMat PARTdataFlds partdataNAME numFldsPARTid  numFldsPARTdata missedFlds hitFlds partDataMISSED partDataHIT missedpartColumns hitpartColumns numFldsPEAKFlds peakColumns spectraColumns spectraFlds PeakFlds

%% get inst info
InstRow = find([INST.InstID]==InstID); %structure search
PartID0 = INST(InstRow).LastPartID;  %get last part ID for current InstID

%% get part data for missed (sem) and hit (set) particles
[missedData, hitData] = read_part(HitPartFile, MissedPartFile); %change code in this script if not using .sem and .set file formats
sizeMissed = size(missedData,1); %get size
sizeHit = size(hitData,1); %get size

%% combine PART data for hit and missed particles
numPart = sizeMissed+sizeHit; %total number of particles
partidTMP = zeros(numPart,numFldsPARTid); %set up temporary partID matrix
partdataTMP = zeros(numPart,numFldsPARTdata,'single'); %set up temporary partData matrix

%add data type1 to PARTdata
if ~isempty(missedData) && ~isempty(hitData)
    partidTMP(:,PARTidFlds.TIME) = [missedData(:,missedFlds.TIME); hitData(:,hitFlds.TIME)];
    partdataTMP(1:sizeMissed,partDataMISSED) = single(missedData(:,missedpartColumns)); 
    partdataTMP((sizeMissed+1):numPart,partDataHIT) = single(hitData(:,hitpartColumns));
elseif isempty(hitData)
    partidTMP(:,PARTidFlds.TIME) = missedData(:,missedFlds.TIME);
    partdataTMP(:,partDataMISSED) = single(missedData(:,missedpartColumns)); 
else
    partidTMP(:,PARTidFlds.TIME) = hitData(:,hitFlds.TIME);
    partdataTMP(:,partDataHIT) = single(hitData(:,hitpartColumns));
end

partidTMP(:,PARTidFlds.INSTID) = InstID; %add instID

%% add data type 2 data to PARTdata (calculated or from .pol, not read in from particle files .sem or .set)

%read in/create polarity/hit data
polfilename = sprintf('%s.pol',HitPartFile(1:end-4));
fip = fopen(polfilename);
if fip ~= -1
    HitPol = textscan(fip,'%f %f','Delimiter',',');
    HitPol = cell2mat(HitPol);
    fclose(fip);
    if size(HitPol,1) ~= sizeHit %if for some weird reason the number of particles in pol file doesn't match number it set file
        HitPol = [4*ones(sizeHit,1) zeros(sizeHit,1)];
    end
else %if no .pol file
    HitPol = [4*ones(sizeHit,1) zeros(sizeHit,1)];
end

%add hit data
if any(strcmpi('HIT',partdataNAME))
    partdataTMP(:,[PARTdataFlds.HIT]) = [zeros(sizeMissed,1); HitPol(:,1)]; 
end

%add ring data
if any(strcmpi('RING',partdataNAME))
    partdataTMP(:,[PARTdataFlds.RING]) = [zeros(sizeMissed,1); HitPol(:,2)]; 
end

%add position in folder 
if any(strcmpi('POSIT',partdataNAME))
    partdataTMP(:,PARTdataFlds.POSIT) = [(1:sizeMissed)';(1:sizeHit)'];
end

%calc da
if any(strcmpi('DA',partdataNAME))
    func2call=char(INST(InstRow).DAcalibFUNCTION{end}); %retrieve function to calculate da
    if ~isempty(func2call)
        %evaluate function using calibration parameters from inst table and
        %particle velocities
        partdataTMP(:,PARTdataFlds.DA)=feval(func2call,INST(InstRow).DAcalibPARAM{end},partdataTMP(:,PARTdataFlds.VELOCITY));   %add DA to PartMat table
    else fprintf('WARN, update da, no da calib function to call');
    end;
end

%%
%sort particles by time
[~,sortTime] = sort(partidTMP(:,PARTidFlds.TIME)); %get sort index
partidTMP = partidTMP(sortTime,:); %reorder using time
partdataTMP = partdataTMP(sortTime,:);
partidTMP(:,PARTidFlds.PARTID) = PartID0 + (1:numPart)'; %add PIDs
INST(InstRow).LastPartID = partidTMP(end,PARTidFlds.PARTID);  %track last PartID for the inst id
INST(InstRow).DAPartID(length(INST(InstRow).DAcalibFUNCTION)) = INST(InstRow).LastPartID; %track what size cal has been applied to which particles


%%

%separate into hit and missed matrices
if any(strcmpi('HIT',partdataNAME))
    partidTMPh = partidTMP(partdataTMP(:,PARTdataFlds.HIT)>0,:);
    partdataTMPh = partdataTMP(partdataTMP(:,PARTdataFlds.HIT)>0,:);
    partidTMPm = partidTMP(partdataTMP(:,PARTdataFlds.HIT)==0,:);
    partdataTMPm = partdataTMP(partdataTMP(:,PARTdataFlds.HIT)==0,:);
    
    %write missed id data to external file
    fiid = fopen(STUDY.PARTidMissed_filename,'a+');
    fwrite(fiid,partidTMPm','double');
    fclose(fiid);
    
    %write missed data to external files
    fdataid = fopen(STUDY.PARTdataMissed_filename,'a+');
    fwrite(fdataid,partdataTMPm','single');  
    fclose(fdataid);
    
else %if for some reason the hit field is not being used, save all particle info to PART table held in matlab
    partidTMPh = partidTMP;
    partdataTMPh = partdataTMP;
    partidTMPm = [];
    partdataTMPm = [];
end

PARTidMat = [PARTidMat; partidTMPh]; %add hit part to table in memory
PARTdataMat = [PARTdataMat; partdataTMPh];

%get hitPART pids
hitPosit = partdataTMPh(:,PARTdataFlds.POSIT);
hitPID = partidTMPh(:,PARTidFlds.PARTID);
[~,sortFolder] = sort(hitPosit);
hitPID = hitPID(sortFolder);

%% get peak data

if ~isempty(PeakFile)
    PeakDataTMP = read_peak(PeakFile); %read in peak file
    PeakData = zeros(size(PeakDataTMP,1), numFldsPEAKFlds);
    PeakData(:,peakColumns) = PeakDataTMP(:,spectraColumns); %add data from peak file
    tmp = PeakDataTMP(:,spectraFlds.PARTID);
    if tmp(1) > 1
        tmp = tmp-tmp(1)+1;
    end
    PeakData(:,PEAKFlds.PARTID) = hitPID(tmp); %add correct PID (not just position in folder)
    PeakData(:,PEAKFlds.INSTID) = InstID*ones(size(PeakData,1),1); %add InstID
    
    % ---------- Write out Peak matrix as a binary file,
    % note that matlab defaults to write out a matrix column wise, so that
    % the first line is the entire first column, but since we want to append
    % peak rows (in the next loop of digesting files) it is best to transpose
    % the PeakMatTemp matrix so that the first row is written out in the frist
    % line. Later, when it is read in, it is read in by sets of lines (or
    % rather by offsets in memory calculated as sets of lines).  Also, that
    % means that the 'single' format here needs to be coordinated with the
    % reading in of Peak Data (see Peakcomm helper)
    fid = fopen(STUDY.PeakMat_filename,'a+');
    fwrite(fid,PeakData','single');
    fclose(fid);
    
    %update PEAK structure in memory to hold indexing information
    %set up 2 rows in PEAK structure to give view of start/end of this
    %part of the peakmat
    if (isempty(PEAK(1).INSTID)) %if first row
        prevrowcnt   =0;
        lastpki      =0;
    else
        lastpki                  =length(PEAK);
        prevrowcnt               = PEAK(lastpki).rowcnt;
    end;
    PEAK(lastpki+1).rowcnt   = prevrowcnt+1;
    PEAK(lastpki+2).rowcnt   = prevrowcnt+size(PeakData,1);
    PEAK(lastpki+1).PARTID   = partidTMP(1,2);
    PEAK(lastpki+2).PARTID   = partidTMP(end,2);
    PEAK(lastpki+1).INSTID   = InstID;
    PEAK(lastpki+2).INSTID   = InstID;
    
    %set Peak Size info for study
    STUDY.NumPkRows=STUDY.NumPkRows+size(PeakData,1);
end
