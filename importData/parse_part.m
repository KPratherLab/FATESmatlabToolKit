function parse_part(HitPartFile,MissedPartFile,PeakFile,InstID)
% PARSE_PART reads files with particle, spectrum, and peak data 
% Call as PARSE_PART(DataDef,PartFile,SpecFile,PeakFile,InstID)
% DataDef is a cell matrix which specifies the data structure.
% Data are contained in the PartFile, SpecFile and PeakFile.  These 
% files are uncommented ytables of numbers which can be read into 
% Matlab quickly.  
%
% The PartFile columns are in the order given by DataDef with 
% these additions:
% 1. the lines start with the particle serial number
% 2. times are expanded into columns for 
%    year month day hour min sec 
%
% The SpecFile columns are in the order given by DataDef with 
% one addition:
% 1. the lines start with the particle serial number
% 
% The PeakFile columns are in the order given by DataDef with
% one addition:
% 1. the lines start with the particle serial number
%
% PARSE_PART adds new data to the PART, SPEC, PEAK ytables in memory.
% It then calls SPLIT_CHUNK to save parts of these ytables to chunk
% files.

% YAADA - Software Toolkit to Analyze Single-Particle Mass Spectral Data
%
% Copyright (C) 1999-2000 California Institute of Technology
% Copyright (C) 2001-2002 Arizona State University

% Jonathan O. Allen  24 Mar 00

% Moved LastPartID to INST ytable to accomodate out-of-order InstID
% Removed LastSpecID and LastPeakID - not needed
% JOA  24 Dec 01

% Minor modifications to make it work in this version (transpose some
% arrays and such)
% Alberto Cazorla 2010-09-30
% PFR 2015-4-01  streamlined version of YAADA, 
%           the largest change is that special 'tables' are replaced
%           by structures, which are filled column wise as data is read in
%           The structures are then converted to matrices, except for peak
%           which is written out to a binary file


global INST PEAK STUDY PEAKFlds PARTidFlds PARTidMat PARTdataMat PARTdataFlds partdataNAME numFldsPARTid  numFldsPARTdata missedFlds hitFlds partDataMISSED partDataHIT missedpartColumns hitpartColumns numFldsPEAKFlds peakColumns spectraColumns spectraFlds PeakFlds
%% get inst info
%PFR change from InstRow = search(INST.InstID,'==',InstID);
InstRow = find([INST.InstID]==InstID); %structure search
PartID0 = INST(InstRow).LastPartID;  %PFR,a sturct reference instead of INST(InstRow,'LastPartID');

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

partidTMP(:,PARTidFlds.INSTID) = InstID; %add inst

%%
%add data type 2 data to PARTdata (calculated or from .pol, not read in from particle files .sem or .set)
%read in/create polarity hit data
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
    partdataTMP(:,[PARTdataFlds.HIT]) = [zeros(sizeMissed,1); HitPol(:,1)]; %add hit
end

%add ring data
if any(strcmpi('RING',partdataNAME))
    partdataTMP(:,[PARTdataFlds.RING]) = [zeros(sizeMissed,1); HitPol(:,2)]; %add polarity
end

%add position in folder if in PARTdataFlds
if any(strcmpi('POSIT',partdataNAME))
    partdataTMP(:,PARTdataFlds.POSIT) = [(1:sizeMissed)';(1:sizeHit)'];
end

%sort by time
[~,sortTime] = sort(partidTMP(:,PARTidFlds.TIME)); %get sort index
partidTMP = partidTMP(sortTime,:); %reorder using time
partdataTMP = partdataTMP(sortTime,:);
partidTMP(:,PARTidFlds.PARTID) = PartID0 + (1:numPart)'; %add PIDs
INST(InstRow).LastPartID = partidTMP(end,PARTidFlds.PARTID);  %track last PartID for the inst id
INST(InstRow).DAPartID(length(INST(InstRow).DAcalibFUNCTION)) = INST(InstRow).LastPartID; %track what size cal has been applied to which particles

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

%separate into hit and missed matrices
if any(strcmpi('HIT',partdataNAME))
    partidTMPh = partidTMP(partdataTMP(:,PARTdataFlds.HIT)>0,:);
    partdataTMPh = partdataTMP(partdataTMP(:,PARTdataFlds.HIT)>0,:);
    partidTMPm = partidTMP(partdataTMP(:,PARTdataFlds.HIT)==0,:);
    partdataTMPm = partdataTMP(partdataTMP(:,PARTdataFlds.HIT)==0,:);
    
    %write missed id data to external file
    id_fname = fullfile(STUDY.ProcDir,sprintf('PARTidMissed_%s.bin',STUDY.Name));
    STUDY.PARTidMissed_filename =id_fname;
    fiid = fopen(id_fname,'a+');
    fwrite(fiid,partidTMPm','double');
    fclose(fiid);
    
    %write missed data to external files
    data_fname =  fullfile(STUDY.ProcDir,sprintf('PARTdataMissed_%s.bin',STUDY.Name));
    STUDY.PARTdataMissed_filename =data_fname;
    fdataid = fopen(data_fname,'a+');
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
    PeakDataTMP = read_peak(PeakFile);
    PeakDataTMP = cell2mat(PeakDataTMP);
    PeakData = zeros(size(PeakDataTMP,1), numFldsPEAKFlds);
    PeakData(:,peakColumns) = PeakDataTMP(:,spectraColumns); %add data from peak file
    tmp = PeakDataTMP(:,spectraFlds.PARTID);
    if tmp(1) > 1
        tmp = tmp-tmp(1)+1;
    end
    PeakData(:,PEAKFlds.PARTID) = hitPID(tmp); %add correct PID (not just position in folder)
    PeakData(:,PEAKFlds.INSTID) = InstID*ones(size(PeakData,1),1); %add InstID
    
    %%
    % PFR: set up output file to hold
    %      peak data, leave PEAK structuer
    %    in memory to hold informatoin
    % ----------------------------
    mat_fname = fullfile(STUDY.ProcDir,sprintf('PEAKMAT_%s.bin',STUDY.Name));
    STUDY.PeakMat_filename =mat_fname;
    
    %set up 2 rows in PEAK structure to give view of start/end of this
    %part of the peakmat
    if (isempty(PEAK(1).INSTID))
        fid          =fopen(mat_fname,'w');
        prevrowcnt   =0;
        lastpki      =0;
    else
        lastpki                  =length(PEAK);
        fid                      = fopen(mat_fname,'a+');
        prevrowcnt               = PEAK(lastpki).rowcnt;
    end;
    PEAK(lastpki+1).rowcnt   = prevrowcnt+1;
    PEAK(lastpki+2).rowcnt   = prevrowcnt+size(PeakData,1);
    PEAK(lastpki+1).PARTID   = partidTMP(1,2);
    PEAK(lastpki+2).PARTID   = partidTMP(end,2);
    PEAK(lastpki+1).INSTID   = InstID;
    PEAK(lastpki+2).INSTID   = InstID;
    
    % ---------- Write out Peak matrix as a binary file,
    % note that matlab defaults to write out a matrix column wise, so that
    % the first line is the entire first column, but since we want to append
    % peak rows (in the next loop of digesting files) it is best to transpose
    % the PeakMatTemp matrix so that the first row is written out in the frist
    % line. Later, when it is read in, it is read in by sets of lines (or
    % rather by offsets in memory calculated as sets of lines).  Also, that
    % means that the 'single' format here needs to be coordinated with the
    % reading in of Peak Data (see Peakcomm helper)
    fwrite(fid,PeakData','single');
    fclose(fid);
    
    %PFR set Peak Size info
    STUDY.NumPkRows=STUDY.NumPkRows+size(PeakData,1);
%     STUDY.NumPkCols=length(PEAKFlds);
end