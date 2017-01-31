function writeHitPARTandPEAKdata_alabama(InstID)
% Call as writeHitPARTandPEAKdata_alabama(InstID)
% InstID is an instrument/experiment generated via parse_inst to help
% identify particle and peak data.  
% writeHitPARTandPEAKdata_alabama parses/writes hit particle data to PARTdataMat and
% PARTidMat of a FATES STUDY as well as PEAK data to the external peakmat
% file. writeHitPARTandPEAKdata_alabama is called by parse_part_alabama.

global hitData missedData hitFlds numFldsPARTid numFldsPARTdata PARTidFlds INST PEAKFlds STUDY PEAK PARTidMat PARTdataMat partDataHIT MZValue Polarity hitpartColumns numFldsPEAKFlds

%% reorganize hitData a little bit
hitData{:,hitFlds.TIME} = (datenum(hitData{:,hitFlds.TIME},'dd.mm.yyyy HH:MM:SS')); %change date string to num
PeakDataTMP = cell2mat(hitData(3:end)); %make matrix of spectra data
hitData = cell2mat(hitData(1:2)); %remove specra data from part data

%% remove hit particles from detected particles list
if ~isempty(missedData)
    missedData = setdiff(missedData, hitData, 'rows');
end

%% write partData for hitData to study
sizeHit = size(hitData,1); %get size
numPart = sizeHit; %total number of particles
partidTMP = zeros(numPart,numFldsPARTid); %set up temporary partID matrix
partdataTMP = zeros(numPart,numFldsPARTdata,'single'); %set up temporary partData matrix

%add data type1 to PARTdata
partidTMP(:,PARTidFlds.TIME) = hitData(:,hitFlds.TIME);
partdataTMP(:,partDataHIT) = single(hitData(:,hitpartColumns));
partidTMP(:,PARTidFlds.INSTID) = InstID; %add instID

%% add data type 2 data to PARTdata (calculated or from .pol, not read in from particle files .sem or .set)
% this is all from ATOFMS stuff, not applicable to ALABAMA currently
% %read in/create polarity/hit data
% polfilename = sprintf('%s.pol',HitPartFile(1:end-4));
% fip = fopen(polfilename);
% if fip ~= -1
%     HitPol = textscan(fip,'%f %f','Delimiter',',');
%     HitPol = cell2mat(HitPol);
%     fclose(fip);
%     if size(HitPol,1) ~= sizeHit %if for some weird reason the number of particles in pol file doesn't match number it set file
%         HitPol = [4*ones(sizeHit,1) zeros(sizeHit,1)];
%     end
% else %if no .pol file
%     HitPol = [4*ones(sizeHit,1) zeros(sizeHit,1)];
% end
% 
% %add hit data
% if any(strcmpi('HIT',partdataNAME))
%     partdataTMP(:,[PARTdataFlds.HIT]) = [zeros(sizeMissed,1); HitPol(:,1)]; 
% end
% 
% %add ring data
% if any(strcmpi('RING',partdataNAME))
%     partdataTMP(:,[PARTdataFlds.RING]) = [zeros(sizeMissed,1); HitPol(:,2)]; 
% end
% 
% %add position in folder 
% if any(strcmpi('POSIT',partdataNAME))
%     partdataTMP(:,PARTdataFlds.POSIT) = [(1:sizeMissed)';(1:sizeHit)'];
% end
% 
% %calc da
% if any(strcmpi('DA',partdataNAME))
%     func2call=char(INST(InstRow).DAcalibFUNCTION{end}); %retrieve function to calculate da
%     if ~isempty(func2call)
%         %evaluate function using calibration parameters from inst table and
%         %particle velocities
%         partdataTMP(:,PARTdataFlds.DA)=feval(func2call,INST(InstRow).DAcalibPARAM{end},partdataTMP(:,PARTdataFlds.VELOCITY));   %add DA to PartMat table
%     else fprintf('WARN, update da, no da calib function to call');
%     end;
% end

%% add unique particle identifier info
PartID0 = INST.LastPartID;  %get last part ID for current InstID
InstRow = find([INST.InstID]==InstID); %structure search
partidTMP(:,PARTidFlds.PARTID) = PartID0 + (1:numPart)'; %add PIDs
INST(InstRow).LastPartID = partidTMP(end,PARTidFlds.PARTID);  %track last PartID for the inst id
%track what size cal has been applied to which particles, right now not
%used for ALABAMA
% INST(InstRow).DAPartID(length(INST(InstRow).DAcalibFUNCTION)) = INST(InstRow).LastPartID; 


%% add data to FATES structures
PARTidMat = [PARTidMat; partidTMP]; %add hit part to table in memory
PARTdataMat = [PARTdataMat; partdataTMP];

%get hitPART pids
% hitPosit = partdataTMPh(:,PARTdataFlds.POSIT);
% hitPID = partidTMPh(:,PARTidFlds.PARTID);
% [~,sortFolder] = sort(hitPosit);
% hitPID = hitPID(sortFolder);

%% get peak data
if ~isempty(PeakDataTMP)
    numMZ = length(MZValue);
    PeakData = single(zeros(numMZ*sizeHit, numFldsPEAKFlds));
    
    %% reorganize/create peak data from HitParticle file
    %peak intensity values
    PeakDataTMP = PeakDataTMP';
    PeakData(:,PEAKFlds.AREA) = PeakDataTMP(:);
    
    %mz values
    PeakData(:,PEAKFlds.MZ) = repmat(MZValue,sizeHit,1);
    
    %spectra polarity values
    PeakData(:,PEAKFlds.SPECID) = repmat(Polarity,sizeHit,1);
    
    %particle index values
    tmp = single(1:sizeHit);
    tmp = repmat(tmp,numMZ,1);
    tmp = tmp(:);
    %set partID and instID
    hitPID = partidTMP(:,PARTidFlds.PARTID);
    PeakData(:,PEAKFlds.PARTID) = hitPID(tmp); %add correct PID (not just position in folder)
    PeakData(:,PEAKFlds.INSTID) = InstID*ones(size(PeakData,1),1); %add InstID
    
    %% ---------- Write out Peak matrix as a binary file,
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
    
    %% update PEAK structure in memory to hold indexing information
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