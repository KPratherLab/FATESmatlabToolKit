function parse_part_alabama(HitPartFile,MissedPartFile,InstID)
% PARSE_PART_alabama reads files with particle and spectral data
% and populates FATES study data structures and external data files
% for the alabama data format HitParticles.txt and DetectedParticles.txt
% Call as parse_part(HitPartFile,MissedPartFile,PeakFile,InstID)
% Particle data are contained in the HitPartFile(HitParticles.txt) and
% MissedPartFile(DetectedParticles.txt). DetectedParticles.txt contains
% data on the both hit and missed particles. Spectra data is in the
% HitPartFile. The InstID provided will be applied to all data imported into FATES supplied in the raw data files
% parse_part_alabama is called by loadDATA_alabama

global INST PEAK STUDY PEAKFlds PARTidFlds PARTidMat PARTdataMat PARTdataFlds partdataNAME numFldsPARTid  numFldsPARTdata missedFlds hitFlds partDataMISSED partDataHIT missedpartColumns hitpartColumns numFldsPEAKFlds peakColumns spectraColumns spectraFlds formatHit PeakFlds hitData missedData

%% get part data for all detected particles (DetectedParticles.txt)
%DetectedParticles.txt actually contains data on both hit and missed
%particles. Redundant data in the hitFile is eventually removed from
%missedData to leave only the missed particle data.

%currently reading all in at once. Don't know if this will work out if data
%files are super large.
if ~isempty(MissedPartFile)
    fim = fopen(MissedPartFile);
    if fim ~= -1 %change code from here to fclose if using a different file format
        missedData = textscan(fim, '%s %f','Delimiter',{'\t'}); %reads in data from all detected particle file
        fclose(fim);
        missedData{:,missedFlds.TIME} = (datenum(missedData{:,missedFlds.TIME},'dd.mm.yyyy HH:MM:SS')); %change string date to num
        missedData = cell2mat(missedData); %make matrix
    else
        missedData = []; %do not change
    end
else
    missedData = []; %do not change
end

%% get part data for hit particles files (HitParticles.txt)
%during this process hit particles in missedData are eliminated
%leaving only missed particles in missedData
%only reading HitPartFile by 50000 blocks at a time to prevent taxing the
%memory too much
if ~isempty(HitPartFile)
    fit = fopen(HitPartFile);
    if fit ~= -1 %change code from here to fclose if using a different file format
        while ~feof(fit) %exit out of loop once hit end of file
            %only reading in data by blocks of 5000 lines to prevent
            %memory overload
            hitData = textscan(fit,formatHit,50000,'Delimiter','\t'); %read in data from hit particles (particle and spectra info)
            if ~isempty(hitData{1}) %if file is exactly multiple of block size, could get an empty output
                writeHitPARTandPEAKdata_alabama(InstID) %write particle and spectra data to various parts of study, also remove hit particles from missedData
            end
        end
        fclose(fit);
    else
        hitData = []; %do not change
    end
else
    hitData = []; %do not change
end

%% get current inst info
InstRow = find([INST.InstID]==InstID); %structure search
PartID0 = INST.LastPartID;  %get last part ID for current InstID

%% at this point have eliminated all hit data from missedData.
%Write missed particle data to study if there is anything left
sizeMissed = size(missedData,1); %get size
if sizeMissed > 0
    %% write missed Part data to study
    numPart = sizeMissed; %total number of particles
    partidTMP = zeros(numPart,numFldsPARTid); %set up temporary partID matrix
    partdataTMP = zeros(numPart,numFldsPARTdata,'single'); %set up temporary partData matrix
    partidTMP(:,PARTidFlds.TIME) = missedData(:,missedFlds.TIME); %add time
    partdataTMP(:,partDataMISSED) = single(missedData(:,missedpartColumns)); %add any other part data
    partidTMP(:,PARTidFlds.INSTID) = InstID; %add instID
    
    %% add data type 2 data to PARTdata (calculated or from .pol, not read in from particle files .sem or .set)
    % ALABAMA data doesn't do any of this, all of this is leftover from
    % ATOFMS
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
    
    %% add unique particle identifiers to matrices
    partidTMP(:,PARTidFlds.PARTID) = PartID0 + (1:numPart)'; %add PIDs
    INST(InstRow).LastPartID = partidTMP(end,PARTidFlds.PARTID);  %track last PartID for the inst id
    %track what size cal has been applied to which particles, used in ATOFMS analysis, not ALABAMA currently
    %   INST(InstRow).DAPartID(length(INST(InstRow).DAcalibFUNCTION)) = INST(InstRow).LastPartID;
    
    
    %% write missed data to external files
    
    %write missed id data to external file
    fiid = fopen(STUDY.PARTidMissed_filename,'a+');
    fwrite(fiid,partidTMP','double');
    fclose(fiid);
    
    %write missed data to external files
    fdataid = fopen(STUDY.PARTdataMissed_filename,'a+');
    fwrite(fdataid,partdataTMP','single');
    fclose(fdataid);
    
end
end