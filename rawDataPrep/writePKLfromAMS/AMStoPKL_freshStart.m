function AMStoPKL_freshStart(fullfileAMS, topDir,PoscalibIntercept,PoscalibSlope,NegcalibIntercept,NegcalibSlope)
% Call as AMStoPKL_freshStart(fullfileAMS, topDir,PoscalibIntercept,PoscalibSlope,NegcalibIntercept,NegcalibSlope)
% AMStoPKL_freshStart is a wrapper function that is used to write calibrated pkl text files (in m/z space) from uncalibrated ams file (in time space) using
% calibration parameters provided.  Also writes supplementary .pol files.
% These .pkl and .pol files are the file format the yaada accepts for
% atofms data.  
% topDir - the full folder path in which all the .ams files are held
% PoscalibIntercept - the intercept for the positive spectra calibration
% PoscalibSlope - the slope for the positive spectra calibration
% NegcalibIntercept - the intercept for the negative spectra calibration
% NegcalibSlope - the slope for the negative spectra calibration
% fullfileAMS - a cell array of .ams fullfile names

if nargin < 6
   error('Too few input arguments.');
elseif nargin > 6  
    error('Too many input arguments.');
end

%set up variables
tic
count = 0;
previousHitParticleNumber = -999;
NumRemoved = 0;
Ring = 0;
RealData = []; %will contain all peak data for all spectra in folder

%determine number of particles in folder
%note # of particles ~= # of ams files
numAMS = length(fullfileAMS);
counter = cell(numAMS,1);
for i = 1:numAMS
    counter{i} = fullfileAMS{i}(end-10:end-7);
end
numPart = length(unique(counter)); 

%create a tracker for the polarity of spectra for each hit particle
%1- both polarity spectra acquired
%2- pos only
%3- neg only
%4- particle is hit but polarity data not recorded for some reason
PolTracker = 4*ones(numPart,1); 
RingTracker = zeros(numPart,1); %zero if no ring detected, 1 if ring detected 

for m=1:numAMS 
%     fullfileAMS{m} %useful to uncomment this is code is failing for some
%     reason
    [HitParticleCounter, IonType, Speed,TimeStampData, LaserPower, Data] = get_spectrumAMS(fullfileAMS{m}); %extract raw data from .ams file
    if IonType == 0 %pos spectra
        %note need to use the appropriate script here if using LVN (oldTOF)
        %or SLY (newTOF) data
        [FinalPeakSave, tempRing] = PeakList_gen_LVNFinal(Data, IonType, PoscalibIntercept, PoscalibSlope,3); % outputs calibrated peak lists and baseline corrected data
    else %neg spectra
        [FinalPeakSave, tempRing] = PeakList_gen_LVNFinal(Data, IonType, NegcalibIntercept, NegcalibSlope,3); % outputs calibrated peak lists and baseline corrected data
    end
    % % info for how data is organized in FinalPeakSave
    % %             FinalPeakSave(k,1) =  Height
    % %             FinalPeakSave(k,2) =  MZ
    % %             FinalPeakSave(k,3) =  Rel Area
    % %             FinalPeakSave(k,4) =  Final Area
    % %             FinalPeakSave(k,5) =  For area graph, min data pointer
    % %             FinalPeakSave(k,6) =  For area graph, max data pointer
    % %             FinalPeakSave(k,7) =  For area graph, peak data pointer
    % %             FinalPeakSave(k,8) =  Blow scale indicator (true if
    %                                     height > 7500)
    
    if isempty(FinalPeakSave) == 1 %no peaks in spectra
        NumRemoved = NumRemoved + 1;
    else
        numPeaks = size(FinalPeakSave,1);
        RealData = [RealData; HitParticleCounter*ones(numPeaks,1) ~IonType*ones(numPeaks,1) (1:numPeaks)'  FinalPeakSave(:,[2,4,1,8,3])];        
        % % info for how data is organized in RealData
        % % RealData(k,1) = location of hit particle in folder
        % % RealData(k,2) = ion type of spectra (1 = pos, 0 = neg)
        % % RealData(k,3) = within each spectra, peak #
        % % RealData(k,4:8) = data for each peak (mz, area, height, blow
        % scale, rel ara)
        if HitParticleCounter ~= previousHitParticleNumber %if this is first spectra for this particle
            if rem(HitParticleCounter,100) == 0 % Just a status update
                disp(sprintf('processing particles %d - %d...',HitParticleCounter, (HitParticleCounter+100)))
            end
            %track polarity of spectra acquired
            if IonType == 0
                PolTracker(HitParticleCounter) = 2; 
                RingTracker(HitParticleCounter) = tempRing*2;
            else
                PolTracker(HitParticleCounter) = 3;
                RingTracker(HitParticleCounter) = tempRing*3;
            end
            %         % particles are pos/neg
            %         count = count+1; % # that are unique
            %         RealData(count).Pos = []; % initialize so that there is a value
            %         RealData(count).Neg = [];
            %         RealData(count).Time = TimeStampData;
            %         RealData(count).LaserPower = LaserPower;
            %         RealData(count).Speed = Speed;
            %         RealData(count).Filename = fullfileAMS{m};
            %
            %         previousHitParticleNumber = HitParticleCounter; %HitParticleCounter;
            %         if IonType == 0 % 0 == positive spectra, 1 = negative
            %             RealData(count).Pos = FinalPeakSave;
            %         else
            %             RealData(count).Neg = FinalPeakSave;
            %         end
        else %both pos and neg spectra acquired 
            PolTracker(HitParticleCounter) = 1;
            if RingTracker > 0 %track whether spectra ringed or not
                RingTracker(HitParticleCounter) = 1;
            end 
        end
        previousHitParticleNumber = HitParticleCounter;
    end
end

%get last folder of topDir
f = filesep; %get file separator
splitFolders = regexp(topDir,f,'split');
botFolder = splitFolders{end};

% with data compiled correctly, can now write m/z spectra to .pkl
pklfilename = fullfile(topDir, sprintf('%s.pkl',botFolder));
fid = fopen(pklfilename, 'w+');
fprintf(fid,'%u, %u, %u, %f, %f, %i, %u, %f\n', RealData');
fclose(fid);

% write data that tracks polarity and ring of spectra for hit particles to
% .pol 
polfilename = fullfile(topDir, sprintf('%s.pol',botFolder));
fid = fopen(polfilename, 'w+');
fprintf(fid,'%u, %u\n',[PolTracker RingTracker]');
fclose(fid);
% StatusReport(TotalFolders).Filename = pklfilename;
% StatusReport(TotalFolders).NumRemoved = NumRemoved;
% StatusReport(TotalFolders).Ring = Ring;
clear RealData;
ElapsedTime=toc;
disp(sprintf('Processed folder in %0.0f seconds. %d spectra removed',ElapsedTime, NumRemoved))

end