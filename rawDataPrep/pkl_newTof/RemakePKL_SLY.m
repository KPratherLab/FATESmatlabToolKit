%% Call file to remake all pkl's in a dataset using new baseline and peak
%% finding algorithms - JC 2014
% input: Top folder path, Pos/Neg Calibration numbers

                    % % % % %Data_FilePath = 'C:\Users\Loki\Desktop\New folder\';
                    % % % % Data_FilePath = 'G:\StudyData\2016_01_CalWater2_BBY\February_26_2016\';
                    % % % % %Data_FilePath = 'G:\scripts\testdata2\raw\';
                    % % % % 
                    % % % % 
                    % % % % 
                    % % % % %% SLY CW2016 BBY 
                    % % % % % February 25, 2016
                    % % % % PoscalibIntercept = 70.51347;  
                    % % % % PoscalibSlope = 8.913774e-06;
                    % % % % NegcalibIntercept = 78.86043;
                    % % % % NegcalibSlope = 8.819549e-06;
                    % % % % %%


TotalFolders = 0;
StatusReport =[];

Top_listing=dir([sprintf('%s*', Data_FilePath)]); % gets all names in directory
if isempty(Top_listing) == 1
else
    for i=3:length(Top_listing(:,1)) % skip first two in listing directory, always '.' and '..'
        if Top_listing(i,1).isdir == 1 % is another folder
            %             disp(sprintf('Date Folder %s',Top_listing(i,1).name))
            Mid_FilePath = sprintf('%s%s\\', Data_FilePath, Top_listing(i,1).name); % this is usually the date folder
            Mid_listing=dir([sprintf('%s*', Mid_FilePath)]); % gets all names in directory
            if isempty(Mid_listing) == 1 % No folders inside
            else
                for j=3:length(Mid_listing(:,1))
                    if Mid_listing(j,1).isdir == 1 % is another folder
                        %                         disp(sprintf('Data Folder %s',Mid_listing(j,1).name))
                        Bot_FilePath = sprintf('%s%s\\', Mid_FilePath, Mid_listing(j,1).name); % this is usually the study folder
                        Bot_listing=dir([sprintf('%s*', Bot_FilePath)]); %
                        if isempty(Bot_listing) == 1 % no folders inside
                        else
                            for k=3:length(Bot_listing(:,1))
                                if Bot_listing(k,1).isdir == 1
                                    VeryBot_FilePath = sprintf('%s%s\\', Bot_FilePath, Bot_listing(k,1).name); % this is usually the letter folders
                                    VeryBot_listing=dir([sprintf('%s*', VeryBot_FilePath) '.ams']); % find spectra
                                    if isempty(VeryBot_listing) == 1
                                    else
                                        TotalFolders = TotalFolders + 1;
                                        count = 0;
%                                         adjcount = 0;
%                                         both = 0;
                                        previousHitParticleNumber = -999;
                                        NumRemoved = 0;
                                        Ring = 0;
                                        
                                        disp(sprintf('Folder %s ',VeryBot_FilePath))
                                        disp(sprintf('processing particles 1-100...'))
                                        tic;
                                        for m=1:length(VeryBot_listing(:,1)) % 1 here bc '.' and '..' are excluded by using '.ams'
                                            MS_filename = sprintf('%s%s', VeryBot_FilePath, VeryBot_listing(m,1).name);
                                            [HitParticleCounter, IonType, Speed,TimeStampData, LaserPower, Data] = get_spectrumAMS_JC(MS_filename);
                                            if IonType == 0
                                                [FinalPeakSave tempRing] = PeakList_gen_SLYFinal(Data, IonType, PoscalibIntercept, PoscalibSlope,3); % outputs peak lists and baseline corrected data
                                            else
                                                [FinalPeakSave tempRing] = PeakList_gen_SLYFinal(Data, IonType, NegcalibIntercept, NegcalibSlope,3); % outputs peak lists and baseline corrected data
                                            end
                                            Ring = Ring + tempRing; % counts the possible spectra with significant ringing
                                            % %             FinalPeakSave(k,1) =  Height
                                            % %             FinalPeakSave(k,2) =  MZ
                                            % %             FinalPeakSave(k,3) =  middle Data pointer
                                            % %             FinalPeakSave(k,4) =  Final Area
                                            % %             FinalPeakSave(k,5) =  For area graph, min data pointer
                                            % %             FinalPeakSave(k,6) =  For area graph, max data pointer         
                                            
                                            % save the ringing trigger (#
                                            % peaks < 5)
                                            if isempty(FinalPeakSave) == 1
                                                NumRemoved = NumRemoved + 1;
                                                Ring = Ring + tempRing;
                                            end
                                            
                                            if HitParticleCounter ~= previousHitParticleNumber
                                                if rem(HitParticleCounter,100) == 0 % Just a status update
                                                    disp(sprintf('processing particles %d - %d...',HitParticleCounter, (HitParticleCounter+100)))
                                                end
                                                
                                                % particles are pos/neg
                                                count = count+1; % # that are unique
                                                RealData(count).Pos = []; % initialize so that there is a value
                                                RealData(count).Neg = [];
                                                RealData(count).Time = TimeStampData;
                                                RealData(count).LaserPower = LaserPower;
                                                RealData(count).Speed = Speed;
                                                RealData(count).Filename = MS_filename;
                                                
                                                previousHitParticleNumber = HitParticleCounter; %HitParticleCounter;
                                                if IonType == 0 % 0 == positive spectra, 1 = negative
                                                    RealData(count).Pos = FinalPeakSave;
                                                else
                                                    RealData(count).Neg = FinalPeakSave;
                                                end
                                            else
%                                                 both = both+1; % # they are dual polarity
                                                if IonType == 0  % still use count as this is the true particle number
                                                    RealData(count).Pos = FinalPeakSave; % puts with the same particle
                                                    RealData(count).Filename = MS_filename; % just to match the way the old pkl was.  Positive ion filepaths would be first, Matlab arranges them second
                                                else % note: have to re-write old values in case it turns out that there is no peaks in the spectrum
                                                    RealData(count).Neg = FinalPeakSave;
                                                    RealData(count).Time = TimeStampData;
                                                end
                                            end
                                        end
                                        % with data compiled correctly, can now write to .pkl
                                        pklfilename = sprintf('%s%s.pkl', Bot_FilePath, Bot_listing(k,1).name); % puts it into the directory above letter numbers (for easy access ;) )
                                        writeToPkl_JC(pklfilename, RealData);
                                        StatusReport(TotalFolders).Filename = pklfilename;
                                        StatusReport(TotalFolders).NumRemoved = NumRemoved;
                                        StatusReport(TotalFolders).Ring = Ring;
                                        clear RealData;
                                        ElapsedTime=toc;
                                        disp(sprintf('Processed folder in %0.0f seconds. %d spectra removed',ElapsedTime, NumRemoved))
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% This is to give some indications of anamolous results from the new peak
% finding algorithm.  It outputs the number of particles in a folder with
% oddly low numbers of peaks (and particles which were removed). Typically
% this is only due to extensive ringing occuring in the mass spectra
fid = fopen(sprintf('%sStatusReport.txt',Data_FilePath), 'w');
for i=1:length(StatusReport(:))
    tempStr = strrep(StatusReport(1,i).Filename, '\', '\\');
    fprintf(fid, sprintf('%s, %d, %d', tempStr, StatusReport(i).NumRemoved, StatusReport(i).Ring));
    fprintf(fid, '\n');
end
fclose(fid);


