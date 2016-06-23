function make_study
% Call as MAKE_STUDY
% A study must already have been initialized by init_study and loaded
% into STUDY.  Make_ydb will use raw data files as specified in STUDY.RawDir
% and populate the fates datafiles, matrices, and structures. 
% Make_study utilizes .set (hit particle), .sem (missed particle), .pkl (spectra)
% .pol (polarity), and .inst (instrument) data files.  .pol files are
% optional
%
% PFR 2015-4-01  streamlined version of FATES, 
% CMS 2016 more streamlined version fo FATES

global STUDY procDATE INST runbatch nameSET namePKL nameSEM nameINST PARTdataFlds PEAK PEAKFlds PARTidMat PARTidFlds PARTdataMat partdataNAME missedFlds hitFlds missedNAME hitNAME partDataMISSED partDataHIT missedpartColumns hitpartColumns peakFldsNAME spectraNAME spectraColumns spectraFlds peakColumns

%% input checks
if (nargin>0)
    fprintf('INFO, make study, argument is unused \n');
    fprintf('NOTE, make study, uses study name loaded into STUDY structure\n')
end;

%PFR for runbatch, 1 is for remaking pk2files (only used if runbatch=1 in startup_fates)
makepk2_batchexecvalue=0;  

%PFR  use global study variable
StudyName = STUDY.NameFull;

%PFR added info msg
fprintf('INFO, makestudy: about to work on %s: \n ',STUDY.Name)

if ~ischar(StudyName) || isempty(StudyName)
  error('Expecting word for StudyName');
end

%%
%STUDY should have been created by init_study
load(StudyName);

%set up data files
STUDY.DataFile = fullfile(STUDY.ProcDir,sprintf('Data_%s.mat',STUDY.Name));
STUDY.PARTidMissed_filename = fullfile(STUDY.ProcDir,sprintf('PARTidMissed_%s.bin',STUDY.Name));
STUDY.PARTdataMissed_filename = fullfile(STUDY.ProcDir,sprintf('PARTdataMissed_%s.bin',STUDY.Name));
STUDY.PeakMat_filename = fullfile(STUDY.ProcDir,sprintf('PEAKMAT_%s.bin',STUDY.Name));
STUDY.datadef = fullfile(STUDY.ProcDir,sprintf('datadef_%s.mat',STUDY.Name));
STUDY.DataList = fullfile(STUDY.ProcDir,sprintf('DataList_%s.mat',STUDY.Name));
STUDY.LastInstID = 0;

%check to see if processed data already in folder
%get info check fo rbatch
if (runbatch==0)
    %check to see if any files in processed directory
    fileList = dir(fullfile(STUDY.ProcDir,sprintf('*%s*',STUDY.Name)));
    fileList = {fileList.name}';
    fileList(ismember(fileList,sprintf('%s.mat',STUDY.Name))) = [];
    if ~isempty(fileList) %found files in processed directory
        RemDataFile=input('WARNING,files found in processed directory. Do you want to overwrite database files? yes = 1, no = 0  ');
        if RemDataFile == 1 %choose to overwrite files
            fclose('all');
            delete(STUDY.DataFile, STUDY.PARTidMissed_filename, STUDY.PARTdataMissed_filename, STUDY.PeakMat_filename,fullfile, STUDY.datadef, STUDY.DataList);
        else %stop make_study if user chooses to not overwrite
            fprintf('INFO, make_study, user selected to not overwrite existing data files\n');
            fprintf('INFO, make_study, user could choose new study or archive existing data files \n');
            error('Execution stopped, make_study, user selected to stop');
        end;
    end;
else
    MakePK2    = makepk2_batchexecvalue;
    RemDataFile=1;
end; %runbatch ck  PFR 

%% set up database structure
%set up fields in study
studyFields
%determine how to read in particle files (.sem and .set) and peak files (.pkl)
dataFileFields
%match PARTdataFlds to columns from particle files
hitpartColumns = [];
missedpartColumns = [];
partDataHIT = [];
partDataMISSED = [];
for i = 1:length(partdataNAME)
    hitTMP = strcmpi(partdataNAME{i}, hitNAME);
    if any(hitTMP)
        partDataHIT = [partDataHIT PARTdataFlds.(partdataNAME{i})];
        hitpartColumns = [hitpartColumns hitFlds.(hitNAME{hitTMP})];
    end
    missedTMP = strcmpi(partdataNAME{i}, hitNAME);
    if any(missedTMP)
        partDataMISSED = [partDataMISSED PARTdataFlds.(partdataNAME{i})];
        missedpartColumns = [missedpartColumns missedFlds.(missedNAME{missedTMP})];
    end    
end
%match PEAKFlds to columns from spectra files
spectraColumns = [];
peakColumns = [];
for i = 1:length(peakFldsNAME)
    peakTMP = strcmpi(peakFldsNAME{i}, spectraNAME);
    if any(peakTMP)
        peakColumns = [peakColumns PEAKFlds.(peakFldsNAME{i})];
        spectraColumns = [spectraColumns spectraFlds.(spectraNAME{peakTMP})];
    end   
end


%% 
%find all files where data is (looking for .sem, .set, .pkl, .inst, .pol files)
fprintf('INFO, make study about to find all data files\n')
findData(STUDY.RawDir);
save(fullfile(STUDY.ProcDir,'DataList.mat'),'nameSEM','nameSET','namePKL','nameINST','procDATE');

%get particle and inst and peak data
single(PARTdataMat);
for i = 1:length(nameSET)
    [InstID,NewInst] = parse_inst(nameINST{i});
    disp(nameSET{i})
    parse_part(nameSET{i},nameSEM{i},namePKL{i},InstID);
end

%save files
save(STUDY.DataFile,'INST','PARTidMat','PARTdataMat','PARTidFlds','PARTdataFlds','PEAK','PEAKFlds','-v7.3');
save(STUDY.NameFull,'STUDY');
clearvars -global runbatch nameSET namePKL nameSEM nameINST procDATE partdataNAME missedFlds hitFlds missedNAME hitNAME partDataMISSED partDataHIT missedpartColumns hitpartColumns numFldsPARTid numFldsPARTdata3 numFldsPEAKFlds peakColumns spectraColumns spectraFlds
return



