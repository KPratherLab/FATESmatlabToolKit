function make_study
% Call as MAKE_STUDY
% A study must already have been initialized by init_study and loaded
% into STUDY.  Make_ydb will use raw data files as specified in STUDY.RawDir
% and populate the yaada datafiles, matrices, and structures. 
% Make_ydb utilizes .set (hit particle), .sem (missed particle), .pkl (spectra)
% .pol (polarity), and .inst (instrument) data files.  .pol files are
% optional
%
% YAADA - Software Toolkit to Analyze Single-Particle Mass Spectral Data
%
% Copyright (C) 1999-2000 California Institute of Technology
% Copyright (C) 2001-2008 Arizona State University
% Copyright (C) 2008 Jonathan O. Allen

% Jonathan O. Allen  2008-02-06

% Recognizes TW04 and call the corresponding digest function
% Alberto Cazorla 2010-09-17
% PFR 2015-4-01  streamlined version of FATES, 
% CMS 2016 more streamlined version fo FATES

global STUDY procDATE INST runbatch nameSET namePKL nameSEM nameINST PARTdataFlds PEAK PEAKFlds PARTidMat PARTidFlds PARTdataMat partdataNAME missedFlds hitFlds missedNAME hitNAME partDataMISSED partDataHIT missedpartColumns hitpartColumns peakFldsNAME spectraNAME spectraColumns spectraFlds peakColumns

if (nargin>0)
    fprintf('INFO, make study, argument is unused \n');
    fprintf('NOTE, make study, uses study name given in init_study script\n')
end;

%PFR for runbatch, 1 is for remaking pk2files (only used if runbatch=1 in startup_fates)
makepk2_batchexecvalue=0;  

%PFR  use global study variable
StudyName = STUDY.NameFull;

%PFR added info msg
fprintf('INFO, makestudy: about to work on %s: \n ',StudyName)

if ~ischar(StudyName) || isempty(StudyName)
  error('Expecting word for StudyName');
end

%STUDY should have been created by init_study
load(StudyName);

%set up data file
STUDY.DataFile = sprintf('%s/Data_%s.mat',STUDY.ProcDir,STUDY.Name);

%PFR get info check fo rbatch
if (runbatch==0)  
  if (exist(STUDY.DataFile,'file'))
      RemDataFile=input(sprintf('WARNING,previous data files,%s, and %s \n,       will be OVERWRITTEN, OK? (y/n) [y]',STUDY.DataFile,STUDY.PeakMat_filename),'s');
      if isempty(RemDataFile) 
          RemDataFile=1;
      end;
      if bool2num(RemDataFile) == 1
          fclose('all');
          delete(STUDY.DataFile, STUDY.PARTidMissed_filename, STUDY.PARTdataMissed_filename, STUDY.PeakMat_filename);
      else bool2num(RemDataFile)==0
          fprintf('INFO, make_study, user selected to not overwrite existing data files\n');
          fprintf('INFO, make_study, user could choose new study or archive existing data files \n');
          error('Execution stopped, make_study, user selected to stop');
      end;
  else RemDataFile=1; %go ahead and continue making new study 
  end;
else
  MakePK2    = makepk2_batchexecvalue;
  RemDataFile=1;
end; %runbatch ck  PFR

STUDY.LastInstID = 0; %PFR 
 
%PFR data def is where struct tables were defined first. 

%set up fields in study
studyFields
%determine how to read in particle files (.sem and .set) and peak files
%(.pkl)
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

%find all files where data is (looking for .sem, .set, .pkl, .inst, .pol
%files)
fprintf('INFO, make study about to find all data files\n')
findData(STUDY.RawDir);
save(fullfile(STUDY.ProcDir,'DataList.mat'),'nameSEM','nameSET','namePKL','nameINST','procDATE');

%get particle and inst data
single(PARTdataMat);
for i = 1:length(nameSET)
    [InstID,NewInst] = parse_inst(nameINST{i});
    nameSET{i}
    parse_part(nameSET{i},nameSEM{i},namePKL{i},InstID);
end

save(STUDY.DataFile,'INST','PARTidMat','PARTdataMat','PARTidFlds','PARTdataFlds','PEAK','PEAKFlds','-v7.3');
save(STUDY.NameFull,'STUDY');
clearvars -global runbatch nameSET namePKL nameSEM nameINST procDATE partdataNAME missedFlds hitFlds missedNAME hitNAME partDataMISSED partDataHIT missedpartColumns hitpartColumns numFldsPARTid numFldsPARTdata3 numFldsPEAKFlds peakColumns spectraColumns spectraFlds
return



