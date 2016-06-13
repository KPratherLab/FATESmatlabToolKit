function init_study(StudyName_Full)
% INIT_YDB interactively initializes YAADA dataset 
% Call as init_study(StudyName) or init_study
%   Note that init study only sets up STUDY information and does NOT
%   actually open or make a yaada data set (see open_ydb or make_ydb)
%   Also slashes are set to forward so that this can run
%   better on unix systems aswell as Windows
% INIT_STUDY interactively updates the global STUDY variable with the fields:
%   STUDY.Name
%   STUDY.RawDir
%   STUDY.ProcDir
%   STUDY.RawFormat  Raw data format, e.g. cmsPKL2016
%   STUDY.Gauge
%   STUDY.LastInstID
%
% STUDY is saved in [STUDY.Name].mat

% YAADA - Software Toolkit to Analyze Single-Particle Mass Spectral Data
%
% Copyright (C) 1999-2000 California Institute of Technology
% Copyright (C) 2001-2008 Arizona State University
% Copyright (C) 2008 Jonathan O. Allen

% Jonathan O. Allen  2008-02-05

% Recognizes mac, unix and pc
% Recognizes TW04 again and it is set as default
% Calls make_ydb to create the data base
% Alberto Cazorla 2010-09-17
% PFR 2015-4-01  streamlined version of FATES, 
% CMS 2016 more streamlined version of FATES


global FATES STUDY DATADEF  runbatch

if nargin > 1
  error('too many arguments, Call as init_study(StudyName)');
end

%-------------------------------------------------
%PFR  for batch exec use a 'runbatch' flag = 1 in startup
%  see also startup_fatesa
%----------------------------------------------
filesep2use='/'; 
if (runbatch==0)
  if nargin == 0
    %change slashes in case whole path is input
    StudyName = input('Study Name? ','s');
    StudyName_Path = strrep(input('Path to save study to? ','s'),'\',filesep2use);
    StudyName_Full = fullfile(StudyName_Path,StudyName); 
    StudyName_Ext = '';
  else %nargin is 1, ie StudyName_Full is the argument
    [StudyName_Path StudyName StudyName_Ext]=fileparts(StudyName_Full);     
  end  
else
  %PFR NOTE: if you are running batch, 
  % pass study name on command line to this function
  %eg init_study('2013_Jan_Intensive_Whole')
  if nargin == 0
     error('batch run missing argument, Call as init_study(StudyName)');
  else %nargin is 1
     %StudyName_Full is the argument
     [StudyName_Path StudyName StudyName_Ext]=fileparts(StudyName_Full);     
  end  
end
% now do some edit check on study name
if (length(StudyName)==0)
     error('study name empty, Call as init_study(StudyName)');
end
if (length(StudyName_Ext)>0 && strcmp(upper(StudyName_Ext),'.MAT')~=1)
    fprintf('WARNING, init study, %s extension for study name ignored \n',StudyName_Ext)
end
if (length(StudyName_Path)>0 && exist(StudyName_Path))
   fprintf('INFO, init study, will use %s for path for data\n',StudyName_Path);
elseif (length(StudyName_Path)>0 && ~exist(StudyName_Path))
    error('Path for study data not found: %s\n',StudyName_Path);
%else  studyname path is 0 so it will be ignored
end

k = strfind(StudyName,' '); %check to make sure no spaces in study name
if ~isempty(k)
  error('Expecting word for StudyName (ie no spaces)');
end

% If the study name mat file exists already load it
if exist([StudyName_Full '.mat'],'file')
  fprintf('INFO, init study, loading prexisting study mat file for %s \n',StudyName)
  load(StudyName_Full);
else  %otherwise prompt user for info
    STUDY.Name = StudyName;
    STUDY.NameFull = StudyName_Full;
   
    % set defaults
    if (length(StudyName_Path)>0)
         path2use=StudyName_Path;
    else path2use=sprintf('%s/data',FATES.fates_data);
    end
    STUDY.RawDir  = sprintf('%s/%s/raw',path2use,StudyName);
    STUDY.ProcDir  = sprintf('%s/%s/ydb',path2use,StudyName);
%     STUDY.PlotDir = sprintf('%s/%s/plot',path2use,StudyName);

    STUDY.RawFormat ='cmsPKL2016';
    STUDY.LastInstID = 0;  %PFR double(instid);
    STUDY.NumPkRows = 0;
    STUDY.Gauge = 1;
end

%-------------------------------------------------
%PFR for batch execution:
%-------------------------------------------------
if (runbatch==0)  %PFR, 0 means get info from user

% RawDir
s = strrep(input(sprintf('Raw data directory? [%s] ',STUDY.RawDir),'s'),'\',filesep2use);
if ~isempty(s)
  STUDY.RawDir = s;
end
if ~exist(STUDY.RawDir,'dir');
  if ~mkdir(STUDY.RawDir);
    error('Cannot create directory %s',STUDY.RawDir);
  end
end

% Processed Data Directory
s = strrep(input(sprintf('Processed data directory? [%s] ',STUDY.ProcDir),'s'),'\',filesep2use);
if ~isempty(s)
  STUDY.ProcDir = s;
end
if ~exist(STUDY.ProcDir,'dir');
  if ~mkdir(STUDY.ProcDir);
    error('Cannot create directory %s',STUDY.ProcDir);
  end
end

% Ask if want to use gauge board data.  If gauge board data is not
% necessary/absent, can eliminate these columns in the PARTdataMat and save
% memory
STUDY.Gauge = input('Do you want to use gauge (sizing PMT) board data? 0-no 1-yes  ');

% get file extensions for data files
STUDY.missedEXT = input('What is the file extension for the files containing missed particle data? (txt, sem, ...)  ','s');
STUDY.hitEXT = input('What is the file extension for the files containing particle data? (txt, set, ...)  ','s');
STUDY.spectraEXT = input('What is the file extension for the files containing spectra data? (txt, pkl, ...)  ','s');
STUDY.instEXT = input('What is the file extension for the files containing instrument/experiment data? (txt, inst, ...)  ','s');
STUDY.missedEXT = sprintf('%s%s','.',STUDY.missedEXT);
STUDY.hitEXT = sprintf('%s%s','.',STUDY.hitEXT);
STUDY.spectraEXT = sprintf('%s%s','.',STUDY.spectraEXT);
STUDY.instEXT = sprintf('%s%s','.',STUDY.instEXT);

end; %if runbatch=1
% ---------------------

fprintf('INFO, saving study info to %s\n',STUDY.Name);
save(StudyName_Full,'STUDY');

return