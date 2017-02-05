function init_study(StudyName_Full)
% INIT_STUDY interactively initializes FATES dataset 
% Call as init_study(StudyName) or init_study
%   Note that init study only sets up STUDY information and does NOT
%   actually open or make a fates data set (see open_study or make_study)
%   Also slashes are set to forward so that this can run
%   better on unix systems aswell as Windows
% INIT_STUDY interactively updates the global STUDY variable with the fields:
%   STUDY.Name
%   STUDY.RawDir
%   STUDY.ProcDir
%   STUDY.RawFormat  Raw spectra data format, e.g. cmsPKL2016
%   STUDY.NumPkRows
%   STUDY.missedEXT 
%   STUDY.hitEXT 
%   STUDY.spectraEXT 
%   STUDY.instEXT
%   %
% STUDY is saved in [STUDY.NameFull].mat

% FATES - Software Toolkit to Analyze Single-Particle Mass Spectral Data
% PFR 2015-4-01  streamlined version of FATES, 
% Camille M Sultana 2016 more streamlined version of FATES


global FATES STUDY runbatch

if nargin > 1
  error('too many arguments, Call as init_study(StudyName)');
end

%-------------------------------------------------
%for batch exec use a 'runbatch' flag = 1 in startup
%  see also startup_fates
runbatch = 0;
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
if (length(StudyName_Path)>0 && ~exist(StudyName_Path))
    error('Path for study data not found: %s\n',StudyName_Path);
%else  studyname path is 0 so it will be ignored
end
%check to make sure no spaces in study name
k = strfind(StudyName,' ');
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
    STUDY.ProcDir  = sprintf('%s/%s/proc',path2use,StudyName);
    STUDY.RawFormat ='cmsPKL2016';
    STUDY.NumPkRows = 0;
end

%-------------------------------------------------
%PFR for batch execution:
%-------------------------------------------------
if (runbatch==0)  %PFR, 0 means get info from user
    
    % RawFormat
    s = strrep(input('Spectra/peak file format? (cmsPKL2016...)  ','s'),'\',filesep2use);
    if ~isempty(s)
        STUDY.RawFormat = s;
    end
    
    % RawDir
    s = strrep(input('Raw data directory?  ','s'),'\',filesep2use);
    if ~isempty(s)
        STUDY.RawDir = s;
    end
    if ~exist(STUDY.RawDir,'dir');
        if ~mkdir(STUDY.RawDir);
            error('Cannot create directory %s',STUDY.RawDir);
        end
    end
    
    % Processed Data Directory
    s = strrep(input('Processed data directory?  ','s'),'\',filesep2use);
    if ~isempty(s)
        STUDY.ProcDir = s;
    end
    if ~exist(STUDY.ProcDir,'dir');
        if ~mkdir(STUDY.ProcDir);
            error('Cannot create directory %s',STUDY.ProcDir);
        end
    end
    
%     t = input('What data type is in the raw directory? \n1:ATOFMS, 2:ALABAMA 3:OTHER    ');
%     STUDY.SPMSData = t;
    
    % get file identifiers for SPMS data
    STUDY.missedEXT = input('What is the unique file identifier for the files containing missed particle data? (sem,...)  ','s');
    STUDY.hitEXT = input('What is the unique file identifier for the files containing particle data? (set,...)  ','s');
    STUDY.spectraEXT = input('What is the unique file identifier for the files containing spectra data? (pkl,...)  ','s');
    STUDY.instEXT = input('What is the unique file identifier for the files containing instrument/experiment data? (inst,...)  ','s');
%     STUDY.missedEXT = sprintf('%s%s','.',STUDY.missedEXT);
%     STUDY.hitEXT = sprintf('%s%s','.',STUDY.hitEXT);
%     STUDY.spectraEXT = sprintf('%s%s','.',STUDY.spectraEXT);
%     STUDY.instEXT = sprintf('%s%s','.',STUDY.instEXT);
        
    
end; %if runbatch=1
% ---------------------

fprintf('INFO, saving study info to %s\n',STUDY.Name);
save(STUDY.NameFull,'STUDY');

return