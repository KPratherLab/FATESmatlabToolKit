function merge_study_load(StudyName1_Full, StudyName2_Full, StudyName3_Full )
%
% Call as merge_study_load(StudyName1,2,3) or merge_study_load
%
% PFR 2015-7-01 
%    This essentially perform an initstudy and open ydb on 2 different
%    studies, . STUDY2 will be merged
%    onto study1, and the new merged study is StudyName3, 
%    The Inst and Part ID updating happens in merge_study_inst or 
%    merge_study_part
%
%    To use these functions the user needs to review the instruments
%    in study1 and study2, then decide if the data for a pair of
%    instruments needs to be merged.  There are 3 choices
%    a) merely append all inst data b/c the inst ids are all different
%         merge_study_inst(1)   %1 signifies append data
%    b) append all data but renumber inst ids in study2 b/c they overlap to
%    inst ids in study 1
%          merge_study_inst(0)  or merge_study_inst()
%          That program will renumber inst2 ids to 
%              start after the last inst1 ids           
%    c)merge data from the same instumnets in study1 and study2
%          merge_study_part({1,2,4},{1,6,3})
%             This means merge Study1,inst id 1 with Study2,inst id 1;
%                              Study1,inst id 2 with Study2,inst id 6;
%                              Study1,inst id 4 with Study2,inst id 3;
%            The final merged data will have inst ids from Study1
%

global FATES STUDY runbatch  STUDY1 STUDY2 STUDY3 PARTdataMat PARTidMat PARTdataMat1 PARTidMat1 PARTdataMat2 PARTidMat2 PARTidMiss2 INST INST1 INST2 INST3 PEAK1 PEAK2 PEAK3 PEAK PARTidFlds PARTdataFlds1 PARTdataFlds; 

%-------------------------------------------------
%PFR  for batch exec use a 'runbatch' flag = 1 in startup
%  see also startup_yaadaa
%----------------------------------------------

%---------------------------------------------
% % first,set up full path of study file names
%   note, if StudyXname_Full is only a name, the path and ext are empty 
%   and the name ends up in StudyNameX variable
% --------------------------------------------
filesep2use='/'; 
if (runbatch==0)
  if nargin ==0
    %change slashes in case whole path is input
    StudyName1_Full = strrep(input('First Study Name? ','s'),'\',filesep2use);
    [StudyName1_Path StudyName1 StudyName1_Ext]=fileparts(StudyName1_Full);
    StudyName2_Full = strrep(input('Second Study Name (will merge into first)? ','s'),'\',filesep2use);
    [StudyName2_Path StudyName2 StudyName2_Ext]=fileparts(StudyName2_Full);
    StudyName3_Full = strrep(input('New Merged Study Name? ','s'),'\',filesep2use);
    [StudyName3_Path StudyName3 StudyName3_Ext]=fileparts(StudyName3_Full);
  elseif nargin==3 %nargin is 1, ie StudyName_Full is the argument
    [StudyName1_Path StudyName1 StudyName1_Ext]=fileparts(StudyName1_Full);     
    [StudyName2_Path StudyName2 StudyName2_Ext]=fileparts(StudyName2_Full);     
    [StudyName3_Path StudyName3 StudyName3_Ext]=fileparts(StudyName3_Full);     
  else 
      error('INFO, merge study, use 3 arguments')
  end  
else
  %PFR NOTE: if you are running batch, 
  % pass study name on command line to this function
  %eg init_study('2013_Jan_Intensive_Whole')
  if nargin == 0
     error('batch run missing argument, Call as init_study(StudyName)');
  elseif nargin==3 %nargin is 1, ie StudyName_Full is the argument
    [StudyName1_Path StudyName1 StudyName1_Ext]=fileparts(StudyName1_Full);     
    [StudyName2_Path StudyName2 StudyName2_Ext]=fileparts(StudyName2_Full);     
    [StudyName3_Path StudyName3 StudyName3_Ext]=fileparts(StudyName3_Full);     
  else 
      error('INFO, merge study, use 3 arguments')
  end  
end

%no edit checks, since we assume study1 and 2 are already set up
% but, now do some edit check on study name3
if isempty(StudyName3)
     error('study name empty, Call as init_study(StudyName)');
end
if (~isempty(StudyName3_Ext) && strcmp(upper(StudyName3_Ext),'.MAT')~=1)
    fprintf('WARNING, init study, %s extension for study name ignored \n',StudyName3_Ext)
end
if (~isempty(StudyName3_Path) && exist(StudyName3_Path,'var'))
   fprintf('INFO, init study, will use %s for path for data\n',StudyName3_Path);
elseif (~isempty(StudyName3_Path)>0 && ~exist(StudyName3_Path,'var'))
    error('Path for study data not found: %s\n',StudyName3_Path);
%else  studyname path is 0 so it will be ignored
end

k = strfind(StudyName3,' '); %check to make sure no spaces in study name
if ~isempty(k)
  error('Expecting word for StudyName3 (ie no spaces)');
end

% If the study name mat file exists already load it
if exist([StudyName3 '.mat'],'file')
  fprintf('INFO, merge study,a preexisting study mat file for %s exists! \n',StudyName3);
  fprintf('      You must delete it before running merge. \n')
  error('Merge study will not overwrite preexisting study files');
end;

% -------------------------------------------------------------------------
%NOW, open study 1, put it aside, and 2, and put it aside (by renaming the data
%matrices), note some data or matrices not renamed are assumed to be the same for
% study 1 and 2 (ie FATES data structure)
% -------------------------------------------------------------------------
if exist([StudyName1 '.mat'],'file')
  fprintf('INFO, merge study, loading study2 %s \n',StudyName1)
  load(StudyName1);   %study info points to data info
  load(STUDY.DataFile);
  STUDY1    = STUDY;clear STUDY;  %NOTE, this is renaming and not using memory
  INST1     = INST; clear INST;
  PARTdataMat1  = PARTdataMat; clear PARTdataMat;
  PARTidMat1 = PARTidMat; clear PARTidMat;
  PEAK1     = PEAK; clear PEAK;
  
  if STUDY1.Gauge
      PARTdataFlds1 = PARTdataFlds;
  end

%   Pk2List1  = Pk2List;            %whee is this again?
  
else error('ERROR, merge study, study1 name not found %s \n,StudyName1');
end;

if exist([StudyName2 '.mat'],'file')
  fprintf('INFO, merge study, loading study2 %s \n',StudyName2)
  load(StudyName2);  %study info points to data info
  load(STUDY.DataFile);
  STUDY2    = STUDY;clear STUDY;  %NOTE, this is renaming and not using memory
  INST2     = INST; clear INST;
  PARTdataMat2  = PARTdataMat; clear PARTdataMat;
  PARTidMat2 = PARTidMat; clear PARTidMat;
  PEAK2     = PEAK; clear PEAK;
  %Pk2List2  = Pk2List;
  %load missed id data files
  fid = fopen(STUDY2.PARTidMissed_filename);
  PARTidMiss2 = cell2mat(textscan(fid,'%f %f %f','Delimiter',','));
  fclose(fid);
  
else error('ERROR, merge study, study2 name not found %s \n',StudyName2);
end;


% ---------------------------------------------------------------
% ---------------------------------------------------------------
%new (merged) study variables are called STUDY3
% these have to be set up, just like in initstudy
% ---------------------------------------------------------------
STUDY3.Name = StudyName3;
   
    % set defaults
if ~isempty(StudyName3_Path)
         path2use=StudyName3_Path;
    else path2use=sprintf('%s/data',FATES.yaada_data);
end
    STUDY3.RawDir  = sprintf('%s/%s/raw',path2use,StudyName3);
%     STUDY3.PK2Dir  = sprintf('%s/%s/pk2',path2use,StudyName3);
    STUDY3.ProcDir  = sprintf('%s/%s/proc',path2use,StudyName3);
%     STUDY3.PlotDir = sprintf('%s/%s/plot',path2use,StudyName3);

    STUDY3.RawFormat ='TW04';
    STUDY3.LastInstID = 0;  %PFR double(instid);
% ---------------------------
%PFR for batch execution:
%------------------------------
if (runbatch==0)  %PFR, 0 means get info from user

% ProcDir
s = strrep(input(sprintf('Processed data directory? [%s] ',STUDY3.ProcDir),'s'),'\',filesep2use);
if ~isempty(s)
  STUDY3.ProcDir = s;
end
if ~exist(STUDY3.ProcDir,'dir');
  if ~mkdir(STUDY3.ProcDir);
    error('Cannot create directory %s',STUDY3.ProcDir); 
  end
end
STUDY3.DataFile = fullfile(STUDY3.ProcDir, sprintf('Data_%s.mat',STUDY3.Name));
fid = fopen(STUDY3.DataFile,'w'); %create file for later use
fclose(fid);
STUDY3.PARTidMissed_filename = fullfile(STUDY3.ProcDir, sprintf('PARTidMissed_%s.bin',STUDY3.Name));
STUDY3.PARTdataMissed_filename = fullfile(STUDY3.ProcDir, sprintf('PARTdataMissed_%s.bin',STUDY3.Name));
STUDY3.PeakMat_filename  = fullfile(STUDY3.ProcDir, sprintf('PEAKMAT_%s.bin',STUDY3.Name));
if STUDY1.Gauge || STUDY2.Gauge
    STUDY3.Gauge = 1;
else
    STUDY3.Gauge = 0;
end

end; %if runbatch=1
% --------------------------------------- end of set up study 3 names

