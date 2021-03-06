function merge_study_load
%
% Call as merge_study_load
%
% PFR 2015-7-01
%    This essentially perform an initstudy and open ydb on 2 different
%    studies, . STUDY2 will be merged
%    onto study1, and the new merged study is StudyName3,
%    Updating the InstID and PartID in the INST structure, the PARTidMat matrix,
%    the external PARTidMissed binary file, 
%    and the PEAKmat external binary fileThe happens in merge_study_inst or
%    merge_study_part.  Merge_study_load merges the PARTdataMat matrix and 
%    the external PARTdataMissed binary file.  
%
%    NOTE: As written merge_study_load requires the PEAKFlds structures in
%    STUDY1 and STUDY2 to be identical.  This means the two external PEAKMat
%    binary files will be easily appended to each other without
%    reformatting.
%
%    To use these functions the user needs to review the instruments
%    in study1 and study2, then decide if the data for a pair of
%    instruments needs to be merged.  There are 3 choices
%    a) merely append all inst data b/c the inst ids are all different
%         merge_study_inst(0)   %0 signifies append data
%    b) append all data but renumber inst ids in study2 b/c they overlap to
%    inst ids in study 1
%          merge_study_inst(1)
%          That program will renumber inst2 ids to
%              start after the last inst1 ids
%    c)merge data from the same instumnets in study1 and study2
%          merge_study_part((1,2,4),(1,6,3))
%             This means merge Study1,inst id 1 with Study2,inst id 1;
%                              Study1,inst id 2 with Study2,inst id 6;
%                              Study1,inst id 4 with Study2,inst id 3;
%            The final merged data will have inst ids from Study1
%
% Updated for final version
% added flexibilty if PART and INST data structures don't match
% Camille Sultana 2016
%
global FATES STUDY runbatch  STUDY1 STUDY2 STUDY3 PARTdataMat PARTidMat PARTidMissed INST INST1 INST2 INST3 PEAK1 PEAK2 PEAK3 PEAK PARTidMat2 PARTidFlds PARTidFlds1 PARTidFlds2 PARTdataFlds  PEAKFlds PEAKFlds1 PEAKFlds2 num1 num2 PARTmisseddataFlds;

% -------------------------------------------------
% PFR  for batch exec use a 'runbatch' flag = 1 in startup
%  see also startup_yaadaa
% ----------------------------------------------
% 
% ---------------------------------------------
% % first,set up full path of study file names
%   note, if StudyXname_Full is only a name, the path and ext are empty
%   and the name ends up in StudyNameX variable
% --------------------------------------------

%% identify studies
filesep2use='/';
if (runbatch==0)
    if nargin ==0
        %change slashes in case whole path is input
        StudyName1_Full = strrep(input('Full Path to First Study file? ','s'),'\',filesep2use);
        [StudyName1_Path StudyName1 StudyName1_Ext]=fileparts(StudyName1_Full);
        StudyName2_Full = strrep(input('Full Path to Second Study file (will merge into first)? ','s'),'\',filesep2use);
        [StudyName2_Path StudyName2 StudyName2_Ext]=fileparts(StudyName2_Full);
        StudyName3 = input('Name for new study? ','s');
        StudyName3_Path = strrep(input('Path to save new study to? ','s'),'\',filesep2use);
        StudyName_Full3 = fullfile(StudyName3_Path,StudyName3); 
        StudyName3_Ext = '.mat';
%     elseif nargin==3 %nargin is 1, ie StudyName_Full is the argument
%         [StudyName1_Path StudyName1 StudyName1_Ext]=fileparts(StudyName1_Full);
%         [StudyName2_Path StudyName2 StudyName2_Ext]=fileparts(StudyName2_Full);
%         [StudyName3_Path StudyName3 StudyName3_Ext]=fileparts(StudyName3_Full);
    else
        error('INFO, merge study does not take arguments, use 3 arguments')
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

%% do some edit check on study name3
if isempty(StudyName3)
    error('study name empty, Call as init_study(StudyName)');
end
if (~isempty(StudyName3_Ext) && strcmp(upper(StudyName3_Ext),'.MAT')~=1)
    fprintf('WARNING, init study, %s extension for study name ignored \n',StudyName3_Ext)
end
if (~isempty(StudyName3_Path) && exist(StudyName3_Path,'dir'))
    fprintf('INFO, init study, will use %s for path for data\n',StudyName3_Path);
elseif (~isempty(StudyName3_Path)>0 && ~exist(StudyName3_Path,'dir'))
    error('Path for study data not found: %s\n',StudyName3_Path);
    %else  studyname path is 0 so it will be ignored
end

k = strfind(StudyName3,' '); %check to make sure no spaces in study name
if ~isempty(k)
    error('Expecting word for StudyName3 (ie no spaces)');
end

STUDY3.Name = StudyName3;
STUDY3.NameFull = StudyName_Full3;

% STUDY 3 Processed Data Directory
s = strrep(input('Processed data directory?  ','s'),'\',filesep2use);
if ~isempty(s)
    STUDY3.ProcDir = s;
end
if ~exist(STUDY3.ProcDir,'dir')
    if ~mkdir(STUDY3.ProcDir)
        error('Cannot create directory %s',STUDY.ProcDir);
    end
end

%check to see if any study files in processed directory
fileList = dir(fullfile(STUDY3.ProcDir,sprintf('*%s*',STUDY3.Name)));
if ~isempty(fileList)
    %found files in processed directory
    fprintf('INFO, merge study,a preexisting study files for %s exists! \n',STUDY3.Name);
    fprintf('      You must delete them before running merge. \n')
    error('Merge study will not overwrite preexisting study files');
end

%% load pre existing studies
% -------------------------------------------------------------------------
%NOW, open study 1, put it aside, and 2, and put it aside (by renaming the data
%matrices), note some data or matrices not renamed are assumed to be the same for
% study 1 and 2 (ie FATES data structure)
% -------------------------------------------------------------------------
if exist([StudyName1_Full '.mat'],'file')
    fprintf('INFO, merge study, loading study2 %s \n',StudyName1)
    load(StudyName1_Full);   %study info points to data info
    open_study
    STUDY1    = STUDY;  %NOTE, this is renaming and not using memory
    INST1     = INST; 
    PARTdataMat1  = PARTdataMat;
    PARTidMat1 = PARTidMat; 
    PARTdataFlds1 = PARTdataFlds; 
    PARTidFlds1 = PARTidFlds; 
    PEAK1     = PEAK; 
    PEAKFlds1 = PEAKFlds;      
    PARTmisseddataFlds1 = PARTmisseddataFlds;
else error('ERROR, merge study, study1 name not found %s \n,StudyName1');
end;

if exist([StudyName2_Full '.mat'],'file')
    fprintf('INFO, merge study, loading study2 %s \n',StudyName2)
    load(StudyName2_Full);  %study info points to data info
    open_study
    STUDY2    = STUDY; STUDY = [];  %NOTE, this is renaming and not using memory
    INST2     = INST; INST = [];
    PARTdataMat2  = PARTdataMat; PARTdataMat = [];
    PARTidMat2 = PARTidMat; PARTidMat = [];
    PARTdataFlds2 = PARTdataFlds; PARTdataFlds = [];
    PARTidFlds2 = PARTidFlds; PARTidFlds = [];
    PEAK2     = PEAK; PEAK = [];
    PEAKFlds2 = PEAKFlds; PEAKFlds = [];
    PARTmisseddataFlds2 = PARTmisseddataFlds; PARTmisseddataFlds = [];
else error('ERROR, merge study, study2 name not found %s \n',StudyName2);
end;

num1 = size(PARTidMat1,1);
num2 = size(PARTidMat2,1);

%% set up study 3 structure

%set up file names and directories
STUDY3.DataFile = fullfile(STUDY3.ProcDir, sprintf('Data_%s.mat',STUDY3.Name));
fid = fopen(STUDY3.DataFile,'w'); %create file for later use
fclose(fid);
STUDY3.PARTidMissed_filename = fullfile(STUDY3.ProcDir, sprintf('PARTidMissed_%s.bin',STUDY3.Name));
STUDY3.PARTdataMissed_filename = fullfile(STUDY3.ProcDir, sprintf('PARTdataMissed_%s.bin',STUDY3.Name));
STUDY3.PeakMat_filename  = fullfile(STUDY3.ProcDir, sprintf('PEAKMAT_%s.bin',STUDY3.Name));

%add raw dirs data to new study
STUDY3.RawDir = {STUDY1.RawDir STUDY2.RawDir};

%set datadef file to study1
STUDY3.datadef = fullfile(STUDY3.ProcDir, sprintf('datadef_%s.mat',STUDY3.Name)); 
copyfile(STUDY1.datadef, STUDY3.datadef);

%combine raw datalists
STUDY3.DataList = fullfile(STUDY3.ProcDir, sprintf('DataList_%s.mat',STUDY3.Name)); 
load(STUDY1.DataList);
nameSEM1 = nameSEM;
nameSET1 = nameSET;
namePKL1 = namePKL;
procDATE1 = procDATE;
load(STUDY2.DataList);
nameSEM = [nameSEM1; nameSEM];
nameSET = [nameSET1; nameSET];
namePKL = [namePKL1; namePKL];
procDATE = [procDATE1; procDATE];
save(STUDY3.DataList,'nameSEM','nameSET','namePKL','nameINST','procDATE');

%% need to do checks to make sure both studies are compatible

%check PEAK
if ~isequal(PEAKFlds1, PEAKFlds2)
    error('%s\n%s\n%s','PEAKFlds must be the same for the two studies to merge the external binary spectra file.',...
        'The number of columns of each PEAKMat table must be equal.',...
        'The data type stored in each column of PEAKMat must match between studies');
else
    PEAKFlds = PEAKFlds1;
    %copy PeakMat data from STUDY1 which will remain unchanged
    copyfile(STUDY1.PeakMat_filename, STUDY3.PeakMat_filename);
end

%update PARTdataMat
if isequal(PARTdataFlds1,PARTdataFlds2) %check to see if fields are equal
    %update hit particle
    PARTdataMat = [PARTdataMat1; PARTdataMat2]; %combine if equal
    PARTdataFlds = PARTdataFlds1;
else %pad then combine if not equal
    warning('%s\n%s\n','Warning: PARTdataFlds between the two studies do not match.',...
        'The new PARTdataMat matrix will be padded or reorganized to account for this mismatch');
    warning('%s\n%s\n','Warning: The datadef file for STUDY 3 is simply a copy of STUDY 1.',...
        'User will need to manually alter STUDY3 datadef file for it to accurately detail database structure');    
    %compare fields and get locations for placing data
    PDATAidx = cell(2,2);
    origFlds1name = fieldnames(PARTdataFlds1);
    origFlds2name = fieldnames(PARTdataFlds2);    
    PARTdataFlds1name = lower(fieldnames(PARTdataFlds1));
    PARTdataFlds2name = lower(fieldnames(PARTdataFlds2));
    PARTdataFlds3name = union(PARTdataFlds1name, PARTdataFlds2name,'stable');
    [~,PDATAidx{1,2}] = ismember(PARTdataFlds1name, PARTdataFlds3name);
    [~,PDATAidx{2,2}] = ismember(PARTdataFlds2name, PARTdataFlds3name);
    PDATAidx{1,2}(PDATAidx{1,2} == 0) = [];
    PDATAidx{2,2}(PDATAidx{2,2} == 0) = [];
    PDATAidx{1,1} = zeros(1,length(origFlds1name));
    PDATAidx{2,1} = zeros(1,length(origFlds2name));    
    for i = 1:length(PARTdataFlds1name)
        PDATAidx{1,1}(i) = PARTdataFlds1.(origFlds1name{i});
    end
    for i = 1:length(PARTdataFlds2name)
        PDATAidx{2,1}(i) = PARTdataFlds2.(origFlds2name{i});
    end
    
    %create PARTdataFlds structure for new study
        %reset PARTdataFlds to uppercase or original
    PARTdataFlds3name([1:length(PARTdataFlds1name)]) = origFlds1name;
    for i = 1:length(PARTdataFlds3name)
        PARTdataFlds.(PARTdataFlds3name{i}) = i;
    end
    
    %do padding for hit particles and combine into new PARTdataMat matrix
    numFlds = length(PARTdataFlds3name);
    PARTdataMat = single(nan(num1+num2,numFlds));
    PARTdataMat(1:num1,PDATAidx{1,2}) = PARTdataMat1(:,PDATAidx{1,1});
    PARTdataMat((num1+1):end,PDATAidx{2,2}) = PARTdataMat2(:,PDATAidx{2,1});
end
clear PARTdataMat1 PARTdataMat2 PARTdataFlds1 PARTdataFlds2

%update missed particle data
if isequal(PARTmisseddataFlds1,PARTmisseddataFlds2) %check to see if fields are equal   
    %note this is probably not the fastest way
    copyfile(STUDY1.PARTdataMissed_filename, STUDY3.PARTdataMissed_filename);
    PARTmisseddataFlds = PARTmisseddataFlds2;
    [~,PARTdataMissedTMP] = load_missed('notpath', STUDY2.PARTdataMissed_filename);
    fdataid = fopen(STUDY3.PARTdataMissed_filename,'a');
    fwrite(fdataid,PARTdataMissedTMP','single');
    fclose(fdataid);
    clear PARTdataMissedTMP
else
    warning('%s\n%s\n','Warning: PARTmisseddataFlds between the two studies do not match.',...
        'The new missed PARTdataMat matrix will be padded or reorganized to account for this mismatch');
    warning('%s\n%s\n','Warning: The datadef file for STUDY 3 is simply a copy of STUDY 1.',...
        'User will need to manually alter STUDY3 datadef file for it to accurately detail database structure');
    %compare fields and get locations for placing data
    PDATAidx = cell(2,2);
    origFlds1name = fieldnames(PARTmisseddataFlds1);
    origFlds2name = fieldnames(PARTmisseddataFlds2);    
    PARTdataFlds1name = lower(fieldnames(PARTmisseddataFlds1));
    PARTdataFlds2name = lower(fieldnames(PARTmisseddataFlds2));
    PARTdataFlds3name = union(PARTdataFlds1name, PARTdataFlds2name,'stable');
    [~,PDATAidx{1,2}] = ismember(PARTdataFlds1name, PARTdataFlds3name);
    [~,PDATAidx{2,2}] = ismember(PARTdataFlds2name, PARTdataFlds3name);
    PDATAidx{1,2}(PDATAidx{1,2} == 0) = [];
    PDATAidx{2,2}(PDATAidx{2,2} == 0) = [];
    PDATAidx{1,1} = zeros(1,length(origFlds1name));
    PDATAidx{2,1} = zeros(1,length(origFlds2name));    
    for i = 1:length(PARTdataFlds1name)
        PDATAidx{1,1}(i) = PARTdataFlds1.(origFlds1name{i});
    end
    for i = 1:length(PARTdataFlds2name)
        PDATAidx{2,1}(i) = PARTdataFlds2.(origFlds2name{i});
    end
    
    %load in missed particle data and combine with padding
    numFlds = length(PARTdataFlds3name);
    PARTmisseddataFlds = PARTmisseddataFlds1;
    [~,PARTdataMissedTMP] = load_missed('notpath', STUDY1.PARTdataMissed_filename);
    num1m = size(PARTdataMissedTMP,1);
    PARTdataMissed = single(nan(num1m,numFlds));
    PARTdataMissed(:,PDATAidx{1,2}) = PARTdataMissedTMP(:,PDATAidx{1,1});
    PARTmisseddataFlds = PARTmisseddataFlds2;
    [~,PARTdataMissedTMP] = load_missed('notpath', STUDY2.PARTdataMissed_filename);
    PARTmisseddataFlds = [];
    num2m = size(PARTdataMissedTMP,1);
    PARTdataMissed((num1m+1):(num1m+num2m),PDATAidx{2,2}) = PARTdataMissedTMP(:,PDATAidx{2,1});
    
    %create PARTdataFlds structure for new study
    %reset PARTdataFlds to uppercase or original
    PARTdataFlds3name([1:length(PARTdataFlds1name)]) = origFlds1name;
    for i = 1:length(PARTdataFlds3name)
        PARTmisseddataFlds.(PARTdataFlds3name{i}) = i;
    end
    
    %write to external file
    fdataid = fopen(STUDY3.PARTdataMissed_filename,'w');
    fwrite(fdataid,PARTdataMissed','single');
    fclose(fdataid);
    clear PARTdataMissed PARTdataMissedTMP PARTmisseddataFlds2 PARTmisseddataFlds1
end

%check PARTidMat structure
if ~isequal(PARTidFlds1,PARTidFlds2)
    fprintf('%s\n%s\n','Warning: PARTidFlds between the two studies do not match.',...
        'The new PARTidMat matrix will be padded or reorganized to account for this mismatch');
    fprintf('%s\n%s\n','Warning: The datadef file for STUDY 3 is simply a copy of STUDY 1.',...
        'User will need to manually alter STUDY3 datadef file for it to accurately detail database structure');  
    %compare fields and get locations
    PIDidx = cell(2,2);
    PARTidFlds1name = lower(fieldnames(PARTidFlds1));
    PARTidFlds2name = lower(fieldnames(PARTidFlds2));
    PARTidFlds3name = union(PARTidFlds1name, PARTidFlds2name,'stable');
    [~,PIDidx{1,2}] = ismember(PARTidFlds1name, PARTidFlds3name);
    [~,PIDidx{2,2}] = ismember(PARTidFlds2name, PARTidFlds3name);
    PIDidx{1,2}(PDATAidx{1,2} == 0) = [];
    PIDidx{2,2}(PDATAidx{2,2} == 0) = [];
    for i = 1:length(PARTidFlds1name)
        PIDidx{1,1} = PARTidFlds1.(PARTidFlds1name{i});
    end
    for i = 1:length(PARTidFlds2name)
        PIDidx{2,1} = PARTidFlds2.(PARTidFlds2name{i});
    end 
    
    %set up PARTidMat variables (will have to change ids later)
    PARTidMat = nan(num1+num2,length(PARTidFlds3name));
    PARTidMat(1:num1,PIDidx{1,2}) = PARTidMat1(:,PIDidx{1,1});
    PARTidMat((num1+1):end,PIDidx{2,2}) = PARTidMat2(:,PIDidx{2,1});
    
    %create new fields structure
    %reset PARTidFlds to uppercase or original
    PARTidFlds3name([1:length(PARTidFlds1name)]) = PARTidFlds1name;
    for i = 1:length(PARTidFlds3name)
        PARTidFlds.(PARTidFlds3name{i}) = i;
    end
    
    %set up PARTidMissed external files
    PARTidMissedTMP = load_missed(STUDY1.PARTidMissed_filename,'notpath');    
    PARTidMissed = nan(size(PARTidMissedTMP,1),length(PARTidFlds3name));
    PARTidMissed(:,PIDidx{1,2}) = PARTidMissedTMP(:,PIDidx{1,1});
    fdataid = fopen(STUDY3.PARTidMissed_filename,'w');
    fwrite(fdataid,PARTidMissed','double');  
    fclose(fdataid);    
    
    %load in study 2 but have to wait til later to rewrite ids
    PARTidMissedTMP = load_missed(STUDY2.PARTidMissed_filename,'notpath');    
    PARTidMissed = nan(size(PARTidMissedTMP,1),length(PARTidFlds3name));
    PARTidMissed(:,PIDidx{2,2}) = PARTidMissedTMP(:,PIDidx{2,1});
    clear PARTidMissedTMP    
else
    PARTidFlds = PARTidFlds1;
    PARTidMat = [PARTidMat1; PARTidMat2];
    copyfile(STUDY1.PARTidMissed_filename, STUDY3.PARTidMissed_filename);
    PARTidMissed = load_missed(STUDY2.PARTidMissed_filename,'notpath');    
end

%check INST
INSTf1 = sort(fieldnames(INST1));
INSTf2 = sort(fieldnames(INST2));
INSTfL1 = lower(INSTf1);
INSTfL2 = lower(INSTf2);
if ~isequal(INSTf1,INSTf2)
    fprintf('%s\n%s\n','Warning: INST field names between the two studies do not match.',...
        'The new INST structure will be padded or reorganized to account for this mismatch');
    fprintf('%s\n%s\n','Warning: The datadef file for STUDY 3 is simply a copy of STUDY 1.',...
        'User will need to manually alter STUDY3 datadef file for it to accurately detail database structure');  
    %find matching field names except for capitalization
    [test3,test4] = ismember(INSTfL2, INSTfL1); %compare lowercase
    [test7,test8] = ismember(INSTf2,INSTf1); %compare original
    fixCap2 = test3 & ~test7;
    fixCap1 = test4(fixCap2);
    fixCap2 = INSTf2(fixCap2); %fields to replace
    fixCap1 = INSTf1(fixCap1); %new field names
    %replace INST2 with INST1 name format
    for i = 1:length(fixCap2)
        [INST2.(fixCap1{i})] = INST2.(fixCap2{i});
        INST2 = rmfield(INST2,fixCap2{i});
    end
    
    %ok everything should match now unless it has an actually different
    %field name
    %find field names unique to each inst
    INSTf2 = fieldnames(INST2);
    INSTf3 = union(INSTf1,INSTf2);
    addF2 = setdiff(INSTf3,INSTf2);
    addF1 = setdiff(INSTf3,INSTf1);
    %add missing field names to each inst
    for i = 1:length(addF2)
        INST2().(addF2{i}) = [];
    end
    for i = 1:length(addF1)
        INST1().(addF1{i}) = [];
    end
end
% 
save(STUDY3.NameFull,'STUDY');
disp('NEW STUDY INITIALIZATION IS NOT COMPLETE');
disp('now run merge_study_part or merge_study_inst to complete!');
