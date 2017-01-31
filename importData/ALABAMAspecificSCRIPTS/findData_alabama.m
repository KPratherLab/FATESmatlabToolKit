function findData_alabama(topDir)
% finds all alabama data files to load into FATES study/structures/datafiles.
% This includes DetectedParticles and HitParticles files if available.
% Filenames are added to the global lists nameHit nameDetected
% All folders and subfolders within topDir are searched.
% findData_alabama is called by loadDATA_alabama

global nameHIT nameDET procDATE STUDY nameINST

if ~ischar(topDir)
    error('Expecting string for topDir');
end

%find all subfolders
DirList = dir(topDir);
findDir = [DirList.isdir]; 
DirNames = {DirList(findDir).name};
DirNames(ismember(DirNames,{'.','..'})) = [];

%look for data
%have to have INST file to recognize as data folder, otherwise will skip
%over
findINST = dir(sprintf('%s%s%s',topDir,'\*',STUDY.instEXT)); %find instrument file
if length(findINST) > 1
    error('%s  %s \n', 'Only one inst file allowed per folder: ', topDir); %only one inst file per folder
elseif ~isempty(findINST) %create list of data files in folder if inst file found
    findDET = dir(fullfile(topDir,STUDY.detName)); %check for all detected particles file
    findHIT = dir(fullfile(topDir,STUDY.hitName)); %check for all hit particles file
    nameINSTtmp = {findINST.name};
    if ~isempty(findDET) && isempty(findHIT) %no hit particles though detected particles found
        nameDETtmp = {findDET.name}; %get names of files
        nameHITtmp = {''};
        fprintf('%s  %s \n','Warning no HitParticle.txt file found, data likely missing: ',topDir); %this can happen if only a few particles in folder most likely
    elseif isempty(findDET) && ~isempty(findHIT) %hit particles found and no detected particle file
        nameHITtmp = {findHIT.name}; %get names of files
        nameDETtmp = {''};
        fprintf('%s  %s \n','Warning no DetectedParticle.txt file found, data likely missing: ',topDir); %this can happen if only a few particles in folder most likely
    else
        nameHITtmp = {findHIT.name}; %get names of files
        nameDETtmp = {findDET.name}; %get names of files        
    end
        %create full file names
        nameHITtmp2 = fullfile(topDir,nameHITtmp); 
        nameDETtmp2 = fullfile(topDir,nameDETtmp);
        %save files names
        nameHIT = [nameHIT; nameHITtmp2']; 
        nameDET = [nameDET; nameDETtmp2'];        
        %save inst files and date processed
        nameINSTtmp = fullfile(topDir, nameINSTtmp);
        ProcDATEtmp = now;
        nameINST = [nameINST; {nameINSTtmp}];
        procDATE = [procDATE; ProcDATEtmp];
end

%iterate into subfolders
if ~isempty(DirNames)
    %find all pkls
    for i = 1:length(DirNames)
        findData_alabama(fullfile(topDir, DirNames{i}));
    end
end 