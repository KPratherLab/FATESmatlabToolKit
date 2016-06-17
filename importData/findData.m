function findData(topDir)
% finds all atofms data files to load into YAADA study/structures/datafiles.
% This includes SET, SEM, PKL,and INST files and POL files if available.
% Filenames are added to the global lists nameSET, nameSEM, namePKL, nameINST .
% All folders and subfolders within topDir are searched.
% STUDY.missedEXT = input('What is the file extension for the files containing missed particle data? (txt, sem, ...)','s');
% STUDY.hitEXT = input('What is the file extension for the files containing particle data? (txt, set, ...)','s');
% STUDY.spectraEXT = input('What is the file extension for the files containing spectra data? (txt, pkl, ...)','s');
% STUDY.instEXT = input('What is the file extension for the files containing instrument/experiment data? (txt, inst, ...)','s');
global nameSET nameSEM namePKL nameINST procDATE STUDY

if ~ischar(topDir)
    error('Expecting string for topDir');
end

%find all subfolders
DirList = dir(topDir);
findDir = [DirList.isdir]; 
DirNames = {DirList(findDir).name};
DirNames(ismember(DirNames,{'.','..'})) = [];
%file type if this is not usually the case.

% if ~isempty(findSEM)
findINST = dir(sprintf('%s%s%s',topDir,'\*',STUDY.instEXT)); %find instrument file
%     if isempty(findINST)
%         error('%s  %s\n', 'No inst file found',topDir); %all folders with data must contain an instrument file
%     else
if length(findINST) > 1
    error('%s  %s\n', 'Only one inst file allowed per folder', topDir); %only one inst file per folder
elseif ~isempty(findINST) %create list of data files in folder %found an inst file
    findSEM = dir(sprintf('%s%s%s',topDir,'\*',STUDY.missedEXT)); %check for sem (missed data) file
    nameINSTtmp = {findINST.name};
    findSET = dir(sprintf('%s%s%s',topDir,'\*',STUDY.hitEXT));
    findPKL = dir(sprintf('%s%s%s',topDir,'\*',STUDY.spectraEXT'));
    if isempty(findSET) %no hit particles though missed particles found
        fprintf('%s  %s\n','Warning no SET file found, data likely missing',topDir); %this can happen if only a few particles in folder
        nameSEMtmp = {findSEM.name};
        nameSEMtmp2 = fullfile(topDir,nameSEMtmp);
        nameSEM = [nameSEM; nameSEMtmp2'];
        if iscell(nameSEMtmp2)
            nameSET = [nameSET; cell(length(nameSEMtmp2),1)];
            namePKL = [namePKL; cell(length(nameSEMtmp2),1)];
        else
            nameSET = [nameSET; cell(1,1)];
            namePKL = [namePKL; cell(1,1)];
        end
    else %hit and missed particles found
        nameSETtmp = {findSET.name}; %get names of files
        nameSEMtmp = {findSEM.name};
        namePKLtmp = {findPKL.name};
        compSET = regexprep(nameSETtmp,STUDY.hitEXT,''); %remove file type to compare
        compSEM = regexprep(nameSEMtmp,STUDY.missedEXT,'');
        compPKL = regexprep(namePKLtmp,STUDY.spectraEXT,'');
        if ~isequal(compSET,compPKL) %names of hit part file and peak files must match
            error('%s\n%s', 'Number and names of pkl and set files must match ', topDir);
        elseif ~isequal(compSEM,compSET) %the names of hit and missed particles probably should match. If they don't, files are likely missing or in wrong folder
            fprintf('%s\n%s','Warning number or names of SEM and SET files do not match. Likely data missing', topDir);
            %procedure to sort and figure out what sem and set files
            %match or do not match
            allname = union(compSEM, compSET);
            [~,matchSEM] = intersect(allname,compSEM);
            [~,matchSET] = intersect(allname,compSET);
            %keep the length of the lists equal only add names where appropiate
            nameSETtmp2 = cell(length(allname),1);
            nameSEMtmp2 = cell(length(allname),1);
            namePKLtmp2 = cell(length(allname),1);
            nameSETtmp2(matchSET) = fullfile(topDir,nameSETtmp);
            nameSEMtmp2(matchSEM) = fullfile(topDir,nameSEMtmp);
            namePKLtmp2(matchSET) = fullfile(topDir,namePKLtmp);
            %save files names
            nameSET = [nameSET; nameSETtmp2];
            nameSEM = [nameSEM; nameSEMtmp2];
            namePKL = [namePKL; namePKLtmp2];
        else %all file names match
            nameSETtmp2 = fullfile(topDir,nameSETtmp); %create full file names
            nameSEMtmp2 = fullfile(topDir,nameSEMtmp);
            namePKLtmp2 = fullfile(topDir,namePKLtmp);
            nameSET = [nameSET; nameSETtmp2']; %save files names
            nameSEM = [nameSEM; nameSEMtmp2'];
            namePKL = [namePKL; namePKLtmp2'];
        end
    end
    
    %save inst files and date processed
    nameINSTtmp = fullfile(topDir, nameINSTtmp);
    ProcDATEtmp = now;
    if iscell(nameSEMtmp2)
        nameINSTtmp2 = cell(length(nameSEMtmp2),1);
        nameINSTtmp2([1:end]) = {nameINSTtmp};
        nameINST = [nameINST; nameINSTtmp2];
        ProcDatetmp2 = zeros(length(nameSEMtmp2),1);
        ProcDatetmp2(1:end,1) = ProcDATEtmp;
        procDATE = [procDATE; ProcDatetmp2];
    else
        nameINST = [nameINST; nameINSTtmp];
        procDATE = [procDATE; ProcDatetmp];
    end
end
% end
%     
%     %save date processed
%     ProcDATEtmp = now;
%     if iscell(nameSEM
%         ProcDatetmp2 = zeros(length(nameSEMtmp2),1);
%         ProcDatetmp2(1:end,1) = ProcDATEtmp;
%         procDATE = [procDATE; ProcDatetmp2];
%     end

%iterate into subfolders
if ~isempty(DirNames)
    %find all pkls
    for i = 1:length(DirNames)
        findData(fullfile(topDir, DirNames{i}));
    end
end 