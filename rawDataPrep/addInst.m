function addInst(rawDir,InstFile)
% addInst adds an .inst file of your choosing (InstFile should include the
% full path to the .inst file) to the lowest subdirectories
% in the path rawDir (ie to all folders which contain no sub-folders).  

%Camille Sultana 2015

%InstFile = 'C:\Users\atofms\Desktop\JKE_CAICE_20130121.inst'; %path for .inst file to be added
%InstFile = 'C:\Users\atofms\Desktop\LVN030413.inst';
%InstFile = 'C:\Users\atofms\Desktop\SLY_CAICE2014.inst';
% InstFile = 'E:\portable_tank_data\test_newYaada\testweirdfolders\d\LVN030413.inst';
% InstFile = 'G:\camilleControlSalt\ELD2015_lp.inst';
dirData = dir(rawDir);
dirIndex = [dirData.isdir];
folderList = {dirData(dirIndex).name}';
folderList(ismember(folderList,{'.','..'})) = [];

if ~isempty(folderList) %identify if folder contains subfolders
    for i = 1:length(folderList) 
            addInst(fullfile(rawDir,folderList{i}),InstFile); %iterate     
    end
else
    findInst = dir(fullfile(rawDir,'*.inst')); %find any inst files there currently. Inst files found will be deleted.
    if ~isempty(findInst)
        nameInst = {findInst.name}; %get pkl names
        fullfileInst = fullfile(rawDir,nameInst);
        for i = 1:length(fullfileInst)
            delete(fullfileInst{i});
        end
    end
    copyfile(InstFile, rawDir); %if there are no subfolders copy the .inst file           
end

end 
