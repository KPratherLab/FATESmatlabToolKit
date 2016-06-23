function addInst(rawDir,InstFile)
% addInst adds an .inst file of your choosing (InstFile should include the
% full path to the .inst file) to the lowest subdirectories
% in the path rawDir (ie to all folders which contain no sub-folders).  

%Camille Sultana 2015

% get directory info
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
