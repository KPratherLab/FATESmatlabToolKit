%move pkl files into their respective folder
function redopkls(topDir,PoscalibIntercept,PoscalibSlope,NegcalibIntercept,NegcalibSlope)
% Call as redopkls(topDir,PoscalibIntercept,PoscalibSlope,NegcalibIntercept,NegcalibSlope)
% redopkls is a wrapper function that is used to write calibrated pkl text files (in m/z space) from uncalibrated ams file (in time space) using
% calibration parameters provided.  Also writes supplementary .pol files
% These .pkl and .pol files are the file format the yaada accepts for
% atofms data.  
% topDir - the topfolder over which redopkls will search all folders and
% subfolders for .ams files
% PoscalibIntercept - the intercept for the positive spectra calibration
% PoscalibSlope - the slope for the positive spectra calibration
% NegcalibIntercept - the intercept for the negative spectra calibration
% NegcalibSlope - the slope for the negative spectra calibration


if nargin < 5
   error('Too few input arguments.');
elseif nargin > 5  
    error('Too many input arguments.');
end

%find all subfolders
DirList = dir(topDir);
findDir = [DirList.isdir]; 
DirNames = {DirList(findDir).name};
DirNames(ismember(DirNames,{'.','..'})) = [];

%find ams files
findAMS = dir(sprintf('%s%s',topDir,'\*.ams'));
if ~isempty(findAMS)
    nameAMS = {findAMS.name}; %get pkl names
    fullfileAMS = fullfile(topDir,nameAMS);
    AMStoPKL_freshStart(fullfileAMS,topDir,PoscalibIntercept,PoscalibSlope,NegcalibIntercept,NegcalibSlope); %make pkls
end

%iterate over all subfolders
if ~isempty(DirNames)
    %find all pkls
    for i = 1:length(DirNames)
        redopkls(fullfile(topDir, DirNames{i}),PoscalibIntercept,PoscalibSlope,NegcalibIntercept,NegcalibSlope);
    end
end    