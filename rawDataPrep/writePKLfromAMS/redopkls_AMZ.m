function redopkls_AMZ(topDir,PoscalibIntercept,PoscalibSlope,NegcalibIntercept,NegcalibSlope)
% writes calibrated pkl text files (in m/z space) from uncalibrated compressed amz file (in time space) using
% calibration parameters provided.  AMZ files first have to be unzipped to
% AMS files.
% topDir - the topfolder over which redopkls will search all folders and
% subfolders for .amz files
% PoscalibIntercept - the intercept for the positive spectra calibration
% PoscalibSlope - the slope for the positive spectra calibration
% NegcalibIntercept - the intercept for the negative spectra calibration
% NegcalibSlope - the slope for the negative spectra calibration
% Camille Sultana 2016

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

findAMZ = dir(fullfile(topDir,'*.amz')); %find .amz files
if ~isempty(findAMZ)
    nameAMZ = {findAMZ.name}; %get pkl names
    nameAMZ(ismember(nameAMZ,{'avg_pos.amz', 'avg_neg.amz'})) = [];
    fullfileAMZ = fullfile(topDir,nameAMZ);
    if iscell(fullfileAMZ)
        fullfileAMS = cell(1,length(fullfileAMZ));
        for i = 1:length(fullfileAMZ)
            tmp = unzip(fullfileAMZ{i},topDir); %uncompress amz file
            fullfileAMS{i} = tmp{1};
        end
    else
        fullfileAMS = unzip(fullfileAMZ,topDir);
    end
    AMStoPKL_freshStart(fullfileAMS,topDir,PoscalibIntercept,PoscalibSlope,NegcalibIntercept,NegcalibSlope); %write pkl file
    delete(fullfile(topDir,'*.ams')); %delete uncompressed ams file
end

if ~isempty(DirNames)
    %find all pkls
    for i = 1:length(DirNames)
        redopkls_AMZ(fullfile(topDir, DirNames{i}),PoscalibIntercept,PoscalibSlope,NegcalibIntercept,NegcalibSlope);
    end
end    