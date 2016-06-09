%% Spectrum Viewer Config File
defaultPath = 'E:\data\cw2_2016\';
defaultPosA = 1.893889E-06;
defaultPosB = -4.170349E+02;
defaultNegA = 1.391047E-06;
defaultNegB = -3.9547796E+02;
%% Get current folder
fu = mfilename('fullpath'); % find folder for spectrumViewer
fid = fopen([fu(1:end-20) 'defaultConfig.cal'],'w'); %load file for writing
outFormat = '%s \n'; % output format
fprintf(fid,outFormat,defaultPath); % write Path
outFormat = '%d \n'; % output format
fprintf(fid,outFormat,defaultPosA); % write cal value
fprintf(fid,outFormat,defaultPosB); % write cal value
fprintf(fid,outFormat,defaultNegA); % write cal value
fprintf(fid,outFormat,defaultNegB); % write cal value
fclose(fid);

fid = fopen([fu(1:end-20) 'defaultConfig.cal'],'r'); % open file for loading cal values
defaultPath = fgetl(fid); % load Path 
defaultPosA = fgetl(fid); % load cal value
defaultPosB = fgetl(fid); % load cal value
defaultNegA = fgetl(fid); % load cal value
defaultNegB = fgetl(fid); % load cal value
fclose(fid);