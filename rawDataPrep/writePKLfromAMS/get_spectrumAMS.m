function [HitParticleCounter,IonType,Speed,TimeStampData,LaserPower,Data] = get_spectrumAMS(MS_filename)
% Call as [HitParticleCounter,IonType,Speed,TimeStampData,LaserPower,Data] = get_spectrumAMS_JC_mod(MS_filename)
% MS_filename - a fullfile path to a .ams file which stores the raw spectra data
% HitParticleCounter - particle location in folder
% IonType - spectra identifier (pos or neg) 
% Speed - speed of particle
% TimeStampData - time spectra acquired
% LaserPower - laser power recorded
% Data - raw spectra data
% opens and stores the data from the .ams files.  Annoying to program as it
% was saved in binary format, likely for TasWare's convience. Pain in my
% ass - Jack Cahill
% originally written Jack Cahill 2014

filestream = fopen(MS_filename);

Version = fread(filestream,1,'short');
NumPoints = fread(filestream,1,'short');
HitParticleCounter = fread(filestream,1,'int');
ScatDelay = fread(filestream,1,'short');
Speed = fread(filestream,1,'float');
Size = fread(filestream,1,'float');
IonType = fread(filestream,1,'short');
TimeStamp_type = fread(filestream,1,'short');
TimeStampData = fread(filestream,1,'double');  % good

TimeText = fread(filestream,20,'*char');


LaserPower = fread(filestream,1,'float');
ParticleHit = fread(filestream,1,'short');
FileType = fread(filestream,1,'short');

FileName = (fread(filestream,80,'*char'));

% I assume all of these outputs are selectable from Tasware, like selecting
% labels, re-caling and such
Cruft = fread(filestream,107,'*char');
TotalIntegral = fread(filestream,1,'int');
Baseline = fread(filestream,1,'short');
Calibrated = fread(filestream,1,'short');
CalibSlope = fread(filestream,1,'double');
CalibIntercept = fread(filestream,1,'double');
CalibData = fread(filestream,32,'double');
Labels = fread(filestream,1,'short');
Label = fread(filestream,480,'char'); % this value is actually a structure, 20 bytes of char, 4 bytes (float). Am not going to seperate
Cruft2 = fread(filestream,158, '*char');

Data = fread(filestream,15000,'short');  % DATA!!
% 
% Cruft3 = fread(filestream,1,2);
fclose(filestream);