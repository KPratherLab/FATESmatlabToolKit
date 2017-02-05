function dataFileFields
% define in what column the data is stored in in the missedData and hitData
% matrices created by read_part function reading in the and missed (.sem) and hit
%(.set) particle files.  also define in what column data is stored in 
% the and PeakDataTMP matrix created by read_peak reading in the spectra (pkl) files.   
% Note currently that missedData and hitData have two less columns than the data in .sem
% and .set files because read_part ignores the first two columns in these files. 
% For example
% hitFlds.LASERPOWER =2 %this means the second column of missedData has laser power data for the particle
% All fieldnames defined here that match a fieldname in the corresponding structure in 
% studyFields will be stored in the database. All other data will be ignored.
% The fieldnames for missedFlds and hitFlds should match
% the fieldnames defined by PARTdataFlds and PARTidFlds.TIME.
% The fieldnames for spectraFlds should match the fieldnames defined by PEAKFlds.

global missedFlds hitFlds spectraFlds missedNAME hitNAME spectraNAME

%% set up missed data fields 
% currently .sem and .set files have identical organization so missedFlds == hitFlds

% Current Prather ATOFMS .sem format
missedFlds.TIME = 3; %TIME ALWAYS NEEDS TO BE THE FIRST FIELD DEFINED DO NOT CHANGE 
missedFlds.VELOCITY = 1; 
missedFlds.LASERPOWER = 2;
% missedFlds.SAHEIGHT = 4;
% missedFlds.SAAREA = 5;
% missedFlds.SAHEIGHT = 6;
% missedFlds.SAAREA = 7;

%Alabama detectedparticles.txt format
%comment out or change if not using Alabama DetectedParticles.txt format
% missedFlds.TIME = 1; %TIME ALWAYS NEEDS TO BE THE FIRST FIELD DEFINED DO NOT CHANGE 
% missedFlds.VELOCITY = 2; 

%Commercial ATOFMS .sem format
%comment out or change if not using Commerical ATOFMS .sem format
% missedFlds.TIME = 2; %TIME ALWAYS NEEDS TO BE THE FIRST FIELD DEFINED DO NOT CHANGE 
% missedFlds.VELOCITY = 1; 

%always do after setting up missedFlds
missedNAME = fieldnames(missedFlds);
missedNAME = missedNAME(2:end);

%% set up hit data fields 
%(currently .sem and .set files have identical organization so missedFlds == hitFlds

% Current Prather ATOFMS .set format
%comment out or change if not using Prather ATOFMS .sem format
hitFlds.TIME = 3; %TIME ALWAYS NEEDS TO BE THE FIRST FIELD DEFINED DO NOT CHANGE
hitFlds.VELOCITY = 1;
hitFlds.LASERPOWER = 2;
% hitFlds.SAHEIGHT = 4;
% hitFlds.SAAREA = 5;
% hitFlds.SAHEIGHT = 6;
% hitFlds.SAAREA = 7;

%Alabama hitparticles.txt format
%comment out or change if not using Alabama HitParticles.txt format
% hitFlds.TIME = 1; %TIME ALWAYS NEEDS TO BE THE FIRST FIELD DEFINED DO NOT CHANGE
% hitFlds.VELOCITY = 2;

%Commercial ATOFMS .set format
%comment out or change if not using Commerical ATOFMS .sem format
% hitFlds.TIME = 3; %TIME ALWAYS NEEDS TO BE THE FIRST FIELD DEFINED DO NOT CHANGE 
% hitFlds.VELOCITY = 1;
% hitFlds.LASERPOWER = 2;

%always do after setting up hitFlds
hitNAME = fieldnames(hitFlds);
hitNAME = hitNAME(2:end);

%% set up spectra data fields 

%Current Prather ATOFMS .pkl format
%comment out or change if not using Prather ATOFMS .pkl format
spectraFlds.PARTID = 1; %PARTID ALWAYS NEEDS TO BE THE FIRST FIELD DEFINED DO NOT CHANGE %particle identifier %the particle identifiers in the pkl file are overwritten by particle identifiers that will be unique to the study. 
spectraFlds.SPECID = 2; %spectrum (polarity) identifier: 1 is pos, 0 is neg
spectraFlds.PEAKID = 3; 
spectraFlds.MZ = 4;  %mz of peak
spectraFlds.AREA = 5; %area of peak
spectraFlds.HEIGHT = 6; %peak height
spectraFlds.BLOWSCALE = 7; %true if peak blows scale
spectraFlds.RELAREA = 8; %peak area relative to total area of spectrum

%Alabama hitparticles.txt format
%comment out or change if not using Alabama HitParticles.txt format
% spectraFlds.PARTID = 4; %PARTID ALWAYS NEEDS TO BE THE FIRST FIELD DEFINED DO NOT CHANGE %particle identifier %the particle identifiers in the pkl file are overwritten by particle identifiers that will be unique to the study. 
% spectraFlds.SPECID = 3; %spectrum (polarity) identifier: 1 is pos, 0 is neg
% spectraFlds.MZ = 2;  %mz of peak
% spectraFlds.AREA = 1; %area of peak

%Current Commerical ATOFMS .pkl format
%comment out or change if not using Commercial ATOFMS .pkl format
% spectraFlds.PARTID = 6; %PARTID ALWAYS NEEDS TO BE THE FIRST FIELD DEFINED DO NOT CHANGE %particle identifier %the particle identifiers in the pkl file are overwritten by particle identifiers that will be unique to the study. 
% spectraFlds.SPECID = 7; %spectrum (polarity) identifier: 1 is pos, 0 is neg
% spectraFlds.MZ = 1;  %mz of peak
% spectraFlds.AREA = 4; %area of peak
% spectraFlds.HEIGHT = 2; %peak height
% spectraFlds.BLOWSCALE = 5; %true if peak blows scale
% spectraFlds.RELAREA = 3; %peak area relative to total area of spectrum

%always do after setting up spectraFlds
spectraNAME = fieldnames(spectraFlds);
spectraNAME = spectraNAME(2:end);