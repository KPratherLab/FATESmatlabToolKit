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

%set up missed data fields (sem data format) (currently .sem and .set files have
%identical organization so missedFlds == hitFlds
%
missedFlds.TIME = 3; %TIME ALWAYS NEEDS TO BE THE FIRST FIELD DEFINED DO NOT CHANGE 
missedFlds.VELOCITY = 1; 
missedFlds.LASERPOWER = 2;
missedFlds.SAHEIGHT = 4;
missedFlds.SAAREA = 5;
missedFlds.SAHEIGHT = 6;
missedFlds.SAAREA = 7;
missedNAME = fieldnames(missedFlds);
missedNAME = missedNAME(2:end);

%set up hit data fields (set data format) (currently .sem and .set files have
%identical organization so missedFlds == hitFlds
hitFlds.TIME = 3; %TIME ALWAYS NEEDS TO BE THE FIRST FIELD DEFINED DO NOT CHANGE
hitFlds.VELOCITY = 1;
hitFlds.LASERPOWER = 2;
hitFlds.SAHEIGHT = 4;
hitFlds.SAAREA = 5;
hitFlds.SAHEIGHT = 6;
hitFlds.SAAREA = 7;
hitNAME = fieldnames(hitFlds);
hitNAME = hitNAME(2:end);

%set up spectra data fields (pkl file format)
spectraFlds.PARTID = 1; %PARTID ALWAYS NEEDS TO BE THE FIRST FIELD DEFINED DO NOT CHANGE %particle identifier %the particle identifiers in the pkl file are overwritten by particle identifiers that will be unique to the study. 
spectraFlds.SPECID = 2; %spectrum (polarity) identifier: 1 is pos, 0 is neg
spectraFlds.PEAKID = 3; 
spectraFlds.MZ = 4;  %mz of peak
spectraFlds.AREA = 5; %area of peak
spectraFlds.HEIGHT = 6; %peak height
spectraFlds.BLOWSCALE = 7; %true if peak blows scale
spectraFlds.RELAREA = 8; %peak area relative to total area of spectrum
spectraNAME = fieldnames(spectraFlds);
spectraNAME = spectraNAME(2:end);