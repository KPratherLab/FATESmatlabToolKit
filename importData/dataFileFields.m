%define in what column the data is stored in in the missed (.sem) and hit
%(.set) particle files and the spectra (pkl) files.  These files are read in during parse_part by
%textscan and the data is used to fill a matrix with each column pertaining
%to a certain variable (speed, laser power, area, etc).  All fieldnames defined here
%that match a fieldname in the corresponding structure in 
%studyFields will be stored in the study variables. All other data will be ignored.
%The fieldnames for missedFlds and hitFlds MUST match
%the fieldnames defined by PARTdataFlds (data type 1 and data type 3) and PARTidFlds.TIME.
%The fieldnames for spectraFlds MUST match the fieldsnames defined by
%PEAKFlds.

global missedFlds hitFlds spectraFlds missedNAME hitNAME spectraNAME

%set up missed data fields (sem data format) (currently .sem and .set files have
%identical organization so missedFlds == hitFlds
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