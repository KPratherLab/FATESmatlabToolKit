%define in what column the data is stored in in the missed (.sem) and hit
%(.set) particle files.  These files are read in during parse_part by
%textscan and the data is used to fill a matrix with each column pertaining
%to a certain variable. The fieldnames for missedFlds and hitFlds MUST match
%the fieldnames defined by PARTdataFlds (data type 1 and data type 3) and PARTidFlds.TIME in studyFields.
%The value assigned to each field is the column where the data is stored in
%the particle file

global missedFlds hitFlds missedNAME hitNAME

%set up missed data fields (currently .sem and .set files have
%identical organization so missedFlds == hitFlds
missedFlds.TIME = 3; %TIME ALWAYS NEEDS TO BE THE FIRST FIELD DEFINED 
missedFlds.VELOCITY = 1;
missedFlds.LASERPOWER = 2;
missedFlds.SAHEIGHT = 4;
missedFlds.SAAREA = 5;
missedFlds.SAHEIGHT = 6;
missedFlds.SAAREA = 7;
missedNAME = fieldnames(missedFlds);
missedNAME = missedNAME(2:end);

%set up hit data fields (currently .sem and .set files have
%identical organization so missedFlds == hitFlds
hitFlds.TIME = 3; %TIME ALWAYS NEEDS TO BE THE FIRST FIELD DEFINED 
hitFlds.VELOCITY = 1;
hitFlds.LASERPOWER = 2;
hitFlds.SAHEIGHT = 4;
hitFlds.SAAREA = 5;
hitFlds.SAHEIGHT = 6;
hitFlds.SAAREA = 7;
hitNAME = fieldnames(hitFlds);
hitNAME = hitNAME(2:end);