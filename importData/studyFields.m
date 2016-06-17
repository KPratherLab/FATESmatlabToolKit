function studyFields
%this sets up all the structures that yaada will rely upon for the study.
%Users can alter this file depending upon the data types/variables in their
%own data sets.  However for each structure (INST, PEAK, etc) read the comments carefully as
%to what fields can be altered/removed and what must remain for yaada to
%continue functioning.

global PEAKFlds PARTidFlds PARTdataFlds PEAK STUDY INST partdataNAME numFldsPARTid numFldsPARTdata peakFldsNAME numFldsPEAKFlds

%set up INST flds
% NOTES ON FLEXIBILITY
%The INST flds to be used are flexible EXCEPT do not remove InstID,
%InstCODE, ExpNAME, DAcalibFUNCTION, or DAcalibPARAM, or LastPartID.  Other scripts/yaada
%depend upon these existing.  You may rename/change/add/delete all other
%INST flds at your discretion
%Note a unique InstID within the study for every inst file that has a unique
%InstCODE and ExpNAME.  Other inst fields do not determine the InstID created (see parse_inst)
INST.InstID = []; %DO NOT DELETE %The unique instid associated with each unique InstCODE and ExpNAME pair. Each particle has an instid as part of its identifying information
INST.InstCODE = []; %DO NOT DELETE %A code associated with the instrument used to acquire the data ex: ELD, LVN, TSI, etc
INST.InstDESC = []; %A description of the instrument used to acquire the data
INST.ExpNAME = []; %DO NOT DELETE %The name of the experiment in which data was collected. 
INST.ExpDESC = []; %A description of the experiment in which data was collected
INST.DAcalibFUNCTION = []; %DO NOT DELETE %a name of the function evaluated which calibrates the particle sizes from velocity
INST.DAcalibPARAM = []; %DO NOT DELETE %the parameters used in the calibration function
INST.LastPartID = []; %The number of particles associated with this InstID

%set up variables to save data definitions file
instFldsNAME = fieldnames(INST);
instTABLE = cell(length(instFldsNAME),1);
instTABLE(1:end) = {'INST'};
instDESC = cell(length(instFldsNAME),1);
%the instDESC variable has to be manually altered if you change the fields
%this variable is just a reference for future users so if you don't change
%it, the study will still load properly but the record of what each field
%is in DATADEF will not be accureate
instDESC(1:end) = {'inst parameters identifier'; 'instrument code'; 'instrument description'; 'experiment name'; 'experimnet description'; 'calibration function'; 'calibration parameters'; 'number of particles with this instid'};

%set up peak fields
% NOTES ON FLEXIBILITY
%The PEAKFlds flds to be used are flexible EXCEPT do not remove or InstID,
%PARTID and do not remove MZ.  Other scripts/yaada
%depend upon these existing. In addition if you remove SPECID (identifies polarity of spectra)
%you will have to alter get_spectrum, get_int_spectrum_SUM, and
%get_NONint_spectrum_SUM, to remove references to SPECID, but this should
%be relatively straightforward.  You would only consider doing this if you only had single
%polarity data. You may rename/change/add/delete all other
%PEAKFlds flds at your discretion.  
PEAKFlds.INSTID = 1; %DO NOT DELETE %instrument identifier
PEAKFlds.PARTID = 2; %DO NOT DELETE %particle identifier
PEAKFlds.SPECID = 3; %probably DO NOT DELETE %spectrum (polarity) identifier: 1 is pos, 0 is neg
PEAKFlds.PEAKID = 4; 
PEAKFlds.MZ = 5;  %DO NOT DELETE %mz of peak
PEAKFlds.AREA = 6; %area of peak
PEAKFlds.HEIGHT = 7; %peak height
PEAKFlds.BLOWSCALE = 8; %true if peak blows scale
PEAKFlds.RELAREA = 9; %peak area relative to total area of spectrum
numFldsPEAKFlds = length(fieldnames(PEAKFlds)); %DON'T DELETE
STUDY.NumPkCols= numFldsPEAKFlds;

%set up variables to save data definitions file
peakFldsNAME = fieldnames(PEAKFlds);
peakFTABLE = cell(length(peakFldsNAME),1);
peakFTABLE(1:end) = {'PEAKFlds'};
peakFDESC = cell(length(peakFldsNAME),1);
%the peakFDESC variable has to be manually altered if you change the fields
%this variable is just a reference for future users so if you don't change
%it, the study will still load properly but the record of what each field
%is in DATADEF will not be accureate
peakFDESC(1:end) = {'inst parameters identifier'; 'particle identifier'; 'spectrum (polarity) identifier'; 'peak identifier'; 'mz of peak'; 'area of peak'; 'peak height'; 'true if peak blows scale'; 'peak area relative to total area of spectrum'};

%set up PEAK
% NOTES ON FLEXIBILITY
%The PEAKFlds flds to be used are flexible EXCEPT do not change/remove rowcnt,
%PARTID, INSTID.  Other scripts/yaada
%depend upon these existing.  You may rename/change/add/delete all other
%PEAKFlds flds at your discretion
PEAK.rowcnt = []; %DO NOT CHANGE/DELETE %nth peak in whole study
PEAK.INSTID = []; %DO NOT CHANGE/DELETE %instrument identifier
PEAK.PARTID = []; %DO NOT CHANGE/DELETE %particle identifier for hit particle

%set up variables to save data definitions file
peakNAME = fieldnames(PEAK);
peakTABLE = cell(length(peakNAME),1);
peakTABLE(1:end) = {'PEAK'};
peakDESC = cell(length(peakNAME),1);
%the peakDESC variable has to be manually altered if you change the fields
%this variable is just a reference for future users so if you don't change
%it, the study will still load properly but the record of what each field
%is in DATADEF will not be accureate
peakDESC(1:end) = {'nth peak in whole study'; 'instrument identifier'; 'particle identifier'};

%set up partID fields
% NOTES ON FLEXIBILITY
%The PARTidFlds flds to be used are flexible EXCEPT do not change TIME,
%PARTID, INSTID.  Other scripts/yaada
%depend upon these existing.  You may rename/change/add/delete all other
%PARTidFlds flds at your discretion
%NOTE PARTICLE DATA has been split between PARTid (double precision) and
%PARTdata (single precision) matrices to save memory.  If you want to add
%fields add to the appropriate matrix/structure keeping memory (using less
%memory is good!) and precision (how important is it to retain the value
%EXACTLY) in mind
PARTidFlds.INSTID = 1; %DO NOT CHANGE/DELETE %instrument identifier
PARTidFlds.PARTID = 2; %DO NOT CHANGE/DELETE %particle identifier
PARTidFlds.TIME = 3; %DO NOT CHANGE/DELETE %time of paritcle detection, read in from particle data file (.sem or .set)

%set up variables to save data definitions file
partidNAME = fieldnames(PARTidFlds);
partidTABLE = cell(length(partidNAME),1);
partidTABLE(1:end) = {'PARTid'};
partidDESC = cell(length(partidNAME),1);
%the partidDESC variable has to be manually altered if you change the fields
%this variable is just a reference for future users so if you don't change
%it, the study will still load properly but the record of what each field
%is in DATADEF will not be accureate
partidDESC(1:end) = {'inst parameters identifier'; 'particle identifier'; 'time of paritcle detection'};

numFldsPARTid = length(partidNAME); %DON'T DELETE

%set up partData fields
%NOTE PARTICLE DATA has been split between PARTid (double precision) and
%PARTdata (single precision) matrices to save memory.  If you want to add
%fields add to the appropriate matric/structure keeping memory (using less
%memory is good!) and precision (how important is it to retain the value
%EXACTLY) in mind.  The function parse_part is what reads in particle data
%and populates the PART matrices.  Make sure changes to the PART flds are
%compatible with parse_part.
% Data has been organized into 3 types in the PARTdata matrix, to simplify
% and make flexible the code.  Pay close attention to the comments if you
% plan on changing any of this.  

% Data type 1: Non-gauge board (sizing PMT) data read in directly from the
% particle file (.sem or .set).  Currently this only consists of Velocity and LaserPower.
% This data type NEEDS to remain the first set defined to keep code in
% parse_part consistent.  However, all PARTdataFlds.DataType1 flds you may 
% rename/change/add/delete at your discretion
PARTdataFlds.VELOCITY = 1; %velocity of particle 
PARTdataFlds.LASERPOWER = 2; %laser power of ldi
% PARTdataFlds.SAHEIGHT = 7; %PMTA Height
% PARTdataFlds.SAAREA = 8; %PMTA Area
% PARTdataFlds.SBHEIGHT = 9; %PMTB HEIGHT
% PARTdataFlds.SBAREA = 10; %PMTB AREA

%Data type 2: Data calculated in parse_part, except .RING which is read in
%from the .pol file if there is one.  HIT is also read in from the .pol
%file for the hit particles if it exists, other is given a value of 1.
%This data type NEEDS to remain the second set defined to keep code in
%parse_part consistent. You may add/delete these flds, but to add
%flds you will manually have to add in the relevant code (parse_part lines 80-127), as this is data
%assumed not to be simply read in from the particle file.
%Note if you remove the HIT fld all data (for hit and missed data)
%will be stored in the PART tables held in memory by MATLAB, significantly 
%increasing memory demands.  In addition there will not be a marker for 
%hit and missed particles.  So it is not advised to remove the HIT field.
PARTdataFlds.DA = 3; %vacuum aerodynamic diameter, calculated using provided calibration function and parameters and particle velocity
PARTdataFlds.HIT = 4; %>0 if hit. 0= not hit, if .pol file exists 1 = both pos and neg spectra generated, 2 = pos spectra only, 3 = neg spectra only. If no pol file 1 = any spectra generated ('hit particle')
PARTdataFlds.RING = 5; %>0 if ring detected in spectra. This is read in from .pol file if it exists.
PARTdataFlds.POSIT = 6; %position of particle in folder, taken from file list

%set up variables to save data definitions file
partdataNAME = fieldnames(PARTdataFlds); %DON'T DELETE
numFldsPARTdata = length(partdataNAME); %DON'T DELETE
partdataTABLE = cell(length(partdataNAME),1);
partdataTABLE(1:end) = {'PARTFlds'};
partdataDESC = cell(length(partdataNAME),1);
%the allDESC variable has to be manually altered if you change the fields
%this variable is just a reference for future users so if you don't change
%it, the study will still load properly but the record of what each field
%is in DATADEF will not be accureate
allDESC = {'velocity of particle'; 'laser power of ldi';'%vacuum aerodynamic diameter';'>0 if hit'; '>0 if ring detected in spectra';'%position of particle in folder'};
%  'PMTA Height'; 'PMTA Area'; 'PMTB HEIGHT'; 'PMTB Area'

%create datadef variable to hold study definitions for future reference
DATADEF.table = 1;
DATADEF.name = 2;
DATADEF.desc = 3;
%fill in table values
[DATADEF((1):(length(partidNAME))).table] = partidTABLE{:};
[DATADEF((length(partidNAME)+1):(length(partidNAME)+length(peakFldsNAME))).table] = peakFTABLE{:};
[DATADEF((length(partidNAME)+length(peakFldsNAME)+1):(length(partidNAME)+length(peakFldsNAME)+length(peakNAME))).table] = peakTABLE{:};
[DATADEF((length(partidNAME)+length(peakFldsNAME)+length(peakNAME)+1):(length(partdataNAME)+length(partidNAME)+length(peakFldsNAME)+length(peakNAME))).table] = partdataTABLE{:};
[DATADEF((length(partdataNAME)+length(partidNAME)+length(peakFldsNAME)+length(peakNAME)+1):(length(partdataNAME)+length(partidNAME)+length(peakFldsNAME)+length(instFldsNAME)+length(peakNAME))).table] = instTABLE{:};


%fill in table names
[DATADEF((1):(length(partidNAME))).name] = partidNAME{:};
[DATADEF((length(partidNAME)+1):(length(partidNAME)+length(peakFldsNAME))).name] = peakFldsNAME{:};
[DATADEF((length(partidNAME)+length(peakFldsNAME)+1):(length(partidNAME)+length(peakFldsNAME)+length(peakNAME))).name] = peakNAME{:};
[DATADEF((length(partidNAME)+length(peakFldsNAME)+length(peakNAME)+1):(length(partdataNAME)+length(partidNAME)+length(peakFldsNAME)+length(peakNAME))).name] = partdataNAME{:};
[DATADEF((length(partdataNAME)+length(partidNAME)+length(peakFldsNAME)+length(peakNAME)+1):(length(partdataNAME)+length(partidNAME)+length(peakFldsNAME)+length(instFldsNAME)+length(peakNAME))).table] = instFldsNAME{:};


%fill in description
[DATADEF((1):(length(partidNAME))).desc] = partidDESC{:};
[DATADEF((length(partidNAME)+1):(length(partidNAME)+length(peakFldsNAME))).desc] = peakFDESC{:};
[DATADEF((length(partidNAME)+length(peakFldsNAME)+1):(length(partidNAME)+length(peakFldsNAME)+length(peakNAME))).desc] = peakDESC{:};
[DATADEF((length(partidNAME)+length(peakFldsNAME)+length(peakNAME)+1):(length(partdataNAME)+length(partidNAME)+length(peakFldsNAME)+length(peakNAME))).desc] = partdataDESC{:};
[DATADEF((length(partdataNAME)+length(partidNAME)+length(peakFldsNAME)+length(peakNAME)+1):(length(partdataNAME)+length(partidNAME)+length(peakFldsNAME)+length(instFldsNAME)+length(peakNAME))).desc] = instDESC{:};

save(fullfile(STUDY.ProcDir,'datadef.mat'),'DATADEF'); %save data definitions file