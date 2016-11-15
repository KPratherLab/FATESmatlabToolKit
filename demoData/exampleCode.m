%%
% example code provided with FATES
% make sure FATES toolkit is in matlab path

%% raw data prep: do this before making a FATES study

%add pkls and pol files
%writing pkls from ams files is most time consuming step currently
% topDir = 'C:\FATESmatlabToolKit\demoData\raw\';
% PoscalibIntercept = 32.98938;
% PoscalibSlope = 8.925975E-6;
% NegcalibIntercept = 39.38443;
% NegcalibSlope = 8.794423E-6;
% redopkls(topDir,PoscalibIntercept,PoscalibSlope,NegcalibIntercept,NegcalibSlope);

% generate psl cal
% would NEED to run this to populate .inst file calibration
% parameters
% pslTime = [startTime1 endTime1; startTime2 endTime2; etc]
% pslVels = [minVel1 maxVel1; minVel2 maxVel2; etc]
% pslSizes = [pslSize1 pslSize2 etc];
% N = 3 (default)
% bins = 100 (default);
% [calParams, speedModes] = pslCal(pslTime, pslVels, pslSizes, N, bins)

%add inst file
% InstFile = 'C:\FATESmatlabToolKit\demoData\ELD2015_lp.inst';
% addInst(topDir,InstFile)

%% startup yaada

startup_yaada
%create startup file if this is the first time

%% make a study
% once you have completed the raw data prep you can now make a study using
% the data provided
init_study %initialize study
make_study %make a study
% open_study %open study that is already loaded in (load(STUDYNAME))

%% merge studies
% merge_study_load %initialize 2 study merge Study1 and Study2 should already have been created.  

%After merge_study_load run one of the three options below to create new
%merged study
% a) merely append all inst data b/c the PIDinst ids are all different
% merge_study_inst(1)   %1 signifies append data

% b) append all data but renumber inst ids in study2 b/c they overlap to
%    inst ids in study 1
%          merge_study_inst(0)  or merge_study_inst()
%          That program will renumber inst2 ids to 
%              start after the last inst1 ids           

% c)merge data from the same instumnets in study1 and study2
%          merge_study_part({1,2,4},{1,6,3})
%             This means merge Study1,inst id 1 with Study2,inst id 1;
%                              Study1,inst id 2 with Study2,inst id 6;
%                              Study1,inst id 4 with Study2,inst id 3;
%            The final merged data will have inst ids from Study1


%% load in missed particle data if desired
%these variables take up a lot of memory, which is why it is not held in the 
%study all the time
%to run this you need to have already created a FATES study and have it
%loaded into the workspace

[PARTidMissed, PARTdataMissed] = load_missed(STUDY.PARTidMissed_filename, STUDY.PARTdataMissed_filename);

%% do some analyses!

%get integer spectrum
partData = PARTidMat(1:500,1:2); %PID now consists of a 2 column matrix, with the full particle identifier made up of each row
MaxMZ = 300; %300;
ResponseType = 'Area';
Polarity = 2;
[negData,posData] = get_int_spectrum_SUM(partData,MaxMZ,ResponseType,Polarity);
%NOTE ON RETRIEVING SPECTRA: Reading the external PEAK file is what takes
%the most time in pretty much any analysis.  Thus you should try to limit
%the number of times you access this file (limit # of calls to
%get_int_spectrum_SUM, get_NONint_spectrum_SUM, get_spectrum, etc).  For
%example calling get_int_spectrum once on a longer PID list is often more 
%efficient than calling it multiple times on shorter PID list

% get non integer spectrum (for a more refined look at peak quality and resolution)
% MZbins = 0.25:0.5:300.5;
%[NegResponse,PosResponse] = get_NONint_spectrum_SUM(PID,MZbins,ResponseType,Polarity)

%run art2a
Polarity = 2;
MaxMass = 300; %350;
LearningRate = 0.05;
VigilanceFactor = 0.8;
MaxIteration = 20;
RegroupFactor = 0.85;

%art2a
[outPartID, outWeightMatrix, partIDCount, trackNeuron, avgPROXALL] = run_art2a(partData, Polarity, MaxMass, LearningRate, VigilanceFactor, MaxIteration, negData, posData); %run art2a initially

[outPartID2, outWeightMatrix2] = regroup_art2a(outPartID, outWeightMatrix, RegroupFactor); %combine output clusters
% [outPartID3] = match_art2a(PARTidMat(301:500,1:2), outWeightMatrix2, Polarity, VigilanceFactor); %match a list of particles to existing weight matrices

outPartIDX2 = zeros(size(partData,1),1);
for i = 1:length(outPartID2)
[~,tmp] = intersect(partData, outPartID2{i}, 'rows');
outPartIDX2(tmp) = i;
end


%% all further analyses shown rely upon using only logical indexing and 2 functions: intersect and histc!!

%logical indexing examples
%working on all particles < 1
submicronIDX = PARTdataMat(:,PARTdataFlds.DA) < 1; 
submicronPID = PARTidMat(submicronIDX,1:2);
submicronVEL = PARTdataMat(submicronIDX,PARTdataFlds.VELOCITY); %get velocities
submicronTIME = PARTidMat(PARTdataMat(:,3) < 1,PARTidFlds.TIME); %get time (combined all work into one line)

%get all particles < 1 and which generated negative and positive spectra
%Note using explicit column #s here rather than column fields. This is a
%little faster to type but less flexible if want to apply code on FATES
%studies with different structures
twoCondINDX = PARTdataMat(:,3) < 1 & PARTdataMat(:,4) == 1; %get index. 
twoCondTIME = PARTidMat(twoCondINDX,3); %get time

% Additionally lists of particles can be compared with the built-in function intersect.  
% For example another way to return the particle ids of submicron particles
% with both pos and neg spectra 
bothSpecIDX = PARTdataMat(:,4) == 1; 
bothSpecPID = PARTidMat(bothSpecIDX,1:2);
% generate particle id list of submicron particles with both pos and neg
% mass spectra
twoCondPID = intersect(bothSpecPID, submicronPID,'rows');

% find the PIDs for particles in partData with -35 > 250
selectPID = partData(negData(35,:)>250,1:2);
% find the particle indexes for the PIDs in the entire study
[~,selectIDX] = intersect(PARTidMat(:,1:2),selectPID,'rows');

% get particles with peak areas of -35 OR -93 greater than 250
peakPID = partData(negData(35,:)>250 | negData(93,:) > 250,1:2); 

%intersect example
%SPECIFYING ROWS WHEN INTERSECTING the 'PIDS' IS IMPORTANT
%get spectra for all particles in cluster 1
[~,idxPID] = intersect(partData,outPartID2{1},'rows'); %find where cluster 1 particles are in PID list used to generate spectra matrix
posC1 = posData(:,idxPID);

%histc example
%get particles binned by size
allSize = PARTdataMat(:,PARTdataFlds.DA); %sizes for all particles 
bins = 0:0.5:3; %size bins
[binCounts,~,binIDX] = histcounts(allSize,bins); %create bins and indexes
sizePIDS = cell(1,length(bins)); %list of pids by bin
for i = 1:length(bins)
    binTmp = binIDX == i;
    sizePIDS{i} = PARTidMat(binTmp,1:2);
end

% A similar technique can also be used to find the distribution of m/z signal for a set of particles.
% get the response at m/z -62 for all partData
neg35 = negData(35,:);
% create bins for m/z 62
neg35bin = min(neg35):100:max(neg35);
% Find number of particles in each response bin and the response bin into which each particle belongs
[neg35BinCounts,~, neg35BinIDX] = histcounts(neg35, neg35bin);

%% GUIfates Initialization Demo
% Set parameters to get peak areas
MZbins = 0.25:0.5:300.5;
ResponseType = 'RelArea';
Polarity = 2;
% Get peak areas using get_NONint_spectrum_SUM
[negResponse, posResponse] = get_NONint_spectrum_SUM(partData,MZbins,ResponseType,Polarity);
% concatenate peak areas for GUIFates initialization--IMPORTANT to do it this way
useMZdata = [fliplr(negResponse') zeros(size(negResponse,2),1) posResponse']; % 
MZvector = -300:0.5:300; %MZ need to be monotonically increasing or decreasing

% particle data inputs
timeData = PARTidMat(1:500,PARTidFlds.TIME); % Time data for each particle--CAN BE SUBSTITUTED FOR A DIFFERENT PARAMETER--one example is to use altitude data for aircraft studies

outPartIDX2;  % Cluster identifier for each particle output by a clustering algorithm, see art2a section
outPartIDX2(outPartIDX2 > 10 | outPartIDX2 == 0) = 10; %simplifying a little bit for example purposes to make display easy
sizeData = PARTdataMat(1:500,PARTdataFlds.DA); % Size data for each particle--CAN BE SUBSTITUTED FOR A DIFFERENT PARAMETER--one example is to use total ion yield
clustRelation = []; % cluster relation statistic

% Call GUIfates
GUIfates(useMZdata,MZvector,partData,timeData,sizeData,outPartIDX2,clustRelation)

