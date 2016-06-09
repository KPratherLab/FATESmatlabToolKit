function usefig = setup_plotclick(PID_list,PID_clusters,XVector,YVector,useSetup,specData,mz)
% setup_plotclick sets up a gscatter plot of particles for the user to interact with
% using plotclick.  
% consider running clearPlotclick to clear global variables created during
% setup_plotclick after analysis with plotclick is complete. This will
% free up memory.
% 
% Call as usefig = setup_plotclick(PID_list,PID_clusters,XVector,YVector,useSetup,specData,mz)
% INPUT
% PID_list = list of PIDS, Nx2 matrix
% PID_clusters = cell array containing lists of PIDS (eg list of clusters).  All PIDS in
% PID_clusters must be in PID_list.  However PID_list can contain PIDS not
% in PID_clusters
% XVector = metric pertaining to particles in PID_list (eg size of particle), a NX1 vector. These 
% values will be plotted on the x axis of usefig.  The nth row of XVector
% should correspond to the nth particle (row) of PID_list.  
% YVector = metric pertaining to particles in PID_list, a NX1 vector (eg laser power). These 
% values will be plotted on the Y axis of usefig.  The nth row of XVector
% should correspond to the nth particle (row) of PID_list and each row is a m/z value.
% useSetup = 0 or 1. 0 if the user wants to give input for setting up plots
% by interacting in the command window. 1 if the user has pre-populated all
% the inputs into the script userSetupScatter
% OPTIONAL INPUT
% specData = MxN matrix containing particle spectra.  The Nth column of
% specData should be the spectra for the nth particle in PID_list.  Note that if
% using only neg spectra the first row should be the greatest absolute m/z
% and the last row should be the smallest absolute m/z (ie -150 (first
% row), -149, -148....-3, -2 -1 (last row)).  If using spectra that has both NEG and POS values
% the matrix should be created so the the first row is the greatest
% absolute neg m/z and the last row is the greatest positive m/z.  In
% addition a row for 0 m/z must exist for plotting purposes. (ie -150 (first
% row), -149...-2, -1, 0, 1, 2...149, 150 (last row))
% mz = A Gx1 vector of the m/z values corresponding to the rows of the spectra
% matrices contained in specData.  See the explanation above. 
% OUTPUT
% usefig = a figure handle to be used in plotclick for user exploration.
% Usefig is a gscatter plot of the particles in PID_clusters with XVector 
% values plotted on the x axis and YVector values plotted on the yaxis.  
% Particles are colored by the cluster they are in in PID_clusters.  
%
% Note: the order 1...N of Particle IDs in PID_list should match the 
%    order 1...N of particles in the data inputs (PosData, NegData, XVector, YVecotr)
%
% There are a wide number of plots that the user is prompted to select
% plotCS: line plot of average mz spectrum for selected particles
% plotPS: line plot of all mz spectra overlaid.  NOTE THIS PLOT IS VERY
% TIME/MEMORY CONSUMING
% plotHM: linear heatmap of mz spectra
% plotHMlog: log10 heatmap of mz spectra
% plotSize: relative number fraction over size of selected particles 
% plotSizeTime: scatter plot of particles over time and size
% plotSim: scatter plot of dot product1 (avg cluster spectra . particle
% spectra) and dot product2 (avg all selected spectra . particle spectra)

global PARTidFlds PARTdataFlds PARTdataMat scatterLogSPEC plotTrack rownum scatterSPEC scatterMZ XVectorUse YVectorUse PARTidMat labels2use trackintersectPART SCATclr origPID scatterTime scatterTimeBin  scatterSIZE scatterSIZEbin scatterAllSize ;

%% check inputs
if nargin < 5 || nargin > 7
  error('Call as usefig = setup_plotclick(PID_list,PID_clusters,XVector,YVector,useSetup,specData,mz)')
end

if nargin ==4
    disp('No spectra data are supplied. Will not be able to create spectra plots using plotclick')
end

if nargin == 5
    error('Must supply both specData (spectra data) and mz (m/z vector to plot against)')
end

if length(XVector) ~= length(YVector) 
    error('XVector and YVector must be same length')
end

if length(XVector) ~= size(PID_list,1)
    error('The number of rows of PID_list must match the length of XVector and YVector')
end

if exist('specData')
    if size(PID_list,1) ~= size(specData,2)
        error('The number of rows of PID_list must match the number of columns of specData');
    end
    if size(specData,1) ~= length(mz)
        error('The number of rows specData must match the length of vector mz');
    end
end

%% get input from user for setting up plots
if useSetup == 1 %read in file from user with everything all setup
    userSetupScatter
else %get user inputs one by one   
    % checks spectra matrix
    if exist('useSPEC')
        plotSpec = input('Do you want to make spectra plots? 0-no, 1-yes  ');
    else
        plotSpec = 0;
    end
    
    if plotSpec
        plotPS = input('Plot all particle spectra overlaid? 0-no, 1-yes  ');
        plotCS = input('Plot average spectra of each cluster overlaid? 0-no, 1-yes  ');
        plotHM = input('Create the heat map of all particle spectra? 0-no, 1-yes  ');
        plotHMlog = input('Create the heat map of all particle spectra with log10 scaling? 0-no, 1-yes  ');
        plotSim = input('Plot the similarity of spectra to averages? 0-no, 1-yes  ');
    else
        plotPS = 0;
        plotHM = 0;
        plotHMlog = 0;
        plotCS = 0;
        plotSim = 0;
    end
    
    plotSizeTime = input('Plot the size and time of all particles? 0-no, 1-yes  ');
    plotSize = input('Plot the size by number fraction of all particles? 0-no, 1-yes  ');
    
    %get info for binning by size
    if plotSize == 1 || plotSizeTime == 1
        minSize = input('\nWhat is the min size (um) to plot?  ');
        maxSize = input('What is the max size (um) to plot?  ');
        numBins = input('What # of size bins is desired?  ');
    end
    
    %get values for time scatter plot
    if plotSizeTime == 1
        minTime = input('What is the starting time to plot? (format ex: 27-Oct-2015 15:45:40)  ','s');
        maxTime = input('What is the ending time to plot? (format ex: 27-Oct-2015 15:45:40)  ','s');
        scatterTimeBin = [datenum(minTime) datenum(maxTime)];
    end
    
    %ask user if want to show titles and legends for plots
    plotRender = input('\nDo you want to do fast rendering (no titles, labels, legends are put onto plots)?\n This may give some speedups of graphics 0-no 1-yes  ');
    
    %ask user if want to show plots separately or make single plots
    plotType = input('\nWould you like all plots separate (0) or combined into a single window (1)?  ');
    plotTrack = [plotPS plotCS plotHM plotHMlog plotSizeTime plotSize plotSim plotRender plotType];
end


%% Find location of all particles in PID_clusters in PID_list and save indexes
numPart = cellfun('size',PID_clusters,1);
endPart = cumsum(numPart);
startPart = [1 endPart(1:(end-1))+1];
totalPart = endPart(end);

labels2use = zeros(totalPart,1); %will be group input into gscatter, also used by plotclick  
trackintersectMZ = zeros(totalPart,1); %index of location of appropriate spectra 
trackintersectPART = zeros(totalPart,1); %index of location of particles in PARTidMat, PARTdataMat in study in YAADA. , Used by plotclick
origPID = zeros(totalPart,2); %PID_clusters reorderedm. Used by plotclick
tmpNum = 1;
for i=1:length(PID_clusters)
    if ~isempty(PID_clusters{i});
        [~, iMZ]=intersect(PID_list,PID_clusters{i},'rows');  %get indx to particle spectra
        [~, iP] = intersect(PARTidMat(:,1:2),PID_clusters{i},'rows'); %get indx to all particle data in study
        if isempty(iMZ) || length(iMZ) < size(PID_clusters{i},1)
            error('particles in PID_cluster %i not found in PID_list',i)
        else
            labels2use(startPart(i):endPart(i)) = repmat(i,1,length(iMZ)); %append data
            origPID(startPart(i):endPart(i),:) = PID_clusters{i};
            trackintersectMZ(startPart(i):endPart(i)) = iMZ;
        end
        if isempty(iP) || length(iP) < size(PID_clusters{i},1)
            error('particles in PID_cluster %i not found in PARTidMat',i)
        elseif plotSize == 1 || plotSizeTime == 1
            trackintersectPART(startPart(i):endPart(i)) = iP; %append data
        end
    end
end;

%get spectra data for particles in PID_clusters only
if exist('specData')
    scatterSPEC = specData(:,trackintersectMZ);
    scatterMZ = mz;
    if size(scatterMZ,2) > 1
        scatterMZ = scatterMZ';
    end
end

if plotHMlog
    scatterLogSPEC = log10(scatterSPEC);
end

%% find size and time values for clusters
if plotSize == 1 || plotSizeTime == 1
    %bin all particles used by size
    binStep = (maxSize-minSize)/numBins; %get width of bin
    scatterSIZEbin = minSize:binStep:maxSize; %create size bin vector
    scatterAllSize = PARTdataMat(trackintersectPART,PARTdataFlds.DA); %get size of particles
    
    %bin all particles used by cluster by size
    scatterSIZE = zeros(length(scatterSIZEbin),length(PID_clusters));
    for i = 1:length(PID_clusters) %bin each cluster by size
        indexClust = labels2use == i;
        clustSize = scatterAllSize(indexClust);
        scatterSIZE(:,i) = histc(clustSize,scatterSIZEbin); %bin
    end
end

if plotSizeTime
    scatterTime = PARTidMat(trackintersectPART,PARTidFlds.TIME); %get size of particles
end

%% if plotting all plots in a single figure set up subplots
%set number of rows and columns in subplot
if plotType == 1
    numSPECplot = plotPS+plotHM+plotHMlog+plotCS+plotSizeTime;
    numPlot = numSPECplot+plotSizeTime+plotSize; %number of plots to combine
    if plotSize || plotSim
        colnum = 1;
        rownum = numSPECplot+1;
    else
        colnum = 1;
        rownum = numSPECplot;
    end
end

%% set up figure to be used by plotclick
usefig = figure;
XVectorUse = XVector(trackintersectMZ);
YVectorUse = YVector(trackintersectMZ);
gs      = gscatter(XVectorUse, YVectorUse,labels2use,[],[],[],'on');  % sym color?
SCATclr = {gs(:).Color};
