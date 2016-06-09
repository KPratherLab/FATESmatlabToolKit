function dendoFig = setup_clickDendo(outWM,useSetup,clustCOUNTS,useSPEC,useMZ,clusterTimeBin,TimeBin,clusterSizeBin,SizeBin)
% setup_clickDendo sets up a dendogram created by grouping spectra using
% hierachichal clustering.  The user can then interact with the dendogram 
% using clickDendo to make helpful plots
% consider running clearDendo to clear global variables created during
% setup_clickDendo after analysis with clickDendo is complete. This will
% free up memory.
% Call as dendoFig = setup_clickDendo(outWM,useSetup,clustCOUNTS,useSPEC,useMZ,clusterTimeBin,TimeBin,clusterSizeBin,SizeBin)
% 
% INPUT
% outWM= a MxN matrix containing N spectra. The spectra are used to group 
% by hierachical clustering (this would usually be the output outWM from art2a)
% and then create the dendogram.  
% useSetup = 0 or 1. 0 if the user wants to give input for setting up plots
% by interacting in the command window. 1 if the user has pre-populated all
% the inputs into the script userSetupDendo
% clustCOUNTS = a 1xN vector containing the number of particles in each
% cluster in the cluster order corresponding to outWM.
%
% OPTIONAL INPUTS
% useSPEC = 1xN cell array.  Each cell should contain an GxH matrix which
% contains the individual spectra for particles belonging to cluster N. Each
% column is a particle spectra and each row is a m/z value.  Note that if
% using only neg spectra the first row should be the greatest absolute m/z
% and the last row should be the smallest absolute m/z (ie -150 (first
% row), -149, -148....-3, -2 -1 (last row)).  If using spectra that has both NEG and POS values
% the matrix should be created so the the first row is the greatest
% absolute neg m/z and the last row is the greatest positive m/z.  In
% addition a row for 0 m/z must exist for plotting purposes. (ie -150 (first
% row), -149...-2, -1, 0, 1, 2...149, 150 (last row))
% useMZ = A Gx1 vector of the m/z values corresponding to the rows of the spectra
% matrices contained in useSPEC.  See the explanation above. 
% clusterTimeBin = a KxN matrix containing the temporal counts for each cluster.
% Each column is the data for a single cluster in the cluster order corresponding to outWM.
% Each row corresponds to a time bin. 
% TimeBin = A Kx1 vector which is the mid point of each time bin to plot
% the temporals in clusterTimeBin against.
% clusterSizeBin = a JxN matrix containing the number fraction by size for each cluster.
% Each column is the data for a single cluster in the cluster order corresponding to outWM.
% Each row corresponds to a size bin. 
% SizeeBin = A Jx1 vector which is the mid point of each size bin to plot
% the number fraction in clusterSizeBin against.
%
% OUTPUT
% dendoFig = a figure handle to be used in clickDendo for user exploration.
% dendoFig is sets up a dendogram created by grouping spectra using
% hierachichal clustering.
%
% There are a wide number of plots that the user is prompted to select
% plotCS: line plot of average mz spectrum for selected particles
% plotPS: line plot of all mz spectra overlaid.  NOTE THIS PLOT IS VERY
% TIME/MEMORY CONSUMING
% plotHM: linear heatmap of mz spectra
% plotHMlog: log10 heatmap of mz spectra
% plotSize: relative number fraction over size of selected particles 
% plotTime: relative temporal of selected particles
% plotSim: scatter plot of dot product1 (avg cluster spectra . particle
% spectra) and dot product2 (avg all selected spectra . particle spectra)

global PARTidMat PARTdataMat dendoSPEC dendoLogSPEC dendoSPECavg dendoSIM dendoMZ dendoTime dendoSize dendoTimeBin dendoSizeBin linkCluster dendoCOUNTS Xlines Ylines finalbaseCluster lineCluster plotTrack rownum 

%% check inputs
if nargin < 2 || nargin > 9 || nargout > 1
  error('Call as dendoFig = setup_clickHist(outWM,useSetup,useSPEC,useMZ,clusterTimeBin,TimeBin,clusterSizeBin,SizeBin)')
end

if exist('useSPEC') && ~isempty(useSPEC)
    if ~exist('useMZ') || isempty(useMZ)
        error('Spectra data are supplied.  Must supply m/z vector for plotting purposes')
    else
        if length(useMZ) ~= size(useSPEC{1},1)
            error('length of m/z vector must equal number of rows in each cell of useSPEC')
        elseif size(useMZ,2) > 1
            dendoMZ = useMZ';
        else
            dendoMZ = useMZ;            
        end
    end
    if length(useSPEC) ~= size(outWM,2)
        error('number of spectra cells (useSPEC) must match number of weight matrices (outWM columns)')
    end
else
    warning('Spectra data are not supplied.  Note will not be able to create spectra plots')    
end

if useSetup ~= 0 && useSetup ~= 1
    error('useSetup must be either 0 or 1')
end

if exist('clusterTimeBin') && ~isempty(clusterTimeBin)
    if ~exist('TimeBin') || isempty(TimeBin)
        error('TimeBin is not supplied.  Must supply TimeBin vector for plotting purposes')
    else
        if length(TimeBin) ~= size(clusterTimeBin,1)
            error('length of TimeBin must equal number of rows in in clusterTimeBin')     
        end
    end
    if size(clusterTimeBin,2) ~= size(outWM,2)
        error('number of columns of clusterTimeBin must match number of weight matrices (outWM columns)')
    end
else
    warning('Time data are not supplied.  Note will not be able to create time plots')    
end

if exist('clusterSizeBin') && ~isempty(clusterSizeBin)
    if ~exist('SizeBin') || isempty(SizeBin)
        error('SizeBin is not supplied.  Must supply SizeBin vector for plotting purposes')
    else
        if length(SizeBin) ~= size(clusterSizeBin,1)
            error('length of SizeBin must equal number of rows in in clusterSizeBin')           
        end
    end
    if size(clusterSizeBin,2) ~= size(outWM,2)
        error('number of columns of clusterSizeBin must match number of weight matrices (outWM columns)')
    end
else
    warning('Size data is not supplied.  Note will not be able to create size plots')    
end

%% do hierarchical clustering analysis
distCluster = pdist(outWM');
linkCluster = linkage(distCluster);

%% get input from user for setting up plots
if useSetup == 1 %read in file from user with everything all setup
    userSetupDendo
else %get user inputs one by one 
    
    % checks spectra matrix
    if exist('useSPEC')
        if ~isempty(useSPEC)
            plotSpec = input('Do you want to make spectra plots? 0-no, 1-yes  ');
        else
            plotSpec = 0;
        end
    else
        plotSpec = 0;
    end
    
    if plotSpec
        plotPS = input('Plot all particle spectra overlaid? 0-no, 1-yes  ');
        plotCS = input('Plot average spectra of each cluster/branch overlaid? 0-no, 1-yes  ');
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
    
    if exist('clusterTimeBin')
        if ~isempty(clusterTimeBin)
            plotTime = input('Do you want to plot clusters over time? 0-no 1-yes   ');
        else
            plotTime = 0;
        end
    else
        plotTime = 0;
    end
    
    if exist('clusterSizeBin')
        if ~isempty(clusterSizeBin)
            plotSize = input('Do you want to plot clusters size number fraction? 0-no 1-yes   ');
        else
            plotSize = 0;
        end
    else
        plotSize = 0;
    end
    
    %ask user if want to show titles and legends for plots
    plotRender = input('\nDo you want to do fast rendering (no titles, labels, legends are put onto plots)?\n This may give some speedups of graphics 0-no 1-yes  ');
    plotTrack = [plotPS plotCS plotHM plotHMlog plotTime plotSize plotSim plotRender];
end

%% create initial histogram plot
dendoFig = figure;
[Dhandle, IDXdendo, XLabel] = dendrogram(linkCluster,0);
set(Dhandle,'Color',[0 0 0]);

%% get data ready for plotting
%if plotting all plots in a single figure set up subplots
%set number of rows and columns in subplot
numSPECplot = plotPS+plotHM+plotHMlog+plotCS+plotTime;
numPlot = numSPECplot+plotTime+plotSize; %number of plots to combine
if plotSize || plotSim
    colnum = 1;
    rownum = numSPECplot+1;
else
    colnum = 1;
    rownum = numSPECplot;
end

%get # of particles in clusters
numCluster = size(outWM,2);
dendoCOUNTS = clustCOUNTS;

%set up spectra data
if plotHMlog
    dendoLogSPEC = cellfun(@log10, useSPEC, 'UniformOutput',0);
end

if plotCS || plotSim
    dendoSPECavg = cellfun(@(x) mean(x,2), useSPEC, 'UniformOutput',0);
end
    
if plotSim
    dendoSIM = cellfun(@(x,y) x'*y, useSPEC, dendoSPECavg, 'UniformOutput',0);
end
    
if ~plotHM && ~plotSim
    clear useSPEC
else
    dendoSPEC = useSPEC;
end

% set up time bins
if plotTime
    dendoTimeBin = TimeBin;
    dendoTime = clusterTimeBin;
end

%set up size bins
if plotSize
    dendoSizeBin = SizeBin;
    dendoSize = clusterSizeBin;
end

%% get indexes to be able to interact with dendogram
%get initial cluster labels on dendogram using data drom dendogram output
maxCluster = numCluster+length(Dhandle);
baseCluster = cell(1,maxCluster);
for i = 1:numCluster
    baseCluster{i} = i;
end

for i = (numCluster+1):(maxCluster)
    baseCluster{i} = linkCluster(i-numCluster,1:2);
end

origbaseCluster = baseCluster;
%rewrite dendo clusters using clusters orignally in art2a
%this is used to reference to particles for more plots later
for i = numCluster+1:maxCluster
    allBase = baseCluster{i} > numCluster; %check if any clusters listed are not base clusters originally in art2a
    if any(allBase)
        tmp = baseCluster{i}(allBase);
        addClust = cell2mat(baseCluster(tmp));
        baseCluster{i} = [baseCluster{i}(~allBase),addClust];
    end
end

%find out the base clusters for the left and right arm of each branch
finalbaseCluster = cell(1,length(Dhandle));
for i = numCluster+1:maxCluster
    finalbaseCluster{i-numCluster} = cell(1,2);
    tmpClust1 = baseCluster{origbaseCluster{i}(1)};
    tmpClust2 = baseCluster{origbaseCluster{i}(2)};
    [~, Idx1, IdxC1] = intersect(XLabel,tmpClust1,'stable');
    [~, Idx2, IdxC2] = intersect(XLabel,tmpClust2,'stable');
    min1 = min(Idx1);
    min2 = min(Idx2);
    if min1 < min2
        finalbaseCluster{i-numCluster}{1} = tmpClust1(IdxC1);
        finalbaseCluster{i-numCluster}{2} = tmpClust2(IdxC2);        
    else
        finalbaseCluster{i-numCluster}{2} = tmpClust1(IdxC1);
        finalbaseCluster{i-numCluster}{1} = tmpClust2(IdxC2);         
    end
end

%this is used to track lines on dendogram to recolor w 
lineCluster = origbaseCluster(1:length(Dhandle));
for i = numCluster+1:maxCluster
    allBase = origbaseCluster{i} > numCluster; %check if any clusters listed are not base clusters originally in art2a
    if any(allBase)
        tmp = origbaseCluster{i}(allBase);
        tmp2 = [i tmp];
        while any(allBase)
            addClust = cell2mat(origbaseCluster(tmp));
            allBase = addClust > numCluster;
            tmp = addClust(allBase);
            tmp2 = [tmp2,tmp];
        end
        lineCluster{i-numCluster} = tmp2-numCluster;
    end
end

%% get location of all lines on dendogram
Xlines = zeros(length(Dhandle),4);
Ylines = zeros(length(Dhandle),4);

for i = 1:length(Dhandle)
    tmp = get(Dhandle(i));
    Xlines(i,:) = tmp.XData;
    Ylines(i,:) = tmp.YData;
end

%reorder X/Ylines (want min X(:,1:2) and max X(:,3:4)
reorderX = Xlines(:,1) > Xlines(:,3);
Xlines(reorderX,[4 3 2 1]) = Xlines(reorderX,[1 2 3 4]);
Ylines(reorderX,[4 3 2 1]) = Ylines(reorderX,[1 2 3 4]);

