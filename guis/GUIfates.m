function GUIfates(useMZdata,MZvector,partData,timeData,sizeData,clustData,clustRelation,inputColors,clustOrdered)
% Call as GUIfates(useMZdata,MZvector,partData,timeData,sizeData,clustData,clustRelation)
% GUIfates is a graphical user interface that allows users to explore
% spectra data based on spectra characteristics, size, time, and cluster
% data. The gui allows users to create new clusters/groupings of particles
% and output the particle ids of any existing clusters.
% The user interacts with the GUI through the uicontrols implemented
% in GUIfates.  The GUI sets up the following plots
% 1)heatmap of all particle spectra (each row is a particle)
% 2)scatter plot of all particle sizes or whatever data was input into sizeData, Panel B
% 3 )scatter plot of all particle times or whatever data was input into
% timeData, Panel A
% 4)scatter plot of particle cluster relation statistics coded by cluster
% identity using color, or whatever data was input into clustRelation,
% Panel C
%
% INPUT 
% useMZdata= a MxN matrix containing M spectra with N data points per
% spectra. The order of the particles generating the M spectra needs to be
% maintained for partData, timeData, sizeData, clustData, and
% clustRelation.  (the xth spectra of useMZdata needs to correspond to the
% same particle as the xth size in sizeData).  
% MZvector= a vector with length N that corresponds to the MZ values in
% useMZdata. The yth entry in MZvector corresponds to the yth column in useMZdata
% The values in MZ vector need to be evenly spaced (either monotonically increasing
% or decreasing).  This means that when going from positive to negative values MZ vector has to
% proceed like ..-4 -3 -2 -1 0 1 2 3....  0 must be included.  Thus
% useMZdata must also have a column corresponding to 0 as well
%
% "OPTIONAL" INPUTS
% NOTE: To be able to use this data for plotting and data analysis purposes it must be
% supplied here.  If you do not supply this data then use [].  Any empty
% variable will be written as a vector of length M of ones to use so that
% the plots don't crash.
% partData: a Mx2 matrix set of paritcle identifiers stored as
%       partData(:,1) = InstID
%       partData(:,2) = PartID
% The yth row (particle) of partData corresponds to the yth row (spectra)
% of useMZdata
% timeData: a vector of length M with any data for the
% particles. Data will be plotted in Panel A. It is suggested 
% to supply particle time data (in numeric form) as Panel A can display the
% x axis using date strings.
% sizeData: a vector of length M with any data for the particles. Data will be 
% plotted in Panel B. It is suggested to supply particle size data or total 
% mass spectra ion intensity as useful metrics.
% clustData: a vector of length M with the cluster/group indetifier for the
% particles in partData.  For example for a list of 10 particles that is
% made up of two groups clustData may look like [1 1 1 2 2 2 2 2 1 1].  In
% this case particles 1-3 and 9-10 belong to cluster 1, and particles 4-8
% belong to cluster 2.
% clustRelation: a vector of length M with any data for the particles. Data will be
% plotted in Panel C and color coded by clustData. It is suggested to supply 
% a cluster relation statistic. A cluster relation statistic is
% computed by the user and could be supplied by a a variety of
% calculations.  It is meant to show the relation of the particles to the 
% cluster it is assigned to. One option is to use the output from silhouette (built-in
% matlab function) or to take the dot product of each particle spectra with the
% average spectra of the cluster to which it belongs. GUIFates allows the
% user to recalculate the cluster statistic after changes to clustering are
% made and replot the values in Panel C.
% inputColors: a Kx3 matrix with K being the number of unique
% clusters/groups in clustData.  Each row is a color vector in Matlab.
% InputColors specifies the colors used to display the different clusters.
% clustOrdered: a vector of length K containing the unique clusters/groups in 
% clustData.  The order that the clusters/groups are listed in clustOrdered
% are the order that they will display in the GUI.  Otherwise the clusters
% default to displaying in numeric order.   
%
% OUTPUT
% "Output PIDs": The user can choose to output the clusters created in the
% current GUI to the workspace by using the "Output PIDs" button.
% The name of the variable is specified by the user in the text box labeled variable name.
% The user can either choose to output a single cluster or all clusters in the GUI. 
% If all clusters are output the output is a cell array with each cell containing the particles
% in the relevant cluster.  The id of the cluster is specified in the second 
% row of the cell array in this case. 

numPart = size(useMZdata,1); %number of particles 
uniqueCluster = unique(clustData);

%% check inputs
if nargin < 2
   error('Too few input arguments.');
elseif nargin > 9  
    error('Too many input arguments.');
end

if isempty(MZvector) || ~exist('MZvector','var')
    error('must input vector for mz values to plot as x axis');
else
    if length(MZvector) < size(useMZdata,2)
        error('length of mz vector must equal number of columns of useMZdata');
    end
end

if isempty(timeData) || ~exist('timeData','var')
    timeData = zeros(numPart,1);
else
    if size(timeData,1) < size(timeData,2)
        timeData = timeData';
    end
    if size(timeData,1) ~= numPart
        error('timeData vector needs length needs same number of rows as spectra data (same # of particles)');
    end
end

if isempty(sizeData) || ~exist('sizeData','var')
    sizeData = zeros(numPart,1);
else
    if size(sizeData,1) < size(sizeData,2)
        sizeData = sizeData';
    end
    if size(sizeData,1) ~= numPart
        error('sizeData vector needs length needs same number of rows as spectra data (same # of particles)');
    end
end

if isempty(clustData) || ~exist('clustData','var')
    clustData = zeros(numPart,1);
else
    if size(clustData,1) < size(clustData,2)
        clustData = clustData';
    end
    if size(clustData,1) ~= numPart
        error('clustData vector needs length needs same number of rows as spectra data (same # of particles)');
    end
end

if isempty(partData) || ~exist('partData','var')
    partData = zeros(numPart,2);
else
    if size(partData,1) < size(partData,2)
        partData = partData';
    end
    if size(partData,1) ~= numPart
        error('partData vector needs length needs same number of rows as spectra data (same # of particles)');
    end
end

if isempty(clustRelation) || ~exist('clustRelation','var')
    clustRelation = zeros(numPart,1);
else
    if size(clustRelation,1) < size(clustRelation,2)
        clustRelation = clustRelation';
    end
    if size(clustRelation,1) ~= numPart
        error('clustRelation vector needs length needs same number of rows as spectra data (same # of particles)');
    end
end

if exist('inputColors','var')
    if ~isempty(inputColors)
        if size(inputColors,1) ~= length(uniqueCluster)
            error('Number of rows of inputColors must equal the number of unique clusters in clustData');
        end
        if size(inputColors,2) ~= 3
            error('inputColors must be a Mx3 matrix. Each row must be a MATLAB rbg triplet');
        end
    end
end

if exist('clustOrdered','var')
    if ~isempty(clustOrdered)
        checkClust = setxor(uniqueCluster,clustOrdered);
        if ~isempty(checkClust)
            error('No clusters can exist in clustData that are not in clustOrdered, and vice versa');
        end
    end
end

%% create plots for GUI
%initialize figure
gf = figure('Visible','off');
set(gf,'units','normalized');
set(gf,'position',[.01 .07 .75 .83]);

%set up variables
tmpIDX = [1:numPart]'; %index for spectra when "sorting" is applied
origClustData = clustData; %save original data to be able to easily reset
origClustRelation = clustRelation; %save original data to be able to easily reset
displaySubset = 0; %variable determines if all clusters (0) or only a subset (1) is to be displayed
currentAxis = 'linear'; %variable determines whether to display spectra imagemap with linear or log scale color
useIDX = tmpIDX; %final index used to determine order of spectra/which spectra to display taking into account all 'sorting', 'filtering', 'cluster selection', etc
filterIDX = true(numPart,1); %variable to filter out spectra based on filter selections
timeLimits = [min(timeData) max(timeData)]; 
sizeLimits = [min(sizeData) max(sizeData)];
clustRLimits = [min(clustRelation) max(clustRelation)];
isZero = useMZdata == 0; %find all locations where data supplied is zero, will need to use this later when switching between log, linear color maps
redoClust = []; %this is a tracker for clusters that need averages recalculated

%set up colors for associating with cluster identity
colormap(hsv) %easier to see differences with this map
cMap = colormap;
if exist('inputColors','var') && ~isempty(inputColors)
    colorIDs = [inputColors; cMap(1:20:end,:); cMap(10:20:end,:); cMap(5:20:end,:); cMap(15:20:end,:)];
else
    colorIDs = [cMap(1:20:end,:); cMap(10:20:end,:); cMap(5:20:end,:); cMap(15:20:end,:)];
end

%create heatmap spectra plot
ax = axes('Parent',gf,'Units','normalized','Position',[0.2 0.17 0.70 0.5]);
colormap(jet) %easier to see differences with this map
%find max and min peak intensity values to limit colormap
[maxVal, maxIdx] = max(useMZdata(:));
[minVal, minIdx] = min(useMZdata(:));
minVal2 = min(useMZdata(~isZero));
%if for some reason all values are the same, then all values are probably
%zero. Need to adjust clim otherwise throws error
if maxVal == minVal
    maxVal = minVal + 1;
    hidelog = 1;
else
    hidelog = 0;
end
%find max and min for MZvector
minMZ = min(MZvector);
maxMZ = max(MZvector);
spectraMap = imagesc([MZvector(1) MZvector(end)],[1 numPart],useMZdata,[minVal maxVal]); %create heatmap plot
%check to make sure MZvector was input correclty
shouldReverseX = minMZ < 0 && maxMZ < 0 && maxMZ == MZvector(1); 
if shouldReverseX %reverse Xdirection if input by user backwards
    ax.XDir = 'reverse';
end
ax.YTick = [];

% create average spectra graph
avgSpectra = nanmean(useMZdata,1); %get average of all m/z data input
cx = axes('Parent',gf,'Units','normalized','Position',[0.2 0.82 0.7 0.15]);
averagePlot = plot(MZvector,avgSpectra);
cx.XTick = [];
if shouldReverseX %reverse x axis if necessary
    cx.XDir = 'reverse';
end

% create PanelA plot, time plot 
dx = axes('Parent',gf,'Units','normalized','Position',[0.05 0.17 0.073 0.5]);
timePlot = plot(timeData,tmpIDX,'.');
dx.XLim = [min(timeData) max(timeData)];
dx.XTickLabelRotation = 90;
dx.XAxisLocation = 'top';
dx.XAxis.TickLabelFormat = '%.6f';
set(dx,'YDir','reverse'); %do this to match heatmap y direction

% create PanelB plot, size plot
bx = axes('Parent',gf,'Units','normalized','Position',[0.125 0.17 0.073 0.5]);
sizePlot = scatter(sizeData,tmpIDX,'.');
axis tight
set(bx,'YDir','reverse'); %do this to match heatmap y direction
bx.YTick = [];

% setup cluster data
% get unique clusters, use an ordered version if supplied
if exist('clustOrdered','var') && ~isempty(clustOrdered)
    uniqueCluster = clustOrdered; 
end
origuniqueCluster = uniqueCluster; %save for resetting values
uniqueLabel = makeClusterLabel(uniqueCluster); %get label/value for each cluster to display in listboxes
origuniqueLabel = uniqueLabel; %save for resetting values
%set up color data for clusters and particles
colorLabel = zeros(length(uniqueCluster),3); %color for each cluster
clusterColor = zeros(size(useMZdata,1),3); %color for each datapoint
averageClust = zeros(length(uniqueCluster),size(useMZdata,2)); %average spectra for each cluster
cnt = 0; %tracker to see if have exhausted color list
for m = 1:length(uniqueCluster)
    cnt = cnt+1;
    if cnt > size(colorIDs,1)
        cnt = 1;
    end
    tmp = clustData == uniqueCluster(m); %find all particles of a certain cluster
    averageClust(m,:) = nanmean(useMZdata(tmp,:),1); %determine average cluster mass spectra
    colorLabel(m,:) = colorIDs(cnt,:); %select color for cluster
    clusterColor(tmp,:) = repmat(colorLabel(m,:),sum(tmp),1); %apply color data for all particles of certain cluster
end
origcolorLabel = colorLabel; %save for resetting values

%make cluster graph, Panel C
ex = axes('Parent',gf,'Units','normalized','Position',[0.9 0.17 0.05 0.5]);
clustPlot = scatter(clustRelation,tmpIDX,[],clusterColor,'.'); %create plot
set(ex,'YDir','reverse');
ex.YAxisLocation = 'right';

% create average cluster mass spectra graph
fx = axes('Parent',gf,'Units','normalized','Position',[0.2 0.67 0.7 0.15]);
ClustAveragePlot = plot(MZvector,averageClust(1,:),'Color',colorLabel(1,:));
fx.XTick = [];
fx.XTick = [];
if shouldReverseX
    fx.XDir = 'reverse';
end

% link axes for plots so everything scales together
linkaxes([ax,bx,dx,ex],'y');
linkaxes([ax,cx,fx],'x');

%axData is a structure that holds info on ax, importantly x and y limits
%that are using when zooming. For some reason these are not being updated
%properly when subsets of the data are displayed, even when xlim/ylim are
%reset (ax.YLim = [x y]) so zooming out doesn't work perfectly. This is a hack to force the
%ylim in the "appData" be correct to get the zooming functionality correct.
axData = localCreateViewInfoCOPY(ax);  
bxData = localCreateViewInfoCOPY(bx);  
dxData = localCreateViewInfoCOPY(dx);  
exData = localCreateViewInfoCOPY(ex);  

%% labels and user input boxes
%label particle data (often time, size, and cluster statistics) plots
ALabel = uicontrol('Parent',gf,'Style','text','String','A','BackgroundColor','black','ForegroundColor','red','FontWeight','bold','Units','normalized','Position',[0.106 0.646 .015 .02],'FontSize',10);
BLabel = uicontrol('Parent',gf,'Style','text','String','B','BackgroundColor','black','ForegroundColor','red','FontWeight','bold','Units','normalized','Position',[0.183 0.646 .015 .02],'FontSize',10);
CLabel = uicontrol('Parent',gf,'Style','text','String','C','BackgroundColor','black','ForegroundColor','red','FontWeight','bold','Units','normalized','Position',[0.933 0.646 .015 .02],'FontSize',10);

%display parameters label
diplayParam = uicontrol('Parent',gf,'Style','text','String','DISPLAY PARAMETERS','FontWeight','bold','ForegroundColor','red','Units','normalized','Position',[0.07 0.11 .15 .02]);

%button group display log or linear color bar
bigGraph = uibuttongroup('Parent',gf,'SelectionChangedFcn',@replotBigGraph,'Units','normalized','Position',[0.01,0.065,0.08,.04]);
linPlot = uicontrol(bigGraph,'Style','radiobutton','String','linear','Units','normalized','Position',[0.45 0.3 .85 .6],'TooltipString','Use linear colormap scheme for peak intensity');
logPlot = uicontrol(bigGraph,'Style','radiobutton','String','log','Units','normalized','Position',[0 0.3 .4 .6],'TooltipString','Use log10 colormap scheme for peak intensity');
if hidelog == 1
    logPlot.Visible = 'off';
end

%create checkbox to sort by cluster or not
clusterGroup = uicontrol('Parent',gf,'Style','checkbox','Callback',@groupByCluster,'String','Group by Cluster','Units','normalized','Position',[0.01,0.04,0.1,.02],'TooltipString','Check to display data grouped by cluster');

%create checkbox to display Panel A axis as time
ATime = uicontrol('Parent',gf,'Style','checkbox','Callback',@panelAdateTick,'String','Use datetick PanelA','Units','normalized','Position',[0.01,0.01,0.1,.02],'TooltipString','Check to display A x-axis as time, using datetick');

%create listbox to select cluster to plot in avg spectrum
spectraLabel = uicontrol('Parent',gf,'Style','text','String','Change Average Cluster Spectra','Units','normalized','Position',[0.195,.085,0.08,.035]);
clusterSelect = uicontrol('Parent',gf,'Style','listbox','String',uniqueLabel,'Callback',@changeCluster,'Units','normalized','Position',[0.21,0.026,0.05,.055],'TooltipString','Select average cluster spectra to display');
recalcAvg = uicontrol('Parent',gf,'Style','pushbutton','String','Recalc Cluster MS Avg','Callback',@redoAverage,'Units','normalized','Position',[0.195,0.005,0.08,.02],'TooltipString','Click to recalculate cluster average mass spectra');

%create control to limit the clusters displayed
selectCluster = uicontrol('Parent',gf,'Style','pushbutton','String','Display Select Clusters','Callback',@selectClusters,'Units','normalized','Position',[0.11,0.085,0.08,.02],'TooltipString','Click to display data only from clusters selected in box below');
clusterSelect5 = uicontrol('Parent',gf,'Style','listbox','String',uniqueLabel,'Max',length(uniqueCluster),'Units','normalized','Position',[0.125,0.026,0.05,.055]);
resetDisplay = uicontrol('Parent',gf,'Style','pushbutton','String','Display all Clusters','Callback',@resetDisplays,'Units','normalized','Position',[0.11,0.005,0.08,.02],'TooltipString','Click to display data from all clusters');

%sort parameters label
sortParam = uicontrol('Parent',gf,'Style','text','String','SORT PARAMETERS','ForegroundColor','red','FontWeight','bold','Units','normalized','Position',[0.32 0.11 .1 .02]);

%button group to select sort method
sortGraph = uibuttongroup('Parent',gf,'SelectionChangedFcn',@sortBigGraph,'Units','normalized','Position',[0.29,0.013,0.15,.1]);
noSort = uicontrol(sortGraph,'Style','radiobutton','String','noSort','Units','normalized','Position',[0.5 0.4 .4 .18],'TooltipString','Select for particles to be ordered as input into GuiFates');
clustSort = uicontrol(sortGraph,'Style','radiobutton','String','C (clusterStat)','Units','normalized','Position',[0.5 0.8 .4 .18],'TooltipString','Select to order particles using data displayed in panel C');
timeSort = uicontrol(sortGraph,'Style','radiobutton','String','A (time)','Units','normalized','Position',[0 .78 .4 .18],'TooltipString','Select to order particles using data displayed in panel A');
sizeSort = uicontrol(sortGraph,'Style','radiobutton','String','B (size)','Units','normalized','Position',[0 .4 .4 .18],'TooltipString','Select to order particles using data displayed in panel B');
mzSort = uicontrol(sortGraph,'Style','radiobutton','String','mz','Units','normalized','Position',[0 0.05 .4 .18],'TooltipString','Select to order particles using m/z value entered in box at right');

%user input to select mz value to sort by
mzValue = uicontrol('Parent',gf,'Style','edit','Units','normalized','Position',[0.33,0.018,0.07,.02],'Callback',@sortbymz,'TooltipString','Must be a value contained in mzVector supplied to GuiFates.');

%filter parameters label
doFilter = uicontrol('Parent',gf,'Style','pushbutton','String','APPLY FILTER','Units','normalized','Position',[0.45 0.12 .075 .02],'Callback',@filterData,'TooltipString','Click to display only data which falls within filters set below');
noFilter = uicontrol('Parent',gf,'Style','pushbutton','String','CLEAR FILTER','Units','normalized','Position',[0.535 0.12 .075 .02],'Callback',@clearFilter,'TooltipString','Click to display data without filters');

%button group to select filter method
% filterGraph = uibuttongroup('Parent',gf,'SelectionChangedFcn',@filterLabels,'Units','normalized','Position',[0.45,0.055,0.13,.05]);
timeFilter = uicontrol('Parent',gf,'Style','checkbox','String','A','Units','normalized','Position',[0.45 0.085 .05 .02],'TooltipString','Check to filter by data displayed in panel A');
sizeFilter = uicontrol('Parent',gf,'Style','checkbox','String','B','Units','normalized','Position',[0.45 0.06 .05 .02],'TooltipString','Check to filter by data displayed in panel B');
clustFilter = uicontrol('Parent',gf,'Style','checkbox','String','C','Units','normalized','Position',[0.45 0.035 .05 .02],'TooltipString','Check to filter by data displayed in panel C');
mzFilter = uicontrol('Parent',gf,'Style','checkbox','String','mz','Units','normalized','Position',[0.45 0.01 .05 .02],'TooltipString','Check to filter by an m/z value');

%all the stuff for filtering
%text boxes to select max and min for filter method
sizeMin = uicontrol('Parent',gf,'Style','edit','Units','normalized','Position',[0.51 0.06, .04 .02]);
sizeMax = uicontrol('Parent',gf,'Style','edit','Units','normalized','Position',[0.6,.06,0.04,.02]);
clustMin = uicontrol('Parent',gf,'Style','edit','Units','normalized','Position',[0.51 0.035, .04 .02]);
clustMax = uicontrol('Parent',gf,'Style','edit','Units','normalized','Position',[0.6,.035, 0.04,.02]);
timeMin = uicontrol('Parent',gf,'Style','edit','Units','normalized','Position',[0.51 0.085, .04 .02]);
timeMax = uicontrol('Parent',gf,'Style','edit','Units','normalized','Position',[0.6,.085,0.04,.02]);
mzMin = uicontrol('Parent',gf,'Style','edit','Units','normalized','Position',[0.51 0.01 .04 .02]);
mzMax = uicontrol('Parent',gf,'Style','edit','Units','normalized','Position',[0.6,.01,0.04,.02]);
%labels for max and min columns
maxString = sprintf('Set max value for filter in boxes below. \n Use value less than max value of data, provided at right.');
minString = sprintf('Set min value for filter in boxes below. \n Use value greater than min value of data, provided at right.');
MaxLabel = uicontrol('Parent',gf,'Style','text','String','max','Units','normalized','Position',[0.61,.107,0.02,.015],'TooltipString',maxString);
MinLabel = uicontrol('Parent',gf,'Style','text','String','min','Units','normalized','Position',[0.52,.107,0.02,.015],'TooltipString',minString);
%data values for max and min
sizeMinV = uicontrol('Parent',gf,'Style','text','String',num2str(sizeLimits(1)),'Units','normalized','Position',[0.552,.06,0.04,.02],'TooltipString','min value of data in panel B');
sizeMaxV = uicontrol('Parent',gf,'Style','text','String',num2str(sizeLimits(2)),'Units','normalized','Position',[0.642,.06,0.04,.02],'TooltipString','max value of data in panel B');
clustMinV = uicontrol('Parent',gf,'Style','text','String',num2str(clustRLimits(1)),'Units','normalized','Position',[0.552,.035,0.04,.02],'TooltipString','min value of data in panel C');
clustMaxV = uicontrol('Parent',gf,'Style','text','String',num2str(clustRLimits(2)),'Units','normalized','Position',[0.642,.035,0.04,.02],'TooltipString','max value of data in panel C');
timeMinV = uicontrol('Parent',gf,'Style','text','String',num2str(timeLimits(1)),'Units','normalized','Position',[0.552,.085,0.04,.02],'TooltipString','min value of data in panel A');
timeMaxV = uicontrol('Parent',gf,'Style','text','String',num2str(timeLimits(2)),'Units','normalized','Position',[0.642,.085,0.04,.02],'TooltipString','max value of data in panel A');
mzMinV = uicontrol('Parent',gf,'Style','text','Visible','off','String','','Units','normalized','Position',[0.552,.01,0.04,.02],'TooltipString','min value of m/z specified');
mzMaxV = uicontrol('Parent',gf,'Style','text','Visible','off','String','','Units','normalized','Position',[0.642,.01,0.04,.02],'TooltipString','max value of m/z specified');
%user input to select mz value to sort by
mzValue2 = uicontrol('Parent',gf,'Style','edit','Units','normalized','Position',[0.475,.01,0.03,.02],'Callback',@filterbymz,'TooltipString','Must be a value contained in mzVector supplied to GuiFates.');

%create controls to combine clusters 
combineCluster = uicontrol('Parent',gf,'Style','pushbutton','String','Combine Clusters','Callback',@combineClusters,'Units','normalized','Position',[0.69,0.11,0.09,.02],'TooltipString','Click to combine clusters listed in box below');
clusterSelect1 = uicontrol('Parent',gf,'Style','edit','Units','normalized','Position',[0.7125,0.045,0.05,.06],'Max',2,'Min',0,'TooltipString','List clusters to combine into single cluster. e.g. 1,3,5..');
resetCluster = uicontrol('Parent',gf,'Style','pushbutton','String','Reset Original Clusters','Callback',@resetClusters,'Units','normalized','Position',[0.69,0.022,0.09,.02],'TooltipString','Click to reset cluster assignments to original values input to GuiFates');

%create controls to split an individual cluster
separateCluster = uicontrol('Parent',gf,'Style','pushbutton','String','Split Single Cluster','Callback',@separateClusters,'Units','normalized','Position',[0.8,0.11,0.07,.02],'TooltipString','Click to select location to split cluster specified in box below.');
clusterSelect3 = uicontrol('Parent',gf,'Style','listbox','String',uniqueLabel,'Units','normalized','Position',[0.81,0.045,0.05,.06],'TooltipString','Select cluster to split into two');
separatePrompt = uicontrol('Parent',gf,'Style','text','Visible','off','String','Select where to split!','Units','normalized','Position',[0.8,.01,0.07,.03]);

%create controls to split all displayed clusters
separateAll = uicontrol('Parent',gf,'Style','pushbutton','String','<html>  Split All Displayed <br> Clusters','Callback',@separateData,'Units','normalized','Position',[0.8,0.001,0.07,.04],'TooltipString','Click to select location to split all clusters displayed. One new cluster will be made.');
separateDir = uibuttongroup('Parent',gf,'Units','normalized','Position',[0.872,0.001,0.053,.04]);
belowLine = uicontrol(separateDir,'Style','radiobutton','String','below','Units','normalized','Position',[0.02 0.02 .95 .4],'TooltipString','Create new cluster grouping everything below cursor');
aboveLine = uicontrol(separateDir,'Style','radiobutton','String','above','Units','normalized','Position',[0.02 0.5 .95 .4],'TooltipString','Create new cluster grouping everything above cursor');

%create control to recalculate cluster statistic
recalcClusterRel = uicontrol('Parent',gf,'Style','pushbutton','String','Recalculate Cluster Stats','Callback',@recalcClust,'Units','normalized','Position',[0.69,0.001,0.09,.02],'TooltipString','Click to recalculate cluster statistics (panel C) using current cluster assignments.');

%create controls to outputPIDs
outputLabel = uniqueLabel;
outputLabel{length(uniqueLabel)+1} = 'ALL';
outputPID = uicontrol('Parent',gf,'Style','pushbutton','String','Output PIDs to Workspace','Callback',@outputPIDs,'Units','normalized','Position',[0.89,0.11,0.1,.02],'TooltipString','Click to ouput pid list to workspace for clusters selected in box below.');
clusterSelect4 = uicontrol('Parent',gf,'Style','listbox','String',outputLabel,'Units','normalized','Position',[0.89,0.045,0.05,.06],'TooltipString','Select clusters to ouput to workspace');
varName = uicontrol('Parent',gf,'Style','edit','Units','normalized','Position',[0.945,.055,0.05,.02],'TooltipString','Provide a name for variable output to workspace');
varLabel = uicontrol('Parent',gf,'Style','text','String','variable name','Units','normalized','Position',[0.945,.076,0.05,.02]);

gf.Visible = 'on';

%% gui functions
    %function to convert PanelA x-axis between datenum (numeric) and
    %datetime(string) tick labels
    function panelAdateTick(source,callbackdata)
        axes(dx)
        delete(timePlot);
        if ATime.Value == 1
            timeData = datetime(timeData,'ConvertFrom','datenum'); %get new time value
            timePlot = plot(timeData(useIDX),1:length(useIDX),'.'); %create plot
            dx.XAxis.TickLabelFormat = 'dd-MM HH:mm'; %set tick label format
        else
            timeData = datenum(timeData); %get new time value
            timePlot = plot(timeData(useIDX),1:length(useIDX),'.'); %create plot
            dx.XAxis.TickLabelFormat = '%.6f'; %set tick label format
        end
        %format axes
        dx.XTickLabelRotation = 90;
        dx.XAxisLocation = 'top';
        set(dx,'YDir','reverse'); %do this to match heatmap y direction
        dx.XLim = [min(timeData(useIDX)) max(timeData(useIDX))];
    end

%----------------------------------------------------%
    %function to filter data based on user selection
    function filterData(source,callbackdata)
       %do filter
        filterIDX = true(numPart,1);
        if mzFilter.Value == 1
            filterValue = mzValue2.String;
            filterValue = str2double(filterValue);
            filterValue = find(MZvector == filterValue);
            filterIDX = useMZdata(:,filterValue) < str2double(mzMax.String) & useMZdata(:,filterValue) > str2double(mzMin.String);
        end
        if sizeFilter.Value == 1
            filterIDX = filterIDX & (sizeData < str2double(sizeMax.String) & sizeData > str2double(sizeMin.String));
        end
        if timeFilter.Value == 1
            filterIDX = filterIDX & (timeData < str2double(timeMax.String) & timeData > str2double(timeMin.String));
        end
        if clustFilter.Value == 1
            filterIDX = filterIDX & (clustRelation < str2double(clustMax.String) & clustRelation > str2double(clustMin.String));
        end
        if mzFilter.Value == 0 && timeFilter.Value == 0 && sizeFilter.Value == 0 && clustFilter.Value == 0
            filterIDX = true(numPart,1);
        end
        %update plots
        plotHelper;
    end

%----------------------------------------------------%
    function clearFilter(source,callbackdata)
        filterIDX = true(numPart,1);
        plotHelper
    end

    %funtion to determine filter limits for mz if mz is reset when mz
    %button is already selected
    function filterbymz(source,callbackdata)
            filterValue = mzValue2.String;
            filterValue = str2double(filterValue);
            filterValue = find(MZvector == filterValue);
            mzMaxV.String = num2str(max(useMZdata(:,filterValue)));
            mzMinV.String = num2str(min(useMZdata(:,filterValue)));
            mzMaxV.Visible = 'on';
            mzMinV.Visible = 'on';
    end

%----------------------------------------------------%
    %function to determine limits for filter type selected and display them
    function filterLabels(source,callbackdata)
        currentButton = filterGraph.SelectedObject.String; %find filter type selected
        %find filter limits
        if strcmp(currentButton,'mz')
            filterValue = mzValue2.String;
            filterValue = str2double(filterValue);
            filterValue = find(MZvector == filterValue);
            MaxValue.String = num2str(max(useMZdata(:,filterValue)));
            MinValue.String = num2str(min(useMZdata(:,filterValue)));
            MaxValue.Visible = 'on';
            MinValue.Visible = 'on';
        elseif strcmp(currentButton,'time');
            MaxValue.String = num2str(timeLimits(2));
            MinValue.String = num2str(timeLimits(1));
            MaxValue.Visible = 'on';
            MinValue.Visible = 'on';
        elseif strcmp(currentButton,'size');
            MaxValue.String = num2str(sizeLimits(2));
            MinValue.String = num2str(sizeLimits(1));
            MaxValue.Visible = 'on';
            MinValue.Visible = 'on';
        elseif strcmp(currentButton,'clusterStat');
            MaxValue.String = num2str(clustRLimits(2));
            MinValue.String = num2str(clustRLimits(1));
            MaxValue.Visible = 'on';
            MinValue.Visible = 'on';        
        else
            MaxValue.Visible = 'off';
            MinValue.Visible = 'off';
        end
    end

%----------------------------------------------------%
    %function to recalculate cluster relation
    function recalcClust(source,callbackdata)
        %note the user could choose to use another method for calculating
        %the particle to cluster relation statistic. The method here is
        %chosen over silhouette because it is much faster
        if ~isempty(redoClust)
            redoAverage
        end
        for i = 1:length(uniqueCluster)
            cidx = clustData == uniqueCluster(i);
            tmp = (averageClust(i,:)/norm(averageClust(i,:)))*normc(useMZdata(cidx,:)');
            clustRelation(cidx) = tmp;
        end
        set(clustPlot,'XData',clustRelation(useIDX));
    end

%----------------------------------------------------%
    %function to display all clusters
    function resetDisplays(source,callbackdata)
        displaySubset = 0;
        plotHelper;
    end

%----------------------------------------------------%    
    %function to select what clusters to display
    function selectClusters(source,callbackdata)
        displaySubset = 1;
        plotHelper;
    end 

%----------------------------------------------------%    
    %function to output PID lists to workspace
    function outputPIDs(source,callbackdata)
        findCluster = clusterSelect4.Value; %get cluster to output
        nameVar = varName.String; %get name for variable in workspace
        if findCluster <= length(uniqueCluster) %for single cluster selection
            outPID = partData(clustData == uniqueCluster(findCluster),:);
        else %to output all clusters
            outPID = cell(2,length(uniqueCluster));
            for i = 1:length(uniqueCluster)
                outPID{1,i} = partData(clustData == uniqueCluster(i),:);
                outPID{2,i} = uniqueCluster(i);
            end
        end
        assignin('base',nameVar,outPID); %output data to workspace
    end

%----------------------------------------------------%
    %function to split a single cluster into two
    function separateClusters(source,callbackdata)
        axes(ex);
        separatePrompt.Visible = 'on'; %indicate to user to make split
        [~,y] = ginput(1); %get user input on where to make split
        separatePrompt.Visible = 'off';
        y = round(y);
        %get indexes for splitting
        tmp = clustData(useIDX) == uniqueCluster(clusterSelect3.Value);
        sum1 = sum(tmp);
        tmp(1:y) = 0;
        if any(tmp) && sum(tmp) ~= sum1 %make sure user input creates a split
            %get new color for new cluster
            [~,indexColor] = intersect(colorIDs,colorLabel(end,:),'rows'); 
            indexColor = indexColor+1; 
            if indexColor > size(colorIDs,1)
                indexColor = 1;
            end
            %add new color to data
            colorLabel = [colorLabel; colorIDs(indexColor,:)];
            newC = max(uniqueCluster)+1; %new cluster index
            clustData(useIDX(tmp)) = newC; %set new cluster ids for each particle
            clusterColor(useIDX(tmp),:) = repmat(colorLabel(end,:),sum(tmp),1); %set new color for each particle
            %update data with new cluster
            uniqueCluster(length(uniqueCluster)+1) = newC;
            uniqueLabel{length(uniqueCluster)} = num2str(newC);
            outputLabel = uniqueLabel;
            outputLabel{length(uniqueLabel)+1} = 'ALL';
            redoClust = [redoClust; uniqueCluster(clusterSelect3.Value); newC];
            averageClust(length(uniqueCluster),:) = zeros(1,size(useMZdata,2)); %just use placeholder for average cluster
            %update display
            clusterSelect3.String = uniqueLabel;
            clusterSelect.String = uniqueLabel;
            clusterSelect5.String = uniqueLabel;
            clusterSelect4.String = outputLabel;
            clusterSelect5.Max = length(uniqueLabel);
            plotHelper; %update plot
        end
    end

%----------------------------------------------------%
    %function to split all displayed data
    function separateData(source,callbackdata)
        axes(ex);
        separatePrompt.Visible = 'on';
        [~,y] = ginput(1); %get user input on where to make split
        separatePrompt.Visible = 'off';
        y = round(y);
        if (y > 1 && y < length(useIDX)) %make sure user clicked on the plot (in y axis)          
            %set up color indexes/labels
            [~,indexColor] = intersect(colorIDs,colorLabel(end,:),'rows'); %new color for new cluster
            indexColor = indexColor+1;
            if indexColor > size(colorIDs,1)
                indexColor = 1;
            end
            colorLabel = [colorLabel; colorIDs(indexColor,:)];
            
            %create index for particles to put into new cluster
            newC = max(uniqueCluster)+1; %new cluster index
            tmp = true(length(useIDX),1); 
            if belowLine.Value
                tmp(1:y) = 0; %leave everything above the line alone
            else
                tmp((y+1):end) = 0; %leave everything below the line alone
            end
            
            %identify all clusters that need averages to be redone
            renameList = unique(clustData(useIDX(tmp))); %all cluster ids that are being affected
            renameIDX = find(ismember(uniqueCluster,renameList)); %find location of clusterIDs to rename in uniqueCluster list
    
            %new cluster/color data for particles
            clustData(useIDX(tmp)) = newC; %set new cluster ids for each particle
            clusterColor(useIDX(tmp),:) = repmat(colorLabel(end,:),sum(tmp),1); %set new color for each particle
            %update data with new cluster
            uniqueCluster(length(uniqueCluster)+1) = newC;
            uniqueLabel{length(uniqueCluster)} = num2str(newC);
            outputLabel = uniqueLabel;
            outputLabel{length(uniqueLabel)+1} = 'ALL';

            %update display
            clusterSelect3.String = uniqueLabel;
            clusterSelect.String = uniqueLabel;
            clusterSelect5.String = uniqueLabel;
            clusterSelect4.String = outputLabel;
            clusterSelect5.Max = length(uniqueLabel);           
            plotHelper; %update plot

            %redo averages for all clusters that were altered
            redoClust = [redoClust; renameList; newC];
            averageClust(length(uniqueCluster),:) = zeros(1,size(useMZdata,2));
        end
    end

%----------------------------------------------------%
    %function to redo average cluster mass spectra
    function redoAverage(source,callbackdata)
        %get list of clusters that need averages to be redone
        redoClust = unique(redoClust); %get unique list, there could be repeats
        clustIDX = find(ismember(uniqueCluster,redoClust)); %get position of clusters to be redone in uniqueCluster list
        if logPlot.Value == 1 %if log data is being used convert back to linear
            for i = 1:length(clustIDX)
                tmpMZ = useMZdata(clustData == uniqueCluster(clustIDX(i)),:);
                tmpMZ = 10.^tmpMZ;
                tmpMZ = tmpMZ*minVal2;
                averageClust(clustIDX(i),:) = nanmean(tmpMZ,1);
            end
        else
            for i = 1:length(clustIDX)
                averageClust(clustIDX(i),:) = nanmean(useMZdata(clustData == uniqueCluster(clustIDX(i)),:),1);
            end
        end
        %reset display if necessary
        if ismember(clusterSelect.Value, clustIDX) 
            changeCluster(clusterSelect);
        end
        redoClust = []; %reset tracker
    end

%----------------------------------------------------%
    %function to return to original input cluster 
    function resetClusters(source,callbackdata)
        %create dialog box to verify the reset
        d = dialog('Position',[300 300 300 150],'Name','My Dialog');
        txt = uicontrol('Parent',d,'Style','text','Units','normalized','Position',[.05 .6 .85 .3],...
            'String',{'Are you sure you want to reset clusters?', 'Consider outputting your clusters before resetting.'});
        
        yesSel = uicontrol('Parent',d,'Style','pushbutton','String','Yes Reset! ','Units','normalized','Position',[0.1 0.35 .4 .2],'Callback',@userYes);
        noSel = uicontrol('Parent',d,'Style','pushbutton','String','Nope','Units','normalized','Position',[0.6 0.35 .3 .2],'Callback',@userNo);
        
        function userYes(yesSel,callbackdata)
            ResetC = 1;
            delete(gcf);
        end
        
        function userNo(noSel,callbackdata)
            ResetC = 0;
            delete(gcf);
        end
        
        uiwait(d)
        %reset data
        if ResetC == 1
            uniqueCluster = origuniqueCluster;
            uniqueLabel = origuniqueLabel;
            colorLabel = origcolorLabel;
            clustData = origClustData;
            clustRelation = origClustRelation;
            displaySubset = 0;
            %recalculate some values
            for q = 1:length(uniqueCluster)
                tmp = clustData == uniqueCluster(q);
                clusterColor(tmp,:) = repmat(colorLabel(q,:),sum(tmp),1);
            end
            averageClust = zeros(length(uniqueCluster),size(useMZdata,2));
            redoAverage(1:length(uniqueCluster));
            %set all selections to 1 in case previous selections are no longer
            %viable are resetting
            clusterSelect2.Value = 1;
            clusterSelect1.Value = 1;
            clusterSelect.Value = 1;
            clusterSelect3.Value = 1;
            clusterSelect4.Value = 1;
            clusterSelect5.Value = 1;
            %update display
            uniqueLabel = makeClusterLabel(uniqueCluster);
            outputLabel = uniqueLabel;
            outputLabel{length(uniqueLabel)+1} = 'ALL';
            clusterSelect3.String = uniqueLabel;
%             clusterSelect1.String = uniqueLabel;
%             clusterSelect2.String = uniqueLabel;
            clusterSelect.String = uniqueLabel;
            clusterSelect4.String = outputLabel;
            clusterSelect5.String = uniqueLabel;
            clusterSelect5.Max = length(uniqueLabel);
            changeCluster(clusterSelect);
            plotHelper; %update plot
        end
    end

%----------------------------------------------------%
    %create cell array of strings from number vector
    function uniqueStr = makeClusterLabel(uniqueVector)
        if size(uniqueVector,1) > size(uniqueVector,2)
            uniqueStr = strtrim(cellstr(num2str(uniqueVector))');
        else
            uniqueStr = strtrim(cellstr(num2str(uniqueVector'))');
        end
    end

%----------------------------------------------------%
    %function to combine two clusters
    function combineClusters(source,callbackdata)
        %get clusters to combine
        clusterList = str2num(clusterSelect1.String); %#ok<ST2NM> %convert cluster list to numeric array
        %make sure all clusters listed are an option to combine
        NOTunqC = setdiff(clusterList, uniqueCluster); %will return any clusters in clusterList(from the user) that are not actual clusters at the time
        
        if isempty(NOTunqC) %everything in clusterList is good, start combining
            %determine how to combine clusters
            minCID = min(clusterList); %find min clusterID to combine
            minCIDX = find(uniqueCluster == minCID); %find location of min clusterID in uniqueCluster list
            renameList = clusterList(clusterList > minCID); %find all other clusterIDs that will be renamed
            renameIDX = ismember(uniqueCluster,renameList); %find location of clusterIDs to rename in uniqueCluster list
            
            %combine clusters
            tmp = ismember(clustData,renameList); %find data points to rename cluster
            clustData(tmp) = minCID; %rename to min cluster
            clusterColor(tmp,:) = repmat(colorLabel(minCIDX,:),sum(tmp),1); %fix color of particles to match min cluster
            redoClust = [redoClust; minCID]; %track clusters that have changed
            
            %get all clusters currently selected by cluster5 (clusters to
            %display
            cSelect = clusterSelect5.Value;
            cSelect = uniqueCluster(cSelect);
            
            %remove old cluster data
            averageClust(renameIDX,:) = [];
            colorLabel(renameIDX,:) = [];
            uniqueCluster(renameIDX) = [];
            
            %set all selections to min cluster in case previous selections are no longer
            %viable after resetting
            clusterSelect.Value = minCIDX;
            clusterSelect3.Value = minCIDX;
            clusterSelect4.Value = minCIDX;
            changeCluster(clusterSelect);
            
            %update display
            uniqueLabel(renameIDX) = [];
            outputLabel(renameIDX) = [];
            clusterSelect3.String = uniqueLabel;
            clusterSelect.String = uniqueLabel;
            clusterSelect4.String = outputLabel;
            clusterSelect5.String = uniqueLabel;
            clusterSelect5.Max = length(uniqueLabel);
            
            %reset cluster5 selection, keeping previous user input
            cSelect = setdiff(cSelect, renameList);
            clusterSelect5.Value = find(ismember(uniqueCluster, cSelect));
            
            plotHelper; %update plots
        else %user listed nonviable cluster
        end
    end

%----------------------------------------------------%
    %function to change the average cluster spectrum displayed
    function changeCluster(source,callbackdata)
        axes(ex);
        newClusterIDX = source.Value; %get new cluster to display
        ClustAveragePlot.YData = averageClust(newClusterIDX,:);
        ClustAveragePlot.Color = colorLabel(newClusterIDX,:);
    end

%----------------------------------------------------%
    %wrapper function to call plotHelper from callback
    function groupByCluster(source,callbackdata)
        plotHelper;
    end

%----------------------------------------------------%
    %function to sort data by mz value when mz value selection is changed
    function sortbymz(source,callbackdata)
        currentButton = sortGraph.SelectedObject.String; %only perform sort if mz sort is selected
        if strcmp(currentButton,'mz')
            %find value to sort by
            sortValue = source.String; 
            sortValue = str2double(sortValue);
            sortValue = find(MZvector == sortValue);
            [~, tmpIDX] = sortrows(useMZdata,sortValue); %sort spectra data
            plotHelper; %update plots
        end
    end

 %----------------------------------------------------%   
    %function to redo plots
    function plotHelper
        isChecked = clusterGroup.Value; %determine if display is grouped by cluster
        useIDX = tmpIDX(filterIDX(tmpIDX)); %index to use for displaying data
        CtmpIDX2 = clustData(useIDX); %index cluster data for sorting
        if ~isChecked
            lengthTMP = length(useIDX);
            if displaySubset == 1 %if not displaying all clusters
                %get index for particles belonging to selected clusters
                selectValues = clusterSelect5.Value;
                selectLogic = zeros(length(useIDX),1);
                for i = 1:length(selectValues)
                    selectLogic = selectLogic | CtmpIDX2 == uniqueCluster(selectValues(i));
                end
                useIDX = useIDX(selectLogic);
            end
        else %if display is grouped by cluster
            numTrack = 1;
            CtmpIDX = [];
            if displaySubset == 0 %if displaying all clusters
                %order indexes by clusters
                for i = 1:length(uniqueCluster);
                    tmp = useIDX(CtmpIDX2 == uniqueCluster(i));
                    lengthClust = size(tmp,1);
                    CtmpIDX(numTrack:(numTrack+lengthClust-1)) = tmp;
                    numTrack = numTrack+lengthClust;
                end
            else
                %order indexes by clusters for only select clusters
                selectValues = clusterSelect5.Value;
                for i = 1:length(selectValues)
                    tmp = useIDX(CtmpIDX2 == uniqueCluster(selectValues(i)));
                    lengthClust = size(tmp,1);
                    CtmpIDX(numTrack:(numTrack+lengthClust-1)) = tmp;
                    numTrack = numTrack+lengthClust;
                end
            end
            useIDX = CtmpIDX;
        end
        lengthTMP = length(useIDX); %determine number of particles to display
        %update plots
        spectraMap.CData = useMZdata(useIDX,:);
        spectraMap.YData = [1:lengthTMP];
        set(sizePlot,'XData',sizeData(useIDX),'YData',1:lengthTMP);
        set(clustPlot,'XData',clustRelation(useIDX),'YData',1:lengthTMP,'CData',clusterColor(useIDX,:));
        set(timePlot,'XData',timeData(useIDX),'YData',1:lengthTMP);
        ax.YLim = [1,lengthTMP];
        %this forces the zoom out function to work correctly. It's a hack.
        %Hopefully works for a while
        axData.YLim = [1, lengthTMP]; 
        bxData.YLim = [1, lengthTMP];
        dxData.YLim = [1, lengthTMP]; 
        exData.YLim = [1, lengthTMP]; 
        setappdata(ax,'matlab_graphics_resetplotview',axData);
        setappdata(bx,'matlab_graphics_resetplotview',bxData);
        setappdata(dx,'matlab_graphics_resetplotview',dxData);
        setappdata(ex,'matlab_graphics_resetplotview',exData);
    end

%----------------------------------------------------%
    %function to sort data based upon selection
    function sortBigGraph(source,callbackdata)
        currentButton = sortGraph.SelectedObject.String; %how to sort
        %do sorts
        if strcmp(currentButton,'mz')
            sortValue = mzValue.String;
            sortValue = str2double(sortValue);
            sortValue = find(MZvector == sortValue);
            [~, tmpIDX] = sortrows(useMZdata,sortValue);
        elseif strcmp(currentButton,'A (time)')
            [~,tmpIDX] = sort(timeData);
        elseif strcmp(currentButton,'B (size)')
            [~,tmpIDX] = sort(sizeData);
        elseif strcmp(currentButton,'C (clusterStat)')
            [~,tmpIDX] = sort(clustRelation);
        else
            tmpIDX = [1:size(useMZdata,1)]';
        end
        plotHelper %update plots
    end

%----------------------------------------------------%
    %function to switch between log and linear color axis for spectra heatmap
    %plot
    function replotBigGraph(source,callbackdata)
        currentAxis = source.SelectedObject.String; %determine if log or linear selected
        if strcmp(currentAxis,'log') %get log data if it hasn't been calculated yet
            useMZdata = useMZdata/(minVal2);
            useMZdata = log10(useMZdata);
            useMZdata(isZero) = 0;
            spectraMap.CData = useMZdata(useIDX,:); %update plot with log data
            caxis(ax,[0 useMZdata(maxIdx)]);
        else
%             useMZdata = linData;
            useMZdata = 10.^useMZdata;
            useMZdata = useMZdata*minVal2;
            useMZdata(isZero) = 0;
            spectraMap.CData = useMZdata(useIDX,:); %updata plot with linear data
            caxis(ax,[0 maxVal]);
        end
    end
end

% %----------------------------------------------------%
% %this is a copy of a subfunction found in resetplotview (internal matlab
% %function) that is used to create a structure to hold info on an axis
%     function viewinfo = localCreateViewInfoCOPY(hAxes)
%         axes_properties = {'DataAspectRatio',...
%             'CameraViewAngle',...
%             'PlotBoxAspectRatio',...
%             'CameraPosition',...
%             'CameraTarget',...
%             'CameraUpVector',...
%             'XLim',...
%             'YLim',...
%             'ZLim'};
%         
%         % Save the value of each axes property and its mode
%         for i = 1:numel(axes_properties)
%             current_prop = axes_properties{i};
%             current_mode = [axes_properties{i} 'Mode'];
%             viewinfo.(current_mode) = get(hAxes,current_mode);
%             % Only get properties in manual since getting them in auto can trigger
%             % an auto-calc
%             if strcmp(get(hAxes,current_mode),'manual')
%                 viewinfo.(current_prop) = get(hAxes,current_prop);
%             end
%         end
%         
%         [az, el] = view(hAxes);
%         viewinfo.View = [az, el];
%     end
