function GUIfates(useMZdata,MZvector,partData,timeData,sizeData,clustData,clustRelation,inputColors,clustOrdered)
% Call as GUIfates(useMZdata,MZvector,partData,timeData,sizeData,clustData,clustRelation)
% GUIfates is a graphical user interface that allows users to explore
% spectra data based on spectra characteristics, size, time, and cluster
% data. The gui allows users to create new clusters/groupings of particles
% and output the particle ids of any existing clusters.
% The user interacts with the GUI through the uicontrols implemented
% in GUIfates.  The GUI sets up the following plots
% 1)heatmap of all particle spectra (each row is a particle)
% 2)scatter plot of all particle sizes
% 3)scatter plot of all particle times
% 4)scatter plot of particle cluster relation statistics coded by cluster
% identity using color
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
% timeData: a vector of length M with time data (in numeric form) for the
% particles in partData
% sizeData: a vector of length M with size data for the particles in
% partData
% clustData: a vecotr of length M with the cluster/group indetifier for the
% particles in partData.  For example for a list of 10 particles that is
% made up of two groups clustData may look like [1 1 1 2 2 2 2 2 1 1].  In
% this case particles 1-3 and 9-10 belong to cluster 1, and particles 4-8
% belong to cluster 2.
% clustRelation: a vector of length M with the cluster relation statistic
% for the particles in partData.  This cluster relation statistic is
% computed by the user and could be supplied by a a variety of
% calculations.  It is meant to show the relation of the particles to the 
% cluster it is assigned to. One option is to use the output from silhouette (built-in
% matlab function) or to take the dot product of each particle spectra with the
% average spectra of the cluster to which it belongs.
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
tmpIDX = [1:numPart]';
origClustData = clustData;
origClustRelation = clustRelation;
logData = [];
displaySubset = 0;
currentAxis = 'linear';
useIDX = tmpIDX;
filterIDX = true(numPart,1);
currentCluster = [];
timeLimits = [min(timeData) max(timeData)];
sizeLimits = [min(sizeData) max(sizeData)];
clustRLimits = [min(clustRelation) max(clustRelation)];

%make heatmap spectra plot
ax = axes('Parent',gf,'Units','normalized','Position',[0.2 0.17 0.70 0.5]);
colormap(jet) %easier to see differences with this map
[maxVal, maxIdx] = max(useMZdata(:));
[minVal, minIdx] = min(useMZdata(:));
%if for some reason all values are the same, then all values are probably
%zero. Need to adjust clim otherwise throws error
if maxVal == minVal
    maxVal = minVal + 1;
    hidelog = 1;
else
    hidelog = 0;
end
minVal2 = min(useMZdata(useMZdata(:)>0));
minMZ = min(MZvector);
maxMZ = max(MZvector);
spectraMap = imagesc([MZvector(1) MZvector(end)],[1 numPart],useMZdata,[minVal maxVal]);
shouldReverseX = minMZ < 0 && maxMZ < 0 && maxMZ == MZvector(1);
if shouldReverseX
    ax.XDir = 'reverse';
end
ax.YTick = [];

%set up colors for associating with cluster identity
cMap = colormap;
if exist('inputColors','var') && ~isempty(inputColors)
    colorIDs = [inputColors; cMap(1:20:end,:); cMap(10:20:end,:); cMap(5:20:end,:); cMap(15:20:end,:)];
else
    colorIDs = [cMap(1:20:end,:); cMap(10:20:end,:); cMap(5:20:end,:); cMap(15:20:end,:)];
end

% create average spectra graph
avgSpectra = nanmean(useMZdata,1);
cx = axes('Parent',gf,'Units','normalized','Position',[0.2 0.82 0.7 0.15]);
averagePlot = plot(MZvector,avgSpectra);
cx.XTick = [];
if shouldReverseX
    cx.XDir = 'reverse';
end

% create time plot 
dx = axes('Parent',gf,'Units','normalized','Position',[0.05 0.17 0.073 0.5]);
timePlot = scatter(timeData,tmpIDX,'.');
datetick('x', 2);
dx.XTick = [min(timeData):(max(timeData)-min(timeData))/3:max(timeData)];
dx.XLim = [min(timeData) max(timeData)];
dx.XTickLabelRotation = 90;
dx.XAxisLocation = 'top';
set(dx,'YDir','reverse');

% create size graph
bx = axes('Parent',gf,'Units','normalized','Position',[0.125 0.17 0.073 0.5]);
sizePlot = scatter(sizeData,tmpIDX,'.');
axis tight
set(bx,'YDir','reverse');
bx.YTick = [];

% add cluster data
% get unique clusters, use an ordered version if supplied
if exist('clustOrdered','var') && ~isempty(clustOrdered)
    uniqueCluster = clustOrdered; 
end
origuniqueCluster = uniqueCluster; %save for resetting values
uniqueLabel = makeClusterLabel(uniqueCluster); %get label/value for each cluster
origuniqueLabel = uniqueLabel; %save for resetting values
colorLabel = zeros(length(uniqueCluster),3); %color for each cluster
clusterColor = zeros(size(useMZdata,1),3); %color for each datapoint
averageClust = zeros(length(uniqueCluster),size(useMZdata,2)); %average spectra for each cluster
% create matrix holding color for each particle
cnt = 0;
for m = 1:length(uniqueCluster)
    cnt = cnt+1;
    if cnt > size(colorIDs,1)
        cnt = 1;
    end
    tmp = clustData == uniqueCluster(m);
    averageClust(m,:) = nanmean(useMZdata(tmp,:),1);
    colorLabel(m,:) = colorIDs(cnt,:);
    clusterColor(tmp,:) = repmat(colorLabel(m,:),sum(tmp),1);
end
origcolorLabel = colorLabel; %save for resetting values
%make cluster graph
ex = axes('Parent',gf,'Units','normalized','Position',[0.9 0.17 0.05 0.5]);
clustPlot = scatter(clustRelation,tmpIDX,[],clusterColor,'.'); %create plot
set(ex,'YDir','reverse');
ex.YAxisLocation = 'right';

% create average spectra graph
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

%%
%display parameters label
diplayParam = uicontrol('Parent',gf,'Style','text','String','DISPLAY PARAMETERS','ForegroundColor','red','Units','normalized','Position',[0.03 0.11 .1 .02]);

%button group to sort spectra by time, size, spectra value
bigGraph = uibuttongroup('Parent',gf,'SelectionChangedFcn',@replotBigGraph,'Units','normalized','Position',[0,0.03,0.08,.05]);
linPlot = uicontrol(bigGraph,'Style','radiobutton','String','linear','Units','normalized','Position',[0.45 0.42 .85 .4]);
logPlot = uicontrol(bigGraph,'Style','radiobutton','String','log','Units','normalized','Position',[0 0.42 .4 .4]);
if hidelog == 1
    logPlot.Visible = 'off';
end

%create checkbox to sort by cluster or not
clusterLabel = uicontrol('Parent',gf,'Style','text','String','Group by Cluster','Units','normalized','Position',[0.08,.06,0.065,.02]);
clusterGroup = uicontrol('Parent',gf,'Style','checkbox','Callback',@groupByCluster,'Units','normalized','Position',[0.09,0.045,0.05,.02]);

%create listbox to select cluster to plot in avg spectrum
spectraLabel = uicontrol('Parent',gf,'Style','text','String','Change Average Cluster Spectra','Units','normalized','Position',[0.22,.11,0.07,.03]);
clusterSelect = uicontrol('Parent',gf,'Style','listbox','String',uniqueLabel,'Callback',@changeCluster,'Units','normalized','Position',[0.23,0.045,0.05,.06]);

%create control to limit the clusters displayed
selectCluster = uicontrol('Parent',gf,'Style','pushbutton','String','Clusters To Display','Callback',@selectClusters,'Units','normalized','Position',[0.14,0.11,0.08,.02]);
clusterSelect5 = uicontrol('Parent',gf,'Style','listbox','String',uniqueLabel,'Max',length(uniqueCluster),'Units','normalized','Position',[0.155,0.045,0.05,.06]);
resetDisplay = uicontrol('Parent',gf,'Style','pushbutton','String','Display all Clusters','Callback',@resetDisplays,'Units','normalized','Position',[0.14,0.022,0.08,.02]);

%sort parameters label
sortParam = uicontrol('Parent',gf,'Style','text','String','SORT PARAMETERS','ForegroundColor','red','Units','normalized','Position',[0.32 0.11 .1 .02]);

%button group to select sort method
sortGraph = uibuttongroup('Parent',gf,'SelectionChangedFcn',@sortBigGraph,'Units','normalized','Position',[0.3,0.055,0.13,.05]);
noSort = uicontrol(sortGraph,'Style','radiobutton','String','noSort','Units','normalized','Position',[0.28 0 .38 .4]);
clustSort = uicontrol(sortGraph,'Style','radiobutton','String','clusterStat','Units','normalized','Position',[0.28 0.42 .38 .82]);
timeSort = uicontrol(sortGraph,'Style','radiobutton','String','time','Units','normalized','Position',[0 .42 .25 .82]);
sizeSort = uicontrol(sortGraph,'Style','radiobutton','String','size','Units','normalized','Position',[0 0 .25 .4]);
mzSort = uicontrol(sortGraph,'Style','radiobutton','String','mz','Units','normalized','Position',[.7 0.42 .3 .82]);

%user input to select mz value to sort by
mzValue = uicontrol('Parent',gf,'Style','edit','Units','normalized','Position',[0.39,.06,0.038,.02],'Callback',@sortbymz);

%filter parameters label
doFilter = uicontrol('Parent',gf,'Style','pushbutton','String','APPLY FILTER','Units','normalized','Position',[0.45 0.11 .06 .02],'Callback',@filterData);

%button group to select filter method
% filterGraph = uibuttongroup('Parent',gf,'SelectionChangedFcn',@filterLabels,'Units','normalized','Position',[0.45,0.055,0.13,.05]);
sizeFilter = uicontrol('Parent',gf,'Style','checkbox','String','size','Units','normalized','Position',[0.45 0.085 .05 .02]);
clustFilter = uicontrol('Parent',gf,'Style','checkbox','String','clusterStat','Units','normalized','Position',[0.45 0.06 .05 .02]);
timeFilter = uicontrol('Parent',gf,'Style','checkbox','String','time','Units','normalized','Position',[0.45 0.035 .05 .02]);
mzFilter = uicontrol('Parent',gf,'Style','checkbox','String','mz','Units','normalized','Position',[0.45 0.01 .05 .02]);

%text boxes to select max and min for filter method
sizeMin = uicontrol('Parent',gf,'Style','edit','Units','normalized','Position',[0.51 0.085 .04 .02]);
sizeMax = uicontrol('Parent',gf,'Style','edit','Units','normalized','Position',[0.6,.085,0.04,.02]);
clustMin = uicontrol('Parent',gf,'Style','edit','Units','normalized','Position',[0.51 0.06 .04 .02]);
clustMax = uicontrol('Parent',gf,'Style','edit','Units','normalized','Position',[0.6,.06,0.04,.02]);
timeMin = uicontrol('Parent',gf,'Style','edit','Units','normalized','Position',[0.51 0.035 .04 .02]);
timeMax = uicontrol('Parent',gf,'Style','edit','Units','normalized','Position',[0.6,.035,0.04,.02]);
mzMin = uicontrol('Parent',gf,'Style','edit','Units','normalized','Position',[0.51 0.01 .04 .02]);
mzMax = uicontrol('Parent',gf,'Style','edit','Units','normalized','Position',[0.6,.01,0.04,.02]);
%labels for max and min
MaxLabel = uicontrol('Parent',gf,'Style','text','String','max','Units','normalized','Position',[0.61,.107,0.02,.015]);
MinLabel = uicontrol('Parent',gf,'Style','text','String','min','Units','normalized','Position',[0.52,.107,0.02,.015]);
%values for max and min
sizeMinV = uicontrol('Parent',gf,'Style','text','String',num2str(sizeLimits(1)),'Units','normalized','Position',[0.552,.085,0.04,.02]);
sizeMaxV = uicontrol('Parent',gf,'Style','text','String',num2str(sizeLimits(2)),'Units','normalized','Position',[0.642,.085,0.04,.02]);
clustMinV = uicontrol('Parent',gf,'Style','text','String',num2str(clustRLimits(1)),'Units','normalized','Position',[0.552,.06,0.04,.02]);
clustMaxV = uicontrol('Parent',gf,'Style','text','String',num2str(clustRLimits(2)),'Units','normalized','Position',[0.642,.06,0.04,.02]);
timeMinV = uicontrol('Parent',gf,'Style','text','String',num2str(timeLimits(1)),'Units','normalized','Position',[0.552,.035,0.04,.02]);
timeMaxV = uicontrol('Parent',gf,'Style','text','String',num2str(timeLimits(2)),'Units','normalized','Position',[0.642,.035,0.04,.02]);
mzMinV = uicontrol('Parent',gf,'Style','text','Visible','off','String','','Units','normalized','Position',[0.552,.01,0.04,.02]);
mzMaxV = uicontrol('Parent',gf,'Style','text','Visible','off','String','','Units','normalized','Position',[0.642,.01,0.04,.02]);
%user input to select mz value to sort by
mzValue2 = uicontrol('Parent',gf,'Style','edit','Units','normalized','Position',[0.475,.01,0.03,.02],'Callback',@filterbymz);

%create controls to combine clusters 
combineCluster = uicontrol('Parent',gf,'Style','pushbutton','String','Combine Clusters','Callback',@combineClusters,'Units','normalized','Position',[0.71,0.11,0.09,.02]);
clusterSelect1 = uicontrol('Parent',gf,'Style','listbox','String',uniqueLabel,'Units','normalized','Position',[0.7,0.045,0.05,.06]);
clusterSelect2 = uicontrol('Parent',gf,'Style','listbox','String',uniqueLabel,'Units','normalized','Position',[0.755,0.045,0.05,.06]);
resetCluster = uicontrol('Parent',gf,'Style','pushbutton','String','Reset Original Clusters','Callback',@resetClusters,'Units','normalized','Position',[0.7,0.022,0.11,.02]);

%create controls to separate clusters
separateCluster = uicontrol('Parent',gf,'Style','pushbutton','String','Split Cluster','Callback',@separateClusters,'Units','normalized','Position',[0.815,0.11,0.05,.02]);
clusterSelect3 = uicontrol('Parent',gf,'Style','listbox','String',uniqueLabel,'Units','normalized','Position',[0.815,0.045,0.05,.06]);
separatePrompt = uicontrol('Parent',gf,'Style','text','Visible','off','String','Select where to split!','Units','normalized','Position',[0.815,.01,0.05,.03]);

%create control to recalculate cluster relation
recalcClusterRel = uicontrol('Parent',gf,'Style','pushbutton','String','Recalculate Cluster Stats','Callback',@recalcClust,'Units','normalized','Position',[0.7,0.001,0.11,.02]);

%create controls to outputPIDs
outputLabel = uniqueLabel;
outputLabel{length(uniqueLabel)+1} = 'ALL';
outputPID = uicontrol('Parent',gf,'Style','pushbutton','String','Output PIDs to Workspace','Callback',@outputPIDs,'Units','normalized','Position',[0.89,0.11,0.1,.02]);
clusterSelect4 = uicontrol('Parent',gf,'Style','listbox','String',outputLabel,'Units','normalized','Position',[0.89,0.045,0.05,.06]);
varName = uicontrol('Parent',gf,'Style','edit','Units','normalized','Position',[0.945,.055,0.05,.02]);
varLabel = uicontrol('Parent',gf,'Style','text','String','variable name','Units','normalized','Position',[0.945,.076,0.05,.02]);

gf.Visible = 'on';

%%
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

    %function to recalculate cluster relation
    function recalcClust(source,callbackdata)
        %note the user could choose to use another method for calculating
        %the particle to cluster relation statistic. The method here is
        %chosen over silhouette because it is much faster
        for i = 1:length(uniqueCluster)
            cidx = clustData == uniqueCluster(i);
            tmp = (averageClust(i,:)/norm(averageClust(i,:)))*normc(useMZdata(cidx,:)');
            clustRelation(cidx) = tmp;
        end
        set(clustPlot,'XData',clustRelation(useIDX));
    end

    %function to display all clusters
    function resetDisplays(source,callbackdata)
        displaySubset = 0;
        plotHelper;
    end
    
    %function to select what clusters to display
    function selectClusters(source,callbackdata)
        displaySubset = 1;
        plotHelper;
    end 
    
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
        assignin('base',nameVar,outPID);
    end

    %function to split a cluster into two
    function separateClusters(source,callbackdata)
        isChecked = clusterGroup.Value;
        newC = max(uniqueCluster)+1; %new cluster index
        axes(ex);
        separatePrompt.Visible = 'on';
        [~,y] = ginput(1); %get user input on where to make split
        separatePrompt.Visible = 'off';
        y = round(y);
        %get indexes for splitting
        tmp = clustData(useIDX) == uniqueCluster(clusterSelect3.Value);
        sum1 = sum(tmp);
        tmp(1:y) = 0;
        if any(tmp) && sum(tmp) ~= sum1 %make user input creates a split
            [~,indexColor] = intersect(colorIDs,colorLabel(end,:),'rows'); %new color for new cluster
            indexColor = indexColor+1;
            if indexColor > size(colorIDs,1)
                indexColor = 1;
            end
            colorLabel = [colorLabel; colorIDs(indexColor,:)];
            clustData(useIDX(tmp)) = newC; %set new cluster ids for each particle
            clusterColor(useIDX(tmp),:) = repmat(colorLabel(end,:),sum(tmp),1); %set new color for each particle
            %update display with new cluster
            uniqueCluster(length(uniqueCluster)+1) = newC;
            uniqueLabel{length(uniqueCluster)} = num2str(newC);
            outputLabel = uniqueLabel;
            outputLabel{length(uniqueLabel)+1} = 'ALL';
            averageClust(clusterSelect3.Value,:) = nanmean(useMZdata(clustData == uniqueCluster(clusterSelect3.Value),:),1);
            averageClust(length(uniqueCluster),:) = nanmean(useMZdata(clustData == newC,:),1);
            clusterSelect3.String = uniqueLabel;
            clusterSelect1.String = uniqueLabel;
            clusterSelect2.String = uniqueLabel;
            clusterSelect.String = uniqueLabel;
            clusterSelect5.String = uniqueLabel;
            clusterSelect4.String = outputLabel;
            clusterSelect5.Max = length(uniqueLabel);
            if clusterSelect.Value == clusterSelect3.Value
                changeCluster(clusterSelect);
            end
            plotHelper; %update plot
        end
    end

    %function to return to original input cluster 
    function resetClusters(source,callbackdata)
        %reset data
        uniqueCluster = origuniqueCluster;
        uniqueLabel = origuniqueLabel;
        colorLabel = origcolorLabel;
        clustData = origClustData;
        clustRelation = origClustRelation;
        %recalculate some values
        averageClust = zeros(length(uniqueCluster),size(useMZdata,2));
        for q = 1:length(uniqueCluster)
            tmp = clustData == uniqueCluster(q);
            clusterColor(tmp,:) = repmat(colorLabel(q,:),sum(tmp),1);
            averageClust(q,:) = nanmean(useMZdata(tmp,:),1);
        end
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
        clusterSelect1.String = uniqueLabel;
        clusterSelect2.String = uniqueLabel;
        clusterSelect.String = uniqueLabel;
        clusterSelect4.String = outputLabel;
        clusterSelect5.String = uniqueLabel;   
        clusterSelect5.Max = length(uniqueLabel);
        changeCluster(clusterSelect);
        plotHelper; %update plot
    end

    %create cell array of strings from number vector
    function uniqueStr = makeClusterLabel(uniqueVector)
        if size(uniqueVector,1) > size(uniqueVector,2)
            uniqueStr = strtrim(cellstr(num2str(uniqueVector))');
        else
            uniqueStr = strtrim(cellstr(num2str(uniqueVector'))');
        end
    end

    %function to combine two clusters
    function combineClusters(source,callbackdata)
        %get clusters to combine
        cluster1 = clusterSelect1.Value;
        cluster2 = clusterSelect2.Value;
        if cluster1 ~= cluster2 %make sure different clusters are selected
            %min cluster is the label always kept
            minCIDX = min([cluster1 cluster2]);
            maxCIDX = max([cluster1 cluster2]);
            maxLabel = uniqueCluster(maxCIDX);
            minLabel = uniqueCluster(minCIDX);
            tmp = clustData == maxLabel; %find clusters to rename
            clustData(tmp) = minLabel; %rename to min cluster
            clusterColor(tmp,:) = repmat(colorLabel(minCIDX,:),sum(tmp),1); %fix color of particles to match min cluster
            averageClust(minCIDX,:) = nanmean(useMZdata(clustData == minLabel,:),1); %recalc average cluster spectrum
            %remove old cluster data
            averageClust(maxCIDX,:) = [];
            colorLabel(maxCIDX,:) = [];
            uniqueCluster(maxCIDX) = [];
            %set all selections to min cluster in case previous selections are no longer
            %viable after resetting
            clusterSelect2.Value = minCIDX;
            clusterSelect1.Value = minCIDX;
            clusterSelect.Value = minCIDX;
            clusterSelect3.Value = minCIDX;
            clusterSelect4.Value = minCIDX;
            clusterSelect5.Value = 1;
            changeCluster(clusterSelect);
            plotHelper; %update plots
            %update display
            uniqueLabel(maxCIDX) = [];
            outputLabel(maxCIDX) = [];
            clusterSelect3.String = uniqueLabel;
            clusterSelect1.String = uniqueLabel;
            clusterSelect2.String = uniqueLabel;
            clusterSelect.String = uniqueLabel;
            clusterSelect4.String = outputLabel;
            clusterSelect5.String = uniqueLabel;
            clusterSelect5.Max = length(uniqueLabel);
        end
    end

    %function to change the average cluster spectrum displayed
    function changeCluster(source,callbackdata)
        axes(ex);
        newClusterIDX = source.Value; %get new cluster to display
        ClustAveragePlot.YData = averageClust(newClusterIDX,:);
        ClustAveragePlot.Color = colorLabel(newClusterIDX,:);
    end

    %wrapper function to call plotHelper from callback
    function groupByCluster(source,callbackdata)
        plotHelper;
    end

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
        if strcmp(currentAxis,'log')
            spectraMap.CData = logData(useIDX,:);
        else
            spectraMap.CData = useMZdata(useIDX,:);
        end
        spectraMap.YData = [1 lengthTMP];
        set(sizePlot,'XData',sizeData(useIDX),'YData',1:lengthTMP);
        set(clustPlot,'XData',clustRelation(useIDX),'YData',1:lengthTMP,'CData',clusterColor(useIDX,:));
        set(timePlot,'XData',timeData(useIDX),'YData',1:lengthTMP);
        bx.YLim = [1,lengthTMP]; %update y limits
    end

    %function to sort data based upon selection
    function sortBigGraph(source,callbackdata)
        currentButton = sortGraph.SelectedObject.String; %how to sort
        %do sorts
        if strcmp(currentButton,'mz')
            sortValue = mzValue.String;
            sortValue = str2double(sortValue);
            sortValue = find(MZvector == sortValue);
            [~, tmpIDX] = sortrows(useMZdata,sortValue);
        elseif strcmp(currentButton,'time');
            [~,tmpIDX] = sort(timeData);
        elseif strcmp(currentButton,'size');
            [~,tmpIDX] = sort(sizeData);
        elseif strcmp(currentButton,'clusterStat');
            [~,tmpIDX] = sort(clustRelation);
        else
            tmpIDX = [1:size(useMZdata,1)]';
        end
        plotHelper %update plots
    end

    %function to switch between log and linear color axis for spectra heatmap
    %plot
    function replotBigGraph(source,callbackdata)
        currentAxis = source.SelectedObject.String; %determine if log or linear selected
        if strcmp(currentAxis,'log') %get log data if it hasn't been calculated yet
            if isempty(logData)
                logData = useMZdata/(minVal2);
                logData = log10(logData);
            end
            spectraMap.CData = logData(useIDX,:); %update plot with log data
            caxis(ax,[0 logData(maxIdx)]);
        else
            spectraMap.CData = useMZdata(useIDX,:); %updata plot with linear data
            caxis(ax,[0 maxVal]);
        end
    end
end
