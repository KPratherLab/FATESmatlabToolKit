function clickDendo_GUI(outWM,clustLabels,clustData,useMZdata,MZvector,partData,timeData,sizeData,clustRelation)
% clickDendo_GUI sets up a dendogram using relations calculated from input
% spectra.  The user can then interact with the GUI and select
% nodes/branches to examine the particles in each cluster in GUIfates or output the particle ids to the
% workspace. 
% Call as clickDendo_GUI(outWM,clustLabels,clustData,useMZdata,MZvector,partData,timeData,sizeData,clustRelation)
% 
% INPUT
% outWM= a JxK matrix containing K spectra. The spectra are used to group 
% by hierachical clustering (this would usually be the output outWM from art2a or 
% some other clustering algorithm (k-means))
% and then create the dendogram.  
% clustLabels = A vector of length K with each entry corresponding to the
% cluster label to be used for the corresponding entry in outWM. If
% clustLabels is empty then the vector labels will just be 1:K.
% clustData: a vector of length M with the cluster/group indetifier for the
% particles in partData. The cluster identifiers in clustData should correspond to the 
% clusters listed in in clustLabels to be able to explore the data in GUIfates.
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
% clustRelation: a vector of length M with the cluster relation statistic
% for the particles in partData.  This cluster relation statistic is
% computed by the user and could be supplied by a a variety of
% calculations.  It is meant to show the relation of the particles to the 
% cluster it is assigned to. One option is to use the output from silhouette (built-in
% matlab function) or to take the dot product of each particle spectra with the
% average spectra of the cluster to which it belongs.
%
% Output: "Output all PIDs to workspace"
% If the user selects this then all particle identifiers for all clusters
% belonging to the node selected will be output to the workspace as well as
% the cluster labels for each cluster.  The name of the variable is
% specified in the text box. 

%% check inputs
numClust = size(outWM,2);
numPart = size(useMZdata,1);

if nargin < 5
   error('Too few input arguments.');
elseif nargin > 9  
    error('Too many input arguments.');
end

if isempty(clustLabels)
    clustLabels = 1:numClust;
elseif length(clustLabels) ~= numClust
        error('The length of clustLabels must equal the number of columns of outWM. This means for each cluster label there must be a corresponding spectra');
end

if isempty(clustData) || ~exist('clustData','var')
    error('Must provide clustData to specify the cluster that each particle belongs to');
else
    if size(clustData,1) < size(clustData,2)
        clustData = clustData';
    end
    if size(clustData,1) ~= numPart
        error('clustData vector needs length needs same number of rows as spectra data (same # of particles)');
    end
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

if isempty(useMZdata) || ~exist('timeData','var')
    error('must input matrix of spectra values to plot in guifates');
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
    clustRelation = zeros(numPart,2);
else
    if size(clustRelation,1) < size(clustRelation,2)
        clustRelation = clustRelation';
    end
    if size(clustRelation,1) ~= numPart
        error('clustRelation vector needs length needs same number of rows as spectra data (same # of particles)');
    end
end

%% setup GUI figure

% do hierarchical clustering analysis
distCluster = pdist(outWM');
linkCluster = linkage(distCluster);

%initialize figure window
dendoFig = figure('Visible','off');
set(dendoFig,'units','normalized');
set(dendoFig,'position',[.1 .1 .5 .65]);

% create axes at right to display selection stats
statx = axes('Parent',dendoFig,'Position', [0.03 0.14 0.025 0.82]);
statx.XTick =[];
statx.YDir = 'reverse';

% create initial histogram plot
dFx = axes('Parent',dendoFig,'Units','normalized','Position',[0.15 0.14 0.82 0.82]);
[Dhandle, IDXdendo, XLabel] = dendrogram(linkCluster,0);
set(Dhandle,'Color',[0 0 0]);
dFx.YLabel = [];

% fix x tick label 
dendoXlabel = clustLabels(XLabel); %reorder cluster labels in the order they appear in the dendogram
dFx.XTickLabel = strtrim(cellstr(num2str(dendoXlabel'))');
hold on

% set up uicontrols
drawSelect = uicontrol('Parent',dendoFig,'Style','pushbutton','String','Select Node for Inspection','Callback',@nodeGUIfates,'Units','normalized','Position',[0.1,0.07,0.4,.03]); %pushbotton to create selection
%button group to select how to display clusters
displayType = uibuttongroup('Parent',dendoFig,'Units','normalized','Position',[0.1,0.01,0.4,.05]);
allclust = uicontrol(displayType,'Style','radiobutton','Units','normalized','Position',[0.05 0.05 .03 .9]);
leftright = uicontrol(displayType,'Style','radiobutton','Units','normalized','Position',[0.35 0.05 .03 .9]);
output = uicontrol(displayType,'Style','radiobutton','Units','normalized','Position',[0.7 0.03 .05 .9]);
%labels for radiobuttons
allclustLabel = uicontrol('Parent',dendoFig,'Style','text','String','display as clusters','Units','normalized','Position',[0.15,.015,0.06,.04]);
leftrightLabel = uicontrol('Parent',dendoFig,'Style','text','String','display as left right branch','Units','normalized','Position',[0.27,.015,0.08,.04]);
outputLabel = uicontrol('Parent',dendoFig,'Style','text','String','output all PIDs to workspace','Units','normalized','Position',[0.4,.015,0.09,.04]);
% text box to put variable name for output
varName = uicontrol('Parent',dendoFig,'Style','edit','Units','normalized','Position',[0.51,.02,0.15,.02]); %name for output varialbe
varLabel = uicontrol('Parent',dendoFig,'Style','text','String','output variable name','Units','normalized','Position',[0.51,.05,0.15,.02]); %label for varName

%pushbotton to clear all drawn lines
clearSelect = uicontrol('Parent',dendoFig,'Style','pushbutton','String','Clear Lines','Callback',@clearLines,'Units','normalized','Position',[0.75,0.02,0.2,.03]); %pushbutton to clear all selections

%text to show % of particles selected
partPercent = uicontrol('Parent',dendoFig,'Style','text','Visible','off','Units','normalized','Position',[0.01,.05,0.05,.03]); %label for varName

%% get data ready for plotting
%get # of particles in clusters
numCluster = size(outWM,2);
dendoCOUNTS = zeros(1,length(clustLabels));
for i = 1:length(clustLabels)
    dendoCOUNTS(i) = sum(clustData == clustLabels(i));
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

%set up variables
clickcnt = 0;
trackcnt = 0;
dendoFig.Visible = 'on';
tmp = colormap(jet);
lineholder = gobjects(0);
usecolors = tmp(4:10:end,:);

    %function to modify display with selection and then initiate GUIfates
    function nodeGUIfates(source,callbackdata)        
        %begin interacting with dendogram
        %get user input
        trackcnt = trackcnt+1;
        [Xall, Yall] = ginput(1); %the first click point sets the location, the second click location is not important but is needed
        checkYCoord = linkCluster(:,3) < Yall(1)+0.01 & linkCluster(:,3) > Yall(1)-0.01; %check for nearby node in y direction
        
        %check for nearby nodes
        if any(checkYCoord)
            checkXCoord = Xlines(:,1) < Xall(1) & Xlines(:,3) > Xall(1) & checkYCoord; %check for nearby node in x direction
            if any(checkXCoord) %nearby node found
                clickcnt = clickcnt+1; %track clicks to set color for highlighter line
                if clickcnt > size(usecolors,1); %reset color cycle
                    clickcnt = 1;
                end
                findRow = find(checkXCoord);
                findXLine = Xlines(findRow,:); %find nearby nodes
                findYLine = Ylines(findRow,:);
                if length(findRow) > 1 %if more than one nearby node select closest in y direction
                    distY = abs(Ylines(findRow,2) - Yall(1)); %calc distance
                    [~, idxmin] = min(distY); %find nearest nodes
                    findRow = findRow(idxmin);
                    findXLine = Xlines(findRow,:); %find nearby nodes
                    findYLine = Ylines(findRow,:);
                end
                
                %if all clusters selected
                if displayType.Children(3).Value == 1
                    
                    %modify display
                    lineholder(trackcnt) = line(findXLine, findYLine, 'Color', usecolors(clickcnt,:),'Parent',dFx); %highlight branch investigating
                    delete(findobj(statx, 'Type','rectangle')); %get rid of population info on the left
                    delete(findobj(statx, 'Type','text')); %get rid of population info on the left
                    
                    %find clusters belonging to node
                    leftIDs = finalbaseCluster{findRow}{1};
                    rightIDs = finalbaseCluster{findRow}{2};
                    allIDs = [leftIDs rightIDs];
                    
                    %get stats on clusters selected
                    leftCounts = dendoCOUNTS(leftIDs);
                    rightCounts = dendoCOUNTS(rightIDs);
                    allCounts = dendoCOUNTS(allIDs);
                    leftFraction = sum(leftCounts)/(sum(leftCounts)+sum(rightCounts));
                    rightFraction = sum(rightCounts)/(sum(leftCounts)+sum(rightCounts));
                    fracCluster = allCounts/(sum(allCounts)); %relative contribution of each cluster to selection
                    
                    %create tmp figure to figure out how to color later plots
                    a = figure('Visible','off');
                    x = [0 1];
                    yL = ones(length(allIDs),2);
                    g = plot(x,yL);
                    useCtmp = {g.Color};
                    close(a);
                    useCtmp = useCtmp';
                    useCtmp = cell2mat(useCtmp);
                    
                    %create legend at right showing percentage of each branch and containing clusters
                    for q = 1:length(allIDs) % percentage of each cluster
                        clr = useCtmp(q,:);
                        rectangle('Parent',statx,'Position',[0,sum(fracCluster(1:q))-fracCluster(q),1,fracCluster(q)],'FaceColor',clr,'EdgeColor',clr);
                        clabel = num2str(clustLabels(allIDs(q)));
                        text(1.1,sum(fracCluster(1:q))-fracCluster(q)/2,clabel,'Parent',statx)
                    end
                    
                    %create boxes framing all left and all right clusters
                    rectangle('Parent',statx,'Position',[0,0,1,leftFraction],'EdgeColor',[0 0 0],'LineWidth',2);
                    rectangle('Parent',statx,'Position',[0,leftFraction,1,rightFraction],'EdgeColor',[0 0 0],'LineWidth',2);
                    statx.YTick = [leftFraction/2, leftFraction+rightFraction/2];
                    statx.YTickLabel = {'left branch' 'right branch'};
                    statx.YTickLabelRotation = 90;
                    statx.YLim = [0 1];
                    
                    %fraction of particles selected out of all in dendogram
                    t = sprintf('%.2f %%',100*sum(allCounts)/sum(dendoCOUNTS));
                    partPercent.String = t;
                    partPercent.Visible = 'on';
                    
                    %get location of data for selected clusters
                    tmpClust = clustLabels(allIDs); %need to get clusters and also preserve the order as displayed in the dendogram
                    idxClust = false(size(useMZdata,1),1);
                    for m = 1:length(allIDs)
                        idxClust = idxClust | clustData == tmpClust(m);
                    end
                    %MAKE PLOT
                    GUIfates(useMZdata(idxClust,:),MZvector,partData(idxClust,:),timeData(idxClust),sizeData(idxClust),clustData(idxClust),clustRelation(idxClust),useCtmp,tmpClust);
                    
                %if left right branch selected
                elseif displayType.Children(2).Value == 1
                    
                    delete(findobj(statx, 'Type','rectangle')); %get rid of population info on the left
                    delete(findobj(statx, 'Type','text')); %get rid of population info on the left
                    
                    %find location of line to highlight
                    leftX = [findXLine(1:2) findXLine(2)+(findXLine(3)-findXLine(2))/2];
                    leftY = findYLine(1:3);
                    rightX = [findXLine(2)+(findXLine(3)-findXLine(2))/2 findXLine(3:4)];
                    rightY = findYLine(2:end);
                    lineholder(trackcnt) = line(leftX, leftY, 'Color', 'red','Parent',dFx); %highlight left branch as red
                    trackcnt = trackcnt + 1;
                    lineholder(trackcnt) = line(rightX, rightY, 'Color', 'blue','Parent',dFx); %highlight right branch as blue
                    
                    %find clusters belonging to the left and right branch
                    leftIDs = finalbaseCluster{findRow}{1};
                    rightIDs = finalbaseCluster{findRow}{2};
                    
                    %get stats on clusters selected
                    leftCounts = dendoCOUNTS(leftIDs);
                    rightCounts = dendoCOUNTS(rightIDs);
                    leftFraction = sum(leftCounts)/(sum(leftCounts)+sum(rightCounts)); %relative contribution of left branch to selection
                    rightFraction = sum(rightCounts)/(sum(leftCounts)+sum(rightCounts)); %relative contribution of right branch ot selection
                    
                    %create boxes framing all left and all right clusters
                    rectangle('Parent',statx,'Position',[0,0,1,leftFraction],'FaceColor','red','EdgeColor','red');
                    rectangle('Parent',statx,'Position',[0,leftFraction,1,rightFraction],'FaceColor','blue','EdgeColor','blue');
                    statx.YTick = [leftFraction/2, leftFraction+rightFraction/2];
                    statx.YTickLabel = {'left branch' 'right branch'};
                    statx.YTickLabelRotation = 90;
                    statx.YLim = [0 1];
                    
                    %fraction of particles selected out of all in dendogram
                    t = sprintf('%.2f %%',100*(sum(leftCounts)+sum(rightCounts))/sum(dendoCOUNTS));
                    partPercent.String = t;
                    partPercent.Visible = 'on';    
                    
                    %get location of data for selected clusters
                    tmpClustData = clustData;
                    %left ID location
                    tmpClust = clustLabels(leftIDs); %need to get clusters and also preserve the order as displayed in the dendogram
                    LidxClust = false(size(useMZdata,1),1);
                    for m = 1:length(tmpClust)
                        LidxClust = LidxClust | clustData == tmpClust(m);
                    end
                    %right ID location
                    tmpClust = clustLabels(rightIDs); %need to get clusters and also preserve the order as displayed in the dendogram
                    RidxClust = false(size(useMZdata,1),1);
                    for m = 1:length(tmpClust)
                        RidxClust = RidxClust | clustData == tmpClust(m);
                    end
                    idxClust = RidxClust | LidxClust; %idx for all data selected
                    
                    %alter clust data to combine clusters into single
                    %left/right groups
                    tmpClustData(LidxClust) = 1;
                    tmpClustData(RidxClust) = 2;           
                    
                    useCtmp = [1 0 0; 0 0 1]; %use red and blue to color points is GUIfates for consistency
                    
                    %MAKE PLOT
                    GUIfates(useMZdata(idxClust,:),MZvector,partData(idxClust,:),timeData(idxClust),sizeData(idxClust),tmpClustData(idxClust),clustRelation(idxClust),useCtmp);
                 
                    
                 %if just want to output PIDS and highlight branch
                else                  
                    delete(findobj(statx, 'Type','rectangle')); %get rid of population info on the left
                    delete(findobj(statx, 'Type','text')); %get rid of population info on the left
                    
                    %find clusters belonging to node
                    leftIDs = finalbaseCluster{findRow}{1};
                    rightIDs = finalbaseCluster{findRow}{2};
                    allIDs = [leftIDs rightIDs];
                    
                    %fraction of particles selected out of all in dendogram
                    allCounts = dendoCOUNTS(allIDs);
                    t = sprintf('%.2f %%',100*sum(allCounts)/sum(dendoCOUNTS));
                    partPercent.String = t;
                    partPercent.Visible = 'on';
                    
                    %get location of data for selected clusters
                    tmpClust = clustLabels(allIDs); %need to get clusters and also preserve the order as displayed in the dendogram
                    idxClust = false(size(useMZdata,1),1);
                    for m = 1:length(allIDs)
                        idxClust = idxClust | clustData == tmpClust(m);
                    end
                    
                    %output data
                    nameVar = varName.String;
                    assignin('base',nameVar,{partData(idxClust,:) clustData(idxClust)});
                    
                    %highlight selection on plot
                    allXcoord = Xlines(lineCluster{findRow},:);
                    allYcoord = Ylines(lineCluster{findRow},:);
                    lineholder(trackcnt:(trackcnt+size(allXcoord,1)-1)) = plot(allXcoord', allYcoord', 'Color', usecolors(clickcnt,:),'LineWidth',2,'Parent',dFx);
                end
            else
                disp('No nearby nodes.  Please try again');
            end
        else
            disp('No nearby nodes.  Please try again');
        end
    end

    function clearLines(source,callbackdata)
        delete(findobj(statx, 'Type','rectangle')); %get rid of population info on the left
        delete(findobj(statx, 'Type','text')); %get rid of population info on the left
        partPercent.Visible = 'off';
        delete(lineholder);
        lineholder = gobjects(0);
        trackcnt = 0;
    end

end