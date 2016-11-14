function scatterFATES(XVector,YVector,CVector,useMZdata,MZvector,partData,timeData,sizeData,clustRelation) 
% scatterFATES sets up a gscatter figure of particles in FATES and allows the
% user to select regions of particles and then inspect those regions of particles
% using GUIfates which automatically makes a number of
% plots (spectra, time, size, etc) pertaining to the selected particles.  
% To select particles the user selects two points on the gscatter plot
% which is used to create a rectangle, all paritcles within the rectangle
% fall within the selection.
% Call as scatterFATES(XVector,YVector,CVector,useMZdata,MZvector,partData,timeData,sizeData,clustRelation) 
%
% INPUT
% XVector = metric pertaining to particles in partData (eg size of particle), a NX1 vector. These 
% values will be plotted on the x axis of the gscatter plot.  The nth row of XVector
% should correspond to the nth particle (row) of partData.  
% YVector = metric pertaining to particles in partData, a NX1 vector (eg laser power). These 
% values will be plotted on the Y axis of of the gscatter plot.  The nth row of YVector
% should correspond to the nth particle (row) of partData.
% CVector = a vector of length M with the cluster/group indetifier for the
% particles in partData. These values will be used to color the points in the gscatter plot.
% For example for a list of 10 particles that is
% made up of two groups clustData may look like [1 1 1 2 2 2 2 2 1 1].  In
% this case particles 1-3 and 9-10 belong to cluster 1, and particles 4-8
% belong to cluster 2.
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
% timeData: a vector of length M with time data (in numeric form) for the
% particles in partData
% sizeData: a vector of length M with size data for the particles in
% partData
% clustData: a vecotr of length M with the cluster/group indetifier for the
% particles in partData.  For example for a list of 10 particles that is
% made up of two groups clustData may look like [1 1 1 2 2 2 2 2 1 1].  In
% this case particles 1-3 and 9-10 belong to cluster 1, and particles 4-8
% belong to cluster 2.
% clustRelation:a vector of length M with the cluster relation statistic
% for the particles in partData.  This cluster relation statistic is
% computed by the user and could be supplied by a a variety of
% calculations.  It is meant to show the relation of the particles to the 
% cluster it is assigned to. One option is to use the output from silhouette (built-in
% matlab function) or to take the dot product of each particle spectra with the
% average spectra of the cluster to which it belongs.
% 
% OUTPUT
% "Output Selection to Workspace": If you click on this button it will
% output a 3xM cell array where M is either equal to 1 or the number of
% selections that currently exist.  The first row is the particle id of all
% the particles that fall within the selection. The second row is the
% xlimits for the rectangle that made up the selection, and the third row
% is the ylimits. The name of the variable is specified by the user in the text box 
% labeled variable name.

%% create figure

% initialize figure
usefig = figure('Visible','off');
set(usefig,'units','normalized');
set(usefig,'position',[.1 .1 .5 .65]);
ux = axes('Parent',usefig,'Units','normalized','Position',[0.11 0.14 0.85 0.82]);
gs      = gscatter(XVector,YVector,CVector,[],[],[],'on');  % sym color?
ux.YLabel = [];
ux.XLabel = [];
hold on

%add user controls
drawSelect = uicontrol('Parent',usefig,'Style','pushbutton','String','Select Particles','Callback',@drawBox,'Units','normalized','Position',[0.05,0.02,0.2,.03]); %pushbotton to create selection
clearSelect = uicontrol('Parent',usefig,'Style','pushbutton','String','Clear Selections','Callback',@clearBox,'Units','normalized','Position',[0.26,0.02,0.2,.03]); %pushbutton to clear all selections
doGUIfates = uicontrol('Parent',usefig,'Style','pushbutton','String','Inspect Selection','Callback',@startupGUIfates,'Units','normalized','Position',[0.55,0.08,0.10,.03]); %pushbutton to activate guifates
listSelect = uicontrol('Parent',usefig,'Style','listbox','Units','normalized','Position',[0.55,0.005,0.10,.07]); %list of selections created, feeds doGUIfates
outputSelect = uicontrol('Parent',usefig,'Style','pushbutton','String','Output Selection to Workspace','Callback',@outputData,'Units','normalized','Position',[0.75,0.08,0.2,.03]); %pushbutton to output selections
listSelect2 = uicontrol('Parent',usefig,'Style','listbox','Units','normalized','Position',[0.75,0.005,0.10,.07]); %list of selections created, feeds outputSelect
varName = uicontrol('Parent',usefig,'Style','edit','Units','normalized','Position',[0.86,.02,0.1,.02]); %name for output varialbe
varLabel = uicontrol('Parent',usefig,'Style','text','String','variable name','Units','normalized','Position',[0.86,.05,0.1,.02]); %label for varName

%get color for each group in plot (will feed into guifates to keep colors
%consistent)
SCATclr = {gs(:).Color};
SCATclr = SCATclr';
SCATclr = cell2mat(SCATclr);
%get order of labels in the scatter plot
SCATlabel = {gs(:).DisplayName};
SCATlabel = cellfun(@str2num, SCATlabel);
[~, labelIDX] = sort(SCATlabel); %want the same order as if the cluster ids were sorted
SCATclr = SCATclr(labelIDX,:); %reorder the colors 
uniqueLabels = unique(CVector); 

%setup variables
clickcount = 0;
trackX = []; %X coordinate output trackers
trackY = []; %Y coordinate output trackers
selectIndex = {};
selectionList = {};
outputList ={};

usefig.Visible = 'on';

%% callback functions

    %function to create selection boxes on the figure
    function drawBox(source,callbackdata)
        clickcount = clickcount + 1; %increment click count, also becomes selection #        
        %draw rectangle to show selected data using input clicks from user
        [Xall, Yall] = ginput(2); %get two click points creating a selected rectangle of points
        minX = min(Xall);
        minY = min(Yall);
        maxX = max(Xall);
        maxY = max(Yall);
        rectangle('Position',[minX minY maxX-minX maxY-minY])
        text(maxX, maxY, sprintf('%i',clickcount));
        
        %add values to output trackers
        trackX = [trackX; minX maxX]; %rectangle info
        trackY = [trackY; minY maxY];
        
        %find particles that fall in rectangle
        findX = XVector > trackX(clickcount,1) & XVector < trackX(clickcount,2);
        findY = YVector > trackY(clickcount,1) & YVector < trackY(clickcount,2);
        selectIndex{clickcount} = findX == 1 & findY == 1;
        
        %update display
        selectionList{clickcount} = sprintf('%i, %.2f%%',clickcount,100*(sum(selectIndex{clickcount})/length(XVector)));
        listSelect.String = selectionList;
        outputList = selectionList;
        outputList{1+clickcount} = 'ALL';
        listSelect2.String = outputList;
    end

    %function to clear selections (rectangles, labels, and lists)
    function clearBox(source,callbackdata)
        delete(findobj(usefig, 'Type','rectangle'));
        delete(findobj(usefig, 'Type','text'));
        clickcount = 0;
        trackX = [];
        trackY = [];
        selectionList = {};
        selectIndex = {};
        outputList = {};
        listSelect.String = selectionList;
        listSelect2.String = outputList;
        listSelect.Value = 1;
        listSelect2.Value = 1;
    end

    %wrapper function to feed into GUIfates
    function startupGUIfates(source,callbackdata)
        findSelection = listSelect.Value; %determine selection
        %find all data points in selection    
        getIDX  = selectIndex{findSelection};
        clustersSelected = unique(CVector(getIDX)); %find the clusters that were in the selection
        [~,colorIDX] = intersect(uniqueLabels, clustersSelected); %get the index of the clusters that were selected out of all clusters
        useColor = SCATclr(colorIDX,:); %get colors corresponding to the selected clusters
        GUIfates(useMZdata(getIDX,:),MZvector,partData(getIDX),timeData(getIDX),sizeData(getIDX),CVector(getIDX),clustRelation(getIDX),useColor);
    end

    %function to output particle lists and selection lists to workspace
    function outputData(source,callbackdata)
        findSelect = listSelect2.Value; %get cluster to output
        nameVar = varName.String; %get name for variable in workspace
        if findSelect <= clickcount %for single cluster selection
            outPID = cell(4,1);
            outPID{1,1} = partData(selectIndex{findSelect},:);
            outPID{2,1} = CVector(selectIndex{findSelect});
            outPID{3,1} = trackX(findSelect,:);
            outPID{4,1} = trackY(findSelect,:);
        else %to output all clusters
            outPID = cell(4,clickcount);
            for i = 1:clickcount
                outPID{1,i} = partData(selectIndex{i},:);
                outPID{2,i} = CVector(selectIndex{i});
                outPID{3,i} = trackX(i,:);
                outPID{4,i} = trackY(i,:);
            end
        end
        assignin('base',nameVar,outPID);
    end
end
