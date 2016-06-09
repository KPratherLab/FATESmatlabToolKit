function groupedClusters = clickDendo(dendoFig)
% setup_clickDendo MUST BE RUN BEFORE RUNNING clickDendo
% consider running clearDendo to clear global variables created during
% setup_clickDendo after analysis with clickDendo is complete. This will
% free up memory.
% Call as groupedClusters = clickDendo(dendoFig)
% clickDendo allows the user to interact with a dendogram of clusters,
% which have been linked by the average cluster spectra
% Branches or leaves of the dendrogram can be selected to produce plots of particle/cluster
% spectra and temporals and size distributions 
% To exit function while interacting with figure press any key on keyboard twice while dendoFig is top window 
%
% When interacting with the plot the user needs to click twice to select
% dendogram branches.  The location of the first click determines what branch is
% selected.  The location of the second click is not important but
% determines what type of plots are created according to the rules below...
% Left click, Left click: creates plots showing all clusters of selected
% branch
% Left click, Right click: creates plots separating data into only the left
% and right branch
% Right click, Right click: no plots are created, but all clusters
% belonging to the selected branch are sent to the groupedClusters output,
% and the branches are highlighted on the plot
%
% dendoFig = a dendogram made in setup_clickDendo
% groupedClusters = a 1xN cell array containing the list of clusters 
% belonging to each highlighted group of branches created by the user
% right clicking twice on the dendogram

global dendoSPEC dendoLogSPEC dendoSPECavg dendoSIM dendoMZ dendoTime dendoSize dendoTimeBin dendoSizeBin linkCluster dendoCOUNTS Xlines Ylines finalbaseCluster lineCluster plotTrack rownum

%% begin interacting with dendogram
%get user input
cont = 1;
clickcnt = 0;
clustsave = 0;
tmp = colormap(jet);
usecolors = tmp(4:10:end,:);
redraw = input('Do you want to redraw dendrogram (clearing all annotations and colors)? 0-no, 1-yes  ');
savePlots = input('Do you want to save all plots created? 0-no, 1-yes \nNote saving plots slows down the GUI. ');
if savePlots
    plotFolder = input('What is the full path to save all plots to?  ','s');
    exportToPPTX('close'); %This will through a warning if tmp PPTX currently open, but that's ok.  Otherwise program crashes on next line if there is a current tmp PPTX open.
    exportToPPTX('new'); %initialize powerpoint to aggregate plots
end

if redraw
    figure(dendoFig)
    [Dhandle] = dendrogram(linkCluster,0);
    set(Dhandle,'Color',[0 0 0]);
end

while cont
    drawnow %have to run this cause of a bug in matlab, otherwise doesn't use user input
    cont = input('Get next selection (1), or stop interacting with figure (0), select single cluster to plot (2)?  '); %this pause allows the user to zoom in and out of the figure before selecting nodes
    if cont == 0
        exportToPPTX('saveandclose',fullfile(plotFolder,'allplots')); %save plots
        %attempt to clear out some memory if won't be working with
        %dendogram again shortly.  Whether the user wants to do this or not
        %will depend on how much memory is available and how long
        %setup_clickHist took
        clearG = input('\nDo you want to clear global variables initiated in setup_clickHist?  0-no, 1-yes \nThis will free up memory, but you will have to rerun setup_clickHist if you want to work with this dendogram again.  ');
        if clearG == 1
            clear global dendoSPEC dendoLogSPEC dendoSPECavg dendoSIM dendoMZ dendoTime dendoSize dendoTimeBin dendoSizeBin linkCluster dendoCOUNTS Xlines Ylines finalbaseCluster lineCluster
        end
        break
    elseif cont == 1
    figure(dendoFig)
    [Xall, Yall, buttonAll] = ginput(2); %the first click point sets the location, the second click location is not important but is needed
    if buttonAll(1) > 3 | buttonAll(2) > 3 %press any key twice to exit interaction
        %attempt to clear out some memory if won't be working with
        %dendogram again shortly.  Whether the user wants to do this or not
        %will depend on how much memory is available and how long
        %setup_clickHist took
        clearG = input('\nDo you want to clear global variables initiated in setup_clickHist?  0-no, 1-yes \nThis will free up memory, but you will have to rerun setup_clickHist if you want to work with this dendogram again.  ');
        if clearG == 1
            clear dendoSPEC dendoLogSPEC dendoSPECavg dendoSIM dendoMZ dendoTime dendoSize dendoTimeBin dendoSizeBin linkCluster dendoCOUNTS Xlines Ylines finalbaseCluster lineCluster plotTrack rownum
        end
        cont = 0;
    end
    checkYCoord = linkCluster(:,3) < Yall(1)+0.01 & linkCluster(:,3) > Yall(1)-0.01; %check for nearby node in y direction
    if any(checkYCoord) 
        checkXCoord = Xlines(:,1) < Xall(1) & Xlines(:,3) > Xall(1) & checkYCoord; %check for nearby node in x direction
        if any(checkXCoord) %nearby node found            
            clickcnt = clickcnt+1; %track clicks to set color for highlighter line
            clrcnt = clickcnt;
            if clrcnt > size(usecolors,1); %reset color cycle
                clrcnt = 1;
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
            
            %if two left clicks then plots will highlight all clusters in
            %the node
            if buttonAll(1) == 1 && buttonAll(2) == 1
                line(findXLine, findYLine, 'Color', usecolors(clrcnt,:)); %highlight branch investigating
                %find clusters belonging to node
                leftIDs = finalbaseCluster{findRow}{1};
                rightIDs = finalbaseCluster{findRow}{2};
                allIDs = [leftIDs rightIDs];
                IDstr = [leftIDs 0 rightIDs];
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
                useCtmpL = {g(1:length(leftIDs)).Color};
                useCtmpR = {g((length(leftIDs)+1):end).Color};
                close(a);
                %MAKE PLOTS
                singlePlotClick_dendo


            %if left click followed by right click then plots will show the
            %left versus right branch of the node
            elseif buttonAll(1) == 1 && buttonAll(2) == 3
                %find location of line to highlight
                leftX = [findXLine(1:2) findXLine(2)+(findXLine(3)-findXLine(2))/2];
                leftY = findYLine(1:3);
                rightX = [findXLine(2)+(findXLine(3)-findXLine(2))/2 findXLine(3:4)];
                rightY = findYLine(2:end);  
                line(leftX, leftY, 'Color', 'red'); %highlight left branch as red
                line(rightX, rightY, 'Color', 'blue'); %highlight right branch as blue
                %find clusters belonging to the left and right branch
                leftIDs = finalbaseCluster{findRow}{1};
                rightIDs = finalbaseCluster{findRow}{2};
                %get stats on clusters selected
                leftCounts = dendoCOUNTS(leftIDs);
                rightCounts = dendoCOUNTS(rightIDs);
                leftFraction = sum(leftCounts)/(sum(leftCounts)+sum(rightCounts)); %relative contribution of left branch to selection
                rightFraction = sum(rightCounts)/(sum(leftCounts)+sum(rightCounts)); %relative contribution of right branch ot selection
                %MAKE PLOTS
                singlePlotClick_dendo_leftRight      
                
            %if two right clicks then no plots are made but clusters are saved to output and all selected branches are highlighted    
            elseif (buttonAll(1) == 3 && buttonAll(2) == 3) || (buttonAll(1) == 3 && buttonAll(2) == 1) 
                clustsave = clustsave + 1;
                allXcoord = Xlines(lineCluster{findRow},:);
                allYcoord = Ylines(lineCluster{findRow},:);
                hold on
                plot(allXcoord', allYcoord', 'Color', usecolors(clrcnt,:),'LineWidth',2);
                text(findXLine(2)+(findXLine(3)-findXLine(2))/2, findYLine(2), sprintf('%i',clustsave));
                groupedClusters{clustsave} = [finalbaseCluster{findRow}{1} finalbaseCluster{findRow}{2}];
                hold off
            end
        else
            disp('No nearby nodes.  Please click again');
        end
    else
        disp('No nearby nodes.  Please click again');
    end
    elseif cont == 2
        clustNUM = input('\nWhat number cluster do you want to plot?  ');
        if clustNUM >0 && clustNUM < length(dendoCOUNTS)
            clickcnt = clickcnt+1; %track clicks
            singlePlotCluster_dendo
        else
        disp('\n Not an available cluster number');
        end          
    end        
end