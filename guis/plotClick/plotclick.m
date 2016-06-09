function [selectPIDs, trackX, trackY] = plotclick(usefig)
% SETUP_PLOTCLICK MUST BE RUN BEFORE RUNNING PLOTCLICK
% consider running clearPlotclick to clear global variables created during
% setup_plotclick after analysis with plotclick is complete. This will
% free up memory.
% PLOTCLICK takes a gscatter figure of particles in yaada and allows the
% user to select N (<100) regions of particles and automatically makes a number of
% plots (spectra, time, size, etc) pertaining to the selected particles.  
% To exit function press any key on keyboard twice while usefig is top window 
%
% Call as [selectPIDs, trackX, trackY] = plotclick(usefig)
% usefig = a gscatter plot made in setup_plotclick    
% selectPIDs = a 1xM cell array with each cell containing a Nx2 matrix of all
% the particle identifiers within each selection box made by the user.
% trackX: a Mx2 matrix noting the location of the x coordinates for the
% lower and upper corners of each selection box made by the user
% trackY: a Mx2 matrix noting the location of the y coordinates for the
% lower and upper corners of each selection box made by the user

global scatterLogSPEC plotTrack rownum scatterSPEC scatterMZ XVectorUse YVectorUse labels2use trackintersectPART SCATclr origPID scatterTime scatterTimeBin  scatterSIZE scatterSIZEbin scatterAllSize ;

uniqueLabels = unique(labels2use); %get cluster ids 

redraw = input('Do you want to redraw scatter plot (clearing all annotations and boxes)? 0-no, 1-yes  ');

%ask user how want to save plots
plotSave = input('\nWould you like to save plots?  0-no, 1-yes \nNote saving plots slows down the GUI. Consider only doing when making final selections.   ');
if plotSave == 1
    plotFolder = input('Full folder path   ','s');
    plotChoose = input('Save all plots automatically (0) or select plots to save (1)?   ');
    exportToPPTX('close'); %This will through a warning if tmp PPTX currently open, but that's ok.  Otherwise program crashes on next line if there is a current tmp PPTX open.
    exportToPPTX('new'); %initialize powerpoint to aggregate plots
end

%redraw scatter
if redraw
    figure(usefig)
    clf('reset')
    gs = gscatter(XVectorUse, YVectorUse,labels2use,[],[],[],'on');  % sym color?
end

%% start interacting with figure and creating plots
%set parameters to interact with figure
cont     = true; %keep interacting
clickcnt = 0 %tracks number of selections, can make up to 100 selections
trackX = []; %X coordinate output trackers
trackY = []; %Y coordinate output trackers

while cont && clickcnt<100
    
   %interact with gscatter plot 
   fprintf('\nplot click processing..%i',clickcnt)
   drawnow %have to run this cause of a bug in matlab, otherwise doesn't use user input 
   figure(usefig)
   hold on %to draw rectangles later
   [Xall, Yall, buttonAll] = ginput(2); %get two click points creating a selected rectangle of points
   if buttonAll(1) > 3 && buttonAll(2) > 3 %press any key twice to exit
       if plotSave == 1 %save ppt before exiting
           exportToPPTX('saveandclose',fullfile(plotFolder,'allplots'));
       end
       clearG = input('\nDo you want to clear global variables initiated in setup_plotclick?  0-no, 1-yes \nThis will free up memory, but you will have to rerun setup_plotclick if you want to work with this scatter plot again.  ');
       if clearG == 1
           clear global scatterLogSPEC plotTrack rownum scatterSPEC scatterMZ XVectorUse YVectorUse labels2use trackintersectPART SCATclr origPID scatterTime scatterTimeBin  scatterSIZE scatterSIZEbin scatterAllSize
       end
       return
   end
   clickcnt = clickcnt + 1; %increment click count, also becomes selection #
   
   %draw rectangle to show selected data using input clicks from user
   minX = min(Xall); 
   minY = min(Yall);
   maxX = max(Xall);
   maxY = max(Yall);
   rectangle('Position',[minX minY maxX-minX maxY-minY])
   text(maxX, maxY, sprintf('%i',clickcnt));
   
   trackX = [trackX; minX maxX]; %add values to output trackers
   trackY = [trackY; minY maxY];
   
   %find all data points in selection
   findX = XVectorUse > minX & XVectorUse < maxX; 
   findY = YVectorUse > minY & YVectorUse < maxY;
   findInt = findX == 1 & findY == 1; 

   %% if any data points in selection
   if any(findInt)
       % add selection to PID list
       selectPIDs{clickcnt} = origPID(findInt,:);
       
       % get index of particles and labels selected
       indexInt = find(findInt); %get index
       tmplabels = labels2use(findInt); %get labels for particles
       [uniqueLabelsTMP, uniqueIDX] = unique(tmplabels); %get list of clusters selected 
       uniqueIDX = [uniqueIDX; length(tmplabels)+1]; %find location of each cluster in list
       numCluster = uniqueIDX(2:end) - uniqueIDX(1:(end-1)); %get number of particles in each cluster out of selection
       fracCluster = numCluster/(sum(numCluster)); %get fraction of particles in each cluster out of selection
       [~,LabelIdx,LabelIdxTMP] = intersect(uniqueLabels,uniqueLabelsTMP);

       if plotTrack(9) == 0
           multiPlotClick; %plotting script for multiple windows
       elseif plotTrack(9) == 1
           singlePlotClick %plotting script for single window
       end
       
   else disp('no row in useVector found for cursor');
       selectPIDs{clickcnt} = []; %no PIDS to add
   end
end; 
