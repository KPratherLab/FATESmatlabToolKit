%% spectra plots
%if useSPEC is provided can set to 0 or 1, if no useSPEC then must set to 0
plotSpec = 1;
    
%if plotSpec is set to 1, can set to 0 or 1, otherwise all must be set to 0
%Plot all particle spectra overlaid? 0-no, 1-yes  
plotPS = 1;
%Plot average spectra of each cluster overlaid? 0-no, 1-yes
plotCS = 1;
%Create the heat map of all particle spectra? 0-no, 1-yes
plotHM = 1;
%Create the heat map of all particle spectra with log10 scaling? 0-no, 1-yes 
plotHMlog = 1;
%Plot the similarity of spectra to averages? 0-no, 1-yes 
plotSim = 1;

%% size and time plots
%Plot the size and time of all particles? 0-no, 1-yes 
plotSizeTime = 0;
%Plot the size by number fraction of all particles? 0-no, 1-yes  
plotSize = 0;

%get info for binning by size if plotSize == 1 or plotSizeTime == 1,
%otherwise don't worry about these
%all values should be um
%What is the min size (um) to plot?
minSize = 0; 
%What is the max size (um) to plot?
maxSize = 3;
%What # of size bins is desired? 
numBins = 6;

%get info for time plot is plotSizeTime == 1
%otherwise don't worry about these
%What is the starting time to plot? (format ex: 27-Oct-2015 15:45:40)
minTime = '14-Jul-2015 10:46:32';
%What is the ending time to plot? (format ex: 27-Oct-2015 15:45:40) 
maxTime = '30-Nov-2015 21:20:22';
%don't mess with this
scatterTimeBin = [datenum(minTime) datenum(maxTime)];

%ask user if want to show titles and legends for plots
%Do you want to do fast rendering (no titles, labels, legends are put onto plots)?\n This may give some speedups of graphics 0-no 1-yes 
plotRender = 0;
    
%ask user if want to show plots separately or make single plots
%Would you like all plots separate (0) or combined into a single window (1)? 
plotType = 1;

%don't mess with this one
plotTrack = [plotPS plotCS plotHM plotHMlog plotSizeTime plotSize plotSim plotRender plotType];
