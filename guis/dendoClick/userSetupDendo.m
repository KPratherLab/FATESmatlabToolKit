%user file to automatically populate all user input variables needed to get started
%with clickHist

%% figure out basic plots to make
% Do you want to make spectra plots? 0-no, 1-yes
%Set to 0 if useSPEC is empty 
plotSpec = 1;
%Do you want to plot clusters over time? 0-no 1-yes
plotTime = 0;
%Do you want to plot clusters size number fraction? 0-no 1-yes
plotSize = 0;

%% if plotSpec is 1 then the following can be set to 0 or 1.
% if plotSpec is 0 then set all of the following to 0
% Plot all particle spectra overlaid? 0-no, 1-yes
plotPS = 0;
%Plot average spectra of each cluster overlaid? 0-no, 1-yes
plotCS = 1;
%Create the heat map of all particle spectra? 0-no, 1-yes  
plotHM = 1;
%Create the heat map of all particle spectra? 0-no, 1-yes 
plotHMlog = 1;
%Plot the similarity of spectra to averages? 0-no, 1-yes
plotSim = 0;

%% How to render plots
%Do you want to do fast rendering (no titles, labels, legends are put onto plots)?\n This may give some speedups of graphics 0-no 1-yes  '
plotRender = 0;

%%
%don't mess with this
plotTrack = [plotPS plotCS plotHM plotHMlog plotTime plotSize plotSim plotRender];
    