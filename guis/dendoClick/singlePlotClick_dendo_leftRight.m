%%this script creates tiled plots for the left and right branches selected during user interation with clickDendo

%% initialize figure
tracker = 0;
q = figure(100+clickcnt);
set(q,'units','normalized');
set(q,'position',[.01 .07 .6 .83]);

%create legend at top showing all clusters selected for each branch
for i = 1:length(leftIDs) %left branch
    rectangle('Position',[i-0.5,0,1,1],'FaceColor','red','EdgeColor','red');
end

for i = 1:length(rightIDs) %right branch
    rectangle('Position',[length(leftIDs)+i+0.5,0,1,1],'FaceColor','blue','EdgeColor','blue');
end

%label x ticks with cluster labels
ax = gca;
ax.XTick = 1:1:(length(leftIDs)+length(rightIDs)+1);
ax.XTickLabel = [cellstr(num2str(leftIDs')); ' '; cellstr(num2str(rightIDs'))];
ax.FontSize = 8;
ax.YTick = [];
ax.XLim = [0.5,(length(leftIDs)+length(rightIDs)+1.5)];
xlabel(ax,'<--left IDs     right IDs-->'); 
set(ax, 'Position', [0.03 0.97 0.96 0.02]);

%create legend at left showing percentage of each branch and containing clusters
axes('Position', [0.03 0.05 0.01 0.90]);
rectangle('Position',[0,0,1,leftFraction],'EdgeColor','red','FaceColor','red');
rectangle('Position',[0,leftFraction,1,rightFraction],'EdgeColor','blue','FaceColor','blue');
ax = gca;
ax.XTick = [];
t = title(sprintf('%.2f %%',100*(sum(leftCounts)+sum(rightCounts))/sum(dendoCOUNTS)));
set(t,'Position',[1 -0.05]);

%% create tiled plots

%% create line plot of all mz spectra overlaid
%note this plot is very time/memory consuming
if plotTrack(1)
    tracker = tracker + 1;
    subplot(rownum,1,tracker,'replace')
    useMZdata = cell2mat(dendoSPEC(leftIDs));
    %this method is much faster than plotting all lines individually, do not revert to that
    useMZdata = [useMZdata; nan(1,size(useMZdata,2))]; 
    useMZtmp = repmat([dendoMZ; 1],size(useMZdata,2),1);
    plot(useMZtmp,useMZdata(:),'Color','red');
    hold on
    useMZdata = cell2mat(dendoSPEC(rightIDs));
    useMZdata = [useMZdata; nan(1,size(useMZdata,2))];
    useMZtmp = repmat([dendoMZ; 1],size(useMZdata,2),1);
    plot(useMZtmp,useMZdata(:),'Color','blue');
    axis tight
    hold off
    if ~plotTrack(8)
        t = title('All Spectra','FontSize',9);
    end
end

%% create plot of spectra averages
if plotTrack(2)
    tracker = tracker + 1;
    subplot(rownum,1,tracker,'replace')
    leftAvg = sum(bsxfun(@times, cell2mat(dendoSPECavg(leftIDs)), leftCounts/sum(leftCounts)),2);
    plot(dendoMZ,leftAvg,'Color','red');
    hold on
    rightAvg = sum(bsxfun(@times, cell2mat(dendoSPECavg(rightIDs)), rightCounts/sum(rightCounts)),2);
    plot(dendoMZ,rightAvg,'Color','blue');
    axis tight
    allAvg = leftAvg*leftFraction+rightAvg*rightFraction;
    plot(dendoMZ,allAvg,'Color','black');
    hold off
    if ~plotTrack(8)
        t = title('Avg Spectra','FontSize',9);
    end
end

%% create heatmap of mz spectra
if plotTrack(3)
    tracker = tracker + 1;
    subplot(rownum,1,tracker,'replace')
    colormap(jet)
    useMZdata = cell2mat(dendoSPEC([leftIDs rightIDs]));
    maxVal = max(useMZdata(:));
    minVal = min(useMZdata(:));
    imagesc([dendoMZ(1) dendoMZ(end)],[1,size(useMZdata,2)], useMZdata',[minVal maxVal]);
    if ~plotTrack(8)
        t = title('Spectra HeatMap','FontSize',9);
    end
    
    %add legend to heat map
    hold on
    rectangle('Position', [0 0 5 sum(leftCounts)], 'FaceColor', 'red', 'EdgeColor', 'red');
    rectangle('Position', [0 sum(leftCounts) 5 sum(rightCounts)+sum(leftCounts)], 'FaceColor', 'blue', 'EdgeColor', 'blue');
    hold off
end

%% create log10 heatmap of mz spectra
if plotTrack(4)
    tracker = tracker + 1;
    subplot(rownum,1,tracker,'replace')
    colormap(jet)
    useMZdata = cell2mat(dendoLogSPEC([leftIDs rightIDs]));
    maxVal = max(useMZdata(:));
    minVal = min(useMZdata(:));
    imagesc([dendoMZ(1) dendoMZ(end)],[1,size(useMZdata,2)], useMZdata',[minVal maxVal]);
    if ~plotTrack(8)
        t = title('Spectra Log HeatMap','FontSize',9);
    end
    
    %add legend to heat map
    hold on
    rectangle('Position', [0 0 5 sum(leftCounts)], 'FaceColor', 'red', 'EdgeColor', 'red');
    rectangle('Position', [0 sum(leftCounts) 5 sum(rightCounts)+sum(leftCounts)], 'FaceColor', 'blue', 'EdgeColor', 'blue');
    hold off
end

%% create relative temporal
if plotTrack(5)
    tracker = tracker + 1;
    subplot(rownum,1,tracker,'replace')
    leftTime = sum(dendoTime(:,leftIDs),2);
    rightTime = sum(dendoTime(:,rightIDs),2);
    a = plot(dendoTimeBin,[leftTime/max(leftTime) rightTime/max(rightTime)],'LineStyle','-','Marker','.');
    a(1).Color = 'red';
    a(2).Color = 'blue';
    xlim([dendoTimeBin(1) dendoTimeBin(end)]);
    if ~plotTrack(8)
        ax.XTick = [dendoTimeBin(1):(dendoTimeBin(end)-dendoTimeBin(1))/5:dendoTimeBin(end)];
        ax.XTickLabelRotation = 45;
        datetick('x','mm/dd/yy HH','keeplimits','keepticks');
        title('Left/Right Relative Temporal','FontSize',9);
    end
end

%% create number fraction over size
if plotTrack(6)
    subplot(rownum,2,(rownum-1)*2+1,'replace')
    leftSize = sum(dendoSize(:,leftIDs),2);
    rightSize = sum(dendoSize(:,rightIDs),2);
    a = plot(dendoSizeBin,[leftSize/max(leftSize) rightSize/max(rightSize)]); %normalize counts so can plot on y scale 0-1
    a(1).Color = 'red';
    a(2).Color = 'blue';
    if ~plotTrack(8)
        title('Left/Right Size Number Fraction','FontSize',9);
    end
end

%% create spectra similiarity plots
if plotTrack(7)
    subplot(rownum,2,(rownum-1)*2+2,'replace')
    if plotTrack(2) == 0
        leftAvg = sum(bsxfun(@times, cell2mat(dendoSPECavg(leftIDs)), leftCounts/sum(leftCounts)),2);
        rightAvg = sum(bsxfun(@times, cell2mat(dendoSPECavg(rightsIDs)), leftCounts/sum(leftCounts)),2);
        allAvg = leftAvg*leftFraction+rightAvg*rightFraction;
    end
    %normalize counts so can plot on y scale 0-1 and see shape of curves
    %easily
    SimLeft = cell2mat(dendoSPEC(leftIDs))'*leftAvg;
    SimRight = cell2mat(dendoSPEC(rightIDs))'*rightAvg;
    SimLavg = cell2mat(dendoSPEC(leftIDs))'*allAvg;
    SimRavg = cell2mat(dendoSPEC(rightIDs))'*allAvg;
    scatter(SimLavg,SimLeft,'.','red');
    hold on
    scatter(SimRavg,SimRight,'.','blue');
    hold off
    if ~plotTrack(8)
        title('Left/Right Spectra Similarity','FontSize',9);
    end
end

%% save plot
if savePlots
    exportToPPTX('addslide');
    exportToPPTX('addpicture', q');
    savefig(fullfile(plotFolder,sprintf('dendo%i',100+clickcnt)));
end
    