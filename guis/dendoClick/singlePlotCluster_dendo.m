%%this script creates tiled plots for the single cluster selected during user interation with clickDendo

%% initialize figure
tracker = 0;
q = figure(100+clickcnt);
set(q,'units','normalized');
set(q,'position',[.01 .07 .6 .83]);

%% create tiled plots

%% create line plot of all mz spectra overlaid
%note this plot is very time/memory consuming
if plotTrack(1)
    tracker = tracker + 1;
    subplot(rownum,1,tracker,'replace')
    useMZdata = dendoSPEC{clustNUM};
    useMZdata = [useMZdata; nan(1,size(useMZdata,2))]; %this method is much faster than plotting all lines individually, do not revert to that
    useMZtmp = repmat([dendoMZ; 1],size(useMZdata,2),1);
    plot(useMZtmp,useMZdata(:));
    axis tight
    if ~plotTrack(8)
        t = title(sprintf('Spectra Cluster %i',clustNUM),'FontSize',9);
    end
end

%% create plot of spectra average
if plotTrack(2)
    tracker = tracker + 1;
    subplot(rownum,1,tracker,'replace')
    plot(dendoMZ,dendoSPECavg{clustNUM});
    if ~plotTrack(8)
        t = title(sprintf('Avg Spectra Cluster %i',clustNUM),'FontSize',9);
    end
end

%% create heatmap of mz spectra
if plotTrack(3)
    tracker = tracker + 1;
    subplot(rownum,1,tracker,'replace')
    colormap(jet)
    useMZdata = dendoSPEC{clustNUM};
    maxVal = max(useMZdata(:));
    minVal = min(useMZdata(:));
    imagesc([dendoMZ(1) dendoMZ(end)],[1,size(useMZdata,2)], useMZdata',[minVal maxVal]);
    if ~plotTrack(8)
        t = title(sprintf('Spectra Heatmap Cluster %i',clustNUM),'FontSize',9);
    end
end

%% create log10 heatmap of mz spectra
if plotTrack(4)
    tracker = tracker + 1;
    subplot(rownum,1,tracker,'replace')
    colormap(jet)
    useMZdata = dendoLogSPEC{clustNUM};
    maxVal = max(useMZdata(:));
    minVal = min(useMZdata(:));
    imagesc([dendoMZ(1) dendoMZ(end)],[1,size(useMZdata,2)], useMZdata',[minVal maxVal]);
    if ~plotTrack(8)
        t = title(sprintf('Spectra Log Heatmap Cluster %i',clustNUM),'FontSize',9);
    end
end

%% create relative temporal
if plotTrack(5)
    tracker = tracker + 1;
    subplot(rownum,1,tracker,'replace')
    a = plot(dendoTimeBin,dendoTime(:,clustNUM),'LineStyle','-','Marker','.');
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
    a = plot(dendoSizeBin,dendoSize(:,clustNUM));
    if ~plotTrack(8)
        title('Left/Right Size Number Fraction','FontSize',9);
    end
end

%% create spectra similiarity plots
% if plotTrack(7)
%     subplot(rownum,2,(rownum-1)*2+2,'replace')
%     if plotTrack(2) == 0
%         leftAvg = sum(bsxfun(@times, cell2mat(dendoSPECavg(leftIDs)), leftCounts/sum(leftCounts)),2);
%         rightAvg = sum(bsxfun(@times, cell2mat(dendoSPECavg(rightsIDs)), leftCounts/sum(leftCounts)),2);
%         allAvg = leftAvg*leftFraction+rightAvg*rightFraction;
%     end
%     SimLeft = cell2mat(dendoSPEC(leftIDs))'*leftAvg;
%     SimRight = cell2mat(dendoSPEC(rightIDs))'*rightAvg;
%     SimLavg = cell2mat(dendoSPEC(leftIDs))'*allAvg;
%     SimRavg = cell2mat(dendoSPEC(rightIDs))'*allAvg;
%     scatter(SimLavg,SimLeft,'.','red');
%     hold on
%     scatter(SimRavg,SimRight,'.','blue');
%     hold off
%     if ~plotTrack(8)
%         title('Left/Right Spectra Similarity','FontSize',9);
%     end
% end

if savePlots
    exportToPPTX('addslide');
    exportToPPTX('addpicture', q');
    savefig(fullfile(plotFolder,sprintf('dendo%i',100+clickcnt)));
end
