%%this script creates tiled plots for all clusters selected during user interation with clickDendo

%% initialize figure
tracker = 0;
q = figure(100+clickcnt);
set(q,'units','normalized');
set(q,'position',[.01 .07 .6 .83]);

%create legend at top showing all clusters selected for each branch
for i = 1:length(leftIDs) %left branch clusters
    clr = useCtmpL{i};
    rectangle('Position',[i-0.5,0,1,1],'FaceColor',clr,'EdgeColor',clr);
end

for i = 1:length(rightIDs) %right branch clusters
    clr = useCtmpR{i};
    rectangle('Position',[i+0.5+length(leftIDs),0,1,1],'FaceColor',clr,'EdgeColor',clr);
end

%label x ticks with cluster labels
ax = gca;
ax.XTick = 1:1:(length(allIDs)+1);
ax.XTickLabel = [cellstr(num2str(IDstr'))];
ax.FontSize = 8;
ax.YTick = [];
ax.XLim = [0.5,length(allIDs)+1.5];
set(ax, 'Position', [0.03 0.97 0.96 0.02]);

%create legend at right showing percentage of each branch and containing clusters
axes('Position', [0.03 0.05 0.01 0.90]);
for i = 1:length(allIDs) % percentage of each cluster
    clr = useCtmp{i};
    rectangle('Position',[0,sum(fracCluster(1:i))-fracCluster(i),1,fracCluster(i)],'FaceColor',clr,'EdgeColor',clr);
end
ax = gca;
ax.XTick = [];

%create boxes framing all left and all right clusters
rectangle('Position',[0,0,1,leftFraction],'EdgeColor',[0 0 0],'LineWidth',2);
rectangle('Position',[0,leftFraction,1,rightFraction],'EdgeColor',[0 0 0],'LineWidth',2);
ax.YTick = [leftFraction/2, leftFraction+rightFraction/2];
ax.YTickLabel = {'left branch' 'right branch'};
ax.YTickLabelRotation = 90;
ax.YLim = [0 1];

%fraction of particles selected out of all in dendogram
t = title(sprintf('%.2f %%',100*sum(allCounts)/sum(dendoCOUNTS)));
set(t,'Position',[1 -0.05]);


%% create tiled plots 

%% create line plot of all mz spectra overlaid
%note this plot is very time/memory consuming
if plotTrack(1)
    %initialize subplot
    tracker = tracker + 1;
    subplot(rownum,1,tracker,'replace')
    Ctrack = 0;
    hold on
    %plot left branch clusters
    for i = 1:length(leftIDs) %this method is much faster than plotting all lines individually, do not revert to that
        useMZdata = [dendoSPEC{leftIDs(i)};nan(1,size(dendoSPEC{leftIDs(i)},2))];
        useMZtmp = repmat([dendoMZ; 1],size(dendoSPEC{leftIDs(i)},2),1);
        plot(useMZtmp,useMZdata(:),'Color',useCtmpL{i});
    end
    
    %plot right branch clusters
    for i = 1:length(rightIDs) %this method is much faster than plotting all lines individually, do not revert to that
        useMZdata = [dendoSPEC{rightIDs(i)};nan(1,size(dendoSPEC{rightIDs(i)},2))];
        useMZtmp = repmat([dendoMZ; 1],size(dendoSPEC{rightIDs(i)},2),1);
        plot(useMZtmp,useMZdata(:),'Color',useCtmpR{i},'LineStyle','--');
    end
    if ~plotTrack(8)
        t = title('All Spectra','FontSize',9);
    end
    hold off
end

%% create plot of spectra averages
if plotTrack(2)
    %initialize subplot
    tracker = tracker + 1;
    subplot(rownum,1,tracker,'replace')
    %get data
    useMZdata = cell2mat(dendoSPECavg(leftIDs));
    %make plot left branch clusters
    plot(dendoMZ,useMZdata);     
    hold on
    %get data
    useMZdata = cell2mat(dendoSPECavg(rightIDs));
    %make plot right branch clusters
    plot(dendoMZ,useMZdata,'LineStyle','--');
    %get average of all clusters
    allAvg = sum(bsxfun(@times, cell2mat(dendoSPECavg(allIDs)), fracCluster),2);
    %plot average of all clusters
    plot(dendoMZ,allAvg,'Color','black','LineStyle',':')     

    if ~plotTrack(8)
        t = title('Avg Spectra','FontSize',9);
    end
    hold off
end
    
%% create heatmap of mz spectra
if plotTrack(3)    
    %initialize subplot
    tracker = tracker + 1;
    subplot(rownum,1,tracker,'replace')
    colormap(jet) %easier to see differences with this map
    %get data
    useMZdata = cell2mat(dendoSPEC(allIDs));
    %set color limits
    maxVal = max(useMZdata(:));
    minVal = min(useMZdata(:));
    %make plot
    imagesc([dendoMZ(1) dendoMZ(end)],[1,size(useMZdata,2)], useMZdata',[minVal maxVal]);
    if ~plotTrack(8)
        t = title('Spectra HeatMap','FontSize',9);
    end
    
    %add legend to heat map
    hold on
    for i = 1:length(allIDs)
        rectangle('Position', [0 sum(allCounts(1:i))-allCounts(i)+0.5 5 allCounts(i)], 'FaceColor', useCtmp{i}, 'EdgeColor', useCtmp{i});
    end
    hold off    
end

%% create log10 heatmap of mz spectra
if plotTrack(4)
    %initialize subplot
    tracker = tracker + 1;
    subplot(rownum,1,tracker,'replace')
    colormap(jet)
    %get data
    useMZdata = cell2mat(dendoLogSPEC(allIDs));
    %set color limits
    maxVal = max(useMZdata(:)); 
    minVal = min(useMZdata(:));
    %make plot
    imagesc([dendoMZ(1) dendoMZ(end)],[1,size(useMZdata,2)], useMZdata',[minVal maxVal]);
    if ~plotTrack(8)
        t = title('Spectra HeatMap','FontSize',9);
    end
    
    %add legend to heat map
    hold on
    for i = 1:length(allIDs)
        rectangle('Position', [0 sum(allCounts(1:i))-allCounts(i)+0.5 5 allCounts(i)], 'FaceColor', useCtmp{i}, 'EdgeColor', useCtmp{i});
    end
    hold off   
end
    
%% create relative temporal
if plotTrack(5)
    %initialize subplot
    tracker = tracker + 1;
    subplot(rownum,1,tracker,'replace')    
    %get time for clusters
    leftTime = dendoTime(:,leftIDs);
    rightTime = dendoTime(:,rightIDs);
    %normalize counts so everything can plot on 0-1 scale
    leftTime = bsxfun(@rdivide,leftTime,max(leftTime,[],1));
    rightTime = bsxfun(@rdivide,rightTime,max(rightTime,[],1));
    %make plot
    plot(dendoTimeBin,leftTime,'LineStyle','-','Marker','.');
    hold on
    plot(dendoTimeBin,rightTime,'LineStyle','--','Marker','.');
    if ~plotTrack(8)
        ax.XTick = [dendoTimeBin(1):(dendoTimeBin(end)-dendoTimeBin(1))/5:dendoTimeBin(end)];
        ax.XTickLabelRotation = 45;
        datetick('x','mm/dd/yy HH','keeplimits','keepticks');
        title('Cluster Relative Temporal','FontSize',9);
    end
end

%% create number fraction over size
if plotTrack(6)
    %initialize subplot
    subplot(rownum,2,(rownum-1)*2+1,'replace')
    %get size for clusters
    leftSize = dendoSize(:,leftIDs);
    rightSize = dendoSize(:,rightIDs);
    %normalize counts so everything can plot on 0-1 scale
    leftSize = bsxfun(@rdivide,leftSize,max(leftSize,[],1));
    rightSize = bsxfun(@rdivide,rightSize,max(rightSize,[],1));
    %make plot
    plot(dendoSizeBin,leftSize);
    hold on
    plot(dendoSizeBin,rightSize,'LineStyle','--');
    if ~plotTrack(8)
        title('Cluster Size Number Fraction','FontSize',9);
    end
    hold off
end

%% create spectra similiarity plots
if plotTrack(7)
    %initialize subplot
    subplot(rownum,2,(rownum-1)*2+2,'replace')
    %get average of all clusters
    if plotTrack(2) == 0
        allAvg = sum(bsxfun(@times, cell2mat(dendoSPECavg(allIDs)), fracCluster),2);
    end
    hold on
    %make plot
    for i = 1:length(leftIDs)
        tmpSim = dendoSPEC{leftIDs(i)}'*allAvg;
        scatter(tmpSim,dendoSIM{leftIDs(i)},'.')
    end
    
    for i = 1:length(rightIDs)
        tmpSim = dendoSPEC{rightIDs(i)}'*allAvg;
        scatter(tmpSim,dendoSIM{rightIDs(i)},15,'+')
    end
    hold off
    
    if ~plotTrack(8)
        title('Cluster Spectra Similarity','FontSize',9);
    end
end

%% save plot
if savePlots
    exportToPPTX('addslide');
    exportToPPTX('addpicture', q');
    savefig(fullfile(plotFolder,sprintf('dendo%i',100+clickcnt)));
end
