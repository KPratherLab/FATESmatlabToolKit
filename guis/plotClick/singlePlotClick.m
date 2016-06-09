%% create tiled spectra/time/size plots based off user interaction with plotclick

%% initialize figure
tracker = 0;
q = figure(100+clickcnt);
set(q,'units','normalized');
set(q,'position',[.01 .07 .6 .83]);

%create legend at top showing all clusters selected
for i = 1:length(uniqueLabelsTMP)
    clr = SCATclr{LabelIdx(i)}; %get appropriate color
    rectangle('Position',[i-0.5,0,1,1],'FaceColor',clr,'EdgeColor',clr); 
end
ax = gca;
ax.XLim = [0.5 i+0.5];
ax.XTick = 1:1:length(uniqueLabelsTMP); %set up ticks
ax.XTickLabel = cellstr(num2str(uniqueLabelsTMP)); %label ticks
ax.FontSize = 8;
ax.YTick = [];
set(ax, 'Position', [0.03 0.97 0.96 0.02]);

%create legend at left showing percentage of each cluster selected
axes('Position', [0.03 0.05 0.01 0.90]);
rectangle('Position',[0 0 1 1]);
for i = 1:length(uniqueLabelsTMP)
    clr = SCATclr{LabelIdx(i)}; %get appropriate color
    rectangle('Position',[0,sum(fracCluster(i:end))-fracCluster(i),1,fracCluster(i)],'FaceColor',clr,'EdgeColor',clr);
end
ax = gca;
ax.XTick = []; %no ticks
t = title(sprintf('%i# %0.1f%%',clickcnt,100*length(indexInt)/length(XVectorUse))); %add selection number and % of particles selected
set(t,'Position',[1 -0.05]);

%% create tiled plots
%% %%%%%%%%%%create line plot of all mz spectra overlaid
if plotTrack(1) == 1
    tracker = tracker + 1;
    subplot(rownum,1,tracker,'replace')
    hold on
    for i = 1:length(uniqueLabelsTMP) %this method is much faster than plotting all lines individually, do not revert to that
        useMZdata = tmplabels == uniqueLabelsTMP(i);
        useMZdata = [scatterSPEC(:,indexInt(useMZdata)); nan(1,sum(useMZdata))];
        useMZ = repmat([scatterMZ; 1],size(useMZdata,2),1);
        plot(useMZ,useMZdata(:),'Color',SCATclr{LabelIdx(i)});
    end
    hold off
    if ~plotTrack(8)
        t = title(sprintf('Spectra Select %i', clickcnt),'FontSize',9);
    end
end

%%   %%%%%%%%%create line plot of average mz spectrum
if plotTrack(2) == 1
    tracker = tracker + 1;
    subplot(rownum,1,tracker,'replace');
    avgSpec = zeros(length(scatterMZ),length(uniqueLabelsTMP));
    hold on
    for i = 1:length(uniqueLabelsTMP) %get average of particles selected grouped by cluster
        useMZdata = tmplabels == uniqueLabelsTMP(i);
        avgSpec(:,i) = mean(scatterSPEC(:,indexInt(useMZdata)),2);
        plot(scatterMZ,avgSpec(:,i),'Color',SCATclr{LabelIdx(i)});
    end
    hold off
    if ~plotTrack(8)
        t = title(sprintf('Avg Spectra Select %i', clickcnt),'FontSize',9);
    end
end


%% %%%%%%%create heatmap of mz spectra
if plotTrack(3)
    tracker = tracker + 1;
    subplot(rownum,1,tracker,'replace');
    colormap(jet)
    useMZdata = scatterSPEC(:,findInt);
    maxVal = max(useMZdata(:)); %create colorbar limits for best image
    minVal = min(useMZdata(:));
    imagesc([scatterMZ(1) scatterMZ(end)],[1,size(useMZdata,2)],useMZdata',[minVal maxVal]);
    if ~plotTrack(8)
        t = title(sprintf('Spectra Heatmap Select %i', clickcnt),'FontSize',9);
    end
    
    %add legend to heat map indicating labels or "clusters"
    hold on
    for i = 1:length(uniqueLabelsTMP)
        rectangle('Position', [0 uniqueIDX(i)-0.5 5 uniqueIDX(i+1)-uniqueIDX(i)], 'FaceColor', SCATclr{LabelIdx(i)}, 'EdgeColor', SCATclr{LabelIdx(i)});
        text(2, floor(uniqueIDX(i)-0.5+(uniqueIDX(i+1)-uniqueIDX(i))/2), sprintf('%i',uniqueLabelsTMP(i)));
    end
    hold off
end

%% %%%%%%%create log10 heatmap of mz spectra
if plotTrack(4)
    tracker = tracker + 1;
    e = subplot(rownum,1,tracker,'replace');
    colormap(jet)
    useMZdata = scatterLogSPEC(:,findInt);
    maxVal = max(useMZdata(:)); %create colorbar limits for best image
    minVal = min(useMZdata(:));
    imagesc([scatterMZ(1) scatterMZ(end)],[1,size(useMZdata,2)],useMZdata',[minVal maxVal]);
    if ~plotTrack(8)
        t = title(sprintf('Spectra Log Heatmap Select %i', clickcnt),'FontSize',9);
    end
    
    %add legend to heat map indicating labels or "clusters"
    hold on
    for i = 1:length(uniqueLabelsTMP)
        rectangle('Position', [0 uniqueIDX(i)-0.5 5 uniqueIDX(i+1)-uniqueIDX(i)], 'FaceColor', SCATclr{LabelIdx(i)}, 'EdgeColor', SCATclr{LabelIdx(i)});
        text(2, floor(uniqueIDX(i)-0.5+(uniqueIDX(i+1)-uniqueIDX(i))/2), sprintf('%i',uniqueLabelsTMP(i)));
    end
    hold off
end

%% %%%%%%%%create scatter of size and time
if plotTrack(5)
    tracker = tracker + 1;
    j = subplot(rownum,1,tracker,'replace');
    gs  = gscatter(scatterTime(indexInt), scatterAllSize(indexInt),tmplabels,[],[],[],'off');
    %get appropriate color
    for i = 1:length(gs)
        gs(i).Color = SCATclr{LabelIdx(i)};
    end
    if ~plotTrack(8)
        ylabel('DA (um)');
        origSize = get(gca, 'Position');
        xlim(scatterTimeBin);
        ax = gca;
        ax.XTick = [scatterTimeBin(1):(scatterTimeBin(2)-scatterTimeBin(1))/4:scatterTimeBin(2)];
        datetick('x','mm/dd/yy HH','keeplimits','keepticks');
        set(gca,'Position', origSize);
        t = title(sprintf('TimeSize Select %i', clickcnt),'FontSize',9);
    end
end

%% %%%%%%%%create number fraction over size
if plotTrack(6)
    subplot(rownum,2,(rownum-1)*2+1,'replace')
    %bin all particles used by cluster by size
    clustSizeBinsTMP = zeros(length(scatterSIZEbin),length(uniqueLabelsTMP));
    for i = 1:length(uniqueLabelsTMP) %bin each cluster by size
        indexClust = tmplabels == uniqueLabelsTMP(i);
        clustSize = scatterAllSize(indexInt(indexClust));
        clustSizeBinsTMP(:,i) = histc(clustSize,scatterSIZEbin); %bin
    end
    
    %get number fraction relative to all particles in plot
    clustAllFraction = bsxfun(@rdivide, clustSizeBinsTMP, sum(scatterSIZE,2));
    clustAllFraction = bsxfun(@rdivide, clustAllFraction, max(clustAllFraction,[],1)); %normalize so that everything easily viewable on y axis from 0 to 1
    %get number fraction relative to particles within cluster
    [~,LabelIdx,LabelIdxTMP] = intersect(uniqueLabels,uniqueLabelsTMP);
    clustClustFraction = zeros(length(scatterSIZEbin),length(uniqueLabelsTMP));
    for i = 1:length(uniqueLabelsTMP)
        clustClustFraction(:,i) = clustSizeBinsTMP(:,LabelIdxTMP(i))./scatterSIZE(:,LabelIdx(i));
    end
    clustClustFraction = bsxfun(@rdivide, clustClustFraction, max(clustClustFraction,[],1)); %normalize so that everything easily viewable on y axis from 0 to 1    
    plotSize = scatterSIZEbin(1:(end-1))+(scatterSIZEbin(1)-scatterSIZEbin(2))/2;
    h = plot(plotSize, clustAllFraction(1:(end-1),:),'--',plotSize, clustClustFraction(1:(end-1),:));
    
    if ~plotTrack(8)
        t = title(sprintf('NumFraction All-Dotted Cluster-Solid %i', clickcnt),'FontSize',9);
        xlabel('DA (um)');
    end
    %get appropriate color
    for i = 1:(length(uniqueLabelsTMP))
        h(i).Color = SCATclr{LabelIdx(i)};
        h(i+length(uniqueLabelsTMP)).Color = SCATclr{LabelIdx(i)};
    end
    
end

%% %%%%%%%%create spectra similiarity plots
if plotTrack(7)
    subplot(rownum,2,(rownum-1)*2+2,'replace')
    avgSpec = zeros(length(scatterMZ),length(uniqueLabelsTMP));
    simTmp = zeros(length(tmplabels),1);
    for i = 1:length(uniqueLabelsTMP)  %multiply each avg spectra (grouped by cluster) by all particles selected in that cluster
        useMZdata = tmplabels == uniqueLabelsTMP(i);
        useSPEC = scatterSPEC(:,indexInt(useMZdata));
        avgSpec = mean(useSPEC,2);
        avgSpec = avgSpec/norm(avgSpec); %normalize so that multiplication yields dot product (max 1)
        simTmp(useMZdata) = useSPEC'*avgSpec;
    end
    %multiply avg spectra of all selected particles by all particle spectra
    useMZdata = scatterSPEC(:,findInt);
    avgAll = mean(useMZdata,2);
    simAll = useMZdata'*(avgAll/norm(avgAll));
    gs  = gscatter(simAll, simTmp,tmplabels,[],[],[],'off'); 
    axis tight
    %get appropriate color
    for i = 1:length(gs)
        gs(i).Color = SCATclr{LabelIdx(i)};
    end
    if ~plotTrack(8)
        title('Spectra Similarity','FontSize',9);
    end
end

%% %%%%%%%%save plots
if plotSave == 1
    if plotChoose == 0 %save plots automatically
        exportToPPTX('addslide');
        exportToPPTX('addpicture', q');
        savefig(fullfile(plotFolder,sprintf('select%i.fig',clickcnt)));
    else
        plotSaveNow = input('\nSave current figure? 0-no, 1-yes  '); %save select plots
        if plotSaveNow == 1
            savefig(fullfile(plotFolder,sprintf('select%i.fig',clickcnt)));
            exportToPPTX('addslide');
            exportToPPTX('addpicture', q');
        end
    end
end
    