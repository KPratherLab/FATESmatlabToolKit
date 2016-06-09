%% create spectra/time/size plots based off user interaction with plotclick

%% %%%%%%%%%%create line plot of all mz spectra overlaid
if plotTrack(1) == 1
    q = figure(100+clickcnt);
    hold on
    for i = 1:length(uniqueLabelsTMP) %this method is much faster than plotting all lines individually, do not revert to that
        useMZdata = tmplabels == uniqueLabelsTMP(i);
        useMZdata = [scatterSPEC(:,indexInt(useMZdata)); nan(1,sum(useMZdata))];
        useMZ = repmat([scatterMZ; 1],size(useMZdata,2),1);
        plot(useMZ,useMZdata(:),'Color',SCATclr{LabelIdx(i)});
    end
    hold off
    Fname = sprintf('allSpectra%i', clickcnt);
    if ~plotTrack(8)
        t = title(Fname,'FontSize',9);
    end
    savePlotClick %save plot
end

%%   %%%%%%%%%create line plot of average mz spectrum
if plotTrack(2) == 1
    q = figure(200+clickcnt);
    avgSpec = zeros(length(scatterMZ),length(uniqueLabelsTMP));
    hold on
    for i = 1:length(uniqueLabelsTMP) %get average of particles selected grouped by cluster
        useMZdata = tmplabels == uniqueLabelsTMP(i);
        avgSpec(:,i) = mean(scatterSPEC(:,indexInt(useMZdata)),2);
        plot(scatterMZ,avgSpec(:,i),'Color',SCATclr{LabelIdx(i)});
    end
    hold off
    Fname = sprintf('avgSpectra%i',clickcnt);
    if ~plotTrack(8)
        t = title(Fname,'FontSize',9);
    end
    savePlotClick %save plot      
end


%% %%%%%%%create heatmap of mz spectra
if plotTrack(3)
    q = figure(300+clickcnt);
    colormap(jet)
    useMZdata = scatterSPEC(:,findInt);
    maxVal = max(useMZdata(:)); %create colorbar limits for best image
    minVal = min(useMZdata(:));
    imagesc([scatterMZ(1) scatterMZ(end)],[1,size(useMZdata,2)],useMZdata',[minVal maxVal]);
    Fname = sprintf('SpectraHeatmap%i',clickcnt);
    if ~plotTrack(8)
        t = title(Fname,'FontSize',9);
    end
    
    %add legend to heat map indicating labels or "clusters"
    hold on
    for i = 1:length(uniqueLabelsTMP)
        rectangle('Position', [0 uniqueIDX(i)-0.5 5 uniqueIDX(i+1)-uniqueIDX(i)], 'FaceColor', SCATclr{LabelIdx(i)}, 'EdgeColor', SCATclr{LabelIdx(i)});
        text(2, floor(uniqueIDX(i)-0.5+(uniqueIDX(i+1)-uniqueIDX(i))/2), sprintf('%i',uniqueLabelsTMP(i)));
    end
    hold off
    savePlotClick %save plot  
end

%% %%%%%%%create log10 heatmap of mz spectra
if plotTrack(4)
    q = figure(400+clickcnt);
    colormap(jet)
    useMZdata = scatterLogSPEC(:,findInt);
    maxVal = max(useMZdata(:)); %create colorbar limits for best image
    minVal = min(useMZdata(:));
    imagesc([scatterMZ(1) scatterMZ(end)],[1,size(useMZdata,2)],useMZdata',[minVal maxVal]);
    Fname = sprintf('SpectraLogHeatmap%i',clickcnt);
    if ~plotTrack(8)
        t = title(Fname,'FontSize',9);
    end
    
    %add legend to heat map indicating labels or "clusters"
    hold on
    for i = 1:length(uniqueLabelsTMP)
        rectangle('Position', [0 uniqueIDX(i)-0.5 5 uniqueIDX(i+1)-uniqueIDX(i)], 'FaceColor', SCATclr{LabelIdx(i)}, 'EdgeColor', SCATclr{LabelIdx(i)});
        text(2, floor(uniqueIDX(i)-0.5+(uniqueIDX(i+1)-uniqueIDX(i))/2), sprintf('%i',uniqueLabelsTMP(i)));
    end
    hold off
    savePlotClick %save plot  
end

%% %%%%%%%%create scatter of size and time
if plotTrack(5)
    q = figure(500+clickcnt);
    gs  = gscatter(scatterTime(indexInt), scatterAllSize(indexInt),tmplabels,[],[],[],'off');
    %get appropriate color
    for i = 1:length(gs)
        gs(i).Color = SCATclr{LabelIdx(i)};
    end
    Fname = sprintf('TimeSize%i',clickcnt);
    if ~plotTrack(8)
        ylabel('DA (um)');
        origSize = get(gca, 'Position');
        xlim(scatterTimeBin);
        ax = gca;
        ax.XTick = [scatterTimeBin(1):(scatterTimeBin(2)-scatterTimeBin(1))/4:scatterTimeBin(2)];
        datetick('x','mm/dd/yy HH','keeplimits','keepticks');
        set(gca,'Position', origSize);
        t = title(Fname,'FontSize',9);
    end
    savePlotClick %save plot      
end

%% %%%%%%%%create number fraction over size
if plotTrack(6)
    q = figure(600+clickcnt);
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
    Fname = sprintf('Size%i',clickcnt);
    if ~plotTrack(8)
        t = title(Fname,'FontSize',9);
        xlabel('DA (um)');
    end
    %get appropriate color
    for i = 1:(length(uniqueLabelsTMP))
        h(i).Color = SCATclr{LabelIdx(i)};
        h(i+length(uniqueLabelsTMP)).Color = SCATclr{LabelIdx(i)};
    end
    savePlotClick %save plot      
end

%% %%%%%%%%create spectra similiarity plots
if plotTrack(7)
    q = figure(700+clickcnt);
    avgSpec = zeros(length(scatterMZ),length(uniqueLabelsTMP));
    simTmp = zeros(length(tmplabels),1);
    for i = 1:length(uniqueLabelsTMP) %multiply each avg spectra (grouped by cluster) by all particles selected in that cluster
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
    Fname = sprintf('Similarity%i',clickcnt);
    if ~plotTrack(8)
        title(Fname,'FontSize',9);
    end
    savePlotClick %save plot  
end

%%
%create scatter of size and time
if plotSizeTime == 1
    q = figure(900+clickcnt);
        gs  = gscatter(PARTidMat(trackintersectPART(indexInt),PARTidFlds.TIME), allSize(indexInt),tmplabels,[],[],[],'off');  % sym color?
        %get appropriate color
        for i = 1:length(gs)
            gs(i).Color = SCATclr{LabelIdx(i)};
        end
        if ~plotRender
            ylabel('DA (um)');
            origSize = get(gca, 'Position');
            xlim([minTime maxTime]);
            ax = gca;
            ax.XTick = [minTime:(maxTime-minTime)/4:maxTime];
            ax.XTickLabelRotation = 45;
            datetick('x','mm/dd/yy HH','keeplimits','keepticks');
            set(gca,'Position', origSize);
            Fname = sprintf('TimeSize Select %i', clickcnt);
            t = title(Fname);
        end
    savePlotClick %save plot  
end

%% 
%create number fraction over size
if plotSize == 1
    
    %bin all particles used by cluster by size
    clustSizeBinsTMP = zeros(length(uniqueLabelsTMP),length(allSizeBins));
    for i = 1:length(uniqueLabelsTMP) %bin each cluster by size
        indexClust = tmplabels == uniqueLabelsTMP(i);
        clustSize = allSize(indexInt(indexClust));
        clustSizeBinsTMP(i,:) = histc(clustSize,sizeBins); %bin
    end
    
    %get number fraction relative to all particles in plot
    clustAllFraction = bsxfun(@rdivide, clustSizeBinsTMP, allSizeBins');
    
    %get number fraction relative to particles within cluster
    [~,LabelIdx,LabelIdxTMP] = intersect(uniqueLabels,uniqueLabelsTMP);
    clustClustFraction = zeros(length(uniqueLabelsTMP),length(allSizeBins));
    for i = 1:length(uniqueLabelsTMP)
        clustClustFraction(i,:) = clustSizeBinsTMP(LabelIdxTMP(i),:)./clustSizeBins(LabelIdx(i),:);
    end
    
    %create plot of number fraction
    q = figure(1000+clickcnt);
    h = plot(sizeBins(1:(end-1))+binStep/2, clustAllFraction(:,1:(end-1)),'--',sizeBins(1:(end-1))+binStep/2, clustClustFraction(:,1:(end-1)));
    %recolor lines
    for i = 1:(length(uniqueLabelsTMP))
        h(i).Color = SCATclr{LabelIdx(i)};
        h(i+length(uniqueLabelsTMP)).Color = SCATclr{LabelIdx(i)};
    end
    if ~plotRender
        %add legend
        legendText = cell(1,length(h));
        for i = 1:(length(uniqueLabelsTMP))
            legendText{i} = sprintf('all%i', uniqueLabels(LabelIdx(i)));
            legendText{i+length(uniqueLabelsTMP)} = sprintf('clust%i', uniqueLabels(LabelIdx(i)));
        end
        legend(legendText,'Location','best');
        xlabel('DA (um)');
        Fname = sprintf('NumberFracSize%i', clickcnt);
        title(Fname);
    end
    savePlotClick %save plot  
end