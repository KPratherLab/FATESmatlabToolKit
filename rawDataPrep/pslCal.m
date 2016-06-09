function [calParams, speedModes] = pslCal(pslTime, pslVels, pslSizes, N, bins)
% Size calibration based on speed and spit out distribution modes in microns
% Call as [calParams, speedModes] = pslCal(pslTime, pslVels, pslSizes, N, bins)
%
% pslTime is a Nx2 matrix specifying the start and stop times for an N
% number of psl sizes run.  Start times are in column 1 and stop times are
% in column 2.
% pslVels is a Nx2 matrix specifying the min (column 1) and max (column 2) velocities to be considered
% times for an N number of psl sizes run.  
% pslSizes is a N vector that specifies the size of the PSL used in the
% calibration
% N indicates the order of the polynomial to be fit (default is 3)
% bins indicates the number or bins (default is 100)
% note depending on the width of the range given by each row may
% want to adjust the bin value (wider range, more bins)
% calParams are the new calibration parameters
% speedModes is the velocity of each psl size

% Modified by Alberto Cazorla. October 4, 2010
% Modified by Camille Sultana 2016
global PARTdataMat PARTidFlds PARTidMat PARTdataFlds

if nargin == 2
    N = 3;
    bins = 100;
elseif nargin == 3
    bins = 100;
end
if N < 1 && N > 5
    error('The order of the polynomial has to be a number between 1 and 5');
end

speedModes = zeros(size(pslTime,1),1);
for i = 1:size(pslTime,1)
    getSpeed = PARTdataMat(:,PARTdataFlds.VELOCITY) > pslVels(i,1) & PARTdataMat(:,PARTdataFlds.VELOCITY) < pslVels(i,2); 
    getTime = PARTidMat(:,PARTidFlds.TIME) > pslTime(i,1) & PARTidMat(:,PARTidFlds.TIME) < pslTime(i,2);
    Speed = PARTdataMat(getSpeed & getTime, PARTdataFlds.VELOCITY);
    % calculates mean and standar deviation of the speed
    avg   = mean(Speed);
    stdev = std(Speed);
    
    % Binning
    edges = 1:(bins+1);
    edges = avg-stdev+(2*stdev*(edges-1))/bins;
    centers = edges+2*stdev/(2*bins);
        
    histogramData = histc(Speed,edges);    % vector of speed histogram counts
    
    % find most frequent speed (mode in the histogram)
    [~, maxIDX] = max(histogramData);
    
    speedModes(i) = centers(maxIDX,1); % speed at the location of the distribution mode

    % speed histogram plots
    subplot(4,3,i), bar(edges,histogramData,'histc') 
    xlim([edges(1),edges(end)])
    xlabel('Speed (m/s)','FontSize',10);
    ylabel('frequency','FontSize',10);
    title(sprintf('%d nm PSL',round(pslSizes(i)*1000)),'FontSize',10);
    set(gca, 'FontSize', 6)
end

p = polyfit(speedModes, pslSizes, N);            % fits data to polynomial
x = linspace(min(speedModes), max(speedModes), 100); % creates data input for the polynomial evaluation
yFit = polyval(p,x);                                 % evaluate the polynomial
Rs = rSquared(pslSizes, polyval(p,speedModes));      % R-squared

% Return new calibration parameters based on the last fit
calParams = zeros(1,10);
for i=1:length(p)
    calParams(i) = p(length(p)-i+1);
end
calParams(9)  = min(speedModes)-50;
calParams(10) = max(speedModes)+25;

% plots
figure, plot(speedModes,pslSizes, 'o', x, yFit);
str = [];
for i=1:length(p)-1
    str = strcat(str, sprintf('%+g*X^%d ', p(i), length(p)-i));
end
str = strcat(str, sprintf('%+g', p(i+1)));
str = strcat(str, sprintf('\nR^2 = %.4f', Rs));
text(min(speedModes)+10,polyval(p,min(speedModes)),str,'FontSize',10); % print the equation and R^2 on the figure
xlabel('speed (m/s)');
ylabel('size (microns)')
title('Size Calibration Curve')
set(gca, 'FontSize', 6)