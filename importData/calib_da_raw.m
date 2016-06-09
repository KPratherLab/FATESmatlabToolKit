% CALIB_DA_RAW
% Calibrate Da - velocity constants from raw data in .SEF files
% Edit this script to update list of calibration files and particle sizes.
% This script generates a quality assurance plots, you should check these 
% to ensure that peak transit time are correctly found for each calibration
% file.
%
% See also CALIB_DA_RAW

% YAADA - Software Toolkit to Analyze Single-Particle Mass Spectral Data
%
% Copyright (C) 2008 Jonathan O. Allen

% Jonathan O. Allen  2008-08-30

% check for stats package

% JC 2013 legacy update
if exist('nlinfit') ~=2
  error('Need statistics package to run CALIB_DA_NOZ');
end

% hardware parameters
% these should match values in *.inst
ScatterLength = 0.06; % distance between timing lasers (m)
TimerResolution = 50e-9; % time for one count (s)

% limit of valid velocities
MinVelocity = 50; % m/s
MaxVelocity = 150;

% transit time count limits
MaxTTCounts = ScatterLength / MinVelocity / TimerResolution;
MinTTCounts = ScatterLength / MaxVelocity / TimerResolution;

% calibration data files with time limits
% Update file location, file names, and calibration particle sizes below
% Start and Stop can be used to limit the particles used in each file by 
% acquisition times.
CalibDir = 'c:/data/d317d/calib';

%%% edit lines like the following ones for your calibration
% SEFFile{fi} = fullfile(CalibDir,'0.269/0.269.sef');
% Da(fi) = 0.269;

fi = 0;
fi = fi + 1;
SEFFile{fi} = fullfile(CalibDir,'0.269/0.269.sef');
Da(fi) = 0.269;
Start(fi) = NaN;
Stop(fi)  = NaN;

fi = fi + 1;
SEFFile{fi} = fullfile(CalibDir,'0.299/0.299.sef');
Da(fi) = 0.299;
Start(fi) = NaN;
Stop(fi)  = NaN;

fi = fi + 1;
SEFFile{fi} = fullfile(CalibDir,'0.350/0.350.sef');
Da(fi) = 0.350;
Start(fi) = NaN;
Stop(fi)  = NaN;

fi = fi + 1;
SEFFile{fi} = fullfile(CalibDir,'0.491/0.491.sef');
Da(fi) = 0.491;
Start(fi) = NaN;
Stop(fi)  = NaN;

fi = fi + 1;
SEFFile{fi} = fullfile(CalibDir,'0.701/0.701.sef');
Da(fi) = 0.701;
Start(fi) = NaN;
Stop(fi)  = NaN;

fi = fi + 1;
SEFFile{fi} = fullfile(CalibDir,'0.993/0.993.sef');
Da(fi) = 0.993;
Start(fi) = NaN;
Stop(fi)  = NaN;

% find peak transit times for each calibration file
for fi = 1:length(SEFFile)
  SEFData = read_sem(SEFFile{fi},NaN,NaN,1);

  SEFTime = SEFData(:,2);
  SEFTTCounts = SEFData(:,3);

  % select for time range
  if ~isnan(Start(fi))
    if ~isnan(Stop(fi))
      idx = range_search(SEFTime,'=',[Start(fi) Stop(fi)]);
    else
      idx = range_search(SEFTime,'=',[Start(fi) SEFTime(end)]);
    end
  else
    if ~isnan(Stop(fi))
      idx = range_search(SEFTime,'=',[SEFTime(1) Stop(fi)]);      
    else
      idx = 1:length(SEFTime);
    end
  end    
  
  idx2 = setdiff(1:length(SEFTime),idx,'rows','legacy'); %JC 2013
  SEFTime1 = SEFTime(idx2);   % keep track of removed data
  SEFTTCounts1 = SEFTTCounts(idx2);

  SEFTime = SEFTime(idx);
  SEFTTCounts = SEFTTCounts(idx);

  % remove particles with velocities outside physically reasonable limits
  idx = find(SEFTTCounts >= MinTTCounts & SEFTTCounts <= MaxTTCounts); 
  idx2 = setdiff(1:length(SEFTime),idx,'rows','legacy'); % JC 2013

  SEFTime2 = SEFTime(idx2);   % keep track of removed data
  SEFTTCounts2 = SEFTTCounts(idx2);

  SEFTime = SEFTime(idx);
  SEFTTCounts = SEFTTCounts(idx);
  
  % remove particles with low count transit times
  N = histc(SEFTTCounts,MinTTCounts:MaxTTCounts);
  MinN = median(N) + 5;
  LowN = MinTTCounts + find(N < MinN);
  idx2 = find(ismember(SEFTTCounts,LowN),'rows','legacy'); % JC 2013
  idx = setdiff(1:length(SEFTime),idx2, 'rows','legacy'); % JC 2013

  SEFTime3 = SEFTime(idx2);   % keep track of removed data
  SEFTTCounts3 = SEFTTCounts(idx2);

  SEFTime = SEFTime(idx);
  SEFTTCounts = SEFTTCounts(idx);
  MedianTTCounts(fi) = median(SEFTTCounts);

  % QA plot of raw data
  figure(2*fi-1);
  clf; 
  plot((SEFTime-floor(SEFTime(1)))*24,SEFTTCounts,'k.','markersize',2);
  LegText = {'Valid Points'};
  hold on;
  if ~isempty(SEFTime1)
    plot((SEFTime1-floor(SEFTime(1)))*24,SEFTTCounts1,'c.','markersize',2);
    LegText = {LegText{:}, 'Outside Time Range'}; 
  end
  if ~isempty(SEFTime2)
    plot((SEFTime2-floor(SEFTime(1)))*24,SEFTTCounts2,'g.','markersize',2);
    LegText = {LegText{:}, 'Outside TT Range'}; 
  end
  if ~isempty(SEFTime3)
    plot((SEFTime3-floor(SEFTime(1)))*24,SEFTTCounts3,'m.','markersize',2);
    LegText = {LegText{:}, 'Too Few TT Counts'}; 
  end

  xlabel('Hour of Day');
  ylabel('Transit Time Counts');
  title(sprintf('Da Calibration (%s)',SEFFile{fi}),'Interpreter','none');
  legend(LegText);
  
  figure(2*fi);
  clf; 
  hist(SEFTTCounts,300);
  n = hist(SEFTTCounts,300);
  hold on;
  plot(MedianTTCounts(fi),max(n),'r*');

  xlabel('Transit Time Counts');
  ylabel('Particle Counts');
  title(sprintf('Da Calibration (%s)',SEFFile{fi}),'Interpreter','none');
end

MedianVelocity = ScatterLength ./ (MedianTTCounts * TimerResolution);

% fit data to Da = D* (vg/v - 1)^1/b (see Jayne et al 2000)
% Beta = [D* vg b]
Beta0 = [0.3 260 0.6];
[Beta,R,J,CovB,MSE] = nlinfit(MedianVelocity,Da,@da_jayne,Beta0);

% plot fit with confidence intervals
V = 10.^(log10(MinVelocity):0.01:log10(MaxVelocity));
[DaPred, DaDelta] = nlpredci(@da_jayne,V,Beta,R,'jacobian',J);

figure(2*fi+1);
clf;
loglog(MedianVelocity,Da,'*');
hold on;
loglog(V,DaPred,'-k');
loglog(V,DaPred+DaDelta','--g');
loglog(V,DaPred-DaDelta','--g');

title(sprintf('Da Calibration Using Data in %s',CalibDir),'Interpreter','none');
xlabel('Velocity (m/s)');
ylabel('Aerodynamic Diameter, Da (\mum)');
lim = axis;
lim(1) = 30;
lim(2) = 300;
axis(lim);

% next step instructions
disp(sprintf('Fit calibration data to Da = D* (vg/v - 1)^1/b'));
disp(sprintf('Constants are D* = %f',Beta(1)));
disp(sprintf('              vg = %f',Beta(2)));
disp(sprintf('              b  = %f',Beta(3)));
disp(sprintf('              vmin  < %f',min(MedianVelocity)));
disp(sprintf('              vmax  > %f',max(MedianVelocity)));
disp('');
disp('You should update INST.DaCalibFunction and INST.DaCalibParam');
disp('then use update_da (see http://www.yaada.org/wiki for details)');

Param = [Beta min(MedianVelocity) max(MedianVelocity)];

return
