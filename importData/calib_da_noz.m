function Param = calib_da_noz(NOZFile)
% Calibrate Da - velocity constants using data from .NOZ file
% Call as Param = CALIB_DA_NOZ(NOZFile)
% where NOZFile - full file name for .noz file
%       Param   - fitting parameters [D* vg b vmin vmax]
%
% See also CALIB_DA_RAW

% YAADA - Software Toolkit to Analyze Single-Particle Mass Spectral Data
%
% Copyright (C) 2008 Jonathan O. Allen

% Jonathan O. Allen  2008-09-04

if ~isword(NOZFile)
  error('Expecting word for NOZFile');
end

% check for stats package
if exist('nlinfit') ~=2
  error('Need statistics package to run CALIB_DA_NOZ');
end

% hardware parameters
% these should match values in *.inst
ScatterLength = 0.06; % distance between timing lasers (m)
TimerResolution = 50e-9; % time for one count (s)

% display velocities
MinVelocity = 50; % m/s
MaxVelocity = 150;

% read .noz file 
% relevant lines have the form
% Data Set 1=C:\SEASAW\D317\ATOFMS\20070326\Size_calibration\0.701\0.701;13569;0.701; 65.798

% get transit time counts (first number) and calibration Da (second number)

fid = fopen(NOZFile,'r');
if fid == -1
  error('Cannot find NOZ file %s',NOZFile);
end

t1 = textscan(fid,'%s','delimiter','\n');
fclose(fid);

% keep only lines which start 'Data Set'
t1 = t1{:};
DataSetLine  = regexp(t1,'^Data Set');
IsDataSetLine = ~cellfun('isempty',DataSetLine);
t1 = t1(IsDataSetLine);

if length(t1) == 0
  error('No Data Set lines in %s',NOZFile);
end

NumLine = length(t1);
MedianTTCounts = zeros(NumLine,1);
Da             = zeros(NumLine,1);
for i = 1:length(t1)
  [str,tt,da,num] = strread(t1{i},'%s%f%f%s','delimiter',';');
  MedianTTCounts(i) = tt;
  Da(i) = da;
end

disp(MedianTTCounts)
disp(Da)

MedianVelocity = ScatterLength ./ (MedianTTCounts * TimerResolution);

% fit data to Da = D* (vg/v - 1)^1/b (see Jayne et al 2000)
% Beta = [D* vg b]
Beta0 = [0.3 260 0.6];
[Beta,R,J,CovB,MSE] = nlinfit(MedianVelocity,Da,@da_jayne,Beta0);

% plot fit with confidence intervals
V = 10.^(log10(MinVelocity):0.01:log10(MaxVelocity));
[DaPred, DaDelta] = nlpredci(@da_jayne,V,Beta,R,'jacobian',J);

figure(1);
clf;
loglog(MedianVelocity,Da,'*');
hold on;
loglog(V,DaPred,'-k');
loglog(V,DaPred+DaDelta','--g');
loglog(V,DaPred-DaDelta','--g');

title(sprintf('Da Calibration Using %s',NOZFile),'Interpreter','none');
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
