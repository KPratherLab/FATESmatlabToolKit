function Da = da_jayne(Coef,Velocity)
% DA_JAYNE calculates Da from velocity using logarithmic fit 
% Call as Da = DA_JAYNE(Coef,Velocity)
% where Da is particle aerodynamic diameter in um
%       Coef are the fitting coefficients and calibration limits
%       Coef = [D* vg b vmin vmax]
%       Velocity is in m/s
% Da is calculated as
% Da = D* * (vg/v - 1)^1/b
%
% Da = NaN for Velocity < Coef(4) or Velocity >= Coef(5)
% Velocity limits are given in m/s.
% 
% See also CALIB_DA_NOZ, CALIB_DA_RAW, DA_TSI

% YAADA - Software Toolkit to Analyze Single-Particle Mass Spectral Data
%
% Copyright (C) 2008 Jonathan O. Allen

% Jonathan O. Allen  2008-08-20

% Updated to not return NaNs when called with 3 parameters 
% as is done in calib_da_*
% JOA  2008-11-14

if nargin ~= 2
  error('Call as Da = DA_JAYNE(Coef,Velocity)');
end
if iscell(Coef)
  Coef = double(Coef{:});
end
if ~isvector(Coef)
  error('Fitting coefficients must be vector');
end

ThreeParam = 0;
if length(Coef) == 3 
  % called for optimization so don't return NaNs
  ThreeParam = 1;
elseif length(Coef) ~= 5
  error('Coefficient vector must be in form [D* vg b vmin vmax]');
end

Velocity = double(Velocity);

if ~ThreeParam
  idx = find(Velocity < Coef(4) | Velocity > Coef(5));
  Velocity(idx) = NaN;

  % replace imaginary results with NaNs
  idx = find((Coef(2)./Velocity - 1) < 0);
  Velocity(idx) = NaN;
end

Da = Coef(1) * (Coef(2)./Velocity - 1).^(1/Coef(3));

return
