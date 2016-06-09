function Da = da_tsi(Coef,Velocity)
% DA_TSI calculates Da from velocity using 3rd order polynomial in transit time
% Call as Da = DA_TSI(Coef,Velocity)
% where Da is vector of particle aerodynamic diameters in um
%       Coef are the logarithmic coefficients and calibration limits 
%       Velocity is a vector of particle velocities in m/s
% Da is calculated as 
%   Da = C(1) + C(2)*tt + C(3)*tt^2 + C(4)*tt^3 
%   for Min >= Velocity > Max
%   Da = NaN for Velocity < Min or Velocity >= Max
%   Coef = [C(1) C(2) C(3) C(4) Min Max]
% 
% See also DA_LOG, DA_POLY, DA_NOZ

% YAADA - Software Toolkit to Analyze Single-Particle Mass Spectral Data
%
% Copyright (C) 1999-2000 California Institute of Technology
% Copyright (C) 2001-2008 Arizona State University
% Copyright (C) 2008 Jonathan O. Allen

% Jonathan O. Allen  2005-10-04

% assumed hardware parameters
% these should match values in *.inst
% ScatterLength = 0.06; % distance between timing lasers (m)
% TimerResolution = 50e-9; % time for one count (s)

% Changed to Velocity = (Min Max)
% JOA  2008-08-11

if nargin ~= 2
  error('Call as Da = DA_TSI(Coef,Velocity)');
end
if iscell(Coef)
  Coef = double(Coef{:});
end
if min(size(Coef)) ~= 1
  error('Coefficients must be vector');
end
if max(size(Coef)) < 6
  error('Coefficient vector must include exponential terms and limit');
end

% calc transit time counts
%TTCount = ScatterLength ./ Velocity ./ TimerResolution;

N = length(Coef)-3;

% reverse order of polynomial coeffs for polyval function
c2 = Coef(N+1:-1:1);
Da = polyval(c2,Velocity);

OutsideLimit = find(Velocity < Coef(end-1) || Velocity > Coef(end));
Da(OutsideLimit) = NaN;

return
