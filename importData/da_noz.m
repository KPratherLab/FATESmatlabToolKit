function Da = da_noz(Coef,Velocity)
% DA_NOZ calculates Da from velocity using mixed log and polynomial
% Call as Da = DA_NOZ(Coef,Velocity)
% where Da is vector of particle aerodynamic diameters in um
%       Coef are the logarithmic coefficients and calibration limits 
%       Velocity is a vector of particle velocities in m/s
% Da is calculated as 
%   Da = C(1) + C(2)*v + C(3)*v^2 ... + C(N+1)*v^N + C(N+2)*e^(C(N+3)*Velocity)
%   for Min >= Velocity > Max
%   Da = NaN for Velocity < Min or Velocity >= Max
%   Coef = [C(1) C(2) C(3) ... C(N+1) C(N+2) C(N+3) Min Max]
% 
% See also DA_LOG, DA_POLY

% YAADA - Software Toolkit to Analyze Single-Particle Mass Spectral Data
%
% Copyright (C) 1999-2000 California Institute of Technology
% Copyright (C) 2001-2008 Arizona State University
% Copyright (C) 2008 Jonathan O. Allen

% Jonathan O. Allen  02 Nov 99

% Changed to Velocity = (Min Max)
% JOA  2008-08-11

if nargin ~= 2
  error('Call as Da = DA_NOZ(Coef,Velocity)');
end
if iscell(Coef)
  Coef = double(Coef{:});
end
if min(size(Coef)) ~= 1
  error('Coefficients must be vector');
end
if max(size(Coef)) < 5 
  error('Coefficient vector must include exponential terms and limit');
end

Velocity = double(Velocity);

N = length(Coef)-5;

% reverse order of polynomial coeffs for polyval function
c2 = zeros(1,N+1);
for i = 1:N+1
  c2(i) = Coef(N+2-i);
end
Da = polyval(c2,Velocity);
Da = Da + Coef(N+2).*exp(Velocity*Coef(N+3));

OutsideLimit = find(Velocity < Coef(N+4) | Velocity > Coef(N+5));
Da(OutsideLimit) = NaN;

return
