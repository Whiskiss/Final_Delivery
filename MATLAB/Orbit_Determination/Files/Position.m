%--------------------------------------------------------------------------
%
% Position.m
%
% Purpose:
%   Position vector (r [m]) from geodetic coordinates (Longitude [rad],
%   latitude [rad], altitude [m])
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
function [r] = Position(lon, lat, h)

global const

R_equ = const.R_Earth;
f     = const.f_Earth;

e2     = f*(2.0-f);   % Square of eccentricity
CosLat = cos(lat);    % (Co)sine of geodetic latitude
SinLat = sin(lat);

% Position vector 
N = R_equ / sqrt(1.0-e2*SinLat*SinLat);

r(1) =  (         N+h)*CosLat*cos(lon);
r(2) =  (         N+h)*CosLat*sin(lon);
r(3) =  ((1.0-e2)*N+h)*SinLat;

