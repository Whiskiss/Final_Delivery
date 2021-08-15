%--------------------------------------------------------------------------
%
% Purpose:
%  Computes azimuth, elevation and partials from local tangent coordinates
%
% Input:
%   s      Topocentric local tangent coordinates (East-North-Zenith frame)
% 
% Outputs:
%   A      Azimuth [rad]
%   E      Elevation [rad]
%   dAds   Partials of azimuth w.r.t. s
%   dEds   Partials of elevation w.r.t. s
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
function [Az, El, dAds, dEds] = AzElPa(s)

pi2 = 2.0*pi;

rho = sqrt(s(1)*s(1)+s(2)*s(2));

% Angles
Az = atan2(s(1),s(2));

if (Az<0.0) 
    Az = Az+pi2;
end

El = atan ( s(3) / rho );

% Partials
dAds = [ s(2)/(rho*rho), -s(1)/(rho*rho), 0.0 ];
dEds = [ -s(1)*s(3)/rho, -s(2)*s(3)/rho , rho ] / dot(s,s);

