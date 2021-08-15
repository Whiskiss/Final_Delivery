%--------------------------------------------------------------------------
%
% G_AccelHarmonic.m
%
% Purpose:
%   Computes the gradient of the Earth's harmonic gravity field 
%
% Inputs:
%   r           Satellite position vector in the true-of-date system
%   U           Transformation matrix to body-fixed system
%   n           Gravity model degree
%   m 			Gravity model order
%
% Output:
%   G    		Gradient (G=da/dr) in the true-of-date system
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
function [G] = G_AccelHarmonic( r, U, n_max, m_max )

d = 1.0;   % Position increment [m]

G = zeros(3);
dr = zeros(3,1);

% Gradient
for i=1:3
    % Set offset in i-th component of the position vector
    dr(:) = 0.0;
    dr(i) = d;
    % Acceleration difference
    da = AccelHarmonic ( r+dr/2,U, n_max, m_max ) - ...
         AccelHarmonic ( r-dr/2,U, n_max, m_max );
    % Derivative with respect to i-th axis
    G(:,i) = da/d;      
end

