%--------------------------------------------------------------------------
%
% Purpose:
%   Computes the eccentric anomaly for elliptic orbits
%
% Inputs:
%   M         Mean anomaly in [rad]
%   e         Eccentricity of the orbit [0,1]
% 
% Output:
%             Eccentric anomaly in [rad]
%
% Last modified:   2015/08/12   M. Mahooti
% 
%--------------------------------------------------------------------------
function [E]  = EccAnom (M, e)

maxit = 15;
i = 1;

% Starting value
M = mod(M, 2.0*pi);

if (e<0.8)
    E = M; 
else
    E = pi;
end

f = E - e*sin(E) - M;
E = E - f / ( 1.0 - e*cos(E) );

% Iteration
while (abs(f) > 1e2*eps)    
    f = E - e*sin(E) - M;
    E = E - f / ( 1.0 - e*cos(E) );
    i = i+1;
    if (i==maxit)
        error(' convergence problems in EccAnom');
    end  
end

