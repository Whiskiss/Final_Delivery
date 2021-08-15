%--------------------------------------------------------------------------
%
% Purpose:
%   Computes the osculating Keplerian elements from the satellite state
%   vector for elliptic orbits
%
% Input:
%    y        State vector (x,y,z,vx,vy,vz)
% 
% Outputs:
%    p        semilatus rectum [m]
%    a        Semimajor axis 
%    e        Eccentricity 
%    i        Inclination [rad]
%    Omega    Longitude of the ascending node [rad]
%    omega    Argument of pericenter [rad]
%    M        Mean anomaly [rad]
%
% Notes:
%   The function cannot be used with state vectors describing a circular
%   or non-inclined orbit.
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
 function [p, a, e, i, Omega, omega, M] = elements(y)
 
 global const
 
 pi2 = 2*pi;
 
 r = y(1:3);                                        % Position
 v = y(4:6);                                        % Velocity
 
 h = cross(r,v);                                    % Areal velocity
 magh = norm(h);
 p = magh*magh/const.GM_Earth;
 H = norm(h);
 
 Omega = atan2 ( h(1), -h(2) );                     % Long. ascend. node 
 Omega = mod(Omega,pi2);
 i     = atan2 ( sqrt(h(1)*h(1)+h(2)*h(2)), h(3) ); % Inclination        
 u     = atan2 ( r(3)*H, -r(1)*h(2)+r(2)*h(1) );    % Arg. of latitude   
 
 R  = norm(r);                                      % Distance           
 
 a = 1/(2/R-dot(v,v)/const.GM_Earth);               % Semi-major axis    
 
 eCosE = 1-R/a;                                     % e*cos(E)           
 eSinE = dot(r,v)/sqrt(const.GM_Earth*a);           % e*sin(E)           
 % test
     eSinE = real(eSinE);
     eCosE = real(eCosE);
 % END test 
 e2 = eCosE*eCosE +eSinE*eSinE;
 e  = sqrt(e2);                                     % Eccentricity 

 E  = atan2(eSinE,eCosE);                           % Eccentric anomaly  
 
 M  = mod(E-eSinE,pi2);                             % Mean anomaly
 
 nu = atan2(sqrt(1.0-e2)*eSinE, eCosE-e2);          % True anomaly
 
 omega = mod(u-nu,pi2);                             % Arg. of perihelion 
 
 