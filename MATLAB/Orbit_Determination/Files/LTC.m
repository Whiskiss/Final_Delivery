%--------------------------------------------------------------------------
% 
% LTC.m
%
% Purpose:
%   Transformation from Greenwich meridian system to 
%   local tangent coordinates
%
% Inputs:
%   lon      -Geodetic East longitude [rad]
%   lat      -Geodetic latitude [rad]
%   
% Output:
%   M        -Rotation matrix from the Earth equator and Greenwich meridian
%             to the local tangent (East-North-Zenith) coordinate system
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
function M = LTC(lon, lat)

M = R_y(-1.0*lat)*R_z(lon);

for j=1:3
    Aux=M(1,j);
    M(1,j)=M(2,j);
    M(2,j)=M(3,j);
    M(3,j)= Aux;
end

end

