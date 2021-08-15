%--------------------------------------------------------------------------
%
% ECI2ECEF: Transforms Earth Centered Inertial (ECI) coordinates to Earth
%           Centered Earth Fixed (ECEF) coordinates
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------
function Y = ECI2ECEF(MJD_UTC, Y0)

global const eopdata

[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,MJD_UTC,'l');
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
[year, month, day, hour, minute, sec] = invjday(MJD_UTC+2400000.5);
[DJMJD0, DATE] = iauCal2jd(year, month, day);
TIME = (60*(60*hour + minute) + sec)/86400;
UTC = DATE + TIME;
TT = UTC + TT_UTC/86400;
TUT = TIME + UT1_UTC/86400;
UT1 = DATE + TUT;

% Form bias-precession-nutation matrix
NPB = iauPnm06a(DJMJD0, TT);
% Form Earth rotation matrix
Theta  = iauRz( iauGst06(DJMJD0, UT1, DJMJD0, TT, NPB),eye(3) );
% Polar motion matrix (TIRS->ITRS, IERS 2003)
Pi = iauPom00(x_pole, y_pole, iauSp00(DJMJD0, TT));

% ICRS to ITRS transformation matrix and derivative
S = zeros(3);
S(1,2) = 1; S(2,1) = -1;              			 % Derivative of Earth rotation
Omega  = 7292115.8553e-11+4.3e-15*( (MJD_UTC-const.MJD_J2000)/36525 ); % [rad/s]
% Omega = const.omega_Earth-0.843994809*1e-9*LOD;  % IERS
dTheta = Omega*S*Theta;           				 % matrix [1/s]
U      = Pi*Theta*NPB;                			 % ICRS to ITRS transformation
dU     = Pi*dTheta*NPB;               			 % Derivative [1/s]

% Transformation from ICRS to WGS
r = U*Y0(1:3)';
v = U*Y0(4:6)' + dU*Y0(1:3)';
Y = [r;v];

