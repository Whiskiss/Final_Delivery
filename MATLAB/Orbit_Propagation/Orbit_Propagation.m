clc
clearvars -except dataMan data
%clear all
format long g
%close all
profile on
tic

global const Cnm Snm AuxParam eopdata swdata SOLdata DTCdata APdata PC

% Initial state for orbit propagation
strFile = 'InitialStateADS.txt';


SAT_Const
constants
load DE436Coeff.mat
PC = DE436Coeff;

% read Earth gravity field coefficients
Cnm = zeros(181,181);
Snm = zeros(181,181);
fid = fopen('GGM03S.txt','r');
for n=0:180
    for m=0:n
        temp = fscanf(fid,'%d %d %f %f %f %f',[6 1]);
        Cnm(n+1,m+1) = temp(3);
        Snm(n+1,m+1) = temp(4);
    end
end
fclose(fid);

% read Earth orientation parameters
fid = fopen('eop19620101.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------
eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
fclose(fid);

% read space weather data
fid = fopen('sw19571001.txt','r');
%  ---------------------------------------------------------------------------------------------------------------------------------
% |                                                                                             Adj     Adj   Adj   Obs   Obs   Obs 
% | yy mm dd BSRN ND Kp Kp Kp Kp Kp Kp Kp Kp Sum Ap  Ap  Ap  Ap  Ap  Ap  Ap  Ap  Avg Cp C9 ISN F10.7 Q Ctr81 Lst81 F10.7 Ctr81 Lst81
%  ---------------------------------------------------------------------------------------------------------------------------------
swdata = fscanf(fid,'%4i %3d %3d %5i %3i %3i %3i %3i %3i %3i %3i %3i %3i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4f %2i %4i %6f %2i %6f %6f %6f %6f %6f',[33 inf]);
fclose(fid);

% read solar storm indices
fid = fopen('SOLFSMY.txt','r');
%  ------------------------------------------------------------------------
% | YYYY DDD   JulianDay  F10   F81c  S10   S81c  M10   M81c  Y10   Y81c
%  ------------------------------------------------------------------------
SOLdata = fscanf(fid,'%d %d %f %f %f %f %f %f %f %f %f',[11 inf]);
fclose(fid);

% read geomagnetic storm indices
fid = fopen('DTCFILE.txt','r');
%  ------------------------------------------------------------------------
% | DTC YYYY DDD   DTC1 to DTC24
%  ------------------------------------------------------------------------
DTCdata = fscanf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d',[26 inf]);
fclose(fid);

% read geomagnetic storm indices
fid = fopen('SOLRESAP.txt','r');
%  ------------------------------------------------------------------------
% | YYYY DDD  F10 F10B Ap1 to Ap8
%  ------------------------------------------------------------------------
APdata = fscanf(fid,'%d %d %f %f %f %f %f %f %f %f %f',[12 inf]);
fclose(fid);

% model parameters
AuxParam = struct('Mjd_UTC',0,'area_solar',0,'area_drag',0,'mass',0,'Cr',0,...
                  'Cd',0,'n',0,'m',0,'sun',0,'moon',0,'sRad',0,'drag',0,...
                  'planets',0,'SolidEarthTides',0,'OceanTides',0,'Relativity',0);

% epoch state (Envisat)
fid = fopen(strFile,'r');
tline = fgetl(fid);
year = str2num(tline(1:4));
mon = str2num(tline(6:7));
day = str2num(tline(9:10));
hour = str2num(tline(12:13));
min = str2num(tline(15:16));
sec = str2num(tline(18:23));
Y0 = zeros(1,6);
for ii=1:6
    tline = fgetl(fid);
    Y0(ii) = str2num(tline);
end
tline = fgetl(fid);
AuxParam.area_solar = str2num(tline(49:end));
tline = fgetl(fid);
AuxParam.area_drag = str2num(tline(38:end));
tline = fgetl(fid);
AuxParam.mass = str2num(tline(19:end));
tline = fgetl(fid);
AuxParam.Cr = str2num(tline(5:end));
tline = fgetl(fid);
AuxParam.Cd = str2num(tline(5:end));
fclose(fid);

% If data was provided by Airbus DS, convert from ECEF to ECI reference
% frame (file name ending by 'ADS') _______________________________________
Mjd_UTC = Mjday(year, mon, day, hour, min, sec);
if strcmp(strFile,'InitialStateADS.txt') ~= 1
    Y0 = ECEF2ECI(Mjd_UTC, Y0);
end
% ________________________________________________________________________/

% Global vriables for the force model
AuxParam.Mjd_UTC = Mjd_UTC;
% Earth harmonics
AuxParam.n       = 10;
AuxParam.m       = 10;
% Sun gravitational force
AuxParam.sun     = 1;
% Moon gravitational force
AuxParam.moon    = 1;
% Solar system planets gravitational force
AuxParam.planets = 0;
% Solar radiation
AuxParam.sRad    = 0;
% Atmospheric drag
AuxParam.drag    = 0;
% Solid tides
AuxParam.SolidEarthTides = 0;
% Ocean tides
AuxParam.OceanTides = 0;
% Relativity
AuxParam.Relativity = 0;
% Maneuvers (type to be modified in Accel.m folder)
AuxParam.Maneuver = 0;
AuxParam.a = [];
AuxParam.t = [];
AuxParam.MJD = [];

Mjd0   = Mjd_UTC;
% Time step for propagation + total step number
Step   = 10;   % [s]
N_Step = 6000;

% Propagation _____________________________________________________________
Eph = Ephemeris(Y0, N_Step, Step);
% Eph = time(1), position(2:4), velocities(5:7)
%_________________________________________________________________________/

Y = Eph(:,2:4);

datetime(Mjd0,'ConvertFrom','modifiedjuliandate')
temp = datetime(year,mon,day,hour,min,sec+Eph(:,1));
obsTime = [temp.Year, temp.Month, temp.Day, temp.Hour, temp.Minute, temp.Second];


%% Creation of ground station observation measurements ____________________
% Kiruna Point station
lat = 67.8790708;  
lon = 21.038 ;
alt = 527.0;
Rs = [lat, lon, alt]';
% Azimuth/Elevation/Range
aerKiruna = eci2aer(Y, obsTime, repmat(Rs',size(Y,1),1));

% California Big Bear station
lat = 34.258; 
lon = -116.921; 
alt = 2067.0;
Rs = [lat, lon, alt]';
% Azimuth/Elevation/Range
aerBigBear = eci2aer(Y, obsTime, repmat(Rs',size(Y,1),1));
% ________________________________________________________________________/

%% Conversion to Keplerian elemnts ________________________________________
kep = [];
for idx =1:size(Eph,1)
	[~, a, e, i, Omega, omega, M] = elements(Eph(idx,2:7)) ;
    kep = [kep; [a, e, i, Omega, omega, M]];
end
% ________________________________________________________________________/

%% Creation of a data structure to save all information (manually) ________
%  Run simulation with and without maneuver, then save both files 'data'
%  and 'dataMan' manually. The loop to run both simulation has been
%  implemented yet.
if AuxParam.Maneuver
    idxKiruna = find(aerKiruna(:,2)>0);
    idxBigBear = find(aerBigBear(:,2)>0);
    dataMan = {Eph, AuxParam, aerKiruna, aerBigBear, ...
                [obsTime,aerKiruna(:,2),aerKiruna(:,1),aerKiruna(:,3)], ...
                [obsTime,aerBigBear(:,2),aerBigBear(:,1),aerBigBear(:,3)], ...
                [obsTime(idxKiruna,:),aerKiruna(idxKiruna,2),aerKiruna(idxKiruna,1),aerKiruna(idxKiruna,3)], ...
                [obsTime(idxBigBear,:),aerBigBear(idxBigBear,2),aerBigBear(idxBigBear,1),aerBigBear(idxBigBear,3)], ...
                kep}; 
    disp('Maneuver');
else
    idxKiruna = find(aerKiruna(:,2)>0);
    idxBigBear = find(aerBigBear(:,2)>0);
    data = {Eph, AuxParam, aerKiruna, aerBigBear, ...
                [obsTime,aerKiruna(:,2),aerKiruna(:,1),aerKiruna(:,3)], ...
                [obsTime,aerBigBear(:,2),aerBigBear(:,1),aerBigBear(:,3)], ...
                [obsTime(idxKiruna,:),aerKiruna(idxKiruna,2),aerKiruna(idxKiruna,1),aerKiruna(idxKiruna,3)], ...
                [obsTime(idxBigBear,:),aerBigBear(idxBigBear,2),aerBigBear(idxBigBear,1),aerBigBear(idxBigBear,3)], ...
                kep}; 
    disp('No Maneuver');
end
% ________________________________________________________________________/

disp('Eos');
