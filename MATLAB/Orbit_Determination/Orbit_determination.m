%%
% EKF for heterogenous measurement.

%%
clc
clearvars -except P
format long g
fclose('all') ;

%% Constant and variables
constantLoader;
global const Cnm Snm AuxParam eopdata n_eqn
%% load data structure .mat from the orbit propagation ( data and dataMan )
load HPOP_6.mat

%% Interval for start and end the EKF processing
iStart = 1;
iEnd = min(300,size(data{1,5},1)); 

% Noise(1) = bNoise 
% Noise(2:4) = signal to noise ratio for Azimuth Elvation Range
Noise = [1,50,50,50]; 

% Time iteration for initial orbit determination
timeIdx = [1,5,10];

% Force model
% Earth harmonics
AuxParam.n       = 20;
AuxParam.m       = 20;
% Sun gravitational force
AuxParam.sun     = 1;
% Moon gravitational force
AuxParam.moon    = 1;
% Solar system planets gravitational force
AuxParam.planets = 1; 

% Initial covariance matrix P0 - initialP(0): Position, initialP(1): Velocities 
initialP = [1e8, 1e3];

% Process noise - diagQ(0): Position, diagQ(1): Velocities 
diagQ = [0.1 0.1];

% Observation measurements to consider in the EKF _________________________
% Data type selection
% bEKF_AER(1:3): boolean to include Azimuth Elevation and Range
bEKF_AER = [1,1,1];
% Ground station selection
bStationKiruna = 1;
bStationBigBear = 0;
% ________________________________________________________________________/

% Boolean of Maneuver - 
% True:  maneuver included in the dataset (dataMan)
% False: No maneuver included in the dataset (data)
bManeuver = 0;

% Boolean to consider only realistic measurement (Elevation >0)
bVisible = 0;

% Option to configure for getCovariance.m
iChoice = 5;

% Definition of a fixed measurement noise variance
sigma_range = 1000;
sigma_az = 0.001*const.Rad; % [rad]
sigma_el = 0.001*const.Rad; % [rad]

K_all = []; % Kalman gain
Cov_all = []; % Variance of the covariance matrix
kep = []; % Keplerian elements
Y_all = []; % Cartesian variables (P,V)


% Do not change
ObsGap = 1;
idxTimeGap = Inf;

% _______________________________________
% data structure for data recording (Debug purpose only)
data_spec.iStart = iStart;
data_spec.iEnd = iEnd;
data_spec.Noise = Noise;
data_spec.timeIdx = timeIdx;
data_spec.AuxParam = AuxParam;
data_spec.initialP = initialP;
data_spec.diagQ = diagQ;
data_spec.bEKF_AER = bEKF_AER;
data_spec.bStationKiruna = bStationKiruna;
data_spec.bStationBigBear = bStationBigBear;
data_spec.bManeuver = bManeuver;


%%

dataAll = {};
for i=5:6
    %No Maneuvers from Kiruna and Big Bear
	dataAll = [dataAll data{1,i}(iStart:iEnd,:)];
end
if exist('dataMan','var')
    %Maneuvers from Kiruna and Big Bear
    for i=5:6
        dataAll = [dataAll dataMan{1,i}(iStart:iEnd,:)];
    end
end


dataVisible = {};
for i=7:8
    %No Maneuvers from Kiruna and Big Bear
    dataVisible = [dataVisible data{1,i}];
end
if exist('dataMan','var')
    %Maneuvers from Kiruna and Big Bear
    for i=7:8
        dataVisible = [dataVisible dataMan{1,i}];
    end
end

%% Noise addition
if Noise(1)
    snrAER = Noise(2:4);  
    for i=1:size(dataAll,2)
        % dataAll -> El Az Ra
        dataAll{1,i}(:,7) = awgn(dataAll{1,i}(:,7),snrAER(2));
        dataAll{1,i}(:,8) = awgn(dataAll{1,i}(:,8),snrAER(1));  
        dataAll{1,i}(:,9) = awgn(dataAll{1,i}(:,9),snrAER(3)); 
    end
    for i=1:size(dataVisible,2)
        % dataAll -> El Az Ra
        dataVisible{1,i}(:,7) = awgn(dataVisible{1,i}(:,7),snrAER(2));
        dataVisible{1,i}(:,8) = awgn(dataVisible{1,i}(:,8),snrAER(1));  
        dataVisible{1,i}(:,9) = awgn(dataVisible{1,i}(:,9),snrAER(3)); 
    end      
end
%% Station loading: location
stationLoader;

%% Selection of the data to be fed to the EKF according visibility (elevation>0)
if bVisible
    dataEKF = dataVisible;
    for i=1:size(dataEKF,2)
        dataEKF{1,i} = dataEKF{1,i}(1:min(iEnd,size(dataEKF{1,i},1)),:);
    end
else
    dataEKF = dataAll;
end




% Data loading for the EKF ________________________________________________
iData = 1;
if bManeuver
    iData = 3;
end
if bStationKiruna
    dataEKF{2,iData} = []; %K
    dataEKF{3,iData} = []; %Var
    dataEKF{4,iData} = []; %Kep
    dataEKF{5,iData} = []; %Y
    dataEKF{6,iData} = []; %Residuals Az
    dataEKF{7,iData} = []; %Residuals El
    dataEKF{8,iData} = []; %Residuals Ra
end
if bStationBigBear
    dataEKF{2,iData+1} = []; %K
    dataEKF{3,iData+1} = []; %Var
    dataEKF{4,iData+1} = []; %Kep
    dataEKF{5,iData+1} = []; %Y
    dataEKF{6,iData+1} = []; %Residuals Az
    dataEKF{7,iData+1} = []; %Residuals El
    dataEKF{8,iData+1} = []; %Residuals Ra
end
    
nobs = size(dataEKF{1,iData},1);
obsSK = zeros(nobs,4);
obsCB = zeros(nobs,4);

% read observations
for idx=1:size(dataEKF{1,iData},1)
    Y = dataEKF{1,iData}(idx,1);
    M = dataEKF{1,iData}(idx,2);
    D = dataEKF{1,iData}(idx,3);
    h = dataEKF{1,iData}(idx,4);
    m = dataEKF{1,iData}(idx,5);
    s = dataEKF{1,iData}(idx,6);
    
    % Sweden Kiruna
    azSK = dataEKF{1,iData}(idx,8);
    elSK = dataEKF{1,iData}(idx,7);
    raSK = dataEKF{1,iData}(idx,9);
    
    obsSK(idx,1) = Mjday(Y,M,D,h,m,s);
    obsSK(idx,2) = const.Rad*azSK;
    obsSK(idx,3) = const.Rad*elSK;
    obsSK(idx,4) = raSK;    
    
    if ~bVisible
        % California Big Bear
        azCB = dataEKF{1,iData+1}(idx,8);
        elCB = dataEKF{1,iData+1}(idx,7);
        raCB = dataEKF{1,iData+1}(idx,9); 
        
        obsCB(idx,1) = Mjday(Y,M,D,h,m,s);
        obsCB(idx,2) = const.Rad*azCB;
        obsCB(idx,3) = const.Rad*elCB;
        obsCB(idx,4) = raCB;          
    end
end
% ________________________________________________________________________/

%% Extended Kalman Filter iterations


% Computation of the Initial Orbit Determination __________________________
if bStationKiruna
    % IOD using Gauss Gibbs
    Mjd1 = obsSK(timeIdx(1),1);
    Mjd2 = obsSK(timeIdx(2),1);
    Mjd3 = obsSK(timeIdx(3),1);
    [r2,v2] = anglesg(obsSK(timeIdx(1),2),obsSK(timeIdx(2),2),obsSK(timeIdx(3),2),obsSK(timeIdx(1),3),obsSK(timeIdx(2),3),obsSK(timeIdx(3),3),...
                     Mjd1,Mjd2,Mjd3,Rs{1,1},Rs{1,1},Rs{1,1});
    LTSK = LTC(Rs{2,1}(1),Rs{2,1}(2)); %lon lat
end
if bStationBigBear
    % IOD using Gauss & Gibbs
    Mjd1 = obsCB(timeIdx(1),1);
    Mjd2 = obsCB(timeIdx(2),1);
    Mjd3 = obsCB(timeIdx(3),1);
    [r2,v2] = anglesg(obsCB(timeIdx(1),2), obsCB(timeIdx(2),2),obsCB(timeIdx(3),2),obsCB(timeIdx(1),3),obsCB(timeIdx(2),3),obsCB(timeIdx(3),3),...
                 Mjd1,Mjd2,Mjd3,Rs{1,2},Rs{1,2},Rs{1,2}); 
    LTCB = LTC(Rs{2,2}(1),Rs{2,2}(2)); %lon lat
end
     
Y0_apr = [r2;v2];
Mjd0 = Mjday(dataEKF{1,iData}(1,1),dataEKF{1,iData}(1,2),dataEKF{1,iData}(1,3),dataEKF{1,iData}(1,4),dataEKF{1,iData}(1,5),dataEKF{1,iData}(1,6));
% ________________________________________________________________________/

Mjd_UTC = obsSK(timeIdx(2),1);
AuxParam.Mjd_UTC = Mjd_UTC;
n_eqn  = 6;
% Pre-computation
Y = DEInteg(@Accel,0,-(obsSK(timeIdx(2),1)-Mjd0)*86400.0,1e-13,1e-6,6,Y0_apr);

% Initial covariance matrix
P = zeros(6);
for idx=1:3
    P(idx,idx)=initialP(1); %1e8;
end
for idx=4:6
    P(idx,idx)=initialP(2); %1e3;
end

yPhi = zeros(42,1);
Phi  = zeros(6);
% Measurement loop
t = 0;
disp(sprintf('Total observations: %i',nobs));

% Data Measurement iterations
for idx=1:ObsGap:min(10000,nobs)
    % Previous step
    t_old = t;
    Y_old = Y;

    % Computation of the new epoch
    Mjd_UTC = obsSK(idx,1)  ;               % Modified Julian Date
    t = (Mjd_UTC-Mjd0)*86400.0;        	% Time since epoch [s]
    
    [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,Mjd_UTC,'l');
    [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_TT;

    % PREDICTION __________________________________________________________
    
    for ii=1:6
        yPhi(ii) = Y_old(ii);
        for j=1:6  
            if (ii==j) 
                yPhi(6*j+ii) = 1; 
            else
                yPhi(6*j+ii) = 0;
            end
        end
    end
    yPhi = DEInteg (@VarEqn,0,t-t_old,1e-13,1e-6,42,yPhi);

    % Extract state transition matrices
    for j=1:6
        Phi(:,j) = yPhi(6*j+1:6*j+6);
    end
    Y = DEInteg(@Accel,0,t-t_old,1e-13,1e-6,6,Y_old);

    % Compute covariance matrix
    getCovariance;

    % ESTIMATION __________________________________________________________
    
    % Computation of the Earth orientation at epoch Mjd_UT1
    theta = gmst(Mjd_UT1);                   
    U = R_z(theta); 
    
    if idx == 1, dataEKF{2,iData}=[0;0;0;0;0;0]; end; % First initeration
    % Computation if Kiruna station is available
    if bStationKiruna
        if bEKF_AER(1) % Azimuth
            % Topocentric coordinates
            r = Y(1:3);
            s = LTSK*(U*r-Rs{1,1});     % ENZ frame     % Topocentric position [m]
            % Azimuth and partials
            [Azim, Elev, dAds, dEds] = AzElPa(s);     % Azimuth, Elevation
            dAdY = [dAds*LTSK*U,zeros(1,3)];
            % Measurement update
            [K, Y, P] = MeasUpdate ( Y, obsSK(idx,2), Azim, sigma_az, dAdY, P, 6); 
            dataEKF{6,iData} = [dataEKF{6,iData} obsSK(idx,2)-Azim ];
        end
        if bEKF_AER(2) % Elevation
            % Elevation and partials
            r = Y(1:3);
            s = LTSK*(U*r-Rs{1,1});                     % Topocentric position [m]
            [Azim, Elev, dAds, dEds] = AzElPa(s);     % Azimuth, Elevation
            dEdY = [dEds*LTSK*U,zeros(1,3)];
            % Measurement update
            [K, Y, P] = MeasUpdate ( Y, obsSK(idx,3), Elev, sigma_el, dEdY, P, 6); 
            dataEKF{7,iData} = [dataEKF{7,iData} obsSK(idx,3)-Elev ];
        end
        if bEKF_AER(3) % Range
            % Range and partials
             r = Y(1:3);
             s = LTSK*(U*r-Rs{1,1});                          % Topocentric position [m]
             Dist = norm(s);
             dDds = (s/Dist)';         % Range
             dDdY = [dDds*LTSK*U,zeros(1,3)];
             % Measurement update
             [K, Y, P] = MeasUpdate ( Y, obsSK(idx,4), Dist, sigma_range, dDdY, P, 6 );
             dataEKF{8,iData} = [dataEKF{8,iData} obsSK(idx,4)-Dist ];
        end
    end
    % Computation if BigBear station is available
    if bStationBigBear
        if bEKF_AER(1) % Azimuth
            % Topocentric coordinates
            r = Y(1:3);
            s = LTCB*(U*r-Rs{1,2});     % ENZ frame     % Topocentric position [m]
            % Azimuth and partials
            [Azim, Elev, dAds, dEds] = AzElPa(s);     % Azimuth, Elevation
            dAdY = [dAds*LTCB*U,zeros(1,3)];
            % Measurement update
            [K, Y, P] = MeasUpdate ( Y, obsCB(idx,2), Azim, sigma_az, dAdY, P, 6 );
            dataEKF{6,iData+1} = [dataEKF{6,iData+1} obsCB(idx,2)-Azim ];
        end
        if bEKF_AER(2) % Elevation
            % Elevation and partials
            r = Y(1:3);
            s = LTCB*(U*r-Rs{1,2});                     % Topocentric position [m]
            [Azim, Elev, dAds, dEds] = AzElPa(s);     % Azimuth, Elevation
            dEdY = [dEds*LTCB*U,zeros(1,3)];
            % Measurement update
            [K, Y, P] = MeasUpdate ( Y, obsCB(idx,3), Elev, sigma_el, dEdY, P, 6); 
            dataEKF{7,iData+1} = [dataEKF{7,iData+1} obsCB(idx,3)-Elev ];
        end
        if bEKF_AER(3) % Range
            % Range and partials
             r = Y(1:3);
             s = LTCB*(U*r-Rs{1,2});                          % Topocentric position [m]
             Dist = norm(s); dDds = (s/Dist)';         % Range
             dDdY = [dDds*LTCB*U,zeros(1,3)];
             % Measurement update
             [K, Y, P] = MeasUpdate ( Y, obsCB(idx,4), Dist, sigma_range, dDdY, P, 6);
             dataEKF{8,iData+1} = [dataEKF{8,iData+1} obsCB(idx,4)-Dist ];
        end        
    end
    
    % Data recording
    if bStationKiruna
        dataEKF{2,iData} = [dataEKF{2,iData} K];
        dataEKF{3,iData} = [dataEKF{3,iData} [P(1,1); P(2,2); P(3,3); P(4,4); P(5,5);P(6,6)]];
        dataEKF{5,iData} = [dataEKF{5,iData} Y];
        [~, a, e, i, Omega, omega, M] = elements(Y);
        dataEKF{4,iData} = [dataEKF{4,iData} [a,e,i,Omega,omega,M]'];
    end
    if bStationBigBear
        dataEKF{2,iData+1} = [dataEKF{2,iData+1} K];
        dataEKF{3,iData+1} = [dataEKF{3,iData+1} [P(1,1); P(2,2); P(3,3); P(4,4); P(5,5);P(6,6)]];
        dataEKF{5,iData+1} = [dataEKF{5,iData+1} Y];
        [~, a, e, i, Omega, omega, M] = elements(Y);
        dataEKF{4,iData+1} = [dataEKF{4,iData+1} [a,e,i,Omega,omega,M]'];     
    end    
    
    % Display idx value for visual follow up
    if mod(idx,25)==0
        disp(sprintf('Number of operation processed: %i', idx));
    end
end

%% Plot

if bManeuver
    inputData = dataMan;
else
    inputData = data;
end


if bVisible
    idx = find(inputData{1,5}(:,7)>0);
    idx = idx(1:min(iEnd,size(dataEKF{1,iData},1)));
else
    idx = [iStart:iEnd];
end


if bStationKiruna 
    figure;
    text = {'K_x','K_y','K_z','K_vx','K_vy','K_vz'};
    for i=1:6
        subplot(6,1,i);
        plot(dataEKF{2,iData}(i,:));
        grid on;
        legend(text{i});
    end
    subplot(6,1,1); title('Kalman gain from Kiruna');    
    subplot(6,1,6); xlabel('time[min]'); 
    
    figure;
    text = {'\sigma^2_x','\sigma^2_y','\sigma^2_z','\sigma^2_vx','\sigma^2_vy','\sigma^2_vz'};
    for i=1:6
        subplot(6,1,i);
        plot(dataEKF{3,iData}(i,:));
        grid on;
        legend(text{i});
    end
    subplot(6,1,1); title('Variance from Kiruna');        
    subplot(6,1,6); xlabel('time[min]');
      
    figure;
    text = {'x[m]','y[m]','z[m]','vx[m/s]','vy[m/s]','vz[m/s]'};
    for i=1:6
        subplot(6,1,i);
        plot(inputData{1,1}(idx,i+1));
        hold on;
        plot(dataEKF{5,iData}(i,:)');
        grid on;
        legend(text{i});
    end
    subplot(6,1,1);
    subplot(6,1,1); title('PV from Kiruna');  
    subplot(6,1,6); xlabel('time[min]');
    
    errX = inputData{1,1}(idx,2) - dataEKF{5,iData}(1,:)';
    errY = inputData{1,1}(idx,3) - dataEKF{5,iData}(2,:)';
    errZ = inputData{1,1}(idx,4) - dataEKF{5,iData}(3,:)';
    err = sqrt(errX.^2 + errY.^2 + errZ.^2);
    figure;
    subplot(3,1,1);
    plot(errX./1000);
    xlabel('time[min]');ylabel('error X [km]');grid on;
    subplot(3,1,2);
    plot(errY./1000);
    xlabel('time[min]');ylabel('error Y[km]');grid on;
    subplot(3,1,3);
    plot(errZ./1000);
    xlabel('time[min]');ylabel('error Z [km]');grid on;
    subplot(3,1,1);
    title('Position error XYZ from Kiruna');    
    
    figure;
    plot(err./1000);
    grid on;
    title('Position error from Kiruna');
    xlabel('time[min]');
    ylabel('distance[km]');
    
    figure;
    text = {'sma[m]','e','i[rad]','Omega[rad]','omega[rad]','ma[rad]'};
    for i=1:6
        subplot(6,1,i);
        plot(dataEKF{4,iData}(i,:));
        hold on;
        plot(inputData{1,9}(idx,i));
        grid on;
        legend(text{i});
    end
    subplot(6,1,1);
    subplot(6,1,1); title('Keplerian Elements from Kiruna');
    subplot(6,1,6); xlabel('time[min]');

    if bEKF_AER(1)
        figure;
        plot(dataEKF{6,iData}(1,:));
        grid on;
        title('Az Residuals from Kiruna');
        xlabel('time[min]');
        ylabel('angle[rad]');
    end
    if bEKF_AER(2)
        figure;
        plot(dataEKF{7,iData}(1,:));
        grid on;
        title('El Residuals from Kiruna');
        xlabel('time[min]');
        ylabel('angle[rad]');
    end
    if bEKF_AER(3)
        figure;
        plot(dataEKF{8,iData}(1,:));
        grid on;
        title('Ra Residuals from Kiruna');
        xlabel('time[min]');
        ylabel('distance[m]');
    end   
end


if bStationBigBear
    
    figure;
    text = {'K_x','K_y','K_z','K_vx','K_vy','K_vz'};
    for i=1:6
        subplot(6,1,i);
        plot(dataEKF{2,iData+1}(i,:));
        grid on;
        legend(text{i});
    end
    subplot(6,1,1); title('Kalman gain from Kiruna');    
    subplot(6,1,6); xlabel('time[min]');
     
    figure;
    text = {'\sigma^2_x','\sigma^2_y','\sigma^2_z','\sigma^2_vx','\sigma^2_vy','\sigma^2_vz'};
    for i=1:6
        subplot(6,1,i);
        plot(dataEKF{3,iData+1}(i,:));
        grid on;
        legend(text{i});
    end
    subplot(6,1,1); title('Variance from Kiruna');        
    subplot(6,1,6); xlabel('time[min]');
     
    figure;
    text = {'x[m]','y[m]','z[m]','vx[m/s]','vy[m/s]','vz[m/s]'};
    for i=1:6
        subplot(6,1,i);
        plot(inputData{1,1}(idx,i+1));
        hold on;
        plot(dataEKF{5,iData+1}(i,:)');
        grid on;
        legend(text{i});
    end
    subplot(6,1,1); title('PV from Big Bear');  
    subplot(6,1,6); xlabel('time[min]');
    
    errX = inputData{1,1}(idx,2) - dataEKF{5,iData+1}(1,:)';
    errY = inputData{1,1}(idx,3) - dataEKF{5,iData+1}(2,:)';
    errZ = inputData{1,1}(idx,4) - dataEKF{5,iData+1}(3,:)';
    err = sqrt(errX.^2 + errY.^2 + errZ.^2);
    figure;
    plot(err./1000);
    grid on;
    title('Position error from Big Bear');
    xlabel('time[min]');
    ylabel('distance[km]');
    figure;
    subplot(3,1,1);
    plot(errX./1000);
    xlabel('time[min]');ylabel('error X [km]');grid on;
    subplot(3,1,2);
    plot(errY./1000);
    xlabel('time[min]');ylabel('error Y[km]');grid on;
    subplot(3,1,3);
    plot(errZ./1000);
    xlabel('time[min]');ylabel('error Z [km]');grid on;
    subplot(3,1,1);
    title('Position error XYZ from Big Bear');
     
    figure;
    text = {'sma[m]','e','i[rad]','Omega[rad]','omega[rad]','ma[rad]'};
    for i=1:6
        subplot(6,1,i);
        plot(dataEKF{4,iData+1}(i,:));
        hold on;
        plot(inputData{1,9}(idx,i));
        grid on;
        legend(text{i});
    end
    subplot(6,1,1); title('Keplerian Elements from Big Bear');
    subplot(6,1,6); xlabel('time[min]');

    if bEKF_AER(1)
        figure;
        plot(dataEKF{6,iData+1}(1,:));
        grid on;
        title('Az Residuals from Big Bear');
        xlabel('time[min]');
        ylabel('angle[rad]');
    end
    if bEKF_AER(2)
        figure;
        plot(dataEKF{7,iData+1}(1,:));
        grid on;
        title('El Residuals from Big Bear');
        xlabel('time[min]');
        ylabel('angle[rad]');
    end
    if bEKF_AER(3)
        figure;
        plot(dataEKF{8,iData+1}(1,:));
        grid on;
        title('Ra Residuals from Big Bear');
        xlabel('time[min]');
        ylabel('distance[m]');
    end
end

%%
disp('EoS')
return;

%% Computation and display of mean and differential values
windowSize = 10;

meanVal = zeros(size(dataEKF{4,iData},1), size(dataEKF{4,iData},2)-windowSize);
diffVal = [];
for i = 1:size(dataEKF{4,iData},1) 
    for j = 1:size(dataEKF{4,iData},2)-windowSize
        meanVal(i,j) = mean(dataEKF{4,iData}(i,j:j+windowSize));

    end
    diffVal = [diffVal; diff(dataEKF{4,iData}(i,:))];
end

diffRes = [];
for i = 6:size(dataEKF(:,iData),1)
    if ~isempty(dataEKF{i,iData})
        diffRes = [diffRes; diff(dataEKF{i,iData}(:))'];
    end
end

figure;
text = {'sma[m]','e','i[rad]','Omega[rad]','omega[rad]','ma[rad]'};
for i=1:6
    subplot(6,1,i);
    plot(meanVal(i,:));
    grid on;
    legend(text{i});
end
subplot(6,1,1); title('Mean Keplerian Elements from Kiruna');
subplot(6,1,6); xlabel('time[min]');

figure;
text = {'sma[m]','e','i[rad]','Omega[rad]','omega[rad]','ma[rad]'};
for i=1:6
    subplot(6,1,i);
    plot(diffVal(i,:));
    grid on;
    legend(text{i});
end
subplot(6,1,1); title('Diff Keplerian Elements from Kiruna');
subplot(6,1,6); xlabel('time[min]');

figure;
text = {'Residual Az[rad]','Residual El[rad]','Residual Ra[m]'};
for i=1:size(diffRes,1)
    subplot(size(diffRes,1),1,i);
    plot(diffRes(i,:));
    grid on;
    legend(text{i});
end
subplot(size(diffRes,1),1,1); title('Diff of residuals');
subplot(size(diffRes,1),1,size(diffRes,1)); xlabel('time[min]');

%%
disp('EoS')
return;
