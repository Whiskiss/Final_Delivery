clc
clear all
format long g
fclose('all') 
close all

%% Constant and variables
constantLoader;
global const Cnm Snm AuxParam eopdata n_eqn
%%
load HPOP_ADS_FineMan.mat
tManeuver = 83; % DO NOT CHANGE
tManeuver = 167; % DO NOT CHANGE
tManeuver = 1000; % DO NOT CHANGE

iLen = 35;
iStart = 950; %50
iEnd = min(1050,size(data{1,4},1)); 
% Noise(1) = bNoise / Noise(2:4) = snrAER
Noise = [0,10,10,10]; 
timeIdx = [1,5,10];

AuxParam.n       = 1; %20
AuxParam.m       = 1; %20
AuxParam.sun     = 1; %1
AuxParam.moon    = 1; %1
AuxParam.planets = 0; %1

initialP = [1e8, 1e3];
diagQ = [0.1 0.1];

% Observation to consider in the EKF
bEKF_AER = [1,1,0];
bStationKiruna = 1;
bStationBigBear = 0;

bManeuver = 1;

idxTimeGap = Inf;
% _______________________________________
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
	dataAll = [dataAll data{1,i}];
end
if exist('dataMan','var')
    %Maneuvers from Kiruna and Big Bear
    for i=5:6
        dataAll = [dataAll dataMan{1,i}];
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


% % Noise added
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

stationLoader;

%%

for iData = 1:2:size(dataAll,2)
    for jData = iStart:iEnd
        clc
        disp(sprintf('%i %i/%i',(iData+1)/2, jData, iEnd));
        
        K_all = [];
        Cov_all = [];
        kep = [];
        Y_all = [];

        if bStationKiruna
            dataAll{2,iData} = []; %K
            dataAll{3,iData} = []; %Var
            dataAll{4,iData} = []; %Kep
            dataAll{5,iData} = []; %Y
            dataAll{6,iData} = []; %Residuals Az
            dataAll{7,iData} = []; %Residuals El
            dataAll{8,iData} = []; %Residuals Ra
        end
        if bStationBigBear
            dataAll{2,iData+1} = []; %K
            dataAll{3,iData+1} = []; %Var
            dataAll{4,iData+1} = []; %Kep
            dataAll{5,iData+1} = []; %Y
            dataAll{6,iData+1} = []; %Residuals Az
            dataAll{7,iData+1} = []; %Residuals El
            dataAll{8,iData+1} = []; %Residuals Ra
        end

        nobs = size(dataAll{1,iData},1);
        obsSK = zeros(nobs,4);
        obsCB = zeros(nobs,4);

        % read observations
        for idx=1:size(dataAll{1,iData},1)
            Y = dataAll{1,iData}(idx,1);
            M = dataAll{1,iData}(idx,2);
            D = dataAll{1,iData}(idx,3);
            h = dataAll{1,iData}(idx,4);
            m = dataAll{1,iData}(idx,5);
            s = dataAll{1,iData}(idx,6);

            % Sweden Kiruna
            azSK = dataAll{1,iData}(idx,8);
            elSK = dataAll{1,iData}(idx,7);
            raSK = dataAll{1,iData}(idx,9);

            % California Big Bear
            azCB = dataAll{1,iData+1}(idx,8);
            elCB = dataAll{1,iData+1}(idx,7);
            raCB = dataAll{1,iData+1}(idx,9);    

            obsSK(idx,1) = Mjday(Y,M,D,h,m,s);
            obsSK(idx,2) = const.Rad*azSK;
            obsSK(idx,3) = const.Rad*elSK;
            obsSK(idx,4) = raSK;

            obsCB(idx,1) = Mjday(Y,M,D,h,m,s);
            obsCB(idx,2) = const.Rad*azCB;
            obsCB(idx,3) = const.Rad*elCB;
            obsCB(idx,4) = raCB;    

        end

        sigma_range = 92.5;
        sigma_az = 0.0224*const.Rad; % [rad]
        sigma_el = 0.0139*const.Rad; % [rad]

        if bStationKiruna
            % IOD using Gibbs
            Mjd1 = obsSK(timeIdx(1)+jData-1,1);
            Mjd2 = obsSK(timeIdx(2)+jData-1,1);
            Mjd3 = obsSK(timeIdx(3)+jData-1,1);
            [r2,v2] = anglesg(obsSK(timeIdx(1)+jData-1,2),obsSK(timeIdx(2)+jData-1,2),obsSK(timeIdx(3)+jData-1,2), ...
                              obsSK(timeIdx(1)+jData-1,3),obsSK(timeIdx(2)+jData-1,3),obsSK(timeIdx(3)+jData-1,3),...
                             Mjd1,Mjd2,Mjd3,Rs{1,1},Rs{1,1},Rs{1,1});
            LTSK = LTC(Rs{2,1}(1),Rs{2,1}(2)); %lon lat
        end
        if bStationBigBear
            % IOD using Gibbs
            Mjd1 = obsCB(timeIdx(1)+jData-1,1);
            Mjd2 = obsCB(timeIdx(2)+jData-1,1);
            Mjd3 = obsCB(timeIdx(3)+jData-1,1);
            [r2,v2] = anglesg(obsCB(timeIdx(1)+jData-1,2), obsCB(timeIdx(2)+jData-1,2),obsCB(timeIdx(3)+jData-1,2), ...
                              obsCB(timeIdx(1)+jData-1,3),obsCB(timeIdx(2)+jData-1,3),obsCB(timeIdx(3)+jData-1,3),...
                         Mjd1,Mjd2,Mjd3,Rs{1,2},Rs{1,2},Rs{1,2}); 
            LTCB = LTC(Rs{2,2}(1),Rs{2,2}(2)); %lon lat
        end

        Y0_apr = [r2;v2];
        Mjd0 = Mjday(dataAll{1,iData}(1,1),dataAll{1,iData}(1,2),dataAll{1,iData}(1,3),dataAll{1,iData}(1,4),dataAll{1,iData}(1,5),dataAll{1,iData}(1,6));


        Mjd_UTC = obsSK(timeIdx(2)+jData-1,1);

        AuxParam.Mjd_UTC = Mjd_UTC;
        % AuxParam.n       = 1; %20
        % AuxParam.m       = 1; %20
        % AuxParam.sun     = 0; %1
        % AuxParam.moon    = 0; %1
        % AuxParam.planets = 0; %1

        n_eqn  = 6;
        Y = DEInteg(@Accel,0,-(obsSK(timeIdx(2)+jData-1,1)-Mjd0)*86400.0,1e-13,1e-6,6,Y0_apr);

        % Covariance matrix initialization
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
        disp(sprintf('%i samples to be processed by EKF', iLen));

        for idx=jData:1:min(jData+iLen-1,nobs)
            % Previous step
            t_old = t;
            Y_old = Y;

            % Time increment and propagation
            Mjd_UTC = obsSK(idx,1)  ;               % Modified Julian Date
            t = (Mjd_UTC-Mjd0)*86400.0;        	% Time since epoch [s]

            [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,Mjd_UTC,'l');
            [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
            Mjd_TT = Mjd_UTC + TT_UTC/86400;
            Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
            AuxParam.Mjd_UTC = Mjd_UTC;
            AuxParam.Mjd_TT = Mjd_TT;

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

            iChoice = 4;
            getCovariance;

            theta = gmst(Mjd_UT1);                    % Earth rotation
            U = R_z(theta); 

            if bStationKiruna
                if bEKF_AER(1)
                    % Topocentric coordinates
                    r = Y(1:3);
                    s = LTSK*(U*r-Rs{1,1});     % ENZ frame     % Topocentric position [m]
                    % Azimuth and partials
                    [Azim, Elev, dAds, dEds] = AzElPa(s);     % Azimuth, Elevation
                    dAdY = [dAds*LTSK*U,zeros(1,3)];
                    % Measurement update
                    [K, Y, P] = MeasUpdate ( Y, obsSK(idx,2), Azim, sigma_az, dAdY, P, 6 );
                    dataAll{6,iData} = [dataAll{6,iData} obsSK(idx,2)-Azim ];
                end

                if bEKF_AER(2)
                    % Elevation and partials
                    r = Y(1:3);
                    s = LTSK*(U*r-Rs{1,1});                     % Topocentric position [m]
                    [Azim, Elev, dAds, dEds] = AzElPa(s);     % Azimuth, Elevation
                    dEdY = [dEds*LTSK*U,zeros(1,3)];
                    % Measurement update
                    [K, Y, P] = MeasUpdate ( Y, obsSK(idx,3), Elev, sigma_el, dEdY, P, 6 );
                    dataAll{7,iData} = [dataAll{7,iData} obsSK(idx,3)-Elev ];
                end

                if bEKF_AER(3)
                    % Range and partials
                     r = Y(1:3);
                     s = LTSK*(U*r-Rs{1,1});                          % Topocentric position [m]
                     Dist = norm(s);
                     dDds = (s/Dist)';         % Range
                     dDdY = [dDds*LTSK*U,zeros(1,3)];
                     % Measurement update
                     [K, Y, P] = MeasUpdate ( Y, obsSK(idx,4), Dist, sigma_range, dDdY, P, 6 ); 
                     dataAll{8,iData} = [dataAll{8,iData} obsSK(idx,4)-Dist ];
                end
            end
            if bStationBigBear
                if bEKF_AER(1)
                    % Topocentric coordinates
                    r = Y(1:3);
                    s = LTCB*(U*r-Rs{1,2});     % ENZ frame     % Topocentric position [m]
                    % Azimuth and partials
                    [Azim, Elev, dAds, dEds] = AzElPa(s);     % Azimuth, Elevation
                    dAdY = [dAds*LTCB*U,zeros(1,3)];
                    % Measurement update
                    [K, Y, P] = MeasUpdate ( Y, obsCB(idx,2), Azim, sigma_az, dAdY, P, 6 );
                    dataAll{6,iData+1} = [dataAll{6,iData+1} obsCB(idx,2)-Azim ];
                end
                if bEKF_AER(2)
                    % Elevation and partials
                    r = Y(1:3);
                    s = LTCB*(U*r-Rs{1,2});                     % Topocentric position [m]
                    [Azim, Elev, dAds, dEds] = AzElPa(s);     % Azimuth, Elevation
                    dEdY = [dEds*LTCB*U,zeros(1,3)];
                    % Measurement update
                    [K, Y, P] = MeasUpdate ( Y, obsCB(idx,3), Elev, sigma_el, dEdY, P, 6 );
                    dataAll{7,iData+1} = [dataAll{7,iData+1} obsCB(idx,3)-Elev ];
                end
                if bEKF_AER(3)
                    % Range and partials
                     r = Y(1:3);
                     s = LTCB*(U*r-Rs{1,2});                          % Topocentric position [m]
                     Dist = norm(s); dDds = (s/Dist)';         % Range
                     dDdY = [dDds*LTCB*U,zeros(1,3)];
                     % Measurement update
                     [K, Y, P] = MeasUpdate ( Y, obsCB(idx,4), Dist, sigma_range, dDdY, P, 6 );
                     dataAll{8,iData+1} = [dataAll{8,iData+1} obsCB(idx,4)-Dist ];
                end        
            end

            if bStationKiruna
            % Data saved
                dataAll{2,iData} = [dataAll{2,iData} K];
                dataAll{3,iData} = [dataAll{3,iData} [P(1,1); P(2,2); P(3,3); P(4,4); P(5,5);P(6,6)]];
                dataAll{5,iData} = [dataAll{5,iData} Y];
                [~, a, e, i, Omega, omega, M] = elements(Y);
                dataAll{4,iData} = [dataAll{4,iData} [a,e,i,Omega,omega,M]'];
            end
            if bStationBigBear
                 % Data saved
                dataAll{2,iData+1} = [dataAll{2,iData+1} K];
                dataAll{3,iData+1} = [dataAll{3,iData+1} [P(1,1); P(2,2); P(3,3); P(4,4); P(5,5);P(6,6)]];
                dataAll{5,iData+1} = [dataAll{5,iData+1} Y];
                [~, a, e, i, Omega, omega, M] = elements(Y);
                dataAll{4,iData+1} = [dataAll{4,iData+1} [a,e,i,Omega,omega,M]'];     
            end    

            if mod(idx-jData,25)==0
                disp(idx-jData);
            end
        end
        GlobalRec{jData-(iStart-1),(iData+1)/2} = dataAll(:,iData);
    end
end

%%
kepler = 4;
cartesian = 5;
stateVar = 2;
element2display = kepler;

for iData = 1:size(GlobalRec,2)
    value{iData} = [];
    finalLine{iData} = [];
    for jData = 1:size(GlobalRec,1)
        value{iData} = [value{iData}; GlobalRec{jData,iData}{element2display,1}(stateVar,:)];
        finalLine{iData} = [finalLine{iData} GlobalRec{jData,iData}{element2display,1}(stateVar,end)];
    end 
    figure;
    imagesc(value{iData}(1:end,:));
    xlabel('x: time');
    ylabel('y: time translation');
    zlabel('z: sma');
    figure; mesh(value{iData});
    xlabel('x: time');
    ylabel('y: time translation');
    zlabel('z: sma');
    if iData == 1 
        title('No maneuver');
    else
        title(sprintf('Delta maneuver %i',tManeuver-iStart));
    end 
end

lMax = 1.05*max([finalLine{1} finalLine{2}]);
lMin = 0.95*min([finalLine{1} finalLine{2}]);

for iData = 1:size(GlobalRec,2)
    figure;
    plot(finalLine{iData});
    xlabel('x: time');
    ylabel('y: sma');  
    grid on;
    if iData == 1 
        title('No maneuver');
        hold on;
        plot(data{1,9}(iStart:end,stateVar));
    else
        title(sprintf('Delta maneuver %i',tManeuver-iStart));
        hold on;
        plot(dataMan{1,9}(iStart:end,stateVar));
    end  
    ylim([lMin lMax]);
end


diff = value{1} - value{2};
figure; mesh(diff);
xlabel('x: time');
ylabel('y: time translation');
zlabel('z: sma');
title(sprintf('Error: Delta maneuver %i',tManeuver-iStart));


%%
disp('EoS')
return;
