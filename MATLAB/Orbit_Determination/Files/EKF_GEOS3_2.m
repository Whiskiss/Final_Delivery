%--------------------------------------------------------------------------
%
%  Initial Orbit Determination using Gauss and Extended Kalman Filter methods
%
% References:
%   O. Montenbruck, E. Gill, "Satellite Orbits - Models, Methods, and
%   Applications", Springer Verlag, Heidelberg, 2000.
%   
%   D. Vallado, "Fundamentals of Astrodynamics and Applications", 
%   4th Edition, 2013.
%
%   G. Seeber, "Satellite Geodesy", 2nd Edition, 2003.
%
% Last modified:   2020/03/16   Meysam Mahooti 
%--------------------------------------------------------------------------
clc
clear all
clearvars -except K_all_Man Cov_all_Man K_all Cov_all
format long g
%close all
%%
global const Cnm Snm AuxParam eopdata n_eqn PC

SAT_Const
load DE430Coeff.mat
PC = DE430Coeff;

% % load LEO_1.mat
%  load MAN_2_all.mat
% %data_Man = data;
% %data_Man = data(find(data(:,7)>0),:);
% a = 2000; %2000
% b = 2100; %2100
% data_Man = data(a:b,:);
% %data_Man = data;
% data = [];
% 
%  load NoMAN_1_all.mat
% %data_NoMan = data;
% %data_NoMan = data(find(data(:,7)>0),:);
% data_NoMan = data(a:b,:);
%%
%  load HPOP_data_3.mat
% %dataAll = {data_Man(1:400,:),data_NoMan(1:400,:)};
% for idx =1:size(dataMan{1,1},1)
% 	[~, a, e, i, Omega, omega, M] = elements(dataMan{1,1}(idx,2:end)') ;
%     dataMan{1,1}(idx,8) = a;
%     dataMan{1,1}(idx,9) = e;
%     dataMan{1,1}(idx,10) = i;
%     dataMan{1,1}(idx,11) = Omega;    
%     dataMan{1,1}(idx,12) = omega;
%     dataMan{1,1}(idx,13) = M;    
% end
% for idx =1:size(data{1,1},1)
% 	[~, a, e, i, Omega, omega, M] = elements(data{1,1}(idx,2:end)') ;
%     data{1,1}(idx,8) = a;
%     data{1,1}(idx,9) = e;
%     data{1,1}(idx,10) = i;
%     data{1,1}(idx,11) = Omega;    
%     data{1,1}(idx,12) = omega;
%     data{1,1}(idx,13) = M;    
% end
% iStart = 1; %50
% iEnd = min(200,size(data{1,column},1)); %iStart+50
% dataAll = {data{1,column}(iStart:iEnd,:),dataMan{1,4}(iStart:iEnd,:)};
% timeIdx = [1 9 18]; %[1,3,5]; %[1,9,18];

%%
% load HPOP_ENVISAT_Nopert.mat
load HPOP_ideal.mat
%dataAll = {data_Man(1:400,:),data_NoMan(1:400,:)};
for idx =1:size(data{1,1},1)
	[~, a, e, i, Omega, omega, M] = elements(data{1,1}(idx,2:end)') ;
    data{1,1}(idx,8) = a;
    data{1,1}(idx,9) = e;
    data{1,1}(idx,10) = i;
    data{1,1}(idx,11) = Omega;    
    data{1,1}(idx,12) = omega;
    data{1,1}(idx,13) = M;    
end
column = 4;
iStart = 29; %50
iEnd = min(127,size(data{1,column},1)); %iStart+50
dataAll = {[data{1,column}(iStart:iEnd,:) ]};
%dataAll{1,1}(find(data{1,column}(:,7)<=0),:) = [];
%dataAll{1,1}=dataAll{1,1}(1:26,:);


timeIdx = [1,5,7]; %[1,3,5]; %[1,9,18];



% % Noise added
% for idx = 1:size(dataAll,2)
%     dataAll{1,1}(:,7) = awgn(dataAll{1,1}(:,7),10);
%     dataAll{1,1}(:,8) = awgn(dataAll{1,1}(:,8),10);
% end

% for i = 1:size(dataAll{1,1},1)-1
%     if abs(dataAll{1,1}(i,8)-dataAll{1,1}(i+1,8))>300
%         dataAll{1,1}(i+1:end,8) = dataAll{1,1}(i+1:end,8)+360;
%     end
% end





%%

% % Noise added
% for idx = 1:size(dataAll,2)
%     dataAll{1,1}(:,7) = awgn(dataAll{1,1}(:,7),10);
%     dataAll{1,1}(:,8) = awgn(dataAll{1,1}(:,8),10);
% end

% load LEO_1.mat
% dataAll = {data};



% load LEO_1.mat
% dataAll = {data(20:30,:)};

K_all = [];
Cov_all = [];
kep = [];
Y_all = [];


for iData = 1:size(dataAll,2)
    dataAll{2,iData} = []; %K
    dataAll{3,iData} = []; %Var
    dataAll{4,iData} = []; %Kep
    dataAll{5,iData} = []; %Y
    dataAll{6,iData} = []; %Residuals Az
    dataAll{7,iData} = []; %Residuals El
    
    Cnm = zeros(181,181);
    Snm = zeros(181,181);
    % Global Gravity Model
    fid = fopen('GGM03S.txt','r');
    for n=0:180
        for m=0:n
            temp = fscanf(fid,'%d %d %f %f %f %f',[6 1]);        
            Cnm(n+1,m+1) = temp(3);
            Snm(n+1,m+1) = temp(4);
        end
    end

    % Model parameters
    AuxParam = struct ('Mjd_UTC',0,'n',0,'m',0);

    % read Earth orientation parameters
    fid = fopen('eop19620101.txt','r');
    %  ----------------------------------------------------------------------------------------------------
    % |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
    % |(0h UTC)           "         "          s          s          "        "          "         "     s 
    %  ----------------------------------------------------------------------------------------------------
    eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
    fclose(fid);

    nobs = size(dataAll{1,iData},1);
    obs = zeros(nobs,4);

    % read observations
    for idx=1:size(dataAll{1,iData},1)
        Y = dataAll{1,iData}(idx,1);
        M = dataAll{1,iData}(idx,2);
        D = dataAll{1,iData}(idx,3);
        h = dataAll{1,iData}(idx,4);
        m = dataAll{1,iData}(idx,5);
        s = dataAll{1,iData}(idx,6);
        az1 = dataAll{1,iData}(idx,8);
        el1 = dataAll{1,iData}(idx,7);

        %Dist = str2num(tline(45:54));
        obs(idx,1) = Mjday(Y,M,D,h,m,s);
        obs(idx,2) = const.Rad*az1;
        obs(idx,3) = const.Rad*el1;
        
        %obs(i,4) = 1e3*Dist;
    end

    sigma_range = 92.5;          % [m]
    sigma_az = 0.0224*const.Rad; % [rad]
    sigma_el = 0.0139*const.Rad; % [rad]
    %sigma_az = 0.0001; % [rad]
    %sigma_el = 0.0001; % [rad]

    % Kiruna Point station
    lat = const.Rad*(67.8790708);     % [rad]
    lon = const.Rad*(21.038); % [rad]
    alt = 527.0;                % [m]
    Rs = Position(lon, lat, alt)';

    lat = const.Rad*34.258; %52.3943; 
    lon = const.Rad*116.921; %16.8746;
    alt = 2067.0; %85;
    Rs = [Rs Position(lon, lat, alt)'];
    
    % IOD using Gibbs

    Mjd1 = obs(timeIdx(1),1);
    Mjd2 = obs(timeIdx(2),1);
    Mjd3 = obs(timeIdx(3),1);

    [r2,v2] = anglesg(obs(timeIdx(1),2),obs(timeIdx(2),2),obs(timeIdx(3),2),obs(timeIdx(1),3),obs(timeIdx(2),3),obs(timeIdx(3),3),...
                      Mjd1,Mjd2,Mjd3,Rs(:,1),Rs(:,1),Rs(:,1));
    %[r2,v2] = anglesdr(obs(timeIdx(1),2),obs(timeIdx(2),2),obs(timeIdx(3),2),obs(timeIdx(1),3),obs(timeIdx(2),3),obs(timeIdx(3),3),...
    %                   Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);

    Y0_apr = [r2;v2];
    Mjd0 = Mjday(dataAll{1,iData}(1,1),dataAll{1,iData}(1,2),dataAll{1,iData}(1,3),dataAll{1,iData}(1,4),dataAll{1,iData}(1,5),dataAll{1,iData}(1,6));
    Mjd_UTC = obs(timeIdx(2),1);

    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.n      = 1; %20
    AuxParam.m      = 1; %20
    AuxParam.sun     = 0; %1
    AuxParam.moon    = 0; %1
    AuxParam.planets = 0; %1

    n_eqn  = 6;

    Y = DEInteg(@Accel,0,-(obs(timeIdx(2),1)-Mjd0)*86400.0,1e-13,1e-6,6,Y0_apr);

    % Covariance matrix initialization
    P = zeros(6);
    for idx=1:3
        P(idx,idx)=1e8; %1e8;
    end
    for idx=4:6
        P(idx,idx)=1e3; %1e3;
    end

    LT = LTC(lon,lat);

    yPhi = zeros(42,1);
    Phi  = zeros(6);

    % Measurement loop
    t = 0;


    disp(nobs);
%     figure;
   
    for idx=1:1:min(10000,nobs)
        % Previous step
        t_old = t;
        Y_old = Y;

        % Time increment and propagation
        Mjd_UTC = obs(idx,1)  ;               % Modified Julian Date
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



        % Time update
        
        Q = 0 * eye(size(Phi));
        val1 = 0; %0.01; %0.01*randi(10);%0.01;
        val2 = val1; %0.001*randi(10);%0.001;
        if mod(idx,1)==0 %10
            Q11 = val1 * eye(3); %1 %0.0000001
            Q12 = zeros(3);
            Q21 = zeros(3);
            Q22 = val2 * eye(3); %1 %0.0000001
            Q = [Q11,Q12;Q21,Q22];
        end
        
        % Predicted error coviarance
        P = TimeUpdate(P, Phi, Q);

       bTest = 0;
       if bTest
            K1 = 0;
            Y1 = 0;
            P1 = 0;

            for iStation=1:2
                % Topocentric coordinates
                theta = gmst(Mjd_UT1);                    % Earth rotation
                U = R_z(theta);
                r = Y(1:3);
                s = LT*(U*r-Rs(:,iStation));                          % Topocentric position [m]

                % Azimuth and partials
                [Azim, Elev, dAds, dEds] = AzElPa(s);     % Azimuth, Elevation
                dAdY = [dAds*LT*U,zeros(1,3)];
                % Measurement update
                [Kaz, Yaz, Paz] = MeasUpdate ( Y, obs(idx,2*iStation), Azim, sigma_az, dAdY, P, 6 );
                %dataAll{6,iData} = [dataAll{6,iData} obs(idx,2*iStation)-Azim ];
                K1 = K1+Kaz;
                Y1 = Y1+Yaz;
                P1 = P1+Paz;

                % Elevation and partials
                r = Y(1:3);
                s = LT*(U*r-Rs(:,iStation));                          % Topocentric position [m]
                [Azim, Elev, dAds, dEds] = AzElPa(s);     % Azimuth, Elevation
                dEdY = [dEds*LT*U,zeros(1,3)];
                % Measurement update
                [Kel, Yel, Pel] = MeasUpdate ( Y, obs(idx,1+2*iStation), Elev, sigma_el, dEdY, P, 6 );
                %dataAll{7,iData} = [dataAll{7,iData} obs(idx,1+2*iStation)-Elev ];
                K1 = K1+Kel;
                Y1 = Y1+Yel;
                P1 = P1+Pel;
            end
            K = K1./4;
            Y = Y1./4;
            P = P1./4;
        else
            % Topocentric coordinates
            theta = gmst(Mjd_UT1);                    % Earth rotation
            U = R_z(theta);
            r = Y(1:3);
            s = LT*(U*r-Rs(:,1));     % ENZ frame     % Topocentric position [m]


            % Azimuth and partials
            [Azim, Elev, dAds, dEds] = AzElPa(s);     % Azimuth, Elevation
            
            
            dAdY = [dAds*LT*U,zeros(1,3)];
            % Measurement update
            [K, Y, P] = MeasUpdate ( Y, obs(idx,2), Azim, sigma_az, dAdY, P, 6 );
            dataAll{6,iData} = [dataAll{6,iData} obs(idx,2*1)-Azim ];
            
            % Elevation and partials
            % r = Y(1:3);
            % s = LT*(U*r-Rs(:,1));                          % Topocentric position [m]
            % [Azim, Elev, dAds, dEds] = AzElPa(s);     % Azimuth, Elevation
            dEdY = [dEds*LT*U,zeros(1,3)];
            % Measurement update
            [K, Y, P] = MeasUpdate ( Y, obs(idx,3), Elev, sigma_el, dEdY, P, 6 );
            dataAll{7,iData} = [dataAll{7,iData} obs(idx,1+2*1)-Elev ];
        end
        
        % Data saved
        dataAll{2,iData} = [dataAll{2,iData} K];
        dataAll{3,iData} = [dataAll{3,iData} [P(1,1); P(2,2); P(3,3); P(4,4); P(5,5);P(6,6)]];
        dataAll{5,iData} = [dataAll{5,iData} Y];
        if mod(idx,25)==0
            disp(idx);
        end
        [~, a, e, idx, Omega, omega, M] = elements(Y);
        dataAll{4,iData} = [dataAll{4,iData} [a,e,idx,Omega,omega,M]'];
        
    end
end
idx = size(obs,1);
[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,obs(idx,1),'l');
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
Mjd_TT = Mjd_UTC + TT_UTC/86400;
AuxParam.Mjd_UTC = Mjd_UTC;
AuxParam.Mjd_TT = Mjd_TT;
Y0 = DEInteg (@Accel,0,-(obs(idx,1)-obs(2,1))*86400.0,1e-13,1e-6,6,Y);

% err = data{1,1}(1,2:7)-Y0';
% correct = dataAll{5,iData} - repmat(err,[size(dataAll{5,iData},2) 1])';
% val = [];
% for i=1:size(correct,2)
%     [~, a, e, idx, Omega, omega, M] = elements(correct(:,i));
%     val = [val [a,e,idx,Omega,omega,M]'];
% end
% figure;plot(val(1,:))

%%
figure;
elt = 1;
plot(dataAll{4,1}(elt,1:end));
hold on;
plot(data{1,1}(iStart:iEnd,elt+7));

figure;
text = {'x','y','z','vx','vy','vz'};
for i=1:6
    subplot(6,1,i);
    plot(data{1,1}(iStart:iEnd,i+1));
    hold on;
    plot(dataAll{5,1}(i,:)');
    grid on;
    legend(text{i});
end

errX = data{1,1}(iStart:iEnd,2) - dataAll{5,1}(1,:)';
errY = data{1,1}(iStart:iEnd,3) - dataAll{5,1}(2,:)';
errZ = data{1,1}(iStart:iEnd,4) - dataAll{5,1}(3,:)';

err = sqrt(errX.^2 + errY.^2 + errZ.^2);
figure;
plot(err);
grid on;
disp('EoS')

return;
%%
dataAll{2,3} = dataAll{2,1} - dataAll{2,2};
dataAll{3,3} = dataAll{3,1} - dataAll{3,2};
dataAll{4,3} = dataAll{4,1} - dataAll{4,2};
dataAll{5,3} = dataAll{5,1} - dataAll{5,2};

idx = size(obs,1);
[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,obs(idx,1),'l');
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
Mjd_TT = Mjd_UTC + TT_UTC/86400;
AuxParam.Mjd_UTC = Mjd_UTC;
AuxParam.Mjd_TT = Mjd_TT;

Y0 = DEInteg (@Accel,0,-(obs(idx,1)-obs(1,1))*86400.0,1e-13,1e-6,6,Y);

%%
plotData( dataAll, 1, 20000 );
figure;
plot(dataAll{4,1}(1,1:end));
hold on;
plot(data{1,1}(:,8));
plot(data{1,column}(:,8)*10^6);

%%
[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,obs(46,1),'l');
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
Mjd_TT = Mjd_UTC + TT_UTC/86400;
AuxParam.Mjd_UTC = Mjd_UTC;
AuxParam.Mjd_TT = Mjd_TT;

Y0 = DEInteg (@Accel,0,-(obs(46,1)-obs(1,1))*86400.0,1e-13,1e-6,6,Y);

Y_true = [5753.173e3, 2673.361e3, 3440.304e3, 4.324207e3, -1.924299e3, -5.728216e3]';

fprintf('\nError of Position Estimation\n');
fprintf('dX%10.1f [m]\n',Y0(1)-Y_true(1));
fprintf('dY%10.1f [m]\n',Y0(2)-Y_true(2));
fprintf('dZ%10.1f [m]\n',Y0(3)-Y_true(3));
fprintf('\nError of Velocity Estimation\n');
fprintf('dVx%8.1f [m/s]\n',Y0(4)-Y_true(4));
fprintf('dVy%8.1f [m/s]\n',Y0(5)-Y_true(5));
fprintf('dVz%8.1f [m/s]\n',Y0(6)-Y_true(6));

