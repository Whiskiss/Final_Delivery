clc
clear all
clearvars -except K_all_Man Cov_all_Man K_all Cov_all
format long g

%%
global const Cnm Snm AuxParam eopdata n_eqn PC

SAT_Const
load DE430Coeff.mat
PC = DE430Coeff;

%%
load HPOP_ENVISAT_Nopert.mat
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
iStart = 35; %50
iEnd = min(127,size(data{1,4},1)); %iStart+50

ival = 1;
bReset = false;
dataAll{1,ival} = [];
for idx =1:size(data{1,4},1)-1
    if (data{1,4}(idx,8)-data{1,4}(idx+1,8))>300
        ival = ival+1;
        dataAll{1,ival} = [];
        bReset = true;
    elseif data{1,4}(idx,8)>24
        bReset = false;
        dataAll{1,ival} = [dataAll{1,ival}; data{1,4}(idx,:)];
    end
end


timeIdx = [1,9,18]; %[1,3,5]; %[1,9,18];

K_all = [];
Cov_all = [];
kep = [];
Y_all = [];


% Covariance matrix initialization
P = zeros(6);
for idx=1:3
    P(idx,idx)=1e4; %1e8;
end
for idx=4:6
    P(idx,idx)=1e2; %1e3;
end
yPhi = zeros(42,1);
Phi  = zeros(6);


%%
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


    %%
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
        az = dataAll{1,iData}(idx,8);
        el = dataAll{1,iData}(idx,7);
        %Dist = str2num(tline(45:54));
        obs(idx,1) = Mjday(Y,M,D,h,m,s);
        obs(idx,2) = const.Rad*az;
        obs(idx,3) = const.Rad*el;
        %obs(i,4) = 1e3*Dist;
    end

    sigma_range = 92.5;          % [m]
    sigma_az = 0.0224*const.Rad; % [rad]
    sigma_el = 0.0139*const.Rad; % [rad]

    % Kiruna Point station
    lat = const.Rad*67.8790708;     % [rad]
    lon = const.Rad*(21.038); % [rad]
    alt = 527.0;                % [m]
    Rs = Position(lon, lat, alt)';

    % IOD using Gibbs
    Mjd1 = obs(timeIdx(1),1);
    Mjd2 = obs(timeIdx(2),1);
    Mjd3 = obs(timeIdx(3),1);

    [r2,v2] = anglesg(obs(timeIdx(1),2),obs(timeIdx(2),2),obs(timeIdx(3),2),obs(timeIdx(1),3),obs(timeIdx(2),3),obs(timeIdx(3),3),...
                      Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);

    Y0_apr = [r2;v2];

    Mjd0 = Mjday(dataAll{1,iData}(1,1),dataAll{1,iData}(1,2),dataAll{1,iData}(1,3),dataAll{1,iData}(1,4),dataAll{1,iData}(1,5),dataAll{1,iData}(1,6));

    Mjd_UTC = obs(timeIdx(2),1);

    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.n      = 70; %20
    AuxParam.m      = 70; %20
    AuxParam.sun     = 1; %1
    AuxParam.moon    = 1; %1
    AuxParam.planets = 0; %1

    n_eqn  = 6;

    Y = DEInteg(@Accel,0,-(obs(timeIdx(2),1)-Mjd0)*86400.0,1e-13,1e-6,6,Y0_apr);

    LT = LTC(lon,lat);



    % Measurement loop
    t = 0;

    disp(nobs);
    for idx=1:1:min(10000,nobs)
        % Previous step
        t_old = t;
        Y_old = Y;

        % Time increment and propagation
        Mjd_UTC = obs(idx,1);                 % Modified Julian Date
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

        % Topocentric coordinates
        theta = gmst(Mjd_UT1);                    % Earth rotation
        U = R_z(theta);
        r = Y(1:3);
        s = LT*(U*r-Rs);                          % Topocentric position [m]

        % Time update

        Q = 0 * eye(size(Phi));
        val1 = 0; %0.01*randi(10);%0.01;
        val2 = 0; %0.001*randi(10);%0.001;
        if mod(idx,5)==0 %10
            Q11 = val1 * eye(3); %1 %0.0000001
            Q12 = zeros(3);
            Q21 = zeros(3);
            Q22 = val2 * eye(3); %1 %0.0000001
            Q = [Q11,Q12;Q21,Q22];
        end

        % Predicted error coviarance
        P = TimeUpdate(P, Phi, Q);

        % Azimuth and partials
        [Azim, Elev, dAds, dEds] = AzElPa(s);     % Azimuth, Elevation
        dAdY = [dAds*LT*U,zeros(1,3)];
        % Measurement update
        [K, Y, P] = MeasUpdate ( Y, obs(idx,2), Azim, sigma_az, dAdY, P, 6 );
        dataAll{6,iData} = [dataAll{6,iData} obs(idx,2)-Azim ];


        % Elevation and partials
        r = Y(1:3);
        s = LT*(U*r-Rs);                          % Topocentric position [m]
        [Azim, Elev, dAds, dEds] = AzElPa(s);     % Azimuth, Elevation
        dEdY = [dEds*LT*U,zeros(1,3)];
        % Measurement update
        [K, Y, P] = MeasUpdate ( Y, obs(idx,3), Elev, sigma_el, dEdY, P, 6 );
        dataAll{7,iData} = [dataAll{7,iData} obs(idx,3)-Elev ];


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


%%
for iData=1:size(dataAll,2)
    figure;
    plot(dataAll{4,iData}(1,1:end));
    hold on;
    plot(data{1,1}(:,8));
end