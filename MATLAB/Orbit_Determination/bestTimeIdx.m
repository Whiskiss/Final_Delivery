clc
clear all
format long g
%%
global const Cnm Snm AuxParam eopdata n_eqn PC

SAT_Const
load DE430Coeff.mat
PC = DE430Coeff;

%%
load HPOP_ideal.mat
for idx =1:size(data{1,1},1)
	[~, a, e, i, Omega, omega, M] = elements(data{1,1}(idx,2:end)') ;
    data{1,1}(idx,8) = a;
    data{1,1}(idx,9) = e;
    data{1,1}(idx,10) = i;
    data{1,1}(idx,11) = Omega;    
    data{1,1}(idx,12) = omega;
    data{1,1}(idx,13) = M;    
end
column = 5;
iStart = 29; %50
iEnd = min(127,size(data{1,column},1)); %iStart+50
dataAll = {[data{1,column}(iStart:iEnd,:) ]};


err = {};

for inc=1:60
    inc
    timeIdx = [1,1+inc,1+2*inc]; %[1,3,5]; %[1,9,18];
    err{inc,2} = timeIdx;
    err{inc,3} = data{1,1}(29+inc,2:7);
    
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

    nobs = size(dataAll{1,1},1);
    obs = zeros(nobs,4);

    % read observations
    for idx=1:size(dataAll{1,1},1)
        Y = dataAll{1,1}(idx,1);
        M = dataAll{1,1}(idx,2);
        D = dataAll{1,1}(idx,3);
        h = dataAll{1,1}(idx,4);
        m = dataAll{1,1}(idx,5);
        s = dataAll{1,1}(idx,6);
        az1 = dataAll{1,1}(idx,8);
        el1 = dataAll{1,1}(idx,7);

        obs(idx,1) = Mjday(Y,M,D,h,m,s);
        obs(idx,2) = const.Rad*az1;
        obs(idx,3) = const.Rad*el1;
    end


    % Kiruna Point station
    lat = const.Rad*(67.8790708);     % [rad]
    lon = const.Rad*(21.038); % [rad]
    alt = 527.0;                % [m]
    Rs = Position(lon, lat, alt)';

    % IOD using Gibbs
    Mjd1 = obs(timeIdx(1),1);
    Mjd2 = obs(timeIdx(2),1);
    Mjd3 = obs(timeIdx(3),1);
    [r2,v2] = anglesg(obs(timeIdx(1),2),obs(timeIdx(2),2),obs(timeIdx(3),2), ...
                      obs(timeIdx(1),3),obs(timeIdx(2),3),obs(timeIdx(3),3),...
                      Mjd1,Mjd2,Mjd3, ...
                      Rs(:,1),Rs(:,1),Rs(:,1));

    Y0_apr = [r2;v2];
    err{inc,4} = Y0_apr;
    err{inc,1} = sqrt(sum((err{inc,3}'-err{inc,4}).^2));
end
