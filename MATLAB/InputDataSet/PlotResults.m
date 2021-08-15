clear all
clc
close all
%%

path = 'C:\Users\antho\OneDrive\PUT\SAAS\Thesis\From_Mathwork\InputDataSet\HPOP_';

idx = 1;
for i = 0:6
    strcat(path,sprintf('%i',i),'.mat')
    load(strcat(path,sprintf('%i',i),'.mat'));
    Gdata{idx} = data;
    GdataMan{idx} = dataMan;
    idx = idx+1;
end

%%
txtKep = {'Semimajor axis [m]', 'Eccentricity', 'Inclination [rad]', ...
    'Longitude of the ascending node [rad]', 'Argument of pericenter [rad]', ...
    'Mean anomaly [rad]'};
txtPV = {'x [m]', 'y [m]', 'z [m]', ...
    'vx [m/s]','vy [m/s]', 'vz [m/s]'};

iMin = 0+1;
iMax = 3+1;

for elt = 2:7
    figure;
    variable = 1;
    %elt = 1;
    for i = iMin:iMax
        plot(Gdata{i}{variable}(:,elt),'DisplayName',sprintf('%i',i-1));
        hold on;
        grid on;

        ylabel(txtPV{elt-1});
        xlabel('Time [Observation]');
    end
    legend show
end

for elt = 1:6
    figure;
    variable = 9;
    %elt = 1;
    for i = iMin:iMax
        plot(Gdata{i}{variable}(:,elt),'DisplayName',sprintf('%i',i-1));
        hold on;
        grid on;

        ylabel(txtKep{elt});
        xlabel('Time [Observation]');
    end
    legend show
end