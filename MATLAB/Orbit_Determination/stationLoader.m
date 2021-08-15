global const Cnm Snm AuxParam eopdata n_eqn

% Swedish Kiruna Point station
lat = const.Rad*(67.8790708);     % [rad]
lon = const.Rad*(21.038);         % [rad]
alt = 527.0;                      % [m]
Rs{1,1} = Position(lon, lat, alt)';
Rs{2,1} = [lon, lat, alt];

% California Big Bear station
lat = const.Rad*(34.258);  
lon = const.Rad*(-116.921); 
alt = 2067.0;
Rs{1,2} = Position(lon, lat, alt)';
Rs{2,2} = [lon, lat, alt];