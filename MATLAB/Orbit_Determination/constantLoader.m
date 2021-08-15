%% Ephemerids loading
global const Cnm Snm AuxParam eopdata n_eqn PC
SAT_Const
load DE430Coeff.mat
PC = DE430Coeff;

% Model parameters
AuxParam = struct ('Mjd_UTC',0,'n',0,'m',0);

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


% read Earth orientation parameters
fid = fopen('eop19620101.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------
eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
fclose(fid);