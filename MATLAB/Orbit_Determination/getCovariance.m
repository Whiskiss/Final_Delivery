% Each choise is a different configuration of covariance matrix creation at
% each iteration
% Default case: iChoice == 99
%%
% Q process noise
% P covariance matrix
% Phi state transition matrix
%%
% Bibliography
% Matthew B Rhudy
% Increasing the Convergence Rate of the Extended Kalman Filter
% https://www.researchgate.net/publication/282223245 
% 2015

%%
if iChoice == 1
    Q = 0 * eye(size(Phi));
    if mod(idx,1)==0 %10
        Q11 = diagQ(1) * eye(3); 
        Q12 = zeros(3);
        Q21 = zeros(3);
        Q22 = diagQ(2) * eye(3);
        Q = [Q11,Q12;Q21,Q22];
    end
    
    % Event detector large time gap, high residuals, ... to be defined by
    % the user
    
    if abs(t-t_old)>100 || (idx>20 && abs(dataEKF{6,iData}(end))>0.001)
        disp('Variance stimulus');
        idxTimeGap = 0;
    end
    
    if idxTimeGap <10
        %P = TimeUpdate(10*P, Phi, Q);
        P = TimeUpdate(3*P, Phi, Q);
        idxTimeGap = idxTimeGap+1;
    else
        P_old = P;
        P = TimeUpdate(1*P, Phi, Q);
    end
    
%%
elseif iChoice == 2
    Q = 0 * eye(size(Phi));
    if mod(idx,1)==0 %10
        Q11 = 10*abs(t-t_old) * eye(3);
        Q12 = zeros(3);
        Q21 = zeros(3);
        Q22 = 1*abs(t-t_old) * eye(3); 
        Q = [Q11,Q12;Q21,Q22];
    end
    
    % Event detector large time gap, high residuals, ... to be defined by
    % the user
    if idx>20 && abs(dataEKF{6,iData}(end)) >= 0.001
        disp('Variance stimulus');
        idxTimeGap = 0;
    end

    if idxTimeGap <10
        P = TimeUpdate(10*P, Phi, Q);
        idxTimeGap = idxTimeGap+1;
    else
        P_old = P;
        P = TimeUpdate(3*P, Phi, Q);
    end    
    
%%
elseif iChoice == 3
    Q = 0 * eye(size(Phi));
    if mod(idx,1)==0 %10
        Q11 = 10*abs(t-t_old) * eye(3);
        Q12 = zeros(3);
        Q21 = zeros(3);
        Q22 = 1*abs(t-t_old) * eye(3); 
        Q = [Q11,Q12;Q21,Q22];
    end
    
    % Event detector large time gap, high residuals, ... to be defined by
    % the user
    if idx>20 && abs(dataEKF{6,iData}(end)) >= 0.001
        disp('Variance stimulus');
        idxTimeGap = 0;
    end

    if idxTimeGap <10
        P = TimeUpdate(10*P, Phi, Q);
        idxTimeGap = idxTimeGap+1;
    else
        P_old = P;
        P = TimeUpdate(3*P, Phi, Q);
    end     
   
%%    
elseif iChoice == 4
    % For Maneuver_detection
    Q = 0 * eye(size(Phi));
    if mod(idx,1)==0 %10
        Q11 = diagQ(1) * eye(3); 
        Q12 = zeros(3);
        Q21 = zeros(3);
        Q22 = diagQ(2) * eye(3);
        Q = [Q11,Q12;Q21,Q22];
    end
    
    % Event detector large time gap, high residuals, ... to be defined by
    % the user
    if abs(t-t_old)>100 
        disp('Variance stimulus');
        idxTimeGap = 0;
    end
    
    if idxTimeGap <10
        %P = TimeUpdate(10*P, Phi, Q);
        P = TimeUpdate(10*P, Phi, Q);
        idxTimeGap = idxTimeGap+1;
    else
        P_old = P;
        P = TimeUpdate(3*P, Phi, Q);
    end
        
%% 
elseif iChoice == 5
    Q = 0 * eye(size(Phi));
    if mod(idx,1)==0 %10
        Q11 = diagQ(1) * eye(3); 
        Q12 = zeros(3);
        Q21 = zeros(3);
        Q22 = diagQ(2) * eye(3);
        Q = [Q11,Q12;Q21,Q22];
    end
    
    % Event detector large time gap, high residuals, ... to be defined by
    % the user
    
    if abs(t-t_old)>100 || (idx>1 && abs(dataEKF{6,iData}(end))>0.001)
        disp('Variance stimulus');
        idxTimeGap = 0;
    end
    
    if idxTimeGap <10
        %P = TimeUpdate(10*P, Phi, Q);
        P = TimeUpdate(10*P, Phi, Q);
        idxTimeGap = idxTimeGap+1;
    else
        P_old = P;
        P = TimeUpdate(3*P, Phi, Q);
    end

%%
elseif iChoice >= 99
    % No change
    Q = 0 * eye(size(Phi));
    if mod(idx,1)==0 %10
        Q11 = diagQ(1) * eye(3); 
        Q12 = zeros(3);
        Q21 = zeros(3);
        Q22 = diagQ(2) * eye(3);
        Q = [Q11,Q12;Q21,Q22];
    end
    
    P = TimeUpdate(1*P, Phi, Q);
end