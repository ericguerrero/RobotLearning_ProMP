function [demoY,Nt, dof, x, mw, Sw] = loadTrajectory(trajpath,trajname, Ndemos, mode)
%This function load the trajectories and adapt the roll, pitch
% and yaw to be coherent (sometimes there is a strange change of sign)
for k=1:Ndemos
    traj = load(sprintf(['%s/',mode ,'%d.txt'], trajpath, k-1));   
    demoY(:,:,k)=traj; % (time, dof, demo)
end

%% post processing for ROLL - PITCH - YAW
% THEY SHOULD BE EULER'S ANGLES - It works anyway

%check if the first value (of roll pitch and yaw) of the first trajectory
%is negative or positive
sign_rpy=[]; %sign=[sign_roll sign_pitch sign_yaw]; 
sign_rpy=[sign(demoY(1,4,1)) sign(demoY(1,5,1)) sign(demoY(1,6,1))];

% check if the roll pitch and yaw start with the same sign or the first
% traj, if not change sign to be the same of the first one
for k=2:Ndemos
    if(sign(demoY(1,4,k))~=sign_rpy(1)) 
        demoY(:,4,k)=-demoY(:,4,k);
    end
    if(sign(demoY(1,5,k))~=sign_rpy(2)) 
        demoY(:,5,k)=-demoY(:,5,k);
    end
    if(sign(demoY(1,6,k))~=sign_rpy(3)) 
        demoY(:,6,k)=-demoY(:,6,k);
    end
end

%Get data
demoY=demoY(1:size(demoY,1),:,:);
Nt=size(demoY,1);   % number of time steps
dof=size(demoY,2);  % dof
x = (0:(Nt-1))/(Nt-1);% Time [0 - 1]

% Load the weights and the covariance
mw=load(sprintf('%smw.txt', trajname)); % mean
Sw=load(sprintf('%sSw.txt', trajname)); % covariance

end

