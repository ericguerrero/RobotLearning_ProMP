function [ demoY] = loadTraj(trajPath, Ndemos)
%This function load the trajectories and adapt the roll, pitch
%and yaw
for k=1:Ndemos
    traj = load(sprintf('%s/Cartesian%d.txt', trajPath, k-1));   
    demoY(:,:,k)=traj; % (time, dof, demo)
end

%% post processing for ROLL - PITCH - YAW

%maybe the orther of roll pitch and yaw is not the same of th eprogram but
%the results is the same

%check if the first value (of roll pitch and yaw) of the first trajectory
%is negative or positive
sign_rpw=[]; %sign=[sign_roll sign_pitch sign_yaw]; 
sign_rpw=[sign(demoY(1,4,1)) sign(demoY(1,5,1)) sign(demoY(1,6,1))]

% check if the roll pitch and yaw start with the same sign or the first
% traj, if not change sign to be the same of the first one

end

