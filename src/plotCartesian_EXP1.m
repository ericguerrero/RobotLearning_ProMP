clc;clear;close all;format compact;

%% Load files (specify Ndemos)
trajpath = '../experiments/narrow_passage_box';
addpath(trajpath);
trajname = 'test';  
Ndemos = 20;   
[demoY,Nt, dof, x, mw, Sw] = loadTrajectory(trajpath, trajname, Ndemos, 'Cartesian');


%% Distributions
Nf=size(mw,1)/dof;  % # of basis functions
C=(0:1:Nf-1)/Nf;    % Centers
D=0.025;            % Amplitud

PHI = obtainPhi(Nf,Nt,dof,C,D);

% Get weights
demos = reshape(demoY,Nt*dof,Ndemos);
w = pinv(PHI')*demos;

% Mean
mu_w = mean(w,2);

% Standard deviation
std_w = std(w,0,2);

% Covariance
cov_w = (w-mu_w*ones(1,Ndemos))*(w-mu_w*ones(1,Ndemos))'/Ndemos;

%Mean Trajectory
mu_y = PHI'*mu_w; 
Mean_Traj = reshape(mu_y,Nt,dof);

% Confidence intervals
upper_ci = PHI'*(mu_w+1.96*std_w/sqrt(Nf)); % Confidence Intervals (95%)
lower_ci = PHI'*(mu_w-1.96*std_w/sqrt(Nf)); % Confidence Intervals (95%)
upper_ci = reshape(upper_ci,Nt,dof);
lower_ci = reshape(lower_ci,Nt,dof);

% Sample trajectories
Sampled_Traj = sampleTrajectories(5,mu_w,cov_w,PHI,Nt,dof);

%% Via Points
t = round([0.8]'*Nt);
y = [0.4]';
via_point_var = 0.003*eye(1);

mu_w_vp = mu_w + cov_w*PHI(:,t)*inv(via_point_var + PHI(:,t)'*cov_w*PHI(:,t))*(y-PHI(:,t)' * mu_w);
cov_w_vp = cov_w - cov_w*PHI(:,t)*inv(via_point_var + PHI(:,t)'*cov_w*PHI(:,t)) * PHI(:,t)'*cov_w;

mu_y_vp = PHI'*mu_w_vp;  
Mean_Traj_vp = reshape(mu_y_vp,Nt,dof);

% Confidence intervals
upper_ci_vp = PHI'*(mu_w_vp+1.96*diag(cov_w_vp)/sqrt(Nf)); % Confidence Intervals (95%)
lower_ci_vp = PHI'*(mu_w_vp-1.96*diag(cov_w_vp)/sqrt(Nf)); % Confidence Intervals (95%)
upper_ci_vp = reshape(upper_ci_vp,Nt,dof);
lower_ci_vp = reshape(lower_ci_vp,Nt,dof);

% Sample trajectories
Sampled_Traj_vp = sampleTrajectories(5,mu_w_vp,cov_w_vp,PHI,Nt,dof);

%% Plot trajectories
figure(1); hold on;
plot2D(x,demoY,'b',1)
plot2D(x,Mean_Traj,'r',2)
% plot2D(x,upper_ci,'--r',2)
% plot2D(x,lower_ci,'--r',2)
% plot2D(x,Sampled_Traj,'r',1)


% Plot3D
% figure(2);hold on;
% plot3D(demoY,'b')
% plot3D(Sampled_Traj,'r')
% plot3(Mean_Traj(:,1),Mean_Traj(:,2),Mean_Traj(:,3),'r','lineWidth',2);

%% Via Points Plots 
figure(1); hold on;
plot2D(x,Mean_Traj_vp,'g',2)
% plot2D(x,upper_ci_vp,'--g',2)
% plot2D(x,lower_ci_vp,'--g',2)
% plot2D(x,Sampled_Traj_vp,'g',1)


% Plot3D
% figure(2);hold on;
% plot3D(Sampled_Traj_vp,'g')
% plot3(Mean_Traj_vp(:,1),Mean_Traj_vp(:,2),Mean_Traj_vp(:,3),'g','lineWidth',2);


%% Comparative
% figure(3)
% % Weights
% mu_w_reshaped = reshape(mu_w,Nf*dof,1);
% plot(mw,'b');hold on;
% plot(mu_w_reshaped,'g')

rmpath(trajpath);
