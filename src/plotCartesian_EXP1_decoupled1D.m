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

phi = zeros(Nf,Nt);
for i=1:Nt      % trajectory points
    for s=1:Nf  % basis functions
        phi(s,i)=evalexp(i/Nt,C(s),D);
    end
end

% Via Points
t = round([0.8, 0.8, 0.8]'*Nt);
y = [0.7, 0.2, -0.25]';
via_point_var = 0.3*eye(3);

for k=1:dof
    % Get weights
    demos = reshape(demoY(:,k,:),Nt,Ndemos);
    w = pinv(phi')*demos;

    % Mean
    mu_w = mean(w,2);

    % Standard deviation
    std_w = std(w,0,2);

    % Covariance
    cov_w = (w-mu_w*ones(1,Ndemos))*(w-mu_w*ones(1,Ndemos))'/Ndemos;

    %Mean Trajectory
    mu_y = phi'*mu_w; 
    

    % Confidence intervals
    upper_ci = phi'*(mu_w+1.96*std_w/sqrt(Nf)); % Confidence Intervals (95%)
    lower_ci = phi'*(mu_w-1.96*std_w/sqrt(Nf)); % Confidence Intervals (95%)


    % Sample trajectories
%     Nsamples = 5;
%     for i = 1:Nsamples
%         w_traj = mvnrnd(mu_w,cov_w)';
%         trajSampled(:,i) = phi'*w_traj;
%     end
    

    %% Plot trajectories
    figure(1); hold on;
    subplot(2,3,k);hold on;grid on;
    for i=1:Ndemos % Demos
        plot(x',demos(:,i),'b','lineWidth',1);
    end
    plot(x',mu_y,'r','lineWidth',2); %Mean
%     for i=1:Nsamples % Sampled
%         plot(x',trajSampled(:,i),'r','lineWidth',1);
%     end
    
    %% Via Point
    if k<=size(t,1)
        mu_w_vp = mu_w + cov_w*phi(:,t(k))*inv(via_point_var + phi(:,t(k))'*cov_w*phi(:,t(k)))*(y(k)-phi(:,t(k))' * mu_w);
        cov_w_vp = cov_w - cov_w*phi(:,t(k))*inv(via_point_var + phi(:,t(k))'*cov_w*phi(:,t(k))) * phi(:,t(k))'*cov_w;
        mu_y_vp = phi'*mu_w_vp;  

        % Confidence intervals
        upper_ci_vp = phi'*(mu_w_vp+1.96*diag(cov_w_vp)/sqrt(Nf)); % Confidence Intervals (95%)
        lower_ci_vp = phi'*(mu_w_vp-1.96*diag(cov_w_vp)/sqrt(Nf)); % Confidence Intervals (95%)
    
        % Sample trajectories
%         Nsamples = 5;
%         for i = 1:Nsamples
%             w_traj = mvnrnd(mu_w_vp,cov_w_vp)';
%             trajSampled_vp(:,i) = phi'*w_traj;
%         end
    
        %% Plot trajectories
        figure(1); hold on;
        subplot(2,3,k);hold on;grid on;
        plot(x',mu_y_vp,'k','lineWidth',2); %Mean VP
%         for i=1:Nsamples % Sampled VP
%             plot(x',trajSampled_vp(:,i),'g','lineWidth',1);
%         end
    end
end

%% Via Points

    % 
%     % plot2D(x,upper_ci,'--r',2)
%     % plot2D(x,lower_ci,'--r',2)
%     plot2D(x,Sampled_Traj,'r',1)
% 
% 
%     % Plot3D
%     figure(2);hold on;
%     plot3D(demoY,'b')
%     plot3D(Sampled_Traj,'r')
%     plot3(Mean_Traj(:,1),Mean_Traj(:,2),Mean_Traj(:,3),'r','lineWidth',2);
% 
%     %% Via Points Plots 
%     figure(1); hold on;
%     plot2D(x,Mean_Traj_vp,'k',2)
%     % plot2D(x,upper_ci_vp,'--g',2)
%     % plot2D(x,lower_ci_vp,'--g',2)
%     plot2D(x,Sampled_Traj_vp,'g',1)
% 
% 
%     % Plot3D
%     figure(2);hold on;
%     plot3D(Sampled_Traj_vp,'g')
%     plot3(Mean_Traj_vp(:,1),Mean_Traj_vp(:,2),Mean_Traj_vp(:,3),'g','lineWidth',2);


%% Comparative
% figure(3)
% % Weights
% mu_w_reshaped = reshape(mu_w,Nf*dof,1);
% plot(mw,'b');hold on;
% plot(mu_w_reshaped,'g')

rmpath(trajpath);
