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

PHI = zeros(Nf,Nt);
for i=1:Nt      % trajectory points
    for s=1:Nf  % basis functions
        PHI(s,i)=evalexp(i/Nt,C(s),D);
    end
end

% Via Points
T = round([0.9,0.9,0.9]'*Nt);
Y = [0.40,0.35,-0.22]';
via_point_var = 0.00001;

for k=1:3
    % Get weights
    demos = reshape(demoY(:,k,:),Nt,Ndemos);
    w = pinv(PHI')*demos;

    % Mean
    mu_w = mean(w,2);

    % Standard deviation
    std_w = std(w,0,2);

    % Covariance
    cov_w = (w-mu_w*ones(1,Ndemos))*(w-mu_w*ones(1,Ndemos))'/Ndemos;

    %Mean Trajectory
    mu_y = PHI'*mu_w; 
    

    % Confidence intervals
    std_y = 2*real(sqrt(diag(PHI'*cov_w*PHI)));
    u_ci = mu_y + std_y;
    l_ci = mu_y - std_y;
    
    %Via points
    t = T(k);
    y = Y(k);

    mu_w_vp = mu_w + cov_w*PHI(:,t)*inv(via_point_var + PHI(:,t)'*cov_w*PHI(:,t))*(y-PHI(:,t)' * mu_w);
    cov_w_vp = cov_w - cov_w*PHI(:,t)*inv(via_point_var + PHI(:,t)'*cov_w*PHI(:,t)) * PHI(:,t)'*cov_w;

    mu_y_vp = PHI'*mu_w_vp; 
    
    
    std_vp = 2*real(sqrt(diag(PHI'*cov_w_vp*PHI)));
    
    u_ci_vp = mu_y_vp + std_vp;
    l_ci_vp = mu_y_vp - std_vp;


    
    % Plots
%     DOF = {'x[m]','y[m]','z[m]'};
%     figure(k);hold on;grid on;
%     plot(x',mu_y,'r','lineWidth',2)
%     plot(x',u_ci,'r')
%     for i=1:Ndemos
%         plot(x',demoY(:,k,i),'b','lineWidth',1)
%     end
%     fill_between_lines(x',u_ci,l_ci,[1 0 0],0.2)
%     plot(x',l_ci,'r')
%     plot(x',mu_y,'r','lineWidth',2)
%     plot(x',u_ci,'r')
%     xlabel('time [0-1]');ylabel(DOF(k))

    DOF = {'x[m]','y[m]','z[m]'};
    figure(k+3);hold on;grid on;
    plot(x',mu_y,'r','lineWidth',2)
    plot(x',mu_y_vp,'g','lineWidth',2)
    plot(x',u_ci,'r')
    plot(x',u_ci_vp,'g')
    fill_between_lines(x',u_ci,l_ci,[1 0 0],0.2)
    fill_between_lines(x',u_ci_vp,l_ci_vp,[0 1 0],0.2)
    plot(x',l_ci,'r')
    plot(x',l_ci_vp,'g')
    plot(x',mu_y,'r','lineWidth',2)
    plot(x',mu_y_vp,'g','lineWidth',2)
    plot(x',u_ci,'r')
    plot(x',u_ci_vp,'g')
    xlabel('time [0-1]');ylabel(DOF(k))
    plot(T(k)/Nt,Y(k),'bX','markerSize',10,'lineWidth',2)
end

rmpath(trajpath);
