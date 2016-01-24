clc;clear;close all;format compact;

%% Load files (specify Ndemos)
% trajpath = './traj1';   % path to the folder where the files are located
% trajname = 'RL_02';     % name of the file to load (with out mw or Sw)
% trajpath = './traj2';   % path to the folder where the files are located
% trajname = 'RL_03';     % name of the file to load (with out mw or Sw)
trajpath = './narrow_passage_box';   % path to the folder where the files are located
trajname = 'test';  

Ndemos = 20;            % number of the demos = number of cartesian trajectories

for k=1:Ndemos
    traj = load(sprintf('%s/Cartesian%d.txt', trajpath, k-1));   
    demoY(:,:,k)=traj; % (time, dof, demo)
end
demoY=demoY(1:size(demoY,1),:,:);
Nt=size(demoY,1);   % number of time steps
dof=size(demoY,2);  % dof
time = 1:Nt;

% Load the weights and the covariance
mw=load(sprintf('%smw.txt', trajname)); % mean
Sw=load(sprintf('%sSw.txt', trajname)); % covariance

%% Plot trajectories
figure(1);
xmin = [0,0,0,0,0,0];xmax = [60,60,60,60,60,60];
ymin = [0.2,-0.5,-0.1,-1,-1,-1];ymax = [0.8,0.5,0.1,1,1,1];
for i=1:dof
    strPlot = ['23',num2str(i)]; % for dof=6
    subplot(strPlot);hold on;
    for k=1:Ndemos % Plot trajectories
        plot(demoY(:,i,k)); axis([xmin(i) xmax(i) ymin(i) ymax(i)]);
    end
end


%% Basis functions
n = Nt;           % number of basis functions
sigma = 1.5;    %variance
mu = linspace(0, Nt, n);

phi = zeros(n,Nt);
for i = 1:n 
    phi(i,:)=gaussBasis(time,mu(i),sigma); % one bf per row
    phi(i,:)=phi(i,:)/sum(phi(i,:));
end

%% Get weights
w = zeros(n,dof,Ndemos);
for k = 1:Ndemos %For each demostration
    y = demoY(:,:,k);
    %% Least squares
    weight=pinv(phi')*y;
    w(:,:,k)=weight;
end

% Mean
mu_w = mean(w,3);
% Standard deviation
std_w= std(w,0,3);
% % Covariance
% cov_w = (w-mu_w*ones(1,Ndemos))*(w-mu_w*ones(1,Ndemos))'/Ndemos;
% lamda = 10^(-3);
% cov_w = (cov_w +cov_w)/2 + eye(size(cov_w))*lamda;



mu_y = phi'*mu_w;  %Mean Traj
upper_ci = phi'*(mu_w+1.96*std_w/sqrt(n)) % Confidence Intervals (95%)
lower_ci = phi'*(mu_w-1.96*std_w/sqrt(n)) % Confidence Intervals (95%)

%% Plots
for i=1:dof
    strPlot = ['23',num2str(i)]; % for dof=6
    subplot(strPlot);hold on;
    
    % Basis functions
    %plot(time,phi,'g');

    % Weights
    %plot(mu,w/40,'r');

    % Trajectory
    plot(time,mu_y(:,i),'r','lineWidth',2);
    
    % CI
%     plot(time,upper_ci(:,i),'--r','lineWidth',2);
%     plot(time,lower_ci(:,i),'--r','lineWidth',2);
    hold off
end

%% Plot3D
figure(2);hold on;
for i=1:Ndemos
    trajX = reshape(demoY(:,1,:),Nt,Ndemos);
    trajY = reshape(demoY(:,2,:),Nt,Ndemos);
    trajZ = reshape(demoY(:,3,:),Nt,Ndemos);
    plot3(trajX,trajY,trajZ,'b');
end
plot3(mu_y(:,1),mu_y(:,2),mu_y(:,3),'r','lineWidth',2)


%% Comparative
% figure(1)
% load EM_traj.mat
% for i=1:dof
%     strPlot = ['23',num2str(i)]; % for dof=6
%     subplot(strPlot);hold on;
%     % Trajectory
%     plot(time,Y(:,i),'y','lineWidth',2);
%     plot(time,Ynew(:,i),'g','lineWidth',2);
% end
% figure(2)
% plot3(Y(:,1),Y(:,2),Y(:,3),'y','lineWidth',2)
% plot3(Ynew(:,1),Ynew(:,2),Ynew(:,3),'g','lineWidth',2)

% %% Plot Weights
% figure(3)
% i = 1;
% for i=1:dof
%     strPlot = ['23',num2str(i)]; % for dof=6
%     subplot(strPlot);hold on;
%     % Weights
%     for k=1:Ndemos
%     plot(mu,w(:,i,k),'r');
%     end
%     plot(mw(i:(i+9)),'b');
%     i = i+10;
% end

