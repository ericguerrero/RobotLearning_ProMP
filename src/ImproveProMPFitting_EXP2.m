close all;
clear;
%% Load files (specify Ndemos)
trajpath = '../experiments/avoidance_ball_01';
addpath(trajpath);
trajname = 'test';  
Ndemos = 5;   
[demoY,Nt, dof, x, mw, Sw] = loadTrajectory(trajpath, trajname, Ndemos,'shifted');



% Basis Functions
Nf=size(mw,1)/dof;  % # of basis functions
C=(0:1:Nf-1)/Nf;    % Centers;
D=0.025;            % Amplitud
Y=[];
SY=[];
Sw=Sw+eye(size(Sw))*0.1;    %1e-3; Regularization

if(min(eig(Sw))<0)
   warning('REGULARITZAR MILLOR!') % Increase value for the regularization
end
% Ndemos=20;
% DEMOS = zeros(Nt,dof,Ndemos);
% DEMOS(:,:,1:5) = demoY0;
% DEMOS(:,:,6:10) = demoY0.*(rand(size(demoY0))*0.01+1);
% DEMOS(:,:,11:15) = demoY0.*(rand(size(demoY0))*0.01+1);
% DEMOS(:,:,16:20) = demoY0.*(rand(size(demoY0))*0.01+1);
% demoY = DEMOS;

%% Expectation Maximization
% Get basis functions
Ymean=[];
SYmean=[];
% SY2mean=[];
GGT=[];
for i=1:Nt      % trajectory points
    for s=1:Nf  % basis functions
        GT(1,s)=evalexp(i/Nt,C(s),D);
    end
    AUX=kron(eye(dof),GT);
    GGT=[GGT;AUX];
end
WK=[];
mwnew=mw;
Swnew=Sw;
maxeigold=1000000;




% Expectation Maximization
 L=1000000
 LL=[];
for iteration=1:100
%      [iteration max(eig(Swnew))]
      [iteration L];
%       if maxeigold<max(eig(Swnew))
%          break 
%       else
	  % compute the weights and the covariance with Expectation Maximization
          [Om_new,mwnew,WK,Swnew,Sy,MK,Sk,GTpos]=ProMP_WMLE(demoY,dof,eye(dof),WK,mwnew,Swnew,eye(dof)*0.5,ones(Ndemos,1),GGT);
          L=Weightedlikelihood(Swnew,Sk,Sy,demoY,eye(dof),MK,GTpos,mwnew,ones(Ndemos,1));
          LL=[LL;L];
%       end
%       maxeigold=max(eig(Swnew));
end
Swnew;
mwnew;


%% Sampling
Nsample = 5
for k = 1:Nsample
    w_traj = mvnrnd(mwnew,Swnew*0.1)';
    trajSampled(:,k) = GGT*w_traj;
end
Sampled_Traj = reshape(trajSampled,Nt,dof,Nsample);

%% Plots
% Compute trajectories
Y=[];
SY=[];
for i=1:Nt
    for s=1:Nf
        GT(1,s)=evalexp(i/Nt,C(s),D);
    end
    AUX=kron(eye(dof),GT);
    Y=[Y;(AUX*mw)'];      
    Syt=AUX*Sw*AUX';
    SY=[SY;2*real(sqrt(diag(Syt)'))];
end


Ynew=[];
SYnew=[];
for i=1:Nt
    for s=1:Nf
        GT(1,s)=evalexp(i/Nt,C(s),D);
    end
    AUX=kron(eye(dof),GT);
    Ynew=[Ynew;(AUX*mwnew)'];      
    Syt=AUX*Swnew*AUX';
    SYnew=[SYnew;2*real(sqrt(diag(Syt)'))];
end


Nsamples = 5;
SampledTraj = zeros(Nt,dof,Nsamples);
for k=1:Nsamples
    Ysampled=[];
    w_traj = mvnrnd(mwnew,Swnew)'
    for i=1:Nt
        for s=1:Nf
            GT(1,s)=evalexp(i/Nt,C(s),D);
        end
        AUX=kron(eye(dof),GT);
        Ysampled = [Ysampled;(AUX*w_traj)'];      
    end
    SampledTraj(:,:,k) = Ysampled;
end


%% Plot
% Demos
for i = 1:dof
    figure(i);
%     figure
    hold on;
    for d = 1:Ndemos
        plot(demoY(:,i,d),'b', 'LineWidth', 1);
    end
    fill_between_lines(x',u_ci(:,k),l_ci(:,k),[1 0 0],0.2)
    plot(Y(:,i),'r', 'LineWidth', 2);
    plot(Y(:,i)+SY(:,i),'r', 'LineWidth', 1);
    plot(Y(:,i)-SY(:,i),'r', 'LineWidth', 1);

end
for i = 1:dof
    figure(i);
%     figure
    hold on;
    for s = 1:Nsamples
        plot(SampledTraj(:,i,s),'g', 'LineWidth', 1);
    end
    plot(Y(:,i),'b', 'LineWidth', 2);
    plot(Y(:,i)+SY(:,i),'b', 'LineWidth', 1);
    plot(Y(:,i)-SY(:,i),'b', 'LineWidth', 1);
    
    plot(Ynew(:,i),'r','LineWidth', 2);
    plot(Ynew(:,i)+SYnew(:,i),'r','LineWidth', 1);
    plot(Ynew(:,i)-SYnew(:,i),'r','LineWidth', 1);

end



csvwrite([trajpath,'/', sprintf('%s_new_mw.txt', trajname)],mwnew);
csvwrite([trajpath,'/', sprintf('%s_new_Sw.txt', trajname)],Swnew);
csvwrite([trajpath,'/', sprintf('%s_old_mw.txt', trajname)],mw);
csvwrite([trajpath,'/', sprintf('%s_old_Sw.txt', trajname)],Sw);



rmpath(trajpath);
 