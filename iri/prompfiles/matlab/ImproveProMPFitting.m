close all;
clear;

%trajpath = './trajs/C_trajectory';
%trajpath = './trajs/gerardProMptest';
%trajpath = './trajs/gerardJustFeed';% this is the path to the folder where the files are located
trajpath = './test2';

trajname = 'test';% name of the file to load (with out mw or Sw)
%trajname = [trajpath '/gerardPROMPparamsTest']; %% C trajectory
%trajname = [trajpath '/gerardProMpfeed'];% full feed
%trajname = [trajpath '/justfeed'];% just feed;

Ndemos = 10; % <<<<< IMPORTANT: number of the demos = number of cartesian trajectories
for k=1:Ndemos
    traj = load(sprintf('%s/Cartesian%d.txt', trajpath, k-1));   
    demoY(:,:,k)=traj; % (time, dof, demo)
end
demoY=demoY(1:size(demoY,1),:,:);
Nt=size(demoY,1);
dof=size(demoY,2);

%here it is loading the weights and the covariance
mw=load(sprintf('%smw.txt', trajname)); % carreguem la mitja de la ProMP
Sw=load(sprintf('%sSw.txt', trajname)); % carreguem la Sw

Nf=size(mw,1)/dof; % AL LORO! el reproduce ProMP agafa aix?? d??un arxiu
C=(0:1:Nf-1)/Nf;%0:0.1:0.9;
D=0.025;
Y=[];
SY=[];
Sw=Sw+eye(size(Sw))*0.1;%1e-3; 

if(min(eig(Sw))<0)
   warning('REGULARITZAR MILLOR!') % Aka sumarli m??s a la Sw.
end

%% PLOT RESULTATS
Ymean=[];
SYmean=[];
% SY2mean=[];
GGT=[];
for i=1:Nt
    for s=1:Nf
        GT(1,s)=evalexp(i/Nt,C(s),D);
    end
    AUX=kron(eye(dof),GT);
    GGT=[GGT;AUX];
end
WK=[];
mwnew=mw;
Swnew=Sw;
maxeigold=1000000;


 L=1000000
 LL=[];
for iteration=1:100
%      [iteration max(eig(Swnew))]
      [iteration L]
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

%% some plots - they show the old fitting and the new one

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
% Plot

figure;
for i = 1:dof
    subplot(2,3,i);
%     figure
    hold on;
    for d = 1:Ndemos
        plot(demoY(:,i,d),'k', 'LineWidth', 1);
    end
    plot(Y(:,i),'b', 'LineWidth', 2);
%     plot(Y(:,i)+SY(:,i),'b', 'LineWidth', 1);
%     plot(Y(:,i)-SY(:,i),'b', 'LineWidth', 1);
    
    plot(Ynew(:,i),'r','LineWidth', 2);
%     plot(Ynew(:,i)+SYnew(:,i),'r','LineWidth', 1);
%     plot(Ynew(:,i)-SYnew(:,i),'r','LineWidth', 1);
    
end
csvwrite(sprintf('%s_new_mw.txt', trajname),mwnew);
csvwrite(sprintf('%s_new_Sw.txt', trajname),Swnew);
