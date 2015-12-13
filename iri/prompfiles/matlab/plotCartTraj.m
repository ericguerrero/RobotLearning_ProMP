% close all
% 
% foo=load(strcat('Cartesian0.txt')); %% just to know the size of the matrices
% Ndemos=3;
% for j=1:6
%     Amean{j}=zeros(size(foo,1),1);
% end
% Nt=size(foo,1); % timestep number
% WAM; % WAM MATLAB object
% for i=1:Ndemos
%          COL=[rand rand rand];
%         % load files Aj, each one is a demonstration
%         Aj=load(strcat('Cartesian',int2str(i-1),'.txt'));
%     
%     for j=1:6
%         figure(j)  
%         hold on
%         plot(Aj(:,j),'LineWidth',2,'Color',COL)
%         hold on
%         %mean
%         Amean{j}=Amean{j}+Aj(:,j)/Ndemos;
%     end    
% end
% 
% %now we get std
% STDev=[];
% 
% for j=1:6
% Aux=[];
%     for i=1:Ndemos
%        Aj=load(strcat('Cartesian',int2str(i-1),'.txt'));
%        Aux=[Aux Aj(:,j)];
%     end
%     STDev=[STDev std(Aux')'];
%     figure(j)
%     hold on
%         shadedplot((1:1:Nt),Amean{j}'-2*STDev(:,j)',Amean{j}'+2*STDev(:,j)','c','b')
%     hold on
%     plot(1:1:Nt,Amean{j},'b','LineWidth',2)
%     
% end
% 
% 
% for i=1:Ndemos
%     COL='k';%[rand rand rand];
%         Aj=load(strcat('Cartesian',int2str(i-1),'.txt'));
%     
%     for j=1:6
%         figure(j)  
%         hold on
%         plot(Aj(:,j),'LineWidth',1,'Color',COL)
%         
%         %mean
%         Amean{j}=Amean{j}+Aj(:,j)/Ndemos;
%     end    
% end
% 
% figure(1)
% title 'X'
% figure(2)
% title 'Y'
% figure(3)
% title 'Z'
% figure(4)
% title 'RX'
% figure(5)
% title 'RY'
% figure(6)
% title 'RZ'
close all;
clear;
dof=6;
Nt=118;

%trajpath = './trajs/C_trajectory';
%trajpath = './trajs/gerardProMptest';
trajpath = './trajs/gerardJustFeed';

trajname = 'newProMpParams';
%trajname = [trajpath '/gerardPROMPparamsTest']; %% C trajectory
%trajname = [trajpath '/gerardProMpfeed'];% full feed
%trajname = [trajpath '/justfeed'];% just feed;

mw=load(sprintf('%smw.txt', trajname)); % carreguem la mitja de la ProMP
Sw=load(sprintf('%sSw.txt', trajname)); % carreguem la Sw

Ndemos = 5;

Nf=size(mw,1)/dof;
C=(0:1:Nf-1)/Nf;%0:0.1:0.9;
D=0.025;
Y=[];
SY=[];
Sw=Sw+eye(size(Sw))*0.1;%1e-3; 

if(min(eig(Sw))<0)
   warning('REGULARITZAR MILLOR!') % Aka sumarli més a la Sw.
end
% change trajectory

%% punt desitjat -- JUGAR AQUI
points = []% [43];[99]; % 99 outpoint%[185, 44]; %% C traj
%points = [102]; %% Feed traj
cov_ratio = ones(size(points))*10;%[ 1 1 1];
ydesinc = [[+0.1 -0.11 +0.1]; %x, y, z
           [0.0 0.51 0.0];
           %[0.0 -0.51 0.0];
          ];
Swnew = Sw;
mwdes = mw;
for idx=1:length(points);
    i = points(idx);
    for s=1:Nf
        GT(1,s)=evalexp(i/Nt,C(s),D);
    end
    AUX=kron(eye(dof),GT);
    ydes=(AUX*mwdes);      
    Sydes_=AUX*Swnew*AUX';%+eye(size(Sydes))*1e-6;
    Sydes_=Sydes_+eye(size(Sydes_))*1e-6;
    
    % Modify position
    ydes(1:3) = ydes(1:3) + ydesinc(idx, :)';
    Sydes=Sydes_/cov_ratio(idx);

    Sw_old = Swnew;
    Swnew=Swnew*(eye(size(Swnew))-AUX'*pinv(Sydes+AUX*Swnew*AUX')*AUX*Swnew );
    mwdes=mwdes+Sw_old*AUX'*pinv(Sydes+AUX*Sw_old*AUX')*(ydes-AUX*mwdes);
end
% %     %Test non iterative
% %     i=44;
% %     for s=1:Nf
% %         GT(1,s)=evalexp(i/Nt,C(s),D);
% %     end
% %     AUX=kron(eye(dof),GT);
% %     ydes=(AUX*mw);     
% %     Sydes_=AUX*Sw*AUX';%+eye(size(Sydes))*1e-6;
% %     Sydes_=Sydes_+eye(size(Sydes_))*1e-6;
% %     ydes(2)=ydes(2)+0.3%+randn(6,1)*0.01;
% %     Sydes=Sydes_/10;
% %     
% %     %Sy_old=(AUX*Sw*AUX'+(AUX*Sw*AUX')')/2+eye(size(Sydes))*1e-6;    
% %     %p=mvnpdf(ydes,AUX*mw,Sy_old);
% %     %Sydes=eye(6)*0.1;
% %     Swnew=Sw*(eye(size(Sw))-AUX'*pinv(Sydes+AUX*Sw*AUX')*AUX*Sw );
% %     mwdes=mw+Sw*AUX'*pinv(Sydes+AUX*Sw*AUX')*(ydes-AUX*mw);
% %     
% %     %% condicionar un segon punt
% %     i=15;
% %     for s=1:Nf
% %         GT(1,s)=evalexp(i/Nt,C(s),D);
% %     end
% %     AUX=kron(eye(dof),GT);
% %     ydes=(AUX*mwdes);      
% %     Sydes_=AUX*Swnew*AUX';
% %     Sydes_=Sydes_+eye(size(Sydes_))*1e-6;
% %     
% %     ydes(1)=ydes(1)-0.1%+randn(6,1)*0.01;
% %     Sydes=Sydes_/4;
% %     Swnew=Swnew*(eye(size(Sw))-AUX'*pinv(Sydes+AUX*Swnew*AUX')*AUX*Swnew );
% %     mwdes=mwdes+Swnew*AUX'*pinv(Sydes+AUX*Swnew*AUX')*(ydes-AUX*mwdes);


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
    Ymean=[Ymean;(AUX*mw)'];      
    Syt=AUX*Sw*AUX';
    SYmean=[SYmean;2*real(sqrt(diag(Syt)'))];
end
WK=[];
mwnew=mw;
Swnew=Sw;
maxeigold=1000000;

for k=1:Ndemos
    traj = load(sprintf('%s/Cartesian%d.txt', trajpath, k-1));   
    demoY(:,:,k)=traj;
end

for iteration=1:100
    [iteration max(eig(Swnew))]
    if maxeigold<max(eig(Swnew))
       break 
    else
        [Om_new,mwnew,WK,Swnew,Sy,MK,Sk,GTpos]=ProMP_WMLE(demoY,dof,eye(dof),WK,mwnew,Swnew,eye(dof)*0.5,ones(5,1),GGT);
        [Om_new,mwnew,WK,Swnew,Sy,MK,Sk,GTpos]=ProMP_WMLE(demoY,dof,eye(dof),WK,mwnew,Swnew,eye(dof)*0.5,ones(5,1),GGT);

    end
end
Sw=Swnew;
mw=mwnew;

for j=1:dof
    %figure(j)
    subplot(2,3,j);
    plot(Ymean(:,j),'r','LineWidth',2);    
    hold on
    plot(Ymean(:,j)+SYmean(:,j),'r');    
    plot(Ymean(:,j)-SYmean(:,j),'r');    
    grid on
    
end

Y=[];
SY=[];
% SY2=[];
%close all
for i=1:Nt
    for s=1:Nf
        GT(1,s)=evalexp(i/Nt,C(s),D);
    end
    AUX=kron(eye(dof),GT);
    Y=[Y;(AUX*mwdes)'];      
    Syt=AUX*Swnew*AUX';
    SY=[SY;2*real(sqrt(diag(Syt)'))];
end

titles = {'Pose x', 'Pose y', 'Pose z', 'Orientation x', 'Orientation y', 'Orientation z'};
for j=1:dof
    %figure(j)
    subplot(2,3,j);
    title(titles(j));
    plot(Y(:,j),'g','LineWidth',2);    
    hold on
    plot(Y(:,j)+SY(:,j),'g');    
    plot(Y(:,j)-SY(:,j),'g');    
    grid on
end

%3D plot
sample = 1:4:Nt;
plotSYmean = 0;
plotOrient = 0;
plotSY = 0;
figure;
hold on;
% Mean trajectory
plot3(Ymean(:,1), Ymean(:,2), Ymean(:,3),'r','LineWidth',2);
scatter3(Ymean(1,1), Ymean(1,2), Ymean(1,3), 55, 'r', 'filled'); % Start point
scatter3(Ymean(points,1), Ymean(points,2), Ymean(points,3), 55, 'r', 'filled', 'MarkerEdgeColor','k'); % Conditioned points

if plotOrient
    quiver3(Ymean(sample,1), Ymean(sample,2), Ymean(sample,3), Ymean(sample,4), Ymean(sample,5), Ymean(sample,6), 0.2,'r','LineWidth',1.25);
end
if plotSYmean
    plot3(Ymean(:,1)+SYmean(:,1), Ymean(:,2)+SYmean(:,2), Ymean(:,3)+SYmean(:,3),'r','LineWidth',1);
    plot3(Ymean(:,1)-SYmean(:,1), Ymean(:,2)-SYmean(:,2), Ymean(:,3)-SYmean(:,3),'r','LineWidth',1);
end

% New trajectory
plot3(Y(:,1), Y(:,2), Y(:,3),'g','LineWidth',2);
scatter3(Y(1,1), Y(1,2), Y(1,3), 55, 'g', 'filled'); % Startpoint
scatter3(Y(points,1), Y(points,2), Y(points,3), 55, 'g', 'filled', 'MarkerEdgeColor','k'); % Conditioned points

if plotOrient
    quiver3(Y(sample,1), Y(sample,2), Y(sample,3), Y(sample,4), Y(sample,5), Y(sample,6), 0.2,'g','LineWidth',1.25);
end

if plotSY
    plot3(Y(:,1)+SY(:,1), Y(:,2)+SY(:,2), Y(:,3)+SY(:,3),'g','LineWidth',1);
    plot3(Y(:,1)-SY(:,1), Y(:,2)-SY(:,2), Y(:,3)-SY(:,3),'g','LineWidth',1);
end
hold off;
xlabel('x');ylabel('y');zlabel('z');
view(-96, 3); % For correct C trajectory view


% dw0=[0.1097 -0.1237 0.07135 0.05655 -0.219 0.3005 -0.1745 0.1109 0.02812 0.03921 0.01726 0.1862 -0.1902 0.2773 -0.2162 0.1028 -0.09207 -0.03394 0.00163 0.0003357 0.04676 -0.1079 -0.1105 0.2178 -0.3071 0.224 -0.1833 0.2698 -0.1546 -0.03853 0.03596 -0.00517 -0.135 0.3001 -0.3073 0.2444 -0.2893 0.2828 -0.08483 0.1381 -0.01872 0.05331 -0.1423 0.1583 -0.3435 0.1828 -0.2334 0.07542 -0.09716 -0.001082 -0.146 0.1711 -0.4182 0.4706 -0.5474 0.5103 -0.3233 0.1627 -0.00289 -0.1105]
% 
% 
%  dw0=[0.1673 -0.2932 0.1563 0.2106 -0.543 0.7152 -0.6371 0.4142 -0.1815 0.03586 -0.2042 0.5115 -0.6244 0.5771 -0.3719 0.1467 0.035 -0.1117 0.06715 -0.005538 0.0892 -0.1391 -0.03693 0.4593 -0.8451 1.076 -1.07 0.8488 -0.4416 0.1179 0.1962 -0.2015 0.2228 -0.1091 0.09359 -0.0219 0.04154 0.02471 -0.05118 0.09232 0.1038 -0.2287 0.4351 -0.4726 0.3268 -0.001686 -0.2207 0.3232 -0.1861 0.07886 0.1239 -0.2277 0.4032 -0.4784 0.535 -0.494 0.3682 -0.2045 0.1016 0.0189]
% 
% 
% % dw0=dw0*3;
% % dw=mvnrnd2(0*mw,Sw+eye(size(Sw))*0.001);
% Y2=[];
% Y=[];
% SY=[]
% SY2=[]
% close all
% for i=1:Nt
%     for s=1:Nf
%         GT(1,s)=evalexp(i/Nt,C(s),D);
%     end
%     AUX=kron(eye(6),GT);
%     Y=[Y;(AUX*mw)'];      
%     Y2=[Y2;(AUX*(mwdes+mw*0+dw0'*0))'];   
%     Syt=AUX*Sw*AUX';
%     Syt2=AUX*Swnew*AUX';
%     SY=[SY;2*real(sqrt(diag(Syt)'))];
%     SY2=[SY2;2*real(sqrt(diag(Syt2)'))];
% end
% for j=1:6
%     figure(j)
%     plot(Y(:,j),'b','LineWidth',2);    
%     hold on
%     plot(Y2(:,j),'g','LineWidth',2);  
%     plot(Y(:,j)+SY(:,j),'r');    
%     plot(Y(:,j)-SY(:,j),'r');    
%     plot(Y2(:,j)+SY2(:,j),'m');    
%     plot(Y2(:,j)-SY2(:,j),'m'); 
%     grid on
% end
% 
% 
% % a=load('fulltest_state_sequence.txt');
% % for j=1:7
% %     figure(j)
% %     hold on
% %     plot(a(:,1),a(:,j+1),'g','LineWidth',3)
% % end

% 
% DE=load('test_state_sequence.txt');
% 
% XDE=convert2cartesian(DE(:,2:8));
% XDE2=convert2cartesian(DE(:,9:15));
% XDE(:,4:6)=XDE(:,4:6);
% for j=1:6
%     figure(j)
%    hold on
%   plot(DE(:,1)/8*80,XDE(:,j),'--g','Linewidth',2)
%   plot(DE(:,1)/8*80,XDE2(:,j),'--y','Linewidth',2)
% end


%% Plot all the trajectories
for j=1:dof 
    figure(j+3)
    for i = 0:Ndemos-1
                traj = load(sprintf('%s/Cartesian%d.txt', trajpath, i));
            plot(traj(:,j),'b', 'LineWidth', 1);
            hold on
            plot(Ymean(:,j),'k', 'LineWidth', 2)
    end
end



alltogether = 0;
if alltogether
    figure; hold on;
    for i = 0:Ndemos-1
        traj = load(sprintf('%s/Cartesian%d.txt', trajpath, i));
        plot3(traj(:,1), traj(:,2), traj(:,3), 'LineWidth', 2);
    end
    hold off;
else
    figure;
    for i = 0:Ndemos-1
        subplot(2, 3, i+1);
        traj = load(sprintf('%s/Cartesian%d.txt', trajpath, i));
        plot3(traj(:,1), traj(:,2), traj(:,3), 'LineWidth', 2);
        axis equal
    end
end