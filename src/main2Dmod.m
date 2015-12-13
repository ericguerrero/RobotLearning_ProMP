clc;close all;clear all; format compact;
%% ProMP Project
% trajectoryGeneration(params,z(t))
% phi(params,z(t))
% weightCalculation(phi,trajectorySet)
% sampleTrajectory(weights,CovWeights)
% trajectoryRepresentation(z(t),trajectorySet, sampledTrajectory,weights, CovWeights)

% Set Demostrations
num = 25; % # of demostrations
coords1 = [0 .4 .6 .8 1 ; 1 2.5 2 0.5 -1]; 
coords2 = [0 .4 .6 .8 1 ; 0.2 1 0.5 1.5 2]; 
coords3 = [0 .4 .6 .8 1 ; rand(1,5)]; 
coords4 = [0 .4 .6 .8 1 ; rand(1,5)]; 
coords5 = [0 .4 .6 .8 1 ; rand(1,5)]; 
coords6 = [0 .4 .6 .8 1 ; rand(1,5)]; 

m = 1000;
time = linspace(0,1,m)';


demoJoint1 = setDemonstrations(num,coords1,time,m);
demoJoint2 = setDemonstrations(num,coords2,time,m);
demoJoint3 = setDemonstrations(num,coords3,time,m);
demoJoint4 = setDemonstrations(num,coords4,time,m);
demoJoint5 = setDemonstrations(num,coords5,time,m);
demoJoint6 = setDemonstrations(num,coords6,time,m);

demoTrajectory = [demoJoint1';demoJoint2';demoJoint3';demoJoint4';demoJoint5';demoJoint6'];

Ys= [demoJoint1;demoJoint2;demoJoint3;demoJoint4;demoJoint5;demoJoint6;];
size(Ys)
nJoint = size(Ys,1)/m;

% Basis functions
n=30; % number of basis functions
sigma = 0.001; %variance
phi = setBasisFunctions(n,sigma,time);

size(phi)
% Kronecker product (multiple joints)
kronPhi = kron(eye(nJoint),phi)

size(kronPhi')
%%
% Get weights
[w, mu_w, cov_w] = getWeights(Ys,kronPhi);

size(mu_w)
%% Sample trajectory
% Sample Trajectory
w_traj = mvnrnd(mu_w,cov_w)';
cov_traj = cov_w;

sampleTrajectory = zeros(nJoint,m)
for iJoint = 1:nJoint
    w_traj_joint =  w_traj((1+n*(iJoint-1)):(n*iJoint));
    sampleTrajectory(iJoint,:) = (phi'*w_traj_joint)';
    

    figure(iJoint)
    hold on
%     subplot(nJoint,1,iJoint)
    %     subplot(nJoint/2,nJoint/(nJoint/2),iJoint)
    plot(time,demoTrajectory((1+(iJoint-1)*num):iJoint*num,:),'b')%,time,sampleTrajectory(iJoint,:),'b')
    title(sprintf('Joint %d',iJoint))
    hold off
end

%%
plotCompleteMod(w_traj,cov_traj,phi,time,sampleTrajectory)

% plot(time,sampleTrajectory(:,)

%%
% Plot Random trajectory with CI's and VP
plotComplete(mu_w,cov_w,mu_w_VP,cov_w_VP,phi,time,viaPoint,trajSampled)


%%
% Via Points
viaPoint =[0.8*m, 0.9];
viaPoint_var = 0.01;
[mu_w_VP, cov_w_VP] = constrainWeights(mu_w, cov_w, phi, viaPoint, viaPoint_var);


% Dynamic Plot
% plotDynamic(mu_w,cov_w,phi,time)
