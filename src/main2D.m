clc;close all;clear all; format compact;
%% ProMP Project
% trajectoryGeneration(params,z(t))
% phi(params,z(t))
% weightCalculation(phi,trajectorySet)
% sampleTrajectory(weights,CovWeights)
% trajectoryRepresentation(z(t),trajectorySet, sampledTrajectory,weights, CovWeights)

% Set Demostrations
num = 20; % # of demostrations
coordsX = [0 .4 .6 .8 1 ; 1 2.5 2 0.5 -1]; 
coordsY = [0 .4 .6 .8 1 ; 0.2 1 0.5 1.5 2]; 
samples = 1000;
time = linspace(0,1,samples)';
coordDemoX = setDemonstrations(num,coordsX,time,samples);
coordDemoY = setDemonstrations(num,coordsY,time,samples);

% Basis functions
n=20; % number of basis functions
sigma = 0.001; %variance
phi = setBasisFunctions(n,sigma,time);

% Get weights
[w, mu_w, cov_w] = getWeights(coordDemoX,phi);

% Via Points
viaPoint =[0.8*samples, 0.9];
viaPoint_var = 0.01;
[mu_w_VP, cov_w_VP] = constrainWeights(mu_w, cov_w, phi, viaPoint, viaPoint_var);

% Sample Trajectory
w_traj = mvnrnd(mu_w_VP,cov_w_VP)';
trajSampled = phi'*w_traj;

% Plot Random trajectory with CI's and VP
plotComplete(mu_w,cov_w,mu_w_VP,cov_w_VP,phi,time,viaPoint,trajSampled)


% Dynamic Plot
% plotDynamic(mu_w,cov_w,phi,time)
