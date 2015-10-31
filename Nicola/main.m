clc;close all;clear all; format compact;

%% Demostrations
num = 100; % # of demostrations
position = [1, 2.5, 3]; 
time = [0 30 60];  
[trajT, trajX] = trajGeneration(num,position,time);

%% Basis functions
n=100; % number of basis functions
sigma = 1.5; %variance
mu = linspace(0,max(time),n);

phi = zeros(n,length(trajT));
for i = 1:n 
    phi(i,:)=gaussBasis(trajT,mu(i),sigma); % one bf per row
    %phi(i,:)=phi(i,:)/sum(phi(i,:));
end

% Demostrations
plot(trajT,trajX,'k');hold on;grid on;

% Basis functions
plot(trajT,phi,'g');
