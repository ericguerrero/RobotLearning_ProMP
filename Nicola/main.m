clc;close all;clear all; format compact;

%% Demostrations
num = 20; % # of demostrations
position = [1, 3.5, 3]; 
time = [0 30 60];  
[trajT, trajX] = trajGeneration(num,position,time);

%% Basis functions
n=60; % number of basis functions
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
 
%% Get weights
w = zeros(n,num);
for i = 1:num %For each demostration
    y = trajX(:,i);
    w(:,i) = pinv(phi')*y;
end

%% mean 
w_mean = mean(w,2);
mean_t = phi'*w_mean;

plot(trajT,mean_t,'LineWidth',2)

%% variance
% Standard deviation
std_w= std(w,0,2);
% Covariance
cov_w = (w-w_mean*ones(1,num))*(w-w_mean*ones(1,num))'/num;
lamda = 10^(-3);
cov_w = (cov_w +cov_w)/2 + eye(size(cov_w))*lamda;

%% Confidence Intervals (95%)
upper_ci = phi'*(w_mean+2*sqrt(diag(cov_w)))
lower_ci = phi'*(w_mean-2*sqrt(diag(cov_w)))

upper_ci_std_w = phi'*(w_mean+1.96*std_w/sqrt(n))
lower_ci_std_w = phi'*(w_mean-1.96*std_w/sqrt(n))

% CI
plot(trajT,upper_ci,'--r','lineWidth',2);
plot(trajT,lower_ci,'--r','lineWidth',2);

plot(trajT,upper_ci_std_w,'--b','lineWidth',2);
plot(trajT,lower_ci_std_w,'--b','lineWidth',2);
hold off