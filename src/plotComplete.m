function plotComplete(mu_w,cov_w,mu_w_VP,cov_w_VP,phi,time,viaPoint,trajSampled)
figure(2);title('Random Trajectory created from the fit of the weights');hold on; grid on;
% Mean Trajectories
trajectory = phi'*mu_w;
plot(time,trajectory,'r');
trajectory = phi'*mu_w_VP;
plot(time,trajectory,'g');

% CI
[upper_ci,lower_ci] = confidenceIntervals(mu_w,cov_w,phi);
[upper_ci_VP,lower_ci_VP] = confidenceIntervals(mu_w_VP,cov_w_VP,phi);

% Sampled Trajectories
plot(time,trajSampled,'b','lineWidth',2);

% Plot VP
plot(time(viaPoint(1)),viaPoint(2),'xg','lineWidth',2,'MarkerSize',20);

% Initial CI's
fill_between_lines(time,upper_ci,lower_ci,[1 0 0],0.3)
plot(time,upper_ci,'r');
plot(time,lower_ci,'r');


% Via point CI's
fill_between_lines(time,upper_ci_VP,lower_ci_VP,[0 1 0],0.3)

plot(time,upper_ci_VP,'g');
plot(time,lower_ci_VP,'g');
