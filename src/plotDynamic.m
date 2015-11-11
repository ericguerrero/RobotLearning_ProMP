function plotDynamic(mu_w,cov_w,phi,trajT)
% CI
[upper_ci,lower_ci] = confidenceIntervals(mu_w,cov_w,phi)

str='Y';
while((strcmp(str,'Y') || strcmp(str,'y')))
    str = input('New trajectory?[Y/N]','s')
    if strcmp(str,'Y') || strcmp(str,'y')
        w_traj = mvnrnd(mu_w,cov_w)';
        trajectory = phi'*w_traj;
        figure(2)
        plot(trajT,trajectory,'r','lineWidth',2);
    end
end
disp('_________EXIT________')

%%
% Dynamic plot
disp('SAMPLING TRAJECTORIES')

N_ITERATIONS = 100; % number of sampled trajectories

figure(3)
hold on; grid on;
plot(trajT,upper_ci,'--r','lineWidth',3);
plot(trajT,lower_ci,'--r','lineWidth',3);
fill_between_lines(trajT,upper_ci,lower_ci,[0 1 0],0.3)
title(sprintf('%d sampled trajectories',N_ITERATIONS))

for i=1:N_ITERATIONS
    i
    w_traj = mvnrnd(mu_w,cov_w)';
    trajectory = phi'*w_traj;
    figure(3)
    plot(trajT,trajectory,'r');
    pause(0.1)
end
