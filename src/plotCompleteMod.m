function plotCompleteMod(mu_w,cov_w,phi,time,sampleTrajectory)
% figure(2);title('Random Trajectory created from the fit of the weights');hold on; grid on;
% Mean Trajectories
n = size(phi,1); % Number of basis function per joint

nJoint = size(sampleTrajectory,1);

for iJoint = 1:nJoint
    mu_w_joint =  mu_w((1+n*(iJoint-1)):(n*iJoint));
    cov_w_joint = cov_w((1+n*(iJoint-1)):(n*iJoint),(1+n*(iJoint-1)):(n*iJoint));
    [upper_ci_joint,lower_ci_joint] = confidenceIntervals(mu_w_joint,cov_w_joint,phi);

    
    figure(iJoint)
    clf
    % Sampled trajectory
    trajectoryMean = phi'*mu_w_joint;
    
    plot(time,trajectoryMean,'r', time,sampleTrajectory(iJoint,:),'g')
    hold on
    % Confidence interval
    fill_between_lines(time,upper_ci_joint,lower_ci_joint,[1 0 0],0.3)
    plot(time,upper_ci_joint,'r');
    plot(time,lower_ci_joint,'r');
    hold off

    title(sprintf('Joint %d',iJoint))
end
