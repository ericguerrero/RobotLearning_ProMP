function [Sampled_Traj]=sampleTrajectories(n,mu_w,cov_w,PHI,Nt,dof)

for k = 1:n
    w_traj = mvnrnd(mu_w,cov_w)';
    trajSampled(:,k) = PHI'*w_traj;
end
Sampled_Traj = reshape(trajSampled,Nt,dof,n);
