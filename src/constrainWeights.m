function [mu_w_VP, cov_w_VP] = constrainWeights(mu_w, cov_w, phi, point, var)

t = point(1);% remake
y = point(2);
% Via Points
mu_w_VP = mu_w + cov_w*phi(:,t)*inv(var + phi(:,t)'*cov_w*phi(:,t))*(y-phi(:,t)' * mu_w);
cov_w_VP = cov_w - cov_w*phi(:,t)*inv(var + phi(:,t)'*cov_w*phi(:,t)) * phi(:,t)'*cov_w;
