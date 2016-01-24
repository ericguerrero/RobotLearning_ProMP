function [upper_ci,lower_ci] = confidenceIntervals(mu,cov,phi)
upper_ci = phi'*(mu+2*sqrt(diag(cov)));
lower_ci = phi'*(mu-2*sqrt(diag(cov)));