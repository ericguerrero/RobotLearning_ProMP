function [w, mu_w, cov_w] = getWeights(trajDemo,phi)

n = size(phi,1);
num = size(trajDemo,2);

w = zeros(n,num);
for i = 1:num %For each demostration
    %% Least squares
    w(:,i) = pinv(phi')*trajDemo(:,i);
end

mu_w = mean(w,2);
cov_w = (w-mu_w*ones(1,num))*(w-mu_w*ones(1,num))'/num;
lamda = 10^(-9);
cov_w = (cov_w +cov_w)/2 + eye(size(cov_w))*lamda;
end