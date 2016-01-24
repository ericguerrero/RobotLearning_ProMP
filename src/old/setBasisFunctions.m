function [phi] = setBasisFunctions(n,sigma,time)
mu = linspace(0,1,n);
phi = zeros(n,length(time));
for i = 1:n 
    phi(i,:)=gaussBasis(time,mu(i),sigma); % one bf per row
    phi(i,:)=phi(i,:)/sum(phi(i,:));
end