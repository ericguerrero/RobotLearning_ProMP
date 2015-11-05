function [ y ] = gaussBasis( x,mu,sigma )
% gaussian basis function: formula 2 paper
y = exp(-((x-mu).^2)./(2*sigma));
end

