function [time, trajDemo] = setDemonstrations(num,points,samples)
time = linspace(0,1,samples)';
n = size(points,2);
trajDemo = zeros(samples,num);
for i=1:num
    % Add noise
    r = randn(1,n); 
    x =  points(2,:) +  r/10;
    
    %Create Demonstration
    trajDemo(:,i) = spline(points(1,:),x,time); % spline
end
end

