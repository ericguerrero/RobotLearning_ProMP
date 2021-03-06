function [coordDemo] = setDemonstrations(num,coord,time,samples)

n = size(coord,2);
coordDemo = zeros(samples,num);
for i=1:num
    % Add noise
    r = randn(1,n); 
    x =  coord(2,:) +  r/10;
    
    %Create Demonstration
   coordDemo(:,i) = spline(coord(1,:),x,time); % spline
end
end