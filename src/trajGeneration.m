function [trajT, trajX] =trajGeneration(num,position,time,delta_time)
tt = 0:delta_time:1;

trajX = [];
for i=1:num
    r = randn(1,length(position)); % add some differences
    x =  position +  r/10;
    s = spline(time,x,tt); % spline
    trajX = [trajX s'];
end
trajT = tt';

end

