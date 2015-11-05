function [trajT, trajX] =trajGeneration(num,position,time)

trajX = [];
for i=1:num
    r = randn(1,3); % add some differences
    coord = [time; position + position .* r/10];
    s = fnplt(cscvn(coord)); % spline
    trajX = [trajX s(2,:)'];
end
trajT = s(1,:)';
end

