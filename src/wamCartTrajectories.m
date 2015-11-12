y = 0.5:-0.005:0-0.3;
x= sqrt((0.7*ones(1,length(y))).^2 - y.^2);
z = y./4;
%x =  0.5*ones(1,length(y));

trajectory =[x;y;z]

close all
plot3(x,y,z)
grid on;
xlabel('x')
ylabel('y')
zlabel('z')

writeCartesianTrajectory('ct.txt',trajectory,25)
