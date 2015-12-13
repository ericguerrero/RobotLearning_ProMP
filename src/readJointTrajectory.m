function [j1,j2,j3,j4,j5,j6] = readJointTrajectory( name )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
name = 'joints_trajectory.txt';
fileID = fopen(name,'r');
tline=fgets(fileID)
J = textscan(fileID,'%f %f %f %f %f %f %f');
figure()
subplot(711);plot(J{1,1})
subplot(712);plot(J{1,2})
subplot(713);plot(J{1,3})
subplot(714);plot(J{1,4})
subplot(715);plot(J{1,5})
subplot(716);plot(J{1,6})
subplot(717);plot(J{1,7})



end

