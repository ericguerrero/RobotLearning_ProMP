load Experiment1.mat



for i=1:size(traj,3)
    name = sprintf('demoTraj%d.txt',i-1);
    writeCartesianTrajectory(name,traj(:,:,i)',25)
end

for i=1:size(sampleTrajectory,3)
    name = sprintf('sampleTraj%d.txt',i-1);
    writeCartesianTrajectory(name,sampleTrajectory(:,:,i)',25)
end

for i=1:size(meanTrajectory,3)
    name = sprintf('mean.txt');
    writeCartesianTrajectory(name,meanTrajectory(:,:,i)',25)
end