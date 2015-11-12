function writeCartesianTrajectory( name,trajectory ,time)
% This function writes to file a cartesian trajectory given as input with
% the name given as input and with the time[seconds] specified as input.
% The trajectory as to be 3xn where each row represents
% x y z respectively[meters].

if size(trajectory,1) ~=3 
    error('The trajectory passed to writeCartesianTrajectory function has not the correct format 3xn')
end

%open file
fileID = fopen(name,'w');

%specification of format
formatTime = '%f\n';
formatCoordinates = '%f %f %f\n';

%write
fprintf(fileID,formatTime,time);

for i=1:size(trajectory,2)-1
    fprintf(fileID,formatCoordinates,trajectory(:,i));
end
fprintf(fileID, '%f %f %f',trajectory(:,i+1)); %we don't want the newline at the last row

end

