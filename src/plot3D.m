function []=plot3D(Y,color)

trajX = reshape(Y(:,1,:),size(Y,1),size(Y,3));
trajY = reshape(Y(:,2,:),size(Y,1),size(Y,3));
trajZ = reshape(Y(:,3,:),size(Y,1),size(Y,3));
plot3(trajX,trajY,trajZ,color);
