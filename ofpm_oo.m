function [ A,rhs,sol ] = ofpm_oo( h, printlevel )
%OFPM_OO Summary of this function goes here
%   Detailed explanation goes here

pointcloud = Pointcloud(h,0,0,1,1);
pointcloud.findNeighbours;
pointcloud.organize;
if ( printlevel > 9 )
    pointcloud.plot;
end

coarsecloud = pointcloud.coarsen;

pointcloud.plot;
figure;
coarsecloud.plot;


end

