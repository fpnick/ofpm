function [ A,rhs,sol,pointcloud ] = ofpm_oo( h, lbx,lby,ubx,uby,printlevel )
%OFPM_OO Summary of this function goes here
%   Detailed explanation goes here

solver=Solver(h,lbx,lby,ubx,uby,printlevel);
solver.advance();

A=solver.matrices{1};
rhs=solver.rhss{1};
sol=solver.sol;
pointcloud = solver.pointcloud;

%pointcloud.plot;
%figure;
%coarsecloud.plot;

% TODO Write a class "solver" that uses the pointcloud(s) to set up the
% matrix(matrices)



end

