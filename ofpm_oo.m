function [ A,rhs,sol,coords ] = ofpm_oo( h, lbx,lby,ubx,uby,printlevel )
%OFPM_OO Summary of this function goes here
%   Detailed explanation goes here

solver=Solver(h,lbx,lby,ubx,uby,printlevel);
solver.advance();

A=solver.matrix;
rhs=solver.rhs;
sol=solver.sol;
coords = solver.pointcloud.coords;

%pointcloud.plot;
%figure;
%coarsecloud.plot;

% TODO Write a class "solver" that uses the pointcloud(s) to set up the
% matrix(matrices)



end

