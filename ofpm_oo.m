function [ A,rhs,sol ] = ofpm_oo( h, printlevel )
%OFPM_OO Summary of this function goes here
%   Detailed explanation goes here

solver=Solver(h,printlevel);
solver.advance();

A=solver.matrix;
rhs=solver.rhs;
sol=solver.sol;

%pointcloud.plot;
%figure;
%coarsecloud.plot;

% TODO Write a class "solver" that uses the pointcloud(s) to set up the
% matrix(matrices)



end

