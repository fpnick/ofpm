function [ A,rhs,sol,pointcloud,rho,iter_needed,solver ] = ofpm_oo( h, lbx,lby,ubx,uby,printlevel )
%OFPM_OO Summary of this function goes here
%   Detailed explanation goes here

addpath('./sort_back')

solver=Solver(h,lbx,lby,ubx,uby,printlevel);
[rho,iter_needed] = solver.advance();

A=solver.matrices{1};
rhs=solver.rhss{1};
sol=solver.sol;
max(sol)
pointcloud = solver.pointcloud;

end

