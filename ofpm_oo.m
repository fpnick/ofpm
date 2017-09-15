function [ A,rhs,sol,pointcloud,rho,solver ] = ofpm_oo( h, lbx,lby,ubx,uby,printlevel )
%OFPM_OO Summary of this function goes here
%   Detailed explanation goes here

solver=Solver(h,lbx,lby,ubx,uby,printlevel);
rho = solver.advance();

A=solver.matrices{1};
rhs=solver.rhss{1};
sol=solver.sol;
max(sol)
pointcloud = solver.pointcloud;

end

