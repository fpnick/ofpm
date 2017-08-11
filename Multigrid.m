classdef Multigrid < handle
    %MULTIGRID represents a multigrid solver for the linear system.
    %   Detailed explanation goes here

    properties
      solver
    end

    methods
      function obj = Multigrid(solver)
         obj.solver = solver;
      end

      function solution = solve(obj,u)
         solution = obj.solver.matrices{1} \ obj.solver.rhss{1};
      end
    end

    methods (Access=private)
    end



end %class
