classdef Multigrid < handle
    %MULTIGRID represents a multigrid solver for the linear system.
    %   Detailed explanation goes here

    properties

      solver % Note that solver contains the complete hierarchy + matrices
      smoother = 1
    end

    methods
      function obj = Multigrid(solver)
         obj.solver = solver;
      end

      function solution = solve(obj,u,tol)
         %solution = obj.solver.matrices{1} \ obj.solver.rhss{1};

         res = norm( obj.solver.matrices{1}*u-obj.solver.rhss{1}, 2);
         while ( res > tol )
            u = obj.smooth( obj.solver.matrices{1}, u,obj.solver.rhss{1}, 1);
            res = norm( obj.solver.matrices{1}*u-obj.solver.rhss{1}, 2);
         end

         solution = u;

      end
    end

    methods (Access=private)

      function u = smooth (obj,A, u, f, iter)
         if ( obj.smoother == 1 ) 
            L=tril(A,0);
            R=triu(A,1);
            u0 = u;
            for i=1:iter
              u = L\(f - R*u0);
              u0=u;
            end
         end
      end
    end



end %class
