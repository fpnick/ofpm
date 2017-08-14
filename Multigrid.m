classdef Multigrid < handle
    %MULTIGRID represents a multigrid solver for the linear system.
    %   Detailed explanation goes here

    properties

      solver % Note that solver contains the complete hierarchy + matrices
      smoother = 1
      nPreSmooth = 1
      nPostSmooth = 1
    end

    methods
      function obj = Multigrid(solver)
         obj.solver = solver;
      end

      function solution = solve(obj,u,tol)

         res = norm( obj.solver.matrices{1}*u-obj.solver.rhss{1}, 2);
         while ( res > tol )
            u = obj.cycle( 1, u, obj.solver.rhss{1});
            res = norm( obj.solver.matrices{1}*u-obj.solver.rhss{1}, 2);
         end

         solution = u;

      end
    end

    methods (Access=private)

      function u = cycle( obj, level, u, f )
         if ( level == obj.solver.hierarchy.depth )
            % Coarsest level => direct solve
            disp('Direct solve');
            u = obj.solver.matrices{level} \ f;
         else
            u = obj.smooth( obj.solver.matrices{level}, u, f, obj.nPreSmooth);
            res = norm( obj.solver.matrices{level}*u-f, 2);
            res = obj.restrict( res, level);
         end
         
      end

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

      function R = restrict( obj, res, level)
         R = zeros( obj.solver.hierarchy.pointclouds{level+1}.N, 1);
      end
    end



end %class
