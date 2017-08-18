classdef Multigrid < handle
    %MULTIGRID represents a multigrid solver for the linear system.
    %   Detailed explanation goes here

    properties

      solver % Note that solver contains the complete hierarchy + matrices
      SMOOTHER = 1    % 1: Gauss-Seidel
      RESTRICTION = 1 % 1: Inclusion
      INTERPOLATION = 1 % 1: Weighed based on distance
      nPreSmooth = 0
      nPostSmooth = 0
    end

    methods
      function obj = Multigrid(solver)
         obj.solver = solver;
      end

      function solution = solve(obj,u,tol)

         res = norm( obj.solver.matrices{1}*u-obj.solver.rhss{1}, 2)
         tol_abs = res * tol
         iterations = 0;
         while ( res > tol_abs )
            u = obj.cycle( 1, u, obj.solver.rhss{1});
            res = norm( obj.solver.matrices{1}*u-obj.solver.rhss{1}, 2)
            iterations = iterations + 1;
         end

         solution = u;
         fprintf('Iterations required: %i\n', iterations);

      end
    end

    methods (Access=private)

      function u = cycle( obj, level, u, f )
         if ( level == obj.solver.hierarchy.depth )
            % Coarsest level => direct solve
            disp('Direct solve');
            u = obj.solver.matrices{level} \ f;
         else
            u              = obj.smooth( obj.solver.matrices{level}, u, f, obj.nPreSmooth);
            resVec         = obj.solver.matrices{level}*u-f;
            resVecCoarse   = obj.restrict( resVec, level);
            correction     = zeros(length(resVecCoarse),1);
            correction     = obj.cycle( level+1, correction, resVecCoarse);
            correctionFine = obj.interpolate( correction, level+1);
            u              = u - correctionFine;
            u              = obj.smooth( obj.solver.matrices{level}, u, f, obj.nPostSmooth);
         end
         
      end

      function u = smooth (obj,A, u, f, iter)
         if ( obj.SMOOTHER == 1 ) 
            L=tril(A,0);
            R=triu(A,1);
            u0 = u;
            for i=1:iter
              u = L\(f - R*u0);
              u0=u;
            end
         end
      end

      function R = restrict( obj, resVec, level)
         NCoarse = obj.solver.hierarchy.pointclouds{level+1}.N;
         R = zeros( NCoarse, 1);

         if ( obj.RESTRICTION == 1 )
            for i=1:NCoarse
               R(i) = resVec( obj.solver.hierarchy.coarse2fine{level+1}(i) );
            end
         end

      end

      function r = interpolate( obj, correction, level)
         % level refers to the level that is interpolated FROM
         % NOTE: Setup of interpolation should be done ONCE, not in every
         % iteration!
         NFine = obj.solver.hierarchy.pointclouds{level-1}.N;
         r = zeros (NFine, 1);

         if ( obj.INTERPOLATION == 1 )
            for i=1:NFine
               if ( obj.solver.hierarchy.fine2coarse{level}(i) == 0 )
                  neighbourList = obj.solver.hierarchy.pointclouds{level-1}.neighbourLists{i};
                  distanceList = obj.solver.hierarchy.pointclouds{level-1}.distanceLists{i};
                  nNeighbours = length(neighbourList);

                  sumDistances = 0;
                  for j=1:nNeighbours
                     if  ( obj.solver.hierarchy.fine2coarse{level}(neighbourList(j)) ~= 0 )
                        sumDistances = sumDistances + distanceList(j);
                     end
                  end
                  for j=1:nNeighbours
                     if  ( obj.solver.hierarchy.fine2coarse{level}(neighbourList(j)) ~= 0 )
                        r(i) = r(i) + correction(neighbourList(j)) * (distanceList(j)/sumDistances);
                     end
                  end
               else
                  r(i) = correction(obj.solver.hierarchy.fine2coarse{level}(i));
               end
            end
         end
      end
    end
end %class
