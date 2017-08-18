classdef Multigrid < handle
    %MULTIGRID represents a multigrid solver for the linear system.
    %   Detailed explanation goes here

    properties

      solver % Note that solver contains the complete hierarchy + matrices

      % Parameters
      SMOOTHER      = 1  % 1: Gauss-Seidel
      RESTRICTION   = 1  % 1: Inclusion
      INTERPOLATION = 1  % 1: Weighed based on distance
      nPreSmooth    = 1  % n: Number of pre-smoothing steps
      nPostSmooth   = 1  % n: Number of post-smoothing steps
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

      function obj = Multigrid(solver)
         obj.solver = solver;
      end

      function solution = solve(obj,u,tol)
         res        = norm( obj.solver.matrices{1}*u-obj.solver.rhss{1}, 2);
         res0       = res;
         tol_abs    = res * tol;
         iterations = 0;

         fprintf('Initial residual: %1.3e\n', res);
         
         while ( res > tol_abs )
            u          = obj.cycle( 1, u, obj.solver.rhss{1});
            res        = norm( obj.solver.matrices{1}*u-obj.solver.rhss{1}, 2);
            iterations = iterations + 1;

            fprintf('Residual after iteration %i: %1.3e\n', iterations,res);
         end

         solution = u;
         rho      = (res/res0)^(1/iterations);

         fprintf('rho = %1.3f\n', rho);
      end

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access=private)

      function u = cycle( obj, level, u, f )
         if ( level == obj.solver.hierarchy.depth )
            % Coarsest level => direct solve
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
         R       = zeros( NCoarse, 1);

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
         NFine       = obj.solver.hierarchy.pointclouds{level-1}.N;
         fine2coarse = obj.solver.hierarchy.fine2coarse{level};
         r           = zeros (NFine, 1);

         if ( obj.INTERPOLATION == 1 )
            for i=1:NFine
               if ( fine2coarse(i) == 0 )
                  neighbourList = obj.solver.hierarchy.pointclouds{level-1}.neighbourLists{i};
                  distanceList  = obj.solver.hierarchy.pointclouds{level-1}.distanceLists{i};
                  nNeighbours   = length(neighbourList);

                  sumDistances = 0;
                  for j=1:nNeighbours
                     if ( fine2coarse(neighbourList(j)) ~= 0 )
                        sumDistances = sumDistances + distanceList(j);
                     end
                  end
                  for j=1:nNeighbours
                     if ( fine2coarse(neighbourList(j)) ~= 0 )
                        r(i) = r(i) + correction(fine2coarse(neighbourList(j))) * (distanceList(j)/sumDistances);
                     end
                  end
               else
                  r(i) = correction(fine2coarse(i));
               end
            end
         end
      end
      
    end

end %class
