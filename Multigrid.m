classdef Multigrid < handle
    % MULTIGRID   represents a multigrid solver for the linear system.
    %    Represents the core multigrid solver. Assumes though, that the hierachy
    %    has already been created in the 'solver' object. 

    properties

      solver % Note that solver contains the complete hierarchy + matrices

      % Parameters
      SMOOTHER      = 1  % 1: Gauss-Seidel
      RESTRICTION   = 1  % 1: Injection
      INTERPOLATION = 1  % 1: Weighed based on distance
      NORMALIZE     = 0  % 1: Normalize matrices on every level
      ENFORCE_DIAGDOM = 0 % 1: Add 5% to every diagonal
      nPreSmooth    = 2  % n: Number of pre-smoothing steps
      nPostSmooth   = 2 % n: Number of post-smoothing steps
      nMaxIter      = 1
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

      function obj = Multigrid(solver)
         obj.solver = solver;
      end

      function [solution,rho] = solve(obj,u,tol)
      % SOLVE  Solve the linear system given by the hierachy in 'solver'
      %     solution = solve(u,tol) solves with a residual reduction of tol
      %                             starting with initial guess u.
         DEBUGLEVEL = 11;

         res        = norm( obj.solver.matrices{1}*u-obj.solver.rhss{1}, 2);
         res0       = res;
         tol_abs    = res * tol;
         iterations = 0;

         if ( obj.NORMALIZE == 1 ) 
            for i=1:obj.solver.hierarchy.depth
               obj.solver.matrices{i} = obj.solver.matrices{i}./diag(obj.solver.matrices{i});
            end
         end
         if ( obj.ENFORCE_DIAGDOM == 1 ) 
            for i=1:obj.solver.hierarchy.depth
               obj.solver.matrices{i} = obj.solver.matrices{i}+0.05*diag(diag(obj.solver.matrices{i}));
            end
         end

         fprintf('Initial residual: %1.3e\n', res);
         
         while ( res > tol_abs && iterations < obj.nMaxIter )
            u          = obj.cycle( 1, u, obj.solver.rhss{1});
            resVec     = obj.solver.matrices{1}*u-obj.solver.rhss{1};
            res        = norm( resVec, 2);
            iterations = iterations + 1;

            if DEBUGLEVEL>=10
               condition = condest(obj.solver.matrices{1});
               fprintf('Condition is %1.3e\n', condition);
               obj.solver.plotSolution(obj.solver.hierarchy.pointclouds{1},resVec,sprintf('Residual on level 1'));
            end

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
      % CYCLE  Performe one cycle from level downards.
      %     u = u(level,u,f)  Perform one cycle from level downwards using the
      %                       RHS f and the initial guess u.
         DEBUGLEVEL = 1;

         if ( level == obj.solver.hierarchy.depth )
            % Coarsest level => direct solve
            u = obj.solver.matrices{level} \ f;
         else
            if DEBUGLEVEL>=10
               solution = obj.solver.matrices{level} \ f;
               err = u-solution;
               obj.solver.plotSolution(obj.solver.hierarchy.pointclouds{level},err,sprintf('Error on level %i', level));
            end
            u              = obj.smooth( obj.solver.matrices{level}, u, f, obj.nPreSmooth);
            if DEBUGLEVEL>=10
               err = u-solution;
               obj.solver.plotSolution(obj.solver.hierarchy.pointclouds{level},err,sprintf('Presmoothed Error on level %i', level));
            end
            resVec         = obj.solver.matrices{level}*u-f;
            resVecCoarse   = obj.restrict( resVec, level);
            correction     = zeros(length(resVecCoarse),1);
            correction     = obj.cycle( level+1, correction, resVecCoarse);
            correctionFine = obj.interpolate( correction, level+1);
            if DEBUGLEVEL>=10
               obj.solver.plotSolution(obj.solver.hierarchy.pointclouds{level},correctionFine,sprintf('Correction on level %i', level));
            end
            u              = u - correctionFine;
            if DEBUGLEVEL>=10
               err = u-solution;
               obj.solver.plotSolution(obj.solver.hierarchy.pointclouds{level},err,sprintf('Error after CGC on level %i', level));
            end
            u              = obj.smooth( obj.solver.matrices{level}, u, f, obj.nPostSmooth);
            if DEBUGLEVEL>=10
               err = u-solution;
               obj.solver.plotSolution(obj.solver.hierarchy.pointclouds{level},err,sprintf('Postsmoothed error on level %i', level));
            end
         end
      end

      function u = smooth (obj, A, u, f, iter)
      % SMOOTH  Smoother
      %     u = smooth(A,u,f,iter)  Apply iter iterations of smoother   
      %                             on Au=f. The smoother type is choosen by
      %                             obj.SMOOTHER.

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
      % RESTRICT  Restriction
      %     R = restrict(resVec,level) Restrict the residual vector resVec from
      %                                level level to level level+1

         NCoarse = obj.solver.hierarchy.pointclouds{level+1}.N;
         R       = zeros( NCoarse, 1);

         if ( obj.RESTRICTION == 1 )
            for i=1:NCoarse
               R(i) = resVec( obj.solver.hierarchy.coarse2fine{level+1}(i) );
            end
         end
      end

      function r = interpolate( obj, correction, level)
      % INTERPOLATE  Interpolation
      %     r = interpolate(correction,level) Interpolate correction from level
      %                                       level to level level-1
      %     TODO: Setup of interpolation should be done ONCE, not in every
      %     iteration!
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
