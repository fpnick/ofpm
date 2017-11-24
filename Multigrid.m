classdef Multigrid < handle
    % MULTIGRID   represents a multigrid solver for the linear system.
    %    Represents the core multigrid solver. Assumes though, that the hierachy
    %    has already been created in the 'solver' object. 

    properties

      solver % Note that solver contains the complete hierarchy + matrices

      % Parameters
      SMOOTHER      = 2  % 1: Gauss-Seidel
                         % 2: Gauss-Seidel with RCM reordering
      RESTRICTION   = 3  % 1: Injection
                         % 2: "Half weighting"
                         % 3: "Half weighting" with row-scaling ResOp
      INTERPOLATION = 1  % 1: Weighted based on distance
                         % 2: like 1 but boundary couplings are taken into
                         % account when computing weights, but not for
                         % interpolation
      NORMALIZE     = 0  % 1: Normalize matrices on every level
      ENFORCE_DIAGDOM = 0 % 1: Add 5% to every diagonal
      nPreSmooth    = 1  % n: Number of pre-smoothing steps
      nPostSmooth   = 1 % n: Number of post-smoothing steps
      nMaxIter      = 100

      %
      restriction_setup_done
      resOp
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

      function obj = Multigrid(solver)
         obj.solver = solver;
         obj.restriction_setup_done = zeros( obj.solver.hierarchy.depth, 1 );
      end

      function [solution,rho] = solve(obj,u,tol)
      % SOLVE  Solve the linear system given by the hierachy in 'solver'
      %     solution = solve(u,tol) solves with a residual reduction of tol
      %                             starting with initial guess u.
         DEBUGLEVEL = 0;

         if ( obj.NORMALIZE == 1 ) 
            for i=1:obj.solver.hierarchy.depth
               obj.solver.matrices{i} = obj.solver.matrices{i}./diag(obj.solver.matrices{i});
               obj.solver.rhss{i} = obj.solver.rhss{i}./diag(obj.solver.matrices{i});
            end
         end
         if ( obj.ENFORCE_DIAGDOM == 1 ) 
            for i=1:obj.solver.hierarchy.depth
               obj.solver.matrices{i} = obj.solver.matrices{i}+0.05*diag(diag(obj.solver.matrices{i}));
            end
         end

         resVec     = obj.solver.matrices{1}*u-obj.solver.rhss{1};
         res        = norm( resVec, 2);
         if DEBUGLEVEL>=10
            condition = condest(obj.solver.matrices{1});
            fprintf('Condition is %1.3e\n', condition);
            % obj.solver.plotSolution(obj.solver.hierarchy.pointclouds{1},resVec,sprintf('Initial Residual'));
            solution_db = obj.solver.matrices{1} \ obj.solver.rhss{1};
         end

         res0       = res;
         tol_abs    = res * tol;
         iterations = 0;

         fprintf('Initial residual: %1.3e\n', res);
         
         while ( res > tol_abs && iterations < obj.nMaxIter )
            u          = obj.cycle( 1, u, obj.solver.rhss{1});
            resVec     = obj.solver.matrices{1}*u-obj.solver.rhss{1};
            res        = norm( resVec, 2);
            iterations = iterations + 1;

            if DEBUGLEVEL>=20
               % obj.solver.plotSolution(obj.solver.hierarchy.pointclouds{1},resVec,sprintf('Residual after iteration %i', iterations));
            end

            if DEBUGLEVEL>=10
               err = u-solution_db;
               % obj.solver.plotSolution(obj.solver.hierarchy.pointclouds{1},err,sprintf('Error after iteration %i', iterations));
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
         DEBUGLEVEL = 0;
         
         if DEBUGLEVEL>=0
            solution = obj.solver.matrices{level} \ f;
            err = u-solution;
            if ( level==1 )
               fprintf(' |Max| error on level %i: %1.3e\n', level,max(abs(err)))
               [~,index] = max(abs(err));
               fprintf('   Index: %i\n', index);
               xcoo = obj.solver.hierarchy.pointclouds{level}.coords(index,1);
               ycoo = obj.solver.hierarchy.pointclouds{level}.coords(index,2);
               fprintf('   Coords: %1.3e , %1.3e\n', xcoo,ycoo)
            end
         end
         if DEBUGLEVEL>=10
            % obj.solver.plotSolution(obj.solver.hierarchy.pointclouds{level},err,sprintf('Error on level %i', level));
         end
         if DEBUGLEVEL>=10
            resVec         = obj.solver.matrices{level}*u-f;
            obj.solver.plotSolution(obj.solver.hierarchy.pointclouds{level},resVec,sprintf('Residual on level %i', level));
         end
         if DEBUGLEVEL>0
            resVec         = obj.solver.matrices{level}*u-f;
            fprintf(' |Max| residual on level %i: %1.3e\n', level,max(abs(resVec)))
         end


         if ( level == obj.solver.hierarchy.depth )
            % Coarsest level => direct solve
            u = obj.solver.matrices{level} \ f;
            if DEBUGLEVEL>=2
               condition = condest(obj.solver.matrices{level});
               fprintf('Condition is %1.3e on CG \n', condition);
               % figure;
               % spy(obj.solver.matrices{level});
            end
            if DEBUGLEVEL>=10
               % obj.solver.plotSolution(obj.solver.hierarchy.pointclouds{level},u,sprintf('CG solution'));
               % obj.solver.hierarchy.pointclouds{level}.plot();
            end
         else
            u              = obj.smooth( obj.solver.matrices{level}, u, f, obj.nPreSmooth);
            if DEBUGLEVEL>=20
               err = u-solution;
               obj.solver.plotSolution(obj.solver.hierarchy.pointclouds{level},err,sprintf('Presmoothed Error on level %i', level));
            end
            resVec         = obj.solver.matrices{level}*u-f;
            if DEBUGLEVEL>=10
               % obj.solver.plotSolution(obj.solver.hierarchy.pointclouds{level},resVec,sprintf('Presmoothed Residual on level %i', level));
            end
            if DEBUGLEVEL>0
               fprintf(' |Max| presmoothed residual on level %i: %1.3e\n',level, max(abs(resVec)))
               fprintf(' Sum presmoothed residual on level %i: %1.3e\n',level,sum(resVec))
            end
            if DEBUGLEVEL>=2
               condition = condest(obj.solver.matrices{level});
               fprintf('Condition is %1.3e on level %i \n', condition, level);
            end
            resVecCoarse   = obj.restrict( resVec, level);
            if DEBUGLEVEL>=10
               % obj.solver.plotSolution(obj.solver.hierarchy.pointclouds{level+1},resVecCoarse,sprintf('Restricted Residual on level %i', level+1));
            end
            if DEBUGLEVEL>0
               fprintf(' |Max| restricted residual on level %i: %1.3e\n', level+1, max(abs(resVecCoarse)))
               fprintf(' Sum restricted residual on level %i: %1.3e\n',level+1, sum(resVecCoarse))
            end
            correction     = zeros(length(resVecCoarse),1);
            correction     = obj.cycle( level+1, correction, resVecCoarse);
            if DEBUGLEVEL>=10
               obj.solver.plotSolution(obj.solver.hierarchy.pointclouds{level+1},correction,sprintf('Correction before interpolated to level %i', level));
            end
            % correction     = obj.cycle( level+1, correction, obj.solver.rhss{level+1});
            correctionFine = obj.interpolate( correction, level+1);
            % correctionFine = 0.7 * correctionFine;
            if DEBUGLEVEL>=20
               err = u-solution;
               fprintf('Max error %1.3e\n',max(err));
               fprintf('Max corr_coarse %1.3e, Max corr_fine %1.3e \n',max(correction),max(correctionFine))
            end
            if DEBUGLEVEL>=20
               err = u-solution;
               obj.solver.plotSolution(obj.solver.hierarchy.pointclouds{level},correctionFine,sprintf('Correction on level %i', level));
               obj.solver.plotSolution(obj.solver.hierarchy.pointclouds{level},correctionFine-err,sprintf('corr-err on level %i', level));
            end
            u              = u - correctionFine;
            if DEBUGLEVEL>=20
               err = u-solution;
               obj.solver.plotSolution(obj.solver.hierarchy.pointclouds{level},err,sprintf('Error after CGC on level %i', level));
            end
            u              = obj.smooth( obj.solver.matrices{level}, u, f, obj.nPostSmooth);
            if DEBUGLEVEL>=20
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

         if ( obj.SMOOTHER == 2 ) 
            r = symrcm(A);
            Ar = A(r,r);
            ur = u(r);
            fr = f(r);
            Lr=tril(Ar,0);
            Rr=triu(Ar,1);
            u0 = ur;
            for i=1:iter
              % fprintf(' Smoothing iter %i\n',iter)
              ur = Lr\(fr - Rr*u0);
              u0=ur;
            end

            u = sort_back(ur,r',1);
         end

      end

      function R = restrict( obj, resVec, level)
      % RESTRICT  Restriction
      %     R = restrict(resVec,level) Restrict the residual vector resVec from
      %                                level level to level level+1

         NCoarse = obj.solver.hierarchy.pointclouds{level+1}.N;
         NFine   = obj.solver.hierarchy.pointclouds{level}.N;
         fine2coarse = obj.solver.hierarchy.fine2coarse{level+1};
         coarse2fine = obj.solver.hierarchy.coarse2fine{level+1};
         ibound_type_fine = obj.solver.hierarchy.pointclouds{level}.ibound_type;

         R       = zeros( NCoarse, 1);

         if ( obj.RESTRICTION == 1 )
            for i=1:NCoarse
               % fprintf('i %i, coarse2fine %i\n',i,obj.solver.hierarchy.coarse2fine{level+1}(i));
%                R(i) = resVec( obj.solver.hierarchy.coarse2fine{level+1}(i) ) * obj.solver.hierarchy.pointclouds{level}.HFACTOR_COARSENING;
               R(i) = resVec( coarse2fine(i) ) * 1.0;
            end

         elseif ( obj.RESTRICTION == 2 || obj.RESTRICTION == 3 || obj.RESTRICTION == 4)

            if ( obj.restriction_setup_done(level) == 0 )

               maxentries = 0;
               for i=1:NFine
                  maxentries = maxentries + length(obj.solver.hierarchy.pointclouds{level}.neighbourLists{i});
               end
               row = zeros(maxentries,1);
               col = zeros(maxentries,1);
               val = zeros(maxentries,1);
               ptr=1;

               for i=1:NFine
                  if ( fine2coarse(i) == 0 ) % F-Point!
                     neighbourList = obj.solver.hierarchy.pointclouds{level}.neighbourLists{i};
                     distanceList_hat  = 1./obj.solver.hierarchy.pointclouds{level}.distanceLists{i};
                     nNeighbours   = length(neighbourList);

                     sumDistances = 0;
                     for j=2:nNeighbours
                        % if ( fine2coarse(neighbourList(j)) ~= 0 && ibound_type_fine(neighbourList(j)) == 0)
                        if ( fine2coarse(neighbourList(j)) ~= 0)
                           sumDistances = sumDistances + distanceList_hat(j);
                        end
                     end
                     sumWeights = 0;
                     for j=2:nNeighbours
                        if ( fine2coarse(neighbourList(j)) ~= 0 && ibound_type_fine(neighbourList(j)) == 0)
                        % if ( fine2coarse(neighbourList(j)) ~= 0 )
                           if ( ibound_type_fine(neighbourList(j))==0 )
                              row(ptr) = fine2coarse(neighbourList(j));
                              col(ptr) = i;
                              val(ptr) = distanceList_hat(j)/sumDistances;
                              ptr = ptr+1;
                              sumWeights = sumWeights+(distanceList_hat(j)/sumDistances);
                           else 
                              % fprintf('C neighbour skipped')
                           end
                        end
                     end
                     if (sumWeights>1.1)
                        fprintf('Sum of weights %1.3e\n', sumWeights);
                     end
                  else % C-Point
                     row(ptr) = fine2coarse(i);
                     col(ptr) = i; 
                     val(ptr) = 1.0;
                     ptr=ptr+1;
                  end
               end

               row = row(1:ptr-1);
               col = col(1:ptr-1);
               val = val(1:ptr-1);
               obj.resOp{level} = sparse(row,col,val);

               if ( obj.RESTRICTION == 3 )
                  sums = sparse(sum(obj.resOp{level},2));
                  [rowidx, colidx, vals] = find(obj.resOp{level});
                  obj.resOp{level} = sparse(rowidx, colidx, vals./sums(rowidx),size(obj.resOp{level},1), size(obj.resOp{level},2));
                  % obj.resOp{level} = diag(1./sums) * obj.resOp{level};
               end

               obj.restriction_setup_done(level) = 1;
            end

            % R = 0.75*obj.resOp{level}*resVec + 0.25*resVec(fine2coarse~=0);
            R = obj.resOp{level}*resVec;
            % R = 0.50*obj.resOp{level}*resVec + 0.50*resVec(fine2coarse~=0);

            % for i=1:NCoarse
               % if ( obj.solver.hierarchy.pointclouds{level+1}.ibound_type(i) ~= 0 )
                  % R(i) = resVec( coarse2fine(i) );
               % end
            % end

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
         ibound_type_fine = obj.solver.hierarchy.pointclouds{level-1}.ibound_type;

         if ( obj.INTERPOLATION == 1 || obj.INTERPOLATION == 2)
            for i=1:NFine
               if ( fine2coarse(i) == 0 ) % F-Point!
                  neighbourList = obj.solver.hierarchy.pointclouds{level-1}.neighbourLists{i};
                  distanceList_hat  = 1./obj.solver.hierarchy.pointclouds{level-1}.distanceLists{i};
                  nNeighbours   = length(neighbourList);

                  sumDistances = 0;
                  for j=2:nNeighbours
                     if ( fine2coarse(neighbourList(j)) ~= 0 ) % C-Neighbour
                        sumDistances = sumDistances + distanceList_hat(j);
                     end
                  end
                  for j=2:nNeighbours
                     if ( fine2coarse(neighbourList(j)) ~= 0 ) % C-Neighbour
                        if ( ibound_type_fine(neighbourList(j)) == 0 || obj.INTERPOLATION == 1)
                           r(i) = r(i) + correction(fine2coarse(neighbourList(j))) * (distanceList_hat(j)/sumDistances);
                        end
                     end
                  end
               else
                  % fprintf('i %i, fine2coarse %i\n',i,fine2coarse(i));
                  r(i) = correction(fine2coarse(i));
               end
            end
         end
      end
      
    end

end %class
