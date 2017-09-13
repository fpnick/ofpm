classdef Solver < handle
    %SOLVER represents a the overall solver
    %   Detailed explanation goes here

    properties
       pointcloud
       matrices
       rhss
       printlevel
       sol
       hierarchy
       useMultigrid=1
    end

    methods
       function obj = Solver(h,lbx,lby,ubx,uby,printlevel)

          obj.pointcloud = Pointcloud(h,lbx,lby,ubx,uby);
          obj.pointcloud.findNeighbours;
          obj.pointcloud.organize;
          obj.printlevel = printlevel;
          if ( printlevel > 9 )
             obj.pointcloud.plot;
          end

          obj.hierarchy = Hierarchy(obj.pointcloud);

       end

       function rho = advance(obj)

           
          for i=1:obj.hierarchy.depth
             obj.setupMatrix(obj.hierarchy.pointclouds{i},i);
             obj.setupRHS(obj.hierarchy.pointclouds{i},i);
          end

          % Solve system
          disp('Solving linear system...')
          tic
          % obj.matrices{1} = obj.matrices{1}./diag(obj.matrices{1});
          if ( ~ obj.useMultigrid )
             obj.sol = obj.matrices{1} \ obj.rhss{1};
             rho = -1;
          else
             mg = Multigrid(obj);
             [obj.sol,rho] = mg.solve(zeros(obj.pointcloud.N,1),10^(-8));
          end
          toc

          if ( obj.printlevel > 4 )
             obj.plotSolution(obj.pointcloud,obj.sol,"Solution")
          end
       end

       function plotSolution(obj,pointcloud,sol,descr)
          % disp('Plotting solution...')
          [X,Y] = meshgrid(pointcloud.lbx:sqrt(1/pointcloud.N):pointcloud.ubx,pointcloud.lby:sqrt(1/pointcloud.N):pointcloud.uby);
          Vq = griddata(pointcloud.coords(:,1),pointcloud.coords(:,2),sol,X,Y);
          figure
          mesh(X,Y,Vq);
          hold on
          plot3(pointcloud.coords(:,1),pointcloud.coords(:,2),sol,'o')
          title(descr)
          hold off
          % figure
          % contour(X,Y,Vq);
       end

    end

    methods (Access=private)
       function rhs = setupRHS(obj,pointcloud,level)
          disp(sprintf('Setting up rhs for level %d',level))
          rhs = zeros(pointcloud.N,1);
          for i=1:pointcloud.N
              if ( pointcloud.ibound_type(i)==0 )
                 % point is not on the boundary
                  rhs(i) = obj.loadFunction(pointcloud.coords(i,:));
              elseif ( pointcloud.ibound_type(i)==1 ) 
                  rhs(i) = obj.bcFunctionDirichlet(pointcloud.coords(i,:));
              elseif ( pointcloud.ibound_type(i)==2 )
                  rhs(i) = obj.bcFunctionNeumann(pointcloud.coords(i,:),pointcloud.ibound_location(i));
              end
          end
          obj.rhss{level} = rhs;
       end
          
       function f = loadFunction(obj,point)
          %LOADFUNCTION Summary of this function goes here
          %   Detailed explanation goes here
          x=point(1,1);
          y=point(1,2);

          % f = 0;
          % f = 1;
          f = -8*pi^2 * sin(2*pi*x) * sin(2*pi*y);
          %f = sin(2*pi*x) * sin(2*pi*y);
       end

       function [ f ] = bcFunctionDirichlet(obj,point)
       %BCFUNCTIONDIRICHLET Summary of this function goes here
       %   Detailed explanation goes here
          x=point(1,1);
          y=point(1,2);

          % f = 100;
          f = 0;
          %f = sin(x) * sin(y);
          % f = sin(2*pi*x) * sin(2*pi*y);
       end

       function [ f ] = bcFunctionNeumann(obj,point,location)
       %BCFUNCTIONNEUMANN Summary of this function goes here
       %   Detailed explanation goes here
          x=point(1,1);
          y=point(1,2);

          if ( location == 4 ) 
             f = -2*pi * sin( 2*pi * y);
          elseif ( location == 3 ) 
             f = 2*pi * sin( 2*pi * y);
          elseif (location == 1)
             f = -2*pi * sin( 2*pi * x);
          elseif (location == 2)
             f = 2*pi * sin( 2*pi * x);
          else
             disp('Error in bcFunctionNeumann')
          end

       end

       function setupMatrix(obj,pointcloud,level)
          disp(sprintf('Setting up matrix for level %d',level))
          tic
          %obj.matrix = sparse(pointcloud.N,pointcloud.N);
          b = 0;

          for i=1:pointcloud.N
             if ( pointcloud.ibound_type(i)==0 ) 
                b = [0;0;0;0;2;2]; 
                stencil{i} = obj.setupStencil(pointcloud,i,b);
                n = max(size(pointcloud.neighbourLists{i}));
                if ( n<20 )
                   % fprintf('Point %i has less than 20 neighbours\n', i);
                end
                ja{i} = pointcloud.neighbourLists{i}(1:n);
             elseif ( pointcloud.ibound_type(i)==2 )
                if ( pointcloud.ibound_location(i)==4 )
                   b = [-1;0;0;0;0;0]; 
                elseif ( pointcloud.ibound_location(i)==3 )
                   b = [1;0;0;0;0;0]; 
                elseif ( pointcloud.ibound_location(i)==2 )
                   b = [0;1;0;0;0;0]; 
                elseif ( pointcloud.ibound_location(i)==1 )
                   b = [0;-1;0;0;0;0]; 
                else
                   disp('Error in setupMatrix')
                end
                stencil{i} = obj.setupStencil(pointcloud,i,b);
                n = max(size(pointcloud.neighbourLists{i}));
                if ( n<20 )
                   % fprintf('Point %i has less than 20 neighbours\n', i);
                end
                ja{i} = pointcloud.neighbourLists{i}(1:n);
             end
          end
          nna=0;
          for i=1:pointcloud.N
             if ( pointcloud.ibound_type(i)==0 ) 
                nna = nna + length(ja{i});
             end
          end
          row = zeros(nna,1);
          col = zeros(nna,1);
          val = zeros(nna,1);
          ptr = 1;
          for i=1:pointcloud.N
             diag = 0.0;
             if ( pointcloud.ibound_type(i)~=1 )
                for j=1:length(stencil{i})
                   row(ptr) = i;
                   col(ptr) = ja{i}(j);
                   val(ptr) = -stencil{i}(j);
                   diag = diag + val(ptr);
                   ptr = ptr+1;
                end 
             else % Dirichlet BCs
                row(ptr) = i;
                col(ptr) = i;
                val(ptr) = 1.0;
                ptr = ptr+1;
             end
          end

          obj.matrices{level}=sparse(row,col,val);
          toc

       end

       function stencil = setupStencil(obj,pointcloud,i,b)
          n = max(size(pointcloud.neighbourLists{i}));

          K = obj.setupK( pointcloud.coords(pointcloud.neighbourLists{i}(1:n),:), pointcloud.coords(i,:));
          W = obj.setupWeightMatrix( pointcloud.distanceLists{i}(1:n), pointcloud.h);

          lambda = (K' * W^2 * K) \ (-b); 
          stencil = -W^2 * K * lambda;
       end

       function K = setupK(obj, coords, p0)
          coords_tmp = coords - p0;
          K = [ ones( size(coords_tmp(:,1))), coords_tmp(:,1), coords_tmp(:,2), coords_tmp(:,1).*coords_tmp(:,2) , coords_tmp(:,1).^2 , coords_tmp(:,2).^2 ];
       end

       function [ W ] = setupWeightMatrix(obj, distances, h)
          %SETUPWEIGHTMATRIX Summary of this function goes here
          %   Detailed explanation goes here
          GAMMA = 4;
          r = distances/h;
          W = diag(exp(-GAMMA * r.^2) - 0.0183);
          % Just to be sure...
          for i=1:length(distances)
             if r >= 1
                W(i,i) = 0;
             end
          end
       end
          

    end



end %class
