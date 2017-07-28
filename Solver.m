classdef Solver < handle
    %SOLVER represents a the overall solver
    %   Detailed explanation goes here

    properties
       pointcloud
       matrix
       rhs

    end

    methods
       function obj = Solver(h,printlevel)

          obj.pointcloud = Pointcloud(h,0,0,1,1);
          obj.pointcloud.findNeighbours;
          obj.pointcloud.organize;
          if ( printlevel > 9 )
             obj.pointcloud.plot;
          end

       end

       function advance(obj)
          obj.setupMatrix(obj.pointcloud);
          obj.setupRHS(obj.pointcloud);
       end

    end

    methods (Access=private)
       function rhs = setupRHS(obj,pointcloud)
          disp('Setting up RHS...')
          rhs = zeros(pointcloud.N,1);
          for i=1:pointcloud.N
              if ( ~pointcloud.ibound(i) )
                  rhs(i) = obj.loadFunction(pointcloud.coords(i,:));
              else
                  rhs(i) = obj.bcFunction(pointcloud.coords(i,:));
              end
          end
       end
          
       function f = loadFunction(obj,point)
          %LOADFUNCTION Summary of this function goes here
          %   Detailed explanation goes here
          x=point(1,1);
          y=point(1,2);

          %f = 0;
          %f = -2*sin(x) * sin(y);
          f = -8*pi^2 * sin(2*pi*x) * sin(2*pi*y);
          %f = sin(2*pi*x) * sin(2*pi*y);
       end

       function [ u ] = bcFunction(obj,point)
       %BCFUNCTION Summary of this function goes here
       %   Detailed explanation goes here
          x=point(1,1);
          y=point(1,2);

          %u = 1;
          %u = sin(x) * sin(y);
          u = sin(2*pi*x) * sin(2*pi*y);
       end

       function setupMatrix(obj,pointcloud)
          disp('Setting up matrix...')
          obj.matrix = sparse(pointcloud.N,pointcloud.N);

          for i=1:pointcloud.N
             if ( pointcloud.ibound(i)==0 ) 
                stencil{i} = obj.setupStencil(pointcloud,i);
             end
          end

          for i=1:pointcloud.N
             if ( pointcloud.ibound(i)==0 )
                n = max(size(pointcloud.neighbourLists{i}));
                obj.matrix(i,pointcloud.neighbourLists{i}(2:n)) = -stencil{i};
                obj.matrix(i,i) = sum(stencil{i});
             else
                obj.matrix(i,i) = 1.0;
             end 
          end

       end

       function stencil = setupStencil(obj,pointcloud,i)
          n = max(size(pointcloud.neighbourLists{i}));

          K = obj.setupK(pointcloud.coords(pointcloud.neighbourLists{i}(2:n),:),pointcloud.coords(i,:));
          W = obj.setupWeightMatrix(pointcloud.distanceLists{i}(2:n));

          b = [0;0;0;0;2;2]; % This determines, what operator is approximated!

          lambda = -b \ (K' * W^2 * K);
          stencil = -W^2 * K * lambda';
       end

       function K = setupK(obj, coords, p0 )
          coords_tmp = coords - p0;
          K = [ ones(size(coords_tmp(:,1))), coords_tmp(:,1), coords_tmp(:,2), coords_tmp(:,1).*coords_tmp(:,2) , coords_tmp(:,1).^2 , coords_tmp(:,2).^2 ];
       end

       function [ W ] = setupWeightMatrix(obj,distances )
          %SETUPWEIGHTMATRIX Summary of this function goes here
          %   Detailed explanation goes here
          GAMMA = 4;
          W = diag(exp(-GAMMA * distances.^2) - 0.0183);
       end
          

    end



end %class
