classdef Hierarchy < handle
    %HIERARCHY represents a hierarchy of pointclouds.
    %   Detailed explanation goes here

    properties
      pointclouds
      depth
      fine2coarse
      coarse2fine
      MAXLEVELS = 2
    end

    methods
      function obj = Hierarchy(finePointcloud)
        obj.pointclouds{1} = finePointcloud;
        obj.pointclouds{1}.stats();
        i=1;
        obj.depth=1;
        while ( obj.pointclouds{i}.N > 1200 && i < obj.MAXLEVELS )
           [ obj.pointclouds{i+1}, obj.fine2coarse{i+1}, obj.coarse2fine{i+1} ] = obj.pointclouds{i}.coarsen;
           if ( obj.pointclouds{i+1}.N < 1200 || obj.pointclouds{i}.N / obj.pointclouds{i+1}.N < 1.5)
              obj.depth = i;
              break;
           else
              obj.pointclouds{i+1}.findNeighbours();
              obj.pointclouds{i+1}.organize;
              obj.pointclouds{i+1}.stats()
              i=i+1;
              obj.depth = i;
           end
        end
      end

      function plot(obj,depth)
        figure;
        hold on;

        for i=1:depth
           plot(obj.pointclouds{i}.coords(:,1),obj.pointclouds{i}.coords(:,2),'.')
        end
        legend('Location','northwest');

        hold off;
      end

    end

    methods (Access=private)
    end



end %class
