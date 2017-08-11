classdef Hierarchy < handle
    %HIERARCHY represents a hierarchy of pointclouds.
    %   Detailed explanation goes here

    properties
      pointclouds
      depth
      fine2coarse
      coarse2fine
    end

    methods
      function obj = Hierarchy(finePointcloud)
        obj.pointclouds{1} = finePointcloud;
        obj.pointclouds{1}.stats();
        i=1;
        while ( obj.pointclouds{i}.N > 200000 )
           [ obj.pointclouds{i+1}, obj.fine2coarse{i+1}, obj.coarse2fine{i+1} ] = obj.pointclouds{i}.coarsen;
           obj.pointclouds{i+1}.findNeighbours();
%            obj.pointclouds{i+1}.organize;
           obj.pointclouds{i+1}.stats();
           i=i+1;
        end
        obj.depth = length(obj.pointclouds);
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
