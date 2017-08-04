classdef Hierarchy < handle
    %HIERARCHY represents a the overall solver
    %   Detailed explanation goes here

    properties
      pointclouds
      depth
    end

    methods
      function obj = Hierarchy(finePointcloud)
        obj.pointclouds{1} = finePointcloud;
        i=1;
        while ( obj.pointclouds{i}.N > 50 )
           obj.pointclouds{i+1} = obj.pointclouds{i}.coarsen;
           obj.pointclouds{i+1}.findNeighbours();
           obj.pointclouds{i+1}.organize;
           i=i+1;
           if ( i>5 ) 
              break;
           end
        end
        obj.depth = length(obj.pointclouds);
      end
    end

    methods (Access=private)
    end



end %class
