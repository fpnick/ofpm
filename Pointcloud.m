classdef Pointcloud < handle
    %POINTCLOUD represents a pointcloud
    %   Detailed explanation goes here
    
    properties
        N
        h
        coords
        lbx
        lby
        ubx
        uby
        neighbourLists
        distanceLists
        ibound
    end
    
    methods
        function obj = Pointcloud(h,lbx,lby,ubx,uby,coords,ibound)
            
            obj.h = h;
            obj.lbx = lbx;
            obj.lby = lby;
            obj.ubx = ubx;
            obj.uby = uby;
            
            if ( nargin == 5 )
                obj.N = round(40/h^2);
                
                if ( 1 )
                   obj.coords = zeros(obj.N,2);
                    for i=1:obj.N
                        obj.coords(i,1) = unifrnd(obj.lbx,obj.ubx);
                    end
                    for i=1:obj.N
                        obj.coords(i,2) = unifrnd(obj.lby,obj.uby);
                    end
                else
                    Ntmp = floor(sqrt(abs(obj.N)));
                    for i=1:Ntmp-1
                        for j=1:Ntmp-1
                            obj.coords((i-1)*Ntmp+j-i+1,1) = 1/Ntmp * j;
                            obj.coords((i-1)*Ntmp+j-i+1,2) = 1/Ntmp * i;
                        end
                    end
                    obj.N = length(obj.coords);
                end
                
                NInterior = obj.N;
                
                hBnd = 1/sqrt(obj.N);
                x_bottom = obj.lbx:hBnd:obj.ubx;
                y_bottom = ones(1,size(x_bottom,2))*obj.lby;
                x_top = x_bottom;
                y_top = ones(1,size(x_bottom,2))*obj.uby;
                y_right = obj.lby+hBnd:hBnd:obj.uby-hBnd;
                x_right = ones(1,size(y_right,2))*obj.ubx;
                y_left = y_right;
                x_left = ones(1,size(y_right,2))*obj.lbx;
                
                obj.N = obj.N + length(x_bottom) + length(x_top) + length(x_right) + length(x_left);
                bc = [ x_bottom' y_bottom'; x_top' y_top'; x_right' y_right'; x_left' y_left' ];
                obj.coords = [ obj.coords; bc ];
                
                obj.ibound = zeros(obj.N,1);
                obj.ibound(NInterior+1:obj.N) = 1;
            elseif ( nargin == 7 )
                 obj.N = length(coords);
                 obj.coords = coords;
                 obj.ibound = ibound;
            end
        end

        function findNeighbours(obj)
            [obj.neighbourLists,obj.distanceLists] = rangesearch(obj.coords,obj.coords,obj.h);
        end

        function organize(obj)
            is_active = ones(obj.N,1);
            for i = 1:obj.N
                if is_active(i) && obj.ibound(i)==0
                    for j = 1:length(obj.neighbourLists{i})
                        if i ~= obj.neighbourLists{i}(j) && is_active(obj.neighbourLists{i}(j))
                            if  obj.distanceLists{i}(j) < obj.h*0.1 
                                is_active(i) = 0;
                                break;
                            end
                        end
                    end
                end
            end
            obj.N = sum(is_active);
            obj.coords = [obj.coords(find(is_active),1),obj.coords(find(is_active),2)];
            obj.findNeighbours;
            for i=1:obj.N
               if ( obj.isBoundary(obj.coords(i,:)) )
                   obj.ibound(i) = 1;
               end
           end
        end

        function bol = isBoundary(obj,point)
           % ISBOUNDARY return 1 iff the given point is on a boundary
           bol = (point(1) == obj.lbx || point(2) == obj.lby ||point(1) == obj.ubx || point(2) == obj.uby);
        end

        function plot(obj)
            figure;
            plot(obj.coords(:,1),obj.coords(:,2),'.');
        end

        function drawStar(obj,index)
            obj.plot();
            hold on;
            A=[obj.coords(index,1)*ones(1,length(obj.neighbourLists{index}));obj.coords(obj.neighbourLists{index},1)'];
            B=[obj.coords(index,2)*ones(1,length(obj.neighbourLists{index}));obj.coords(obj.neighbourLists{index},2)'];
            line(A,B);
            hold off;
        end

        function stats(obj)
            avgNeighbours = 0;
            for i=1:obj.N
               avgNeighbours = avgNeighbours + length(obj.neighbourLists{i});
            end
            avgNeighbours = avgNeighbours / obj.N;
            fprintf('Pointcloud size: %d\n', obj.N);
            fprintf('Average neighbours: %f\n', avgNeighbours);
        end
        
        function [ coarsePointcloud, fine2coarse, coarse2fine ] = coarsen(obj)
            H_FACTOR=0.3; % This is the factor that deremines the coarsening rate
            level = zeros(obj.N,1);
            nC = 0;
            nF = 0;            
            fine2coarse = zeros(obj.N);
            coarse2fine = zeros(nC);
            for i=1:obj.N
                if ( level(i) == 0 )
                    level(i) = 2;
                    nC = nC + 1;
                    coarse2fine(nC) = i;
                    fine2coarse(i) = nC;
                    for j=2:length(obj.neighbourLists{i})
                        if ( level(obj.neighbourLists{i}(j)) == 0 && obj.distanceLists{i}(j) <= obj.h*H_FACTOR )
                            level(obj.neighbourLists{i}(j)) = 1;
                            nF = nF +1;
                        end
                    end
                end
            end             
            
            % h->H isn't really correct here
            coarsePointcloud = Pointcloud(sqrt((nC+nF)/nC)*obj.h,obj.lbx,obj.lby,obj.ubx,obj.uby,obj.coords(find(level==2),:),obj.ibound(find(level==2)));
            
        end
    end % METHODS
    
end % CLASS

