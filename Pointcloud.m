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
                obj.coords = zeros(obj.N,2);
                
                if ( obj.N > 0 )
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
                end
                
                NInterior = obj.N;
                
                x_bottom = obj.lbx:h:obj.ubx;
                y_bottom = ones(1,size(x_bottom,2))*obj.lby;
                x_top = x_bottom;
                y_top = ones(1,size(x_bottom,2))*obj.uby;
                y_right = obj.lby+obj.h:obj.h:obj.uby-obj.h;
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
                if is_active(i)
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
            ptr=0;
            index2activeIndex=zeros(obj.N,1);
            for i = 1:obj.N
                if ( is_active(i) )
                    ptr=ptr+1;
                end
                index2activeIndex(i) = ptr;
            end
            obj.N = sum(is_active);
            obj.coords = [obj.coords(find(is_active),1),obj.coords(find(is_active),2)];
            obj.findNeighbours;
        end
        function plot(obj)
            plot(obj.coords(:,1),obj.coords(:,2),'.');
        end
        
        function coarsePointcloud = coarsen(obj)
            H_FACTOR=2;
            level = zeros(obj.N,1);
            for i=1:obj.N
                if ( level(i) == 0 )
                    level(i) = 2;
                    for j=2:length(obj.neighbourLists{i})
                        if ( level(obj.neighbourLists{i}(j)) == 0 && obj.distanceLists{i}(j) <= obj.h/H_FACTOR )
                            level(obj.neighbourLists{i}(j)) = 1;
                        end
                    end
                end
            end
            
            coarsePointcloud = Pointcloud(obj.h,obj.lbx,obj.lby,obj.ubx,obj.uby,obj.coords(find(level==2),:),obj.ibound(find(level==2)));
            
        end
    end % METHODS
    
end % CLASS

