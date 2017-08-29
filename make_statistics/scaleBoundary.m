function [ B ] = scaleBoundary( A, ibound )
%SCALEBONDART Summary of this function goes here
%   Detailed explanation goes here

    innerAvg=0;
    bndAvg=0;
    inner=0;
    bndDirichlet=0;

    for i=1:max(size(A))
        if ( ibound(i) == 0 )
            innerAvg=innerAvg+A(i,i);
            inner = inner + 1;
        elseif ( ibound(i) ~= 1)
            bndAvg=bndAvg+A(i,i);
            bndDirichlet = bndDirichlet + 1;
        end
    end
    
    innerAvg = innerAvg / inner;
    bndAvg = bndAvg / bndDirichlet;
    
    scalFactor = innerAvg/bndAvg;
    
    tmp = ibound*scalFactor;
    for i=1:length(tmp)
        if ( tmp(i) == 0 )
            tmp(i) = 1;
        end
    end
    
    B = bsxfun(@times,A,tmp);
%     for i=1:max(size(B))
%         if ( ibound(i)>0 )
%             B(i,:) = B(i,:) * scalFactor;
%         end
%     end
            

end

