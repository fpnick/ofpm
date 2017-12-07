%> Finds strong (negative) couplings based on the definition in Stueben
%>
%> @param[in]   A   Connectivity matrix to find strong couplings
%> 
%> @param[out]  S   S(i,j) = 1 <=> i is strongly connected to j.


% PARAMETERS:
%   coupling_mode   1   -   Standard. Only consider negative couplings.
%                   2   -   Treat positive and negative couplings equally.
function S = strongCouplings(A)

    coupling_mode = 1;
    positive_couplings = 0;

    m = size(A,1);
    maxvals = zeros(1,m);
    
    estr = 0.25;

    [ii,jj,ss] = find(A);
    Srow = zeros(numel(ii),1);
    Scol = zeros(numel(ii),1);
    Sval = ones(numel(ii),1);

    for i = 1:numel(ii)
       row = ii(i);
       col = jj(i);
       val = ss(i);
       if ( row ~= col )
          if ( coupling_mode == 1 )
             if ( val < 0 && abs(val) > maxvals(row) )
                maxvals(row) = abs(val);
             end
          elseif ( coupling_mode == 2 )
             if ( abs(val) > maxvals(row) )
                maxvals(row) = abs(val);
             end
          end
       end
    end

    ptr=1;
    for i = 1:numel(ii)
       row = ii(i);
       col = jj(i);
       val = ss(i);
       if ( coupling_mode == 1 )
          if ( val < 0.0 && -val >= estr * maxvals(row))
             Srow(ptr) = row;
             Scol(ptr) = col;
             ptr = ptr+1;
          end
       elseif ( coupling_mode == 2 ) 
          if ( abs(val) >= estr * maxvals(row) )
             Srow(ptr) = row;
             Scol(ptr) = col;
             ptr = ptr+1;
             if ( val > 0 ) 
                positive_couplings = positive_couplings + 1;
             end
          end
       end
    end
    
    for i=1:m
        Srow(ptr)=i;
        Scol(ptr)=i;
        Sval(ptr)=0;
        ptr=ptr+1;
    end
    
    S = sparse(Srow(1:ptr-1),Scol(1:ptr-1),Sval(1:ptr-1));

    if ( positive_couplings > 0 )
        disp(sprintf('Used %i positive couplings\n', positive_couplings));
    end

end
