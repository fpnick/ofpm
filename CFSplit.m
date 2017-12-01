%>  Based on the Ruge-Stueben idea, this function computes 
%>  a splitting of the variables in fine- and coarse-grid variables.
%>
%>  @param[in]  S   Matrix indicating strong connections. Strong connections
%>                  should be indicated by S(i,j) = 1.
%>  
%>  @param[out] splitting   Vector of length nnu, whereas splitting(i)=1 indicates
%>                          a coarse grid variable and =-1 indicates a fine grid
%>                          variable.
%>  @param[out] num_cg_vars Number of variables on coarse grid.
%>
%>  @todo Write routine for updating lambda locally instead of re-initializing it.
function [splitting, num_cg_vars] = CFSplit(S,A)

    n = size(S,1);
    splitting = zeros(n,1);   % splitting(i) =  0 : Undecided
                              %                 1 : Coarse grid
                              %                -1 : Fine grid

    num_U = n;
    num_cg_vars=0;

    lambda = initLambda(S, splitting);

    mx = 1;
    while ( mx > 0 )
        % Decide on next coarse grid variable.
        [mx, next_coarse] = max(lambda);
        if ( mx > 0 )
            splitting(next_coarse) = 1;
            num_U = num_U - 1;
            num_cg_vars = num_cg_vars+1;

            % Put strongly connected variables in F.
            [row,col,val] = find(S(:,next_coarse));
            for j=1:numel(row)
                if ( S(row(j),next_coarse) == 1 && splitting(row(j)) == 0 )
                    splitting(row(j)) = -1;
                    num_U = num_U - 1;
                end
            end

            % Update lambda
            lambda = updateLambda(lambda, next_coarse, A);
        end
    end

    left_over_f = 0;
    for i=1:n
        if ( splitting(i) == 0 )
            splitting(i) = -1;
            left_over_f = left_over_f + 1;
        end
    end
    if ( left_over_f > 0 )
        disp(sprintf('%i variables were left over in the coarsening process.\n',left_over_f));
    end

end

function lambda = updateLambda(lambda, iChanged, A)
    
    n = size(A,1);
    lambda(iChanged) = 0;
    
    [row,col,val] = find(A(iChanged,:));
    for j = 1:numel(col)
        lambda(col(j)) = lambda(col(j)) - 1;
    end
end
    

%>  Initialize lambda, which is being used to decide which variable becomes a C variable next.
%>
%>  @param[in]  S           Matrix indicating strong connections. Strong connections
%>                          should be indicated by S(i,j) = 1
%>  @param[in]  splitting   Vector of length nnu, whereas splitting(i)=1 indicates
%>                          a coarse grid variable and =-1 indicates a fine grid
%>                          variable. splitting(i)=0 indicates the variable is
%>                          still undecided.
%>  @param[out] lambda      lambda as being used in the Ruge-Stueben coarsening process. 
function lambda = initLambda(S, splitting)

    n = size(S,1);
    lambda = zeros(n,1);

    for i=1:n
        if ( splitting(i) == 0 )
            [~,col,~] = find(S(i,:));
            sum1 = 0;
            sum2 = 0;
            for j=1:numel(col)
                if ( S(col(j),i) == 1 )
                   if ( splitting(col(j)) == 0 )
                      sum1 = sum1 + 1;
                   elseif ( splitting(col(j)) == -1 )
                      sum2 = sum2 + 1;
                   end
                end
            end
            lambda(i) = sum1 + 2*sum2;
        else
            lambda(i) = 0;
        end
    end
end

