
x=0.03:0.01:0.1; 
condition = zeros(length(x),1);
size = zeros(length(x),1);
eigmax = zeros(length(x),1);
eigmin = zeros(length(x),1);


parfor i=1:length(x)
    x(i)
    [A,rhs,sol]=ofpm_oo(x(i),0);
    condition(i) = condest(A);
    size(i) = length(A);
    eigmin(i) = eigs(A,1,'sm');
    eigmax(i) = eigs(A,1);
    
    A_norm = A./diag(A);
    condition_norm(i) = condest(A_norm);
    eigmin_norm(i) = eigs(A,1,'sm');
    eigmax_norm(i) = eigs(A,1);
    
    matrizen{i} = A;
end
