
addpath('..')

X=1:1:3
H=0.05

condition = zeros(length(X),1);
size = zeros(length(X),1);
eigmax = zeros(length(X),1);
eigmin = zeros(length(X),1);
diagdom = zeros(length(X),1);


parfor i=1:length(X)
    X(i)
    [A,rhs,sol]=ofpm_oo(H,0,0,X(i),X(i),0);
    condition(i) = condest(A);
    size(i) = length(A);
    eigmin(i) = eigs(A,1,'sm');
    eigmax(i) = eigs(A,1);
    diagdom(i) = measureDiagDom(A);
    
    A_norm = A./diag(A);
    condition_norm(i) = condest(A_norm);
    eigmin_norm(i) = eigs(A,1,'sm');
    eigmax_norm(i) = eigs(A,1);
    diagdom(i) = measureDiagDom(A_norm);
    
    matrizen{i} = A;
end

figure;
plot(X,condition,'s',X,condition_norm,'s');
title("Condition Numbers");
legend("Unscaled","Normalized");
legend('Location','northwest');
xlabel("Geometry Length");
figure;
semilogy(X,eigmin,'v',X,eigmax,'^');
title("Eigenvalues of Unscaled Systems (absolute values)");
legend("Min","Max");
legend('Location','southwest');
xlabel("Geometry Length");
figure;
semilogy(X,eigmin_norm,'v',X,eigmax_norm,'^');
title("Eigenvalues of Normalized Systems (absolute values)");
legend("Min","Max");
legend('Location','southwest');
xlabel("Geometry Length");
figure;
plot(X,diagdom,'x');
title("Diagonal Dominance");
xlabel("Geometry Length");

