
addpath('..')

X=1:1:11
H=0.05

condition = zeros(length(X),1);
size = zeros(length(X),1);
eigmax = zeros(length(X),1);
eigmin = zeros(length(X),1);
diagdom = zeros(length(X),1);


parfor i=1:length(X)
    X(i)
    [A,rhs,sol,pointcloud]=ofpm_oo(H,0,0,X(i),X(i),0);
    condition(i) = condest(A);
    size(i) = length(A);
    eigmin(i) = eigs(A,1,'sm');
    eigmax(i) = eigs(A,1);
    diagdom(i) = measureDiagDom(A);
    
    A_norm = A./diag(A);
    condition_norm(i) = condest(A_norm);
    eigmin_norm(i) = eigs(A_norm,1,'sm');
    eigmax_norm(i) = eigs(A_norm,1);
    diagdom(i) = measureDiagDom(A_norm);
    
    A_scaledBoundary = scaleBoundary(A,pointcloud.ibound);
    condition_scaledBoundary(i) = condest(A_scaledBoundary);
    eigmin_scaledBoundary(i) = eigs(A_scaledBoundary,1,'sm');
    eigmax_scaledBoundary(i) = eigs(A_scaledBoundary,1);
    diagdom_scaledBoundary(i) = measureDiagDom(A_scaledBoundary);
    
    matrizen{i} = A;
end

figure;
plot(X,condition,'s',X,condition_norm,'s',X,condition_scaledBoundary,'s');
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
% figure;
% plot(X,diagdom,'x');
% title("Diagonal Dominance");
% xlabel("Geometry Length");
figure;
semilogy(size,eigmin_scaledBoundary,'v',size,eigmax_scaledBoundary,'^');
title("Eigenvalues of System with Scaled Boundaries (absolute values)");
legend("Min","Max");
legend('Location','southwest');
xlabel("Matrix Size");