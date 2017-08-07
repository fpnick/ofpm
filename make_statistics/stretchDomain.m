
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
    [A,rhs,sol,pointcloud]=ofpm_oo(H,0,0,X(i),1,0);
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
plot(size,condition,'o',size,condition_norm,'s',size,condition_scaledBoundary,'x');
title(sprintf('Condition Numbers; H=%3.2f',H));
legend('Unscaled','Normalized','Scaled Boundary');
legend('Location','northwest');
xlabel('Matrix Size');
figure;
semilogy(size,eigmin,'v',size,eigmax,'^');
title(sprintf('Eigenvalues of Unscaled Systems (absolute values); H=%3.2f',H));
legend("Min","Max");
legend('Location','southwest');
xlabel("Matrix Size");
figure;
semilogy(size,eigmin_norm,'v',size,eigmax_norm,'^');
title(sprintf('Eigenvalues of Normalized Systems (absolute values); H=%3.2f',H));
legend('Min','Max');
legend('Location','southwest');
xlabel('Matrix Size');
% figure;
% plot(X,diagdom,'x');
% title("Diagonal Dominance");
% xlabel("Geometry Length");
figure;
semilogy(size,eigmin_scaledBoundary,'v',size,eigmax_scaledBoundary,'^');
title(sprintf('Eigenvalues of System with Scaled Boundaries (absolute values); H=%3.2f',H));
legend('Min','Max');
legend('Location','southwest');
xlabel('Matrix Size');
