
addpath('..')

x=0.005:0.001:0.07; 
% For Darmstadt paper:
% x=0.005:0.001:0.07; 
 % x=0.3:0.1:0.5; 
condition = zeros(length(x),1);
size = zeros(length(x),1);
eigmax = zeros(length(x),1);
eigmin = zeros(length(x),1);
diagdom = zeros(length(x),1);
solmax = zeros(length(x),1);
rho = zeros(length(x),1);
iter_needed = zeros(length(x),1);


for i=1:length(x)
    x(i)
    [A,rhs,sol,pointcloud,rho(i),iter_needed(i)]=ofpm_oo(x(i),0,0,1,1,0);
    solmax(i) = max(sol);

    save(sprintf('/home/fabian/papers/paper-darmstadt/paper/figures/matrix%i',i),'A');
    
%     condition(i) = condest(A);
    size(i) = length(A);
%     eigmin(i) = eigs(A,1,'sm');
%     eigmax(i) = eigs(A,1);
%     diagdom(i) = measureDiagDom(A);
%     
%     A_norm = A./diag(A);
%     condition_norm(i) = condest(A_norm);
%     eigmin_norm(i) = eigs(A_norm,1,'sm');
%     eigmax_norm(i) = eigs(A_norm,1);
%     diagdom(i) = measureDiagDom(A_norm);
    
    % A_scaledBoundary = scaleBoundary(A,pointcloud.ibound_location);
    % condition_scaledBoundary(i) = condest(A_scaledBoundary);
    % eigmin_scaledBoundary(i) = eigs(A_scaledBoundary,1,'sm');
    % eigmax_scaledBoundary(i) = eigs(A_scaledBoundary,1);
    % diagdom_scaledBoundary(i) = measureDiagDom(A_scaledBoundary);
    
%     matrizen{i} = A;
end

% figure;
% % plot(size,condition,'o',size,condition_norm,'s',size,condition_scaledBoundary,'x');
% plot(size,condition,'o',size,condition_norm,'s');
% title("Condition Numbers");
% % legend("Unscaled","Normalized","Scaled Boundary");
% legend("Original","Normalized");
% legend('Location','northwest');
% xlabel("Matrix Size");
% figure;
% semilogy(size,eigmin,'v',size,eigmax,'^');
% title("Eigenvalues of Original Systems (absolute values)");
% legend("Min","Max");
% legend('Location','southwest');
% xlabel("Matrix Size");
% figure;
% semilogy(size,eigmin_norm,'v',size,eigmax_norm,'^');
% title("Eigenvalues of Normalized Systems (absolute values)");
% legend("Min","Max");
% legend('Location','southwest');
% xlabel("Matrix Size");
% figure;
% plot(size,diagdom,'x');
% title("Diagonal Dominance");
% xlabel("Matrix Size");
% figure;
% semilogy(size,eigmin_scaledBoundary,'v',size,eigmax_scaledBoundary,'^');
% title("Eigenvalues of System with Scaled Boundaries (absolute values)");
% legend("Min","Max");
% legend('Location','southwest');
% xlabel("Matrix Size");
% figure;
% semilogx(size,rho,'ro');
% title("Convergence rates");
% xlabel("Matrix Size");
figure;
loglog(size,iter_needed,'ro');
title("Iterations to Converge");
xlabel("Matrix Size");
