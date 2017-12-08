for i=66:-1:1
   load(sprintf('/home/fabian/privatecloud/figs/matrix_notnormalized_hrefine_%i.mat',i));
   % A./diag(A);
   ibound = findBoundary(A);
   B = scaleBoundary(A,ibound);
   A=B;

   condition_refine(i) = condest(A);
   eigmin_refine(i) = eigs(A,1,'sm');
   eigmax_refine(i) = eigs(A,1);
   save('/home/fabian/papers/paper-darmstadt/paper/figures/condition_refine_bndscaled','condition_refine');
   save('/home/fabian/papers/paper-darmstadt/paper/figures/eigmin_refine_bndscaled','eigmin_refine');
   save('/home/fabian/papers/paper-darmstadt/paper/figures/eigmax_refine_bndscaled','eigmax_refine');
end
% for i=1:9
   % load(sprintf('/media/data2/matrices_stretch/stretchmatrix%i.mat',i));
   % A./diag(A);
   % condition_stretch(i) = condest(A);
   % eigmin_stretch(i) = eigs(A,1,'sm');
   % eigmax_stretch(i) = eigs(A,1);
   % save('/home/fabian/papers/paper-darmstadt/paper/figures/condition_stretch_norm','condition_stretch');
   % save('/home/fabian/papers/paper-darmstadt/paper/figures/eigmin_stretch_norm','eigmin_stretch');
   % save('/home/fabian/papers/paper-darmstadt/paper/figures/eigmax_stetch_norm','eigmax_stretch');
% end
% for i=1:20
   % load(sprintf('/media/data2/matrices_stretch1/stretchmatrix_1direction_%i.mat',i));
   % A./diag(A);
   % condition_stretch1(i) = condest(A);
   % eigmin_stretch1(i) = eigs(A,1,'sm');
   % eigmax_stretch1(i) = eigs(A,1);
   % save('/home/fabian/papers/paper-darmstadt/paper/figures/condition_stretch1_norm','condition_stretch1');
   % save('/home/fabian/papers/paper-darmstadt/paper/figures/eigmin_stretch1_norm','eigmin_stretch1');
   % save('/home/fabian/papers/paper-darmstadt/paper/figures/eigmax_stretch1_norm','eigmax_stretch1');
% end
