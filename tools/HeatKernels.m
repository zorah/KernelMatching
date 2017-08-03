function [ Q,beta ] = HeatKernels( surface ,sample, n_evecs)
%% Calculate surface laplacian

% compute cotan Laplacian
nv = length(surface.VERT);
% G1 = metric_scale(surface.VERT(:,1), surface.VERT(:,2), surface.VERT(:,2), surface.TRIV,0);
% [W1, A1]=laplace_beltrami(surface.TRIV,G1);
[W1, A1]=geomProcessing.laplacian(surface.VERT,surface.TRIV);
% A1 = spdiags(A1,0,nv,nv);

[surface.evecs,surface.evals]=eigs(W1,A1,n_evecs,'sm'); %solves W*Phi = A*Phi*Lambda
surface.evals = -diag(surface.evals);

%[surface.evecs,surface.evals,~,~] = geomProcessing.cotan_LB(surface.VERT,surface.TRIV,n_evecs);

[~,sorted_idx] = sort(surface.evals,'ascend');
surface.evals = surface.evals(sorted_idx);
surface.evecs = surface.evecs(:,sorted_idx);

Q = surface.evecs(sample,:);
beta = surface.evals;

% H = 1000 * Q * beta * Q';


end

