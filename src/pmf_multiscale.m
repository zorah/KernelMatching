function [ matches_out ] = pmf_multiscale( X, Y, seed_X, seed_Y, ...
    ind_forbidden_X, ind_forbidden_Y, params)
%MULTISCALE Solving the product manifold filter in a multiscale approach. 
% Inputs: 
%   X, Y - mesh struct with .desc
%   seed_X, seed_Y - indices of initial matches on X, Y, ordered in
%                    corresponding way
%   ind_forbidden_X/Y - indices of forbidden vertices
%   params - see example file
%
% Efficient Deformable Shape Correspondence via Kernel Matching
% Copyright (C) 2017  Zorah Lähner (laehner@in.tum.de)
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% version 07.2017

warning('off')

if ~exist('params', 'var')
   params = struct;
end
if ~isfield(params, 'problemSize')
    params.problemSize = 500;
end
if ~isfield(params, 'k')
    params.k = 2;
end
if ~isfield(params, 'max_num_cells')
    params.max_num_cells = ceil(0.1 * X.n);
end
if ~isfield(params, 'filter')
    params.filter = 'hk';
end
if ~isfield(params, 'mu')
    params.mu = 1;
end
if ~isfield(params, 'initSeeds')
    params.initSeeds = 1000;
end

if params.initSeeds > length(seed_X)
    params.initSeeds = length(seed_X);
end

max_num_iter = ceil(log(X.n / params.problemSize) / log(params.k)) + 10;

%% initialize and add samples by fps

% 
init_seeds = fps_euclidean(X.VERT(seed_X), params.initSeeds, 1);

seed_X_orig = seed_X(init_seeds);
seed_Y_orig = seed_Y(init_seeds);

% indicators for seeds
is_seed_X = false(X.n,1);
is_seed_X(seed_X) = true;
is_seed_Y = false(Y.n,1);
is_seed_Y(seed_Y) = true;

% indicators for forbidden points
is_forbidden_X = false(X.n,1);
is_forbidden_X(ind_forbidden_X) = true;
is_forbidden_Y = false(Y.n,1);
is_forbidden_Y(ind_forbidden_Y) = true;

% fps sample new points from high res
sample_X = [];
sample_Y = [];

% indicator for samples, seeds and samples are disjunct sets
is_sampled_X = false(X.n,1);
is_sampled_X(sample_X) = true;
is_sampled_Y = false(Y.n,1);
is_sampled_Y(sample_Y) = true;

% init matches 
matches_out = zeros(X.n,1);
matches_out(sample_X) = sample_Y;
matches_out(seed_X) = seed_Y; 

%%
step = 1;

% dividing into voronoi cells, init everything to one cell
cell_index_X = zeros(X.n,1);
cell_index_Y = zeros(Y.n,1);
cell_index_X(sample_X) = 1;
cell_index_X(seed_X) = 1;
cell_index_Y(sample_Y) = 1;
cell_index_Y(seed_Y) = 1;

num_cells = 1;

% dividing into voronoi cells, the work
while step < max_num_iter && any(~(is_forbidden_X | is_seed_X)) ...
        && any(~(is_forbidden_Y | is_seed_Y))
    fprintf('starting step: %d \n',step)
    % solve problem at the current scale
    
    is_updated_matches = false;
    for i_cell = 1:num_cells  
        % find sample points in this cell
        sample_X_i = find( (cell_index_X == i_cell) & is_sampled_X);
        sample_Y_i = find( (cell_index_Y == i_cell) & is_sampled_Y);  
        
        % all possible vertices in this cell have been matched already
        if isempty(sample_X_i) || isempty(sample_Y_i)
            continue;
        end
        
        % check for convergence
        is_updated_matches = true;
        
        % get all seeds in this cell
        seed_X_i = find( (cell_index_X == i_cell) & is_seed_X);
        seed_Y_i = matches_out(seed_X_i);
        
        % prestore evecs or distances
        if strcmp(params.filter, 'hk')
           evecsD_X = X.evecs([seed_X_i; sample_X_i; seed_X_orig], :);
           evecsD_Y = Y.evecs([seed_Y_i; sample_Y_i; seed_Y_orig], :);
        elseif strcmp(params.filter, 'gauss')
           evecsD_X = X.D([seed_X_i; sample_X_i], [seed_X_i; sample_X_i]);
           evecsD_Y = Y.D([seed_Y_i; sample_Y_i], [seed_Y_i; sample_Y_i]);
        end
        
        % solve matching in this cell on {seeds, samples}
        tic,
        [matches, non_matched_X_i, non_matched_Y_i] = ...
            mfilter(X.desc([seed_X_i; sample_X_i], :), Y.desc([seed_Y_i; sample_Y_i], :), ...
            evecsD_X, evecsD_Y,...
            X.evals, Y.evals, seed_X_i, seed_Y_i, sample_X_i, sample_Y_i, params, step, params.mu);
        time_elapsed = toc;
        fprintf('step: %d, i_cell: %d of %d, num_seeds: %d, num_samples_X: %d, num_samples_Y: %d, time: %.2f \n',...
                      step, i_cell, num_cells, numel(seed_X_i), numel(sample_X_i) , numel(sample_Y_i), time_elapsed);
 
        matches_out([seed_X_i; sample_X_i]) = matches; % mind the ordering
        
        % make sample points the new seed points
        is_sampled_X(sample_X_i) = false;
        is_sampled_Y(sample_Y_i) = false;
        sample_X_i = setdiff(sample_X_i, non_matched_X_i);
        sample_Y_i = setdiff(sample_Y_i, non_matched_Y_i);
        is_seed_X(sample_X_i) = true;
        is_seed_Y(sample_Y_i) = true;
        is_seed_X(non_matched_X_i) = false;
        is_seed_Y(non_matched_Y_i) = false;
        
    end
    
    assert( any(~is_sampled_X) )
    
    % show current result
    seed_X = find(is_seed_X);
    seed_Y = matches_out(seed_X);
    if isfield(params, 'vis') && params.vis
        vis_partial_matches(X, Y, seed_X, seed_Y),
    end
    
    if step > 1 && ~is_updated_matches
        break;
    end
    
    num_cells = min(num_cells*params.k, params.max_num_cells);    
    
    % sample more points to increase scale by k_factor
    [sample_X, sample_Y] = sample_hires_pts(X, Y, is_seed_X, is_seed_Y, ...
        is_forbidden_X, is_forbidden_Y, params.k);
    is_sampled_X(sample_X) = true;
    is_sampled_Y(sample_Y) = true;
    
    % split sampled points on X to cells
    max_problem_size = Inf;
    if ~isfield(X.lores, 'D')
        try
            mesh = geodesic_new_mesh(X.lores.VERT, X.lores.TRIV);
        catch ex
            geodesic_delete;
            rethrow(ex);
        end
    else
        mesh = [];
    end
    fac = 1;
    while max_problem_size > params.problemSize
        num_cells = ceil(num_cells * fac);
        fac = 1.05;
        [fx, indsX] = fast_voronoi(X, X.lores, num_cells, seed_X, sample_X, mesh);
        max_problem_size = max(accumarray(fx, 1));
    end
    if ~isfield(X.lores, 'D')
        geodesic_delete;
    end
    
    cell_index_X(indsX) = fx;
        
    % transfer voronois of X to Y via seeds
    [fy, indsY] = fast_voronoi(Y, Y.lores, numel(seed_Y), seed_Y, sample_Y);
    y2x_cell_ind = cell_index_X(seed_X);
    cell_index_Y(indsY) = y2x_cell_ind(fy);
    

    step = step+1;
    
end





end


function [sample_X, sample_Y] = sample_hires_pts(X, Y, is_seed_X, is_seed_Y, is_forbidden_X, is_forbidden_Y, k)
    unused_pts = find(~is_seed_X & ~is_forbidden_X);   
    if numel(unused_pts) == 0
       sample_X = [];
       sample_Y = [];
       return;
    end
    num_pts_to_sample = min(numel(unused_pts), (k-1)*sum(is_seed_X));
    sample_X = fps_euclidean(X.VERT(unused_pts,:), num_pts_to_sample, randi(numel(unused_pts)));    
    sample_X = unused_pts(sample_X);
    unused_pts = find(~is_seed_Y & ~is_forbidden_Y); 
    num_pts_to_sample = min(numel(unused_pts), (k-1)*sum(is_seed_Y));
    
    if numel(unused_pts) == 0
       sample_X = [];
       sample_Y = [];
       return;
    end
    sample_Y = fps_euclidean(Y.VERT(unused_pts,:), num_pts_to_sample, randi(numel(unused_pts)));
    sample_Y = unused_pts(sample_Y);
end