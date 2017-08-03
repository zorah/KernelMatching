function [ matches_out ] = pmf_multiscale( X, Y, seed_X, seed_Y, ...
    ind_forbidden_X, ind_forbidden_Y, params)
%MULTISCALE Solving the product manifold filter in a multiscale approach 
% using multiple threads. 
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

if strcmp(params.filter, 'gauss')
   error('Gauss filters and parallel multiscale do not work together at the moment.'), 
end

max_num_iter = ceil(log(X.n / params.problemSize) / log(params.k)) + 5;

%% initialize and add samples by fps

seed_X_orig = seed_X;
seed_Y_orig = seed_Y;

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
while step <= max_num_iter && any(~(is_forbidden_X | is_seed_X)) ...
        && any(~(is_forbidden_Y | is_seed_Y))
    fprintf('starting step: %d \n',step)
    % solve problem at the current scale
    
    is_updated_matches_i = {};
    for i_cell = 1:num_cells  
        % find sample points in this cell
        sample_X_i{i_cell} = find( (cell_index_X == i_cell) & is_sampled_X);
        sample_Y_i{i_cell} = find( (cell_index_Y == i_cell) & is_sampled_Y);  
        
        % check for convergence
        
        % get all seeds in this cell
        seed_X_i{i_cell} = find( (cell_index_X == i_cell) & is_seed_X);
        seed_Y_i{i_cell} = matches_out(seed_X_i{i_cell});
        
        % check for convergence
        if isempty(sample_X_i{i_cell}) || isempty(sample_Y_i{i_cell})
            is_updated_matches_i{i_cell} = false;
        else
            is_updated_matches_i{i_cell} = true;
        end
    end
    
    % pre assemble input for the filter
    X_desc = cell(num_cells, 1);
    Y_desc = cell(num_cells, 1);
    X_evecs = cell(num_cells, 1);
    Y_evecs = cell(num_cells, 1);
    X_evals = cell(num_cells, 1);
    Y_evals = cell(num_cells, 1);
    for k=1:num_cells
       X_desc{k} = X.desc([seed_X_i{k}; sample_X_i{k}], :);
       Y_desc{k} = Y.desc([seed_Y_i{k}; sample_Y_i{k}], :);
       X_evecs{k} = X.evecs([seed_X_i{k}; sample_X_i{k}; seed_X_orig], :);
       Y_evecs{k} = Y.evecs([seed_Y_i{k}; sample_Y_i{k}; seed_Y_orig], :);
       X_evals{k} = X.evals;
       Y_evals{k} = Y.evals;
    end
    
    parfor i_cell = 1:num_cells  
        % solve matching in this cell on {seeds, samples}
        
        if is_updated_matches_i{i_cell} == false
            continue
        end
        
        tic,
        [matches{i_cell}, non_matched_X_i{i_cell}, non_matched_Y_i{i_cell}] = ...
            mfilter(X_desc{i_cell}, Y_desc{i_cell}, X_evecs{i_cell}, Y_evecs{i_cell}, ...
            X_evals{i_cell}, Y_evals{i_cell}, seed_X_i{i_cell}, seed_Y_i{i_cell}, ...
            sample_X_i{i_cell}, sample_Y_i{i_cell}, params, step, params.mu);
        time_elapsed = toc;
        
        fprintf('step: %d, i_cell: %d of %d, num_seeds: %d, num_samples_X: %d, num_samples_Y: %d, time: %.2f \n',...
                      step, i_cell, num_cells, numel(seed_X_i{i_cell}), numel(sample_X_i{i_cell}) , numel(sample_Y_i{i_cell}), time_elapsed);
    end
    
    is_updated_matches = false;
    for i_cell = 1:num_cells 
        
        if is_updated_matches_i{i_cell} == false
            continue;
        end

        matches_out([seed_X_i{i_cell}; sample_X_i{i_cell}]) = matches{i_cell}; % mind the ordering
        is_updated_matches = is_updated_matches || is_updated_matches_i{i_cell};
        % make sample points the new seed points
        is_sampled_X(sample_X_i{i_cell}) = false;
        is_sampled_Y(sample_Y_i{i_cell}) = false;
        sample_X_i{i_cell} = setdiff(sample_X_i{i_cell}, non_matched_X_i{i_cell});
        sample_Y_i{i_cell} = setdiff(sample_Y_i{i_cell}, non_matched_Y_i{i_cell});
        is_seed_X(sample_X_i{i_cell}) = true;
        is_seed_Y(sample_Y_i{i_cell}) = true;
        is_seed_X(non_matched_X_i{i_cell}) = false;
        is_seed_Y(non_matched_Y_i{i_cell}) = false;        
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