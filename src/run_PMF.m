function [ correspondences ] = run_PMF( X, Y, opts)
% This is a wrapper for all possible PMF configurations. Use opt to select
% which configuration you need. 
% inputs:
%   X, Y - struct with .VERT, .TRIV, .evecs, .evals, and .desc containing descriptors
%   opts - struct with .lambda, .tX, .tY, .maxIter, .n_seeds
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

%% depenencies
global geodesic_library;
geodesic_library = 'geodesic_matlab_api';

%% params
if ~isfield(opts, 'use_par')
    opts.use_par = false;
end
if ~isfield(opts, 'lambda')
    opts.lambda = 1e0;
end
if ~isfield(opts, 'maxIter')
    opts.maxIter = 10;
end
if ~isfield(opts, 'problemSize')
    opts.problemSize = max(X.n, Y.n);
end
if ~isfield(opts, 'problemSizeInit')
    opts.problemSizeInit = opts.problemSize;
end
if ~isfield(opts, 'tX')
    opts.tX = 5 * X.area;
end
if ~isfield(opts, 'tY')
    opts.tY = 5 * X.area;
end
if ~isfield(opts, 'sigma')
    opts.sigma = 0.5 * X.area;
end
if ~isfield(opts, 'uniformSampling')
    opts.uniformSampling = false;
end

problemSize = opts.problemSizeInit;

%% Start with solving a coarse (or full) problem
% Choose samples
if isfield(opts, 'partial') && opts.partial
    facX = X.area / max(X.area, Y.area);
    facY = Y.area / max(X.area, Y.area);
    samples_X = fps_euclidean(X.VERT, ceil(facX * problemSize), randi(X.n));
    
    samples_Y = fps_euclidean(Y.VERT, ceil(facY * problemSize), randi(Y.n));
else
    samples_X = fps_euclidean(X.VERT, problemSize, randi(X.n));
    if opts.uniformSampling
        samples_Y = samples_X;
    else
        samples_Y = fps_euclidean(Y.VERT, problemSize, randi(Y.n));
    end
end

if strcmp(opts.filter, 'hk')
    evecsD_X = X.evecs(samples_X, :);
    evecsD_Y = Y.evecs(samples_Y, :);
elseif strcmp(opts.filter, 'gauss')
    evecsD_X = X.D(samples_X, samples_X);
    evecsD_Y = Y.D(samples_Y, samples_Y);
end


[matches, non_matched_X, non_matched_Y] = mfilter(X.desc(samples_X, :), ...
    Y.desc(samples_Y, :), evecsD_X, evecsD_Y, ...
    X.evals, Y.evals, [], [], samples_X, samples_Y, opts);

zero_matches = matches == 0;
matches(zero_matches) = 1;
correspondences = [samples_X, matches];
correspondences(zero_matches, :) = [];

if isfield(opts, 'vis') && opts.vis
    vis_partial_matches(X, Y, correspondences(:,1), correspondences(:,2)),
end

%% Running multi scale

if (isfield(opts, 'oneIteration') && opts.oneIteration) || length(correspondences) >= min(X.n, Y.n)
    return; 
end

% downsampled version of X, Y, only works for windows
if ~isfield(X, 'lores')
    X.lores = remesh(X, struct('vertices', 5e3, 'placement', 0));
end
if ~isfield(Y, 'lores')
    Y.lores = remesh(Y, struct('vertices', 5e3, 'placement', 0));
end

if isfield(opts, 'partial') && opts.partial
   non_sampled_X = (1:X.n)';
   non_sampled_X(samples_X) = [];
   [fx, indsX] = fast_voronoi(X, X.lores, size(correspondences, 1) + numel(non_matched_X), [non_matched_X; correspondences(:,1)], non_sampled_X);
   forbidden_spread = indsX(fx <= numel(non_matched_X));
   non_matched_X = [non_matched_X; forbidden_spread];
   f1 = zeros(X.n, 1);
   f1(non_matched_X) = 1;
   
   non_sampled_Y = (1:Y.n)';
   non_sampled_Y(samples_Y) = [];
   [fy, indsY] = fast_voronoi(Y, Y.lores, numel(correspondences) + numel(non_matched_Y), [non_matched_Y; correspondences(:,2)], non_sampled_Y);
   forbidden_spread = indsY(fy <= numel(non_matched_Y));
   non_matched_Y = [non_matched_Y; forbidden_spread];
   f2 = zeros(Y.n, 1);
   f2(non_matched_Y) = 1;
end

if opts.use_par
    matches = pmf_multiscale_par(X, Y, correspondences(:,1), correspondences(:,2), non_matched_X, non_matched_Y, opts);
else
    matches = pmf_multiscale(X, Y, correspondences(:,1), correspondences(:,2), non_matched_X, non_matched_Y, opts);
end

zero_matches = matches == 0;
matches(zero_matches) = 1;
correspondences = [(1:X.n)', matches];
correspondences(zero_matches, :) = [];

end
