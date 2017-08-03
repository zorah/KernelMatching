clc, clear all, close all,
addpath(genpath('.')),

%% load high res

file1 = 'cat0';
file2 = 'cat6';

X = load_off(strcat('./data/', file1, '.off'));
X.n = size(X.VERT, 1);
X.m = size(X.TRIV, 1);
X.area = sum(calc_tri_areas(X));

Y = load_off(strcat('./data/', file2, '.off'));
Y.n = size(Y.VERT, 1);
Y.m = size(Y.TRIV, 1);
Y.area = sum(calc_tri_areas(Y));

%% add low res
% X.lores, Y.lores should be triangle meshes with around 5000 vertices
% these are not strictly necessary and can be replaced with the original
% meshes but it will slow down the code considerably.
% We calculate geodesics on these which is sensitive to faulty meshes
% and might be the cause of crashes. 

X.lores = load_off(strcat('./data/', file1, '_5k.off'));
Y.lores = load_off(strcat('./data/', file2, '_5k.off'));

%% options
opts = struct;

opts.lambda = 1e7; % lambda weighting the regularizer, alpha = 1/lambda
opts.maxIter = 10; % max iterations in each filter step
opts.k = 3; % factor increasing the sampling in each multi scale iteration
opts.problemSize = 1500; % max LAP size solved
opts.problemSizeInit = 3000; % max LAP size for first iteration (this might 
                             % be bigger than problemSize, if too low bad 
                             % results will propagate too much through the multi scale, but too big will be slow)                         
opts.tX = logspace(log10(500),log10(10),opts.maxIter); % time parameter for the heat kernels on X
opts.tY = logspace(log10(500),log10(10),opts.maxIter); % time parameter for the heat kernels on Y
opts.n_evecs = 500; % number of eigenfunctions computed

opts.partial = false; % allows entire parts of the shape to be not matched, based on difference of area
opts.vis = false; % shows (or not) the result of each multiscale iteration
opts.oneIteration = false; % only samples problemSizeInit many points and
                           % returns without multiscale

% these options can be changed but were fixed in all our experiments
opts.shot_num_bins = 10; % number of bins for shot
opts.shot_radius = 5; % percentage of the diameter used for shot

opts.mu = 1; % mu weighting the data term after the first iteration
opts.t_it = @(t, k, i) (t(i) ./ 1); % function applied to tX/tY in each iteration

opts.filter = 'hk'; % 'gauss' is the other option but it will not work out-of-the-box
opts.use_par = false; % whether to use parfor while working on each scale of the problem
                      % leads to memory problems rather often

%% precalc

% any descriptors you use should be concatinated in X.desc/Y.desc
X.desc = calc_shot(X.VERT', X.TRIV', 1:X.n, opts.shot_num_bins, opts.shot_radius*sqrt(X.area)/100, 3)';
Y.desc = calc_shot(Y.VERT', Y.TRIV', 1:Y.n, opts.shot_num_bins, opts.shot_radius*sqrt(Y.area)/100, 3)';

[X.evecs,X.evals] = HeatKernels(X,1:X.n,opts.n_evecs);
[Y.evecs,Y.evals] = HeatKernels(Y,1:Y.n,opts.n_evecs);


%% run code
a=tic;
matches = run_PMF(X, Y, opts);
toc(a)

%% visualize

figure,
vis_partial_matches(X, Y, matches(:,1), matches(:,2));
% grey dots on either side indicate unmatched points

