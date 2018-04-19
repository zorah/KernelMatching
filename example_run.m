clc, clear all, close all,
addpath(genpath('.')),

%% load high res
% about fixing the meshes: not really necessary but the exactgeodesic code
% is really picky, it might even fail after fixing. Feel free to suggest
% a different code (we actually only use dijkstra but it should be really
% fast)

file1 = 'cat0';
file2 = 'cat6';

% kids
% X = load_off(strcat('./data/', file1, '.off'));
% %X = remesh_and_fix(X);
% X.n = size(X.VERT, 1);
% X.m = size(X.TRIV, 1);
% Y = load_off(strcat('./data/', file2, '.off'));
% %Y = remesh_and_fix(Y);
% Y.n = size(Y.VERT, 1);
% Y.m = size(Y.TRIV, 1);

% tosca
% X = load_off('~/Shapes/shrec2016_Partial/null/cat.off');
% load(strcat('~/Shapes/toscahires-mat/', file1, '.mat'));
% X.lores = struct;
% X.lores.VERT = [surface.X, surface.Y, surface.Z];
% X.lores.TRIV = surface.TRIV;
% Y = load_off(strcat('~/Shapes/shrec2016_Partial/holes/', file2, '.off'));
% load(strcat('/work/laehner/SGP2017/data/SHREC16 remeshed/holes/', file2, '.mat'));
% Y.lores = N_remesh;
% Y.lores = remesh_and_fix(Y.lores);
X = load_off(strcat('./data/', file1, '.off'));
X.n = size(X.VERT, 1);
X.m = size(X.TRIV, 1);
X.lores = load_off(strcat('./data/', file1, '_5k.off'));

Y = load_off(strcat('./data/', file2, '.off'));
Y.n = size(Y.VERT, 1);
Y.m = size(Y.TRIV, 1);
Y.lores = load_off(strcat('./data/', file2, '_5k.off'));

% load(strcat('~/Shapes/toscahires-mat/', file1, '.mat'));
% X = struct;
% X.VERT = [surface.X, surface.Y, surface.Z];
% X.TRIV = surface.TRIV;
% X.lores = X;
% load(strcat('~/Shapes/toscahires-mat/', file2, '.mat'));
% Y = struct;
% Y.VERT = [surface.X, surface.Y, surface.Z];
% Y.TRIV = surface.TRIV;
% Y.lores = Y;

% faust
% X = struct;
% [X.VERT, X.TRIV] = read_ply(strcat('~/Shapes/MPI-FAUST/training/registrations/tr_reg_', file1, '.ply'));
% load(strcat('/work/laehner/SGP2017/data/faust_reg_models_combined/tr_reg_', file1, '.mat'));
% X.lores = X;
% Y = struct;
% [Y.VERT, Y.TRIV] = read_ply(strcat('~/Shapes/MPI-FAUST/training/registrations/tr_reg_', file2, '.ply'));
% load(strcat('/work/laehner/SGP2017/data/faust_reg_models_combined/tr_reg_', file2, '.mat'));
% Y.lores = Y;

X.n = size(X.VERT, 1);
X.m = size(X.TRIV, 1);
Y.n = size(Y.VERT, 1);
Y.m = size(Y.TRIV, 1);

% FOR SHAPES NOT IN TOSCA
% refarea = 1.93e+04; % time parameters are adjusted to a certain surface area
%                     % this is used to rescale the input shapes
% X.VERT = X.VERT ./ sqrt(sum(calc_tri_areas(X))) .* sqrt(refarea);
% Y.VERT = Y.VERT ./ sqrt(sum(calc_tri_areas(Y))) .* sqrt(refarea);


%%
% X.VERT = 100 * X.VERT;
% X.lores.VERT = 100 * X.lores.VERT;
% Y.VERT = 100 * Y.VERT;
% Y.lores.VERT = 100 * Y.lores.VERT;

%% add low res
% X.lores, Y.lores should be triangle meshes with around 5000 vertices
% either precompute them or comment in qslim in remesh_and_fix


% load(strcat('./data/', file1, '.mat'));
% X.lores = model;
% X.lores = remesh_and_fix(X.lores);
% load(strcat('./data/', file2, '.mat'));
% Y.lores = model;
% Y.lores = remesh_and_fix(Y.lores);

X.area = sum(calc_tri_areas(X));
Y.area = sum(calc_tri_areas(Y));

%% options
opts = struct;
opts.shot_num_bins = 10; % number of bins for shot
opts.shot_radius = 5; % percentage of the diameter used for shot

opts.lambda = 1e7; % lambda weighting the regularizer
opts.mu = 1; % mu weighting the data term after the first iteration
opts.maxIter = 10; % max iterations in each filter step
opts.k = 3; % factor increasing the sampling in each multi scale iteration
opts.problemSize = 1000; % max LAP size solved
opts.problemSizeInit = 3000; % max LAP size for first iteration (this might 
                             % be bigger than problemSize, if too low bad 
                             % results will propagate too much through the multi scale, but too big will be slow)
opts.use_par = false; % whether to use parfor while working on each scale of the problem.                             
opts.tX = logspace(log10(500),log10(10),opts.maxIter); % time parameter for the heat kernels on X
opts.tY = logspace(log10(500),log10(10),opts.maxIter); % time parameter for the heat kernels on Y
opts.t_it = @(t, k, i) (t(i) ./ 1);
% opts.tX = logspace(log10(500),log10(50),opts.maxIter); % time parameter for the heat kernels on X

opts.partial = false; % allows entire parts of the shape to be not matched, based on difference of area
opts.vis = false; % shows (or not) the result of each multiscale iteration
opts.n_evecs = 500; % number of eigenfunctions computed
opts.oneIteration = false;

opts.convexify = false;
opts.rho = 0.5;

opts.filter = 'hk';

%% precalc

% any descriptors you use should be in X.desc/Y.desc, normalized
X.desc = calc_shot(X.VERT', X.TRIV', 1:X.n, opts.shot_num_bins, opts.shot_radius*sqrt(X.area)/100, 3)';
Y.desc = calc_shot(Y.VERT', Y.TRIV', 1:Y.n, opts.shot_num_bins, opts.shot_radius*sqrt(Y.area)/100, 3)';

[X.evecs,X.evals] = HeatKernels(X,1:X.n,opts.n_evecs);
[Y.evecs,Y.evals] = HeatKernels(Y,1:Y.n,opts.n_evecs);


%% run code
a=tic;
matches = run_PMF(X, Y, opts);
toc(a)

%% visualize
% figure
vis_partial_matches(X, Y, matches(:,1), matches(:,2))
% figure
% thresh = 0:.001:0.2;
% Y.diameter = sqrt(Y.area);
% gt_matches = load('data/kid02_ref.txt');
% gt_matches = gt_matches(:,[2 1]);
% disp('Calculating errs')
% errs = geoerr(Y,matches,gt_matches,Y.diameter);
% disp('done')
% plot(thresh,calc_err_curve(errs,thresh));
