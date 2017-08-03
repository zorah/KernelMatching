function [matches, non_matched_X, non_matched_Y] = mfilter(desc_X, desc_Y, ...
    evecsD_X, evecsD_Y, evals_X, evals_Y, seed_X, seed_Y,...
    sample_X, sample_Y, params, n, mu)
% MFILTER Calculates a (possible partial) matching between X and Y.  
% Inputs:
%   X, Y - mesh struct with .desc
%   seed_X, seed_Y - indices of sparse starting correspondence, or []
%   sample_X, sample_Y - indices of points to be matched, determines the
%                        complexity of the LAP
%   params - struct with .n_iter, .lambda
% Output:
%   matches - k x 2, indices on X (first column) matching to indices on Y
%             (second column)
%   non_matched_X, non_matched_Y - if partial, indices of unmatched points
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
    
    if ~exist('mu', 'var')
        mu = 1;
    end
    if ~exist('n', 'var')
        n = 1;
    end
    if ~exist('params', 'var')
        params = struct;
    end
    if ~isfield(params, 'maxIter')
        params.maxIter = 10;
    end
    if ~isfield(params, 'lambda')
        params.lambda = 1e0;
    end
    if ~isfield(params, 't_it')
        params.t_it = @(t,k,i) (t(i));
    end

    matches = zeros(numel([seed_X; sample_X]),1);
    non_matched_X = [];
    non_matched_Y = [];
          
    n_seeds = numel(seed_X);    
    sample_X = [seed_X; sample_X];
    sample_Y = [seed_Y; sample_Y];
    
    if strcmp(params.filter, 'hk')
       seed_X_out = (length(sample_X)+1):size(evecsD_X, 1);
       seed_Y_out = (length(sample_Y)+1):size(evecsD_Y, 1);
       kernelfun = @(xin, yin, tX, tY)kde_HK(evecsD_X, evecsD_Y, evals_X, ...
           evals_Y, xin, yin, seed_X_out, seed_Y_out, tX, tY);   
    elseif strcmp(params.filter, 'gauss')
       kernelfun = @(xin, yin, varX, varY)kde_Gauss(evecsD_X, evecsD_Y, xin, yin, params.sigma); 
    else
       error('No kernel method chosen.'); 
    end
    
    xin = 1:n_seeds;
    yin = 1:n_seeds;
    
    F1 = desc_X*desc_Y';
    desc_dim = size(desc_X,2);
    
    F2 = kernelfun(xin, yin, params.t_it(params.tX,n,params.maxIter), params.t_it(params.tY,n,params.maxIter));
    F1 = F1 + params.lambda*(desc_dim/numel(sample_X))*F2;
    
    for i_iter = 1:params.maxIter        
        F2 = kernelfun(xin, yin, params.t_it(params.tX,n,i_iter), params.t_it(params.tY,n,i_iter));
        
        F = mu * F1 + params.lambda*(desc_dim/numel(sample_X))*F2;
        
        if isfield(params, 'convexify') && params.convexify 
           Pi=sparse(xin,yin,params.rho,size(F,1),size(F,2));
           F = F + double(Pi);
        end
        
        Fmin = min(min(F));
        
        
        % if sample_x < sample_y
        if size(F,1) <= size(F,2)
            G = [F; 0.3 * Fmin * rand(size(F,2) - size(F,1), size(F,2)) + eps]; 
            assignment = assignmentAlgs(G,'auction');
            xin = 1:size(F,1);
            yin = assignment(xin);            
            non_matched_Y = assignment(size(F,1)+1:end);
            if ~isempty(non_matched_Y)
                non_matched_Y = sample_Y(non_matched_Y);            
            end
        else % sample_x > sample_y
            G = [F, 0.3 * Fmin * rand(size(F,1), size(F,1) - size(F,2)) + eps];  
            assignment = assignmentAlgs(G,'auction');
            yin = (assignment(assignment <= size(F,2)))';
            xin = (find(assignment <= size(F,2)))';
            non_matched_X = setdiff(1:size(F,1), xin);
            non_matched_X = sample_X(non_matched_X);         
        end
        
            
    end    

    matches(xin) = sample_Y(yin);
     

end

% calculates the inner product of heat kernels
function [ F ] = kde_HK(Qx, Qy, evals_X, evals_Y, xin, yin, seed_X_out, ...
    seed_Y_out, tX, tY)
    nu = 1;
    if size(xin, 2) == 1
        xin = xin';
    end
    if size(yin, 2) == 1
        yin = yin';
    end
    F = Qx * diag(exp(-tX*abs(evals_X))) * ([Qx(xin,:); nu * ...
        Qx(seed_X_out, :)])' * [Qy(yin,:); nu * Qy(seed_Y_out, :)] * ...
        diag(exp(-tY*abs(evals_Y))) * Qy';
    F(seed_X_out, :) = [];
    F(:, seed_Y_out) = [];
    F = F - min(min(F)) + 1e-3;
end

% calculates the inner product of gauss kernels
function [ F ] = kde_Gauss(XD,YD, xin, yin, sigma)
    Kx = kernel(XD, xin, sigma);
    Ky = kernel(YD, yin, sigma);
    F = Kx*Ky'; 
    F = max(0,F) + 1e-5;
end

