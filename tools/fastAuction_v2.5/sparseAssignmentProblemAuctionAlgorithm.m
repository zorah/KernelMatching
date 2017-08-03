% 
% Mex implementation of Bertsekas' auction algorithm [1] for a very fast
% solution of the linear assignment problem.
%
% The implementation is optimised for sparse matrices where an element
% A(i,j) = 0 indicates that the pair (i,j) is not possible as assignment.
% Solving a sparse problem of size 950,000 by 950,000 with around
% 40,000,000 non-zero elements takes less than 8 mins. The method is also
% efficient for dense matrices, e.g. it can solve a 20,000 by 20,000
% problem in less than 3.5 mins.
%
% In practice, with respect to running time, the auction algorithm
% outperforms the Kuhn-Munkres (or Hungarian) algorithm significantly. The
% auction algorithm has average run-time complexity O(N^2*log(N)) and
% worst-case complexity O(N^3), whereas the Kuhn-Munkres algorithm has
% average- and worst-case complexity O(N^3).
%
% Note that only global optima are found for integer-valued benefit
% matrices. For real-valued benefit matrices a scaling of the values needs
% to be applied by multiplication with a large number. This scaling factor
% depends on the desired accuracy, as the global solution is found for the
% integral part of the benefit matrix, whilst there is no guarantee that
% the fractional part of the benefits are properly taken into account.
% However, in practical cases it seems advantageous to not round the
% resulting benefit matrix and retain the fractional parts in the benefit
% matrix. Also note that larger scalings of the benefit matrix increase the
% run-time, so a problem-specific trade-off between runtime and accuracy
% must be chosen.
%
% Input:
%               A                       N-by-N sparse benefit matrix
%                                       (higher values indicate a better 
%                                       match, the value 0 indicates 
%                                       inadmissable assignments)
%               [epsilon]               Initial value of epsilon (optional)
%               [epsilonDecreaseFactor] Decrease factor of epsilon
%                                       (optional)
%               [verbosity]             level of verbosity (0: quiet, 1:
%                                       general infos, 2: full infos)
%                                       (optional)
%
% Output:
%               assignments             Resulting assignments. If there is
%                                       no feasible solution, all
%                                       assignments are -1.
%               [P]                     Permutation matrix output, such
%                                       that trace(P*A') gives the optimal
%                                       value
%               [prices]                Prices used during auctions
%                                       (optional)
%
% Example:      See the function test() below for a usage example.
%               Typically only the benefit matrix A is given as input and
%               the first or second output argument is relevant. epsilon
%               and epsilonDecreaseFactor can be used to heuristically
%               adapt runtime.
%
% Compilation:  mex -largeArrayDims auctionAlgorithmSparseMex.cpp -lut
%           
% When using this implementation in your work, in addition to [1], please
% cite our paper [2].
%
% [1]	Bertsekas, D.P. 1998. Network Optimization: Continuous and Discrete
%       Models. Athena Scientific.
%
% [2]   Bernard, F., Vlassis, N., Gemmar, P., Husch, A., Thunberg, J.,
% 	    Goncalves, J. and Hertel, F. 2016. Fast correspondences for
% 	    statistical shape models of brain structures. SPIE Medical Imaging,
% 	    San Diego, CA, 2016.
%
% Implementation by Florian Bernard ( f.bernardpi [at] gmail [dot] com ).
%
% The author would like to thank Gary Guangning Tan for helpful feedback.
%
% Last modified on 15/12/2016
%

function [assignments, P, prices] = ...
	sparseAssignmentProblemAuctionAlgorithm(A, epsilon, ...
	epsilonDecreaseFactor, verbosity)

	N = size(A,1);

	if ( any(A(:)<0) )
		error('Only non-negative benefits allowed!');
	end

	if ( ~issparse(A) )
% 		warning('Converting A to sparse matrix!');
		A = sparse(A);
	end

	% heuristic for setting epsilon
	A = A*(N+1);
	maxAbsA = full(max(abs(A(:))));
	if ( ~exist('epsilon', 'var') || isempty(epsilon) )
		epsilon = max(maxAbsA/50, 1);
% 		epsilon = 0.5*((N*maxAbsA)/5 + N*maxAbsA); % see page 260 in [1]
	end

	if ( ~exist('epsilonDecreaseFactor', 'var') || isempty(epsilonDecreaseFactor) )
		epsilonDecreaseFactor = 0.2;
	end
	
	if ( ~exist('verbosity', 'var') )
		verbosity = 0;
	end
	
	% if the diagonal element of A is infeasible, we augment it to make
	% sure there exists a feasible solution...
	infeasibleDiags = find(diag(A)==0);
	linIdx = sub2ind([N,N], infeasibleDiags, infeasibleDiags);
	A(linIdx) = -(2*N-1)*maxAbsA-1; % see page 266 in [1]
	
	[assignments, prices] = ...
		auctionAlgorithmSparseMex(A', epsilon, epsilonDecreaseFactor, ...
		maxAbsA, verbosity);
	
	if ( all(assignments<0) )
		warning('No feasible solution exists');
		P = [];
		return;
	end

	P = sparse(1:N, assignments', ones(1,N), N,N);
	
	% now, we check if the solution uses any of the infeasible elements on
	% the diagonal - if that is the case, no feasible solution exists
	diagP = diag(P);
	if ( any(diagP(infeasibleDiags)) )
		warning('No feasible solution exists');
		assignments = -ones(N,1);
		P = [];
	end

end


function test()
%% DEMO
	% compile mex file
	mex -largeArrayDims auctionAlgorithmSparseMex.cpp -lut
	
	% create sample data
	N = 2000;
	
	A = rand(N,N);
	
	% create sparse matrix, since the mex implementation uses the Matlab
	% sparse matrix data structure
	A = sparse(A);

	% scale A such that round(Ascaled) has sufficient accuracy
	scalingFactor = 10^6;
	Ascaled = A*scalingFactor;
	
	% solve assignment problem
	tic
	[assignments,P] = sparseAssignmentProblemAuctionAlgorithm(Ascaled);
	toc
end
