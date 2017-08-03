function [ S ] = fps_euclidean(V, n, seed)
%EUCLIDEANFPS Samples K vertices from V by using farthest point sampling.
% The farthest point sampling starts with vertex v1 and uses the euclidean
% metric of the 3-dimensional embedding space.
% -  V is a n-by-3 matrix storing the positions of n vertices
% -  K is the number of samples
% -  v1 is the index of the first vertex in the sample set. (1<=v1<=n)
% Returns
% -  S is a K-dimensional vector that includes the indeces of the K sample
%    vertices.

%Hint: matlab function pdist2 could be helpful

if n > size(V, 1)
    n = size(V, 1);
end

S = zeros(n,1);
S(1) = seed;
d = pdist2(V,V(seed,:));

for i=2:n
    [~,m] = max(d);
    S(i) = m(1);
    d = min(pdist2(V,V(S(i),:)),d);
end

end

