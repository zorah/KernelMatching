function [HKS,t] = calc_HKS(M, n, t)

if nargin==2
    % based on J.Sun's code
    tmin = abs(4*log(10) / -abs(M.evals(n)));
    tmax = abs(4*log(10) / -abs(M.evals(2)));
    nstep = 100;
    stepsize = (log(tmax) - log(tmin)) / nstep;
    logts = log(tmin):stepsize:log(tmax);
    t = exp(logts);
    %t = 2.^(1: 1/16 : 18-1/16);
end

HKS = zeros(M.n, length(t));

% skipping the first freq. as the authors do in their code
for i=1:length(t)
    HKS(:,i) = sum(...
        (M.evecs(:,2:n)).^2 .* repmat(exp(t(i)*(-abs(M.evals(2:n))))',M.n,1), ...
        2);
end

end
