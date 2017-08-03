function [Pi,cost]=assignmentAlgs(C,alg,auctionScale,opts)
%solves max trace(C.'Pi) = -min trace(-C.'Pi)

if nargin<3 || isempty(auctionScale)
    max_payoff_value = 10^7;
    auctionScale = max_payoff_value/max(C(:));
%     auctionScale=10^4;
    
end

n=size(C,1);

switch alg
    case 'munkres'
        
        [assignment,cost] = munkres(-(C));
        Pi=sparse(1:n,assignment,1,n,n);
        
    case 'auction'
        verbose = 0;
        C_benefits=sparse(C);
        C_benefits=C_benefits-min(C_benefits(:))+10^-3;
        C_benefits=C_benefits*auctionScale;
        [assignment,~] = sparseAssignmentProblemAuctionAlgorithm(C_benefits, [], [], verbose);
        
        Pi=assignment;
        
    case 'linprog'
        %C_costs=sparse(-C);
        %C_costs=C_costs-min(C_costs(:));
        [Pi,cost]=linProgAssignmentProblem(-C);
        
    case 'lapjv'
        [assignment,cost,v,u,costMat] = lapjv(-C);
        Pi=sparse(1:n,assignment,1,n,n);

    case 'Sinkhorn'
        done = false;
        while ~done            
            [Pi,exitflag] = optimalTransportWithEntropicReg(-C,ones(size(C,1),1),ones(size(C,2),1),opts);
            if ~exitflag
                opts.lambda = opts.lambda*opts.betta;
                disp(['Numerical inaccuracies. Increasing regularization parameter to lambda = ',num2str(opts.lambda)]);
            else
                done = true;
            end
        end
end