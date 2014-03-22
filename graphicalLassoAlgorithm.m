function [theta]  = graphicalLassoAlgorithm(S, rho, maxLoops, minPrecision)

%dimension of the covariance matrix
dimS = size(S,1);
allRows = 1:dimS;
% Step 1
currentW = S + rho * eye(dimS);

% Step 2
previousW = currentW;
i = 1;
while true
    for targetRow = dimS:-1:1
        allButTargetRow = allRows(allRows~=targetRow);
        S_11 = currentW(allButTargetRow,allButTargetRow);
        s_12 = S(allButTargetRow,targetRow);
        [V D] = eig(S_11);
        sqrtD = D.^0.5;
        % compute W_11^(1/2)
        sqrtW_11 = V * sqrtD * V';        
        % compute W_11^(-1/2)
        invSqrtD = 1./sqrtD;
        invSqrtD(isinf(invSqrtD)) = 0;
        invSqrtW_11 = V * invSqrtD * V';
        % compute b = W_11^(-1/2) * s_12
        b = invSqrtW_11 * s_12;
        % use Lasso to solve
        %beta = lassoAlgorithm(sqrtW_11, b, rho, maxLoops, minPrecision);
        beta = lasso(sqrtW_11, b, rho);
        
        
        % replace result on matrix
        new_w_12 = currentW(allButTargetRow,allButTargetRow) * beta;
        currentW(allButTargetRow,targetRow) = new_w_12;
        currentW(targetRow,allButTargetRow) = new_w_12';
    end
    
    % stop criterion of precision
    if norm(currentW-previousW,1) < minPrecision, 
        break; 
    end
    
    % stop criterion of max loops
    if i >= maxLoops
        fprintf('Max loops reached\n');
        break;
    end
    
    % prepare for next step
    previousW = currentW;
    i = i+1;
end

fprintf('Total loops: %i\n',i);

theta = inv(currentW);

