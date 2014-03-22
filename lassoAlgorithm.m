function beta = lassoAlgorithm(X, Y, rho, maxLoops, minPrecision)

% Initialization
p = size(X,2);
beta = zeros(p,1);
b_old = beta;
i = 0;

% Precompute X'X and X'Y
XTX = X'*X;
XTY = X'*Y;

% Shooting loop
while i < maxLoops,
    i = i+1;
    for j = 1:p,
        jminus = setdiff(1:p,j);
        S0 = XTX(j,jminus)*beta(jminus) - XTY(j);  % S0 = X(:,j)'*(X(:,jminus)*b(jminus)-Y)
        if S0 > rho,
            beta(j) = (rho-S0) / norm(X(:,j),2)^2;
        elseif S0 < -rho,
            beta(j) = -(rho+S0) / norm(X(:,j),2)^2;
        else
            beta(j) = 0;
        end
    end
    delta = norm(beta-b_old,1);    % Norm change during successive iterations
    if delta < minPrecision, break; end
    b_old = beta;
end
if i == maxLoops,
    fprintf('%s\n', 'Maximum number of iteration reached, shooting may not converge.');
end