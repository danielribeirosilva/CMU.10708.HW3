% Function to run lasso by minimizing the following objective:
% f(beta) = g(beta) + h(beta)
% where g(beta) = 0.5*||y - X*beta||_2^2 
% and h(beta) = lambda*||beta||_1

% Inputs:
%   X: n-by-p matrix of input features
%   y: n-by-1 vector of output values
%   lambda: weight for l1 norm penalty

% Outputs:
%   beta: estimate of regression parameters

function beta = lasso(X,y,lambda)

% set options for optimization
miniter = 10;
maxiter = 5000;
tol = 1e-5;

% set fixed step size
t = 0.0001;

% initialize beta
p = size(X,2);
beta = zeros(p,1);

% optimize lasso objective 
for iter = 1:maxiter
    % compute gradient of g(beta) with respect to beta
    grad = X'*X*beta - X'*y;
    % compute new estimate of beta by appling soft thresholding operator
    beta_new = softthresh(beta-t*grad,t*lambda);
    
    % compute new objective value
    obj_new = sum((y-X*beta_new).^2)/2 + lambda*sum(abs(beta_new));
    % check convergence
    if (iter >= miniter && abs(obj_new-obj) < tol)
        break;
    end
    % store new variables
    beta = beta_new;
    obj = obj_new;
    
end

%disp(beta)

end

% Soft thresholding
function alpha = softthresh(beta,eta)
    
    beta(beta>eta) = beta(beta>eta) - eta;
    beta(beta<-eta) = beta(beta<-eta) + eta;
    beta(beta >= -eta & beta <= eta) = 0;
    
    alpha = beta;
end

