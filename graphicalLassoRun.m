% QUESTION 2 - Graphical Lasso


%params
%rhoList = [0 0.2 0.5 0.8];
rhoList = [0,0.2,0.38,0.4,0.41,0.8];
maxLoops = 100;
minPrecision = 0.000001;

%load data
load('ggm_data.mat');

[n,p] = size(X);

%center distribution
Xcentered = bsxfun(@minus, X, mean(X));

%compute S
S = (1/n)*Xcentered'*Xcentered;

%run graphical lasso
i=1;
for rho = rhoList
    fprintf('current rho = %f\n',rho);
    theta  = graphicalLassoAlgorithm(S, rho, maxLoops, minPrecision);
    
    %adjust colors
    theta = theta & theta;
    theta = ~ theta;
    %disp(theta)
    
    %plot
    figure(i);
    imshow(theta,'InitialMagnification','fit');
    
    i=i+1;
end