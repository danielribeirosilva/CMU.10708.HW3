% QUESTION 1 - Neighborhood Selection


%params
%rhoList = [0,35,37,38,40,41,42];
%rhoList = [0, 20, 30, 40];
rhoList = [0,1,5,10,30];

%load data
load('ggm_data.mat');
n = size(X,1);

%center distribution
Xcentered = bsxfun(@minus, X, mean(X));

%compute S
S = (1/n)*Xcentered'*Xcentered;
theta = inv(S);

p = size(S,1);


i = 1;
for rho = rhoList
    
    thetaNew = eye(p);
    
    for c = 1:p
       A = theta(:,1:p~=c);
       B = theta(:,c);
       beta = lasso(A,B,rho);
       thetaNew(1:p~=c,c)=beta;
    end
    
    %symmetrize
    thetaNew = thetaNew & thetaNew';

    %adjust colors
    thetaNew = ~ thetaNew;
    
    %plot
    figure(i);
    imshow(thetaNew,'InitialMagnification','fit');
    
    i=i+1;
end


