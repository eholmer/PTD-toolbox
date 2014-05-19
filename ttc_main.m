function [CoxLambda,CoxTheta,Cox] = ttc_main(startProb,varargin)%Lambda,Theta,startProb,C)
% Converts a triangular intensity matrix into a coxian representation
% First argument must be the startvector startprob
% If the whole intensity matrix is an argument it must be the second one
% If a generator matrix and absorbing intensities vector is the argument
% the matrix must be the second one and the vector the third one
% The function returns the generator CoxLambda, the absoring intensities
% CoxTheta and the whole intensity matrix Cox

if nargin == 3
    Lambda=varargin{1};
    Theta=varargin{2};
    Q=[Lambda Theta;zeros(1,length(Lambda)+1)];
else if nargin ==2
        Q=varargin{1};
    end
end

n=length(Q)-1; % Number of non-absorbing states
diags=diag(Q); % Sojurn intensities

diagsSort=sort(diags);
diagsSort=-diagsSort; % Sorted sojurn intensities

Cox=diag(-diagsSort); % Matrix of coxian represenation

Qk=zeros(1,n-1); % Equation 9 in Hjelmgren
P=zeros(1,n-1); % Equation 8 in Hjelmgren

% Computes the elements in coxian matrix
for i=1:n-1
    P(i)=ttc_absorbing_jump_prob(i,Q,startProb);
    Qk(i)=P(i)/(1-sum(P(1:i-1)));
    Cox(i,i+1)=-Cox(i,i)*(1-Qk(i));
    Cox(i,end)=-Cox(i,i)*Qk(i);
    
end
Cox(n,end)=-Cox(n,n);
CoxLambda=Cox(1:end-1,1:end-1);
CoxTheta=Cox(1:end-1,end);
