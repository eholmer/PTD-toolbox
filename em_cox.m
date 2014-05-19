function [Lambda, theta, Pi] = em_cox(data, Lambda, theta, Pi, n_iterations)
% The EM-algorithm for a coxian phase-type distribution.
%
% Approximate the maximum likelihood estimates for a general phase-type
% distribution via an iterative algorithm.
%
% This algorithm follows the development in the paper
% "Fitting phase-type distributions via the EM-algorithm"
% by M. Olsson et. al.
%
% Parameters:
%
% data: Vector of data points (p x 1)
% Lambda: Generator matrix (p x p)
% Pi: Initial distribution (1 x p)
% theta: Absorbtion intensities (p x 1)
% n_iterations: Number of iterations in the EM-algorithm

n_samples = length(data);

% Dimension of the phase-type distribution 
p = size(Lambda,1);

% Unit matrix used to pick out columns of matrices
unit = eye(p);

% Disallow negative numbers
options = odeset('NonNegative',1:(p*p+2*p));

% Sorted in order to help ode45 calculate all points simultaneously.
data = sort(data);

% EM-algorithm
Lambda = sparse(Lambda);

for k=1:n_iterations
    Z = zeros(1,p);
    N = zeros(1,p-1);
    super_diag = zeros(1,p-1);
    main_diag = zeros(1,p);
    T = zeros(1,p);
    
    [~,U] = ode45(@(t,y)ph_ode_generator(Lambda,theta',p,y),[0 data],[Pi';theta';zeros(p*p,1)],options);

    % E-step
    for v=1:n_samples
        aa = U(v+1,1:p);
        ba = U(v+1,(p+1):2*p);
        ca = reshape(U(v+1,(2*p+1):end),p,p)';
        bidiag = spdiags(Lambda,1);
        d = dot(Pi, ba);
        
        Z = Z + diag(ca)' / d;
        T = T + (theta .* aa) / d;
        N = N + (bidiag(2:end) .* diag(ca(:,2:end)))' / d;
    end
    
    % M-step
    theta = T ./ Z;
    super_diag = N ./Z(1:(end-1));
    main_diag = [-(theta(1:(end-1))+super_diag) -theta(end)];
    
    
    Lambda = spdiags([main_diag' [0 super_diag]'], [0 1], Lambda);

end
Lambda = full(Lambda);
