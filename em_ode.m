function [Lambda, theta, Pi] = em_ode(data, Lambda, theta, Pi, n_iterations)
% The EM-algorithm for a phase-type distribution.
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
for k=1:n_iterations
    B = zeros(1,p);
    Z = zeros(1,p);
    N = zeros(p,p);
    T = zeros(1,p);
    [~,U] = ode45(@(t,y)ph_ode_generator(Lambda,theta',p,y),[0 data],[Pi';theta';zeros(p*p,1)],options);

    % E-step
    for v=1:n_samples        
        aa = U(v+1,1:p);
        ba = U(v+1,(p+1):2*p);
        ca = reshape(U(v+1,(2*p+1):end),p,p)';
        d = dot(Pi, ba);
        
        assert(sum(aa < 0) == 0)
        assert(sum(ba < 0) == 0)
        assert(sum(ca(:) < 0) == 0)
        assert(sum(d < 0) == 0)

        for i=1:p
            B(i) = B(i) + (Pi(i) * ba(i) / d);
            Z(i) = Z(i) + ca(i,:) * unit(:,i) / d;
            
            % theta update
            T(i) = T(i) + theta(i) * aa(i) / d;

            for j=1:p
                if ~(i==j)
                    N(i,j) = N(i,j) + Lambda(i,j) * ca(i,:) * unit(:,j) / d;
                end
            end

        end
    end
    
    % M-step
    Pi = B / n_samples;
    for i=1:p
        Lambda(i,:) = N(i,:) / Z(i);
    end
    theta = T ./ Z;
    for i=1:p
        Lambda(i,i) = -(theta(i) + sum(Lambda(i,:)));
    end
end
