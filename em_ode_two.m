function [Lambda, theta, Pi, mu] = em_ode_two(data1,data2, Lambda, theta, Pi, mu, n_iterations)
% The EM-algorithm for a phase-type distribution using two samples assumed
% to be different only by a factor (mu) applied to the absorption intensities.
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

n_samples1 = length(data1);
n_samples2 = length(data2)

% Dimension of the phase-type distribution 
p = size(Lambda,1);

% Unit matrix used to pick out columns of matrices
unit = eye(p);

% Set lower tolerance to increase speed
options = odeset('NonNegative',1:(p*p+2*p));

% Sorted in order to help ode45 calculate all points simultaneously.
data1 = sort(data1);
data2 = sort(data2);

% EM-algorithm
for k=1:n_iterations
    %Update 
    
    % Update everything but `mu`
    if mod(k,2) == 1  
        B = zeros(1,p);
        BO = zeros(1,p);
        Z = zeros(1,p);
        ZO = zeros(1,p);
        N = zeros(p,p);
        NO = zeros(p,p);
        T = zeros(1,p);
        TO = zeros(1,p);
        [~,U] = ode45(@(t,y)ph_ode_generator(Lambda,theta',p,y),[0 data1],[Pi';theta';zeros(p*p,1)],options);
        [~,UO] = ode45(@(t,y)ph_ode_generator(Lambda,mu*theta',p,y),[0 data2],[Pi';mu*theta';zeros(p*p,1)],options);

        % E-step
        for v=1:n_samples1        
            aa = U(v+1,1:p);
            ba = U(v+1,(p+1):2*p);
            ca = reshape(U(v+1,(2*p+1):end),p,p)';
            d = dot(Pi, ba);            

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
        for v=1:n_samples2        
            aaO = UO(v+1,1:p);
            baO = UO(v+1,(p+1):2*p);
            caO = reshape(UO(v+1,(2*p+1):end),p,p)';
            dO = dot(Pi, baO);            

            for i=1:p
                BO(i) = BO(i) + (Pi(i) * baO(i) / dO);
                ZO(i) = ZO(i) + caO(i,:) * unit(:,i) / dO;

                % theta update
                TO(i) = TO(i) + mu*theta(i) * aaO(i) / dO;

                for j=1:p
                    if ~(i==j)
                        NO(i,j) = NO(i,j) + Lambda(i,j) * caO(i,:) * unit(:,j) / dO;
                    end
                end
            end
        end
        
        % M-step
        Pi = (B + BO) / (n_samples1+n_samples2);
        for i=1:p
            Lambda(i,:) = (N(i,:) + NO(i,:)) / (Z(i) + ZO(i));
        end
        theta = (T + TO) ./ (Z + mu*ZO);
        for i=1:p
            Lambda(i,i) = -(theta(i) + sum(Lambda(i,:)));
        end

    % Update `mu`
    else
        ZO = zeros(1,p);
        TO = zeros(1,p);
        [~,UO] = ode45(@(t,y)ph_ode_generator(Lambda,mu*theta',p,y),[0 data2],[Pi';mu*theta';zeros(p*p,1)],options);

        % E-step
        for v=1:n_samples2 
            aaO = UO(v+1,1:p);
            baO = UO(v+1,(p+1):2*p);
            caO = reshape(UO(v+1,(2*p+1):end),p,p)';
            dO = dot(Pi, baO);

            for i=1:p
                ZO(i) = ZO(i) + caO(i,:) * unit(:,i) / dO;

                % theta update
                TO(i) = TO(i) + mu*theta(i) * aaO(i) / dO;

            end
        end
        % M-step        
        mu = sum(TO) / sum(theta.*ZO);
    end
end
