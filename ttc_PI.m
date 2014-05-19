function p = ttc_PI(ik, L, lk1, lambdas)
% ik: State in the elementary series we are jumping from
% L: State where we are jumping from
% lk1: State from which the previous event jump occured
% lambdas: vector of sorted sojourn intensities

p=1; % Probability for absorbing jump

% Absorbing jump cannot happen later than (ik) and not before the previous
% event jump
if(L>ik || L<=lk1)
    p=0;
else
    for i= lk1+1 : L-1 % Do not absorb before L
        p=p.*(1-lambdas(ik)./lambdas(i));
    end
    p=p.*(lambdas(ik)./lambdas(L)); % Absorb from L
end
