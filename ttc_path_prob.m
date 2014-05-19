function p = ttc_path_prob(v, Q, Pi)
% Computes the probability of a certain path v
% v is a vector of transition states in a path
% Q is the transition matrix of intesities including absorbing state
% Pi is the start probability vector

p=Pi(v(1)); % Start probability ggggggg

for i=1:length(v)-1 % Transition jumps
    p=p*(Q(v(i),v(i+1))/(-Q(v(i),v(i))));
end
p=p*(Q(v(end),end)/(-Q(v(end),v(end)))); % Absorbing jump
