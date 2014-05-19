function p = ttc_absorbing_jump_prob(k, Q, Pi)
% computes the probability that an absorbing jump occurs from state k
% k is the state where the event jump will occur
% Q is the transition matrix
% Pi is the start probability vector
p=0;

M=Q(1:end-1,1:end-1); % Generator matrix
diags=diag(M); % Diagonal sojurn intensities

% sigma is a function from index of sorted intesities to index in tringular matrix
[diagsSort, sigma]=sort(diags);
diagsSort=-diagsSort;

% sigmaInv is a function from index in tri to index in sorted
[dummy, sigmaInv]=sort(sigma);

% Indicies of all non-absorbing states
states=linspace(1,length(M),length(M));

for i=1:length(M)
    V=combnk(states,i); % Every cominations of i states
    for j=1:length(V(:,1)) % Loops every possible path with i states
        u=V(j,:); % Path in ES indexes
        v=sort(sigmaInv(u))'; % Path in sorted indexes
        if (k <= v(end) && k>=length(v)) % Absorbing jump cannot happen later than
            % the last state and not before the number of states on path
            % Adds the probability of taking a path and the absorbing jump
            % Occuring from k
            prob=ttc_path_prob(u,Q,Pi)*ttc_recurse(1,0,length(v)-1,v,diagsSort,k,0);
            p=p+prob; % ttc_path_prob(u,Q,Pi)*ttc_recurse(1,0,length(v)-1,v,diagsSort,k,0);
        end
    end
end
