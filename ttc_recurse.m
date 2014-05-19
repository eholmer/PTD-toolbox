function y = ttc_recurse(start, result, depth, v, lambda, L, prev)
% Computes the probability that the absorbing jump occurs from state L
% given that the path is v.

% start must always start with value 1
% result must always start with value 0
% depth must always start with length(v)-1
% v must be a vector of increasing integers
% lambda is a vector of sorted intensities high to low
% l is the state where the absorbing jump will occur from
% prev should always start with 0

k=length(v);

if depth==0
    result=ttc_PI(v(end),L,prev,lambda);
else
    for i=start:min(v(k-depth),L-k+(k-depth))
        result = result+ttc_PI(v(k-depth),i,prev,lambda)*ttc_recurse(i+1,0,depth-1,v,lambda,L,i);
    end
end
y=result;
