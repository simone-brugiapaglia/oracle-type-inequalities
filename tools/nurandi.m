function X = nurandi(IMAX,N,p,cdist)
% X = NURANDI(IMAX,N,P) returns a 1-by-N array containing pseudorandom
% integer values drawn from the discrete distribution P / sum(P) on 1:IMAX.
% In particular, 
%
%       Prob{ X(i) = k } = P(k) / sum(P),      i = 1,...,M,   k = 1,...N. 
%
% X = NURANDI(IMAX,N,P,CDIST) returns the same output as NURANDI(IMAX,N,P), 
% but the cumulative distribution CDIST is given as input to accelerate the 
% computation.

% Simone Brugiapaglia, SFU
% simone_brugiapaglia@sfu.ca
% Last Update: December 2017

%% Preliminary checks
% p non-negative vector
if ~ (p > 0) 
    error('p has to be non-negative')
end

if N ~= round(N)
    error('M must be an integer');
elseif IMAX ~= round(IMAX)
    error('N must be an integer');
end

% initialization
URN = 1:IMAX;

p = p(:);

if sum(p) ~=1
    p = p/(sum(p));
end

% build cumulative distribution function (if not already given in input)
if nargin == 3  
    cdist = cumsum(p);
end
if cdist(end) ~= 1
    cdist = cdist/cdist(end);
end
cdist = cdist(:);

% drawing strategy: draw a pseudorandom number Y ~ Unif(0,1) and then find 
% the smallest k such that F(h-1) < Y <= F(h), where F is the cumulative 
% distribution of P.
indices = zeros(1,N);
for i = 1:N
    Y = rand;
    indices(i) = find(Y < cdist,1,'first'); % OLD: IMAX - sum(rnd_num <= cdist) + 1;
end

X = URN(indices);



