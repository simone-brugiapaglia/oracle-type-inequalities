function [pi L]=pi_opt(n,d,options,level)
%% Computes the drawing probability minimizing the coherence, i.e.
% pi_k proportional to || a_k ||_inf^2 with a_k the measurements vectors
% used with the Fourier-Wavelet transform.
%
% INPUT :
%   n = size of the square, cube ...
%   d = space dimension
%   options, level (optional) = parameters for the wavelets
%
% OUTPUT :
%   pi = output probability of size n x ... x n (d times)
%   L  = renormalizing factor used to scale the probability pi
%
% Authors : Pierre Weiss & Nicolas Chauffert
%

if nargin==2
    options.h= compute_wavelet_filter('Symmlet',10);
    options.ti=1;
    level=3;
end
Pi=zeros(n,1);
vect=zeros(n,1);
J=log2(n);
Jmin=max(1,J-level);
for i=1:n
    vect(i)=1;
    FW=perform_wavortho_transf(ifft(vect),Jmin,1,options);
    FW=FW/sqrt(sum(abs(FW).^2));
    vect(i)=0;
    Pi(i)=max(abs(FW))^2;
end
Pi=fftshift(Pi);
if d == 1
    pi = Pi/sum(Pi);
else
    D=ones(n*ones(1,d));
    for i=1:d
        dim=ones(1,d);
        dim2=n*ones(1,d);
        dim2(i)=1;
        dim(i)=n;
        D=D.*repmat(reshape(Pi,dim),dim2);
    end
    L=sum(reshape(D,n^d,1));
    pi=D/L;
end
end
