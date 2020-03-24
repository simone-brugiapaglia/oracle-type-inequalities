function [pi, L]=pi_opt_theta(n,d,options,level,sparsities)
%% Computes the drawing probability minimizing the theta quantity, i.e.
% pi_k proportional to a uniform bound on || B_k B_{k,S}^* ||_{inf -> inf}
% with B_k the measurements block corresponding to a 2D horizontal line
% in the Fourier-Wavelet transform,
% and with the support S structured by block given by 'sparsities'.
%
% INPUT :
%   n = size of the square, cube ...
%   d = space dimension
%   options, level (optional) = parameters for the wavelets
%   sparsities = vector of sparsity degrees in wavelet levels 
%
% OUTPUT :
%   pi = output probability of size n x ... x n (d times)
%   L  = renormalizing factor used to scale the probability pi
%
% Author : Claire Boyer
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
aux = Jmin:J ;
ilevel = 2.^(aux) ;
ilevel= [0, ilevel];

for i=1:n
    vect(i)=1;
    FW=perform_wavortho_transf(ifft(vect),Jmin,1,options);
    FW=FW/sqrt(sum(abs(FW).^2));
    norm_FW = zeros(length(ilevel)-1,1) ;
    for j=1:length(ilevel)-1
        norm_FW(j) = max(abs(FW(ilevel(j)+1:ilevel(j+1)))) ;
    end
    vect(i)=0;
    Pi(i)=sum(norm_FW.*sparsities')*max(norm_FW);
end
Pi=fftshift(Pi);

if d==1
    L = sum(Pi) ;
    pi=Pi/L ;

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
