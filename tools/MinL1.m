function rec = MinL1(scheme,measure,options,level)
%% function rec = MinL1(scheme,measure,options,level)
%
% Solve the l1-problem:
%
%    min ||x||_1 s.t Ax=y
%
% with Douglas-Rachford algorithm,
% when A is a matrix extracted from a 2D Fourier-Wavelet transform.
%
% INPUT :
%   - scheme : sampling mask of size nxn
%   - measure : measurement vector
%   - option, level : parameters for the wavelet transform
% OUTPUT :
%   - rec : reconstructed image 


n=size(scheme,1);
sigma2=circshift(scheme,[n/2 n/2]); % Put low-frequencies in corners
J=log2(n);
J_min=max(1,J-level);


multA = @(x) sigma2.*fft2(perform_wavortho_transf(x,J_min,-1,options))/n;
Prox_N1 = @(x,tau)max(0,1-tau./max(1e-15,abs(x))).*x;
Prox_set = @(x)x+perform_wavortho_transf(n*ifft2(measure-multA(x)),J_min,1,options);

% Parameters of Douglas Rachford
lambda0=1.5;
lambda=lambda0;
eps=1e-5;
gamma=1;
y=0;
N=100;

for i=1:N
    x=Prox_set(y);
    y=y+lambda*(Prox_N1(2*x-y,gamma)-x);
  % err(i)=norm(x_s-x)/norm(x_s);
end

%figure;plot(err)
rec=x;
end
