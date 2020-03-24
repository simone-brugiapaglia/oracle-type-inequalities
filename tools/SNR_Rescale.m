function snr=SNR_Rescale(im,im_ref)
%% function snr=SNR_Rescale(im,im_ref)
%
% Computes SNR by first finding best linear scaling of im
% Author : Pierre Weiss

n=numel(im);
M=[[sum(im(:)),-n];[norm(im(:))^2,-sum(im(:))]];

coeffs=M\[sum(im_ref(:));sum(im_ref(:).*im(:))];
im=im*coeffs(1)-coeffs(2);

Psignal=norm(im_ref(:))^2;
Pnoise=norm(im(:)-im_ref(:))^2;
snr=10*log10(Psignal/Pnoise);
