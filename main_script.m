%% Main script
%
% By running this Matlab script you will reproduce the plots in the paper 
%
%   "On oracle-type local recovery guarantees in compressed sensing" 
%   Information and Inference: A Journal of the IMA (Accepted, 2020)
%   by Ben Adcock, Claire Boyer, and Simone Brugiapaglia
%
% A preprint of the paper can be found at https://arxiv.org/abs/1806.03789
%
%
% DEPENDENCIES: To run this script, you need to add to your Matlab path
%
%  1. Numerical tours by Gabriel Peyré 
%     https://www.numerical-tours.com/matlab/
%  2. export_fig 
%     https://github.com/altmany/export_fig
%  3. SPGL1 
%     https://github.com/mpf/spgl1
%
%
% DISCLAIMER: The following experiments involve randomization. Hence, the 
% plots obtained by running this script might slightly differ from those 
% in the paper.
%
% Code by Claire Boyer and Simone Brugiapaglia

plotting_probabilities_1D_2D;  % Reproduces plots in Figure 2 
Haar_Fourier_test_1D;          % Reproduces plots in Figure 3
function_approximation_test_1; % Reproduces plots in Figure 6 and Figure 7
function_approximation_test_2; % Reproduces plots in Figure 8 and data in Table 1