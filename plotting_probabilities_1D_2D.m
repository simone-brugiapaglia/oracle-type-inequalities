%% Script for generating Figure 2
% Plotting drawing probabilities in 2D for the MRI setting
% (acquisition of Fourier samples and reconstruction of
% wavelet coeff)

clear all
close all


%addpath(genpath('~/Documents/MATLAB')) ;

addpath tools;

n = 2^7 ;


level=3;
J=log2(n);
Jmin=max(1,J-level);

options.h= compute_wavelet_filter('Symmlet',10);
options.ti=1;
d = 1;   % dimension 1D

L = J-level ;

%% Optimal drawing probability, "uniform" for s-sparse vectors

pi = pi_opt(n,d,options,level);
f1 = figure(1);
plot(pi,'LineWidth',2);
set(gca,'fontsize',16)
axis tight
title({'$ \pi_k \propto \|a_k \|_\infty^2 $'},'interpreter','latex','FontSize',20) ;
% tit = strcat('optimal_proba_inf_1D_',num2str(n),'.pdf') ;
% print(f1,'-dpdf', tit) ;

export_fig fig/proba_coherence_1D_2048 -pdf -m3


f11 = figure(11);
Z = pi*pi' ;
Z = Z/sum(Z(:)) ;
surf(1:length(pi),1:length(pi),Z,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.5)
%surf(peaks,'EdgeColor','none','FaceColor','interp','FaceLighting','phong') 
%daspect([1,1,.3]);axis tight; 

export_fig fig/proba_coherence_3D -pdf -m3


%% Defining structured sparsity - defined by tensor product of 1D structured sparsity


% % Structured sparsity - choice 1
aux = Jmin:J ;
nb_level_1D =  length(Jmin:J);
ilevel = 2.^(aux) ;
size_levels = [2^Jmin , ilevel(2:end)-ilevel(1:end-1)] ;
prop = [0.1 0.1 0.075 0.05] ;
struct_spars_1D = round(prop.*size_levels) ;


% % Structured sparsity - choice 2
% nb_level_1D =  length(Jmin:J);
% global_spars_1D = round(0.06*n);
% lev = 1:nb_level_1D ;
% coef = (1/2).^(lev-1) ;
% s1 = round(global_spars_1D*2^nb_level_1D/(2*(2^nb_level_1D-1))) ;
% struct_spars_1D = round(s1*coef) ;



% % Plotting the sparsity structure accross the levels
% f5 = figure(5) ;
% stem(1:nb_level_1D,struct_spars_1D,'color',[0.48 0.41 0.93], 'LineWidth',2) ;
% set(gca,'fontsize',16)
% axis tight
% title('Structured sparsity : sparsity by level','Interpreter','latex') ;
% % tit = 'struct_sparsity_1D.pdf' ;
% % print(f2,'-dpdf', tit) ;


%% Sampling for structured sparsity optimizing the lambda quantity

[pi2, L2]=pi_opt_lambda(n,d,options,level,struct_spars_1D) ;
f2 = figure(2);
plot(pi2,'LineWidth',2);
set(gca,'fontsize',16)
axis tight
title({'$ \pi_k \propto \sum_{\ell} s_\ell \|a_{k,\Omega_\ell} \|_\infty^2  $'},'interpreter','latex','FontSize',20) ;
% tit = strcat('optimal_proba_lambda_1D_',num2str(n),'.pdf') ;
% print(f3,'-dpdf', tit) ;

export_fig fig/proba_lambda_1D_2048 -pdf -m3

f21 = figure(21);
Z2 = pi2*pi2' ;
Z2 = Z2/sum(Z2(:)) ;
surf(1:length(pi2),1:length(pi2),Z2,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.5)

export_fig fig/proba_lambda_3D -pdf -m3


%% Sampling for structured sparsity optimizing the theta quantity

[pi3, L3]=pi_opt_theta(n,d,options,level,struct_spars_1D) ;
f3 = figure(3);
plot(pi3,'LineWidth',2);
set(gca,'fontsize',16)
axis tight
title({'$ \pi_k \propto \|a_{k} \|_\infty \sum_{\ell} s_\ell \|a_{k,\Omega_\ell} \|_\infty  $'},'interpreter','latex','FontSize',20) ;
% tit = strcat('optimal_proba_theta_1D_',num2str(n),'.pdf') ;
% print(f4,'-dpdf', tit) ;


export_fig fig/proba_theta_1D_2048 -pdf -m3

f31 = figure(31);
Z3 = pi3*pi3' ;
Z3 = Z3/sum(Z3(:)) ;
surf(1:length(pi3),1:length(pi3),Z3,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.5)

export_fig fig/proba_theta_3D -pdf -m3



