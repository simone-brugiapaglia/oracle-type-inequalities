%% test 1D
% MRI setting : acquisition of Fourier samples and reconstruction of
% wavelet coefficients with a sparsity-by-level structure
% Different drawing probabilities to construct the sensing matrix are compared.
%
% Dependencies :
% - the MATLAB numerical-tour of Gabriel Peyr\'e is needed.
% - export_fig is also necessary to export nice figures.
%
% Author : Claire Boyer

clear all
close all

rng(1) ;

%addpath(genpath('~/Documents/MATLAB/numerical-tour')) ;
addpath tools ; 

is_save_on = 1; % boolean = 1 to save pictures and .mat
if is_save_on 
    addpath(genpath('../export_fig')) ;
end


n = 2^11 ;



level=3;
J=log2(n);
Jmin=max(1,J-level);

options.h= compute_wavelet_filter('Symmlet',10);
options.ti=1;
d = 1;   % dimension 1D

L = J-level ;

%% Optimal drawing probability propotional to the squared inf-norm of measurements vectors

pi = pi_opt(n,d,options,level);
f1 = figure(1);
plot(pi,'LineWidth',2);
set(gca,'fontsize',16)
axis tight
% title({'$ \pi_k \propto \|a_k \|_\infty^2 $'},'interpreter','latex','FontSize',20) ;
% tit = strcat('optimal_proba_inf_1D_',num2str(n),'.pdf') ;
% print(f1,'-dpdf', tit) ;

if is_save_on 
    export_fig fig/proba_coherence_1D_2048 -pdf -m3
end


%% Structured sparsity

aux = Jmin:J ;
ilevel = 2.^(aux) ;
size_levels = [2^Jmin , ilevel(2:end)-ilevel(1:end-1)] ;
prop = [0.1 0.1 0.075 0.05] ;
struct_spars_1D = round(prop.*size_levels) ;
f2 = figure(2) ;
f2 = bar(1:length(aux),struct_spars_1D/2048,'LineWidth',2) ;
set(gca,'FontSize',18)

if is_save_on 
    %export_fig fig/struct_sparsity_1D -pdf -m3
end

%% Drawing probabilities adapted to structured sparsity

[pi2, L2]=pi_opt_lambda(n,d,options,level,struct_spars_1D) ;
f3 = figure(3);
plot(pi2,'LineWidth',2);
set(gca,'fontsize',16)
axis tight
%title({'$ \pi_k \propto \sum_{\ell} s_\ell \|a_{k,\Omega_\ell} \|_\infty^2  $'},'interpreter','latex','FontSize',20) ;
% tit = strcat('optimal_proba_lambda_1D_',num2str(n),'.pdf') ;
% print(f3,'-dpdf', tit) ;
if is_save_on 
    %export_fig fig/proba_lambda_1D_2048 -pdf -m3
end

[pi3, L3]=pi_opt_theta(n,d,options,level,struct_spars_1D) ;
f4 = figure(4);
plot(pi3,'LineWidth',2);
set(gca,'fontsize',16)
axis tight
%title({'$ \pi_k \propto \|a_{k} \|_\infty \sum_{\ell} s_\ell \|a_{k,\Omega_\ell} \|_\infty  $'},'interpreter','latex','FontSize',20) ;
% tit = strcat('optimal_proba_theta_1D_',num2str(n),'.pdf') ;
% print(f4,'-dpdf', tit) ;
if is_save_on 
    %export_fig fig/proba_theta_1D_2048 -pdf -m3
end 



%% Reconstruction of 1D signals
sampling_prop = 0.25 ;
ntest = 100;
SNR_unif_rec = zeros(ntest,1) ;
SNR_theta_rec = zeros(ntest,1) ;
SNR_lambda_rec = zeros(ntest,1) ;

for t=1:ntest
    % Signal generation with structured sparsity
    aux = Jmin:J ;
    ilevel = 2.^(aux) ;
    aux=[0, aux] ;
    ilevel= [0, ilevel];
    x = zeros(n,1) ;
    for j=1:length(ilevel)-1
        ind = ilevel(j)+randperm(ilevel(j+1)-ilevel(j),struct_spars_1D(j)) ;
        x(ind,1) = 10*rand(struct_spars_1D(j),1)/log(j+1) ;
    end

    % Measurements w/ unif sampling + rec
    samples = Draw_iid_Points_i(pi,round(sampling_prop*n));
    scheme_unif = zeros(n,1) ;
    scheme_unif(samples)=1 ;
    data_unif = @(u)(fftshift(scheme_unif).*(fft(perform_wavortho_transf(u,Jmin,-1,options)))/sqrt(n));
    y_unif = data_unif(x); % measurements

    x_unif_rec = MinL1_1D(scheme_unif,y_unif,options,level) ;
    SNR_unif_rec(t) = SNR_Rescale(x_unif_rec,x) ;

    % Measurements according to theta sampling + rec
    samples = Draw_iid_Points_i(pi3,round(sampling_prop*n));
    scheme_theta = zeros(n,1) ;
    scheme_theta(samples)=1 ;
    data_theta = @(u)(fftshift(scheme_theta).*(fft(perform_wavortho_transf(u,Jmin,-1,options)))/sqrt(n));
    y_theta = data_theta(x); % measurements

    x_theta_rec = MinL1_1D(scheme_theta,y_theta,options,level) ;
    SNR_theta_rec(t) = SNR_Rescale(x_theta_rec,x) ;

    % Measurements according to lambda sampling + rec
    samples = Draw_iid_Points_i(pi2,round(sampling_prop*n));
    scheme_lambda = zeros(n,1) ;
    scheme_lambda(samples)=1 ;
    data_lambda = @(u)(fftshift(scheme_lambda).*(fft(perform_wavortho_transf(u,Jmin,-1,options)))/sqrt(n));
    y_lambda = data_lambda(x); % measurements

    x_lambda_rec = MinL1_1D(scheme_lambda,y_lambda,options,level) ;
    SNR_lambda_rec(t) = SNR_Rescale(x_lambda_rec,x) ;
    if mod(t,10)==0
        fprintf('%d - ',t)
        %pause
    end
end

if is_save_on
    save data/SNR_unif_rec SNR_unif_rec
    save data/SNR_theta_rec SNR_theta_rec
    save data/SNR_lambda_rec SNR_lambda_rec
end


%% Boxplot display

load_data = 0 ; %boolean to load data or not
if load_data
    load data/SNR_unif_rec SNR_unif_rec
    load data/SNR_theta_rec SNR_theta_rec
    load data/SNR_lambda_rec SNR_lambda_rec
end

f10 = figure(10) ;
h = {SNR_unif_rec ; SNR_theta_rec ; SNR_lambda_rec}; % Create a cell array with the data for each group
aboxplot_mine(h,'labels',[],'colorgrad','mine'); % Advanced box plot
l10 = legend('$\pi^\infty$','$\pi^\Theta$','$\pi^\Lambda $') ;
set(l10,'FontSize',20,'Interpreter','Latex','Location','southeast');
tit = strcat('boxplot_1D_',num2str(n),'_',num2str(sampling_prop),'p_pos_sign.pdf') ;
if is_save_on 
    export_fig fig/boxplot_1D -pdf -m3
end 