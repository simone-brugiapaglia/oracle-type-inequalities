%% Test 2: Convergence of the measure via adptive sampling
%
% The goal of this code is to understand if adaptive sampling strategies 
% are able to converge to a Chebyshev-like distribution.
% We consider n Legendre polynomials on (-1,1) and an equispaced grid of 
% n^2 points on (-1,1). The corresponding sampling matrix is approximately
% an isometry. Then, we study the evolution of the sampling measure as a
% function of the iterative step in the outer loop of the adaptive
% procedure.
%
% Dependencies:
% - SPGL1 [https://github.com/mpf/spgl1]

% Author: Simone Brugiapaglia
% Concordia University
% simone.brugiapaglia@concordia.ca

clear all; close all

addpath tools

%% Parameters
n = 100;        % polynomial space dimension
s_signal = 5; % sparsity of the signal
m1 = 10;       % number of adaptive measurements added at each iteration
K = 5;         % number of adaptive iterations (outer loop)


s_adapt_1 = s_signal; % sparsity parameter for (Adapt I)
s_adapt_2 = s_signal; % sparsity parameter for (Adapt II)
n_grid = n^2;  % number of sampling points in the uniform grid 
               % the scaling n^2 guarantees that ||Ax||_2 ~ ||x||_2



%% Approximate isometry
Unif_grid = linspace(-1,1,n_grid+2);
Unif_grid = Unif_grid(2:end-1);
A0 = (n/n_grid) * LegMat(Unif_grid,n);


fprintf('Well-conditioning check: cond(A0) = %1.1e\n',cond(A0))
fprintf('Isometry check: norm(A0''*A0-I) = %1.1e\n',norm(A0'*A0 - eye(n)))



%% Generate sparse random vector (i.e., sparse combination of Legendre polynomials)
x_exact = zeros(n,1);
supp = randperm(n);
supp = supp(1:s_signal);
x_exact(supp) = randn(s_signal,1);



%% (Adapt I) Adaptive sampling Variant I - fixed sparsity
J_1  = zeros(m1,K);  % k-th group of samples in k-th column
xh_1 = zeros(n,K);   % k-th approximation in k-th column
S_1  = zeros(s_adapt_1,K-1); % k-th support in k-th column
pi_1 = zeros(n_grid,K);   % k-th adapted measure in k-th column

% first random sampling (from uniform measure)
J_1(:,1) = randi(n_grid,m1,1);
pi_1(:,1) = 1/n_grid;
A = sqrt(n_grid/m1) * A0(J_1(:,1),:);
y = A * x_exact; % noiseless measurements

opts = spgSetParms('verbosity',0);
xh_1(:,1) = spg_bp(A,y,opts);


% itrations 2,...,K with adaptive sampling
for k = 2:K
    fprintf('(Adapt I) Iteration %d out of %d\n',k,K);
    
    % Update support information
    [xh_sort, I_sort] = sort(abs(xh_1(:,k-1)),1,'descend');
    S_1(:,k-1) = I_sort(1:s_adapt_1);
    
    % Define sampling measure based on the support information
    for j = 1:n_grid
        pi_1(j,k) = norm(A0(j,S_1(:,k-1)),2)^2;
    end
    pi_1(:,k) = pi_1(:,k)/sum(pi_1(:,k));
    cum_dist = cumsum(pi_1(:,k)); % precomputation to make nurandi faster
    
    % Draw random samples using nurandi routine
    J_1(:,k) = nurandi(n_grid, m1, pi_1(:,k), cum_dist);
    
    % Compute current samples and sampling matrix
    Ck = diag(1./sqrt(m1 * pi_1(J_1(:,k),k))); % diagonal preconditioner (due to nonuniform sampling)
    Ak = Ck * A0(J_1(:,k),:);
    yk = Ak * x_exact;
    
    % Stack previous and current samples
    A = 1/sqrt(k)*[sqrt(k-1)*A; Ak];
    y = 1/sqrt(k)*[sqrt(k-1)*y; yk];
    
    % l1 minimization
    xh_1(:,k) = spg_bp(A,y,opts);
end

x_adapt_1 = xh_1(:,K);
err_adapt_1 = norm(x_exact-x_adapt_1);


%% (Adapt II) Adaptive sampling Variant II - increasing sparsity
J_2  = zeros(m1,K);  % k-th group of samples in k-th column
xh_2 = zeros(n,K);   % k-th approximation in k-th column
S_2  = zeros((K-1)*s_adapt_2,K-1); % k-th support in k-th column
pi_2 = zeros(n_grid,K);   % k-th adapted measure in k-th column

% first random sampling (uniform measure)
J_2(:,1) = randi(n,m1,1);
pi_2(:,1) = 1/n_grid;
A = sqrt(n_grid/m1) * A0(J_2(:,1),:);
y = A * x_exact;

opts = spgSetParms('verbosity',0);
xh_2(:,1) = spg_bp(A,y,opts);


% itrations 2,...,K with adaptive sampling
for k = 2:K
    fprintf('(Adapt II) Iteration %d out of %d\n',k,K);
    
    % Update support information
    [xh_sort, I_sort] = sort(abs(xh_2(:,k-1)),1,'descend');
    S_2(1:(k-1)*s_adapt_2,k-1) = I_sort(1:(k-1)*s_adapt_2); % SUPPORT UPDATE (different from Variant I)
    
    % Define sampling measure based on the support information
    for j = 1:n_grid
        pi_2(j,k) = norm(A0(j,S_2(1:(k-1)*s_adapt_2,k-1)),2)^2;
    end
    pi_2(:,k) = pi_2(:,k)/sum(pi_2(:,k));
    cum_dist = cumsum(pi_2(:,k)); % precomputation to make nurandi faster
    
    % Draw random samples
    J_2(:,k) = nurandi(n_grid,m1,pi_2(:,k),cum_dist);
    
    % Compute current samples and sampling matrix
    Ck = diag(1./sqrt(m1 * pi_2(J_2(:,k),k)));
    Ak = Ck * A0(J_2(:,k),:);
    yk = Ak * x_exact;
    
    % Stack previous and current samples
    A = 1/sqrt(k)*[sqrt(k-1)*A; Ak];
    y = 1/sqrt(k)*[sqrt(k-1)*y; yk];
    
    % l1 minimization
    xh_2(:,k) = spg_bp(A,y,opts);
end

x_adapt_2 = xh_2(:,K);
err_adapt_2 = norm(x_exact-x_adapt_2);


% Compute distance from Chebyshev measure
pi_Cheby = 1./sqrt(1-Unif_grid.^2) / sum(1./sqrt(1-Unif_grid.^2));
for k = 1:K
    err_meas_1(k) = norm(pi_1(:,k) - pi_Cheby(:),1);
    err_meas_2(k) = norm(pi_2(:,k) - pi_Cheby(:),1);
end


clear A0;

save('data/fun_approx_test_2')



%% Visualize results

clear all; close all;

load data/fun_approx_test_2.mat

%% Printf data
% Compute distance from Chebyshev measure
pi_Cheby = 1./sqrt(1-Unif_grid.^2) / sum(1./sqrt(1-Unif_grid.^2));
for k = 1:K
    err_meas_1(k) = norm(pi_1(:,k) - pi_Cheby(:),2);
    err_meas_2(k) = norm(pi_2(:,k) - pi_Cheby(:),2);
end

%err_adapt_1
%err_adapt_2




for k = 1:K
    figure;
    semilogy(Unif_grid,pi_1(:,k),'r-',Unif_grid, pi_Cheby,'--k')
    
    ht = title(['$\pi^{(',num2str(k),')}$']);
    hl = legend(['$\pi^{(',num2str(k),')}$'],'Chebyshev');
    set(hl,'interpreter','latex');
    set(ht,'interpreter','latex');
    set(gca,'TickLabelInterpreter', 'latex');
    ylim([min(min(min(pi_1)),min(pi_Cheby)),max(max(max(pi_1)),max(pi_Cheby))])
    
    pbaspect([2.5 1 1])
    
    grid on
    set(gca,'fontsize',20)
    saveas(gca,['fig/Conv_meas_Adapt_I_k_',num2str(k)],'epsc')
end



for k = 1:K
    figure;
    semilogy(Unif_grid,pi_2(:,k),'r-',Unif_grid, pi_Cheby,'--k')
    
    ht = title(['$\pi^{(',num2str(k),')}$']);
    hl = legend(['$\pi^{(',num2str(k),')}$'],'Chebyshev');
    set(hl,'interpreter','latex');
    set(ht,'interpreter','latex');
    set(gca,'TickLabelInterpreter', 'latex');
    ylim([min(min(min(pi_2)),min(pi_Cheby)),max(max(max(pi_2)),max(pi_Cheby))])
    
    pbaspect([2.5 1 1])
    
    grid on
    set(gca,'fontsize',20)
    saveas(gca,['fig/Conv_meas_Adapt_II_k_',num2str(k)],'epsc')
end

disp(' ')
disp('Table 1 data:')

[err_meas_1',err_meas_2']