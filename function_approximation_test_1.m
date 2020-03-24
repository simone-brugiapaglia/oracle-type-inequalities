%% Test 1: Adaptive sampling for 1D function approximation 
%
% Comparison of four stategies for function approximation via l1
% minimization
% (Adapt 1) Adaptive sampling based on a support update with fixed sparsity
% (Adapt 2) Adaptive sampling based on a support update with increasing sparsity
% (Unif I)  Uniform random sampling over [-1,1]
% (Unif II) Uniform random sampling from the discrete isometry A0
% (Cheby)   Random sampling from the Chebyshev measure with preconditioning
%           (as in Rauhut-Ward)
%
% These methods are compared to recover a sparse combination of Legendre
% polynomials
%
% Dependencies:
% - SPGL1 [https://github.com/mpf/spgl1]

% Author: Simone Brugiapaglia
% Concordia University
% simone.brugiapaglia@concordia.ca


clear all
close all

addpath tools

n  = 150;
s_adapt  = 5; % sparsity parameter for the adptive procedure
m1 = 15;
K  = 5;


s_signal = s_adapt; % we recover a (s_signal)-sparse combination of Legendre functions

N_test = 100;


m = K * m1;



% Build A0 and y0
[g,w]=lgwt(n,-1,1);
D = diag(sqrt(w/2));
A0 = D * LegMat(g,n);


fprintf('Isometry check: ||A0''*A0-I|| = %1.1e\n',norm(A0'*A0 - eye(n)))


% error matrix
Err = zeros(5,N_test);

for i_test = 1:N_test
    
    if mod(i_test,10)==0
        fprintf('Test %d of %d\n',i_test,N_test)
    end
    
    %% Generate (s_signal)-sparse vector of coefficients x_exact
    x_exact = zeros(n,1);
    supp = randperm(n);
    supp = supp(1:s_signal);
    x_exact(supp) = randn(s_signal,1);
    
    
    %% (Adapt I) Adaptive sampling Variant I - fixed sparsity
    J_1  = zeros(m1,K);  % k-th group of samples in k-th column
    xh_1 = zeros(n,K);   % k-th approximation in k-th column
    S_1  = zeros(s_adapt,K-1); % k-th support in k-th column
    pi_1 = zeros(n,K);   % k-th adapted measure in k-th column
    
    % first random sampling
    J_1(:,1) = randi(n,m1,1);
    pi_1(:,1) = 1/n;
    A = sqrt(n/m1) * A0(J_1(:,1),:);
    y = A * x_exact; % noiseless measurements
    
    opts = spgSetParms('verbosity',0);
    xh_1(:,1) = spg_bp(A,y,opts);
    
    
    % itrations 2,...,K with adaptive sampling
    for k = 2:K
        %fprintf('Iteration %d out of %d\n',k,K);
        
        % Update support information
        [xh_sort, I_sort] = sort(abs(xh_1(:,k-1)),1,'descend');
        S_1(:,k-1) = I_sort(1:s_adapt);
        
        % Define sampling measure based on the support information
        for j = 1:n
            pi_1(j,k) = norm(A0(j,S_1(:,k-1)),2)^2;
        end
        pi_1(:,k) = pi_1(:,k)/sum(pi_1(:,k));
        cum_dist = cumsum(pi_1(:,k)); % precomputation to make nurandi faster
        
        % Draw random samples using nurandi routine
        J_1(:,k) = nurandi(n,m1,pi_1(:,k),cum_dist);
        
        % Compute current samples and sampling matrix
        Ck = diag(1./sqrt(m1 * pi_1(J_1(:,k),k))); % diagonal preconditioner (nonuniform sampling)
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
    S_2  = zeros((K-1)*s_adapt,K-1); % k-th support in k-th column
    pi_2 = zeros(n,K);   % k-th adapted measure in k-th column
    
    % first random sampling
    J_2(:,1) = randi(n,m1,1);
    pi_2(:,1) = 1/n;
    A = sqrt(n/m1) * A0(J_2(:,1),:);
    y = A * x_exact;
    
    opts = spgSetParms('verbosity',0);
    xh_2(:,1) = spg_bp(A,y,opts);
    
    
    % itrations 2,...,K with adaptive sampling
    for k = 2:K
        %fprintf('Iteration %d out of %d\n',k,K);
        
        % Update support information
        [xh_sort, I_sort] = sort(abs(xh_2(:,k-1)),1,'descend');
        S_2(1:(k-1)*s_adapt,k-1) = I_sort(1:(k-1)*s_adapt); % SUPPORT UPDATE (different from Variant I)
        
        % Define sampling measure based on the support information
        for j = 1:n
            pi_2(j,k) = norm(A0(j,S_2(1:(k-1)*s_adapt,k-1)),2)^2;
        end
        pi_2(:,k) = pi_2(:,k)/sum(pi_2(:,k));
        cum_dist = cumsum(pi_2(:,k)); % precomputation to make nurandi faster
        
        % Draw random samples
        J_2(:,k) = nurandi(n,m1,pi_2(:,k),cum_dist);
        
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
    
    
    %% (Unif I) Uniform random sampling (continuous measure)
    pts_unif_cont = 2 * rand(m,1)-1;
    A = 1/sqrt(m) * LegMat(pts_unif_cont,n);
    y = A * x_exact;
    x_unif_cont = spg_bp(A,y,opts);
    err_unif_cont = norm(x_exact-x_unif_cont);
    
    
    %% (Unif II) Uniform random sampling (discrete measure)
    J_unif = randi(n,m,1);
    A = sqrt(n/m) * A0(J_unif,:);
    y = A * x_exact;
    x_unif_disc = spg_bp(A,y,opts);
    err_unif_disc = norm(x_exact-x_unif_disc);
    
    %% (Cheby) Random sampling from Chebyshev measure w/preconditioning (Rauhut-Ward)
    pts_Cheby = cos(pi * rand(m,1)); % samples from Chebyshev measure
    Preco = diag(sqrt(pi/2) * (1 - pts_Cheby.^2).^(1/4)); % diagonal preconditioner
    A = 1/sqrt(m) * Preco * LegMat(pts_Cheby,n);
    y = A * x_exact;
    x_Cheby = spg_bp(A,y,opts);
    err_Cheby = norm(x_exact-x_Cheby);
    
    %% store data
    Err(:,i_test) = [err_adapt_1; err_adapt_2; err_unif_cont; err_unif_disc; err_Cheby];
    
end

clear A0;

save('data/fun_approx_test_1')




%% Visualize results



clear all;
close all;


load data/fun_approx_test_1.mat

%% Visualize adapted measure (Variant I)
figure;
imagesc(pi_1)
xlabel('Iteration $k$','interpreter','latex')
ylabel('Index $j$','interpreter','latex')
set(gca,'ydir','normal')
hc = colorbar;
colormap(winter)
xticks(1:K)
title('Sampling measure $\pi^{(k)}_j$ (Adapt I)','interpreter','latex')
set(gca,'TickLabelInterpreter', 'latex');
set(hc,'TickLabelInterpreter', 'latex')
set(gca,'fontsize',15)
saveas(gca,'fig/Measure_Adapt_I','png')

%% Visualize adapted measure (Variant II)
figure;
imagesc(pi_2)
xlabel('Iteration $k$','interpreter','latex')
ylabel('Index $j$','interpreter','latex')
set(gca,'ydir','normal')
hc = colorbar;
colormap(winter)
xticks(1:K)
title('Sampling measure $\pi^{(k)}_j$ (Adapt II)','interpreter','latex')
set(gca,'TickLabelInterpreter', 'latex');
set(hc,'TickLabelInterpreter', 'latex')
set(gca,'fontsize',15)
saveas(gca,'fig/Measure_Adapt_II','png')

%% Boxplot
figure;

boxplot(log10(Err'), 'labels',{'(Adapt I)','(Adapt II)','(Unif I)','(Unif II)','(Cheby)'});
title('Adaptive {\it vs.} nonadaptive sampling','interpreter','latex')
%xlabel('Sampling strategy','interpreter','latex')
ylabel('$\log_{10}(\|\hat{x}-x\|_2)$','interpreter','latex')
set(gca,'TickLabelInterpreter', 'latex');
grid on
set(gca,'fontsize',15)
saveas(gca,'fig/Comparison','epsc')


%% plot of f(x)
figure;
t_grid = linspace(-1,1,1000);
plot(t_grid,LegMat(t_grid,n)*x_exact);
title('last randomly generated function $f$','interpreter','latex')
set(gca,'TickLabelInterpreter', 'latex');
grid on
set(gca,'fontsize',15)


% %% plot of LS solution
% figure(5);
% A_grid = LegMat(t_grid,n);
% plot(t_grid, A_grid * x_exact);
% hold on
% plot(t_grid, A_grid * x_adapt_1);
% plot(t_grid, A_grid * x_adapt_2);
% plot(t_grid, A_grid * x_unif_cont);
% plot(t_grid, A_grid * x_unif_disc);
% plot(t_grid, A_grid * x_Cheby);
% hold off
% legend('exact','(Adapt I)','(Apdapt II)','(Unif I)','(Unif II)','(Cheb)','location','southeast','interpreter','latex')
% title('Various approximations to $f$','interpreter','latex')
% set(gca,'TickLabelInterpreter', 'latex');
% grid on
% set(gca,'fontsize',15)

%% comparison of coefficients (last run)
figure;
stem(x_exact,'o');
hold on
stem(x_adapt_1,'x');
stem(x_adapt_2,'x');
stem(x_unif_cont,'x');
stem(x_unif_disc,'x');
stem(x_Cheby,'x');
hold off
legend('exact','Adapt I','Apdapt II','Unif I','Unif II','Cheby')




