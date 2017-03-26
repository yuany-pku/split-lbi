% DESCRIPTION 
% This file compare IRR(\nu) of split lbi and IC condition of genlasso
% under setting of fused lasso 
%% Setting %%
run ../Linear_setting_sim.m;
D = [eye(p - 1), zeros(p - 1, 1); eye(p)] - [zeros(p - 1, 1), eye(p - 1); zeros(p)];
d_true = D * beta_true;
I = find(d_true~=0);
%% Compute IRR %%
log_nu = -20:0.1:20;
nulist = exp(log_nu);
n_nu = length(nulist);
irr_nu_mean = zeros(n_nu,1);
irr_nu_std = zeros(n_nu,1);
fprintf('It will take a for a while to finish computation of IRR...\n');
tic;
for nu_num = 1:n_nu
    irr_nu = zeros(100,1);
    nu = nulist(nu_num);
    for j = 1:100
        X = load(strcat('../Data/X',num2str(j),'_lbi_linear.txt'));
        irr_nu(j) = Check_irr(X,D,I,nu);
    end
    if ismember(nu_num, round(n_nu ./ [100 50 20 10 5 2 1]))
        fprintf('Process: %0.2f%%. Time: %f\n', nu_num / n_nu * 100, toc);
    end
    irr_nu_mean(nu_num) = mean(irr_nu);
    irr_nu_std(nu_num) = std(irr_nu);
end
%% Compute IC_0 and IC_1 %%
ic_0 = zeros(100,1);
ic_1 = zeros(100,1);
for i=1:100
    X = load(strcat('../Data/X',num2str(i),'_lbi_linear.txt'));
    ic_0(i) = Check_IC_0(X, D, I);
    ic_1(i) = Check_IC_1(X, D, I, d_true);
end
ic_0_mean = mean(ic_0);
ic_0_std = std(ic_0);
ic_1_mean = mean(ic_1);
ic_1_std = std(ic_1);

%% Plot: IRR($$\nu$$) v.s. \nu and IC_0,IC_1
[~,ind_change] = min(abs(irr_nu_mean - 1));
figure;
plot(log_nu,irr_nu_mean,'Color','b');
hold on;
plot([min(log_nu),max(log_nu) + 1],[ic_0_mean,ic_0_mean],'Color','r','LineStyle','-');
text(min(log_nu) + 0.5,ic_0_mean - 0.2,'$$IC_0$$','interpreter','latex');
hold on;
plot([min(log_nu),max(log_nu) + 1],[ic_1_mean,ic_1_mean],'Color','r','LineStyle','--');
text(min(log_nu) + 0.5,ic_1_mean - 0.2,'$$IC_1$$','interpreter','latex');
hold on;
plot([min(log_nu),max(log_nu) + 1],[1,1],'Color','r','LineStyle','-.');
hold on;
plot([log_nu(ind_change),log_nu(ind_change)],[0,max([irr_nu_mean;ic_0_mean;ic_1_mean]) + 0.5],'Color',[0.5 0 0.5],'LineStyle','-.');
hold on;
xlabel('$$log\nu$$','Interpreter','latex');
ylabel('$$IRR(\nu)$$','Interpreter','latex');
xlim([min(log_nu)-0.02,max(log_nu) + 0.02]);
ylim([0 - 0.02,max([irr_nu_mean;ic_0_mean;ic_1_mean]) + 0.5]);
title('$$IRR(\nu)$$ v.s $\nu$','Interpreter','latex');
saveas(gcf, [strcat('./','IRR_IC'), '.png']);

    
    