% DESCRIPTION 
% This file provides a simple example to start with split lbi. The split
% lbi is implemented to recover signal vector with group sparsity
% sparsity under linear model. Including computation of auc, visualization of roc 
% and regularization solution path.
%% Setting %%
run ../exp_simulation/Linear_setting_group_sim.m;
X = load(strcat('../exp_simulation/Data/X',num2str(1),'_lbi_linear.txt'));
rng(1);
y = X * beta_true + randn(n,1)*sigma*0.5;
D = eye(p);
m = size(D,1);
d_true = D * beta_true;
I = find(d_true~=0);
%% Split LBI %%
clear opt;
opt.kappa = 200;
opt.nu = 1;
opt.t_ratio = 20;
opt.fast_init = false;
opt.intercept = false;
opt.t_num = 100;
opt.auc = true;
data.y = y;
data.X = X;
data.D = D;
group_type = 'grouped';
lbi = splitlbi.splitlbi(data,opt,[],group_type,group_index);
%% Computation of AUC %%
var_hist = lbi.var_hist;
var_order = lbi.var_order;
plot_roc = true;
name = 'Linear_group';
auc = splitlbi.ROC_AUC(var_hist,var_order,I,m,plot_roc,name);
    
%% Cross validation %%
obj_cross = splitlbi.cv_splitlbi(data,opt,[],[],group_type,group_index);
fprintf('The l2 norm of beta after cross-validation is %0.2f\n',norm(obj_cross.beta - beta_true));
fprintf('The l2 norm of beta_tilde after cross-validation is %0.2f\n',norm(obj_cross.beta_tilde - beta_true));