% DESCRIPTION 
% This file provides a simple example to start with split lbi. The split
% lbi is implemented to recover signal vector with 1-d fused lasso
% sparsity under linear model. Including computation of auc, visualization of roc 
% and regularization solution path.
%% Setting %%
run ../exp_simulation/Linear_setting_sim.m;
X = load(strcat('../exp_simulation/Data/X',num2str(1),'_lbi_linear.txt'));
y = load(strcat('../exp_simulation/Data/y',num2str(1),'_lbi_linear.txt'));
D = [eye(p - 1), zeros(p - 1, 1); eye(p)] - [zeros(p - 1, 1), eye(p - 1); zeros(p)];
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
obj = splitlbi.splitlbi(data,opt);
%% Computation of AUC %%
var_hist = obj.var_hist;
var_order = obj.var_order;
plot_roc = true;
name = 'Linear';
auc = splitlbi.ROC_AUC(var_hist,var_order,I,m,plot_roc,name);
    
%% Cross validation %%
obj_cross = splitlbi.cv_splitlbi(data,opt);
fprintf('The l2 norm of beta after cross-validation is %0.2f\n',norm(obj_cross.beta - beta_true));
fprintf('The l2 norm of beta_tilde after cross-validation is %0.2f\n',norm(obj_cross.beta_tilde - beta_true));