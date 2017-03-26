% DESCRIPTION 
% This file provides a simple example to start with split lbi. The split
% lbi is implemented to recover a linear model with group sparsity
% under logit model. Including computation of auc, visualization of roc 
% and regularization solution path.
%% Setting %%
run ../exp_simulation/Logit_setting_group_sim.m;
rng(1);
X = randn(n,p);
prob = exp(X * beta_true) ./ (1 + exp(X * beta_true));
y = binornd(1,prob,n,1);
y(y == 0) = -1;
D = eye(p);
m = size(D,1);
d_true = D * beta_true;
I = find(d_true~=0);
%% Split LBI %%
clear opt;
opt.kappa = 50;
opt.nu = 1;
opt.t_ratio = 40;
opt.fast_init = false;
opt.intercept = false;
opt.t_num = 100;
opt.auc = true;
data.y = y;
data.X = X;
data.D = D;
family = 'binomial';
group_type = 'grouped';
obj = splitlbi.splitlbi(data,opt,family,group_type,group_index);
%% Computation of AUC %%
var_hist = obj.var_hist;
var_order = obj.var_order;
plot_roc = true;
name = 'Logit_group';
auc = splitlbi.ROC_AUC(var_hist,var_order,I,m,plot_roc,name);

%% Cross validation %%
obj_cross = splitlbi.cv_splitlbi(data,opt,[],family,group_type,group_index);