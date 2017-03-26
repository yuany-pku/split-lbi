% DESCRIPTION 
% Compute average auc for every \nu 
%% Setting %%
run ../Linear_setting_sim.m;
clear opt;
opt.kappa = 200;
opt.nu = nu0;
opt.t_ratio = 20;
opt.fast_init = false;
opt.intercept = false;
plot_roc = false;
opt.auc = true;
%% Computation of AUC %%
auc_nu_fused = zeros(100,1);
auc_nu_identity = zeros(100,1);
for i = 1:100
    X = load(strcat('../Data/X',num2str(i),'_lbi_linear.txt'));
    y = load(strcat('../Data/y',num2str(i),'_lbi_linear.txt'));
    % when D is fused lasso %
    D = [eye(p - 1), zeros(p - 1, 1); eye(p)] - [zeros(p - 1, 1), eye(p - 1); zeros(p)];
    d_true = D * beta_true;
    I_fused = find(d_true ~= 0);
    data.y = y;
    data.X = X;
    data.D = D;
    obj = splitlbi.splitlbi(data,opt);
    var_hist_fused_tmp = obj.var_hist;
    var_order_fused_tmp = obj.var_order;
    auc_nu_fused(i) = splitlbi.ROC_AUC(var_hist_fused_tmp,var_order_fused_tmp,I_fused,size(D,1),plot_roc);
    % when D = I %
    D = eye(p);
    d_true = D * beta_true;
    I_identity = find(d_true ~= 0);
    data.y = y;
    data.X = X;
    data.D = D;
    obj = splitlbi.splitlbi(data,opt);
    var_hist_identity_tmp = obj.var_hist;
    var_order_identity_tmp = obj.var_order;
    auc_nu_identity(i) = splitlbi.ROC_AUC(var_hist_identity_tmp,var_order_identity_tmp,I_identity,size(D,1),plot_roc);
end