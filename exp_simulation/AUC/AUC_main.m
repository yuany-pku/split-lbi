% DESCRIPTION
% This file compute average auc of split lbi
% with setting of fused lasso and standard lasso 
% for different \nu under linear model
%% Computation of AUC v.s \nu %%
nulist = [1,5,10]';
n_nu = length(nulist);
auc_mean_fused = zeros(n_nu,1);
auc_std_fused = zeros(n_nu,1);
auc_mean_identity = zeros(n_nu,1);
auc_std_identity = zeros(n_nu,1);
for nu_num = 1:n_nu
    nu0 = nulist(nu_num);
    run auc_each;
    auc_mean_fused(nu_num) = mean(auc_nu_fused);
    auc_std_fused(nu_num) = std(auc_nu_fused);
    auc_mean_identity(nu_num) = mean(auc_nu_identity);
    auc_std_identity(nu_num) = std(auc_nu_identity);
end

%% Table %%
T = table(nulist,auc_mean_fused,auc_std_fused,auc_mean_identity,auc_std_identity...
            ,'VariableNames',{'nu','Fused_Average_AUC','Fused_Std_AUC'...
            ,'Identity_Average_AUC','Identity_Std_AUC'});
writetable(T,'AUC.csv','Delimiter',',');