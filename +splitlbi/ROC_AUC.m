function AUC = ROC_AUC(var_select,var_order,var_true,m,plot_roc,name)
% DESCRIOTION
% This file is used to compute AUC and plot ROC for splitlbi and lbi
if nargin < 3, error('The input for computaion of auc is incomplete!'); end
var_len = length(var_select);
if var_len ~= length(var_order)
    error('The length of variable selected and that of order vector should be equal!');
end
if sum(ismember(var_true,var_select)) < length(var_true)
    warning('The selected variable does not contain all true variables, the computation is inaccurate!');
end
if var_len > m || max(var_select) > m
    error('The length of selected variable vector should not larger than all variables!');
end
labels = zeros(m,1);
labels(var_true) = 1;
scores = max(min(var_order) / max(var_order) - 0.01,0) * ones(m,1);
scores(var_select) = sort(var_order,'descend') / max(var_order);
[x,y,~,AUC] = perfcurve(labels,scores,1);
if nargin < 5 || isempty(plot_roc), plot_roc = false; end
if plot_roc
    figure;
    plot(x,y);
    if nargin < 6 || isempty(name)
        saveas(gcf, [strcat('./','ROC_CURVE'), '.png']);
    else
        saveas(gcf, [strcat('./',name,'_ROC_CURVE'), '.png']);
    end
end
end
