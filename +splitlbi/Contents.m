% split lbi filter
%
%   splitlbi.splitlbi                - main entry point for package(start here)
%
%   splitlbi.cv_splitlbi             - cross validation for split lbi
%   splitlbi.linear_split            - split lbi for linear model with ungrouped sparsity
%   splitlbi.logistic_split          - split lbi for logit model with ungrouped sparsity
%   splitlbi.linear_grouped_split    - split lbi for linear model with grouped sparsity
%   splitlbi.logistic_split          - split lbi for logit model with grouped sparsity
%   splitlbi.initial                 - Initialize parameters for split lbi.
%   splitlbi.linear_minimize         - Return MLE estimator \beta under linear model
%   splitlbi.logistic_minimize       - Return MLE estimator \beta under logit model
%   splitlbi.proj_Dt                 - return \tilde{\beta} by project \beta on D_{S}
%   splitlbi.plot_splitlbi           - plot regularization solution path of split lbi and lbi
%   splitlbi.ROC_AUC                 - plot ROC curve in simulation and compute AUC
%   splitlbi.table_grouping          - Return Grouping results of candidates
