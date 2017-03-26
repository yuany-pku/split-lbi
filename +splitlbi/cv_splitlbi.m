function obj_cross = cv_splitlbi(data,opt,kfolds,family,group_type,group_index,verbose)
%---------------------------------------------------------------------------------------
% cv_splitlbi.m: Cross-validation method to tuning the parameter t for splitlbi.
%---------------------------------------------------------------------------------------
%
% DESCRIPTION:
%    K-fold cross-validation method is used to tuning the parameter t for
%    split lbi. Mean square error is used for linear model. Miss-classification
%    error is used for binomial model.
%
% USAGE:
%    obj_cross = splitlbi.cv_splitlbi(data)
%    obj_cross = splitlbi.cv_splitlbi(data, opt)
%    obj_cross = splitlbi.cv_splitlbi(data, opt, family)
%    obj_cross = splitlbi.cv_splitlbi(data, opt, family, group_type, group_index)
%    obj_cross = splitlbi.cv_splitlbi(data, opt, family, group_type, group_index,verbose)
%    (Use empty matrix [] to apply the default value for 'family', eg.
%    obj_cross = splitlbi.cv_splitlbi(data, opt, [], group_type, group_index)
%
% INPUT ARGUMENTS:
% data        data.X, dimension nobs x nvars; each row is an
%             observation vector. Can be in sparse matrix format.
%             data.y, dimension nobs x 1; response variable, can
%             be in sparse matrix format.
%             data.D, dimension m x nvars; matrix corresponding to
%             sparsity form, can be in sparse matrix format.
%
% opt         The paramter setting for split lbi, such as kappa,alpha,
%             nu, intercept, fast_init, auc. If defined is empty, we
%             initialize it using +splitlbi/initial.m, you can reference this file
%             for details about the meaning of these parameters.
%
% kfolds      K Folds number for CV. Default is 5.
%
% family      Either 'gaussian' or 'binomial'. Default is 'gaussian'. When
%             family = 'gaussian', the loss function is:
%
%              1/2 ||y - X\beta ||^2 / nobs + 1/2 ||D\beta - \gmma||^2
%
%             When family = 'binomial', the loss function is:
%
%                  -loglik  / nobs + 1/2 ||D\beta - \gmma||^2
%
% group_type  Either 'grouped' or 'ungrouped', representing grouped
%             sparsity or not. Default is 'ungrouped'.
%
% group_index Only be used When group_type = 'grouped'. It is an index
%             vector with each element representing the number of group
%             this element belongs to.
%
% verbose     Output the process or not. Default is true.
%
%
% OUTPUT ARGUMENTS:
% obj_cross         A structure.
% obj_cross.t       The early stopping time t chosen by cross validation.
% obj_cross.beta0   Intercep  at time t, if any
% obj_cross.beta    \beta(t)
% obj_cross.gamma   \gamma(t)
% obj_cross.z       z(t) = rho(t) + \gamma(t) / kappa.
% obj_cross.delta   Step Size.
%
%
% DETAILS:
%  The Split Linearized Bregman solver computes the whole regularization path
%  with general form of structure sparsity under gaussian, binomial and
%  models through iterations. For binomial models, the response variable y
%  is assumed to be a vector of two classes which is transformed in to \{1,-1\}.
%  The definitions of damping parameter kappa, step size alpha are the
%  same, and nu in variable splitting term are the same with those defined
%  in the reference paper.
%
% SEE ALSO:
%   splitlbi.m  linear_split.m logistic_split.m
%   lbi.m cv_lbi.m linear.m logistic.m
%
% EXAMPLES:
% % Gaussian
%    X = randn(50,50);
%    y = randn(100,1);
%    data.X = X; data.y = y;
%    data.D = eye(p);
%    kfolds = 5;
%    obj_cross = splitlbi.cv_splitlbi(data, [], kfolds);
%
% % Binomial:
%    y = binomial(1,0.5,100,1);
%    y(y==0) = -1;
%    data.y = y;
%    family = 'binomial';
%    obj_cross = splitlbi.cv_splitlbi(data, [], kfolds, family);

if ~exist('data','var'), error('No input data!'); end
if ~isfield(data, 'X') || ~isfield(data, 'y')
    error('The input data is not complete,missing X or y!');
else
    X = data.X;
    y = data.y;
    [n,p] = size(X);
    if n ~= length(y)
        error('Number of rows of X must be equal to the length of y!');
    end
    if size(y,2) > 1
        error('y should be an vector!');
    end
end
if ~isfield(data, 'D') || isempty(data.D)
    D = eye(p);
else
    D = data.D;
    if size(D,2) ~= p, error('Number of columns of D and X must be equal!'); end
end
data.D = D;
m = size(D,1);
%% opt %%
if nargin < 2 || isempty(opt), opt = [];end
%% Number of folds %%
if nargin < 3 || isempty(kfolds)
    kfolds = 5;
else if kfolds < 2
        error('The number of folds should greater than 1');
    end
end
%% match the family %%
if nargin < 4 || isempty(family), family = 'gaussian';end
fambase = {'gaussian','binomial','poisson','multinomial','cox','mgaussian'};
famind = find(strncmp(family,fambase,length(family)),1);
if isempty(famind)
    error('family should be one of ''gaussian'', ''binomial''');
else
    family = fambase{famind};
    if strcmp(family,'binomial') && sum(abs(data.y)~=1) > 0, error('y must be in {1,-1}');end
end

%% match the group_type %%
if nargin < 5 || isempty(group_type), group_type = 'ungrouped';end
grobase = {'ungrouped','grouped'};
groind = find(strncmp(group_type,grobase,length(group_type)),1);
if isempty(groind)
    error('group type should be either ''grouped'' or ''ungrouped''');
else
    group_type = grobase{groind};
end
%% group index %%
if nargin < 6 || isempty(group_index)
    group_index = 1:m;
else
    if length(group_index) ~= m
        error('The length of unique elements of group_index should be equal to the number of rows of D!');
    end
    if max(group_index) > m
        error('The group number should less than number of rows of D!');
    end
end
if nargin < 7, verbose = true; end
%% Fit model %%
cp = cvpartition(y,'k',kfolds);
test_error_vec = zeros(opt.t_num,1);
if strcmp(group_type,'ungrouped')
    switch family
        case 'gaussian'
            obj = splitlbi.linear_split(data,opt,verbose);
            opt.t_seq = obj.t_seq;
            opt.delta = obj.delta;
            for k = 1:kfolds
                ind_train = find(training(cp,k));
                ind_test = find(test(cp,k));
                X_train = X(ind_train,:);
                y_train = y(ind_train,:);
                data_k.X = X_train;
                data_k.y = y_train;
                data_k.D = D;
                obj_k = splitlbi.linear_split(data_k,opt,verbose);
                X_test = X(ind_test,:);
                y_test = y(ind_test,:);
                y_pre =  X_test * obj_k.beta - repmat(obj_k.beta0,size(X_test,1),1);
                test_error_vec = test_error_vec + norms(repmat(y_test,1,opt.t_num) - y_pre,[],1)';
            end
            [~,ind_cross] = min(test_error_vec);
            
        case 'binomial' 
            obj = splitlbi.logistic_split(data,opt,verbose);
            opt.t_seq = obj.t_seq;
            opt.delta = obj.delta;
            for k = 1:kfolds
                ind_train = find(training(cp,k));
                ind_test = find(test(cp,k));
                X_train = X(ind_train,:);
                y_train = y(ind_train,:);
                data_k.X = X_train;
                data_k.y = y_train;
                data_k.D = D;
                obj_k = splitlbi.logistic_split(data_k,opt,verbose);
                X_test = X(ind_test,:);
                y_test = y(ind_test,:);
                Xtest_beta = X_test * obj_k.beta + repmat(obj_k.beta0,size(X_test,1),1);
                y_pred = exp(Xtest_beta) ./ (exp(Xtest_beta) + 1);
                test_error_vec = test_error_vec + sum(repmat(y_test,1,opt.t_num) ~= y_pred,1)';
            end
            [~,ind_cross] = min(test_error_vec);
    end
else
    switch family
        case 'gaussian'
            obj = splitlbi.linear_split_grouped(data,opt,group_index,verbose);
            opt.t_seq = obj.t_seq;
            opt.delta = obj.delta;
            for k = 1:kfolds
                ind_train = find(training(cp,k));
                ind_test = find(test(cp,k));
                X_train = X(ind_train,:);
                y_train = y(ind_train,:);
                data_k.X = X_train;
                data_k.y = y_train;
                data_k.D = D;
                obj_k = splitlbi.linear_split(data_k,opt,verbose);
                X_test = X(ind_test,:);
                y_test = y(ind_test,:);
                y_pre =  X_test * obj_k.beta - repmat(obj_k.beta0,size(X_test,1),1);
                test_error_vec = test_error_vec + norms(repmat(y_test,1,opt.t_num) - y_pre,[],1)';
            end
            [~,ind_cross] = min(test_error_vec);
        case 'binomial' 
            obj = splitlbi.logistic_split_grouped(data,opt,group_index,verbose);
            opt.t_seq = obj.t_seq;
            opt.delta = obj.delta;
            for k = 1:kfolds
                ind_train = find(training(cp,k));
                ind_test = find(test(cp,k));
                X_train = X(ind_train,:);
                y_train = y(ind_train,:);
                data_k.X = X_train;
                data_k.y = y_train;
                data_k.D = D;
                obj_k = splitlbi.logistic_split(data_k,opt,verbose);
                X_test = X(ind_test,:);
                y_test = y(ind_test,:);
                Xtest_beta = X_test * obj_k.beta + repmat(obj_k.beta0,size(X_test,1),1);
                y_pred = exp(Xtest_beta) ./ (exp(Xtest_beta) + 1);
                test_error_vec = test_error_vec + sum(repmat(y_test,1,opt.t_num) ~= y_pred,1)';
            end
            [~,ind_cross] = min(test_error_vec);
    end
end
obj_cross.beta0 = obj.beta0(ind_cross);
obj_cross.beta = obj.beta(:,ind_cross);
obj_cross.beta_tilde = obj.beta_tilde(:,ind_cross);
obj_cross.gamma = obj.gamma(:,ind_cross);
obj_cross.t = obj.t_seq(ind_cross);
obj_cross.delta = obj.delta;
end
%------------------------------------------------------------------
% End function cv_splitlbi
%------------------------------------------------------------------