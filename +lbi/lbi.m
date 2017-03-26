function obj = lbi(data,opt,family,group_type,group_index,verbose)
%--------------------------------------------------------------------------
% lbi.m: fit an GLM with l1 sparsity
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Fit a generalized linear model via penalized maximum likelihood. The
%    regularization path can return estimation with l1 sparsity.
%    Fits linear and logistic models.
%
% USAGE:
%    obj = lbi.lbi(data)
%    obj = lbi.lbi(data, opt)
%    obj = lbi.lbi(data, opt, family)
%    obj = lbi.lbi(data, opt, family, group_type, group_index)
%    obj = lbi.lbi(x, y, family, group_type, group_index, verbose)
%    (Use empty matrix [] to apply the default value for 'family', eg.
%    obj = lbi.lbi(data, opt, [], group_type, group_index))
%
% INPUT ARGUMENTS:
% data        data.X, dimension nobs x nvars; each row is an
%             observation vector. Can be in sparse matrix format.
%             data.y: dimension nobs x 1; response variable, can
%             be in sparse matrix format.
%
% opt         The paramter setting for lbi, such as kappa,alpha,
%             intercept, fast_init, auc. If defined is empty, we initialize
%             it using +lbi/initial.m, you can reference this file for details
%             about the meaning of these parameters.
%
% family      Either 'gaussian' or 'binomial'. Default is 'gaussian'. When
%             family = 'gaussian', the loss function is:
%
%                            1/2 RSS / nobs
%
%             When family = 'binomial', the loss function is:
%
%                           -loglik  / nobs
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
% obj         A structure.
% obj.beta0   Intercept, if any.
% obj.beta    The regularization path of \beta(t)
% obj.z       The regularization path of z(t) = rho(t) + \gamma(t) / kappa.
% obj.delta   Step Size.
% obj.K       Length of opt.t_seq;
% obj.class   {'linear','logistic'}
%
%
% DETAILS:
%  The Linearized Bregman solver computes the whole regularization path
%  with general form of l1 sparsity under gaussian, binomial and
%  models through iterations. For binomial models, the response variable y
%  is assumed to be a vector of two classes which is transformed in to \{1,-1\}.
%  The definitions of damping parameter kappa, step size alpha are the
%  same, are the same with those defined in the reference paper.
%
% SEE ALSO:
%   splitlbi.m cv_splitlbi.m linear_split.m
%   logistic_split.m cv_lbi.m linear.m logistic.m
%
% EXAMPLES:
% % Gaussian
%    X = randn(50,50);
%    y = randn(100,1);
%    data.X = X; data.y = y;
%    obj = lbi.lbi(data);
%
% % Binomial:
%    y = binomial(1,0.5,100,1);
%    y(y==0) = -1;
%    data.y = y;
%    family = 'binomial';
%    obj = lbi.lbi(data, [], family);


if ~exist('data','var'), error('No input data!'); end
if ~isfield(data, 'X') || ~isfield(data, 'y')
    error('The input data is not complete,missing X or y!');
else
    if size(data.X,1) ~= length(data.y)
        error('Number of rows of X must be equal to the length of y!');
    end
    if size(data.y,2) > 1
        error('y should be an vector!');
    end
end
p = size(data.X,2);
%% opt %%
if nargin < 2 || isempty(opt), opt = [];end
% match the family
if nargin < 3 || isempty(family), family = 'gaussian';end
fambase = {'gaussian','binomial','poisson','multinomial','cox','mgaussian'};
famind = find(strncmp(family,fambase,length(family)),1);
if isempty(famind)
    error('family should be one of ''gaussian'', ''binomial''');
else
    family = fambase{famind};
    if strcmp(family,'binomial') && sum(abs(data.y)~=1) > 0, error('y must be in {1,-1}');end
end

% match the group_type
if nargin < 4 || isempty(group_type), group_type = 'ungrouped';end
grobase = {'ungrouped','grouped'};
groind = find(strncmp(group_type,grobase,length(group_type)),1);
if isempty(groind)
    error('group type should be either ''grouped'' or ''ungrouped''');
else
    group_type = grobase{groind};
end

% group index
if nargin < 5 || isempty(group_index)
    group_index = 1:p;
else
    if length(group_index) ~= p
        error('The length of unique elements of group_index should be equal to the number of rows of D!');
    end
    if max(group_index) > p
        error('The group number should less than number of rows of D!');
    end
end
if nargin < 6, verbose = true; end
% Fit model
if strcmp(group_type,'ungrouped')
    switch family
        case 'gaussian'
            obj = lbi.linear(data,opt,verbose);
        case 'binomial' 
            obj = lbi.logistic(data,opt,verbose);
    end
else
    switch family
        case 'gaussian'
            obj = lbi.linear_grouped(data,opt,group_index,verbose);
        case 'binomial' 
            obj = lbi.logistic_grouped(data,opt,group_index,verbose);
    end
end
%------------------------------------------------------------------
% End function lbi
%------------------------------------------------------------------