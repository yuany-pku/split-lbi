function obj = splitlbi(data,opt,family,group_type,group_index,verbose)
%--------------------------------------------------------------------------
% splitlbi.m: fit an GLM with structured sparsity
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Fit a generalized linear model via penalized maximum likelihood. The 
%    regularization path can return estimation with structured sparsity.
%    Fits linear and logistic models.
%
% USAGE:
%    obj = splitlbi(data)
%    obj = splitlbi(data, opt)
%    obj = splitlbi(data, opt, family)
%    obj = splitlbi(x, y, family, group_type,group_index)
%    obj = splitlbi(x, y, family, group_type,group_index,verbose)
%    (Use empty matrix [] to apply the default value for 'family', eg. 
%    obj = splitlbi(data, opt, [], group_type,group_index))
%
% INPUT ARGUMENTS:
% data        data.X, dimension nobs x nvars; each row is an
%             observation vector. Can be in sparse matrix format. 
%             data.y: dimension nobs x 1; response variable, can
%             be in sparse matrix format.
%             data.D, dimension m x nvars; matrix corresponding to
%             sparsity form, can be in sparse matrix format.
%
% opt         The paramter setting for split lbi, such as kappa,alpha,
%             nu, intercept, fast_init, auc. If defined is empty, we
%             initialize it using +splitlbi/initial.m, you can reference
%             this file for details about the meaning of these parameters.
%             
% family      Either 'gaussian' or 'binomial'. Default is 'gaussian'. When
%             family = 'gaussian', the loss function is:
%
%              1/2 ||y - X\beta ||^2 / nobs + 1/2 ||D\beta - \gmma||^2
%              
%             When family = 'binomial', the loss function is:
%
%                    -loglik  / nobs + 1/2 ||D\beta - \gmma||^2
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
% obj.gamma   The regularization path of \gamma(t)  
% obj.nu      The paramter for variable splitting term.
% obj.delta   Step Size.
% obj.K       Length of opt.t_seq;
% obj.class   {'linear split','logistic split'}
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
% LICENSE: GPL-2
%
% DATE: 26 Mar 2017
%
% AUTHORS:
%    Chendi Huang, Xinwei Sun, Jiechao Xiong, Yuan Yao
%
% REFERENCES:
%    Huang C, Sun X, Xiong J, et al. 
%    Split LBI: An Iterative Regularization Path with Structural Sparsity.
%    Advances In Neural Information Processing Systems. 2016: 3369-3377.
%
%    Vaiter S, Peyré G, Dossal C, et al. Robust sparse analysis regularization.
%    IEEE Transactions on information theory, 2013, 59(4): 2001-2016.
%
%    Osher S, Ruan F, Xiong J, et al.
%    Sparse recovery via differential inclusions.
%    Applied and Computational %Harmonic Analysis, 2016, 41(2): 436-469.
%
% SEE ALSO:
%    cv_splitlbi.m linear_split.m logistic_split.m 
%    lbi.m cv_lbi.m linear.m logistic.m
%
% EXAMPLES:
% % Gaussian
%    X = randn(50,50);
%    y = randn(100,1);
%    data.X = X; data.y = y;
%    data.D = eye(50);
%    obj = splitlbi.splitlbi(data);
%    
% % Binomial:
%    y = binomial(1,0.5,100,1);
%    y(y==0) = -1;
%    data.y = y;
%    family = 'binomial';
%    obj = splitlbi.splitlbi(data, [], family);
%
% DEVELOPMENT:
%    26 Mar 2017: Original version of splitlbi.m written.

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
if ~isfield(data, 'D') || isempty(data.D)
    D = eye(size(X,2));
else
    D = data.D;
    if size(D,2) ~= size(data.X,2), error('Number of columns of D and X must be equal!'); end
end
data.D = D;
m = size(D,1);
%% opt %%
if nargin < 2 || isempty(opt), opt = [];end
%% match the family %%
if nargin < 3 || isempty(family), family = 'gaussian';end
fambase = {'gaussian','binomial','poisson','multinomial','cox','mgaussian'};
famind = find(strncmp(family,fambase,length(family)),1);
if isempty(famind)
    error('family should be one of ''gaussian'', ''binomial''');
else
    family = fambase{famind};
    if strcmp(family,'binomial') && sum(abs(data.y)~=1) > 0, error('y must be in {1,-1}');end
end
%% match the group_type %%
if nargin < 4 || isempty(group_type), group_type = 'ungrouped';end
grobase = {'ungrouped','grouped'};
groind = find(strncmp(group_type,grobase,length(group_type)),1);
if isempty(groind)
    error('group type should be either ''grouped'' or ''ungrouped''');
else
    group_type = grobase{groind};
end
%% group index %%
if nargin < 5 || isempty(group_index)
    group_index = 1:m;
else
    if length(group_index) ~= m
        error('The length of unique elements of group index should be equal to the number of rows of D!');
    end
    if max(group_index) > m
        error('The group number should less than number of rows of D!');
    end
end
if nargin < 6, verbose = true; end
%% Fit model %%
if strcmp(group_type,'ungrouped')
    switch family
        case 'gaussian'
            obj = splitlbi.linear_split(data,opt,verbose);
        case 'binomial' 
            obj = splitlbi.logistic_split(data,opt,verbose);
    end
else
    switch family
        case 'gaussian'
            obj = splitlbi.linear_split_grouped(data,opt,group_index,verbose);
        case 'binomial' 
            obj = splitlbi.logistic_split_grouped(data,opt,group_index,verbose);
    end
end
%------------------------------------------------------------------
% End function splitlbi
%------------------------------------------------------------------