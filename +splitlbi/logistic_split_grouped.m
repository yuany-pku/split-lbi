function obj = logistic_split_grouped(data, opt,group_index,verbose)
% DESCRIPTION
% Return regularization path with group sparsity under logit model.
%% Initialization of Data %%
obj.class = 'logistic split';
y = data.y;
X = data.X;
D = data.D;
[n, p] = size(X);
m = size(D, 1);
%% Initialization of Parameter %%
opt = splitlbi.initial(opt);
if opt.normalize == true, X = normc(X);end
kappa = opt.kappa;
%% Initialization of \nu %%
if isempty(opt.nu)
    nu = n * norm(full(D), 2)^2 / norm(full(bsxfun(@minus, X, mean(X))), 2)^2; 
else
    nu = opt.nu;
end
%% Initialization of \delta %%
if isempty(opt.delta)
    delta = opt.c * nu / kappa / (1 + nu * norm(full(bsxfun(@minus, X, mean(X))), 2)^2 / n + norm(full(D), 2)^2);
else
    delta = opt.delta;
end
%% Initialize t_seq, t_ratio and t_num %%
if isempty(opt.t_seq)
    if isempty(opt.t_ratio)
        if n < p, opt.t_ratio = 10; else opt.t_ratio = 100; end
    elseif opt.t_ratio <= 1, error('t_max/t_min should be larger than 1.');
    end
else
    opt.t_seq = sort(opt.t_seq);
    if opt.t_seq(1) < 0, error('Time should be non-negative.'); end
    opt.t_num = length(opt.t_seq);
end
if nargin < 3, verbose = true; end
%% Initialize \gamma(0),z(0),\beta(0) and \beta_tilde(0) to zeros %%
if opt.intercept
    X_tilde = [ones(n, 1), X]; D_tilde = [zeros(m, 1), D];
    beta_tilde = splitlbi.logistic_minimize(X_tilde,y,D_tilde,nu);
    beta0 = beta_tilde(1);
    beta = beta_tilde(2:end);
else
    beta0 = 0;
    beta = splitlbi.logistic_minimize(X,y,D,nu);
end
z = zeros(m, 1);
gamma = zeros(m, 1);
obj.beta0 = repmat(beta0, 1, opt.t_num);
obj.beta = repmat(beta, 1, opt.t_num);
obj.z = zeros(m, opt.t_num);
obj.gamma = zeros(m, opt.t_num);
obj.beta_tilde = repmat(beta, 1, opt.t_num);
obj.cost = zeros(1, opt.t_num);
obj.var_hist = [];
obj.var_order = [];


%% The regularization path from 0 to t0 %%
if isempty(opt.t_seq)
    d_gamma = - D * beta / nu;
    if ~isfield(opt,'t0')
        t0 = 1 / max(max(abs(d_gamma)));
    else
        t0 = opt.t0;
    end
    opt.t_seq = logspace(log10(t0), log10(t0 * opt.t_ratio), opt.t_num);
    if opt.fast_init
        z = z - t0 * d_gamma;
    else
        beta = zeros(p,1);
        t0 = 0;
    end
else
    t0 = 0;
end
rec_cur = sum(opt.t_seq <= t0) + 1;
steps_remain = ceil((opt.t_seq(end) - t0) / delta);
fprintf('The number of whole iteration %d\n',steps_remain);
%% Starting Iteration %%
G = length(unique(group_index));
if verbose, fprintf(['Linearized Bregman Iteration (', obj.class, '):\n']); end
tic
var_hist = [];
var_order = [];
for step_cur = 1:steps_remain
    if rec_cur > opt.t_num, break; end
    %% update \beta,z and \gamma %%
    Xbeta = X * beta;
    tmp = - y ./ (1 + exp((Xbeta + repmat(beta0, n, 1)) .* y));
    if opt.intercept
        d_beta0 = mean(tmp);
    end
    d_gamma = (gamma - D * beta) / nu;
    z = z - delta * d_gamma;
    for g = 1:G
        g_ind = find(group_index == g);
        gamma(g_ind) = kappa * max(0, 1 - 1 / norm(z(g_ind))) * z(g_ind);
    end
    if opt.intercept
        beta0 = beta0 - kappa * delta * d_beta0;
    end
    tmp = - y ./ (1 + exp(Xbeta .* y));
    d_beta = X' * tmp / n - D' * d_gamma;
    beta = beta - kappa * delta * d_beta;
    
    %% update var_hist and var_order %%
    if opt.auc 
        gamma_index = find(gamma);
        index_add = ~ismember(gamma_index,var_hist);
        if sum(index_add) > 0
            var_order = [var_order;step_cur * ones(sum(index_add),1)];
            var_hist = [var_hist;gamma_index(index_add)];
        end
    end
    %% Recording some of estimations in the regularization path %%
    
    while true
        dt = step_cur * delta + t0 - opt.t_seq(rec_cur);
        if dt < 0, break; end
        % update \beta(0),\beta,z and \gamma %
        if opt.intercept
            obj.beta0(:, rec_cur) = beta0 + kappa * dt * d_beta0;
        end
        obj.z(:, rec_cur) = z + dt * d_gamma;
        for g = 1:G
            g_ind = find(group_index == g);
            obj.gamma(g_ind,rec_cur) = kappa * max(0, 1 - 1 / norm(obj.z(g_ind,rec_cur))) * obj.z(g_ind,rec_cur);
        end
        obj.beta(:, rec_cur) = beta + kappa * dt * d_beta;
        obj.beta_tilde(:, rec_cur) = splitlbi.proj_Dt(obj.gamma(:,rec_cur),obj.beta(:,rec_cur),D);
        obj.cost(rec_cur) = toc * 1;
        rec_cur = rec_cur + 1;
        if rec_cur > opt.t_num, break; end
    end
    if verbose && ismember(step_cur, round(steps_remain ./ [100 50 20 10 5 2 1]))
        fprintf('Process: %0.2f%%. Time: %f\n', step_cur / steps_remain * 100, toc);
    end
end
fprintf('\n');
obj.nu = nu;
obj.delta = delta;
obj.t_seq = opt.t_seq;
obj.K = length(opt.t_seq);
obj.var_hist = var_hist;
obj.var_order = var_order;
end