function obj = logistic(data, opt,verbose)
% DESCRIPTION:
% Return regularization path of logit model.
%% Initialization of Data %%
obj.class = 'logistic split';
y = data.y;
X = data.X;
[n, p] = size(X);
%% Initialization of Parameter %%
opt = lbi.initial(opt);
if opt.normalize == true, X = normc(X);end
kappa = opt.kappa;
%% Initialization of \delta %%
if isempty(opt.delta)
    delta = opt.c * n / kappa / (1 + norm(X, 2)^2);
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
    X_tilde = [ones(n,1),X];
    beta_tilde = lbi.logistic_minimize(X_tilde,y);
    beta0 = beta_tilde(1);
    beta = beta_tilde(2:end);
else
    beta0 = 0;
    beta = lbi.logistic_minimize(X,y);
end
z = zeros(p, 1);
obj.beta0 = repmat(beta0, 1, opt.t_num);
obj.beta = repmat(beta, 1, opt.t_num);
obj.z = zeros(p, opt.t_num);
obj.cost = zeros(1, opt.t_num);
obj.var_hist = [];
obj.var_order = [];

%% The regularization path from 0 to t0 %%
if isempty(opt.t_seq)
    tmp = - y ./ (1 + exp(repmat(beta0, n, 1) .* y));
    d_beta = X' * tmp / n;
    if ~isfield(opt,'t0')
        t0 = 1 / max(max(abs(d_beta)));
    else
        t0 = opt.t0;
    end
    opt.t_seq = logspace(log10(t0), log10(t0 * opt.t_ratio), opt.t_num);
    if opt.fast_init
        z = z - t0 * d_beta;
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
if verbose, fprintf(['Linearized Bregman Iteration (', obj.class, '):\n']); end
tic
var_hist = [];
var_order = [];
for step_cur = 1:steps_remain
    if rec_cur > opt.t_num, break; end
    %% update \beta,z %%
    tmp = - y ./ (1 + exp((X * beta + repmat(beta0, n, 1)) .* y));
    if opt.intercept
        d_beta0 = mean(tmp);
    end
    d_beta = X' * tmp / n;
    if opt.intercept
        beta0 = beta0 - kappa * delta * d_beta0;
    end
    z = z - delta * d_beta;
    beta = kappa * sign(z) .* max(abs(z) - 1, 0);
    %% update var_hist and var_order %%
    if opt.auc
        beta_index = find(beta);
        if sum(~ismember(beta_index,var_hist)) > 0
            var_order = [var_order;step_cur * ones(sum(~ismember(beta_index,var_hist)),1)];
            var_hist = [var_hist;beta_index(~ismember(beta_index,var_hist))];
        end
    end
    %% Recording some of estimations in the regularization path %%
    while true
        dt = step_cur * delta + t0 - opt.t_seq(rec_cur);
        if dt < 0, break; end
        %% update \beta(0),\beta,z and \gamma %%
        if opt.intercept
            obj.beta0(rec_cur) = beta0 + kappa * dt * d_beta0;
        end
        obj.z(:, rec_cur) = z + dt * d_beta;
        obj.beta(:, rec_cur) = kappa * sign(z + dt * d_beta) .* max(abs(z + dt * d_beta) - 1, 0);
        obj.cost(rec_cur) = toc * 1;
        rec_cur = rec_cur + 1;
        if rec_cur > opt.t_num, break; end
    end
    if verbose && ismember(step_cur, round(steps_remain ./ [100 50 20 10 5 2 1]))
        fprintf('Process: %0.2f%%. Time: %f\n', step_cur / steps_remain * 100, toc);
    end
end
fprintf('\n');
obj.delta = delta;
obj.t_seq = opt.t_seq;
obj.K = length(opt.t_seq);
obj.var_hist = var_hist;
obj.var_order = var_order;
end

