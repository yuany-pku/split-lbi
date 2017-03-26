function obj = linear_split(data, opt, verbose)
% DESCRIPTION
% Return regularization path of linear model.
%% Initialization of Data %%
obj.class = 'linear split';
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
    nu = n * norm(full(D), 2)^2 / norm(full(X), 2)^2; 
else
    nu = opt.nu;
end
%% Initialization of \delta %%
if isempty(opt.delta)
    delta = opt.c * nu / kappa / (1 + nu * norm(full(X), 2)^2 / n + norm(full(D), 2)^2); 
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
    beta_tilde = splitlbi.linear_minimize(X_tilde,y,D_tilde,nu);
    beta0 = beta_tilde(1);
    beta = beta_tilde(2:end);
else
    beta0 = 0;
    beta = splitlbi.linear_minimize(X,y,D,nu);
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
fprintf('The number of whole iteration: %d\n',steps_remain);
%% Starting Iteration %%
if verbose, fprintf(['Linearized Bregman Iteration (', obj.class, '):\n']); end
tic
var_hist = [];
var_order = [];
for step_cur = 1:steps_remain
    if rec_cur > opt.t_num, break; end
    %% update \beta(0),\beta,z and \gamma %%
    Xbeta = X * beta;
    if opt.intercept
        d_beta0 = beta0 + mean(Xbeta - y);
    end
    d_gamma = (gamma - D * beta) / nu;
    z = z - delta * d_gamma;
    gamma = kappa * sign(z) .* max(abs(z) - 1, 0);
    if opt.intercept
        beta0 = beta0 - kappa * delta * d_beta0;
    end
    d_beta = X' * (Xbeta + beta0 - y) / n - D' * d_gamma;
    beta = beta - kappa * delta * d_beta;
    %% update var_hist and var_order %%
    if opt.auc 
        gamma_index = find(gamma);
        if sum(~ismember(gamma_index,var_hist)) > 0
            var_order = [var_order;step_cur * ones(sum(~ismember(gamma_index,var_hist)),1)];
            var_hist = [var_hist;gamma_index(~ismember(gamma_index,var_hist))];
        end
    end
    %% Recording some of estimations in the regularization path %%
    while true
        dt = step_cur * delta + t0 - opt.t_seq(rec_cur);
        if dt < 0, break; end
        %% update \beta(0),\beta,z and \gamma %%
        if opt.intercept
            obj.beta0(:, rec_cur) = beta0 + kappa * dt * d_beta0;
        end
        obj.z(:, rec_cur) = z + dt * d_gamma;
        obj.gamma(:, rec_cur) = kappa * sign(z + dt * d_gamma) .* max(abs(z + dt * d_gamma) - 1, 0);
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