% DESCRIPTION
% This file is the application of split lbi to university dataset.
% For details please see the paper.
%% Data preprocessing %%
load('Data_bas.mat');
Compare = Data_bas.Compare;
Team = Data_bas.Team;
p = length(Team);
n = size(Compare,1);
X = zeros(n,p);
index_left = (Compare(:,1) - 1) * n + [1:n]';
index_right = (Compare(:,2) - 1) * n + [1:n]';
X(index_left) = 1;
X(index_right) = -1;
y = Compare(:,3);

%% Implementing split lbi %%
data.X = X;
data.y = y;
eig_X = eig(X'*X);
data.D = X / sqrt(min(eig_X(abs(eig_X)>=0.001)));
clear opt;
opt.nu = 1;
opt.kappa = 100;
opt.intercept = false;
opt.t_ratio = 20;
opt.fast_init = false;
family = 'gaussian';
group_type = 'ungrouped';
obj = splitlbi.splitlbi(data,opt,family,group_type);

%% plot %%
beta_tilde = obj.beta_tilde;
beta_tilde = beta_tilde(:,1:100);
t_seq = lbi.t_seq;
t_seq = t_seq(1:100);
t_tick = [];
ind_t = 55;
logscale = 1;
fig_option.fig_scale = 1;
fig_option.xlabel = 't';
fig_option.ylabel = 'Regularization Path';
fig_option.title = 'Basketball';
fig_option.legend_list = Team;
name = 'basketball';
splitlbi.plot_splitlbi(beta_tilde,t_seq,logscale,t_tick,fig_option,name)

%% See grouping result %%
score = beta_tilde(:,ind_t);
splitlbi.table_grouping(score,Team')
