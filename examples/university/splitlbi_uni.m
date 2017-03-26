% DESCRIPTION
% This file is the application of split lbi to university dataset.
% For details please see the paper.
%% Load Data %%
load('univer.mat');
compare = univer.uni_robust;
uni_list = unique(compare(:,[1,2]));
n = size(compare,1);
p = length(uni_list);
load('UniName.mat');
ind_t = 313;

%% Define X,y and D for shrinked dataset %%
X = zeros(n,p);
row_X = 1:n;
X((compare(:,1) - 1) * n + row_X') = 1;
X((compare(:,2) - 1) * n + row_X') = -1;
y = compare(:,3);
[index_l,index_r] = vec2pair(p);
row_D = 1:length(index_l);
D = zeros(length(index_l),p);
D((index_l - 1) * length(index_l) + row_D') = 1;
D((index_r - 1) * length(index_r) + row_D') = -1;
X = sparse(X);
D = sparse(D);
D = D / sqrt(p);


%% Implementing split lbi using shrinked dataset  %%
data.D = D;
data.X = X;
data.y = y;
clear opt;
opt.kappa = 20;
opt.t_ratio = 100;
opt.fast_init = false;
opt.intercept = false;
opt.t_num = 400;
family = 'gaussian';
obj = splitlbi.splitlbi(data,opt,family);
beta_tilde = obj.beta_tilde;
t_seq = obj.t_seq;

%% plot sub-graph %%
uni_index = [1,9,16,33,259];
beta_sub = beta_tilde(uni_index,:);
t_tick = [];
logscale = 0;
fig_option.fig_scale = 1;
fig_option.xlabel = 't';
fig_option.ylabel = 'Regularization Path';
fig_option.title = 'University';
fig_option.legend_list = uni_name(uni_index);
name = 'university';
splitlbi.plot_splitlbi(beta_sub,t_seq,logscale,t_tick,fig_option,name);

%% See grouping result %%
score = beta_tilde(:,ind_t);
splitlbi.table_grouping(score,uni_name)

