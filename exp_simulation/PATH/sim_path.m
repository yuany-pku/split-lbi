% DESCRIPTION
% This file plot regularization solution path of split lbi
% for different \nu
%% Setting %%
run ../Linear_setting_sim.m;
X = load(strcat('../Data/X',num2str(1),'_lbi_linear.txt'));
y = load(strcat('../Data/y',num2str(1),'_lbi_linear.txt'));
D = [eye(p - 1), zeros(p - 1, 1); eye(p)] - [zeros(p - 1, 1), eye(p - 1); zeros(p)];
d_true = D * beta_true;
data.y = y;
data.X = X;
data.D = D;
clear opt;
%% Plot Path of \gamma vs t for split lbi %%
nulist = [1,5,10];
for nu_num = 1:length(nulist)
    nu = nulist(nu_num);
    opt.nu = nu;
    opt.kappa = 200;
    opt.intercept = false;
    opt.t_ratio = 20;
    opt.fast_init = false;
    opt.t_num = 1000;
    obj = splitlbi.splitlbi(data,opt);
    gamma = obj.gamma;
    t_seq = obj.t_seq;
    gamma_0 = [zeros(size(gamma,1),1),gamma];
    t_seq_0 = [0,t_seq];
    [~,ind_stop] = min(abs(t_seq - 30));
    gamma_path = gamma_0(:,1:ind_stop);
    t_seq_path = t_seq_0(:,1:ind_stop);
    t_tick = [];
    logscale = 0;
    fig_option_path.fig_scale = 1;
    fig_option_path.xlabel = 't';
    fig_option_path.ylabel = 'Coordinates of $$\gamma$$';
    fig_option_path.title = strcat('$$\nu = $$',num2str(nu));
    name = strcat('path_','nu_',num2str(nu));
    splitlbi.plot_splitlbi(gamma_path,t_seq_path,logscale,t_tick,fig_option_path,name)
end
