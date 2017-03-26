function [] = plot_splitlbi(beta,t_seq,logscale,t_tick,fig_option,name,fig_path,varargin)
% DESCRIPTION
% This file plots regularization solution path of split lbi and lbi
if nargin < 1 || isempty(beta)
    error('Estimators in the path is missing!'); 
else
    [p,t_num] = size(beta);
    beta = beta - repmat(mean(beta),p,1);
end
if nargin < 2 || isempty(t_seq)
    t_seq = logspace(log10(1),log10(1*100),t_num);
end
if nargin < 3 || isempty(logscale), logscale = false; end
if logscale == true
    t_seq = log10(t_seq);
    t_tick = log10(t_tick);
end
if nargin < 4 || isempty(t_tick), t_tick = [];end
if nargin < 5 || isempty(fig_option) 
    fig_scale = 1; 
else
    if isfield(fig_option,'fig_scale')
        fig_scale = fig_option.fig_scale;
    end
end
    
        
figure;
h = plot(repmat(t_seq,p,1)',beta',varargin{:});
xlim_range = [1.02, -0.02] * min(t_seq) + [-0.02, 1.02] * max(t_seq);
xtick = t_tick;
set(gca, 'FontSize', 10 * fig_scale);
set(h, 'MarkerSize', 6 * fig_scale);
set(h,'LineStyle','--');
if length(unique(xlim_range)) > 1, set(gca, 'xlim', xlim_range); end
set(gca, 'xtick', xtick);
if isfield(fig_option,'xlabel')
    xlab_name = fig_option.xlabel;
    xlabel(xlab_name,'Interpreter','latex');
end
if isfield(fig_option,'ylabel')
    ylab_name = fig_option.ylabel;
    ylabel(ylab_name,'Interpreter','latex');
end
if isfield(fig_option,'title')
    tit_name = fig_option.title;
    title(tit_name,'Interpreter','latex');
end
if isfield(fig_option,'legend_list')
    legend_list = fig_option.legend_list;
    legend(h,legend_list,'Location','Best');
end
if nargin >= 6 && ~isempty(name)
    if nargin < 7 || isempty(fig_path) 
        fig_path = './'; 
    end
    saveas(gcf, [strcat(fig_path,name), '.png']);
end

    