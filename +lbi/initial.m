function opt_out = initial(opt_in)
% Damping factor (kappa)
opt_out.kappa = 10;
% Step size (delta)
opt_out.delta = [];
% Factor when calculating step size (delta) from data
opt_out.c = 1;
% Sequence of recorded times
opt_out.t_seq = [];
% Number of recorded times
opt_out.t_num = 100;
% t_max/t_min
opt_out.t_ratio = 100;
% Normalizing the predictors?
opt_out.normalize = false;
% Having intercept?
opt_out.intercept = true;
% Directly setting the intercept to be the MLE before t0?
opt_out.fast_init = false;
% Compute AUC or not?
opt_out.auc = false;
if nargin == 0
    if nargout == 0
        disp('Default options:');
        disp(opt_out);
    end
    return;
end

% Substitute opt_in fields for the default ones
fields = fieldnames(opt_in);
for i = 1:length(fields)
    field = fields{i};
    if isfield(opt_out, field);
        opt_out.(field) = opt_in.(field);
    end
end
end