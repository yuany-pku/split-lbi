% DESCRIPTION
% This file is setting for n,p that will be used for 
% simulation under logit model
n = 100; p = 80;
S1 = (1:floor(p / 20));
S2 = ((floor(p / 20) + 1):floor(p / 10));
S = [S1, S2];
beta_true = zeros(p, 1);
beta_true(S1) = 2; beta_true(S2) = -2;
T = find(~beta_true);
s = length(S); t = length(T);

