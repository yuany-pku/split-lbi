% DESCRIPTION
% This file is setting for n,p that will be used for 
% simulation with group sparisty under logit model
n = 100; p = 80;
S1 = (1:floor(p / 20));
S2 = ((floor(p / 20) + 1):floor(p / 10));
S = [S1, S2];
beta_true = zeros(p, 1);
beta_true(S1) = 2; beta_true(S2) = -2;
T = find(~beta_true);
s = length(S); t = length(T);
group_index = zeros(p,1);
group_index(S1) = 1;
S_t = S1;
for i=1:p/length(S1) - 1
    S_t = S_t + length(S1);
    group_index(S_t) = i + 1;
end