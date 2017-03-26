% DESCRIPTION
% This file is setting for n,p that will be used for 
% simulation under linear model
n = 50; p = 50; 
s = 15;
S1 = floor(p/4) - (0:(s/3-1));
S2 = floor(p/2) - (0:(s/3-1));
S3 = floor(3*p/4) - (0:(s/3-1));
S = [S1, S2, S3];
beta_true = zeros(p, 1);
beta_true(S1) = - 2; beta_true(S2) = 2;beta_true(S3) = - 2;
T = find(~beta_true);
s = length(S); t = length(T);
sigma = 1;

