function beta_tilde = proj_Dt(gamma, beta, D)
% DESCRIPTION
%  This file is used to compute \beta_tilde by projecting
%  \beta onto the subspace of signal set of D.
p = size(D,2);
gamma(abs(gamma) < 1e-12) = 0;
tmp = D(gamma==0, :); 
tmp = sparse(tmp);
tmp = tmp' * tmp;
beta_tilde = (eye(p) - pinv(full(tmp)) * tmp) * beta;
beta_tilde(abs(beta_tilde)<=1e-10) = 0;
end
