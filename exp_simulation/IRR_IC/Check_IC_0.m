function IC = Check_IC_0(X, D, I)
% DESCRIPTION
% Compute IC_0 of genlasso
if nargin < 3
    error('The input is incomplete!');
else
    [~,p] = size(X);
    if size(D,2) ~= p, error('The number of columns of D and X should be equal!'); end
end

XTX = X' * X;
DI = D(I,:);
DJ = D;
DJ(I,:) = [];

[V,S,U] = svd(DJ);
ind = find(diag(S)>1e-12);
U(:,ind) = [];
V(:,ind) = [];
if (length(ind)==p)
    error('Full rank D_J');
end
A = U * pinv(U' * XTX * U) * U';
Omega_s = pinv(DJ')*((XTX * A - eye(p))*(DI'));
IC = norm(Omega_s,'inf');
end